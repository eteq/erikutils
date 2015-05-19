#!/usr/bin/env python
"""
Various utility functions for working with spectra from DEIMOS on Keck II
"""

import os

import numpy as np
from astropy import units as u
from astropy import constants as cnst
from astropy.io import fits
from astropy.coordinates import SkyCoord


DEIMOS_pixscale = .1185 * u.arcsec/u.pixel
DEIMOS_sizescale = .727 * u.mm/u.arcsec


def get_deimos_slits_inorder(maskfn):
    """
    Takes slits from a DEIMOS dsimulator fits file and returns them sorted in
    order of increasing x.  This should be the same orientation as they appear
    on the DEIMOS science images (so the far-left/lowest x object is the 0th
    entry in the returned lists)

    return scs, desislits, objcats
    """
    f = fits.open(maskfn)
    xorder = np.argsort(f['BluSlits'].data['slitX1'])

    ddesi = f['DesiSlits'].data
    scsinorder = SkyCoord(ddesi['slitRa']*u.deg, ddesi['slitDec']*u.deg)[xorder]

    somdata = f['SlitObjMap'].data
    cdata = f['ObjectCat'].data
    cats = cdata[np.searchsorted(cdata['ObjectId'], somdata['ObjectId'][xorder])]

    return scsinorder, ddesi[xorder], cats


def plot_deimos_spec1d(fn, horne=False, smoothing=False, catv=None, mady=False):
    """
    Plots a reduced DEIMOS spec1d file with matplotlib

    Parameters
    ----------
    fn : str
        The spec1d file to load
    horne : bool
        If True, use the horne extraction, otherwise boxcar
    smoothing : float or False
        If not False/0, smooth the spectrum.  If positive, the value gives the
        number of pixels of boxcar smoothing to use.  If negative, gives the
        size of a gaussian kernel to use for smoothing.
    catv : None or astropy Quantity with velocity units
        If not None, shows Halpha and the calcium triplet assuming the star is
        moving at the provided velocity (use E.g. ``0*u.km/u.s`` for rest
        wavelength).
    mady : float or False
        If not False/0, re-scale the part of the spectrum shown to ``mady``
        median absolute deviations from the median



    Returns
    -------
    xb
        wl array for the blue side
    bspec
        spectrum for the blue side
    xr
        wl array for the red side
    rspec
        spectrum for the red side

    """
    from astropy.stats import median_absolute_deviation
    from matplotlib import pyplot as plt
    from scipy import signal

    # Halpha + CaT
    lineswl = [6562.801, 8498.03, 8542.09, 8662.14]*u.angstrom

    with fits.open(fn) as f:
        db = f[1+int(horne)*3].data
        dr = f[2+int(horne)*3].data

        bspec = db['SPEC'][0]
        rspec = dr['SPEC'][0]

        if smoothing:
            if smoothing < 0:
                kernelb = signal.gaussian(len(bspec), -smoothing)
                kernelr = signal.gaussian(len(rspec), -smoothing)
            else:
                kernelr = kernelb = [1/float(smoothing)] * int(smoothing)
            bspec = np.convolve(bspec, kernelb, 'spec')
            rspec = np.convolve(rspec, kernelr, 'spec')

        plt.step(db['LAMBDA'][0], bspec, color='b', where='mid')
        plt.step(dr['LAMBDA'][0], rspec, color='r', where='mid')

        plt.plot(db['LAMBDA'][0], db['IVAR'][0]**-0.5, color='k', ls=':')
        plt.plot(dr['LAMBDA'][0], dr['IVAR'][0]**-0.5, color='k', ls=':')

        if catv is not None:
            for wl in lineswl:
                plt.axvline((wl*(1+catv/cnst.c)).to(u.angstrom).value, c='k')

        plt.xlabel(r'$\lambda [{\rm \AA}]$ ')
        plt.title(fn)

        if mady:
            specbr = np.concatenate((db['SPEC'][0], dr['SPEC'][0]))
            mad = median_absolute_deviation(specbr)
            med = np.median(specbr)

            uppery = med + mad*mady
            lowery = med - mad*mady
            if lowery > 0:
                lowery = 0
            plt.ylim(lowery, uppery)

        return db['LAMBDA'][0], bspec, dr['LAMBDA'][0], rspec


def plot_deimos_slit(fn, madcolorscale=None, scalekwargs=None):
    from matplotlib import pyplot as plt
    from astropy.stats import median_absolute_deviation
    from astropy.visualization import scale_image

    with fits.open(fn) as f:
        h = f[1].header
        d = f[1].data

        xsz, ysz = map(int, h['TDIM1'][1:-1].strip().split(','))

        dflux = d['FLUX'].reshape((ysz, xsz))

    if scalekwargs:
        dflux = scale_image(dflux, **scalekwargs)

    if madcolorscale is None:
        vmin = vmax = None
    else:
        mad = median_absolute_deviation(dflux)
        med = np.median(dflux)
        vmin = med-mad*float(madcolorscale)
        vmax = med+mad*float(madcolorscale)

    plt.title(fn)
    plt.imshow(dflux, interpolation='nearest', vmin=vmin, vmax=vmax)
    plt.colorbar(orientation='horizontal')


_mmperpx = (DEIMOS_pixscale * DEIMOS_sizescale).to(u.mm/u.pixel).value
def deimos_slits_to_xpx(mskorbintabfn):
    """
    Takes a DEIMOS image or bintab (output from the spec2d pipeline) and gives
    the x-pixel of each of the BluSlits in order
    """
    if mskorbintabfn.endswith('.fits'):
        fn = mskorbintabfn
    else:
        fn = '%s.bintabs.fits' % mskorbintabfn
        if not os.path.exists(fn):
            fn = mskorbintabfn+os.sep+fn
    blu = fits.getdata(fn, 'BluSlits')
    xmidmm = (blu.SLITX1+blu.SLITX2+blu.SLITX3+blu.SLITX4)/4.
    # centered on 0

    px = xmidmm/_mmperpx + 1024*4
    chipoffsets = [87, -13, -122, 150]  # crappily empirical - must be chip gap
    chipmsks = px<2048, (2048<px)&(px<4096), (4096<px)&(px<6144), 6144<px

    for i in range(4):
        if sum(chipmsks[i])>0:
            px[chipmsks[i]] = px[chipmsks[i]]+chipoffsets[i]

    return px
