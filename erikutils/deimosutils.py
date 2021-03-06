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

    return sky_coords, desislits, bluslits, objcats
    """
    with fits.open(maskfn) as f:
        dblu = f['BluSlits'].data
        xorder = np.argsort(dblu['slitX1'])

        ddesi = f['DesiSlits'].data
        scsinorder = SkyCoord(ddesi['slitRa']*u.deg, ddesi['slitDec']*u.deg)[xorder]

        somdata = f['SlitObjMap'].data
        cdata = f['ObjectCat'].data
        cats = cdata[np.searchsorted(cdata['ObjectId'], somdata['ObjectId'][xorder])]

        return scsinorder, ddesi[xorder], dblu[xorder], cats


def plot_deimos_slits(maskfn, ax=None, plotkwargs={}, onsky=True):
    """
    Plot all slits from a DEIMOS bintab or design file.

    Parameters
    ----------
    maskfn : str
        The name of the file to load from
    ax : matplotlib axes or None
        The axes to plot into or None to use the current axes
    plotkwargs : dict
        keywords to pass into the plot command used to make the boxes
    onsky : bool
        If True, plot the *intended* slits (in on-sky degrees RA/Dec).  If
        False, plots the slits as dsimulator outputted them on the mask.

    Returns
    -------
    same as `get_deimos_slits_inorder`
    """
    from matplotlib import pyplot as plt

    if ax is None:
        ax = plt.gca()

    res = get_deimos_slits_inorder(maskfn)
    scs, ddesi, dblu, cats = res

    if onsky:
        for drow in ddesi:
            ra0 = drow['SLITRA'] * u.deg
            dec0 = drow['SLITDEC'] * u.deg
            w = drow['SLITWID'] * u.arcsec
            l = drow['SLITLEN'] * u.arcsec
            pa = drow['SLITLPA'] * u.deg
            spa = np.sin(pa)
            cpa = np.cos(pa)

            xs_unrot = [-l/2, -l/2, +l/2, +l/2]
            xs_unrot.append(xs_unrot[0])
            xs_unrot = u.Quantity(xs_unrot)
            ys_unrot = [-w/2, +w/2, +w/2, -w/2]
            ys_unrot.append(ys_unrot[0])
            ys_unrot = u.Quantity(ys_unrot)

            xs = xs_unrot*cpa - ys_unrot*spa
            ys = xs_unrot*spa + ys_unrot*cpa

            plt.plot(xs + ra0, ys + dec0, **plotkwargs)

    else:
        for brow in dblu:
            xs = [brow['SLITX1'], brow['SLITX2'], brow['SLITX3'], brow['SLITX4'], brow['SLITX1']]
            ys = [brow['SLITY1'], brow['SLITY2'], brow['SLITY3'], brow['SLITY4'], brow['SLITY1']]
            plt.plot(xs, ys, **plotkwargs)


    return res


def get_deimos_spec1d(fn, horne=False, smoothing=False, fixfile=True):
    """
    Gets a reduced DEIMOS spec1d file's contents

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
    fixfile : bool
        Silently fix any fits errors encountered



    Returns
    -------
    hdus
        hdus for the selected spectrum type
    xb
        wl array for the blue side
    bspec
        spectrum for the blue side
    xr
        wl array for the red side
    rspec
        spectrum for the red side
    bivar
        inverse variance for the blue side
    rivar
        inverse variance for the red side

    """
    from scipy import signal


    with fits.open(fn) as f:
        if fixfile:
            f.verify('silentfix')

        hdub = f[1+int(horne)*2]
        hdur = f[2+int(horne)*2]
        db = hdub.data
        dr = hdur.data

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

        return [hdub, hdur], db['LAMBDA'][0], bspec, dr['LAMBDA'][0], rspec, db['IVAR'][0], dr['IVAR'][0]

def plot_deimos_spec1d(fn, horne=False, smoothing=False, catv=None, mady=False, fixfile=True):
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
    fixfile : bool
        Silently fix any fits errors encountered



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


    hdus, xb, bspec, xr, rspec, bivar, rivar = get_deimos_spec1d(fn, horne, smoothing, fixfile)

    plt.step(xb, bspec, color='b', where='mid')
    plt.step(xr, rspec, color='r', where='mid')

    plt.plot(xb, bivar**-0.5, color='k', ls=':')
    plt.plot(xr, rivar**-0.5, color='k', ls=':')

    if catv is not None:
        for wl in lineswl:
            plt.axvline((wl*(1+catv/cnst.c)).to(u.angstrom).value, c='k')

    plt.xlabel(r'$\lambda [{\rm \AA}]$ ')
    plt.title(fn)

    if mady:
        specbr = np.concatenate((bspec, rspec))
        mad = median_absolute_deviation(specbr)
        med = np.median(specbr)

        uppery = med + mad*mady
        lowery = med - mad*mady
        if lowery > 0:
            lowery = 0
        plt.ylim(lowery, uppery)

    return xb, bspec, xr, rspec


def show_deimos_slit(spec2dfn, madcolorscale=None, scalekwargs=None, ax=None, fixfile=True):
    """
    Show a deimos 2d slitlet from the spec2d file
    """
    from matplotlib import pyplot as plt
    from astropy.stats import median_absolute_deviation
    from astropy.visualization import scale_image

    if ax is None:
        ax = plt.gca()

    with fits.open(spec2dfn) as f:
        if fixfile:
            f.verify('silentfix')

        h = f[1].header
        d = f[1].data

        xsz, ysz = map(int, h['TDIM1'][1:-1].strip().split(','))

        dflux = d['FLUX'].reshape((ysz, xsz))

    if scalekwargs:
        dflux = scale_image(dflux, **scalekwargs)

    if madcolorscale is None:
        vmin = vmax = None
    else:
        if madcolorscale is not None:
            raise ValueError('cannot rescale if scale_image used')
        mad = median_absolute_deviation(dflux)
        med = np.median(dflux)
        vmin = med-mad*float(madcolorscale)
        vmax = med+mad*float(madcolorscale)

    plt.title(spec2dfn)
    plt.imshow(dflux, interpolation='nearest', vmin=vmin, vmax=vmax)
    if scalekwargs is None:
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
