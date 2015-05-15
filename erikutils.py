from __future__ import division, print_function
"""
Various helpful snippets by Erik Tollerud (erik.tollerud@gmail.com)
"""

import numpy as np
from astropy import units as u
from astropy import constants as cnst

from astropy.io import fits, ascii
from astropy.coordinates import SkyCoord, Angle


def write_to_mac_clipboard(s):
    from subprocess import Popen, PIPE

    p = Popen('pbcopy', stdin=PIPE)
    return p.communicate(s)


SDSS_IMAGE_LIST_URL = 'http://skyserver.sdss.org/dr12/en/tools/chart/list.aspx'
_imglist_post_templ = """
<html>
<head>
<title>SDSS sampled_imagelist form</title>
</head>
<body>
<h1> SDSS sampled_imagelist form </h1>
<p>Using URL {url}</p>

<form action="{url}"
method="post">
<TEXTAREA name="paste">
{text}
</TEXTAREA>
<br>
Scale: <input class="in" type="text" value="0.4" name="scale">
<br>
Opt: <input class="in" type="text" value="" name="opt">
<br>
<input type="submit">
</form>
</body>
</html>
"""
def show_imagelist(scs, names=None, url=SDSS_IMAGE_LIST_URL, posttoimglist=3.):
    """
    Returns the text to be pasted into the sdss image list page.  Also opens
    the page (if `url` is not None) and copies the text to the clipboard if on
    a mac or linux.

    Parameters
    ----------
    ras : SkyCoord
        SkyCoord of objects to show
    url : str or None
        The URL to the SDSS image list page or None to not open in a web
        browser.
    copytoclipboard : bool
        If True, copies the list of images to the clipboard for use on the SDSS
        web site
    posttoimglist : bool or float
        If True, makes a form to post to the URL site. If a float, gives the
        number of seconds to wait until deleting the temporary file (to gives
        the browser time to load).

    Returns
    -------
    text : str
        The table to be pasted into the image list text box

    """
    import webbrowser
    import tempfile
    import time
    import os

    if hasattr(scs, 'ra') and hasattr(scs, 'dec'):
        decs = scs.dec
        ras = scs.ra
    else:
        # assume it's a Table
        decs = ras['dec']
        ras = ras['ra']

    if len(ras) != len(decs):
        raise ValueError('ras and decs not the same size!')

    ras = np.array(ras, copy=False)
    decs = np.array(decs, copy=False)

    if names is None:
        names = [str(i) for i in range(len(ras))]

    text = ['name ra dec']
    for nmi, rai, deci in zip(names, ras, decs):
        text.append('{0} {1} {2}'.format(nmi, rai, deci))
    text = '\n'.join(text)

    if url:
        if posttoimglist:
            page = _imglist_post_templ.format(url=url, text=text)
            tf = tempfile.NamedTemporaryFile(delete=False)
            tf.write(page)
            tf.flush()
            fiurl = 'file://' + os.path.abspath(tf.name)
            webbrowser.open(fiurl)
            if isinstance(posttoimglist, float):
                time.sleep(posttoimglist)
        else:
            webbrowser.open(url)

    return text


def show_skycoord_sdss(sc, openinbrowser=True,
                       baseurl='http://skyserver.sdss.org/dr12/en/tools/chart/navi.aspx?ra={ra}&dec={dec}',
                       **kwargs):
    """
    Shows a `SkyCoord` in the SDSS Navigate interface.

    Parameters
    ----------
    sc : SkyCoord
        The coordinate to view
    openinbrowser : bool
        If True, opens the page in a web browser
    baseurl : str
        The URL to use, with '{ra}' and '{dec}' indicating the coordinates in
        degrees
    kwargs
        Additional string keywords to add to the query

    Returns
    -------
    url : str
        The url to show
    """
    from urllib import urlencode
    import webbrowser

    url = baseurl.format(ra=sc.ra.degree, dec=sc.dec.degree)
    if kwargs:
        url = url + '&' + urlencode(kwargs)
    if openinbrowser:
        webbrowser.open(url)
    return url


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


def to_clipboard(data):
    import subprocess

    p = subprocess.Popen('pbcopy', stdin=subprocess.PIPE)
    p.communicate(data)
    return p.wait()


def keola_html_to_ascii_log(htmlfn, outfn=None):
    """
    Takes an html output like what comes out when you hit the "save to disk"
    button on the Keck e-log tool, and outputs an ASCII-formatted file with the
    observing log info.

    Currently this doesn't output the night comments or weather data, because
    that requires more complex parsing of the html.  Maybe later I'll do that.
    """
    if not htmlfn.endswith('.html'):
        raise ValueError('need the html log from KeOLA')
    if outfn is None:
        outfn = htmlfn[:-5] + '.log'

    header = ascii.read(htmlfn, format='html', htmldict={'table_id': 1})
    observations = ascii.read(htmlfn, format='html', htmldict={'table_id': 2})

    with open(outfn, 'w') as f:
        header.write(f, format='ascii.no_header')
        f.write('\n\n')
        observations.write(f, format='ascii')

    return outfn


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


def regions_to_skycoord(ds9=None):
    """
    Gets the regions from ds9 and returns their locations as a skycoord

    Returns sc, names, others
    """
    import pyds9
    import pyregion

    if ds9 is None:
        ds9 = pyds9.DS9()

    oldsys = ds9.get('regions system')
    oldsf = ds9.get('regions skyformat')
    try:
        ds9.set('regions system wcs')
        ds9.set('regions skyformat degrees')
        regstr = ds9.get('regions')
    finally:
        ds9.set('regions system ' + oldsys)
        ds9.set('regions skyformat ' + oldsf)

    frame = ds9.get('regions sky')
    regs = pyregion.parse(regstr)

    ras = []
    decs = []
    names = []
    others = []
    for s in regs:
        ras.append(s.params[0].v*u.deg)
        decs.append(s.params[1].v*u.deg)

        names.append(s.name)

        other = []
        for p in s.params[2:]:
            if hasattr(p, 'degree'):
                other.append(Angle(p.degree, u.deg))
            else:
                other.append(p.v)
        others.append(tuple(other))

    scs = SkyCoord(ras, decs, unit=u.deg, frame=frame)
    return scs, names, others


def skycoord_to_regions(scs, shape, otherparams, ds9=None):
    """
    Add regions to an open ds9 from a SkyCoord

    Parameters
    ----------
    scs : SkyCoord
        The coordinate(s) to place a region from
    shape : str
        The name of the type of shape (goes straight into the region file)
    otherparams : str or list of str
        The params after the coordinates.  If a list, must match the length of
        `scs`.  If a str, the same string will be used for all.
    ds9 : pyds9.DS9 or None
        The ds9 instance to use.  If None, will try to create one.
    """
    import pyds9

    if ds9 is None:
        ds9 = pyds9.DS9()

    if scs.isscalar:
        scs = SkyCoord([scs])

    reglines = []
    if isinstance(otherparams, basestring):
        for ra, dec in zip(scs.ra.deg, scs.dec.deg):
            reglines.append('icrs; {shape} {ra}d {dec}d '.format(**locals()) + otherparams)
    else:
        if len(otherparams) != len(scs):
            raise ValueError('otherparams must match scs')
        for ra, dec, otherparam in zip(scs.ra.deg, scs.dec.deg, otherparams):
            reglines.append('icrs; {shape} {ra}d {dec}d '.format(**locals()) + otherparam)

    ds9.set('regions', '\n'.join(reglines))
