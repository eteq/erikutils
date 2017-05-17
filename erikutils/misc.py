from __future__ import division, print_function
"""
Various random snippets by Erik Tollerud (erik.tollerud@gmail.com)
"""

import numpy as np
from astropy import units as u

from astropy.io import ascii
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


def get_sdss_cutout_url(coord, size=204.8*u.arcsec, pixelscale=0.4*u.arcsec/u.pixel, opts='',
                        baseurl='http://skyservice.pha.jhu.edu/DR12/ImgCutout/getjpeg.aspx'):
    """
    Gets the URL to get an SDSS image cutout

    Parameters
    ----------
    coord : SkyCoord
        The coordinate to view
    size : astropy quantity (scalar or length-2) with units of arcsec or pixels
        The size of the image.  If a single value, assumes square.  If two
        values, they are (width, height)
    pixelscale : quantity in angle/pixel or pixel/angle
        The pixel scale for the output images.
    opts : str
        The 'opt' to supply to the cutout service (see SDSS docs for what they mean).
    baseurl : str
        The URL to use for the cutout service (without the query part).


    Returns
    -------
    url : str
        The url for the image cutout.

    """
    from six.moves.urllib_parse import urlencode

    if pixelscale.unit.is_equivalent(u.pixel/u.arcsec):
        pixelscale = 1/pixelscale
    pixscaleapp = pixelscale.to(u.arcsec/u.pixel).value

    if size.unit.is_equivalent(u.arcsec):
        size = size / pixelscale

    if size.isscalar:
        wpix = hpix = size.to(u.pixel).value
    elif len(size) == 2:
        wpix, hpix = size.to(u.pixel).value
    else:
        raise ValueError('size must be a scalar or length-2 quantity')

    qrydct = {'ra': coord.ra.deg, 'dec': coord.dec.deg,
              'width': int(wpix), 'height': int(hpix),
              'scale': pixscaleapp, 'opt': opts}
    url = baseurl + '?' + urlencode(qrydct)

    return url


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
    from six.moves.urllib_parse import urlencode
    import webbrowser

    url = baseurl.format(ra=sc.ra.degree, dec=sc.dec.degree)
    if kwargs:
        url = url + '&' + urlencode(kwargs)
    if openinbrowser:
        webbrowser.open(url)
    return url


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


def skycoord_to_regions(scs, shape='point', otherparams='', ds9=None):
    """
    Add regions to an open ds9 from a SkyCoord.

    For example, to get red 2" circles:

        skycoord_to_regions(scs, 'circle', '2" # color = red')

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

def get_distutils_compiler_info():
    """
    Like all good things, this is from Stack Overflow:
    http://stackoverflow.com/questions/7924420/how-to-get-the-actual-command-that-the-distutils-compile-function-uses/28849344#28849344
    """

    import distutils.sysconfig
    import distutils.ccompiler
    compiler = distutils.ccompiler.new_compiler()
    distutils.sysconfig.customize_compiler(compiler)
    return compiler, compiler.compiler_so

def offset_info(fromsc, tosc):
    """
    Compute offsets to get from one SkyCoord to another

    Parameters
    ----------
    fromsc : SkyCoord
        The "from" coordinate
    tosc : SkyCoord
        The "to" coordinate

    Returns
    -------
    str
        A string with the offset information
    """
    dra = (tosc.ra-fromsc.ra)*np.cos((tosc.dec + fromsc.dec)/2)
    ddec = tosc.dec-fromsc.dec
    eorw = 'E' if dra>0 else 'W'
    nors = 'N' if ddec>0 else 'S'
    return '{0} {1} {2} {3}'.format(dra.arcsec, eorw, ddec.arcsec, nors)
