from __future__ import division, print_function
"""
Various helpful snippets by Erik Tollerud (erik.tollerud@gmail.com)
"""

import numpy as np

__all__ = ['write_to_mac_clipboard', 'show_imagelist']


def write_to_mac_clipboard(s):
    from subprocess import Popen, PIPE

    p = Popen('pbcopy', stdin=PIPE)
    return p.communicate(s)


SDSS_IMAGE_LIST_URL = 'http://skyserver.sdss3.org/dr10/en/tools/chart/list.aspx'
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
                       baseurl='http://skyserver.sdss3.org/dr12/en/tools/chart/navi.aspx?ra={ra}&dec={dec}',
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
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astropy.io import fits

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
    from astropy.io import ascii

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
