from __future__ import division, print_function
"""
A variety of scripts for orchestrating reduction of HST/ACS imaging.
"""

import os
import shutil
from warnings import warn

import numpy as np

from astropy.io import fits
from astropy import units as u

from dolphot_runner import DolphotRunner

class Association(object):
    def __init__(self, asnfn):
        self.datadir = os.path.abspath(os.path.split(asnfn)[0])
        asn = fits.open(asnfn)

        self.asn_header = asn[0].header
        self.table = asn[1].data

        self._populate_asn_attributes()
        self._set_flt_info()

    def _populate_asn_attributes(self):
        """
        Use the header and table to figure out what
        """
        from astropy.coordinates import SkyCoord

        hdr = self.asn_header

        self.target_name = hdr['TARGNAME']
        self.target_coords = SkyCoord(hdr['RA_TARG']*u.deg, hdr['RA_TARG']*u.deg)
        self.instrument = hdr['INSTRUME']
        self.detector = hdr['DETECTOR']
        self.propid = hdr['PROPOSID']
        self.data = hdr['DATE']

        self.exposure_names = []
        self.product_name = None

        for nm, typ, cal in self.table:
            if not cal:
                raise ValueError('File {0} was not calibrated!'.format(nm))
            if typ == 'EXP-DTH':
                self.exposure_names.append(nm.lower())
            elif typ == 'PROD-DTH':
                if self.product_name:
                    raise ValueError('Found *two* products: "{0}" and '
                                     '"{1}"'.format(self.product_name, nm))
                self.product_name = nm.lower()
            else:
                raise ValueError('Unrecognized type "{0}" for file {1}'.format(typ, nm))


    def _set_flt_info(self):
        """
        Determines extra info from the flt of the first exposure
        """
        try:
            flfn = self.flts[0]
            self.fl_header = hdr = fits.getheader(flfn, 0)
        except IOError:
            flfn = self.flcs[0]
            self.fl_header = hdr = fits.getheader(flfn, 0)


        f1, f2 = (hdr['FILTER1'], hdr['FILTER2'])
        if f1.lstrip().lower().startswith('clear'):
            self.filter = f2.strip()
        elif f2.lstrip().lower().startswith('clear'):
            self.filter = f1.strip()
        else:
            self.filter = f1.strip() + '-' + f2.strip()

    @property
    def flts(self):
        return [os.path.join(self.datadir, basefn + '_flt.fits') for basefn in self.exposure_names]

    @property
    def flcs(self):
        return [os.path.join(self.datadir, basefn + '_flc.fits') for basefn in self.exposure_names]

    @property
    def drz(self):
        return os.path.join(self.datadir, self.product_name + '_drz.fits')

    @property
    def drc(self):
        return os.path.join(self.datadir, self.product_name + '_drc.fits')


def find_target_associations(datadir, targrenames=None):
    from glob import glob

    asnfns = glob(os.path.join(datadir, '*asn*.fits'))
    asns = [Association(fn) for  fn in asnfns]
    targdct = dict()
    for asn in asns:
        targ_asns = targdct.setdefault(asn.target_name, [])
        targ_asns.append(asn)

    if targrenames:
        for old, new in targrenames.items():
            targdct[new] = targdct.pop(old)

    return targdct, asns

default_dolphot_params = {
'img_apsky': '15 25',
'UseWCS': '1',
'RAper': '4',
'RChi': '2.0',
'RSky0': '15',
'RSky1': '35',
'SkipSky': '2',
'SkySig': '2.25',
'SecondPass': '5',
'SigFindMult': '0.85',
'MaxIT': '25',
'NoiseMult': '0.10',
'FSat': '0.999',
'ApCor': '1',
'RCentroid': '2',
'PosStep': '0.25',
'dPosMax': '2.5',
'RCombine': '1.5',
'RPSF': '10',
'SigPSF': '5.0',
'PSFres': '1',
'PSFPhot': '1',
'FitSky': '1',
'Force1': '0',
'Align': '2',
'Rotate': '1',
#'ACSuseCTE': '1',  #done in the runner function
'FlagMask': '4',
'ACSpsfType': '0',
'DiagPlotType': 'PS',
'VerboseData': '1'
}

default_calcsky_args = ['15', '35', '-128', '2.25', '2.00']

def copy_files(asns, dest_dir, allowexistingdata=False, cte=True, incldrz=True,
               extrafns=[]):
    from warnings import warn

    #first check that the destination dir exists and copy over the data
    if os.path.exists(dest_dir):
        if allowexistingdata == 'clobber':
            print("Deleting directory", dest_dir)
            shutil.rmtree(dest_dir)
        elif allowexistingdata:
            warn('Destination directory "{0}" already exists.'.format(dest_dir))
        else:
            raise IOError('Destination directory "{0}" already exists!'.format(dest_dir))
    if not os.path.exists(dest_dir):
        os.mkdir(dest_dir)

    toprocess_fns = []
    for asn in asns:
        tocpy = []
        if incldrz:
            tocpy.append(asn.drc if cte else asn.drz)
        tocpy.extend(asn.flcs if cte else asn.flts)
        for fn in tocpy:
            targfn = os.path.join(dest_dir, os.path.split(fn)[-1])
            toprocess_fns.append(os.path.split(targfn)[-1].strip())
            if os.path.exists(targfn):
                print(targfn, 'already exists, not copying')
            else:
                print('Copying', fn, '->', targfn)
                shutil.copy(fn, targfn)

    for fn in extrafns:
        targfn = os.path.join(dest_dir, os.path.split(fn)[-1])
        if os.path.exists(targfn):
            print(targfn, 'already exists, not copying')
        else:
            print('Copying', fn, '->', targfn)
            shutil.copy(fn, targfn)

    return toprocess_fns


def do_prep(working_dir, toprocess_fns, dolphotpath=None, calcskyoverride=None):
    outputs = {}
    #replace the 'acsmask' command with various acs commands later
    dptool_runner = DolphotRunner('acsmask', workingdir=working_dir, execpathordirs=dolphotpath,
                                             paramfile=None, logfile=None)

    #first we check if the drz/drc files are missing the FILETYPE header.  If so
    #that means we ran acsmask already, so we skip them
    acsmask_fns = []
    for fn in toprocess_fns:
        if '_drc.fits' in fn or '_drz.fits' in fn:
            fullfn = os.path.join(working_dir, fn)
            if 'DOL_ACS' in fits.getheader(fullfn, 1):
                print('File', fn, 'has already been acsmasked.')
                continue  # don't run acsmask on it
        acsmask_fns.append(fn)

    if len(acsmask_fns) > 0:
        print('\n...Running acsmask...\n')
        outputs[dptool_runner.cmd] = dptool_runner(*acsmask_fns)

    print('\n...Running splitgroups...\n')
    dptool_runner.cmd = 'splitgroups'
    outputs[dptool_runner.cmd] = splout = dptool_runner(*toprocess_fns)

    chip_fns = []
    for line in splout.split('\n'):
        if line.startswith('Writing FITS file '):
            fn = line.split(':')[0].replace('Writing FITS file ', '').strip()
            chip_fns.append(fn)

    print('\n...Running calcsky...\n')
    if calcskyoverride:
        args = list(calcskyoverride)
    else:
        args = list(default_calcsky_args)
    print('calcsky args:', args)
    args.insert(0, '')  # replaced in the for loop

    dptool_runner.cmd = 'calcsky'
    outputs[dptool_runner.cmd] = calsky_outputs = []
    for fn in chip_fns:
        args[0] = fn.split('.fits')[0]
        calsky_outputs.append(dptool_runner(*args))

    return outputs


def do_dolphot(working_dir, reffn, imgfns, outbase, dolphotpath=None, paramoverrides={}):
    params = dict(default_dolphot_params)
    params.update(paramoverrides)

    #first check cte status
    flts = flcs = drc = drz = False
    for fn in imgfns:
        if '_flt.fits' in fn:
            flts = True
        if '_flc.fits' in fn:
            flcs = True
    if '_drc.fits' in reffn:
        drc = True
    if '_drz.fits' in reffn:
        drz = True

    cte = None
    if flts and flcs:
        warn('Mixed flts and flcs - unclear if CTE correction should be applied!')
    elif flts and drc:
        warn('Mixed flts and CTEd drizzle - unclear if CTE correction should be applied!')
    elif flcs and drz:
        warn('Mixed flcs and CTEd not drizzle - unclear if CTE correction should be applied!')
    elif flcs:
        cte = True
    elif flts:
        cte = False

    if cte:
        if params.get('ACSuseCTE', False):
            warn("You're using CTE corrected files but also asked for the correction... mistake?")
        else:
            print('Will not do CTE correction in dolphot')
            params['ACSuseCTE'] = '0'

    else:
        if not params.get('ACSuseCTE', False):
            warn("You're using CTE uncorrected files but also didn't ask for the correction... mistake?")
        else:
            print('Will do CTE correction in dolphot')
            params['ACSuseCTE'] = '1'

    params['Nimg'] = len(imgfns)
    if reffn is not None:
        params['img0_file'] = reffn.split('.fits')[0]
    for i, fn in enumerate(imgfns):
        params['img' + str(i+1) + '_file'] = fn.split('.fits')[0]

    dolphot_runner = DolphotRunner('dolphot', workingdir=working_dir,
                                   execpathordirs=dolphotpath, params=params)

    return dolphot_runner(outbase)


def do_all_dolphot(working_dir, asns, outbase, chipnum='all', cte=True,
                   allowexistingdata=True, dolphotpath=None,
                   calcskyoverride=None, dolphotparamoverrides={},
                   reffn='guess'):

    if dolphotparamoverrides.get('FakeStars', False):
        dolphotparamoverrides.setdefault('DiagPlotType', '')

    if reffn == 'guess':
        toprocess_fns = copy_files(asns, working_dir, cte=cte, incldrz=True,
                                   allowexistingdata=allowexistingdata)
    else:
        toprocess_fns = copy_files(asns, working_dir, cte=cte, incldrz=False,
                                   allowexistingdata=allowexistingdata,
                                   extrafns=[reffn])
        toprocess_fns.append(os.path.split(reffn)[-1])

    output = do_prep(working_dir, toprocess_fns, dolphotpath=dolphotpath,
                     calcskyoverride=calcskyoverride)
    output['prep_fns'] = toprocess_fns

    imgfns = []
    drfns = []
    for line in output['splitgroups'].split('\n'):
        if line.startswith('Writing FITS file '):
            fn = line[18:].split(':')[0].strip()
            if (cte and '_flc' in fn) or (not cte and '_flt' in fn):
                imgfns.append(fn)
            if (cte and '_drc' in fn) or (not cte and '_drt' in fn):
                drfns.append(fn)
    if reffn == 'guess':
        if len(drfns) == 0:
            warn('Could not find any ref files - not using ref')
            reffn = None
        else:
            reffn = drfns[0]
            if len(drfns) > 1:
                warn('Found multiple possible drizzle refs.  Using first one.')
    else:
        # should have been copied and split, so we use the working_dir version
        # that has already gotten splitgroup'd
        reffn = os.path.split(reffn)[-1].replace('.fits', '.chip1.fits')

    if chipnum != 'all':
        chipstr = '.chip' + str(chipnum)
        imgfns = [fn for fn in imgfns if chipstr in fn]

    print('Running dolphot on', imgfns, 'with ref', reffn)

    output['dolphot_reffn'] = reffn
    output['dolphot_imgfns'] = imgfns
    output['dolphot'] = do_dolphot(working_dir, reffn, imgfns, outbase,
                                   dolphotpath=dolphotpath,
                                   paramoverrides=dolphotparamoverrides)

    return output


def do_drizzle(asns, outname, working_dir='sameasoutname', cte=True, **options):
    """
    Note that this requires Ureka to be activated

    A good baseline of options for 2-image ACS/WFC is:
    ``final_pixfrac=0.6,final_scale=.03, final_wcs=True``

    `working_dir` can have the special value "sameasoutname", or anything ending
    with an "_" will have `outname` added at the end
    """
    from stsci.tools import teal
    from drizzlepac import astrodrizzle

    if working_dir == 'sameasoutname':
        working_dir = outname
    if working_dir.endswith('_'):
        working_dir = working_dir + outname

    #toprocess_fns are *relative* to the working_dir, not absolute paths
    toprocess_fns = copy_files(asns, working_dir, cte=cte, incldrz=False,
                               allowexistingdata=True)

    #reset to defaults
    teal.unlearn('astrodrizzle')

    options.setdefault('preserve', False)
    options.setdefault('restore', False)

    if len(set([asn.filter for asn in asns])) > 1:
        if 'final_wht_type' not in options:
            print("You are using multiple filters at the same time but didn't "
                  "set a wheight type.  You probably don't want EXP, so we'll "
                  "try IVM")
            options['final_wht_type'] = 'IVM'

    olddir = os.path.abspath(os.curdir)
    try:
        os.chdir(working_dir)
        astrodrizzle.AstroDrizzle(toprocess_fns, output=outname, **options)
    finally:
        os.chdir(olddir)

    return working_dir




