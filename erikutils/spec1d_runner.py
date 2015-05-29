from __future__ import division, print_function
"""
Scripts to automate running the DEIMOS spec1d pipeline (in IDL)
"""

import os
import sys
import time


invoke_spec1d_templ = "et_spec1d,'{maskname}'\nexit\n"
invoke_mcerrs_templ = "et_mcerrs,'{maskname}',{nmc}\nexit\n"


def _do_invoke(maskname, datadir, idlsrc, logsuffix):
    import subprocess
    from warnings import warn

    if datadir is None:
        if 'DEIMOS_DATA' not in os.environ:
            raise ValueError('datadir not given and DEIMOS_DATA not set')
        path = os.path.join(os.environ['DEIMOS_DATA'], maskname)
    else:
        path = os.path.join(datadir, maskname)

    planfn = os.path.join(path, maskname + '.plan')
    if not os.path.exists(planfn):
        warn('Plan file "{0} not found - is something wrong? continuing anyway'.format(planfn))

    logfn = os.path.abspath(os.path.join(path, maskname + '_{0}.log'.format(logsuffix)))
    logf = open(logfn, 'w')

    newenv = os.environ.copy()
    if datadir is not None:
        if not os.path.isdir(datadir):
            raise ValueError('datadir "{0}" is not a directory'.format(datadir))
        newenv['DEIMOS_DATA'] = datadir + ('' if datadir.endswith(os.sep) else os.sep)

    proc = subprocess.Popen('idl', cwd=path, stdin=subprocess.PIPE, stdout=logf,
                            stderr=subprocess.STDOUT, env=newenv)
    proc.stdin.write(idlsrc)
    proc.maskname = maskname
    return proc


def invoke_spec1d(maskname, datadir=None):
    """
    Runs spec1d on the given mask dir, which should be given as a subdir of
    DEIMOS_DATA or `datadir` can replace DEIMOS_DATA

    Note that you have to manually close the returned proc.stdout!
    """
    idlsrc = invoke_spec1d_templ.format(**locals())
    proc = _do_invoke(maskname, datadir, idlsrc, '1d')
    return proc


def invoke_mcerrs(maskname, nmc=500, datadir=None):
    """
    Runs mcerrs on the given mask dir, which should be given as a subdir of
    DEIMOS_DATA or `datadir` can replace DEIMOS_DATA

    Note that you have to manually close the returned proc.stdout!
    """
    idlsrc = invoke_mcerrs_templ.format(**locals())
    proc = _do_invoke(maskname, datadir, idlsrc, 'mcerr')
    return proc


def try_finish_spec1d(proc):
    if proc.poll() is None:
        return False
    else:
        if proc.returncode != 0:
            print('The process for plan file "{0}" returned {1}... '
                  'possible problem? Check logs.'.format(proc.maskname, proc.returncode))
        if proc.stdout is not None and not proc.stdout.closed:
            proc.stdout.close()
        if proc.stderr is not None and not proc.stderr.closed:
            proc.stderr.close()
        return True


def find_finished_masks(msknames):
    finished_masks = []
    for nm in msknames:
        if os.path.isdir(nm):
            if os.path.isfile(os.path.join(nm, 'doneprocessing.txt')):
                print("Found doneprocessing for", nm)
                finished_masks.append(nm)
            else:
                print("doneprocessing was not found for", nm, 'skipping!')

    return finished_masks


def scatter_spec1ds(dirs, maxtorun=2, waittime=1, verbose=True, nmcerrs=None):
    """
    `dirs` is list of directories with 2d reduced DEIMOS data
    `maxtorun` is the number of simultaneous processes to run
    `waittime` is the time in sec to wait between polling

    `nmcerrs` determines which part gets run: if None, it means run spec1d,
    otherwise `mcerrs`, with `nmc` given by the value of nmcerrs
    """
    procsdone = []
    procsrunning = []

    toinvoke = []
    dirs = [os.path.abspath(d) for d in dirs]
    cpfx, cnm = os.path.split(os.path.commonprefix(dirs))
    cpfx = cpfx + os.sep
    for d in dirs:
        name = d[len(cpfx):]
        if os.path.split(name)[0] != '':
            raise ValueError('Maskname "{0}" has directory components'.format(name))
        toinvoke.append((name, cpfx))

    sleepsdone = 0
    while len(toinvoke) > 0 or len(procsrunning) > 0:
        # first check if any are running that have finished
        for i, p in reversed(list(enumerate(procsrunning))):
            if try_finish_spec1d(p):  # True -> proc done
                if verbose:
                    print('\nFinished running for', p.maskname)
                del procsrunning[i]
                procsdone.append(p)
                sleepsdone = 0

        # now try to invoke any that remain to be invoked
        rem_from_toinvoke = []
        for i, (name, dd) in enumerate(toinvoke):
            if len(procsrunning) < maxtorun:
                if nmcerrs is None:
                    if verbose:
                        print('\nInvoking spec1d for', name, dd)
                    procsrunning.append(invoke_spec1d(name, dd))
                else:
                    if verbose:
                        print('\nInvoking nmcerrs for', name, dd)
                    procsrunning.append(invoke_mcerrs(name, nmcerrs, dd))
                rem_from_toinvoke.append(i)
                sleepsdone = 0
        for i in reversed(sorted(rem_from_toinvoke)):
            del toinvoke[i]

        if verbose:
            sys.stdout.write('Sleeping for {0} sec\r'.format(waittime*sleepsdone))
            sys.stdout.flush()

        time.sleep(waittime)
        sleepsdone += 1

    return dict([(p.maskname, p) for p in procsdone])
