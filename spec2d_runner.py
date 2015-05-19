from __future__ import division, print_function
"""
Scripts to automate running the DEIMOS spec1d pipeline (in IDL)
"""

import os
import sys
import time


invoke_spec2d_templ = "et_domask,'{planfn}'\nexit\n"


def invoke_spec2d(path, maskname):
    """
    Runs spec2d in the given path, assuming there's a {maskname}.plan file

    Note that you have to manually close the returned proc.stdout!
    """
    import subprocess

    planfn = os.path.abspath(os.path.join(path, maskname + '.plan'))
    logfn = os.path.abspath(os.path.join(path, maskname + '.log'))

    if not os.path.isdir(path):
        raise IOError('"{0}" is not a directory!'.format(path))

    if not os.path.isfile(planfn):
        raise IOError('Plan file "{0}" does not exist!'.format(planfn))

    logf = open(logfn, 'w')

    proc = subprocess.Popen('idl', cwd=path, stdin=subprocess.PIPE,
                            stdout=logf, stderr=subprocess.STDOUT)
    proc.stdin.write(invoke_spec2d_templ.format(**locals()))
    # proc = subprocess.Popen('ls', cwd=path, stdin=None,
    #                         stdout=logf, stderr=subprocess.STDOUT)
    proc.maskname = maskname

    return proc


def try_finish_spec2d(proc):
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


def find_unfinished_planfiles(msknames):
    planfiles = []
    for nm in msknames:
        if os.path.isfile(nm):
            planfiles.append(nm)
        elif os.path.isdir(nm):
            path, name = os.path.split(nm)
            if name == '':
                nm = path
                path, name = os.path.split(nm)
            planfiles.append(os.path.join(path, name, name + '.plan'))

    for i, pf in reversed(list(enumerate(planfiles))):
        path, name = os.path.split(pf)
        if os.path.isfile(os.path.join(path, 'doneprocessing.txt')):
            print("doneprocessing was found for", name, 'skipping!')
            del planfiles[i]

    return planfiles


def scatter_spec2ds(planpaths, maxtorun=2, waittime=1, verbose=True):
    """
    `planpaths` is list of planfiles
    `maxtorun` is the number of simultaneous processes to run
    `waittime` is the time in sec to wait between polling
    """
    procsdone = []
    procsrunning = []

    toinvoke = []
    for plp in planpaths:
        if plp.endswith('.plan'):
            plp = plp[:-5]

        path, name = os.path.split(plp)
        toinvoke.append((path, name))

    sleepsdone = 0
    while len(toinvoke) > 0 or len(procsrunning) > 0:
        #first check if any are running that have finished
        for i, p in reversed(list(enumerate(procsrunning))):
            if try_finish_spec2d(p):  # True -> proc done
                if verbose:
                    print('\nFinished spec2d for', p.maskname)
                del procsrunning[i]
                procsdone.append(p)
                sleepsdone = 0

        #now try to invoke any that remain to be invoked
        rem_from_toinvoke = []
        for i, (path, name) in enumerate(toinvoke):
            if len(procsrunning) < maxtorun:
                if verbose:
                    print('\nInvoking spec2d for', path, name)
                procsrunning.append(invoke_spec2d(path, name))
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
