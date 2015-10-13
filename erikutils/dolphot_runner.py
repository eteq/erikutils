from __future__ import division

import os
import sys
import pty
import select
import subprocess


__all__ = ['DolphotRunner']


def find_dolphot(dirstocheck=None):
    if dirstocheck is None:
        return subprocess.check_output(['which', 'dolphot'])
    else:
        for dr in dirstocheck:
            path = os.path.join(dr, 'dolphot')
            if os.path.isfile(os.path.join(dr, 'dolphot')):
                return path
        raise ValueError('Could not find dolphot in {0}'.format(dirstocheck))


class DolphotRunner(object):
    """
    This class is a convenience tool for running the Dolphot photometry code.
    Use it something like this::

        runner = DolphotRunner()
        params = {'Align':1, 'UseWCS':1, ... }
        runner('my_awesome_hst_data1.flc', 'my_awesome_hst_data2.flc', params=params)
    """
    def __init__(self, cmd='dolphot', logfile='auto', paramfile='auto',
                       execpathordirs=None, workingdir='.', params={}):
        """
        if `logfile` is "auto", it means use the first argument
        """
        if isinstance(execpathordirs, basestring):
            dolphot_path = execpathordirs
            if not os.path.isfile(dolphot_path):
                raise ValueError('could not find dolphot path "{0}"'.format(dolphot_path))
        else:
            dolphot_path = find_dolphot(execpathordirs)
        self.dolphot_bin_dir = os.path.abspath(os.path.split(dolphot_path)[0])

        self.workingdir = workingdir
        self.cmd = cmd
        self.logfile = logfile
        self.paramfile = paramfile
        self.params = params

    @property
    def workingdir(self):
        return self._workingdir
    @workingdir.setter
    def workingdir(self, value):
        self._workingdir = os.path.abspath(value)

    def __call__(self, *args):

        exec_path = os.path.join(self.dolphot_bin_dir, self.cmd)
        if self.logfile == 'auto':
            if params.get('FakeStars', False):
                logfile = args[0] + '_' + self.cmd + '_fakestars.log'
            else:
                logfile = args[0] + '_' + self.cmd + '.log'
        else:
            logfile = self.logfile

        if self.paramfile == 'auto':
            paramfile = args[0] + '_' + self.cmd + '.param'
        else:
            paramfile = self.paramfile

        args = list(args)
        args.insert(0, exec_path)

        if paramfile:
            if self.params:
                with open(os.path.join(self.workingdir, paramfile), 'w') as f:
                    for k, v in self.params.items():
                        f.write('{0} = {1}\n'.format(k, v))
                args.append('-p' + paramfile)
        else:
            for k, v in self.params.items():
                args.append('{0}={1}'.format(k, v))

        output = []
        # now actually do the work with the process itself.  We use a pseudo-tty
        # mostly just to allow line-buffering, which doesn't happen with pipes.
        master_fd, slave_fd = pty.openpty()
        try:

            p = subprocess.Popen(args, cwd=self.workingdir, bufsize=1,
                                       stdout=slave_fd, stderr=subprocess.STDOUT,
                                       close_fds=True)

            pty_timeout = .05 # 50 ms


            donereading = False
            while not donereading:
                ready = select.select([master_fd], [], [], pty_timeout)[0]
                if ready:
                    data = os.read(master_fd, 512)
                    if data:
                        output.append(data)
                        sys.stdout.write(data)
                        sys.stdout.flush()
                    else:
                        donereading = True
                elif p.poll() is not None:
                    donereading = True
        finally:
            os.close(slave_fd) # can't do it sooner: it leads to errno.EIO error
            os.close(master_fd)

        # normalize line endings because it seems os.read doesn't?
        self.last_out = ''.join(output).replace('\r\n','\n').replace('\r','\n')

        if logfile is not None:
            with open(os.path.join(self.workingdir, logfile), 'w') as f:
                f.write(self.last_out)

        if p.returncode != 0:
            raise ValueError('command "{0}" returned {1} instead of '
                             '0'.format(self.cmd, p.returncode), self.last_out)
        return self.last_out
