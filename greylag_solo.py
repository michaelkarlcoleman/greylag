#!/usr/bin/env greylag-python

'''Run a greylag job split across multiple local processes.  This is not for
use on a cluster, but useful if you have just one multi-CPU machine, or for
debugging.

If you specify more processes than the total number of spectra in the input
files, you will get harmless "no input spectra" warnings, which may be
ignored.
'''

__copyright__ = '''
    greylag, Copyright (C) 2006-2007, Stowers Institute for Medical Research

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
'''


# This might better be implemented as a shell script under Unix.  It's done in
# Python here as a demonstration, and so that greylag can be easily run on
# multiple CPUs under Windows.

# This is not intended to be used on cluster nodes.  In particular, it just
# naively divides the spectra into N parts and processes them separately,
# making no further attempt at load balancing.


import os
import os.path
from socket import gethostname
import subprocess
import sys


def warn(s):
    print >> sys.stderr, 'warning:', s
def error(s):
    sys.exit('error: ' + s)

def usage():
    print >> sys.stderr, ('Usage: %s <processes>'
                          ' <greylag-grind-options-and-args>...'
                          '\n\n'
                          '%s\n(see "greylag-grind --help" for more'
                          ' information)'
                          % (os.path.basename(sys.argv[0]), __doc__))
    sys.exit()


def slices(N):
    """Generate N equal slices.

    >>> list(slices(4))
    [(0.0, 0.25), (0.25, 0.5), (0.5, 0.75), (0.75, 1.0)]

    """

    s = 1.0 / N
    for i in range(N):
        w0 = i * s
        w1 = (i+1) * s
        if i == N-1:
            w1 = 1.0                    # forestall rounding issues
        yield (round(w0,3), round(w1,3))


def run_parts_and_merge(processes, job_id, args):
    merge_fn = 'greylag-merge-%s.glw' % job_id
    work_fns = [ 'grind_%s_%s-%s.glw' % (job_id, w0, w1)
                 for w0, w1 in slices(processes) ]

    try:
        # list of currently running subprocess.Popen objects
        subprocs = []
        try:
            for w0, w1 in slices(processes):
                p = subprocess.Popen(['greylag-grind', '--job-id='+job_id,
                                      '-w', str(w0), str(w1)] + args)
                subprocs.append(p)
            while subprocs:
                p = subprocs.pop()
                if p.wait():
                    raise EnvironmentError("error status")
        except EnvironmentError, e:
            error("greylag-grind process failed [%s]" % e)
        finally:
            # upon error, try to kill any remaining processes
            try:
                import signal
                for p in subprocs:
                    os.kill(p.pid, signal.SIGINT)
            except Exception:
                pass

        ret = subprocess.call(['greylag-merge'] + work_fns + [merge_fn])
        if ret:
            error("greylag-merge failed")

        ret = subprocess.call(['greylag-sqt', merge_fn])
        if ret:
            error("greylag-sqt failed")
    finally:
        try:
            os.remove(merge_fn)
            for fn in work_fns:
                os.remove(fn)
        except OSError:
            pass


def main(args=sys.argv[1:]):
    if len(args) < 3 or '-h' in args or '--help' in args:
        usage()
    try:
        processes = int(args[0])
        if processes < 1:
            raise ValueError
    except ValueError:
        error("<processes> must be a positive integer")

    # Keep naive users out of trouble.  This limit can be raised by setting
    # environment variable.
    GREYLAGPROCESSES = os.environ.get('GREYLAGPROCESSES', 8)
    if processes > GREYLAGPROCESSES:
        error("current limit is %s processes (see source code for more)"
              % GREYLAGPROCESSES)

    for a in args:
        if a.startswith('--job-id'):
            error('--job-id may not be specified')
        if a.startswith(('-w', '--work-slice')):
            error('--work-slice may not be specified')

    # try to generate a unique prefix, to avoid (unlikely) collision
    job_id = 'solo-%s-%s' % (gethostname(), os.getpid())

    run_parts_and_merge(processes, job_id, args[1:])


if __name__ == '__main__':
    main()
