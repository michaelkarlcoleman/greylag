#!/usr/bin/env python

'''Run a greylag job split across multiple local processes.  This is not for
use on a cluster, but useful if you have just one multi-CPU machine, or for
debugging.
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
# multiple CPUs under Windows, as X!Tandem can.

# This is not intended to be used on cluster nodes.  In particular, like
# X!Tandem, it just divides the spectra into N parts and then processes them
# separately, making no further attempt at load balancing.


import os
import os.path
from socket import gethostname
import subprocess
import sys


GREYLAG_PROGRAM = 'greylag'


def warn(s):
    print >> sys.stderr, 'warning:', s
def error(s):
    sys.exit('error: ' + s)

def usage():
    print >> sys.stderr, ('Usage: %s <processes> <greylag-options-and-args>...'
                          '\n\n'
                          '%s\n(see "greylag.py --help" for more information)'
                          % (os.path.basename(sys.argv[0]), __doc__))
    sys.exit()


def run_parts(processes, partprefix, args):
    ret = subprocess.call([GREYLAG_PROGRAM, '--part-split=%s' % processes,
                           '--part-prefix=%s' % partprefix] + args)
    if ret:
        error("part split failed")

    # list of currently running subprocess.Popen objects
    subprocs = []
    try:
        for n in range(processes):
            p = subprocess.Popen([GREYLAG_PROGRAM,
                                  '--part=%sof%s' % (n+1, processes),
                                  '--part-prefix=%s' % partprefix] + args)
            subprocs.append(p)
        while subprocs:
            p = subprocs.pop()
            if p.wait():
                raise EnvironmentError("error status")
    except EnvironmentError, e:
        error("part process failed [%s]" % e)
    finally:
        # upon error, try to kill any remaining processes
        try:
            import signal
            for p in subprocs:
                os.kill(p.pid, signal.SIGINT)
        except Exception:
            pass

    ret = subprocess.call([GREYLAG_PROGRAM, '--part-merge=%s' % processes,
                           '--part-prefix=%s' % partprefix] + args)
    if ret:
        error("part merge failed")


def main():
    args = sys.argv[1:]
    if len(args) < 2 or '-h' in args or '--help' in args:
        usage()
    try:
        processes = int(args[0])
        if processes < 1:
            raise ValueError
    except ValueError:
        error("<processes> must be a positive integer")
    # This limit can be raised by setting environment variable.
    GREYLAGPROCESSES = os.environ.get('GREYLAGPROCESSES', 8)
    if processes > GREYLAGPROCESSES:
        error("current limit is %s processes (see source code for more)"
              % GREYLAGPROCESSES)
    if any(x for x in args if x.startswith('--part')):
        error("no --part options may be specified")

    # try to generate a unique prefix, to avoid (unlikely) collision
    partprefix = '#greylag-part-%s-%s' % (gethostname(), os.getpid())

    try:
        run_parts(processes, partprefix, args[1:])
    finally:
        # try to clean up temporary part files
        for fn in os.listdir('.'):
            if fn.startswith(partprefix):
                try:
                    os.remove(fn)
                except IOError, e:
                    warn("could not remove '%s' [%s]" % (fn, e))


if __name__ == '__main__':
    main()
