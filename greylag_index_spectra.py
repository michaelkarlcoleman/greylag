#!/usr/bin/env greylag-python

'''Create a trivial index giving the starting point of each spectrum, as a
   byte offset from the file beginning.  (The index is stored in Python pickle
   format, compressed with gzip.)  Also checks that spectra names are unique
   and that spectra are ordered by name, which other greylag programs assume.
'''

from __future__ import with_statement

__copyright__ = '''
    greylag, a collection of programs for MS/MS protein analysis
    Copyright (C) 2006-2008  Stowers Institute for Medical Research

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Mike Coleman
             Stowers Institute for Medical Research
             1000 East 50th Street
             Kansas City, Missouri  64110
             USA
'''

__version__ = "0.0"


import contextlib
import cPickle
import gzip
import optparse
import os.path
import re
import sys


def error(s):
    sys.exit('error: ' + s)


def main(args=sys.argv[1:]):
    parser = optparse.OptionParser(usage=
                                   "usage: %prog [options] <ms2-file>...",
                                   description=__doc__)
    pa = parser.add_option
    pa("--copyright", action="store_true", dest="copyright",
       help="print copyright and exit")
    pa("--version", action="store_true", dest="version",
       help="print version and exit")
    (options, args) = parser.parse_args(args=args)

    if options.copyright:
        print __copyright__
        sys.exit(0)
    if options.version:
        print __version__
        sys.exit(0)

    if (len(args) < 1
        or any(True for f in args if not f.endswith('.ms2'))):
        parser.print_help()
        sys.exit(1)

    for fn in args:
        with open(fn) as specfile:
            contents = specfile.read()
        specnames = set()
        with contextlib.closing(gzip.open(fn + '.idx', 'w')) as idx:
            ms = [ m for m in re.finditer(r'^:.*$', contents, re.MULTILINE) ]
            specnames = [ m.group() for m in ms ]
            if len(set(specnames)) < len(ms):
                error("duplicate spectrum names not allowed")
            if specnames != sorted(specnames):
                error("spectra must be ordered by name")
            offsets = [ m.start() for m in ms ]
            cPickle.dump({ 'offsets' : offsets,
                           'file size' : os.path.getsize(fn) },
                         idx, cPickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    main()
