#!/usr/bin/env python

'''Create a trivial index giving the starting point of each spectrum, as a
   byte offset from the previous starting point.  Also checks that spectra
   names are unique and that spectra are ordered by name.
'''

from __future__ import with_statement

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

__version__ = "0.0"


import optparse
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
        prevname = ''
        offset = 0
        with open(fn + '.idx', 'w') as idx:
            for m in re.finditer('^:.*$', contents, re.MULTILINE):
                specname = m.group()
                if specname in specnames:
                    error("duplicate spectrum names not allowed [%s]"
                          % specname)
                specnames.add(specname)
                if not prevname < specname:
                    error("spectra must be ordered by name [%s]" % specname)
                prevname = specname
                print >> idx, m.start() - offset
                offset = m.start()


if __name__ == '__main__':
    main()
