#!/usr/bin/env python

"""
Convert between FASTA-formatted input and one-sequence-per-line output so that
the sequences can be easily manipulated with UNIX text tools (e.g., grep,
head, wc, split, sort, etc.).

In order for '--inverse' to work correctly, the same flags must be supplied as
were supplied during the forward conversion (the script does not try to
guess).  With '--defline=after', the conversion should be perfectly
invertible, modulo whitespace and wrapping.  For '--defline=omit', an
artificial defline will be constructed based on the filename and line number.

"""

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


import fileinput
import optparse
import re
import sys

import greylag


# no backtrace on SIGPIPE
try:
    import signal
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
except Exception:
    pass


def error(s, *args):
    "fatal error"
    # if we're unit testing, just throw an exception
    if __name__ != "__main__":
        raise Exception((s + " (fatal error)") % args)
    print >> sys.stderr, ("error: " + s) % args
    sys.exit(1)

# errors are fatal
greylag.chase_error = error


def warn(message):
    print >> sys.stderr, "warning: %s [at %s:%s]" \
          % (message, fileinput.filename(), fileinput.filelineno())


def write_flattened_locus(options, defline, sequence):
    if options.defline == 'after':
        print '%s%s>%s' % (sequence, options.delimiter, defline)
    elif options.defline == 'before':
        if options.delimiter in defline:
            warn("delimiter present in defline")
        print '>%s%s%s' % (defline, options.delimiter, sequence)
    else:
        print sequence


def _main():
    parser = optparse.OptionParser(usage="usage: %prog [options] [<file>...]",
                                   description=__doc__)
    parser.add_option("-d", "--delimiter", dest="delimiter", default='\t',
                      help="delimiter between defline and sequence"
                      " [default TAB]", metavar="STRING")
    parser.add_option("-D", "--defline", dest="defline",
                      choices=('before', 'after', 'omit'), default="after",
                      help="position of defline with respect to sequence, one"
                      " of 'before', 'after' [default], or 'omit'",
                      metavar="POSITION")
    parser.add_option("-i", "--inverse", dest="inverse", action="store_true",
                      help="do the inverse transformation (flat to FASTA)")
    DEFAULT_WRAP = 80
    parser.add_option("-w", "--wrap", dest="wrap", type="int",
                      default=DEFAULT_WRAP,
                      help="for --inverse, wrap sequence to specified width"
                      " [default %s, 0 means don't wrap at all]" % DEFAULT_WRAP,
                      metavar="COLUMNS")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
                      help="be verbose")
    parser.add_option("--copyright", action="store_true", dest="copyright",
                      help="print copyright and exit")
    options, args = parser.parse_args()

    if options.wrap < 0:
        parser.print_help()
        sys.exit(1)

    if not options.inverse:
        if not args:
            files = [ sys.stdin ]
        else:
            files = [ open(fn) for fn in args ]

        for f in files:
            for locusname, defline, sequence in greylag.read_fasta_file(f):
                write_flattened_locus(options, defline, sequence)
    else:
        for line in fileinput.input(args):
            if options.defline != 'omit':
                parts = line.split(options.delimiter, 1)
                if len(parts) < 2:
                    error("input line lacks delimiter")
                if options.defline == 'before':
                    defline, sequence = parts
                else:
                    sequence, defline = parts
            else:
                sequence = line
                defline = "%s:%s" % (fileinput.filename(),
                                     fileinput.filelineno())
            sequence = sequence.strip()
            print defline.strip()
            if options.wrap:
                for start in range(0, len(sequence), options.wrap):
                    print sequence[start:start+options.wrap]
            else:
                print sequence


if __name__ == '__main__':
    _main()
