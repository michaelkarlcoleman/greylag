#!/usr/bin/env python

'''Generate a set of sqt files from the specified greylag search result
file.  The sqt files will be created in the current directory, with names
corresponding to the searched ms2 files.
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


import contextlib
import cPickle
import gzip
import math
import optparse
import os.path
from pprint import pprint
import sys


def warn(s):
    print >> sys.stderr, 'warning:', s
def error(s):
    sys.exit('error: ' + s)


def print_header(f, r):
    print >> f, "H\tSQTGenerator\tgreylag"
    print >> f, "H\tSQTGeneratorVersion\t%s" % __version__
    for db in r['databases']:
        print >> f, "H\tDatabase\t%s" % db
    for pk, pv in r['parameters'].items():
        print >> f, "H\tParameter\t%s\t%s" % (pk, pv)
    # How should these be specified?
    # What does DTASelect do with these?
    # H       StaticMod       C=160.1388
    # H       DiffMod TNA*=+2.0


def print_spectrum(f, spectrum, sp_best_matches):
    scan_low, scan_high, _rest = spectrum['name'].split('.', 2)
    print >> f, '\t'.join(str(v) for v in
                          ["S", scan_low, scan_high, spectrum['charge'], 0,
                           'honk', spectrum['mass'],
                           round(math.log(spectrum['total_ion_current']), 4),
                           0, spectrum['comparisons']])

    best_scores = sorted(list(set(m['score'] for m in sp_best_matches)))
    # score -> rank
    rank_map = dict(zip(best_scores, range(1,len(best_scores)+1)))

    best_score = sp_best_matches[0]['score']
    prev_score, prev_sequence = None, None
    for match in sp_best_matches:
        score = match['score']
        if score == 0:
            continue                    # null entry
        rank = rank_map[score]
        score_delta = (best_score - score) / best_score

        # This is a temporary way of stripping duplicates.  It only works if
        # they are adjacent (which will only rarely be false?).
        if score != prev_score or match['peptide_sequence'] != prev_sequence:
            # Note: scores, being log probability, are non-positive, but we
            # flip the sign in the SQT output
            print >> f, '\t'.join(str(v) for v in
                                  ["M", rank, rank,
                                   round(match['predicted_parent_mass'], 5),
                                   round(score_delta, 4), round(-score, 4), 0,
                                   # 1 of 2 ions found--keep DTASelect happy
                                   1, 2,
                                   "-.%s.-" % match['peptide_sequence'],
                                   'U'])
        prev_score, prev_sequence = score, match['peptide_sequence']

        print >> f, 'L\t%s\t%s' % (match['sequence_name'],
                                   match['peptide_begin'])

def main(args=sys.argv[1:]):
    parser = optparse.OptionParser(usage=
                                   "usage: %prog [options] <result-file>",
                                   description=__doc__, version=__version__)
    pa = parser.add_option
    pa("-d", "--output-directory", dest="output_directory",
       help="directory where output files are written [default is '.']",
       metavar="DIR")
    pa("-v", "--verbose", action="store_true", dest="verbose",
       help="be verbose")
    pa("--dump", action="store_true", dest="dump",
       help="dump the result file (for debugging)")
    pa("--copyright", action="store_true", dest="copyright",
       help="print copyright and exit")
    (options, args) = parser.parse_args(args=args)

    if options.copyright:
        print __copyright__
        sys.exit(0)

    if len(args) != 1:
        parser.print_help()
        sys.exit(1)

    with contextlib.closing(gzip.open(args[0])) as r_file:
        r = cPickle.load(r_file)

    if options.dump:
        pprint(r)
        return

    spectrum_fns = r['spectrum files']
    assert len(spectrum_fns) == len(set(spectrum_fns)) # check uniqueness

    spectra = r['spectra']
    best_matches = r['matches']['best_matches']
    assert len(spectra) == len(best_matches)

    # order spectra and best_matches by (filename, spectrum name)
    def filename_specname_order(sm):
        return r['spectrum files'][sm[0]['file_id']], sm[0]['name']

    spec_match = zip(spectra, best_matches)
    spec_match.sort(key=filename_specname_order)

    for spectrum_n, spectrum_fn in enumerate(spectrum_fns):
        assert os.path.dirname(spectrum_fn) == ''
        sqt_fn = os.path.splitext(spectrum_fn)[0] + '.sqt'
        with open(sqt_fn, 'w') as sqtf:
            print_header(sqtf, r)
            for spectrum, best_match in spec_match:
                if spectrum['file_id'] == spectrum_n:
                    print_spectrum(sqtf, spectrum, best_match)


if __name__ == '__main__':
    main()
