#!/usr/bin/env python

'''Merge a set of greylag search result files, producing a search result file
that summarizes the best matches found.

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


import contextlib
import cPickle as pickle
import optparse
import os.path
from pprint import pprint
import sys

from greylag import VERSION


# gc kills our performance, so disable it.  gc only matters for cycles, which
# we (hope we) don't create.  See the gc module docs.
import gc; gc.disable()


def warn(s):
    print >> sys.stderr, 'warning:', s
def error(s):
    sys.exit('error: ' + s)


def check_consistency(r0, r1):
    """Warn on any inconsistencies between the two results that are about to
    be merged."""

    if r0['version'] != r1['version']:
        warn("merging runs searched with different greylag versions!"
             " (rerunning with same version strongly recommended)")
    k0 = set(r0.keys())
    k1 = set(r1.keys())
    if k0 != k1:
        warn("mismatched keys: %s" % (k0 ^ k1))

    # FIX: this restriction goes away later
    if r0['spectrum files'] != r1['spectrum files']:
        error('merging across differing sets of ms2 files not yet implemented'
              ' (%s vs %s)' % (r0['spectrum files'], r1['spectrum files']))

    varying_keys = set(['matches', 'total comparisons', 'argv'])
    for k in k0 - varying_keys:
        if r0[k] != r1.get(k):
            warn("search context differs: %s (%s vs %s)" % (k, r0[k], r1[k]))


def score_equal(s1, s2):
    return abs(s1-s2) < 1e-6


def merge_match(m0, m1, keep):
    """Merge a particular spectrum match from m1 into m0."""
    m0[0]['comparisons'] += m1[0]['comparisons']

    # NB: This merge code, and score_equal above, must do the same merge as
    # the merge code in cgreylag.cpp:evaluate_peptide!

    assert keep >= 1
    matches = m0[1] + m1[1]
    matches.sort(key=lambda x: x['score'], reverse=True)

    merged_matches = [matches.pop()]
    while len(merged_matches) < keep and matches:
        c0, c1 = merged_matches[-1], matches.pop()
        if (score_equal(c0['score'], c1['score'])
            and c0['peptide_sequence'] == c1['peptide_sequence']
            and c0.get('mass_trace') == c1.get('mass_trace')):
            continue
        merged_matches.append(c1)


def merge_matches(m0, m1, keep):
    """Merge the match information in m1 into m0."""
    for m1k, m1v in m1.iteritems():
        if m1k not in m0:
            m0[m1k] = m1v
        else:
            merge_match(m0[m1k], m1v, keep)


def main(args=sys.argv[1:]):
    parser = optparse.OptionParser(usage=
                                   "usage: %prog [options] [<result-file>...]"
                                   " <output-result-file>",
                                   description=__doc__, version=VERSION)
    pa = parser.add_option
    pa("-f", "--files-on-stdin", action="store_true", dest="files_on_stdin",
       help="read filenames to be merged from stdin, one per line, instead of"
       " from the command-line.")
    pa("-v", "--verbose", action="store_true", dest="verbose",
       help="be verbose")
    pa("--copyright", action="store_true", dest="copyright",
       help="print copyright and exit")
    (options, args) = parser.parse_args(args=args)

    if options.copyright:
        print __copyright__
        sys.exit(0)

    if options.files_on_stdin:
        if len(args) != 1:
            parser.print_help()
            sys.exit(1)
        filenames = [ l.strip() for l in sys.stdin ]
        if len(filenames) < 1:
            error("--files-on-stdin given, but stdin was empty")
        result_fn_0 = filenames[0]
        result_fn_1_N = filenames[1:]
    else:
        if len(args) < 2:
            parser.print_help()
            sys.exit(1)
        result_fn_0 = args[0]
        result_fn_1_N = args[1:-1]

    output_fn = args[-1]

    # spectrum name -> match list/spectrum info
    matches = {}
    with contextlib.closing(open(result_fn_0, 'rb')) as r_file:
        r0 = pickle.load(r_file)
    total_comparisons = r0['total comparisons']
    keep = r0['parameters']['best_result_count']
    matches = r0['matches']
    if options.verbose:
        print >> sys.stderr, "loaded", result_fn_0

    for additional_result_fn in result_fn_1_N:
        with contextlib.closing(open(additional_result_fn, 'rb')) as r1_file:
            r1 = pickle.load(r1_file)
        check_consistency(r0, r1)
        total_comparisons += r1['total comparisons']
        keep = max(keep, r1['parameters']['best_result_count'])
        merge_matches(matches, r1['matches'], keep)
        if options.verbose:
            print >> sys.stderr, "merged", additional_result_fn
        del r1                          # free some memory

    r0['matches'] = matches
    r0['total comparisons'] = total_comparisons
    with contextlib.closing(open(output_fn, 'wb')) as output_file:
        pk = pickle.Pickler(output_file, pickle.HIGHEST_PROTOCOL)
        pk.fast = 1                     # stipulate no circular references
        pk.dump(r0)
    if options.verbose:
        print >> sys.stderr, "dumped", output_fn


if __name__ == '__main__':
    main()
