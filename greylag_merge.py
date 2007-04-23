#!/usr/bin/env python2.5

'''Merge a set of greylag search result files, producing a search result file
that summarizes the best matches found.

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
import optparse
import os.path
from pprint import pprint
import sys


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

    varying_keys = set(['matches', 'total comparisons', 'argv'])
    for k in k0 - varying_keys:
        if r0[k] != r1.get(k):
            warn("search context differs: %s (%s vs %s)" % (k, r0[k], r1[k]))


def score_equal(s1, s2):
    return abs(s1-s2) < 1e-6


def merge_match(m0, m1):
    """Merge a particular spectrum match from m1 into m0."""
    m0[0]['comparisons'] += m1[0]['comparisons']

    # NB: This merge code, and score_equal above, must functionally match the
    # merge code in cgreylag.cpp:evaluate_peptide!

    keep = max(len(m0), len(m1))
    assert keep >= 1
    matches = m0[1] + m1[1]
    matches.sort(key=lambda x: x['score'], reverse=True)

    merged_matches = [matches.pop()]
    while len(merged_matches) < keep and matches:
        c0, c1 = merged_matches[-1], matches.pop()
        if (score_equal(c0['score'], c1['score'])
            and c0['peptide_sequence'] == c1['peptide_sequence']
            and c0['mass_trace'] == c1['mass_trace']):
            continue
        merged_matches.append(c1)


def merge_matches(m0, m1):
    """Merge the match information in m1 into m0."""
    for m1k, m1v in m1.iteritems():
        if m1k not in m0:
            m0[m1k] = m1v
        else:
            merge_match(m0[m1k], m1v)


def main(args=sys.argv[1:]):
    parser = optparse.OptionParser(usage=
                                   "usage: %prog [options] <result-file>..."
                                   " <output-result-file>",
                                   description=__doc__, version=__version__)
    pa = parser.add_option
    pa("-v", "--verbose", action="store_true", dest="verbose",
       help="be verbose")
    pa("--copyright", action="store_true", dest="copyright",
       help="print copyright and exit")
    (options, args) = parser.parse_args(args=args)

    if options.copyright:
        print __copyright__
        sys.exit(0)

    if len(args) < 2:
        parser.print_help()
        sys.exit(1)

    result_fn_0 = args[0]
    result_fn_1_N = args[1:-1]
    output_fn = args[-1]

    # spectrum name -> match list/spectrum info
    matches = {}
    with contextlib.closing(gzip.open(result_fn_0)) as r_file:
        r0 = cPickle.load(r_file)
    matches = r0['matches']
    total_comparisons = r0['total comparisons']
    if options.verbose:
        print >> sys.stderr, "loaded", result_fn_0

    for additional_result_fn in result_fn_1_N:
        with contextlib.closing(gzip.open(additional_result_fn)) as r1_file:
            r1 = cPickle.load(r1_file)
        check_consistency(r0, r1)
        merge_matches(matches, r1['matches'])
        total_comparisons += r1['total comparisons']
        if options.verbose:
            print >> sys.stderr, "merged", additional_result_fn

    r0['matches'] = matches
    r0['total comparisons'] = total_comparisons
    with contextlib.closing(gzip.open(output_fn, 'w')) as output_file:
        cPickle.dump(r0, output_file, cPickle.HIGHEST_PROTOCOL)
    if options.verbose:
        print >> sys.stderr, "dumped", output_fn


if __name__ == '__main__':
    main()
