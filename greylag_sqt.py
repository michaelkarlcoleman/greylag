#!/usr/bin/env python

'''Generate a set of sqt files from the specified greylag search result
file.  The sqt files will be created in the current directory, with names
corresponding to the searched ms2 files.
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


from collections import defaultdict
import contextlib
import cPickle
import math
import optparse
import os.path
from pprint import pprint
import sys

import greylag


def warn(s):
    print >> sys.stderr, 'warning:', s
def error(s):
    sys.exit('error: ' + s)


def generate_leaves(tree):
    """Given a tree made of lists, yield the leaves.

    >>> list(generate_leaves([(1,2,3)]))
    [(1, 2, 3)]
    >>> list(generate_leaves([[[(1,2,3)]], (4,5)]))
    [(1, 2, 3), (4, 5)]

    """
    if not isinstance(tree, list):
        yield tree
        return
    for b in tree:
        for l in generate_leaves(b):
            yield l


# Example legacy headers:
# H       Database        /foo/bar/db.fasta
# H       PrecursorMasses AVG
# H       FragmentMasses  MONO
# H       StaticMod       C=160.1388
# H       DiffMod M*=+16.0
# H       DiffMod STY#=+80.0

def print_legacy_headers(f, r):
    """Output legacy H lines that readers of traditional SQT files may
    expect.  Warns if multiple databases, or multiple or isotopic mass regimes
    are present.
    """

    if len(r['databases']) > 1:
        warn("multiple databases present")
    for db in r['databases']:
        print >> f, "H\tDatabase\t%s" % db

    precursor_masses, fragment_masses = 'MONO', 'MONO'
    mass_regimes = r['parameters']['mass_regimes']
    if len(mass_regimes) > 1:
        warn("multiple mass regimes present--legacy SQT headers will be based"
             "on just the first")
    parent_regime, fragment_regime = mass_regimes[0]
    if parent_regime[1] or fragment_regime[1]:
        warn("isotopic mass regimes present--ignored for legacy SQT headers")
    assert parent_regime[0] in ('MONO', 'AVG')
    assert fragment_regime[0] in ('MONO', 'AVG')
    print >> f, "H\tPrecursorMasses\t%s" % parent_regime[0]
    print >> f, "H\tFragmentMasses\t%s" % fragment_regime[0]

    # FIX: merge this with mass calculation code in
    #      greylag.initialize_spectrum_parameters?
    # FIX: currently assuming MONO only (what would be better?)
    static = dict((res, greylag.formula_mass(greylag.RESIDUE_FORMULA[res]))
                  for res in greylag.RESIDUES)
    for m in r['parameters']['pervasive_mods']:
        if isinstance(m[1], float):
            static[m[3]] += m[0] * m[1]
        else:
            static[m[3]] += m[0] * greylag.formula_mass(m[1])
    # all residues get a StaticMod because our figures are more exact than the
    # "default" ones
    for res in static:
        print >> f, "H\tStaticMod\t%s=%s" % (res, static[res])

    diffs = set(generate_leaves(r['parameters']['potential_mods']))
    for m in diffs:
        symbol = m[5]
        if not symbol:
            warn("unmarked potential mod present--omitting")
            continue
        if isinstance(m[1], float):
            delta = m[0] * m[1]
        else:
            delta = m[0] * greylag.formula_mass(m[1])
        print >> f, "H\tDiffMod\t%s%s=%+f" % (m[3], m[5], delta)


def print_header(f, r):
    """Output H lines."""
    print >> f, "H\tSQTGenerator\tgreylag"
    print >> f, "H\tSQTGeneratorVersion\t%s" % greylag.VERSION

    print_legacy_headers(f, r)

    for pk, pv in sorted(r['parameters'].items()):
        print >> f, "H\tParameter\t%s\t%s" % (pk, pv)


def print_regime_manifest(f, regime_manifest):
    """Output R lines."""
    for regime_no, residue, mass in regime_manifest:
        print >> f, "R\t%s\t%s\t%s" % (regime_no, residue, mass)


def generate_marked_sequence(match_name_map, mass_trace, peptide_sequence):
    """Yield the characters in a marked version of peptide_sequence.

    >>> ''.join(generate_marked_sequence({}, [], 'ASDF'))
    'ASDF'
    >>> ''.join(generate_marked_sequence({('S', 80) : ('phosphorylation', '*')},
    ...                                  [{'position' : 1, 'delta' : 80}],
    ...                                  'ASDF'))
    'AS*DF'

    """

    # FIX: handle [ ]
    trace = sorted(mass_trace,
                   key=lambda x: (x['position'], x['delta']), reverse=True)
    for n, r in enumerate(peptide_sequence):
        yield r
        while trace and trace[-1]['position'] == n:
            mark = match_name_map[(r, trace[-1]['delta'])][1]
            if mark:                    # None if no mark char for this mod
                yield mark
            trace.pop()


def print_spectrum(f, mod_name_map, sp_name, sp_matches, enhanced=False):
    """Print the lines (S/M/L/A*) associated with the given spectrum."""
    spectrum, sp_best_matches = sp_matches

    scan_low, scan_high, charge = sp_name.rsplit('.', 2)
    print >> f, '\t'.join(str(v) for v in
                          ["S", scan_low, scan_high, charge, 0, 'honk',
                           spectrum['mass'],
                           round(math.log(spectrum['total_ion_current']), 4),
                           0, spectrum['comparisons']])

    best_scores = sorted(list(set(m['score'] for m in sp_best_matches)))
    if not best_scores:
        return

    # score -> rank
    rank_map = dict(zip(best_scores, range(1,len(best_scores)+1)))

    best_score = best_scores[0]
    for match in sp_best_matches:
        score = match['score']
        if score == 0:
            continue                    # null entry
        rank = rank_map[score]
        score_delta = (best_score - score) / best_score

        match_name_map = mod_name_map[(match.get('mass_regime_index', 0),
                                       match.get('conjunct_index', 0))]
        marked_sequence \
            = ''.join(generate_marked_sequence(match_name_map,
                                               match.get('mass_trace', []),
                                               match['peptide_sequence']))

        # FIX: also need M lines for fixed mods, terminal mods, isotope mods
        # (N15), etc.
        # Note: scores, being log probability, are non-positive, but we
        # flip the sign in the SQT output
        print >> f, '\t'.join(str(v) for v in
                              ["M", rank, rank,
                               round(match['predicted_parent_mass'], 5),
                               round(score_delta, 4), round(-score, 4), 0,
                               # FIX: 1 of 2 ions found--keep DTASelect happy
                               1, 2,
                               '.'.join((match['N_peptide_flank'],
                                         marked_sequence,
                                         match['C_peptide_flank'])),
                               'U'])

        if enhanced:
            if match.get('mass_regime_index', 0) != 0:
                print >> f, "AR\t%s" % match['mass_regime_index']
            if match.get('pca_delta', 0) != 0:
                print >> f, "APCA\t%s" % match['pca_delta']
            for mt in match.get('mass_trace', []):
                name = (match_name_map[(match['peptide_sequence'][mt['position']],
                                        mt['delta'])][0])
                fs = ["AM", mt['position'], round(mt['delta'], 5),
                      name if name else '']
                print >> f, '\t'.join(str(v) for v in fs)

        assert len(match['sequence_name']) == len(match['peptide_begin'])
        for sn, pb in zip(match['sequence_name'], match['peptide_begin']):
            print >> f, 'L\t%s\t%s' % (sn, pb)


def main(args=sys.argv[1:]):
    parser = optparse.OptionParser(usage=
                                   "usage: %prog [options] <result-file>",
                                   description=__doc__, version=greylag.VERSION)
    pa = parser.add_option
    pa("-d", "--output-directory", dest="output_directory",
       help="directory where output files are written [default is '.']",
       metavar="DIR")
    pa("-e", "--enhanced-output", action="store_true", dest="enhanced_output",
       help="include extra DTASelect-incompatible information in output")
    pa("-v", "--verbose", action="store_true", dest="verbose",
       help="be verbose")
    pa("--dump", action="store_true", dest="dump",
       help="just dump the result file to stdout (for debugging)")
    pa("--copyright", action="store_true", dest="copyright",
       help="print copyright and exit")
    (options, args) = parser.parse_args(args=args)

    if options.copyright:
        print __copyright__
        sys.exit(0)

    if len(args) != 1:
        parser.print_help()
        sys.exit(1)

    with contextlib.closing(open(args[0])) as r_file:
        r = cPickle.load(r_file)

    if options.dump:
        pprint(r)
        return

    spectrum_fns = r['spectrum files']
    assert len(spectrum_fns) == len(set(spectrum_fns)) # check uniqueness

    matches = r['matches'].items()
    matches.sort()

    # (regime index, conjunct index) -> (residue, delta)
    #                                                -> (mod name, mod marker)
    mod_name_map = defaultdict(dict)
    for regime in range(len(r['mass regime atomic masses'])):
        for cj_n, (N_cj, C_cj, R_cj) in enumerate(r['modification conjuncts']):
            assert not N_cj and not C_cj, "not yet implemented"
            for cj in R_cj:
                for res in cj[3]:
                    mod_name_map[(regime, cj_n)][(res, cj[6][regime][1])] \
                        = (cj[4], cj[5])
    #pprint(mod_name_map.items())

    for spectrum_n, spectrum_fn in enumerate(spectrum_fns):
        assert os.path.dirname(spectrum_fn) == ''
        sqt_fn = os.path.splitext(spectrum_fn)[0] + '.sqt'
        with open(sqt_fn, 'w') as sqtf:
            print_header(sqtf, r)
            print_regime_manifest(sqtf, r['mass regime manifest'])
            for match in matches:
                if match[0][0] == spectrum_n:
                    print_spectrum(sqtf, mod_name_map, match[0][1], match[1],
                                   options.enhanced_output)


if __name__ == '__main__':
    main()
