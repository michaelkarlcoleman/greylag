#!/usr/bin/env python

"""
Given a set of sqt files, determine a 'valid' set of identifications that
satisfy the specified FDR (false discovery rate).  Optionally rewrite the sqt
files with validation marks set to 'N' for filtered-out spectra.

"""

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


from collections import defaultdict
import fileinput
import optparse
from pprint import pprint
import string
import sys


def warn(s):
    print >> sys.stderr, 'warning:', s
def error(s):
    sys.exit('error: ' + s)
def fileerror(s, *args):
    error(s + (", at line %s of file '%s'"
               % (fileinput.filelineno(), fileinput.filename())),
          *args)

def inplace_warning():
    warn("!!!\nan error occurred while modifying .sqt files in-place--it may"
         " be necessary to recover some or all of the .sqt files from the"
         " corresponding .sqt.bak files.\n!!!")


def reset_marks(options, sqt_fns):
    """Rewrite all evaluation marks to 'U', in-place."""

    try:
        for line in fileinput.input(sqt_fns, inplace=1, backup='.bak'):
            if line.startswith("M\t"):
                fs = line.split("\t")
                if len(fs) >= 11:
                    fs[10] = 'U' + fs[10][1:]
                line = '\t'.join(fs)
            sys.stdout.write(line)
    except:
        inplace_warning()
        raise


def mark(options, thresholds, sp_scores, sqt_fns):
    """Rewrite evaluation marks to 'N', in-place, for spectra not meeting
    score and delta thresholds."""

    spectrum_no = -1
    mark_spectrum = False

    try:
        for line in fileinput.input(sqt_fns, inplace=1, backup='.bak'):
            if line.startswith("S\t"):
                spectrum_no += 1
                charge, score, delta = sp_scores[spectrum_no]
                mark_spectrum = (charge in thresholds
                                 and score != None
                                 and (score < thresholds[charge][0]
                                      or delta < thresholds[charge][1]))
            elif line.startswith("M\t") and mark_spectrum:
                fs = line.split("\t")
                if len(fs) >= 11:
                    fs[10] = 'N' + fs[10][1:]
                line = '\t'.join(fs)
            sys.stdout.write(line)
    except:
        inplace_warning()
        raise


nulltrans = string.maketrans('','')
non_aa_chars = string.digits + string.punctuation

def stripmods(s):
    """Given a peptide, return just the unmodified peptide seen (upcased)."""
    # assuming no control chars present
    return s.upper().translate(nulltrans, non_aa_chars)


def read_sqt_info(decoy_prefix, minimum_trypticity, sqt_fns):
    """Return a pair, the first a dict mapping each charge to a list of
    (score, delta, state), where state is 'real' or 'decoy', for all the
    spectra in sqt_fns, and the second a list of (charge, score, delta), one
    for each spectrum.  Spectra lacking minimum_trypticity are dropped.
    """

    # charge -> [ (score, delta, state, spectrum_name, peptide), ... ]
    #   where state is either 'real' or 'decoy'
    #     and peptide is the peptide as first seen, with mods, no flanks
    z_scores = defaultdict(list)

    # [ (charge, score, delta), ... ]
    sp_scores = []

    current_charge = None
    current_score = None
    current_delta = 0
    current_peptide_trypticity = None   # or one of [ 0, 1, 2 ]
    current_state = set()               # 'real', 'decoy' or both
    highest_rank_seen = 0
    best_peptide_seen = None
    current_peptide = None
    spectrum_name = None

    for line in fileinput.input(sqt_fns):
        fs = line.split('\t')
        if fs[0] == 'S':
            if (current_charge != None and
                current_peptide_trypticity >= minimum_trypticity):
                if current_score != None:
                    if len(current_state) == 1:
                        z_scores[current_charge].append((current_score,
                                                         current_delta,
                                                         current_state.pop(),
                                                         spectrum_name,
                                                         current_peptide))
                    sp_scores.append((current_charge, current_score,
                                      current_delta))
            current_charge = int(fs[3])
            current_score = None
            current_delta = 0
            current_peptide_trypticity = None
            current_state = set()
            highest_rank_seen = 0
            best_peptide_seen = None
            current_peptide = None
            spectrum_name = '.'.join((fileinput.filename(), fs[1], fs[2], fs[3]))
        elif fs[0] == 'M':
            rank = int(fs[1])
            delta = float(fs[4])
            score = float(fs[5])
            sp_rank = int(fs[2])
            peptide = fs[9].strip()     # e.g., A.B@CD*.-
            mass = float(fs[3])

            if rank > highest_rank_seen + 1:
                # ignore aux top SpRank hits, because delta confounds
                continue
            highest_rank_seen = rank

            assert peptide[1] == '.' and peptide[-2] == '.'
            peptide_flanks = (peptide[0], peptide[-1]) # ('A', '-')
            peptide = peptide[2:-2]     # 'B@CD*'
            peptide_stripped = stripmods(peptide) # 'BCD'

            if best_peptide_seen == None:
                best_peptide_seen = peptide_stripped
            if current_peptide == None:
                current_peptide = peptide

            # MyriMatch gives delta=0 when score=0, even if there are better
            # scores!
            if score > 0:
                if delta == 0:
                    assert current_delta == 0
                    current_score = score
                elif (current_delta == 0 or best_peptide_seen == None
                      or best_peptide_seen == peptide_stripped):
                    current_delta = delta
                    if best_peptide_seen != peptide_stripped:
                        # once we've seen a non-equivalent peptide, don't set
                        # delta if we see further equivalent peptides
                        best_peptide_seen = '####'
            if current_peptide_trypticity == None:
                current_peptide_trypticity = 0
                if (peptide_flanks[0] in ('K', 'R')
                    and peptide_stripped[0] != 'P'):
                    current_peptide_trypticity += 1
                if (peptide_stripped[-1] in ('K', 'R')
                    and peptide_flanks[1] != 'P'):
                    current_peptide_trypticity += 1
        elif fs[0] == 'L':
            if current_delta == 0:
                if fs[1].startswith(decoy_prefix):
                    current_state.add('decoy')
                else:
                    current_state.add('real')

    # handle final spectrum, as above
    if (current_charge != None and
        current_peptide_trypticity >= minimum_trypticity):
        if current_score != None:
            if len(current_state) == 1:
                z_scores[current_charge].append((current_score,
                                                 current_delta,
                                                 current_state.pop(),
                                                 spectrum_name,
                                                 current_peptide))
            sp_scores.append((current_charge, current_score,
                              current_delta))

    return (z_scores, sp_scores)


def PPV(reals, decoys):
    """Returns the estimated Positive Predictive Value (== 1 - FDR), given
    counts of reals and decoys."""

    # We know that the decoys are false positives, and we estimate that an
    # equal number of the "reals" are actually false, too.
    false_positives = 2*decoys
    true_positives = reals - decoys

    return float(true_positives) / (true_positives + false_positives)


def calculate_inner_threshold(fdr, charge, spinfo):
    ppv_goal = 1 - fdr

    spinfo = sorted(spinfo, key=lambda x: x[0])

    epsilon = +1e-6

    real_count = sum(1 for x in spinfo if x[2] == 'real')
    decoy_count = len(spinfo) - real_count

    current_threshold = -1e100      # allow all spectra
    for n, sp in enumerate(spinfo):
        if real_count == 0:
            return (None, real_count, decoy_count) # give up

        ppv_est = PPV(real_count, decoy_count)
        if ppv_est >= ppv_goal:
            break
        if sp[2] == 'real':
            real_count -= 1
        else:
            decoy_count -= 1
        # set threshold just high enough to exclude this spectrum
        current_threshold = sp[0] + epsilon
    else:
        current_threshold = spinfo[-1][0] + epsilon # couldn't meet goal

    return (current_threshold, real_count, decoy_count)


def calculate_combined_thresholds(options, z_scores):
    """Find best score/delta thresholds for each charge."""

    ppv_goal = 1 - options.fdr

    # Rather than search every possible value of delta, we're only going to
    # "sample" at this granularity.  This cuts search time dramatically (and
    # makes it O(n) instead of O(n**2).  Extra precision wouldn't really be
    # useful in any case.
    SEARCH_GRANULARITY = 0.001

    # charge -> (score, delta, passing_reals, passing_decoys)
    thresholds = {}

    total_reals, total_decoys = 0, 0

    for charge, spinfo in z_scores.iteritems():
        spinfo0 = sorted(spinfo, key=lambda x: x[1], reverse=True)

        last_value = None

        while spinfo0:
            this_value = spinfo0[-1][1] # current delta
            if (last_value == None
                or abs(this_value - last_value) >= SEARCH_GRANULARITY):

                # "inner" is score
                r = calculate_inner_threshold(options.fdr, charge, spinfo0)
                if r[0] != None:
                    if options.debug:
                        print '#', charge, r[0], this_value, r[1], r[2]
                    if (charge not in thresholds
                        or r[1] > thresholds[charge][2]):
                        thresholds[charge] = (r[0], this_value, r[1], r[2])

                last_value = this_value
            spinfo0.pop()

        if charge in thresholds:
            score_threshold, delta_threshold = thresholds[charge][0:2]

            if options.list_peptides:
                for (score, delta, state,
                     spectrum_name, peptide) in z_scores[charge]:
                    if score >= score_threshold and delta >= delta_threshold:
                        print "%+d %s %s %s" % (charge, spectrum_name, state,
                                                peptide)
            else:
                reals, decoys = thresholds[charge][2], thresholds[charge][3]
                total_reals += reals
                total_decoys += decoys
                print ("%+d: score %s, delta %s -> %s real ids, %s decoys"
                       " (fdr %.4f)"
                       % (charge, score_threshold, delta_threshold,
                          reals, decoys,
                          1 - PPV(thresholds[charge][2],
                                  thresholds[charge][3])))
        else:
            warn("could not calculate thresholds for %+d" % charge)

    if not options.list_peptides:
        print "# total: %s real ids, %s decoys" % (total_reals, total_decoys)

    return thresholds


def main(args=sys.argv[1:]):
    parser = optparse.OptionParser(usage=
                                   "usage: %prog [options] <sqt-file>...",
                                   description=__doc__, version=__version__)
    pa = parser.add_option
    DEFAULT_DECOY_PREFIX = "SHUFFLED_"
    pa("--decoy-prefix", dest="decoy_prefix", default=DEFAULT_DECOY_PREFIX,
       help='prefix given to locus name of decoy (e.g., shuffled) database'
       ' sequences [default=%r]' % DEFAULT_DECOY_PREFIX, metavar="PREFIX")
    DEFAULT_FDR = 0.02
    pa("--fdr", dest="fdr", type="float", default=DEFAULT_FDR,
       help="false discovery rate [default=%s]" % DEFAULT_FDR,
       metavar="PROPORTION")
    pa("-v", "--verbose", action="store_true", dest="verbose",
       help="be verbose")
    pa("--list-peptides", action="store_true", dest="list_peptides",
       help="print information about passing peptides")
    pa("-t", "--minimum-trypticity", type="choice", choices=('0', '1', '2'),
       default='0', dest="minimum_trypticity",
       help="drop peptides with too few tryptic ends"
       " (2=fully tryptic, 1=semi-tryptic, 0=any)", metavar="TRYPTICITY")
    #pa("--graph", dest="graph",
    #   help='create distribution graphs, using the specified prefix',
    #   metavar="PATH PREFIX")
    pa("--debug", action="store_true", dest="debug",
       help="show debug output")
    pa("-m", "--mark", action="store_true", dest="mark",
       help="rewrite the input files, changing some validation marks to 'N',"
       " according to filtering")
    pa("--reset-marks", action="store_true", dest="reset_marks",
       help="rewrite the input files, changing all validation marks to 'U'")
    pa("--copyright", action="store_true", dest="copyright",
       help="print copyright and exit")
    (options, args) = parser.parse_args(args=args)

    if options.copyright:
        print __copyright__
        sys.exit(0)

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    if not (0.0 <= options.fdr <= 1.0):
        error("--fdr must be within range [0.0, 1.0]")
    if options.mark and options.reset_marks:
        error("only one of --mark and --reset-marks may be specified")

    if options.reset_marks:
        reset_marks(options, args)
        return

    z_scores, spectrum_scores = read_sqt_info(options.decoy_prefix,
                                              int(options.minimum_trypticity),
                                              args)

    if options.debug:
        print ('%s passing spectra, of which %s are not decoys'
               % (len(spectrum_scores), sum(len(z_scores[x])
                                            for x in z_scores)))

    thresholds = calculate_combined_thresholds(options, z_scores)

    if options.debug:
        pprint(thresholds)

        for charge, score, delta in spectrum_scores:
            print '#S', charge, score, delta

    if options.mark:
        mark(options, thresholds, spectrum_scores, args)


if __name__ == '__main__':
    main()
