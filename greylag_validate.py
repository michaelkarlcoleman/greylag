#!/usr/bin/env greylag-python

"""
Filter a set of sqt files according to FDR and other criteria, optionally
rewriting them with the validation marks set to 'N' for filtered-out spectra.

"""

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


from collections import defaultdict
import fileinput
import optparse
from pprint import pprint
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


def read_sqt_info(decoy_prefix, is_laf, sqt_fns):
    """Return a pair, the first a dict mapping each charge to a list of
    (score, delta, state), where state is 'real' or 'decoy', for all the
    spectra in sqt_fns, and the second a list of (charge, score, delta), one
    for each spectrum.
    """

    # charge -> [ (score, delta, state), ... ]
    #   where state is either 'real' or 'decoy'
    z_scores = defaultdict(list)

    # [ (charge, score, delta), ... ]
    sp_scores = []

    current_charge = None
    current_score = None
    current_delta = 0
    current_sp_rank = None
    current_peptide_length = None
    current_peptide_trypticity = None   # or one of [ 0, 1, 2 ]
    current_state = set()

    for line in fileinput.input(sqt_fns):
        fs = line.split('\t')
        if fs[0] == 'S':
            if (current_charge != None and
                (not is_laf or (current_sp_rank <= 10
                                and current_peptide_trypticity == 2))):
                if current_score != None and len(current_state) == 1:
                    z_scores[current_charge].append((current_score,
                                                     current_delta,
                                                     current_state.pop()))
                if current_charge != None:
                    sp_scores.append((current_charge, current_score,
                                      current_delta))
            current_charge = int(fs[3])
            current_score = None
            current_delta = 0
            current_sp_rank = None
            current_peptide_length = None
            current_peptide_trypticity = None
            current_state = set()
        elif fs[0] == 'M':
            delta, score, sp_rank, peptide = (float(fs[4]), float(fs[5]),
                                              int(fs[2]), fs[9].strip())
            pep_length = len(peptide) - 4

            if delta == 0:
                current_score = score
            elif current_delta == 0:
                current_delta = delta
            if current_sp_rank == None:
                current_sp_rank = sp_rank
            if current_peptide_length == None:
                current_peptide_length = pep_length
            if current_peptide_trypticity == None:
                current_peptide_trypticity = 0
                if (peptide[0] in ('K', 'R')
                    and peptide[2] != 'P'):
                    current_peptide_trypticity += 1
                if (peptide[-3] in ('K', 'R')
                    and peptide[-1] != 'P'):
                    current_peptide_trypticity += 1
        elif fs[0] == 'L':
            if current_delta == 0:
                if fs[1].startswith(decoy_prefix):
                    current_state.add('decoy')
                else:
                    current_state.add('real')
    # handle final spectrum, as above
    if (current_charge != None and
        (not is_laf or (current_sp_rank <= 10
                        and current_peptide_trypticity == 2))):
        if current_score != None and len(current_state) == 1:
            z_scores[current_charge].append((current_score, current_delta,
                                             current_state.pop()))
        if current_charge != None:
            sp_scores.append((current_charge, current_score, current_delta))

    return (z_scores, sp_scores)


def specificity(positives, negatives):
    return (float(positives - negatives)
            / (positives + negatives))


def calculate_inner_threshold(specificity_goal, charge, spinfo):
    spinfo = sorted(spinfo, key=lambda x: x[0])

    real_count = sum(1 for x in spinfo if x[-1] == 'real')
    decoy_count = len(spinfo) - real_count

    if real_count == 0:
        return (None, real_count, decoy_count) # give up

    current_threshold = -1e100      # allow all spectra
    for n, sp in enumerate(spinfo):
        specificity_est = specificity(real_count, decoy_count)
        if specificity_est >= specificity_goal:
            break
        if sp[-1] == 'real':
            real_count -= 1
        else:
            decoy_count -= 1
        # set threshold just high enough to exclude this spectrum
        current_threshold = sp[0] + 1e-6
    else:
        current_threshold = spinfo[-1][0] + 1e-6 # couldn't meet goal

    return (current_threshold, real_count, decoy_count)


def calculate_combined_thresholds(options, z_scores):
    """Find best score/delta thresholds for each charge."""

    specificity_goal = 1 - options.fdr

    # Rather than search every possible value of delta, we're only going to
    # "sample" at this granularity.  This cuts search time dramatically (and
    # making it O(n) instead of O(n**2).  Extra precision wouldn't really be
    # useful in any case.
    SEARCH_GRANULARITY = 0.001

    # charge -> (score, delta, passing_reals, passing_decoys)
    thresholds = {}

    for charge, spinfo in z_scores.iteritems():
        spinfo0 = sorted(spinfo, key=lambda x: x[1], reverse=True)

        last_value = None

        while spinfo0:
            this_value = spinfo0[-1][1] # current delta
            if (last_value == None
                or abs(this_value - last_value) >= SEARCH_GRANULARITY):

                # "inner" is score
                r = calculate_inner_threshold(specificity_goal, charge,
                                              spinfo0)
                if r[0] != None:
                    if options.debug:
                        print '#', charge, r[0], this_value, r[1], r[2]
                    if (charge not in thresholds
                        or r[1] > thresholds[charge][2]):
                        thresholds[charge] = (r[0], this_value, r[1], r[2])

                last_value = this_value
            spinfo0.pop()

        if charge in thresholds:
            if options.verbose:
                print ("%+d: score %s, delta %s -> %s real ids (fdr %.4f)"
                       % (charge, thresholds[charge][0], thresholds[charge][1],
                          thresholds[charge][2],
                          1 - specificity(thresholds[charge][2],
                                          thresholds[charge][3])))
        else:
            warn("could not calculate thresholds for %+d" % charge)

    return thresholds


def main(args=sys.argv[1:]):
    parser = optparse.OptionParser(usage=
                                   "usage: %prog [options] <sqt-file>...",
                                   description=__doc__, version=__version__)
    pa = parser.add_option
    pa("--decoy-prefix", dest="decoy_prefix", default="SHUFFLED_",
       help='prefix given to locus name of decoy (e.g., shuffled) database'
       ' sequences [default="SHUFFLED_"]', metavar="PREFIX")
    pa("--fdr", dest="fdr", type="float", default="0.02",
       help="false discovery rate [default=0.02]", metavar="PROPORTION")
    pa("-v", "--verbose", action="store_true", dest="verbose",
       help="be verbose")
    pa("--laf", action="store_true", dest="laf",
       help="drop peptides of length < 7 or Sp rank > 10")
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
                                              options.laf, args)

    thresholds = calculate_combined_thresholds(options, z_scores)

    if options.debug:
        pprint(thresholds)

        for charge, score, delta in spectrum_scores:
            print '#S', charge, score, delta

    if options.mark:
        mark(options, thresholds, spectrum_scores, args)


if __name__ == '__main__':
    main()
