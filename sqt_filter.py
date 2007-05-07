#!/usr/bin/env greylag-python

"""
Filter a set of sqt files according to FPR and other criteria, optionally
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


def read_sqt_info(decoy_prefix, sqt_fns):
    """Return a dict mapping each charge to a list of (score, delta, state),
    where state is 'real' or 'decoy', for all the spectra in sqt_fns.
    """

    # charge -> [ (score, delta, state), ... ]
    #   where state is either 'real' or 'decoy'
    z_scores = defaultdict(list)

    # reset fileinput first, in case it's been called before
    try:
        fileinput.close()
    except RuntimeError:
        pass

    current_charge = None
    current_score = None
    current_delta = 0
    current_state = set()

    for line in fileinput.input(sqt_fns):
        fs = line.split('\t')
        if fs[0] == 'S':
            if current_score and len(current_state) == 1:
                z_scores[current_charge].append((current_score, current_delta,
                                                 current_state.pop()))
            current_charge = int(fs[3])
            current_score = None
            current_delta = 0
            current_state = set()
        elif fs[0] == 'M':
            delta, score = float(fs[4]), float(fs[5])
            if delta == 0:
                current_score = score
            elif current_delta == 0:
                current_delta = delta
        elif fs[0] == 'L':
            if current_delta == 0:
                if fs[1].startswith(decoy_prefix):
                    current_state.add('decoy')
                else:
                    current_state.add('real')
    # handle final spectrum
    if current_score and len(current_state) == 1:
        z_scores[current_charge].append((current_score, current_delta,
                                         current_state.pop()))

    #pprint(dict(z_scores))
    return z_scores


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


def calculate_combined_thresholds(options, sqt_fns):
    """Find best score/delta thresholds for each charge.
    """

    z_scores = read_sqt_info(options.decoy_prefix, sqt_fns)

    specificity_goal = 1 - options.fpr

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

        print ("%+d: score %s, delta %s -> %s real ids (fdr %.4f)"
               % (charge, thresholds[charge][0], thresholds[charge][1],
                  thresholds[charge][2],
                  1 - specificity(thresholds[charge][2],
                                  thresholds[charge][3])))

    return thresholds


def main(args=sys.argv[1:]):
    parser = optparse.OptionParser(usage=
                                   "usage: %prog [options] <sqt-file>...",
                                   description=__doc__, version=__version__)
    pa = parser.add_option
    pa("--decoy-prefix", dest="decoy_prefix", default="SHUFFLED_",
       help='prefix given to locus name of decoy (e.g., shuffled) database'
       ' sequences [default="SHUFFLED_"]', metavar="PREFIX")
    pa("--fpr", dest="fpr", type="float", default="0.02",
       help="false positive rate [default=0.02]", metavar="PROPORTION")
    pa("-v", "--verbose", action="store_true", dest="verbose",
       help="be verbose")
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

    if not (0.0 <= options.fpr <= 1.0):
        error("--fpr must be within range [0.0, 1.0]")
    if options.mark and options.reset_marks:
        error("only one of --mark and --reset-marks may be specified")

    if options.reset_marks:
        # do it
        error("not yet implemented")
        return

    thresholds = calculate_combined_thresholds(options, args)

    if options.debug:
        pprint(thresholds)

    if options.mark:
        # do it
        error("not yet implemented")


if __name__ == '__main__':
    main()
