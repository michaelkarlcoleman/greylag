#!/usr/bin/env python

"""
Given a set of sqt files, determine a 'valid' set of identifications that
satisfy the specified FDR (false discovery rate) and other specified criteria.
Optionally rewrite the sqt files with invalid spectra marked 'N' or deleted
altogether.

Note that a peptide end with an invalid flanking residue (currently any of
BJOUXZ) is not considered tryptic if it's an N-terminal end and is considered
tryptic if it's a C-terminal end and otherwise qualifies.  (Flanking '-' is
always considered tryptic.)

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


from collections import defaultdict
import optparse
import os
from pprint import pprint
import re
import string
import sys
import time


# allow this script to be run even if not "installed"
try:
    from greylag import VERSION
except:
    VERSION = 'unknown'


def warn(s):
    print >> sys.stderr, 'warning:', s
def error(s):
    sys.exit('error: ' + s)


def generate_spectra_from_files(sqt_filenames):
    """Yield tuples (filename, H+_lines, S_line, [(M_line, mod_lines, L_lines),
    ...]) for the given SQT files.  H+_lines, for example, is a list
    containing H lines and possibly other kinds of lines up to the next S
    line.  (mod_lines are non-L lines following an M-line.)
    """
    for sqt_filename in sqt_filenames:
        with open(sqt_filename) as sqtf:
            ###FIX###sqt_filename = os.path.splitext(sqt_filename)[0]
            H_lines = []
            S_line, ML_lines = None, []
            for line in sqtf:
                if line.startswith('L\t'):
                    ML_lines[-1][2].append(line)
                elif line.startswith('M\t'):
                    ML_lines.append((line, [], []))
                elif line.startswith('S\t'):
                    if S_line:
                        yield (sqt_filename, H_lines, S_line, ML_lines)
                        S_line, ML_lines = None, []
                    S_line = line
                else:
                    if not S_line:
                        # read headers up to first S line
                        H_lines.append(line)
                    else:
                        # other kind of line following M line
                        ML_lines[-1][1].append(line)
            if S_line:
                yield (sqt_filename, H_lines, S_line, ML_lines)


def write_spectra_to_files(spectra):
    """Given a generator of tuples of the form generated by
    'generate_spectra_from_files', write them to the implied files.  If any of
    those files already exist, they are renamed with a '.bak' suffix.  Temp
    files are written first, then renamed at the end.
    """

    tmp_suffix = '.tmp.%s' % time.time()

    # filename -> file
    filemap = {}

    for sqt_filename, H_lines, S_line, ML_lines in spectra:
        sqtf = filemap.get(sqt_filename)
        if sqtf == None:
            sqtf = filemap[sqt_filename] = open(sqt_filename + tmp_suffix, 'w')
            for h in H_lines:
                sqtf.write(h)
        sqtf.write(S_line)
        for m_line, mod_lines, l_lines in ML_lines:
            sqtf.write(m_line)
            for m in mod_lines:
                sqtf.write(m)
            for l in l_lines:
                sqtf.write(l)

    for fd in filemap.values():
        fd.close()
    for fn in filemap:
        if os.path.exists(fn):
            os.rename(fn, fn + '.bak')
        os.rename(fn + tmp_suffix, fn)


def reset_marks(spectra):
    """Given a generator of tuples of the form generated by
    'generate_spectra_from_files', yield them with the M line marks
    (destructively) reset to U.
    """

    def reset_mark(m_line):
        return re.sub(r'\S(\s*)$', r'U\1', m_line, 1)

    for sqt_filename, H_lines, S_line, ML_lines in spectra:
        ML_lines_1 = [ (reset_mark(M_line), mod_lines, L_lines)
                       for M_line, mod_lines, L_lines in ML_lines ]
        yield (sqt_filename, H_lines, S_line, ML_lines_1)


def set_marks(options, thresholds, sp_scores, spectra, kill=False):
    """Given a generator of tuples 'spectra' of the form generated by
    'generate_spectra_from_files', yield them with the M line marks
    (destructively) set to 'N', for spectra not meeting score and delta
    thresholds.  If 'kill', omit marked spectra altogether.
    """

    def set_mark(m_line):
        return re.sub(r'\S(\s*)$', r'N\1', m_line, 1)

    for spectrum_no, (sqt_filename, H_lines, S_line, ML_lines) \
            in enumerate(spectra):
        mark_spectrum = True
        if sp_scores[spectrum_no]:
            charge, score, delta = sp_scores[spectrum_no]
            mark_spectrum = (charge not in thresholds
                             or (score < thresholds[charge][0]
                                 or delta < thresholds[charge][1]))
        if mark_spectrum and not kill:
            ML_lines = [ (set_mark(M_line), mod_lines, L_lines)
                         for M_line, mod_lines, L_lines in ML_lines ]
        if not (kill and mark_spectrum):
            yield (sqt_filename, H_lines, S_line, ML_lines)


def get_spectrum_info(decoy_prefix, minimum_trypticity, spectra):
    """Return a pair, the first a dict mapping each charge to a list of
    (score, delta, state), where state is 'real' or 'decoy', for all the
    spectra in sqt_fns, and the second a list of (charge, score, delta), one
    for each spectrum, given a generator of tuples 'spectra' of the form
    generated by 'generate_spectra_from_files'.  Spectra lacking
    minimum_trypticity are dropped.
    """

    # charge -> [ (score, delta, state, spectrum_name, peptide), ... ]
    #   where state is either 'real' or 'decoy'
    #     and peptide is the peptide as first seen, with mods, no flanks
    z_scores = defaultdict(list)

    # information about spectra, in file order--None iff spectrum was filtered
    # out
    # [ (charge, score, delta) or None, ... ]
    sp_scores = []

    for sqt_filename, H_lines, S_line, ML_lines in spectra:
        S_fields = S_line.split('\t')
        current_charge = int(S_fields[3])
        spectrum_name = '.'.join((sqt_filename,
                                  S_fields[1], S_fields[2], S_fields[3]))

        current_score = None
        current_delta = None
        current_peptide_trypticity = None
        current_state = set()
        highest_rank_seen = 0
        best_peptide_seen = None
        current_peptide = None

        for M_line, mod_lines, L_lines in ML_lines:
            M_fields = M_line.split('\t')

            rank = int(M_fields[1])
            score = float(M_fields[5])
            peptide = M_fields[9].strip() # e.g., A.B@CD*.-
            #delta = float(M_fields[4])
            #sp_rank = int(M_fields[2])
            #mass = float(M_fields[3])

            if rank > highest_rank_seen + 1:
                # ignore aux top SpRank hits, because delta confounds
                break
            highest_rank_seen = rank

            if current_score == None:
                if score <= 0:
                    break               # spectrum is total crap
                current_score = score

            # Don't trust delta from input file.
            delta = (current_score - score) / current_score

            assert peptide[1] == '.' and peptide[-2] == '.'
            peptide_flanks = (peptide[0], peptide[-1]) # ('A', '-')
            peptide = peptide[2:-2]                    # 'B@CD*'

            if current_peptide == None:
                current_peptide = peptide
                current_peptide_trypticity = 0
                # NB: K*, R* not tryptic!
                if (peptide_flanks[0] in ('K', 'R') and peptide[0] != 'P'
                    or peptide_flanks[0] == '-'):
                    current_peptide_trypticity += 1
                if (peptide[-1] in ('K', 'R') and peptide_flanks[1] != 'P'
                    or peptide_flanks[1] == '-'):
                    current_peptide_trypticity += 1

            # Choose delta from first M line such that peptide differs from
            # peptide in initial (top) M line.  We consider differently
            # modified versions of the same peptide to be different (at least
            # for now).
            if peptide != current_peptide:
                current_delta = delta
                break

            for L_line in L_lines:
                if L_line.split('\t')[1].startswith(decoy_prefix):
                    current_state.add('decoy')
                else:
                    current_state.add('real')

        sps = None
        if None not in (current_score, current_delta):
            if current_peptide_trypticity >= minimum_trypticity:
                if len(current_state) == 1: # real xor decoy
                    z_scores[current_charge].append((current_score,
                                                     current_delta,
                                                     current_state.pop(),
                                                     spectrum_name,
                                                     current_peptide))
                    #print >> sys.stderr, z_scores[current_charge][-1]
                sps = (current_charge, current_score, current_delta)
        sp_scores.append(sps)

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


def calculate_combined_thresholds(fdr, options, z_scores):
    """Find best score/delta thresholds for each charge."""

    ppv_goal = 1 - fdr

    # Rather than search every possible value of delta, we're only going to
    # "sample" at this granularity.  This cuts search time dramatically (and
    # makes it O(n) instead of O(n**2).  Extra precision wouldn't really be
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
                r = calculate_inner_threshold(fdr, charge, spinfo0)
                if r[0] != None:
                    if options.debug:
                        print '#', charge, r[0], this_value, r[1], r[2]
                    if (charge not in thresholds
                        or r[1] > thresholds[charge][2]):
                        thresholds[charge] = (r[0], this_value, r[1], r[2])

                last_value = this_value
            spinfo0.pop()

    return thresholds


def search_adjusting_fdr(options, z_scores):
    """Repeatedly calculate thresholds while adjusting FDR, looking for an FDR
    guess that results in an appropriate goal FDR.  This is needed because the
    thresholds found for per-charge FDRs will be saturated or have
    quantization errors, resulting in overly conservative thresholds.
    """

    def try_fdr(fdr_guess):
        thresholds = calculate_combined_thresholds(fdr_guess, options,
                                                   z_scores)
        charges = sorted(thresholds.keys())
        total_reals = sum(thresholds[charge][2] for charge in charges)
        total_decoys = sum(thresholds[charge][3] for charge in charges)
        fdr_result = 1 - PPV(total_reals, total_decoys)
        return fdr_result, thresholds, charges, total_reals, total_decoys

    low_fdr, high_fdr = options.fdr, min(2*options.fdr, 0.4999999999999)
    low_results = try_fdr(low_fdr)
    high_results = try_fdr(high_fdr)
    initial_total_reals = low_results[3]

    #    #
    #    l        g            h
    #   l       g         h
    #    #

    for i in range(32):
        if abs(high_fdr - low_fdr) <= 0.01*options.fdr:
            break
        guess_fdr = (high_fdr + low_fdr) / 2.0
        guess_results = try_fdr(guess_fdr)
        if options.debug:
            print ("%s %s %s -> %s"
                   % (low_fdr, guess_fdr, high_fdr, guess_results[0]))
        if guess_results[0] > options.fdr:
            high_fdr, high_results = guess_fdr, guess_results
        else:
            low_fdr, low_results = guess_fdr, guess_results
    else:
        warn("FDR adjustment convergence failed")

    fdr_result, thresholds, charges, total_reals, total_decoys = low_results

    for charge in charges:
        score_threshold, delta_threshold = thresholds[charge][0:2]

        if options.list_peptides:
            for (score, delta, state,
                 spectrum_name, peptide) in z_scores[charge]:
                if score >= score_threshold and delta >= delta_threshold:
                    print "%+d %s %s %s" % (charge, spectrum_name, state,
                                            peptide)
        else:
            reals, decoys = thresholds[charge][2], thresholds[charge][3]
            print ("%+d: score %f, delta %.6f -> %s real ids, %s decoys"
                   " (fdr %.4f)"
                   % (charge, score_threshold, delta_threshold,
                      reals, decoys, 0.5 * (1 - PPV(reals, decoys))))

    if not options.list_peptides:
        print ("# total: %s real ids, %s decoys (fdr %.4f)"
               % (total_reals, total_decoys, fdr_result / 2.0))
        print ("# (FDR adjustment found %s more real ids)"
               % (total_reals - initial_total_reals))

    return thresholds


def main(args=sys.argv[1:]):
    parser = optparse.OptionParser(usage=
                                   "usage: %prog [options] <sqt-file>...",
                                   description=__doc__, version=VERSION)
    pa = parser.add_option
    DEFAULT_DECOY_PREFIX = "SHUFFLED_"
    pa("--decoy-prefix", dest="decoy_prefix", default=DEFAULT_DECOY_PREFIX,
       help='prefix given to locus name of decoy (e.g., shuffled) database'
       ' sequences [default=%r]' % DEFAULT_DECOY_PREFIX, metavar="PREFIX")
    DEFAULT_FDR = 0.01
    pa("--fdr", dest="fdr", type="float", default=DEFAULT_FDR,
       help="false discovery rate [default=%s] (Note that this is the"
       " final resulting FDR after the decoys have been removed" % DEFAULT_FDR,
       metavar="PROPORTION")
    pa("-v", "--verbose", action="store_true", dest="verbose",
       help="be verbose")
    pa("--list-peptides", action="store_true", dest="list_peptides",
       help="print information about passing peptides")
    pa("-t", "--minimum-trypticity", type="choice", choices=('0', '1', '2'),
       default='0', dest="minimum_trypticity",
       help="drop peptides with too few tryptic ends"
       " (2=fully tryptic, 1=semi-tryptic, 0=any)", metavar="TRYPTICITY")
    DEFAULT_MINIMUM_PEPTIDES_PER_LOCUS = 1
    pa("-p", "--minimum-peptides-per-locus", dest="minimum_peptides_per_locus",
       type="int", default=DEFAULT_MINIMUM_PEPTIDES_PER_LOCUS,
       help="only loci with this level of peptide coverage are included, and"
       " peptides not associated with such a locus are excluded"
       " [default=%s]" % DEFAULT_MINIMUM_PEPTIDES_PER_LOCUS,
       metavar="COUNT")
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
    pa("-k", "--kill", action="store_true", dest="kill",
       help="rewrite the input files, removing spectra that don't pass"
       " validation, according to filtering")
    pa("--copyright", action="store_true", dest="copyright",
       help="print copyright and exit")
    (options, args) = parser.parse_args(args=args)

    if options.copyright:
        print __copyright__
        sys.exit(0)

    if len(args) < 1 or options.minimum_peptides_per_locus < 0:
        parser.print_help()
        sys.exit(1)

    if not (0.0 <= options.fdr < 0.5):
        error("--fdr must be within range [0.0, 0.5)")
    if (sum(1 for x in (options.mark, options.reset_marks, options.kill) if x)
        > 1):
        error("only one of --mark, --reset-marks, --kill may be specified")

    if options.reset_marks:
        write_spectra_to_files(reset_marks(generate_spectra_from_files(args)))
        return

    # This is to translate the "effective" FDR, which is the rate after decoys
    # have been stripped out of the results (which is what users really care
    # about), into our internal FDR, which includes decoys.
    options.fdr *= 2

    z_scores, spectrum_scores = \
              get_spectrum_info(options.decoy_prefix,
                                int(options.minimum_trypticity),
                                generate_spectra_from_files(args))

    if options.debug:
        print ('%s passing spectra, of which %s are not decoys'
               % (sum(1 for x in spectrum_scores if x),
                  sum(len(z_scores[x]) for x in z_scores)))

    thresholds = search_adjusting_fdr(options, z_scores)

    if options.debug:
        pprint(thresholds)

        for ss in spectrum_scores:
            if ss:
                charge, score, delta = ss
                print '#S', charge, score, delta

    if options.mark or options.kill:
        # Note that we are intentionally reading the input files again here,
        # to avoid having to keep them all in memory.
        write_spectra_to_files(set_marks(options, thresholds, spectrum_scores,
                                         generate_spectra_from_files(args),
                                         options.kill))


if __name__ == '__main__':
    main()
