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
import itertools
import logging; from logging import debug, info, warning
import optparse
import os
from pprint import pprint
import re
import string
import sys
import time

import greylag


def error(s, *args):
    "fatal error"
    # if we're unit testing, just throw an exception
    if __name__ != "__main__":
        raise Exception((s + " (fatal error)") % args)
    logging.error(s, *args)
    sys.exit(1)

# errors are fatal
greylag.chase_error = error


def generate_spectra_from_files(sqt_filenames):
    """Yield tuples (filename, H+_lines, S_line, [(M_line, mod_lines, L_lines),
    ...]) for the given SQT files.  H+_lines, for example, is a list
    containing H lines and possibly other kinds of lines up to the next S
    line.  (mod_lines are non-L lines following an M-line.)
    """
    for sqt_filename in sqt_filenames:
        with open(sqt_filename) as sqtf:
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


def write_spectra_to_files(spectra, filenames):
    """Given a generator of tuples of the form generated by
    'generate_spectra_from_files', write them to the implied files.  If any of
    filenames already exist, they are renamed with a '.bak' suffix.  Temp
    files are written first, then renamed at the end.  (The filenames
    parameter is necessary because spectra may not contain any entries for a
    given filename, thus necessitating its removal.)
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

    assert set(filemap.keys()) <= set(filenames)
    for fn in filenames:
        if os.path.exists(fn):
            os.rename(fn, fn + '.bak')
        if os.path.exists(fn + tmp_suffix):
            os.rename(fn + tmp_suffix, fn)


def reset_marks(spectra):
    """Given a generator of tuples of the form generated by
    'generate_spectra_from_files', yield them with the M line marks
    (destructively) reset to U.

    >>> list(reset_marks([('fn', [], 'S...', [('M... Y', [], [])])]))
    [('fn', [], 'S...', [('M... U', [], [])])]

    """

    def reset_mark(m_line):
        return re.sub(r'(\s+)\S(\s*)$', r'\1U\2', m_line, 1)

    for sqt_filename, H_lines, S_line, ML_lines in spectra:
        ML_lines_1 = [ (reset_mark(M_line), mod_lines, L_lines)
                       for M_line, mod_lines, L_lines in ML_lines ]
        yield (sqt_filename, H_lines, S_line, ML_lines_1)


def set_marks(valid_spectrum_info, spectra, kill=False):
    """Given a generator of tuples 'spectra' of the form generated by
    'generate_spectra_from_files', yield them with the M line marks
    (destructively) set to 'N', for spectra not listed in
    valid_spectrum_info.  If 'kill', omit marked spectra altogether.

    >>> list(set_marks([(0,)], [('fn', [], 'S 000...', [('M... U', [], [])]),
    ...                         ('fn', [], 'S 111...', [('M... U', [], [])])]))
    [('fn', [], 'S 000...', [('M... U', [], [])]), ('fn', [], 'S 111...', [('M... N', [], [])])]
    >>> list(set_marks([(0,)], [('fn', [], 'S 000...', [('M... U', [], [])]),
    ...                         ('fn', [], 'S 111...', [('M... U', [], [])])],
    ...                kill=True))
    [('fn', [], 'S 000...', [('M... U', [], [])])]

    """

    def set_mark(m_line):
        return re.sub(r'\S(\s*)$', r'N\1', m_line, 1)

    passing_spectrum_numbers = frozenset(si[0] for si in valid_spectrum_info)

    for spectrum_no, (sqt_filename, H_lines, S_line, ML_lines) \
            in enumerate(spectra):
        mark_spectrum = spectrum_no not in passing_spectrum_numbers
        if mark_spectrum and not kill:
            ML_lines = [ (set_mark(M_line), mod_lines, L_lines)
                         for M_line, mod_lines, L_lines in ML_lines ]
        if not (kill and mark_spectrum):
            yield (sqt_filename, H_lines, S_line, ML_lines)


nulltrans = string.maketrans('','')
non_aa_chars = string.digits + string.punctuation

def strip_mods(s):
    """Given a peptide, return just the unmodified peptide seen.

    >>> strip_mods('G*REY@LAG')
    'GREYLAG'
    >>> strip_mods('GREYLAG')
    'GREYLAG'

    """
    # assuming no control chars present
    return s.translate(nulltrans, non_aa_chars)


def get_spectrum_info(options, spectra):
    """Given a generator of tuples 'spectra' of the form generated by
    'generate_spectra_from_files', yield an info tuple for each spectrum that
    passes the criteria specified in options (e.g., options.minimum_trypticity
    and options.maximum_sp_rank).
    """

    for spectrum_no, (sqt_filename,
                      H_lines, S_line, ML_lines) in enumerate(spectra):
        S_fields = S_line.split('\t')
        current_charge = int(S_fields[3])
        actual_mass = float(S_fields[6])
        spectrum_name = '.'.join([os.path.splitext(sqt_filename)[0]]
                                 + S_fields[1:4])

        current_score = None
        current_delta = None
        current_peptide_trypticity = None
        current_sp_rank = None
        current_peptide_mass = None
        current_state = set()
        current_loci = set()
        highest_rank_seen = 0
        best_peptide_seen = None
        current_peptide = None

        for M_line, mod_lines, L_lines in ML_lines:
            M_fields = M_line.split('\t')

            rank = int(M_fields[1])
            sp_rank = int(M_fields[2])
            peptide_mass = float(M_fields[3])
            score = float(M_fields[5])
            peptide = M_fields[9].strip().upper() # e.g., A.B@CD*.-
            assert score >= 0 and rank >= 1 and sp_rank >= 1

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
                current_peptide_mass = peptide_mass
                current_peptide_trypticity = 0
                # NB: K*, R* not tryptic!
                if (peptide_flanks[0] in ('K', 'R') and peptide[0] != 'P'
                    or peptide_flanks[0] == '-'):
                    current_peptide_trypticity += 1
                if (peptide[-1] in ('K', 'R') and peptide_flanks[1] != 'P'
                    or peptide_flanks[1] == '-'):
                    current_peptide_trypticity += 1
                current_sp_rank = sp_rank

            # Choose delta from first M line such that peptide differs from
            # peptide in initial (top) M line.  We consider differently
            # modified versions of the same peptide to be different (at least
            # for now).  But, if there are two peptides that differ only in
            # their choice of I (isoleucine) or L (leucine), consider them
            # identical and keep looking for a delta.
            if peptide.replace('I', 'L') != current_peptide.replace('I', 'L'):
                current_delta = delta
                break

            for L_line in L_lines:
                locus = L_line.split('\t')[1]
                current_loci.add(locus)
                if locus.startswith(options.decoy_prefix):
                    current_state.add('decoy')
                else:
                    current_state.add('real')

        if (current_score != None and current_delta != None
            and current_peptide_trypticity >= options.minimum_trypticity
            and current_sp_rank <= options.maximum_sp_rank):

            if len(current_state) == 2:
                current_state = 'both'
            else:
                assert len(current_state) == 1
                current_state = current_state.pop()

            yield (spectrum_no, spectrum_name, current_charge, current_score,
                   current_delta, current_state, current_peptide,
                   strip_mods(current_peptide), actual_mass,
                   actual_mass-current_peptide_mass, frozenset(current_loci))


# FIX: can this be more succinct?
def remove_charge_aliases(spectrum_generator):
    """Given a generator of spectra like get_spectrum_info, for spectra that
    received multiple id's at different charges, yield only the best id.
    (Probably unnecessary unless the precursor tolerance is wide.)
    """

    sp_info_0 = None
    for sp_info_1 in spectrum_generator:
        assert sp_info_1[2] < 10, "single-digit charge"
        if not sp_info_0:
            sp_info_0 = sp_info_1
            continue
        if sp_info_0[1][:-1] != sp_info_1[1][:-1]:
            yield sp_info_0
            sp_info_0 = sp_info_1
            continue
        # if at least one has sufficiently small mass delta, choose the
        # one with the smaller
        if abs(sp_info_0[9]) < abs(sp_info_1[9]) < 10:
            continue
        if abs(sp_info_1[9]) < abs(sp_info_0[9]) < 10:
            sp_info_0 = sp_info_1
            continue
        # otherwise choose the better-scoring one
        if sp_info_0[3] < sp_info_1[3]:
            sp_info_0 = sp_info_1

    if sp_info_0:
        yield sp_info_0


# General algorithm:
# - repeat, adjusting fdr upwards (or downwards?) to meet user goal
#     repeat until set of valid spectra converges:
#         - do pr filter of spectra
#             - if not first time, and no spectra filtered out, break
#         - separately for each charge, choose valid spectra to maximize real ids subject to FDR
#     - saturate (locally or completely)


def PPV(reals, decoys):
    """Returns the estimated Positive Predictive Value (== 1 - FDR), given
    counts of reals and decoys.

    >>> PPV(100, 100)
    0.0
    >>> PPV(99, 100)
    0.0
    >>> PPV(100, 0)
    1.0
    >>> str(PPV(95, 5))
    '0.9'

    """

    # We know that the decoys are false positives, and we estimate that an
    # equal number of the "reals" are actually false, too.
    false_positives = 2*decoys
    true_positives = reals - decoys

    if true_positives <= 0:
        return 0.0
    return float(true_positives) / (true_positives + false_positives)


def calculate_inner_threshold(fdr, spinfo):
    ppv_goal = 1 - fdr

    spinfo = sorted(spinfo, key=lambda x: x[3])

    epsilon = +1e-6

    real_count = sum(1 for x in spinfo if x[5] == 'real')
    decoy_count = sum(1 for x in spinfo if x[5] == 'decoy')

    current_threshold = -1e100      # allow all spectra
    for n, sp in enumerate(spinfo):
        if real_count == 0:
            return (None, real_count, decoy_count) # give up

        ppv_est = PPV(real_count, decoy_count)
        if ppv_est >= ppv_goal:
            break
        if sp[5] == 'real':
            real_count -= 1
        elif sp[5] == 'decoy':
            decoy_count -= 1
        # set threshold just high enough to exclude this spectrum
        current_threshold = sp[3] + epsilon
    else:
        current_threshold = spinfo[-1][3] + epsilon # couldn't meet goal

    return (current_threshold, real_count, decoy_count)


def calculate_combined_thresholds(fdr, spectrum_info):
    """Find best score/delta thresholds for each charge."""

    ppv_goal = 1 - fdr

    # Rather than search every possible value of delta, we're only going to
    # "sample" at this granularity.  This cuts search time dramatically (and
    # makes it O(n) instead of O(n**2).  Extra precision wouldn't really be
    # useful in any case.
    SEARCH_GRANULARITY = 0.001

    # charge -> (score, delta, passing_reals, passing_decoys)
    thresholds = {}

    by_charge = lambda x: x[2]
    spectrum_info.sort(key=by_charge)
    for charge, spinfo in itertools.groupby(spectrum_info, key=by_charge):
        spinfo0 = sorted(spinfo, key=lambda x: x[4], reverse=True)

        last_value = None

        while spinfo0:
            this_value = spinfo0[-1][4] # current delta
            if (last_value == None
                or abs(this_value - last_value) >= SEARCH_GRANULARITY):

                # "inner" is score
                r = calculate_inner_threshold(fdr, spinfo0)
                threshold, real_count, decoy_count = r
                if threshold != None:
                    debug('# %s %s %s %s %s', charge, threshold, this_value,
                          real_count, decoy_count)
                    if (charge not in thresholds
                        or real_count > thresholds[charge][2]):
                        thresholds[charge] = (threshold, this_value,
                                              real_count, decoy_count)
                last_value = this_value
            spinfo0.pop()

    return thresholds


def flatten(gen_of_gens):
    """Yield items from a generator of generators.

    >>> list(flatten([(1,2,3), (4, (5, 'b'), 6)]))
    [1, 2, 3, 4, (5, 'b'), 6]

    """

    for gen in gen_of_gens:
        for item in gen:
            yield item


def redundancy_filter(options, spectrum_info):
    if options.minimum_spectra_per_locus < 2:
        return [ si for si in spectrum_info ]

    # si = (spectrum_no, spectrum_name, charge, score, delta, state,
    #       peptide, stripped_peptide, actual_mass, mass_delta, loci)

    Rs = sorted(flatten(((locus, si) for locus in si[10])
                        for si in spectrum_info))

    result = set()
    for locus, locus_Rs in itertools.groupby(Rs, key=lambda x: x[0]):
        list_locus_si = [ si for locus, si in locus_Rs ]
        if len(list_locus_si) < options.minimum_spectra_per_locus:
            continue
        result |= set(list_locus_si)
    return sorted(result)


def saturate_by_prefix_suffix(spectrum_info, saturation_spectrum_info):
    PREFIX_LENGTH = 16          # empirically seems good
    peptide_prefixes = set(si[7][:PREFIX_LENGTH] for si in spectrum_info)
    sat_si = sorted(ssi for ssi in saturation_spectrum_info
                    if ssi[7][:PREFIX_LENGTH] in peptide_prefixes)

    info("saturation by prefix (w/%s): %s -> %s",
         len(saturation_spectrum_info), len(spectrum_info), len(sat_si))
    return sat_si


def saturate(spectrum_info, saturation_spectrum_info):
    peptides = set(si[7] for si in spectrum_info)
    # FIX: does this actually need to be sorted?
    sat_si = sorted(ssi for ssi in saturation_spectrum_info
                    if ssi[7] in peptides)
    info("saturation (w/%s): %s -> %s", len(saturation_spectrum_info),
         len(spectrum_info), len(sat_si))
    return saturate_by_prefix_suffix(sat_si, saturation_spectrum_info)


def fdr_stats(spectrum_info):
    reals = sum(1 for si in spectrum_info if si[5] == 'real')
    decoys = sum(1 for si in spectrum_info if si[5] == 'decoy')
    fdr = 1 - PPV(reals, decoys)
    return fdr, reals, decoys

def debug_fdr(options, msg, si):
    if options.verbose:
        fdr, reals, decoys = fdr_stats(si)
        info('%s %s %s %s', msg, fdr, reals, decoys)

def try_fdr(fdr_guess, options, remaining_spectrum_info,
            saturation_spectrum_info=None):
    info("trying FDR = %s", fdr_guess)
    debug_fdr(options, "initial", remaining_spectrum_info)

    remaining_spectrum_info_1 = None
    while True:
        thresholds = calculate_combined_thresholds(fdr_guess,
                                                   remaining_spectrum_info)
        remaining_spectrum_info = [ si for si in remaining_spectrum_info
                                    if (si[2] in thresholds
                                        and si[3] >= thresholds[si[2]][0]
                                        and si[4] >= thresholds[si[2]][1]) ]
        debug_fdr(options, "after thresholding", remaining_spectrum_info)
        if remaining_spectrum_info == remaining_spectrum_info_1:
            break
        remaining_spectrum_info_1 = redundancy_filter(options,
                                                      remaining_spectrum_info)
        debug_fdr(options, "after redundancy", remaining_spectrum_info_1)
        if len(remaining_spectrum_info_1) == len(remaining_spectrum_info):
            break
        remaining_spectrum_info = remaining_spectrum_info_1

    fdr_result, total_reals, total_decoys = fdr_stats(remaining_spectrum_info_1)
    if not options.no_saturation and saturation_spectrum_info:
        remaining_spectrum_info_1 = saturate(remaining_spectrum_info_1,
                                             saturation_spectrum_info)
        debug_fdr(options, "after saturation", remaining_spectrum_info_1)
        (fdr_result,
         total_reals, total_decoys) = fdr_stats(remaining_spectrum_info_1)

    return (fdr_result, thresholds, total_reals, total_decoys,
            remaining_spectrum_info_1)


# FIX: at various points here we may have a goal of a minimum number of real
# ids, based on goal FDR or whatever.  We can pass this inward to do
# branch-and-bound, which may save time.

# FIX: this is now somewhat slower than before--any way to speed it up?

def search_adjusting_fdr(options, spectrum_info_0):
    """Repeatedly calculate thresholds while adjusting FDR, looking for an FDR
    guess that results in an appropriate goal FDR.  This might be needed
    because the thresholds found for per-charge FDRs will be saturated or have
    quantization errors, resulting in overly conservative thresholds.  Returns
    a list of items from spectrum_info that are considered valid.
    """

    FDR_INFLATION_FACTOR = 3
    assert FDR_INFLATION_FACTOR > 1

    spectrum_info = redundancy_filter(options, spectrum_info_0)

    low_fdr = options.fdr
    # FIX: spectrum_info_0 has already been filtered for tryptic status
    # Should we use an unfiltered set instead?
    low_results = try_fdr(low_fdr, options, spectrum_info, spectrum_info_0)

    if not options.no_adjust_fdr:
        high_fdr = min(1.0, options.fdr * FDR_INFLATION_FACTOR)
        high_results = try_fdr(high_fdr, options, spectrum_info,
                               spectrum_info_0)
        initial_total_reals = low_results[2]

        for i in range(32):
            debug("adjust infl fdr %s %s -> %s", low_fdr, high_fdr,
                 high_results[0])
            if high_results[0] >= options.fdr:
                break
            low_fdr, low_results = high_fdr, high_results
            high_fdr = min(1.0, high_fdr * FDR_INFLATION_FACTOR)
            if high_fdr == 1.0:
                break
            high_results = try_fdr(high_fdr, options, spectrum_info,
                                   spectrum_info_0)
        else:
            warning("FDR adjustment inflation failed")

        #    #
        #    l        g            h
        #   l       g         h
        #    #

        CONVERGENCE_FACTOR = 0.01
        assert CONVERGENCE_FACTOR > 0
        for i in range(32):
            assert high_fdr >= low_fdr
            if (high_fdr - low_fdr) <= CONVERGENCE_FACTOR*options.fdr:
                break

            # - seems to be slower than the naive method?
            # Bias is empirically determined.  It's supposed to speed
            # convergence, but should not affect the results.
            # FIX: better method possible?
            #GUESS_BIAS = 1.15
            #prop = ((options.fdr - low_results[0])
            #        / (high_results[0] - low_results[0]))
            #prop = min(1, max(0, prop))
            #guess_fdr = low_fdr + prop**GUESS_BIAS * (high_fdr - low_fdr)

            guess_fdr = (high_fdr + low_fdr) / 2.0

            guess_results = try_fdr(guess_fdr, options, spectrum_info,
                                    spectrum_info_0)
            debug("adjust fdr %s %s %s -> %s", low_fdr, guess_fdr, high_fdr,
                 guess_results[0])
            if guess_results[0] > options.fdr:
                high_fdr, high_results = guess_fdr, guess_results
            else:
                low_fdr, low_results = guess_fdr, guess_results
        else:
            warning("FDR adjustment convergence failed")

    if not options.no_adjust_fdr:
        info("FDR adjustment found %s extra real ids",
             low_results[2] - initial_total_reals)

    return low_results


def search_with_possible_saturation(options, spectrum_info):

    for saturate in (True, False):
        if saturate and options.no_saturation:
            continue
        results = search_adjusting_fdr(options, spectrum_info)
        (fdr_result, thresholds,
         total_reals, total_decoys, valid_spectrum_info) = results
        if fdr_result <= options.fdr:
            break
        if saturate:
            info('could not achieve FDR with saturation, retrying without')
            options.no_saturation = True

    if options.output_peptides:
        valid_spectrum_info.sort(key=lambda x: x[1])
        for (spectrum_no, spectrum_name, charge, score, delta, state,
             peptide, stripped_peptide, actual_mass,
             mass_delta, loci) in valid_spectrum_info:
            print >> options.output_peptides, \
                ("%+d %s %s %s %.8f %.8f"
                 % (charge, spectrum_name, state, peptide, actual_mass,
                    mass_delta))

    for charge in sorted(thresholds.keys()):
        score_threshold, delta_threshold = thresholds[charge][0:2]
        reals, decoys = thresholds[charge][2:4]
        info("%+d: score %f, delta %.6f -> %s real ids, %s decoys (fdr %.4f)",
             charge, score_threshold, delta_threshold, reals, decoys,
             0.5 * (1 - PPV(reals, decoys)))

    if not options.quiet:
        print ("%s real ids, %s decoys (FDR = %.4f) [p%s]"
               % (total_reals, total_decoys, fdr_result / 2.0,
                  options.minimum_spectra_per_locus))

    return valid_spectrum_info


def print_mass_error_histogram(valid_spectrum_info):
    hist = defaultdict(int)

    for (spectrum_no, spectrum_name, charge, score, delta, state, peptide,
         stripped_peptide, actual_mass, mass_delta,
         loci) in valid_spectrum_info:
        hist[round(mass_delta)] += 1

    pairs_most_first = sorted(((bin, hist[bin]) for bin in hist),
                              key=lambda x: x[1], reverse=True)

    total = sum(count for bin, count in pairs_most_first)

    cumulative_total = 0
    for bin, count in pairs_most_first:
        cumulative_total += count
        print '# %+d\t%s\t%.3f' % (bin, count,
                                   (float(cumulative_total) / total))


def main(args=sys.argv[1:]):
    parser = optparse.OptionParser(usage=
                                   "usage: %prog [options] <sqt-file>...",
                                   description=__doc__, version=greylag.VERSION)
    pa = parser.add_option
    pa("--decoy-prefix", dest="decoy_prefix",
       default=greylag.DEFAULT_DECOY_PREFIX,
       help='prefix given to locus name of decoy (e.g., shuffled) database'
       ' sequences [default=%r]' % greylag.DEFAULT_DECOY_PREFIX,
       metavar="PREFIX")
    DEFAULT_FDR = 0.01
    pa("--fdr", dest="fdr", type="float", default=DEFAULT_FDR,
       help="false discovery rate [default=%s] (Note that this is the"
       " final resulting FDR after the decoys have been removed" % DEFAULT_FDR,
       metavar="PROPORTION")
    pa("--no-adjust-fdr", action="store_true", dest="no_adjust_fdr",
       help="skip adjustment of internal FDR to achieve a final FDR closer to"
       " requested FDR (which takes longer, but may find more real ids in"
       " some cases)")
    pa("--no-saturation", action="store_true", dest="no_saturation",
       help="skip peptide saturation")
    pa("-v", "--verbose", action="store_true", dest="verbose",
       help="be verbose")
    pa("-q", "--quiet", action="store_true", dest="quiet", help="no warnings")
    pa("-l", "--logfile", dest="logfile",
       help="log to FILE instead of stderr", metavar="FILE")
    pa("--output-peptides", dest="output_peptides",
       help="output information about passing peptides to specified file"
       " ('-' for stdout)", metavar="FILENAME")
    pa("-t", "--minimum-trypticity", type="choice", choices=('0', '1', '2'),
       default='0', dest="minimum_trypticity",
       help="drop peptides with too few tryptic ends"
       " (2=fully tryptic, 1=semi-tryptic, 0=any)", metavar="TRYPTICITY")
    DEFAULT_MINIMUM_SPECTRA_PER_LOCUS = 1
    pa("-p", "--minimum-spectra-per-locus", dest="minimum_spectra_per_locus",
       type="int", default=DEFAULT_MINIMUM_SPECTRA_PER_LOCUS,
       help="only loci with at least this many spectra are included, and"
       " spectra not associated with such a locus are excluded"
       " [default=%s]" % DEFAULT_MINIMUM_SPECTRA_PER_LOCUS,
       metavar="COUNT")
    DEFAULT_MAXIMUM_SP_RANK = 1000000
    pa("--maximum-sp-rank", dest="maximum_sp_rank", type="int",
       default=DEFAULT_MAXIMUM_SP_RANK,
       help="drop peptides with greater Sp rank (the secondary rank in sqt"
       " files, for compatibility with other search programs)"
       " [default=%s]" % DEFAULT_MAXIMUM_SP_RANK, metavar="RANK")
    #pa("--graph", dest="graph",
    #   help='create distribution graphs, using the specified prefix',
    #   metavar="PATH PREFIX")
    pa("-M", "--mass-error-histogram", action="store_true",
       dest="mass_error_histogram",
       help="print a histogram of mass errors (1-Da-wide bins)")
    pa("--debug", action="store_true", dest="debug",
       help="show debug output")
    pa("--mark", action="store_true", dest="mark",
       help="rewrite the input files, changing some validation marks to 'N',"
       " according to filtering")
    pa("--reset-marks", action="store_true", dest="reset_marks",
       help="rewrite the input files, changing all validation marks to 'U'")
    pa("--kill", action="store_true", dest="kill",
       help="rewrite the input files, removing spectra that don't pass"
       " validation, according to filtering")
    pa("--copyright", action="store_true", dest="copyright",
       help="print copyright and exit")
    (options, args) = parser.parse_args(args=args)

    if options.copyright:
        print __copyright__
        sys.exit(0)

    if (len(args) < 1 or options.minimum_spectra_per_locus < 0
        or options.maximum_sp_rank < 1):
        parser.print_help()
        sys.exit(1)

    if not (0.0 <= options.fdr < 0.5):
        error("--fdr must be within range [0.0, 0.5)")
    if (sum(1 for x in (options.mark, options.reset_marks, options.kill) if x)
        > 1):
        error("only one of --mark, --reset-marks, --kill may be specified")

    options.minimum_trypticity = int(options.minimum_trypticity)

    greylag.set_logging(options)

    if options.reset_marks:
        write_spectra_to_files(reset_marks(generate_spectra_from_files(args)),
                               args)
        return

    if options.output_peptides:
        if options.output_peptides == '-':
            options.output_peptides = sys.stdout
        else:
            options.output_peptides = open(options.output_peptides, 'w')

    # This is to translate the "effective" FDR, which is the rate after decoys
    # have been stripped out of the results (which is what users really care
    # about), into our internal FDR, which includes decoys.
    options.fdr *= 2

    info("reading spectra")
    spectrum_info = list(remove_charge_aliases(get_spectrum_info(options,
                                                                 generate_spectra_from_files(args))))

    # valid_spectrum_info is the list of spectrum_info elements that have been
    # chosen as "valid"
    valid_spectrum_info = search_with_possible_saturation(options, spectrum_info)

    if options.mass_error_histogram:
        print_mass_error_histogram(valid_spectrum_info)

    if options.debug:
        valid_spectrum_info.sort()
        for (spectrum_no, spectrum_name, charge, score, delta, state,
             peptide, stripped_peptide, actual_mass,
             mass_delta, loci) in valid_spectrum_info:
            debug('#S %s %s %s', charge, score, delta)

    if options.mark or options.kill:
        # Note that we are intentionally reading the input files again here,
        # to avoid having to keep them all in memory.
        info("rewriting spectrum files")
        write_spectra_to_files(set_marks(valid_spectrum_info,
                                         generate_spectra_from_files(args),
                                         options.kill),
                               args)
        info("finished rewriting")


if __name__ == '__main__':
    main()
