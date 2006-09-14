#!/usr/bin/env python2.4

'''
Analyze mass spectra and assign peptides and proteins.  (Partial
re-implementation of X!Tandem, with extra features.)

'''

__version__ = "$Id$"


import fileinput
import logging
from logging import debug, info, warning
import math
import optparse
import os
from pprint import pprint
import re
import sys

# FIX: this goes away in Python 2.5
try:
    import cElementTree as ElementTree
    print >> sys.stderr, 'imported celementtree'
except ImportError:
    import elementtree.ElementTree as ElementTree

import cxtpy

import xtpy_search


def error(s, *args):
    logging.error(s, *args)
    sys.exit(1)
def fileerror(s, *args):
    error(s + (", at line %s of file '%s'"
               % (fileinput.filelineno(), fileinput.filename())),
          *args)


# name -> value map of processed XML parameters
XTP = {}

# handle to the singleton parameter object shared with the C++ module
CP = cxtpy.cvar.parameters_the

PROTON_MASS = 1.007276

# atom -> (monoisotopic mass, average mass)
ATOMIC_MASS = {
    'H' : (  1.007825035,  1.00794  ),
    'O' : ( 15.99491463,  15.9994   ),
    'N' : ( 14.003074,    14.0067   ),
    'C' : ( 12.0,         12.0107   ),
    'S' : ( 31.9720707,   32.065    ),
    'P' : ( 30.973762,    30.973761 ),
    }

def formula_mass(formula, monoisotopic=True):
    """Return the mass of formula."""
    parts = [ p or '1' for p in re.split(r'([A-Z][a-z]*)', formula)[1:] ]
    # parts for alanine is ('C', '3', 'H', '5', 'O', '1', 'N', '1')
    return sum(ATOMIC_MASS[parts[i]][int(not monoisotopic)] * int(parts[i+1])
               for i in range(0, len(parts), 2))

# residue -> formula
RESIDUE_FORMULA = {
    'A' : "C3H5ON",
    'C' : "C3H5ONS",
    'D' : "C4H5O3N",
    'E' : "C5H7O3N",
    'F' : "C9H9ON",
    'G' : "C2H3ON",
    'H' : "C6H7ON3",
    'I' : "C6H11ON",
    'K' : "C6H12ON2",
    'L' : "C6H11ON",
    'M' : "C5H9ONS",
    'N' : "C4H6O2N2",
    'P' : "C5H7ON",
    'Q' : "C5H8O2N2",
    'R' : "C6H12ON4",
    'S' : "C3H5O2N",
    'T' : "C4H7O2N",
    'V' : "C5H9ON",
    'W' : "C11H10ON2",
    'Y' : "C9H9O2N",
    }

# FIX: selenocysteine (U), etc
# residue -> (monoisotopic mass, average mass)
RESIDUE_MASS = {
    'A' : [ 0.0,   71.0788 ],
    'C' : [ 0.0 , 103.1388 ],
    'D' : [ 0.0 , 115.0886 ],
    'E' : [ 0.0 , 129.1155 ],
    'F' : [ 0.0 , 147.1766 ],
    'G' : [ 0.0 ,  57.0519 ],
    'H' : [ 0.0 , 137.1411 ],
    'I' : [ 0.0 , 113.1594 ],
    'K' : [ 0.0 , 128.1741 ],
    'L' : [ 0.0 , 113.1594 ],
    'M' : [ 0.0 , 131.1926 ],
    'N' : [ 0.0 , 114.1038 ],
    'P' : [ 0.0 ,  97.1167 ],
    'Q' : [ 0.0 , 128.1307 ],
    'R' : [ 0.0 , 156.1875 ],
    'S' : [ 0.0 ,  87.0782 ],
    'T' : [ 0.0 , 101.1051 ],
    'V' : [ 0.0 ,  99.1326 ],
    'W' : [ 0.0 , 186.2132 ], 
    'Y' : [ 0.0 , 163.1760 ],
    }

for residue in RESIDUE_MASS:
    RESIDUE_MASS[residue][0] = formula_mass(RESIDUE_FORMULA[residue])
    

def initialize_spectrum_parameters(quirks_mode):
    """Initialize parameters known to the spectrum module."""

    CP.monoisotopic_atomic_mass.resize(128, 0.0)
    CP.average_atomic_mass.resize(128, 0.0)
    for atom, (monomass, avgmass) in ATOMIC_MASS.iteritems():
        CP.monoisotopic_atomic_mass[ord(atom)] = monomass
        CP.average_atomic_mass[ord(atom)] = avgmass

    CP.monoisotopic_residue_mass.resize(128, 0.0)
    CP.average_residue_mass.resize(128, 0.0)
    for residue, (monomass, avgmass) in RESIDUE_MASS.iteritems():
        CP.monoisotopic_residue_mass[ord(residue)] = monomass
        CP.average_residue_mass[ord(residue)] = avgmass

    CP.modification_mass.resize(128)
    for residue, modvalue in XTP["residue, modification mass"]:
        CP.modification_mass[ord(residue)] = modvalue

    CP.potential_modification_mass.resize(128)
    for residue, modvalue in XTP["residue, potential modification mass"]:
        # NB: SWIG currently only exposes inner vector as tuple
        # CP.potential_modification_mass[ord(residue)].append(modvalue)
        CP.potential_modification_mass[ord(residue)] \
            = CP.potential_modification_mass[ord(residue)] + (modvalue,)
    CP.potential_modification_mass_refine.resize(128)
    for residue, modvalue in XTP["refine, potential modification mass"]:
        # CP.potential_modification_mass_refine[ord(residue)].append(modvalue)
        CP.potential_modification_mass_refine[ord(residue)] \
            = CP.potential_modification_mass_refine[ord(residue)] + (modvalue,)

    CP.cleavage_N_terminal_mass_change \
        = XTP["protein, cleavage N-terminal mass change"]
    CP.cleavage_C_terminal_mass_change \
        = XTP["protein, cleavage C-terminal mass change"]

    CP.parent_monoisotopic_mass_error_plus \
        = XTP["spectrum, parent monoisotopic mass error plus"]
    CP.parent_monoisotopic_mass_error_minus \
        = XTP["spectrum, parent monoisotopic mass error minus"]
    CP.fragment_mass_error = XTP["spectrum, fragment mass error"]

    CP.minimum_ion_count = XTP["scoring, minimum ion count"]

    CP.proton_mass = PROTON_MASS
    CP.hydrogen = formula_mass("H")
    CP.hydroxide = formula_mass("OH")
    CP.water_mass = formula_mass("H2O")

    # CP.factorial[n] == (double) n!
    CP.factorial.resize(100, 1.0)
    for n in range(2, len(CP.factorial)):
        CP.factorial[n] = CP.factorial[n-1] * n

    CP.quirks_mode = bool(quirks_mode)
    

def cleavage_motif_re(motif):
    """Return (regexp, pos), where regexp is a regular expression that will
    match a cleavage motif, and pos is the position of the cleavage with
    respect to the match (e.g., 1 for '[KR]|{P}', or None if absent).
    """
    # FIX: use look-ahead groups to keep actual length == 1
    assert False, "before using, think about overlapping matches"
    cleavage_pos = None
    re_parts = []
    motif_re = re.compile(r'(\||(\[(X|[A-WYZ]+)\])|({([A-WYZ]+)}))')
    parts = [ p[0] for p in motif_re.findall(motif) ]
    if ''.join(parts) != motif:
        raise ValueError('invalid cleavage motif pattern')
    i = 0
    for part in parts:
        if part == '|':
            if cleavage_pos != None:
                raise ValueError("invalid cleavage motif pattern (multiple '|'s)")
            cleavage_pos = i
            continue
        if part == '[X]':
            re_parts.append('.')
        elif part[0] == '[':
            re_parts.append('[%s]' % part[1:-1])
        elif part[0] == '{':
            re_parts.append('[^%s]' % part[1:-1])
        else:
            assert False, "unknown cleavage motif syntax"
        i += 1
    return (''.join(re_parts), cleavage_pos)


def get_prefix_sequence(begin_pos, run_offset, sequence):
    prefix_start = run_offset + begin_pos - 4
    s = sequence[max(0, prefix_start):prefix_start+4]
    if len(s) < 4:
        s = '[' + s
    return s

def get_suffix_sequence(end_pos, run_offset, sequence):
    suffix_start = run_offset + end_pos
    s = sequence[suffix_start:suffix_start+4]
    if len(s) < 4:
        s = s + ']'
    return s


def generate_peptides(seq, cleavage_points, maximum_missed_cleavage_sites):
    """Yield (begin, end, missed_cleavage_count) for each apt peptide in
    sequence.""" 
    len_cp = len(cleavage_points)
    for begin_i in range(len_cp-1):
        for end_i in range(begin_i+1,
                           min(len_cp,
                               begin_i+2+maximum_missed_cleavage_sites)):
            begin, end = cleavage_points[begin_i], cleavage_points[end_i]
            if end - begin >= 4:        # make 4 a param!
                yield begin, end, end_i-begin_i-1


def generate_cleavage_points(sequence):
    """Yields the offsets of the cleavages in sequence."""
    yield 0
    # if pattern length > 1, will miss overlapping matches!
    for m in re.finditer(r'[KR](?=[^P])', sequence):
        yield m.start() + 1
    yield len(sequence)

aa_sequence = re.compile(r'[ARNDCQEGHILKMFPSTWYV]+')

def split_sequence_into_aa_runs(idno, defline, sequence, filename):
    """Returns a tuple (idno, start, defline, seq, is_N, is_C) for each
    contiguous run of residues in sequence, where 'start' is the position of
    'seq' in 'sequence', and 'is_N' and 'is_C' are flags indicating whether
    the the beginning/end of 'seq' is at the beginning/end of 'sequence'."""
    matches = list(aa_sequence.finditer(sequence))
    return [ (idno, m.start(), defline, m.group(),
              n==0 and m.start()==0,
              n==len(matches)-1 and m.end()==len(sequence), filename)
             for n, m in enumerate(matches) ]


def clean_defline(s):
    """Return the given string with control characters removed."""
    return re.sub(r'[^ -~]', '', s)


def abbrev_defline(s):
    """Return an abbreviated version of the defline--about 80 characters."""
    ab = re.match(r'.{,80}\S{,170}', s).group(0)
    if len(ab) < len(s):
        ab += '...'
    return ab
    

def read_fasta_files(filenames):
    """Yield (defline, sequence, filename) tuples as read from FASTA files
    (uppercasing sequence)."""
    defline = None
    seqs = []
    for line in fileinput.input(filenames):
        line = line.strip()
        if line[:1] == '>':
            if defline:
                yield (defline, ''.join(seqs), fileinput.filename())
            elif seqs:
                fileerror("bad format: line precedes initial defline")
            defline = line[1:]
            seqs = []
        else:
            seqs.append(line.upper())
    if defline:
        yield (defline, ''.join(seqs), fileinput.filename())


def read_taxonomy(filename):
    """Return a map of taxa to lists of filenames."""
    root = ElementTree.ElementTree(file=filename).getroot()
    return dict( (taxon.get("label"),
                  [ f.get("URL") for f in taxon.findall('./file') ])
                 for taxon in root.findall('taxon') )


def read_xml_parameters(filename):
    """Return a map of parameters to values, per parameter file fn."""
    root = ElementTree.ElementTree(file=filename).getroot()
    return dict((e.get("label"), e.text)
                for e in root.findall('note')
                if e.get("type") == "input")


# XML parameter file processing

def mod_list(modification_list_specification, unique_mods=True):
    """Check and return a sorted list of (residue, modification) tuples.  If
    unique_mods is True, also check that each residue is specified only once.
    """
    if not modification_list_specification:
        return []
    speclist = []
    for modspec in modification_list_specification.split(','):
        parts = modspec.split('@')
        if len(parts) != 2:
            raise ValueError("invalid modification list specification"
                             " (missing or extra '@'?)")
        value = float(parts[0])
        residue = parts[1]
        if len(residue) != 1:
            raise ValueError("invalid modification list specification"
                             " (single-letter residue expected, got '%s')"
                             % residue)
        if residue not in RESIDUE_FORMULA and residue not in '[]':
            raise ValueError("invalid modification list specification"
                             " (invalid residue '%s')" % residue)
        speclist.append((residue, value))
    if unique_mods:
        if len(set(residue for residue, value in speclist)) != len(speclist):
            raise ValueError("invalid modification list specification"
                             " (residue specified multiple times)")
    return sorted(speclist)

def potential_mod_list(modification_list_specification):
    """Check and return a list of (residue, modification) tuples."""
    return mod_list(modification_list_specification, unique_mods=False)
        

# "ni" means check verifies that "not implemented" functionality is not
# specified
def p_ni_equal(x): return lambda y: y == x
def p_ni_empty(x): return len(x) == 0
def p_positive(x): return x > 0
def p_negative(x): return x < 0
def p_nonnegative(x): return x >= 0

# name -> (type, default, check_fn)
# type may be bool, int, float, str, modlist, or a tuple of values
# default, as a string value, or None if value must be explicitly specified
# check_fn is an optional function that returns True iff the value is valid
# *1: param required, as has differing defaults based on other params, in
# xtandem

XML_PARAMETER_INFO = {
    "list path, default parameters" : (str, ""),
    "list path, taxonomy information" : (str, None),
    "output, histogram column width" : (int, 30), # currently ignored
    "output, histograms" : (bool, "no"),
    "output, log path" : (str, ""),     # ignored
    "output, maximum valid expectation value" : (float, None),
    "output, message" : (str, "."),     # ignored
    "output, one sequence copy" : (bool, "no", p_ni_equal(False)),
    "output, parameters" : (bool, "no"),
    "output, path hashing" : (bool, "no", p_ni_equal(False)),
    "output, path" : (str, ""),
    "output, performance" : (bool, "no"),
    "output, proteins" : (bool, "no"),
    "output, results" : (('all', 'valid', 'stochastic'), "all", p_ni_equal("valid")),
    "output, sequence path" : (str, "", p_ni_empty),
    "output, sequences" : (bool, "no"),
    "output, sort results by" : (('protein', 'spectrum'), "spectrum", p_ni_equal("protein")),
    "output, spectra" : (bool, "no"),
    "output, xsl path" : (str, ""),
    "protein, C-terminal residue modification mass" : (float, "0.0", p_ni_equal(0)), # eliminate, obsolete
    "protein, N-terminal residue modification mass" : (float, "0.0", p_ni_equal(0)), # eliminate, obsolete
    "protein, cleavage C-terminal mass change" : (float, formula_mass("OH")),
    "protein, cleavage N-terminal limit" : (int, "100000000", p_positive),
    "protein, cleavage N-terminal mass change" : (float, formula_mass("H")),
    "protein, cleavage semi" : (bool, "no", p_ni_equal(False)),
    "protein, cleavage site" : (("[RK]|{P}",), None, p_ni_equal("[RK]|{P}")),
    "protein, homolog management" : (bool, "no", p_ni_equal(False)),
    "protein, modified residue mass file" : (str, "", p_ni_empty),
    "protein, taxon" : (str, None),
    "refine" : (bool, "no", p_ni_equal(False)),
    "refine, cleavage semi" : (bool, "no", p_ni_equal(False)),
    "refine, maximum valid expectation value" : (float, None),
    "refine, modification mass" : (mod_list, "", p_ni_empty),
    "refine, point mutations" : (bool, "no", p_ni_equal(False)),
    "refine, potential C-terminus modifications": (potential_mod_list, "", p_ni_empty),
    "refine, potential N-terminus modifications": (potential_mod_list, "", p_ni_empty),
    "refine, potential N-terminus modification position limit" : (int, 50), # nyi
    "refine, potential modification mass" : (potential_mod_list, "", p_ni_empty),
    "refine, potential modification motif" : (str, "", p_ni_empty),
    "refine, sequence path" : (str, "", p_ni_empty),
    "refine, spectrum synthesis" : (bool, "no", p_ni_equal(False)),
    "refine, tic percent" : (float, 20.0, p_nonnegative), # ignored
    "refine, unanticipated cleavage" : (bool, p_ni_equal(False)),
    "refine, use potential modifications for full refinement" : (bool, "no", p_ni_equal(False)),
    "residue, modification mass" : (mod_list, ""),
    "residue, potential modification mass" : (potential_mod_list, "", p_ni_empty),
    "residue, potential modification motif" : (str, "", p_ni_empty),
    "scoring, a ions" : (bool, "no", p_ni_equal(False)),
    "scoring, b ions" : (bool, "yes", p_ni_equal(True)),
    "scoring, c ions" : (bool, "no", p_ni_equal(False)),
    "scoring, cyclic permutation" : (bool, "no", p_ni_equal(False)),
    "scoring, include reverse" : (bool, "no", p_ni_equal(False)),
    "scoring, maximum missed cleavage sites" : (int, None, p_nonnegative),
    "scoring, minimum ion count" : (int, None, p_positive),
    "scoring, pluggable scoring" : (bool, "no"), # ignored
    "scoring, x ions" : (bool, "no", p_ni_equal(False)),
    "scoring, y ions" : (bool, "yes", p_ni_equal(True)),
    "scoring, z ions" : (bool, "no", p_ni_equal(False)),
    "spectrum, check all charges" : (bool, "no", p_ni_equal(False)),
    "spectrum, dynamic range" : (float, "100", p_positive),
    "spectrum, fragment mass error units" : (("Daltons", "ppm"), "Daltons", p_ni_equal("Daltons")),
    "spectrum, fragment mass error" : (float, "0.45", p_positive),
    "spectrum, fragment mass type" : (("average", "monoisotopic"), "monoisotopic", p_ni_equal("monoisotopic")),
    "spectrum, homology error" : (float, "4.5", p_positive), # nyi
    "spectrum, maximum parent charge" : (int, "4", p_positive), # nyi
    "spectrum, minimum fragment mz" : (float, "200.0", p_positive),
    "spectrum, minimum parent m+h" : (float, "850.0", p_positive), # nyi
    "spectrum, minimum peaks" : (int, "5", p_positive),
    "spectrum, neutral loss mass" : (float, "0.0", p_nonnegative), # nyi
    "spectrum, neutral loss window" : (float, "0.0", p_nonnegative), # nyi
    "spectrum, parent monoisotopic mass error minus" : (float, None, p_negative), # *1
    "spectrum, parent monoisotopic mass error plus" : (float, None, p_positive), # *1
    "spectrum, parent monoisotopic mass error units" : (("Daltons", "ppm"), "Daltons", p_ni_equal("Daltons")),
    "spectrum, parent monoisotopic mass isotope error" : (bool, "no", p_ni_equal(False)), # parent error should be <0.5Da
    "spectrum, path" : (str, None),
    "spectrum, sequence batch size" : (int, 1000), # ignored
    "spectrum, threads" : (int, 1),     # ignored
    "spectrum, total peaks" : (int, "50", p_positive),
    "spectrum, use conditioning" : (bool, "yes"),
    "spectrum, use contrast angle" : (bool, "no", p_ni_equal(False)),
    "spectrum, use neutral loss window" : (bool, "no", p_ni_equal(False)),
    "spectrum, use noise suppression" : (bool, "yes", p_ni_equal(False)),
}

def validate_parameters(parameters):
    """Verify that parameters are valid, have valid values, and correspond to
    currently implemented functionality.  Default values are filled in, and a
    name/value dict returned."""

    pmap = {}
    for p_name, p_info in sorted(XML_PARAMETER_INFO.items()):
        type_ = p_info[0]
        default = p_info[1]
        check_fn = len(p_info) > 2 and p_info[2] or None

        v = parameters.get(p_name)
        if isinstance(v, str):
            v = v.strip()
        if v == None:
            if default != None:
                info("parameter '%s' defaulting to '%s'", p_name, default)
                v = default
            else:
                error("missing required parameter '%s'" % p_name)
        if isinstance(type_, tuple):
            if not v in type_:
                error("parameter '%s' value '%s' not in %s (feature not implemented?)"
                      % (p_name, v, type_))
        elif type_ == bool:
            v = { 'yes' : True, 'no' : False }.get(v)
            if v == None:
                error("parameter '%s' requires a value of 'yes' or 'no'")
        else:
            try:
                v = type_(v)
            except ValueError, e:
                error("parameter '%s' has value '%s' with invalid format [%s]"
                      % (p_name, v, e))
        if check_fn and not check_fn(v):
            error("parameter '%s' has invalid value '%s' (or feature not implemented)"
                  % (p_name, v))
        pmap[p_name] = v

    unknown_parameters = set(parameters) - set(XML_PARAMETER_INFO)
    if unknown_parameters:
        warning("%s unknown parameters:\n %s"
             % (len(unknown_parameters), pformat(sorted(list(unknown_parameters)))))

    # FIX: this would trip the "cyclic" param (assume Daltons)
    assert (abs(pmap["spectrum, parent monoisotopic mass error plus"]) > 0.095
            and abs(pmap["spectrum, parent monoisotopic mass error minus"]) > 0.095), \
            "feature not implemented (cyclic param)"

    return pmap


def get_spectrum_expectation(hyper_score, hyperscore_histogram):
    """Return (expectation, survival curve, (a0, a1)) for this hyperscore and
    histogram, where a0 and a1 are the intercept and slope, respectively, for
    the least-squares fit to the survival curve, used to predict the
    expectation.
    """
    scaled_hyper_score = cxtpy.scale_hyperscore(hyper_score)
    max_histogram_index = max(hyperscore_histogram)
    assert max_histogram_index < 1000, "assume these aren't huge"
    # histogram as a list
    histogram = [0] * (max_histogram_index+1)
    for score, count in hyperscore_histogram.iteritems():
        histogram[score] = count

    counts = sum(histogram)
    c = counts
    # survival is just a reversed cumulative distribution
    survival = [0] * len(histogram)
    for i in range(len(histogram)):
        survival[i] = c
        c -= histogram[i]
    # note: neither histogram nor survival have 0 as rightmost element here

    # this next looks bogus, but use for now to replicate xtandem
    # "remove potentially valid scores from the stochastic distribution"
    # makes the survival function non-monotonic?
    lPos = survival[0] / 5
    for lMid in range(len(survival)):
        if survival[lMid] <= lPos:
            break
    a = len(survival) - 1
    assert survival[a] != 0
    lSum = 0
    while a > 0:
        if (survival[a] == survival[a-1] and survival[a] != survival[0]
            and a > lMid):
            lSum = survival[a]
            a -= 1
            while survival[a] == lSum:
                assert a >= 0
                survival[a] -= lSum
                a -= 1
        else:
            survival[a] -= lSum
            a -= 1
    survival[a] -= lSum
    # end bogus

    if counts < 200:
        # use default line
        a0, a1 = 3.5, -0.18
        return 10.0 ** (a0 + a1 * scaled_hyper_score), survival, (a0, a1)
    
    min_limit = 10
    max_limit = int(round(survival[0]/2.0))
    max_i = min(i for i in range(len(survival)) if survival[i] <= max_limit)
    min_i = min(i for i in range(max_i, len(survival)) if survival[i] <= min_limit)
    data_X = range(max_i, min_i)
    data_Y = [ math.log10(survival[x]) for x in data_X ]
    assert data_Y[0] == max(data_Y), "impossible xtandem case?"

    # fit least-squares line
    n = len(data_X)
    sum_X, sum_Y = sum(data_X), sum(data_Y)
    sum_XX = sum(x**2 for x in data_X)
    sum_XY = sum(x*y for x,y in zip(data_X, data_Y))

    m = (n*sum_XY - sum_X*sum_Y) / (n*sum_XX - sum_X**2)
    b = (sum_Y - m*sum_X) / n
    return 10.0 ** (b + m * scaled_hyper_score), survival, (b, m)


def get_final_protein_expect(sp_count, valid_spectra, match_ratio, raw_expect,
                             db_length, candidate_spectra):
    debug('pe: %s', (sp_count, valid_spectra, match_ratio, raw_expect,
                     db_length, candidate_spectra))
    # (From xtandem) Compensate for multiple peptides supporting a protein.
    # The expectation values for the peptides are combined with a simple
    # Bayesian model for the probability of having two peptides from the same
    # protein having the best score in different spectra.
    # CHECK THIS
    if sp_count == 1:
        if raw_expect < 0.0:
            return raw_expect
        else:
            return 1.0
    assert sp_count > 0
    r = raw_expect + math.log10(db_length)
    for a in range(sp_count):
        r += math.log10(float(valid_spectra - a)/(sp_count - a))
    r -= math.log10(valid_spectra) + (sp_count-1) * math.log10(match_ratio)
    debug('x: %s', (match_ratio, candidate_spectra))
    p = min(float(match_ratio) / candidate_spectra, 0.9999999)
    debug("p: %s", p)
    r += (sp_count * math.log10(p)
          + (valid_spectra - sp_count) * math.log10(1 - p))
    return r


def process_results(best_match, hyperscore_histogram, fasta_db,
                    candidate_spectrum_count, spectra, db_residue_count):
    
    # spectrum index -> [ (protein id, [domain info, ...]), ... ]
    # domain info == (...)
    spec_prot_info = {}

    for sp_n, match_info_list in best_match.iteritems():
        prot_info = {}
        for m in match_info_list:
            protein_id = m[11]
            prot_info.setdefault(protein_id, []).append(m)
        spec_prot_info[sp_n] = prot_info.items()

    #debug("spec_prot_info: %s" % spec_prot_info)

    # calculate exp for best match for each spectrum

    # spectrum index -> expectation value
    expect = {}
    # spectrum index -> (scaled, binned hyperscore -> count)
    best_histogram = {}
    # spectrum index -> list of values
    survival_curve = {}
    # spectrum index -> (a0, a1)  [regression line intercept/slope]
    line_parameters = {}
    for sp_n, match_info_list in best_match.iteritems():
        # all hyperscores the same, so choose the first
        # all charges the same (by assumption above)
        sp_hyper_score = match_info_list[0][0]
        hh = hyperscore_histogram.get(sp_n, {})
        best_histogram[sp_n] = hh
        expect[sp_n], survival_curve[sp_n], line_parameters[sp_n] \
                      = get_spectrum_expectation(sp_hyper_score, hh)

    #debug("expect: %s" % expect)

    # these are the ids of spectra considered "good enough"
    passing_spectra = set(sp_n for sp_n, e in expect.iteritems()
                          if e <= XTP["refine, maximum valid expectation value"])

    #debug("passing_spectra: %s" % sorted(list(passing_spectra)))

    # protein id -> list of spectra match info
    best_protein_matches = {}
    for sp_n, match_info_list in best_match.iteritems():
        for match in match_info_list:
            best_protein_matches.setdefault(match[11], []).append(match)

    #debug("best_protein_matches: %s" % best_protein_matches)

    # calculate exp for each protein
    valid_spectra_count = sum(1 for sp_n in passing_spectra
                              if (expect[sp_n]
                                  <= XTP["output, maximum valid expectation value"]))
    best_histogram_sum = sum(sum(h.itervalues())
                             for sp_n, h in best_histogram.iteritems()                             if sp_n in passing_spectra)
    match_ratio = 0
    debug("best_histogram_sum: %s passing_spectra: %s", best_histogram_sum,
          passing_spectra)
    if passing_spectra:
        match_ratio = float(best_histogram_sum) / len(passing_spectra)

    # drop all but the best match for each physical spectrum
    # physical id -> best spectrum id
    best_by_physical = {}
    for sp_n in passing_spectra:
        pid = spectra[sp_n].physical_id
        if (pid not in best_by_physical
            or expect[sp_n] < expect[best_by_physical[pid]]):
            best_by_physical[pid] = sp_n
    passing_spectra = set(best_by_physical.itervalues())

    # why do this?
    if XTP["output, results"] == "valid":
        # spectra sorted by e value here--does it matter???
        passing_spectra = set(sp_n for sp_n in passing_spectra
                              if (expect[sp_n]
                                  <= 0.95 * XTP["output, maximum valid expectation value"]))

    #debug("passing_spectra: %s" % sorted(list(passing_spectra)))

    # (passing) spectrum id's for spectra that are "repeats"--a repeat is a
    # spectrum for which a better corresponding spectrum exists having the
    # same domain 0 match (start, end, protein id).
    repeats = set()
    for pid, matches in best_protein_matches.iteritems():
        domain_0_matches = {}           # spectrum id -> match
        for m in matches:
            if m[6] not in domain_0_matches:
                domain_0_matches[m[6]] = m
        for sid_x in domain_0_matches:
            if sid_x not in passing_spectra:
                continue
            for sid_y in domain_0_matches:
                if sid_y not in passing_spectra:
                    continue
                if sid_x < sid_y:
                    if (domain_0_matches[sid_x][7:9]
                        == domain_0_matches[sid_y][7:9]):
                        if expect[sid_x] > expect[sid_y]:
                            repeats.add(sid_x)
                        else:
                            repeats.add(sid_y)

    #debug("repeats: %s" % repeats)

    # protein id -> protein expectation value (log10 of expectation)
    raw_protein_expect = {}
    # protein id -> set of spectrum ids used in expectation value
    raw_protein_expect_spectra = {}

    # FIX: rework this and other uses of best_protein_matches to use a protein
    # -> supporting spectra map instead?
    for protein_id, matches in best_protein_matches.iteritems():
        for m in matches:
            spectrum_id = m[6]
            if spectrum_id not in passing_spectra or spectrum_id in repeats:
                continue
            log_expect = math.log10(expect[spectrum_id])
            if log_expect > -1.0:
                continue
            # this avoids adding matches for multiple domains
            if True or spectrum_id not in raw_protein_expect_spectra.get(protein_id, set()):
                raw_protein_expect[protein_id] = (raw_protein_expect.get(protein_id, 0)
                                                  + log_expect)
                raw_protein_expect_spectra.setdefault(protein_id,
                                                      set()).add(spectrum_id)

    debug("raw_protein_expect: %s", raw_protein_expect)
    debug("raw_protein_expect_spectra: %s", raw_protein_expect_spectra)
    
    bias = float(candidate_spectrum_count) / db_residue_count

    # protein id -> protein expectation value (log10 of expectation)
    protein_expect = {}
    for protein_id in raw_protein_expect:
        sp_count = len(raw_protein_expect_spectra[protein_id])
        protein_expect[protein_id] = get_final_protein_expect(sp_count,
                                                              len(passing_spectra),
                                                              int(round(match_ratio)),
                                                              raw_protein_expect[protein_id],
                                                              len(fasta_db),
                                                              candidate_spectrum_count)
        if sp_count > 1:
            protein_length = len(fasta_db[protein_id][1])
            if protein_length * bias < 1.0:
                protein_expect[protein_id] += sp_count * math.log10(protein_length*bias)

    #debug("protein_expect: %s" % protein_expect)

    # protein id -> sum of peak intensities for supporting spectra
    intensity = {}
    # protein id -> set of spectrum ids used in intensity value
    intensity_spectra = {}
    # rework?
    for pid, matches in best_protein_matches.iteritems():
        for m in matches:
            sp_n = m[6]
            if sp_n not in passing_spectra:
                continue
            if sp_n not in intensity_spectra.get(pid, set()):
                intensity[pid] = (intensity.get(pid, 0.0)
                                  + spectra[sp_n].sum_peak_intensity)
                intensity_spectra.setdefault(pid, set()).add(sp_n)

    #debug("intensity: %s" % intensity)

    # set output order
    spec_prot_info_items = spec_prot_info.items()
    if XTP["output, sort results by"] == 'protein':
        # sort spectra by protein expect of spectra's best protein expect
        # for proteins grouped together, sort spectra by peptide start pos
        # (BUG?  instead, we only pos sort for equal protein expects)
        def protein_order_less_than(item):
            return (min(protein_expect.get(protein_id, 1000)
                        for protein_id, domains in item[1]),
                    item[1][0][1][0][7], # domain 0 starting position (ouch)
                    expect[item[0]])    # spectrum expect
        spec_prot_info_items.sort(key=protein_order_less_than)
    else:
        # sort each spectra's proteins by (protein expect, protein id)
        assert False

    return (spec_prot_info_items, expect, protein_expect, intensity,
            survival_curve, line_parameters, passing_spectra)


def print_histogram_XML(hlabel, htype, histogram, a0a1=None):
    #hlength = max(histogram)+1+1        # final bin is 0
    hlength = len(histogram)+1          # final bin is 0
    h0 = histogram + [0]

    print '<GAML:trace label="%s" type="%s">' % (hlabel, htype)
    if a0a1:
        print '<GAML:attribute type="a0">%s</GAML:attribute>' % a0a1[0]
        print '<GAML:attribute type="a1">%s</GAML:attribute>' % a0a1[1]
    print '<GAML:Xdata label="%s" units="score">' % hlabel
    print ('<GAML:values byteorder="INTEL" format="ASCII" numvalues="%s">'
           % hlength)
    print ' '.join([ str(x) for x in range(hlength) ]) + ' '
    print '</GAML:values>'
    print '</GAML:Xdata>'
    print '<GAML:Ydata label="%s" units="counts">' % hlabel
    print ('<GAML:values byteorder="INTEL" format="ASCII" numvalues="%s">'
           % hlength)
    #print ' '.join([ str(histogram.get(x, 0)) for x in range(hlength) ])
    print ' '.join([ str(h0[x]) for x in range(hlength) ]) + ' '
    print '</GAML:values>'
    print '</GAML:Ydata>'
    print '</GAML:trace>'


def print_results_XML(options, XTP, db_info, spectrum_fn,
                      spec_prot_info_items, spectra, expect, protein_expect,
                      intensity, second_best_score, survival_curve,
                      line_parameters, passing_spectra):
    """Output the XTandem-style XML results file, to stdout."""

    print '<?xml version="1.0"?>'
    xslpath = XTP.get("output, xsl path") 
    if xslpath:
        print '<?xml-stylesheet type="text/xsl" href="%s"?>' % xslpath
    print '<bioml xmlns:GAML="http://www.bioml.com/gaml/"',
    title = XTP.get("output, title")
    if title:
        print 'label="%s">' % title
    else:
        print '''label="models from '%s'">''' % spectrum_fn

    #debug("spec_prot_info_items: %s" % spec_prot_info_items)

    # there is an apparently redundant check against
    # "output, maximum valid expectation value" here

    for spectrum_id, spectrum_info in spec_prot_info_items:
        if spectrum_id not in passing_spectra:
            continue

        sp = spectra[spectrum_id]
        assert spectrum_info, "only 'valid' implemented"
        si0 = spectrum_info[0][1][0]
        defline, run_seq, seq_filename = db_info[(si0[11], si0[12])]

        if (XTP["output, proteins"] or XTP["output, histograms"]
            or XTP["output, spectra"]):
            print ('<group id="%s" mh="%.6f" z="%s" expect="%.1e"'
                   ' label="%s" type="model" sumI="%.2f" maxI="%.5e"'
                   ' fI="%s" >'
                   % (sp.id, sp.mass, sp.charge, expect[spectrum_id],
                      abbrev_defline(clean_defline(defline)),
                      math.log10(sp.sum_peak_intensity),
                      sp.max_peak_intensity, sp.normalization_factor))

        if XTP["output, proteins"]:
            for pn, (protein_id, domains) in enumerate(spectrum_info):
                d0_defline, d0_run_seq, d0_seq_filename = db_info[(domains[0][11],
                                                                   domains[0][12])]
                print ('<protein expect="%.1f" id="%s.%s" uid="%s" label="%s"'
                       ' sumI="%.2f" >'
                       % (protein_expect[protein_id], sp.id, pn+1,
                          protein_id+1,
                          abbrev_defline(clean_defline(d0_defline)),
                          math.log10(intensity[protein_id])))
                print ('<note label="description">%s'
                       % clean_defline(d0_defline))
                print '</note>'
                print '<file type="peptide" URL="%s"/>' % d0_seq_filename
                print '<peptide start="1" end="%s">' % len(d0_run_seq)
                if XTP["output, sequences"]:
                    seq = d0_run_seq
                    for rowbegin in range(0, len(seq), 50):
                        rowseq = seq[rowbegin:rowbegin+50]
                        if rowseq:
                            sys.stdout.write('\t')
                        for blockbegin in range(0, len(rowseq), 10):
                            if blockbegin:
                                sys.stdout.write(' ')
                            sys.stdout.write(rowseq[blockbegin:blockbegin+10])
                        print
                for dn, dom in enumerate(domains):
                    dom_defline, dom_run_seq, dom_seq_filename \
                                 = db_info[(dom[11], dom[12])]
                    print ('<domain id="%s.%s.%s" start="%s" end="%s"'
                           ' expect="%.1e" mh="%.4f" delta="%.4f"'
                           ' hyperscore="%.1f" nextscore="%.1f" y_score="%.1f"'
                           ' y_ions="%s" b_score="%.1f" b_ions="%s" pre="%s"'
                           ' post="%s" seq="%s" missed_cleavages="%s">'
                           % (sp.id, pn+1, dn+1, dom[7]+1,
                              dom[7]+len(dom[8]), expect[spectrum_id],
                              dom[10], dom[4]-dom[10],
                              cxtpy.scale_hyperscore(dom[0]),
                              cxtpy.scale_hyperscore(second_best_score.get(spectrum_id, 100)),
                              cxtpy.scale_hyperscore(dict(dom[2])['Y']),
                              dict(dom[3])['Y'],
                              cxtpy.scale_hyperscore(dict(dom[2])['B']),
                              dict(dom[3])['B'],
                              get_prefix_sequence(dom[7], dom[12], dom_run_seq),
                              get_suffix_sequence(dom[7]+len(dom[8]), dom[12],
                                                  dom_run_seq),
                              dom[8], dom[16]))
                    # FIX: print '<aa type="C" modified="42" />'s here (mods)
                    print '</domain>'
                print '</peptide>'
                print '</protein>'

        if XTP["output, histograms"]:

            print '<group label="supporting data" type="support">'
            print_histogram_XML("%s.hyper" % (sp.id),
                                "hyperscore expectation function",
                                survival_curve[spectrum_id],
                                line_parameters[spectrum_id])
            # other histograms nyi (data not collected)
            #print_histogram_XML("%s.convolute" % (sp.id),
            #                    "convolution survival function",
            #                    convolution_histogram[spectrum_id])
            #print_histogram_XML("%s.b" % (sp.id),
            #                    "b ion histogram",
            #                    b_histogram)
            #print_histogram_XML("%s.y" % (sp.id),
            #                    "y ion histogram",
            #                    y_histogram)
            print
            print '</group>'

        if XTP["output, spectra"]:
            print '<group type="support" label="fragment ion mass spectrum">'
            if options.quirks_mode:
                print '<note label="Description">no description</note>'
            else:
                print ('<note label="Description">%s:%s</note>'
                       % (sp.file_id, sp.name))
            print ('<GAML:trace id="%s" label="%s.spectrum"'
                   ' type="tandem mass spectrum">'
                   % (sp.id, sp.id))
            print '<GAML:attribute type="M+H">%s</GAML:attribute>' % sp.mass
            print '<GAML:attribute type="charge">%s</GAML:attribute>' % sp.charge
            print ('<GAML:Xdata label="%s.spectrum" units="MASSTOCHARGERATIO">'
                   % (sp.id))
            print ('<GAML:values byteorder="INTEL" format="ASCII" numvalues="%s">'
                   % len(sp.peaks))
            def rstrip_zeros(s):
                if '.' in s:
                    return s.rstrip('0').rstrip('.')
                return s
            for p in sp.peaks:
                print rstrip_zeros('%.2f' % p.mz),
            print ''
            print '</GAML:values>'
            print '</GAML:Xdata>'
            print ('<GAML:Ydata label="%s.spectrum" units="UNKNOWN">'
                   % (sp.id))
            print ('<GAML:values byteorder="INTEL" format="ASCII" numvalues="%s">'
                   % len(sp.peaks))
            for p in sp.peaks:
                print str(int(round(p.intensity))),
            print ''
            print '</GAML:values>'
            print '</GAML:Ydata>'
            print '</GAML:trace>'
            sys.stdout.write('</group>')

        if (XTP["output, proteins"] or XTP["output, histograms"]
            or XTP["output, spectra"]):
            print '</group>'

    # input parameters
    if XTP["output, parameters"]:
        print '<group label="input parameters" type="parameters">'
        for k in sorted(XTP.keys()):
            v = XTP[k]
            if isinstance(v, list):         # modlist
                v = ','.join([ '%s@%s' % (mod, res) for res, mod in v ])
            elif isinstance(v, bool):
                v = { True : 'yes', False : 'no' }[v]
            print '	<note type="input" label="%s">%s</note>' % (k, v)
        if options.quirks_mode:
            print '	<note type="input" label="xtandem quirks mode">yes</note>'
        print '</group>'

    # performance not implemented

    # if we're using anything other than standard monoisotopic masses, output
    # a group as in mreport::masses (not implemented)

    print '</bioml>'                    # end of XML output
    

def main():
    parser = optparse.OptionParser(usage=
                                   "usage: %prog [options] <parameter-file>",
                                   description=__doc__)
    parser.add_option("-o", "--output", dest="output",
                      help="destination file [default as given in parameter"
                      " file, '-' for stdout]",
                      metavar="FILE")
    parser.add_option("--quirks-mode", action="store_true",
                      dest="quirks_mode",
                      help="try to generate results as close as possible to"
                      " those of XTandem (possibly at the expense of"
                      " accuracy)") 
    parser.add_option("-q", "--quiet", action="store_true",
                      dest="quiet", help="no warnings")
    parser.add_option("-v", "--verbose", action="store_true",
                      dest="verbose", help="be verbose")
    parser.add_option("--debug", action="store_true",
                      dest="debug", help="output debugging info")
    parser.add_option("--profile", action="store_true",
                      dest="profile",
                      help="dump profiling output to './xtpy.prof'")
    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.print_help()
        sys.exit(1)

    log_level = logging.WARNING
    if options.quiet:
        log_level = logging.ERROR
    if options.verbose:
        log_level = logging.INFO
    if options.debug:
        log_level = logging.DEBUG
    logging.basicConfig(level=log_level, datefmt='%b %e %H:%M:%S',
                        format='%(asctime)s %(levelname)s: %(message)s')

    # read params
    parameters = read_xml_parameters(args[0])
    default_parameter_fn = parameters.get("list path, default parameters")
    if default_parameter_fn:
        default_parameters = read_xml_parameters(default_parameter_fn)
        default_parameters.update(parameters)
        parameters = default_parameters
    global XTP
    XTP = validate_parameters(parameters)

    if options.output:
        if options.output != '-':
            sys.stdout = open(options.output, 'w')
    else:
        output_path = XTP["output, path"]
        if output_path:
            sys.stdout = open(output_path, 'w')

    taxonomy = read_taxonomy(XTP["list path, taxonomy information"])

    initialize_spectrum_parameters(options.quirks_mode)
    xtpy_search.CP = CP

    # read sequence dbs
    # [(defline, seq, filename), ...]
    fasta_db = list(read_fasta_files(taxonomy[XTP["protein, taxon"]]))
    # [(idno, offset, defline, seq, is_N, is_C, seq_filename), ...]
    db = []
    for idno, (defline, sequence, filename) in enumerate(fasta_db):
        db.extend(split_sequence_into_aa_runs(idno, defline, sequence,
                                              filename))
    db_residue_count = sum(len(dbi[3]) for dbi in db)

    # (idno, offset) -> (defline, run_seq, filename)
    db_info = dict(((idno, offset), (defline, seq, filename))
                   for idno, offset, defline, seq, is_N, is_C, filename in db)

    info("read %s sequences (%s runs, %s residues)", len(fasta_db), len(db),
         db_residue_count)
    if not db:
        error("no database sequences")

    # read spectra
    spectrum_fn = XTP["spectrum, path"]
    spectra = list(cxtpy.spectrum.read_spectra(open(spectrum_fn), 0))
    if not spectra:
        error("no input spectra")
    info("read %s spectra", len(spectra))

    spectra.sort(key=lambda x: x.mass)

    # filter and normalize spectra
    # NI: optionally using "contrast angle" (mprocess::subtract)
    if XTP["spectrum, use conditioning"]:
        spectra = [ s for s in spectra
                    if s.filter_and_normalize(XTP["spectrum, minimum fragment mz"],
                                              XTP["spectrum, dynamic range"],
                                              XTP["spectrum, minimum peaks"],
                                              XTP["spectrum, total peaks"]) ]
    info("     %s spectra after filtering", len(spectra))

    cxtpy.spectrum.set_searchable_spectra(spectra)
    #sys.exit('len = %s' % len(cxtpy.cvar.spectrum_searchable_spectra))
    #print >> sys.stderr, cxtpy.cvar.spectrum_spectrum_mass_index.keys()

    spectrum_mass_index = [ (sp.mass, n) for n, sp in enumerate(spectra) ]
    spectrum_mass_index.sort()

    # spectrum index is wrt 'spectra' above
    # spectrum index -> [ <match info>, ... ]
    best_match = {}
    # spectrum index -> best_hyperscore
    best_score = {}
    # spectrum index -> 2nd-best hyperscore
    second_best_score = {}
    # spectrum index -> (scaled, binned hyperscore -> count)
    hyperscore_histogram = {}

    candidate_spectrum_count = 0
    for idno, offset, defline, seq, is_N, is_C, seq_filename in db:
        if options.verbose:
            sys.stderr.write('p')

        cleavage_points = list(generate_cleavage_points(seq))
        for begin, end, missed_cleavage_count \
                in generate_peptides(seq, cleavage_points,
                                     XTP["scoring, maximum missed cleavage sites"]):
            peptide_seq = seq[begin:end]
            debug('generated peptide: %s', peptide_seq)
            candidate_spectrum_count \
                += xtpy_search.search_peptide(spectrum_mass_index,
                                              spectra, idno, offset, begin,
                                              peptide_seq, is_N, is_C,
                                              missed_cleavage_count,
                                              hyperscore_histogram,
                                              best_score, best_match,
                                              second_best_score)
    if options.verbose:
        print >> sys.stderr

    info('%s candidate spectra examined', candidate_spectrum_count)
    info('processing results')

    # filter results and calculate statistics
    (spec_prot_info_items, expect, protein_expect, intensity, survival_curve,
     line_parameters, passing_spectra) \
         = process_results(best_match, hyperscore_histogram, fasta_db,
                           candidate_spectrum_count, spectra, db_residue_count)

    info('writing results')

    print_results_XML(options, XTP, db_info, spectrum_fn,
                      spec_prot_info_items, spectra, expect, protein_expect,
                      intensity, second_best_score, survival_curve,
                      line_parameters, passing_spectra)

    info('finished')
    logging.shutdown()


if __name__ == '__main__':
#     try:
#         import psyco
#         psyco.full()
#         ###psyco.bind(formula_mass)
#         warning('psyco enabled')
#     except ImportError:
#         pass

    try:
        if '--profile' in sys.argv:
            import hotshot
            import hotshot.stats
            data_fn = "xtpy.prof.tmp"
            prof = hotshot.Profile(data_fn)
            prof.runcall(main)
            prof.close()

            sys.stdout = open("xtpy.prof", 'w')

            stats = hotshot.stats.load(data_fn)
            stats.strip_dirs()
            stats.sort_stats('cumulative')
            stats.print_stats(100)
            stats.sort_stats('time', 'calls')
            stats.print_stats(100)
            os.remove(data_fn)
        else:
            main()
    except SystemExit:
        raise
    except:
        logging.exception("unhandled exception")


# FIXES:
# - need to rigorously check for bad input in any file (fasta, spectrum, xml)
# - escape special XML chars


# NOTES:

# pyro_check: if no N-terminal mod otherwise specified, look for potential
# mods for these N-terminal residues: Q or C+57 -> loss of ammonia, E -> loss
# of water
