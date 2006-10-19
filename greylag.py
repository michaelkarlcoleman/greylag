#!/usr/bin/env python

'''Analyze mass spectra and assign peptides and proteins.  (This is a partial
   re-implementation of the X!Tandem algorithm, with many extra features.)
'''


### "Simplicity is prerequisite for reliability" - Edsger W. Dijkstra ###


__version__ = "$Id$"


import cPickle
import fileinput
import itertools
import logging
from logging import debug, info, warning
import math
import optparse
import os
from pprint import pprint
import re
import sys

import elementtree.ElementTree

import cgreylag


# Try to drop dead immediately on interrupt (control-C), instead of normal
# Python KeyboardInterrupt processing, since we may spend long periods of time
# in uninterruptible C++ calls.
try:
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)
except:
    pass


def error(s, *args):
    logging.error(s, *args)
    sys.exit(1)
def fileerror(s, *args):
    error(s + (", at line %s of file '%s'"
               % (fileinput.filelineno(), fileinput.filename())),
          *args)


# name -> value map of processed XML input parameters
XTP = {}

# handle to the singleton parameter object shared with the C++ module
CP = cgreylag.cvar.parameters_the

# FIX: is this from H1 or H-avg??  (possible error here is ~0.0007 amu)
PROTON_MASS =   1.007276
ELECTRON_MASS = 0.000549                # ?

# Reference values from NIST (http://physics.nist.gov/PhysRefData/)
ATOMIC_MASS = {
    'H'     :  1.00782503214,
    'H-avg' :  1.007947,
    'C'     : 12.00000000,
    'C-avg' : 12.01078,
    'N'     : 14.00307400529,
    'N-avg' : 14.00672,
    'N15'   : 15.00010889849,
    'O'     : 15.994914622115,
    'O-avg' : 15.99943,
    'P'     : 30.9737615120,
    'P-avg' : 30.9737612,
    'S'     : 31.9720706912,
    'S-avg' : 32.0655,
    }

ATOMIC_MASS_XT = {
    'H'     :  1.007825035,
    'H-avg' :  1.00794,
    'C'     : 12.0,
    'C-avg' : 12.0107,
    'N'     : 14.003074,
    'N-avg' : 14.0067,
    'O'     : 15.99491463,
    'O-avg' : 15.9994,
    'P'     : 30.973762,
    'P-avg' : 30.973761,
    'S'     : 31.9720707,
    'S-avg' : 32.065,
    }

# FIX!
ATOMIC_MASS = ATOMIC_MASS_XT


def formula_mass(formula, atomic_mass=ATOMIC_MASS):
    """Return the mass of formula, using the given mass regime (monoisotopic
    by default)."""
    parts = [ p or '1' for p in re.split(r'([A-Z][a-z]*)', formula)[1:] ]
    # parts for alanine is ('C', '3', 'H', '5', 'O', '1', 'N', '1')
    return sum(atomic_mass[parts[i]] * int(parts[i+1])
               for i in range(0, len(parts), 2))

# FIX: selenocysteine (U), etc
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

STANDARD_RESIDUE_MASS = dict((residue, formula_mass(RESIDUE_FORMULA[residue]))
                             for residue in RESIDUE_FORMULA)

# Why aren'te these two sets of average masses equal?
STANDARD_XT_AVG_RESIDUE_MASS = {
    'A' :  71.0788,
    'C' : 103.1388,
    'D' : 115.0886,
    'E' : 129.1155,
    'F' : 147.1766,
    'G' :  57.0519,
    'H' : 137.1411,
    'I' : 113.1594,
    'K' : 128.1741,
    'L' : 113.1594,
    'M' : 131.1926,
    'N' : 114.1038,
    'P' :  97.1167,
    'Q' : 128.1307,
    'R' : 156.1875,
    'S' :  87.0782,
    'T' : 101.1051,
    'V' :  99.1326,
    'W' : 186.2132, 
    'Y' : 163.1760,
    }

STANDARD_AVG_RESIDUE_MASS = dict((residue,
                                  formula_mass(RESIDUE_FORMULA[residue],
                                               { 'H' : ATOMIC_MASS['H-avg'],
                                                 'C' : ATOMIC_MASS['C-avg'],
                                                 'N' : ATOMIC_MASS['N-avg'],
                                                 'O' : ATOMIC_MASS['O-avg'],
                                                 'P' : ATOMIC_MASS['P-avg'],
                                                 'S' : ATOMIC_MASS['S-avg'],
                                                 }))
                                 for residue in RESIDUE_FORMULA)


def initialize_spectrum_parameters(quirks_mode):
    """Initialize parameters known to the spectrum module."""

    # This is the size of vectors that are indexed by residues (A-Z) or
    # special characters ('[]').
    RESIDUE_LIMIT = max(ord(c) for c in 'Z[]') + 1

    # These two aren't currently part of the regime, so we're not handling
    # deuterium yet (pointless?)
    CP.proton_mass = PROTON_MASS
    CP.hydrogen_mass = formula_mass("H")

    # FIX: only a mono/mono regime 0 implemented for now
    mono_regime = cgreylag.mass_regime_parameters()

    mono_regime.hydroxyl_mass = formula_mass("OH")
    mono_regime.water_mass = formula_mass("H2O")
    mono_regime.ammonia_mass = formula_mass("NH3")

    mono_regime.residue_mass.resize(128)
    for residue in RESIDUE_FORMULA:
        mono_regime.residue_mass[ord(residue)] = STANDARD_RESIDUE_MASS[residue]

    # FIX: for the moment we don't differentiate the parent/fragment cases
    mono_regime.modification_mass.resize(128)
    for residue, modvalue in XTP["residue, modification mass"]:
        mono_regime.modification_mass[ord(residue)] = modvalue

    # SWIG currently only exposes the outermost vector as a modifiable object.
    # Inner vectors appear as tuples, and are thus unmodifiable.  They must
    # therefore be assigned all at once.  This shortcoming will probably be
    # fixed in a future version of SWIG.

    rpmm = XTP["residue, potential modification mass"]
    mono_regime.potential_modification_mass.resize(len(rpmm))
    info("searching %s potential mod alternative set(s)", len(rpmm))
    for altn, alternative in enumerate(rpmm):
        v = [ [] for i in range(128) ]
        for residue, modvalue in alternative:
            v[ord(residue)].append(modvalue)
        mono_regime.potential_modification_mass[altn] = v

    #for residue, modvalue in XTP["refine, potential modification mass"]:
    #    mono_regime.potential_modification_mass_refine[0][ord(residue)] \
    #        = mono_regime.potential_modification_mass_refine[0][ord(residue)] + (modvalue,)

    # regime 0 is mono/mono
    CP.parent_mass_regime.append(mono_regime);
    CP.fragment_mass_regime.append(mono_regime);

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
    CP.spectrum_synthesis = XTP["refine, spectrum synthesis"]

    # CP.factorial[n] == (double) n!
    CP.factorial.resize(100, 1.0)
    for n in range(2, len(CP.factorial)):
        CP.factorial[n] = CP.factorial[n-1] * n

    CP.quirks_mode = bool(quirks_mode)
    CP.hyper_score_epsilon_ratio = 0.999 # must be <1

    #debug("CP: %s", pythonize_swig_object(CP, ['the']))
    

def cleavage_motif_re(motif):
    """Return (regexp, pos), where regexp is a regular expression that will
    match a cleavage motif, and pos is the position of the cleavage with
    respect to the match (e.g., 1 for '[KR]|{P}', or None if absent).  (The RE
    actually matches one character, the rest matching as lookahead, which
    means that re.finditer will find all overlapping matches.)
    """
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
    if len(re_parts) == 0:
        re_pattern = '.'
    elif len(re_parts) == 1:
        re_pattern = re_parts[0]
    else:
        re_pattern = re_parts[0] + '(?=' + ''.join(re_parts[1:]) + ')'
    return (re_pattern, cleavage_pos)


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
    for begin_i in xrange(len_cp-1):
        for end_i in xrange(begin_i+1,
                            min(len_cp,
                                begin_i+2+maximum_missed_cleavage_sites)):
            begin, end = cleavage_points[begin_i], cleavage_points[end_i]
            if end - begin >= 5:        # XT says 4, but means 5
                yield begin, end, end_i-begin_i-1


def generate_cleavage_points(cleavage_re, cleavage_pos, sequence):
    """Yields the offsets of the cleavages in sequence."""
    yield 0
    for m in cleavage_re.finditer(sequence, 1, len(sequence)-2):
        yield m.start() + cleavage_pos
    yield len(sequence)


aa_sequence = re.compile(r'[ARNDCQEGHILKMFPSTWYV]+')

def split_sequence_into_aa_runs(idno, defline, sequence, filename):
    """Returns a tuple (idno, start, defline, seq, filename) for each
    contiguous run of residues in sequence, where 'start' is the position of
    'seq' in 'sequence'."""
    matches = list(aa_sequence.finditer(sequence))
    return [ (idno, m.start(), defline, m.group(), filename)
             for n, m in enumerate(matches) ]


def clean_defline(s):
    """Return the given string with tags replaced by spaces and control
    characters removed, then stripped."""
    return re.sub(r'[^ -~]', '', s.replace('\t', ' ')).strip()


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
    root = elementtree.ElementTree.ElementTree(file=filename).getroot()
    return dict( (taxon.get("label"),
                  [ f.get("URL") for f in taxon.findall('./file') ])
                 for taxon in root.findall('taxon') )


def read_xml_parameters(filename):
    """Return a map of parameters to values, per parameter file fn."""
    root = elementtree.ElementTree.ElementTree(file=filename).getroot()
    return dict((e.get("label"), e.text)
                for e in root.findall('note')
                if e.get("type") == "input")


# XML parameter file processing

def mod_list(modification_list_specification, unique_mods=True):
    """Check and return a sorted list of (residue, modification) tuples.  If
    unique_mods is True, also check that each residue is specified only once.
    """
    modification_list_specification = modification_list_specification.strip()
    if not modification_list_specification:
        return []
    speclist = []
    for modspec in modification_list_specification.split(','):
        parts = modspec.strip().split('@')
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

def alternative_mod_list(modification_list_specification):
    """Check and return a list of lists of (residue, modification) tuples.
    The top-level list corresponds to the ';'-separated alternatives.
    """
    alist = [ mod_list(alternative_specification, unique_mods=False)
              for alternative_specification
              in modification_list_specification.split(';') ]
    return [ a for a in alist if a ]    # omit empty alternatives
        

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
    "protein, C-terminal residue modification mass" : (float, "0.0", p_ni_equal(0)), # eliminate, tandem ni
    "protein, N-terminal residue modification mass" : (float, "0.0", p_ni_equal(0)), # eliminate, tandem ni
    "protein, cleavage C-terminal mass change" : (float, formula_mass("OH")),
    "protein, cleavage N-terminal limit" : (int, "100000000", p_positive),
    "protein, cleavage N-terminal mass change" : (float, formula_mass("H")),
    "protein, cleavage semi" : (bool, "no", p_ni_equal(False)),
    "protein, cleavage site" : (str, None),
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
    "refine, potential modification mass" : (alternative_mod_list, "", p_ni_empty),
    "refine, potential modification motif" : (str, "", p_ni_empty),
    "refine, sequence path" : (str, "", p_ni_empty),
    "refine, spectrum synthesis" : (bool, "no"),
    "refine, tic percent" : (float, 20.0, p_nonnegative), # ignored
    "refine, unanticipated cleavage" : (bool, p_ni_equal(False)),
    "refine, use potential modifications for full refinement" : (bool, "no", p_ni_equal(False)),
    "residue, modification mass" : (mod_list, ""),
    "residue, potential modification mass" : (alternative_mod_list, ""),
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
                debug("parameter '%s' defaulting to '%s'", p_name, default)
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

    # FIX: where should this go?
    if pmap["spectrum, fragment mass error"] > 0.5:
        warning("'spectrum, fragment mass error' is %s",
                pmap["spectrum, fragment mass error"])

    return pmap


# FIX: This code really seems iffy.  It may diverge from xtandem and/or
# xtandem's implementation may simply be buggy or ill-conceived.
def get_spectrum_expectation(hyper_score, histogram):
    """Return (expectation, survival curve, (a0, a1)) for this hyperscore and
    histogram, where a0 and a1 are the intercept and slope, respectively, for
    the least-squares fit to the survival curve, used to predict the
    expectation.
    """
    scaled_hyper_score = cgreylag.scale_hyperscore(hyper_score)

    # trim zeros from right end of histogram (which must contain a non-zero)
    histogram = list(histogram)
    while histogram and histogram[-1] == 0:
        histogram.pop()

    counts = sum(histogram)
    c = counts
    # survival is just a reversed cumulative distribution
    survival = [0] * len(histogram)
    for i in xrange(len(histogram)):
        survival[i] = c
        c -= histogram[i]
    # note: neither histogram nor survival have 0 as rightmost element here

    # this next looks bogus, but use for now to replicate xtandem
    # "remove potentially valid scores from the stochastic distribution"
    # makes the survival function non-monotonic?
    if survival:
        lPos = survival[0] / 5
        for lMid in xrange(len(survival)):
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

    if counts < 200 or len(survival) < 3:
        # use default line
        a0, a1 = 3.5, -0.18
        return 10.0 ** (a0 + a1 * scaled_hyper_score), survival, (a0, a1)
    
    survival += [0]                     # FIX
    min_limit = 10
    max_limit = int(round(survival[0]/2.0))
    try:
        max_i = min(i for i in xrange(len(survival))
                    if survival[i] <= max_limit)
        # this assumes rightmost element is < min_limit!
        min_i = min(i for i in xrange(max_i, len(survival))
                    if survival[i] <= min_limit)
    except ValueError:
        # FIX: how would this ever happen?
        warning('bad survival curve? %s', survival)
        a0, a1 = 3.5, -0.18
        return 10.0 ** (a0 + a1 * scaled_hyper_score), survival, (a0, a1)
        
    data_X = range(max_i, min_i)
    data_Y = [ math.log10(survival[x]) for x in data_X ]
    # FIX!
    if len(data_Y) < 2 or data_Y[0] != max(data_Y):
        warning('bad survival curve? (2) [%s] %s', len(data_Y), survival)
        a0, a1 = 3.5, -0.18
        return 10.0 ** (a0 + a1 * scaled_hyper_score), survival, (a0, a1)

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
    #debug('pe: %s', (sp_count, valid_spectra, match_ratio, raw_expect,
    #                 db_length, candidate_spectra))

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
    #debug('valid_spectra: %s, match_ratio: %s', valid_spectra, match_ratio)
    r -= math.log10(valid_spectra) + (sp_count-1) * math.log10(match_ratio)
    #debug('x: %s', (match_ratio, candidate_spectra))
    p = min(float(match_ratio) / candidate_spectra, 0.9999999)
    #debug("p: %s", p)
    r += (sp_count * math.log10(p)
          + (valid_spectra - sp_count) * math.log10(1 - p))
    return r


# def read_ms2_spectrum_masses(spectrum_fns):
#     """Return a list of all parent masses present in the given ms2 files."""
#     masses = []
#     for fn in spectrum_fns:
#         header_line = False
#         for line in open(fn):
#             if line.startswith(':'):
#                 header_line = True
#             elif header_line:
#                 masses.append(float(line.split(None, 1)[0]))
#     return masses
    

def generate_mass_bands(band_count, mass_list):
    """Yield (n, mass_lb, mass_ub) for each mass band, where n ranges from 1
    to band_count.  To generate the bands, the mass list (which is assumed
    already sorted by mass) is evenly partitioned into bands with masses in
    the range [mass_lb, mass_ub).
    """
    assert mass_list and band_count > 0
    band_size = len(mass_list) / band_count + 1
    lb = mass_list[0]
    for bn in range(1, band_count):
        i = min(bn*band_size, len(mass_list)-1)
        ub = mass_list[i]
        yield bn, lb, ub
        lb = ub
    # Since the upper bound is exclusive, we add a little slop to make sure
    # the last spectrum is included.
    yield band_count, lb, round(mass_list[-1]+1)


def filter_ms2_by_mass(f, lb, ub):
    """Yield lines from f (an open ms2 file), zeroing out any spectra with
    mass outside [lb, ub).  Specifically, a zeroed spectrum will have an empty
    name and a mass and charge of zero.  If all charges for a physical
    spectrum are zeroed, its peaklist will be replaced with a single peak
    having mass and intensity zero.  Thus the zeroed spectra are validly
    formatted placeholders.
    """
    line = f.readline()
    while line:
        if not line.startswith(':'):
            error("bad ms2 format: missing header lines?")
        zero_peaks = True
        headers = []
        while line.startswith(':'):
            name = line
            mass_charge = f.readline()
            try:
                mass = float(mass_charge.split()[0])
            except:
                error("bad ms2 format: bad mass")
            if lb <= mass < ub:
                zero_peaks = False
                yield name
                yield mass_charge
            else:
                yield ':\n'
                yield '0 0\n'
            line = f.readline()
        if not zero_peaks:
            while line and not line.startswith(':'):
                yield line
                line = f.readline()
        else:
            yield '0 0\n'
            while line and not line.startswith(':'):
                line = f.readline()
    

class struct:
    "generic struct class used by pythonize_swig_object"
    def __repr__(self):
        return 'struct(%s)' % self.__dict__


def pythonize_swig_object(o, skip_methods=[]):
    """Creates a pure Python copy of a SWIG object, so that it can be
    easily pickled, or printed (for debugging purposes).  If provided,
    'skip_methods' is a list of methods not to include--this is a hack to
    avoid infinite recursion.
    """

    if isinstance(o, str):
        return o
    try:
        len(o)
    except TypeError:
        pass
    else:
        return list(pythonize_swig_object(x) for x in o)
    if hasattr(o, '__swig_getmethods__'):
        s = struct()
        for a in o.__swig_getmethods__:
            if a not in skip_methods:
                setattr(s, a, pythonize_swig_object(getattr(o, a)))
        return s
    return o


def filter_matches(score_statistics):
    """Filter out any close-but-not-quite matches."""
    if not CP.quirks_mode:
        for sp_n in xrange(len(score_statistics.best_match)):
            bs = score_statistics.best_score[sp_n]
            score_statistics.best_match[sp_n] \
                = [ m for m in score_statistics.best_match[sp_n]
                    if m.hyper_score/bs > CP.hyper_score_epsilon_ratio ]
    # else there shouldn't be anything to filter


def merge_score_statistics(ss0, ss1, offset):
    """Merge ss1 into ss0, keeping the best of both.  (Currently this is
    trivial since the statistics are for distinct (adjacent) sets of spectra.)
    Also, offset is added to the spectrum_index member of all matches.
    """
    for sp_n in xrange(len(ss1.best_match)):
        for m in ss1.best_match[sp_n]:
            m.spectrum_index += offset
    ss0.candidate_spectrum_count += ss1.candidate_spectrum_count
    ss0.hyperscore_histogram.extend(ss1.hyperscore_histogram)
    ss0.second_best_score.extend(ss1.second_best_score)
    ss0.best_score.extend(ss1.best_score)
    ss0.best_match.extend(ss1.best_match)


def find_repeats(best_protein_matches, passing_spectra, expect):
    """Find (passing) spectrum id's for spectra that are repeats.  A repeat is
    a spectrum for which a better corresponding spectrum exists having the
    same domain 0 match (start, end, protein id).
    """
    repeats = set()
    for pid, matches in best_protein_matches.iteritems():
        domain_0_matches = {}           # spectrum id -> match
        for m in matches:
            if m.spectrum_index not in domain_0_matches:
                domain_0_matches[m.spectrum_index] = m
        for sid_x in domain_0_matches:
            if sid_x not in passing_spectra:
                continue
            for sid_y in domain_0_matches:
                if sid_y not in passing_spectra:
                    continue
                if sid_x < sid_y:
                    if ((domain_0_matches[sid_x].peptide_begin
                        == domain_0_matches[sid_y].peptide_begin)
                        and (domain_0_matches[sid_x].peptide_sequence
                             == domain_0_matches[sid_y].peptide_sequence)):
                        if expect[sid_x] > expect[sid_y]:
                            repeats.add(sid_x)
                        else:
                            repeats.add(sid_y)
    return repeats


def process_results(score_statistics, fasta_db, spectra, db_residue_count):

    # spectrum index -> [ (protein id, [domain info, ...]), ... ]
    # domain info == (...)
    spec_prot_info = {}

    for sp_n in xrange(len(score_statistics.best_match)):
        prot_info = {}
        for m in score_statistics.best_match[sp_n]:
            protein_id = m.sequence_index
            prot_info.setdefault(protein_id, []).append(m)
        if prot_info:
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

    for sp_n in xrange(len(score_statistics.best_match)):
        # all hyperscores the same, so choose the first
        # all charges the same (by assumption above)
        if not score_statistics.best_match[sp_n]:
            continue
        sp_hyper_score = score_statistics.best_score[sp_n]
        hh = score_statistics.hyperscore_histogram[sp_n]
        best_histogram[sp_n] = hh
        expect[sp_n], survival_curve[sp_n], line_parameters[sp_n] \
                      = get_spectrum_expectation(sp_hyper_score, hh)

    #debug("expect: %s" % expect)

    # these are the ids of spectra considered "good enough"
    passing_spectra = set(sp_n for sp_n, e in expect.iteritems()
                          if e <= XTP["refine, maximum valid expectation value"])

    #debug("passing_spectra: %s" % sorted(list(passing_spectra)))

    # FIX: is this redundant with spec_prot_info??
    # protein id -> list of spectra match info
    best_protein_matches = {}
    for sp_n in xrange(len(score_statistics.best_match)):
        for m in score_statistics.best_match[sp_n]:
            best_protein_matches.setdefault(m.sequence_index, []).append(m)

    #debug("best_protein_matches: %s" % best_protein_matches)

    # calculate exp for each protein
    valid_spectra_count = sum(1 for sp_n in passing_spectra
                              if (expect[sp_n]
                                  <= XTP["output, maximum valid expectation value"]))
    best_histogram_sum = sum(sum(h) for sp_n, h in best_histogram.iteritems()
                             if sp_n in passing_spectra)
    match_ratio = 0
    #debug("best_histogram_sum: %s passing_spectra: %s", best_histogram_sum,
    #      passing_spectra)
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

    repeats = find_repeats(best_protein_matches, passing_spectra, expect)
    #debug("repeats: %s" % repeats)

    # protein id -> protein expectation value (log10 of expectation)
    raw_protein_expect = {}
    # protein id -> set of spectrum ids used in expectation value
    raw_protein_expect_spectra = {}

    # FIX: rework this and other uses of best_protein_matches to use a protein
    # -> supporting spectra map instead?
    for protein_id, matches in best_protein_matches.iteritems():
        for m in matches:
            spectrum_id = m.spectrum_index
            if spectrum_id not in passing_spectra or spectrum_id in repeats:
                continue
            log_expect = math.log10(expect[spectrum_id])
            # FIX: uncomment!!!
            #if log_expect > -1.0:
            #    continue
            # this avoids adding matches for multiple domains
            # FIX: remove 'True or'!!!
            if True or spectrum_id not in raw_protein_expect_spectra.get(protein_id, set()):
                raw_protein_expect[protein_id] = (raw_protein_expect.get(protein_id, 0)
                                                  + log_expect)
                raw_protein_expect_spectra.setdefault(protein_id,
                                                      set()).add(spectrum_id)

    #debug("raw_protein_expect: %s", raw_protein_expect)
    #debug("raw_protein_expect_spectra: %s", raw_protein_expect_spectra)
    
    bias = float(score_statistics.candidate_spectrum_count) / db_residue_count

    # protein id -> protein expectation value (log10 of expectation)
    protein_expect = {}
    for protein_id in raw_protein_expect:
        sp_count = len(raw_protein_expect_spectra[protein_id])
        protein_expect[protein_id] \
            = get_final_protein_expect(sp_count, len(passing_spectra),
                                       int(round(match_ratio)),
                                       raw_protein_expect[protein_id],
                                       len(fasta_db),
                                       score_statistics.candidate_spectrum_count)
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
            sp_n = m.spectrum_index
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
                    item[1][0][1][0].peptide_begin, # domain 0 starting position (ouch)
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


def print_results_XML(options, XTP, db_info, spectrum_fns,
                      spec_prot_info_items, spectra, expect, protein_expect,
                      intensity, survival_curve, line_parameters,
                      passing_spectra, score_statistics):
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
        print 'label="models from %s">' % (spectrum_fns,)

    #debug("spec_prot_info_items: %s" % spec_prot_info_items)

    # there is an apparently redundant check against
    # "output, maximum valid expectation value" here

    for spectrum_id, spectrum_info in spec_prot_info_items:
        if spectrum_id not in passing_spectra:
            continue

        sp = spectra[spectrum_id]
        assert spectrum_info, "only 'valid' implemented"
        si0 = spectrum_info[0][1][0]
        defline, run_seq, seq_filename = db_info[(si0.sequence_index, si0.sequence_offset)]

        if (XTP["output, proteins"] or XTP["output, histograms"]
            or XTP["output, spectra"]):
            print ('<group id="%s" mh="%.6f" z="%s" expect="%.1e"'
                   ' label="%s" type="model" sumI="%.2f" maxI="%g"'
                   ' fI="%g" >'
                   % (sp.id, sp.mass, sp.charge, expect[spectrum_id],
                      abbrev_defline(clean_defline(defline)),
                      math.log10(sp.sum_peak_intensity),
                      sp.max_peak_intensity, sp.normalization_factor))

        if XTP["output, proteins"]:
            for pn, (protein_id, domains) in enumerate(spectrum_info):
                d0_defline, d0_run_seq, d0_seq_filename = db_info[(domains[0].sequence_index,
                                                                   domains[0].sequence_offset)]
                if protein_id not in protein_expect:
                    warning("protein id '%s' not in protein_expect",
                            protein_id)
                    continue
                if protein_id not in intensity:
                    warning("protein id '%s' not in intensity", protein_id)
                    continue
                            
                print ('<protein expect="%.1f" id="%s.%s" uid="%s" label="%s"'
                       ' sumI="%.2f" >'
                       % (protein_expect[protein_id], sp.id, pn+1,
                          protein_id+1,
                          abbrev_defline(clean_defline(d0_defline)),
                          math.log10(intensity[protein_id])))
                print ('<note label="description">%s</note>'
                       % clean_defline(d0_defline))
                print '<file type="peptide" URL="%s"/>' % d0_seq_filename
                print '<peptide start="1" end="%s">' % len(d0_run_seq)
                if XTP["output, sequences"]:
                    seq = d0_run_seq
                    for rowbegin in xrange(0, len(seq), 50):
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
                                 = db_info[(dom.sequence_index, dom.sequence_offset)]

                    delta_precision = 4
                    delta = sp.mass - dom.peptide_mass
                    if CP.quirks_mode and abs(delta) > 0.01:
                        delta_precision = 3
                    print ('<domain id="%s.%s.%s" start="%s" end="%s"'
                           ' expect="%.1e" mh="%.*f" delta="%.*f"'
                           ' hyperscore="%.1f" nextscore="%.1f" y_score="%.1f"'
                           ' y_ions="%s" b_score="%.1f" b_ions="%s" pre="%s"'
                           ' post="%s" seq="%s" missed_cleavages="%s">'
                           % (sp.id, pn+1, dn+1, dom.peptide_begin+1,
                              dom.peptide_begin+len(dom.peptide_sequence), expect[spectrum_id],
                              delta_precision, dom.peptide_mass,
                              delta_precision, delta,
                              cgreylag.scale_hyperscore(score_statistics.best_score[spectrum_id]),
                              cgreylag.scale_hyperscore(score_statistics.second_best_score[spectrum_id]),
                              cgreylag.scale_hyperscore(dom.ion_scores[cgreylag.ION_Y]),
                              dom.ion_peaks[cgreylag.ION_Y],
                              cgreylag.scale_hyperscore(dom.ion_scores[cgreylag.ION_B]),
                              dom.ion_peaks[cgreylag.ION_B],
                              get_prefix_sequence(dom.peptide_begin, dom.sequence_offset, dom_run_seq),
                              get_suffix_sequence(dom.peptide_begin+len(dom.peptide_sequence), dom.sequence_offset,
                                                  dom_run_seq),
                              dom.peptide_sequence, dom.missed_cleavage_count))
                    # print static mods
                    # FIX: handle '[', ']' too
                    for i, c in enumerate(dom.peptide_sequence):
                        delta = CP.fragment_mass_regime[dom.mass_regime].modification_mass[ord(c)]
                        if delta:
                            print ('<aa type="%s" at="%s" modified="%s" />'
                                   % (c, dom.peptide_begin+i+1, delta))
                    # print potential mods
                    mt_items = list(dom.mass_trace)
                    mt_items.reverse()
                    for mt_item in mt_items:
                        mt_item_pos = mt_item.position
                        if mt_item_pos == cgreylag.POSITION_NTERM:
                            mt_item_pos = 0
                        elif mt_item_pos == cgreylag.POSITION_CTERM:
                            mt_item_pos = len(dom.peptide_sequence)-1
                        if mt_item_pos >= 0:
                            print ('<aa type="%s" at="%s" modified="%s" />'
                                   % (dom.peptide_sequence[mt_item_pos],
                                      dom.peptide_begin+mt_item_pos+1,
                                      mt_item.delta))
                        #    print ('<!-- modified="%s" description="%s" -->'
                        #           % (mt_item.delta, mt_item.description))
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
                       % (spectrum_fns[sp.file_id], sp.name))
            print ('<GAML:trace id="%s" label="%s.spectrum"'
                   ' type="tandem mass spectrum">'
                   % (sp.id, sp.id))
            print '<GAML:attribute type="M+H">%s</GAML:attribute>' % sp.mass
            print '<GAML:attribute type="charge">%s</GAML:attribute>' % sp.charge
            print ('<GAML:Xdata label="%s.spectrum" units="MASSTOCHARGERATIO">'
                   % (sp.id))
            print ('<GAML:values byteorder="INTEL" format="ASCII" numvalues="%s">'
                   % len(sp.peaks))
            #def rstrip_zeros(s):
            #    if '.' in s:
            #        return s.rstrip('0').rstrip('.')
            #    return s
            #for p in sp.peaks:
            #    print rstrip_zeros('%.2f' % p.mz),
            for p in sp.peaks:
                print '%g' % p.mz,
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


def zopen(filename, mode='r', compresslevel=None):
    """Open a filename as with 'open', but using compression if indicated by
    the filename suffix.  The compression level defaults to 1 for .gz files
    and 9 for .bz2 files (presumably bzip2 will only be used if size matters
    much more than compression time)."""
    if filename.endswith('.gz'):
        import gzip
        return gzip.GzipFile(filename, mode, compresslevel or 1)
    elif filename.endswith('.bz2'):
        import bz2
        return bz2.BZ2File(filename, mode, compresslevel or 9)
    else:
        return open(filename, mode)


def main():
    parser = optparse.OptionParser(usage=
                                   "usage: %prog [options] <parameter-file>"
                                   " [<ms2-file>...]",
                                   description=__doc__)
    parser.add_option("-o", "--output", dest="output",
                      help="destination file [default as given in parameter"
                      " file, '-' for stdout]", metavar="FILE")
    parser.add_option("--quirks-mode", action="store_true",
                      dest="quirks_mode",
                      help="try to generate results as close as possible to"
                      " those of X!Tandem (possibly at the expense of"
                      " accuracy)") 
    parser.add_option("--part-split", dest="part_split", type="int",
                      help="split input into M parts, to prepare for"
                      " --part runs [NOTE: the same parameter file and same"
                      " spectrum files (in the same order) must be specified"
                      " for all --part* steps]", metavar="M")  
    parser.add_option("--part", dest="part",
                      help="search one part, previously created with"
                      " --part-split; e.g. '1of2' and '2of2'", metavar="NofM")
    parser.add_option("--part-merge", dest="part_merge", type="int",
                      help="merge the previously searched M results and"
                      " continue", metavar="M")  
    default_prefix = 'greylag'
    parser.add_option("--part-prefix", dest="part_prefix",
                      default=default_prefix,
                      help="prefix to use for temporary part files"
                      " [default='%s']" % default_prefix, metavar="PREFIX")
    parser.add_option("--compress-level", dest="compress_level", type="int",
                      help="compression level to use for compressed files"
                      " created [default=1 for *.gz, 9 for *.bz2]",
                      metavar="N")
    parser.add_option("-q", "--quiet", action="store_true",
                      dest="quiet", help="no warnings")
    parser.add_option("-p", "--show-progress", action="store_true",
                      dest="show_progress", help="show running progress")
    parser.add_option("-v", "--verbose", action="store_true",
                      dest="verbose", help="be verbose")
    parser.add_option("--debug", action="store_true",
                      dest="debug", help="output debugging info")
    parser.add_option("--profile", action="store_true",
                      dest="profile",
                      help="dump Python profiling output to './greylag.prof'")
    (options, args) = parser.parse_args()

    if (len(args) < 1
        or not args[0].endswith('.xml')
        or sum(1 for f in args[1:]
               if not (f.endswith('.ms2') or f.endswith('.ms2.gz') or
                       f.endswith('.ms2.bz2')))
        or (options.part_split and options.part_split < 1)
        or (options.part_merge and options.part_merge < 1)):
        parser.print_help()
        sys.exit(1)

    if sum(1 for x in (options.part_split, options.part, options.part_merge)
           if x != None) > 1:
        error("specify only one of --part-split, --part and --part-merge")

    part_fn_pattern = '%s.0.%sof%s.part.%s'
    part_fn_in_suffix = 'ms2+'
    part_fn_out_suffix = 'out'
    part = None                         # "2of10" -> (2, 10)

    if options.part_split:
        part_infn_pattern = (part_fn_pattern % (options.part_prefix, '%s',
                                                options.part_split,
                                                part_fn_in_suffix))
    elif options.part:
        try:
            part = tuple(int(x) for x in options.part.split('of', 1))
        except:
            parser.print_help()
            sys.exit(1)
        if not 1 <= part[0] <= part[1]:
            error("bad --part parameter")
        part_infn_pattern = (part_fn_pattern % (options.part_prefix, '%s',
                                                part[1], part_fn_in_suffix))
        part_outfn_pattern = (part_fn_pattern % (options.part_prefix, part[0],
                                                 part[1], part_fn_out_suffix))
    elif options.part_merge:
        part_outfn_pattern = (part_fn_pattern % (options.part_prefix, '%s',
                                                 options.part_merge,
                                                 part_fn_out_suffix))
    log_level = logging.WARNING
    if options.quiet:
        log_level = logging.ERROR
    if options.verbose:
        log_level = logging.INFO
    if options.debug:
        log_level = logging.DEBUG
    logging.basicConfig(level=log_level, datefmt='%b %e %H:%M:%S',
                        format='%(asctime)s %(levelname)s: %(message)s')
    info("starting")

    # read params
    parameters = read_xml_parameters(args[0])
    default_parameter_fn = parameters.get("list path, default parameters")
    if default_parameter_fn:
        default_parameters = read_xml_parameters(default_parameter_fn)
        default_parameters.update(parameters)
        parameters = default_parameters
    global XTP
    XTP = validate_parameters(parameters)

    if len(args) > 1:
        spectrum_fns = args[1:]
    else:
        spectrum_fns = [XTP["spectrum, path"]]

    taxonomy = read_taxonomy(XTP["list path, taxonomy information"])

    initialize_spectrum_parameters(options.quirks_mode)
    #greylag_search.CP = CP

    if options.part_split:
        # FIX: clean this up
        info("reading spectrum masses")
        sp_files = [ open(fn) for fn in spectrum_fns ]
        masses = cgreylag.spectrum.read_ms2_spectrum_masses([ f.fileno()
                                                              for f in sp_files ])
        for f in sp_files:
            f.close()
        info("writing %s sets of input files", options.part_split)
        mass_bands = list(generate_mass_bands(options.part_split, masses))
        #info("mass bands: %s", mass_bands)
        mass_band_ubs = [ x[2] for x in mass_bands ]
        mass_band_files = [ open(part_infn_pattern % n, 'w')
                            for n in range(1, options.part_split+1) ]
        mass_band_fds = [ f.fileno() for f in mass_band_files ]
        for n, fn in enumerate(spectrum_fns):
            inf = open(fn)
            cgreylag.spectrum.split_ms2_by_mass_band(inf, mass_band_fds, n,
                                                     mass_band_ubs)
            inf.close()
        for f in mass_band_files:
            f.close()

        info("finished, wrote %s sets of input files", options.part_split)
        logging.shutdown()
        return

    # read spectra
    if options.part:
        # read from a --part-split file
        spectra = list(cgreylag.spectrum.read_spectra_from_ms2(open(part_infn_pattern
                                                                    % part[0]),
                                                               -1))
    else:
        spectra = list(itertools.chain(
            *[ cgreylag.spectrum.read_spectra_from_ms2(open(fn), n)
               for n, fn in enumerate(spectrum_fns) ]))
    spectra.sort(key=lambda x: x.mass)

    if not spectra:
        warning("no input spectra")
    info("read %s spectra (mass range %s - %s)", len(spectra),
         spectra[0].mass, spectra[-1].mass)

    # filter and normalize spectra
    # NI: optionally using "contrast angle" (mprocess::subtract)
    if XTP["spectrum, use conditioning"]:
        spectra = [ s for s in spectra
                    if s.filter_and_normalize(XTP["spectrum, minimum fragment mz"],
                                              XTP["spectrum, dynamic range"],
                                              XTP["spectrum, minimum peaks"],
                                              XTP["spectrum, total peaks"]) ]
    cgreylag.spectrum.set_searchable_spectra(spectra)
    info("     %s spectra after filtering", len(spectra))

    # read sequence dbs
    # [(defline, seq, filename), ...]
    fasta_db = list(read_fasta_files(taxonomy[XTP["protein, taxon"]]))
    # [(idno, offset, defline, seq, seq_filename), ...]
    db = []
    for idno, (defline, sequence, filename) in enumerate(fasta_db):
        db.extend(split_sequence_into_aa_runs(idno, defline, sequence,
                                              filename))
    db_residue_count = sum(len(dbi[3]) for dbi in db)

    # (idno, offset) -> (defline, run_seq, filename)
    db_info = dict(((idno, offset), (defline, seq, filename))
                   for idno, offset, defline, seq, filename in db)

    info("read %s sequences (%s runs, %s residues)", len(fasta_db), len(db),
         db_residue_count)
    if not db:
        error("no database sequences")

    score_statistics = cgreylag.score_stats(len(spectra))

    # (cleavage_re, position of cleavage in cleavage_re)
    cleavage_pattern, cleavage_pos \
                      = cleavage_motif_re(XTP["protein, cleavage site"])
    if cleavage_pos == None:
        error("cleavage site '%s' is missing '|'",
              XTP["protein, cleavage site"])
    cleavage_pattern = re.compile(cleavage_pattern)

    if not options.part_merge:
        if part:
            del spectra                 # try to release memory

        for idno, offset, defline, seq, seq_filename in db:
            if options.show_progress:
                sys.stderr.write("\r%s of %s sequences, %s candidates"
                                 % (idno, len(fasta_db),
                                    score_statistics.candidate_spectrum_count))
            cleavage_points = list(generate_cleavage_points(cleavage_pattern,
                                                            cleavage_pos, seq))
            for begin, end, missed_cleavage_count \
                    in generate_peptides(seq, cleavage_points,
                                         XTP["scoring, maximum missed cleavage sites"]):
                peptide_seq = seq[begin:end]
                #debug('generated peptide: %s', peptide_seq)
                cgreylag.spectrum.search_peptide_all_mods(idno, offset, begin,
                                                          peptide_seq,
                                                          missed_cleavage_count,
                                                          score_statistics)
        if options.show_progress:
            sys.stderr.write("\r%60s\r" % ' ') 

        info('%s candidate spectra examined',
             score_statistics.candidate_spectrum_count)
        filter_matches(score_statistics)

        if part:
            # try to release memory
            del db
            del fasta_db
            cgreylag.spectrum.set_searchable_spectra([])

            partfile = zopen(part_outfn_pattern, 'w')
            cPickle.dump((part, pythonize_swig_object(score_statistics)),
                         partfile, cPickle.HIGHEST_PROTOCOL)
            partfile.close()
            info("finished, part file written to '%s'", part_outfn_pattern)
            logging.shutdown()
            return
    else:
        info('loading/merging %s parts' % options.part_merge)
        part0, score_statistics \
               = cPickle.load(zopen(part_outfn_pattern % 1))
        offset = len(score_statistics.best_score)
        for p in range(2, options.part_merge+1):
            part0, score_statistics0 \
                   = cPickle.load(zopen(part_outfn_pattern % p))
            merge_score_statistics(score_statistics, score_statistics0,
                                   offset)
            offset += len(score_statistics0.best_score)
        if len(score_statistics.best_score) != len(spectra):
            error("error during part merge (expecting %s, got %s)",
                  len(spectra), len(score_statistics.best_score))

    info('processing results')

    #debug('best score: %s', tuple(score_statistics.best_score))

    # filter results and calculate statistics
    (spec_prot_info_items, expect, protein_expect, intensity, survival_curve,
     line_parameters, passing_spectra) \
         = process_results(score_statistics, fasta_db, spectra,
                           db_residue_count)

    info('writing results')

    if options.output:
        if options.output != '-':
            sys.stdout = zopen(options.output, 'w')
    else:
        output_path = XTP["output, path"]
        if output_path:
            sys.stdout = zopen(output_path, 'w')

    print_results_XML(options, XTP, db_info, spectrum_fns,
                      spec_prot_info_items, spectra, expect, protein_expect,
                      intensity, survival_curve, line_parameters,
                      passing_spectra, score_statistics)

    sys.stdout.flush()                  # finish writing before
    info('finished')                    # saying 'finished'
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
            data_fn = "greylag.prof.tmp"
            prof = hotshot.Profile(data_fn)
            prof.runcall(main)
            prof.close()

            sys.stdout = open("greylag.prof", 'w')

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
