#!/usr/bin/env python

'''Analyze mass spectra and assign peptides and proteins.  (This is a partial
   re-implementation of the X!Tandem algorithm, with many extra features.)
'''

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

##############################################################################
##    "Simplicity is prerequisite for reliability" - Edsger W. Dijkstra     ##
##############################################################################


__version__ = "$Id$"


import cPickle
import fileinput
import itertools
import logging
from logging import debug, info, warning
import math
import optparse
import os
from pprint import pprint, pformat
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

# FIX: is this mono or avg?  (possible error here is ~0.0007 amu)
PROTON_MASS =   1.007276
ELECTRON_MASS = 0.000549                # ?

# Reference values from NIST (http://physics.nist.gov/PhysRefData/)
MONOISOTOPIC_ATOMIC_MASS = {
    'H' :  1.00782503214,
    'C' : 12.00000000,
    'N' : 14.00307400529,
    'O' : 15.994914622115,
    'P' : 30.9737615120,
    'S' : 31.9720706912,
    }

ISOTOPIC_ATOMIC_MASS = {
    'N15' : 15.00010889849,
    }

AVERAGE_ATOMIC_MASS = {
    'H' :  1.007947,
    'C' : 12.01078,
    'N' : 14.00672,
    'O' : 15.99943,
    'P' : 30.9737612,
    'S' : 32.0655,
    }

# The xtandem average residue masses are about 0.002 amu higher than those
# calculated directly from the above average atomic masses.  None of the
# chemists consulted knew of any reason why, aside from lack of precision in
# the average atomic mass estimates.  This shouldn't matter very much, as
# fragmentation calculations should all be monoisotopic, and we can always
# widen the parent tolerance window a bit.


def formula_mass(formula, atomic_mass=MONOISOTOPIC_ATOMIC_MASS):
    """Return the mass of formula, using the given mass regime (monoisotopic
    by default)."""
    parts = [ p or '1' for p in re.split(r'([A-Z][a-z]*)', formula)[1:] ]
    # e.g., parts for alanine is ('C', '3', 'H', '5', 'O', '1', 'N', '1')
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

RESIDUES = sorted(RESIDUE_FORMULA.keys())
RESIDUES_W_BRACKETS = RESIDUES + ['[', ']']


def mass_regime_atomic_masses(spec):
    """Given a regime spec like ('MONO', [('N15', 0.9)]), return a map of atom
    names to masses.
    """
    name, isotopes = spec
    assert name in ['MONO', 'AVG'] and len(isotopes) <= 1
    if name == 'MONO':
        r = MONOISOTOPIC_ATOMIC_MASS.copy()
    else:
        r = AVERAGE_ATOMIC_MASS.copy()
    if isotopes:
        iname, prevalence = isotopes[0]
        assert iname == 'N15' and 0 <= prevalence <= 1
        # this is a simplification, but additional accuracy pointless?
        if name == 'MONO':
            r['N'] = ISOTOPIC_ATOMIC_MASS['N15']
        else:
            r['N'] += (ISOTOPIC_ATOMIC_MASS['N15'] - r['N']) * prevalence
    return r


def initialize_spectrum_parameters(mass_regimes, fixed_mod_map, quirks_mode):
    """Initialize parameters known to the spectrum module.
    fixed_mod_map maps, for example, 'M' to (1, 'O', False, 'M', 'oxidation').
    """

    debug('fixed_mod_map: %s', fixed_mod_map)
    # This is the size of vectors that are indexed by residues (A-Z) or
    # special characters ('[]').
    RESIDUE_LIMIT = max(ord(c) for c in 'Z[]') + 1

    # These are currently monoisotopic.  (deuterium pointless?)
    CP.proton_mass = PROTON_MASS
    CP.hydrogen_mass = formula_mass("H")

    for regime_pair in mass_regimes:
        assert len(regime_pair) == 2    # parent and fragment
        debug('rp: %s', regime_pair)
        for n, regime in enumerate(regime_pair):
            atmass = mass_regime_atomic_masses(regime)
            creg = cgreylag.mass_regime_parameters()

            creg.hydroxyl_mass = formula_mass("OH", atmass)
            creg.water_mass = formula_mass("H2O", atmass)
            creg.ammonia_mass = formula_mass("NH3", atmass)

            creg.fixed_residue_mass.resize(RESIDUE_LIMIT)
            for r in RESIDUES_W_BRACKETS:
                m = 0
                if r in RESIDUES:
                    m = formula_mass(RESIDUE_FORMULA[r], atmass)
                rmod = fixed_mod_map.get(r)
                if rmod:
                    if isinstance(rmod[1], str):
                        if rmod[2]:
                            m += rmod[0] * formula_mass(rmod[1])
                        else:
                            m += rmod[0] * formula_mass(rmod[1], atmass)
                    else:
                        m += rmod[0] * rmod[1]
                creg.fixed_residue_mass[ord(r)] = m
            if not n:
                CP.parent_mass_regime.append(creg);
            else:
                CP.fragment_mass_regime.append(creg);
    for r in RESIDUES_W_BRACKETS:
        debug('fixed_mass %s: %s', r,
              [ "%.6f/%.6f"
                % (CP.parent_mass_regime[rn].fixed_residue_mass[ord(r)],
                   CP.fragment_mass_regime[rn].fixed_residue_mass[ord(r)])
                for rn in range(len(mass_regimes)) ])
    
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

    CP.maximum_modification_combinations_searched \
        = XTP["scoring, maximum modification combinations searched"]
    CP.maximum_simultaneous_modifications_searched \
        = XTP["scoring, maximum simultaneous modifications searched"]

    # CP.factorial[n] == (double) n!
    CP.factorial.resize(100, 1.0)
    for n in range(2, len(CP.factorial)):
        CP.factorial[n] = CP.factorial[n-1] * n

    CP.quirks_mode = bool(quirks_mode)
    CP.hyper_score_epsilon_ratio = 0.999 # must be slightly less than 1

    CP.check_all_fragment_charges = XTP["spectrum, check all fragment charges"]

    if CP.quirks_mode and CP.maximum_modification_combinations_searched == 0:
        CP.maximum_modification_combinations_searched = 1 << 12 # FIX?

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
                raise ValueError("invalid cleavage motif pattern"
                                 " (multiple '|'s)")
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


def generate_cleavage_points(cleavage_re, cleavage_pos, sequence):
    """Yields the offsets of the cleavages in sequence.  The endpoints are
    always included, by convention."""
    yield 0
    for m in cleavage_re.finditer(sequence):
        p = m.start() + cleavage_pos
        if 0 < p < len(sequence):
            yield p
    yield len(sequence)


aa_sequence = re.compile(r'[ARNDCQEGHILKMFPSTWYV]+')

def split_sequence_into_aa_runs(idno, defline, sequence, filename):
    """Returns a tuple (idno, start, defline, seq, filename) for each
    contiguous run of residues in sequence, where 'start' is the position of
    'seq' in 'sequence'."""
    matches = list(aa_sequence.finditer(sequence))
    return [ (idno, m.start(), defline, m.group(), filename)
             for n, m in enumerate(matches) ]


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
    # try to catch misspellings/etc in the config file
    bad_nodes = [ e for e in root.findall('*')
                  if (e.tag != 'note'
                      or set(e.keys()) not in [ set(), set(['type']),
                                                set(['type', 'label']) ]
                      or e.get("type") not in [ None, "input", "description" ]
                      or e.get("type") == "input" and e.get("label") == None
                      or (e.get("type") == "description"
                          and e.get("label") != None)) ]
    if bad_nodes:
        for bn in bad_nodes:
            warning("%s: invalid '%s' element (type=%s, label=%s)"
                    % (filename, bn.tag, bn.get('type'), bn.get('label')))
    # return just the good nodes
    return dict((e.get("label"), e.text)
                for e in root.findall('note')
                if e.get("type") == "input")


# XML parameter file processing

# FIX: This parsing is way too complex.  How to simplify?

def mass_regime_part(part_specification):
    """Parse a single mass regime specification part (e.g., 'MONO(N15@90%)').
    """
    ps = [ x.strip() for x in part_specification.split('(', 1) ]
    if ps[0] not in ('MONO', 'AVG'):
        raise ValueError("invalid mass regime list specification"
                         " (regime id must be 'MONO' or 'AVG')")
    if len(ps) == 1:
        return (ps[0], [])
    if ps[1][-1] != ')':
        raise ValueError("invalid mass regime list specification"
                         " (expected ')')")
    pps = [ x.strip() for x in ps[1][:-1].split(',') ]
    if len(pps) > 1:
        raise ValueError("invalid mass regime list specification"
                         " (multiple isotopes not yet implemented)")
    ppps = [ x.strip() for x in pps[0].split('@', 1) ]
    if len(ppps) != 2:
        raise ValueError("invalid mass regime list specification"
                         " (expected '@')")
    if ppps[0] not in ('N15',):
        raise ValueError("invalid mass regime list specification"
                         " (isotope id must currently be 'N15')")
    if ppps[1][-1] != '%':
        raise ValueError("invalid mass regime list specification"
                         " (expected '%')")
    prevalence = float(ppps[1][:-1]) / 100
    if not (0 <= prevalence <= 1):
        raise ValueError("invalid mass regime list specification"
                         " (prevalence must be in range 0-100%)")
    return (ps[0], [(ppps[0], prevalence)])

def mass_regime_list(mass_regime_list_specification):
    """Check and return a list of regime tuples (parent_regime,
    fragment_regime), where each regime is a tuple (id, [(isotope_id,
    prevalence), ...]).  So, for example, 'AVG/MONO;MONO;MONO(N15@90%)' would
    return

    [[('AVG', []), ('MONO', [])],
     [('MONO', []), ('MONO', [])],
     [('MONO', [('N15', 0.90000000000000002)]),
      ('MONO', [('N15', 0.90000000000000002)])]]

    Multiple isotopes (when implemented) would be comma-separated.
    """
    result = []
    for regspec in mass_regime_list_specification.split(';'):
        halves = [ x.strip() for x in regspec.split('/') ]
        if len(halves) > 2:
            raise ValueError("invalid mass regime list specification"
                             " (too many '/'?)")
        pr = [ mass_regime_part(h) for h in halves ]
        if len(pr) == 1:
            pr = [ pr[0], pr[0] ]
        result.append(pr)
    debug("mass regime list:\n%s", pformat(result))
    return result


def parse_mod_term(s, is_potential=False):
    """Parse a modification term, returning a tuple (sign, mod, fixed_regime,
    residues, description).  For example:
    
    '-C2H3ON!@C' --> (-1, 'C2H3ON', True, 'C', None)
    '42@STY phosphorylation' --> (1, 42.0, False, 'STY', 'phosphorylation')

    """
    
    m = re.match(r'^\s*(-|\+)?(([1-9][0-9.]*)|([A-Z][A-Z0-9]*))(!?)'
                 r'@([A-Z]+|\[|\])(\s+([A-Za-z0-9_]+))?\s*$', s)
    if not m:
        raise ValueError("invalid modification term specification"
                         " '%s'" % s)
    mg = m.groups()
    invalid_residues = set(mg[5]) - set(RESIDUES_W_BRACKETS)
    if invalid_residues:
        raise ValueError("invalid modification list specification"
                         " (invalid residues %s)" % list(invalid_residues))
    delta = mg[1]
    if mg[2]:
        delta = float(mg[1])
        if is_potential and abs(delta) < 0.0001:
            raise ValueError("invalid modification list specification"
                             " (delta '%s' is too small)" % delta)
    residues = mg[5]
    if not is_potential and len(residues) != 1:
        raise ValueError("invalid modification list specification '%s' (only"
                         " potential modifications may have multiple residues)"
                         % residues)
    if len(residues) != len(set(residues)):
        raise ValueError("invalid modification list specification"
                         " '%s' (duplicate residues prohibited)"
                         % residues)
    return (mg[0] == '-' and -1 or 1, delta, mg[4] == '!', residues, mg[7])


def fixed_mod_list(specification):
    """Check and return a list of modification tuples."""
    if not specification.strip():
        return []
    result = [ parse_mod_term(s) for s in specification.split(',') ]
    residues = [ x[3] for x in result ]
    # this check is across terms; the one above is within terms
    if len(residues) != len(set(residues)):
        raise ValueError("invalid modification list specification"
                         " '%s' (duplicate residues prohibited)"
                         % residues)
    debug("fixed_mod_list:\n%s", pformat(result))
    return result


def parse_mod_basic_expression(s):
    s = s.strip()
    if s[0] == '(':
        tree, rest =  parse_mod_disjunction(s[1:])
        rest = rest.lstrip()
        if rest[:1] != ')':
            raise ValueError("invalid modification list specification"
                             " (expected matching ')')")
        return tree, rest[1:]
    parts = re.split(r'([;,()])', s, 1)
    if len(parts) == 1:
        term, rest = s, ''
    else:
        assert len(parts) == 3
        term, rest = parts[0], parts[1]+parts[2]
    return parse_mod_term(term, is_potential=True), rest

def parse_mod_conjunction(s):
    result = []
    while s:
        tree, s = parse_mod_basic_expression(s)
        result.append(tree)
        s = s.lstrip()
        if s[:1] != ',':
            break
        s = s[1:]
    return result, s

def parse_mod_disjunction(s):
    result = []
    while s:
        tree, s = parse_mod_conjunction(s)
        result.append(tree)
        s = s.lstrip()
        if s[:1] != ';':
            break
        s = s[1:]
    return result, s

def potential_mod_list(specification):
    """Check and return a tree of potential modification tuples.  Nodes at
    even (odd) levels are disjunctions (conjunctions).  (The top list is a
    disjunction.)
    """
    specification = specification.strip()
    if not specification:
        return []
    tree, remainder = parse_mod_disjunction(specification)
    if remainder:
        raise ValueError("invalid modification list specification"
                         " (unexpected '%s')" % remainder)
    debug("potential_mod_list:\n%s", pformat(tree))
    return tree


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
    "refine, modification mass" : (fixed_mod_list, "", p_ni_empty),
    "refine, point mutations" : (bool, "no", p_ni_equal(False)),
    "refine, potential C-terminus modifications": (potential_mod_list, "", p_ni_empty),
    "refine, potential N-terminus modifications": (potential_mod_list, "", p_ni_empty),
    "refine, potential N-terminus modification position limit" : (int, 50), # nyi
    "refine, potential modification mass" : (potential_mod_list, "", p_ni_empty),
    "refine, potential modification motif" : (str, "", p_ni_empty),
    "refine, sequence path" : (str, "", p_ni_empty),
    "refine, spectrum synthesis" : (bool, "no"),
    "refine, tic percent" : (float, 20.0, p_nonnegative), # ignored
    "refine, unanticipated cleavage" : (bool, p_ni_equal(False)),
    "refine, use potential modifications for full refinement" : (bool, "no", p_ni_equal(False)),
    "residue, mass regimes" : (mass_regime_list, "MONO"),
    "residue, modification mass" : (fixed_mod_list, ""),
    "residue, potential modification mass" : (potential_mod_list, ""),
    "residue, potential modification motif" : (str, "", p_ni_empty),
    "scoring, a ions" : (bool, "no", p_ni_equal(False)),
    "scoring, b ions" : (bool, "yes", p_ni_equal(True)),
    "scoring, c ions" : (bool, "no", p_ni_equal(False)),
    "scoring, cyclic permutation" : (bool, "no", p_ni_equal(False)),
    "scoring, include reverse" : (bool, "no", p_ni_equal(False)),
    "scoring, maximum missed cleavage sites" : (int, None, p_nonnegative),
    "scoring, maximum modification combinations searched" : (int, 0, p_nonnegative), # 0 -> no limit # FIX
    "scoring, maximum simultaneous modifications searched" : (int, None, p_nonnegative),
    "scoring, minimum ion count" : (int, None, p_positive),
    "scoring, pluggable scoring" : (bool, "no"), # ignored
    "scoring, x ions" : (bool, "no", p_ni_equal(False)),
    "scoring, y ions" : (bool, "yes", p_ni_equal(True)),
    "scoring, z ions" : (bool, "no", p_ni_equal(False)),
    "spectrum, check all charges" : (bool, "no", p_ni_equal(False)),
    "spectrum, check all fragment charges" : (bool, "no"),
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
    "spectrum, path" : (str, ""),
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
                error("parameter '%s' value '%s' not in %s (feature not"
                      " implemented?)" % (p_name, v, type_))
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
            error("parameter '%s' has invalid value '%s' (or feature not"
                  " implemented)" % (p_name, v))
        pmap[p_name] = v

    unknown_parameters = set(parameters) - set(XML_PARAMETER_INFO)
    if unknown_parameters:
        warning("%s unknown parameters:\n %s"
                % (len(unknown_parameters),
                   pformat(sorted(list(unknown_parameters)))))

    # FIX: this would trip the "cyclic" param (assume Daltons)
    assert (abs(pmap["spectrum, parent monoisotopic mass error plus"]) > 0.095
            and (abs(pmap["spectrum, parent monoisotopic mass error minus"])
                 > 0.095)), "feature not implemented (cyclic param)"

    # FIX: where should this go?
    if pmap["spectrum, fragment mass error"] > 0.5:
        warning("'spectrum, fragment mass error' is %s",
                pmap["spectrum, fragment mass error"])

    return pmap


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


def enumerate_conjunction(mod_tree, limit, conjuncts=[], empty=True):
    if not mod_tree:
        if not empty and len(conjuncts) <= limit:
            yield conjuncts
        return
    first, rest = mod_tree[0], mod_tree[1:]
    if isinstance(first, list):
        for x in enumerate_disjunction(first, limit):
            if x:
                empty = False
            for y in enumerate_conjunction(rest, limit, conjuncts + x, empty):
                yield y
    else:
        for y in enumerate_conjunction(rest, limit, conjuncts, empty):
            yield y
        for y in enumerate_conjunction(rest, limit, conjuncts + [first],
                                       False):
            yield y

def enumerate_disjunction(mod_tree, limit):
    """Generates the conjuncts for mod_tree that are no longer than limit."""
    assert isinstance(mod_tree, list)
    yield []
    for b in mod_tree:
        for s in enumerate_conjunction(b, limit):
            yield s

def get_mod_conjunct_info(mod_tree, limit):
    def has_N_mod(c):
        return bool(sum(1 for t in c if t[3] == '['))
    
    return [ (conjunct, has_N_mod(conjunct))
             for conjunct in enumerate_disjunction(mod_tree, limit) ]


def gen_delta_bag(i, remainder, bag):
    if i < 1:
        assert i == 0
        bag[0] = remainder
        yield tuple(bag)
        return
    for delta in range(1, remainder-i+1):
        bag[i] = delta
        for x in gen_delta_bag(i-1, remainder-delta, bag):
            yield x

def generate_delta_bags(mod_count, conjunct_length):
    """Generate all tuples of positive integers having length conjunct_length
    and sum mod_count.  As a special case, () is such a tuple having length 0
    and sum 0.
    """
    if conjunct_length == 0:
        return mod_count == 0 and [()] or []
    if mod_count < conjunct_length:
        return []
    return gen_delta_bag(conjunct_length - 1, mod_count, [0] * conjunct_length)


def search_all(options, fasta_db, db, cleavage_pattern, cleavage_pos,
               score_statistics):
    """Search sequence database against searchable spectra."""
    warning("assuming no N-term mods (???)")

    mod_limit = XTP["scoring, maximum simultaneous modifications searched"]
    combination_limit \
        = XTP["scoring, maximum modification combinations searched"]

    mod_conjunct_info = get_mod_conjunct_info(
        XTP["residue, potential modification mass"], mod_limit)
    info("%s potential modification conjuncts to search",
         len(mod_conjunct_info))
    debug("mod_conjunct_info:\n%s", pformat(mod_conjunct_info))

    min_peptide_length = 5
    for idno, offset, defline, seq, seq_filename in db:
        if options.show_progress:
            sys.stderr.write("\r%s of %s sequences, %s candidates"
                             % (idno, len(fasta_db),
                                score_statistics.candidate_spectrum_count))
        cleavage_points = list(generate_cleavage_points(cleavage_pattern,
                                                        cleavage_pos, seq))

        for mr_index, mass_regime in enumerate(XTP["residue, mass regimes"]):
            score_statistics.combinations_searched = 0
            for mod_count in range(mod_limit + 1):
                for has_pca in (False, True):
                    for mod_conjunct, has_N_mod in mod_conjunct_info:
                        if has_pca and has_N_mod:
                            continue    # mutually exclusive, for now
                        debug("mod_count: %s", mod_count)
                        debug("mod_conjunct: %s", mod_conjunct)
                        for delta_bag in generate_delta_bags(mod_count,
                                                             len(mod_conjunct)):
                            debug("delta_bag: %s", delta_bag)
                            if (combination_limit
                                and (combination_limit
                                     < score_statistics.combinations_searched)):
                                break

                            base_N_mass = 0
                            base_C_mass = 0
                            delta_mass = 0
#                             cgreylag.spectrum.search_run(
#                                 XTP["scoring, maximum missed cleavage sites"],
#                                 min_peptide_length, idno, offset, seq,
#                                 cleavage_points, mr_index, score_statistics)




    if options.show_progress:
        sys.stderr.write("\r%60s\r" % ' ')

    info('%s candidate spectra examined',
         score_statistics.candidate_spectrum_count)


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
    info("%s spectra with valid models", valid_spectra_count)
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
                    item[1][0][1][0].peptide_begin, # domain 0 starting
                                                    # position (ouch)
                    expect[item[0]])    # spectrum expect
        spec_prot_info_items.sort(key=protein_order_less_than)
    else:
        # sort each spectra's proteins by (protein expect, protein id)
        assert False

    return (spec_prot_info_items, expect, protein_expect, intensity,
            survival_curve, line_parameters, passing_spectra)


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


def print_results_XML(options, db_info, spectrum_fns, spec_prot_info_items,
                      spectra, expect, protein_expect, intensity,
                      survival_curve, line_parameters, passing_spectra,
                      score_statistics):
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
        defline, run_seq, seq_filename = db_info[(si0.sequence_index,
                                                  si0.sequence_offset)]

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
                d0_defline, d0_run_seq, d0_seq_filename \
                            = db_info[(domains[0].sequence_index,
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
                                 = db_info[(dom.sequence_index,
                                            dom.sequence_offset)]

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
                              dom.peptide_begin+len(dom.peptide_sequence),
                              expect[spectrum_id], delta_precision,
                              dom.peptide_mass, delta_precision, delta,
                              cgreylag.scale_hyperscore(score_statistics.best_score[spectrum_id]),
                              cgreylag.scale_hyperscore(score_statistics.second_best_score[spectrum_id]),
                              cgreylag.scale_hyperscore(dom.ion_scores[cgreylag.ION_Y]),
                              dom.ion_peaks[cgreylag.ION_Y],
                              cgreylag.scale_hyperscore(dom.ion_scores[cgreylag.ION_B]),
                              dom.ion_peaks[cgreylag.ION_B],
                              get_prefix_sequence(dom.peptide_begin,
                                                  dom.sequence_offset,
                                                  dom_run_seq),
                              get_suffix_sequence((dom.peptide_begin
                                                   + len(dom.peptide_sequence)),
                                                  dom.sequence_offset,
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
            print ('<GAML:attribute type="charge">%s</GAML:attribute>'
                   % sp.charge)
            print ('<GAML:Xdata label="%s.spectrum" units="MASSTOCHARGERATIO">'
                   % (sp.id))
            print ('<GAML:values byteorder="INTEL" format="ASCII"'
                   ' numvalues="%s">' % len(sp.peaks))
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
            print ('<GAML:values byteorder="INTEL" format="ASCII"'
                   ' numvalues="%s">' % len(sp.peaks))
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
            print '	<note type="input"' \
                  ' label="xtandem quirks mode">yes</note>'
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
    pa = parser.add_option
    pa("-o", "--output", dest="output", help="destination file [default as"
       " given in parameter file, '-' for stdout]", metavar="FILE")
    pa("--quirks-mode", action="store_true", dest="quirks_mode",
       help="try to generate results as close as possible to those of X!Tandem"
       " (possibly at the expense of accuracy)")
    pa("--part-split", dest="part_split", type="int", help="split input into M"
       " parts, to prepare for --part runs [NOTE: the same parameter file and"
       " same spectrum files (in the same order) must be specified for all"
       " --part* steps]", metavar="M") 
    pa("--part", dest="part", help="search one part, previously created with"
       " --part-split; e.g. '1of2' and '2of2'", metavar="NofM")
    pa("--part-merge", dest="part_merge", type="int", help="merge the"
       " previously searched M results and continue", metavar="M")
    default_prefix = 'greylag'
    pa("--part-prefix", dest="part_prefix", default=default_prefix,
       help="prefix to use for temporary part files [default='%s']"
       % default_prefix, metavar="PREFIX") 
    pa("--compress-level", dest="compress_level", type="int",
       help="compression level to use for compressed files created [default=1"
       " for *.gz, 9 for *.bz2]", metavar="N")
    pa("-q", "--quiet", action="store_true", dest="quiet", help="no warnings")
    pa("-p", "--show-progress", action="store_true", dest="show_progress",
    help="show running progress")
    pa("-v", "--verbose", action="store_true", dest="verbose",
       help="be verbose") 
    pa("--copyright", action="store_true", dest="copyright",
       help="print copyright and exit")
    pa("--debug", action="store_true", dest="debug",
       help="output debugging info")
    pa("--profile", action="store_true", dest="profile",
       help="dump Python profiling output to './greylag.prof'") 
    (options, args) = parser.parse_args()

    if options.copyright:
        print __copyright__
        sys.exit(0)

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
        if not spectrum_fns:
            error("input spectrum files not specified")

    taxonomy = read_taxonomy(XTP["list path, taxonomy information"])

    fixed_mod_map = dict((r[3], r) for r in XTP["residue, modification mass"])
    initialize_spectrum_parameters(XTP["residue, mass regimes"], fixed_mod_map,
                                   options.quirks_mode)

    if options.part_split:
        # FIX: clean this up
        info("reading spectrum masses")
        sp_files = [ open(fn) for fn in spectrum_fns ]
        masses = cgreylag.spectrum.read_ms2_spectrum_masses([ f.fileno()
                                                              for f
                                                              in sp_files ])
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
        spectra = cgreylag.spectrum.read_spectra_from_ms2(open(part_infn_pattern
                                                               % part[0]),
                                                          -1)
    else:
        spectra = itertools.chain(
            *[ cgreylag.spectrum.read_spectra_from_ms2(open(fn), n)
               for n, fn in enumerate(spectrum_fns) ])
    spectra = list(spectra)
    spectra.sort(key=lambda x: x.mass)

    if not spectra:
        warning("no input spectra")
    else:
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

    # (cleavage_re, position of cleavage in cleavage_re)
    cleavage_pattern, cleavage_pos \
                      = cleavage_motif_re(XTP["protein, cleavage site"])
    if cleavage_pos == None:
        error("cleavage site '%s' is missing '|'",
              XTP["protein, cleavage site"])
    cleavage_pattern = re.compile(cleavage_pattern)

    if not options.part_merge:
        cgreylag.spectrum.set_searchable_spectra(spectra)
        score_statistics = cgreylag.score_stats(len(spectra))

        if part:
            del spectra                 # try to release memory

        search_all(options, fasta_db, db, cleavage_pattern, cleavage_pos,
                   score_statistics)

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
        info('%s candidate spectra were examined',
             score_statistics.candidate_spectrum_count)

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

    print_results_XML(options, db_info, spectrum_fns, spec_prot_info_items,
                      spectra, expect, protein_expect, intensity,
                      survival_curve, line_parameters, passing_spectra,
                      score_statistics)

    # finish writing before saying 'finished' (this should be a flush(), but
    # BZFile doesn't have that method)
    sys.stdout.close()
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
        logging.shutdown()
        sys.exit(1)


# FIXES:
# - need to rigorously check for bad input in any file (fasta, spectrum, xml)
# - escape special XML chars
