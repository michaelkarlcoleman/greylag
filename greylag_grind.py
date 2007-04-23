#!/usr/bin/env python

'''This program does the actual work to search mass spectra against a sequence
database.  <job-id> is a unique identifier used as a prefix for work and
output directories.  The <configuration-file> contains program options.  The
spectra are in <ms2-file>s; greylag-index-spectra must already have been run
on them.

This program can operate in standalone mode (useful for testing) or command
mode.  If the --work-slice option is given, that slice will be searched in
standalone mode, then the program will exit.  Otherwise, the program will
process commands until an exit command is received.

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

##############################################################################
##    "Simplicity is prerequisite for reliability" - Edsger W. Dijkstra     ##
##############################################################################


__version__ = "0.0"


import ConfigParser
from collections import defaultdict
import contextlib
import cPickle
import fileinput
import gzip
import itertools
import logging
from logging import debug, info, warning
import math                             #??
import optparse
import os
from pprint import pprint, pformat
import re
from socket import gethostname
import sys

import cgreylag


# Try to drop dead immediately on SIGINT (control-C), instead of normal Python
# KeyboardInterrupt processing, since we may spend long periods of time in
# uninterruptible C++ calls.  Also die immediately on SIGPIPE.
try:
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
except:
    pass


def error(s, *args):
    "fatal error"
    logging.error(s, *args)
    if __name__ != "__main__":
        raise Exception("(fatal error, but unit testing, so not exiting)")
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

# most prevalent only (1 in 1000)
ISOTOPIC_ATOMIC_MASS = {                # prevalence (in %)
    'C13' : 13.003354837810,            # 1.078
    'N15' : 15.00010889849,             # 0.3687
    'O18' : 17.99916049,                # 0.20514
    'S33' : 32.9714585012,              # 0.762
    'S34' : 33.9678668311,              # 4.2928
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
    by default).

    >>> formula_mass('H2O', { 'H':1, 'O':16 })
    18
    >>> # monoisotopic mass of glycine
    >>> str(round(formula_mass('C2H3ON'), 4))
    '57.0215'

    """
    parts = [ p or '1' for p in re.split(r'([A-Z][a-z]*)', formula)[1:] ]
    # parts for glycine = ['C', '2', 'H', '3', 'O', '1', 'N', '1']
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


# accessed as a global variable (FIX?)
# [0][1] -> 'H' -> fragment mass of H for regime 0
MASS_REGIME_ATOMIC_MASSES = []


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


def initialize_spectrum_parameters(options, mass_regimes, fixed_mod_map):
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

    global MASS_REGIME_ATOMIC_MASSES
    for regime_pair in mass_regimes:
        assert len(regime_pair) == 2    # parent and fragment
        info('mass regime: %s', regime_pair)
        MASS_REGIME_ATOMIC_MASSES.append([])
        for n, regime in enumerate(regime_pair):
            atmass = mass_regime_atomic_masses(regime)
            MASS_REGIME_ATOMIC_MASSES[-1].append(atmass)
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
            # assuming these are monoisotopic (not regime)
            creg.fixed_residue_mass[ord('[')] += formula_mass("H")
            creg.fixed_residue_mass[ord(']')] += formula_mass("OH")
            if not n:
                CP.parent_mass_regime.append(creg);
            else:
                creg.fixed_residue_mass[ord('[')] -= CP.hydrogen_mass
                creg.fixed_residue_mass[ord(']')] -= creg.hydroxyl_mass
                CP.fragment_mass_regime.append(creg);
    for r in RESIDUES_W_BRACKETS:
        info('fixed_mass %s: %s', r,
             [ "%.6f/%.6f"
               % (CP.parent_mass_regime[rn].fixed_residue_mass[ord(r)],
                  CP.fragment_mass_regime[rn].fixed_residue_mass[ord(r)])
               for rn in range(len(mass_regimes)) ])
    for r in RESIDUES:
        for rn in range(len(mass_regimes)):
            # physically impossible (and results would be garbage)
            if CP.parent_mass_regime[rn].fixed_residue_mass[ord(r)] < 1.0:
                raise ValueError('bogus parent mass specification for %s' % r)
            if CP.fragment_mass_regime[rn].fixed_residue_mass[ord(r)] < 1.0:
                raise ValueError('bogus parent mass specification for %s' % r)

    CP.parent_mass_tolerance_1 = XTP["parent_mz_tolerance"]
    CP.parent_mass_tolerance_max = (XTP["parent_mz_tolerance"]
                                    * XTP["charge_limit"])

    CP.fragment_mass_tolerance = XTP["fragment_mass_tolerance"]
    CP.intensity_class_count = XTP["intensity_class_count"]

    CP.minimum_peptide_length = XTP["min_peptide_length"]

    # CP.ln_factorial[n] == ln(n!)
    CP.ln_factorial.resize(int(XTP["max_parent_spectrum_mass"]
                               / XTP["fragment_mass_tolerance"] + 100), 0.0)
    for n in range(2, len(CP.ln_factorial)):
        CP.ln_factorial[n] = CP.ln_factorial[n-1] + math.log(n)

    CP.estimate_only = bool(options.estimate_only)
    CP.show_progress = bool(options.show_progress)


def cleavage_motif_re(motif):
    """Return (regexp, pos), where regexp is a regular expression that will
    match a cleavage motif, and pos is the position of the cleavage with
    respect to the match (or None).  (The RE actually matches one character,
    the rest matching as lookahead, so that re.finditer will find all
    overlapping matches.)

    >>> cleavage_motif_re('[KR]|{P}')
    ('[KR](?=[^P])', 1)
    >>> cleavage_motif_re('[X]|[X]')
    ('.(?=.)', 1)

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
    always included, by convention.

    >>> list(generate_cleavage_points(re.compile('[KR](?=[^P])'), 1,
    ...                               'ARBCKDEKPF'))
    [0, 2, 5, 10]

    """
    yield 0
    for m in cleavage_re.finditer(sequence):
        p = m.start() + cleavage_pos
        if 0 < p < len(sequence):
            yield p
    yield len(sequence)


AA_SEQUENCE = re.compile(r'[ARNDCQEGHILKMFPSTWYV]+')

def split_sequence_into_aa_runs(idno, locusname, defline, sequence, filename):
    """Returns a tuple (idno, start, locusname, defline, seq, filename) for
    each contiguous run of residues in sequence, where 'start' is the position
    of 'seq' in 'sequence'.

    >>> pprint(split_sequence_into_aa_runs(123, 'ln', 'ln defline',
    ...                                    'STSS*DEFABA', 'filename'))
    [(123, 0, 'ln', 'ln defline', 'STSS', 'filename'),
     (123, 5, 'ln', 'ln defline', 'DEFA', 'filename'),
     (123, 10, 'ln', 'ln defline', 'A', 'filename')]

    """
    return [ (idno, m.start(), locusname, defline, m.group(), filename)
             for n, m in enumerate(AA_SEQUENCE.finditer(sequence)) ]


def read_fasta_files(filenames):
    """Yield (locusname, defline, sequence, filename) tuples as read from
    FASTA files (uppercasing sequence)."""

    loci_seen = set()
    locusname, defline = None, None
    seqs = []
    for line in fileinput.input(filenames):
        line = line.strip()
        if line[:1] == '>':
            if defline != None:
                yield (locusname, defline, ''.join(seqs), fileinput.filename())
            elif seqs:
                fileerror("bad format: line precedes initial defline")
            defline = line[1:]
            locusname = defline.split(None, 1)[0]
            if locusname in loci_seen:
                error("locus name '%s' is not unique in the search database(s)"
                      % locusname)
            loci_seen.add(locusname)
            seqs = []
        else:
            seqs.append(line.upper())
    if defline:
        yield (locusname, defline, ''.join(seqs), fileinput.filename())


def read_spectra_slice(spectrum_fns, offset_indices, slice):
    s_l, s_u = slice
    assert 0 <= s_l <= s_u <= 1

    total_spectra = sum(len(oi) for oi in offset_indices)
    sp_l = int(float(s_l) * total_spectra)
    sp_u = int(float(s_u) * total_spectra)

    seeing_l = 0
    spectra = []

    for f_no, (fn, oi) in enumerate(zip(spectrum_fns, offset_indices)):
        seeing_u = seeing_l + len(oi)
        f_l = max(sp_l, seeing_l)
        f_u = min(sp_u, seeing_u)

        if f_l < f_u:
            b_l = oi[f_l - seeing_l]
            b_u = oi[f_u - seeing_l] if f_u - seeing_l < len(oi) else -1 # oo

            with open(fn) as f:
                fsp = cgreylag.spectrum.read_spectra_from_ms2(f, f_no,
                                                              b_l, b_u)
                spectra.extend(fsp)
        seeing_l = seeing_u

    return spectra


# XML parameter file processing

# FIX: This parsing is way too complex.  How to simplify?

def mass_regime_part(part_specification):
    """Parse a single mass regime specification part.

    >>> mass_regime_part('MONO')
    ('MONO', [])
    >>> mass_regime_part('AVG')
    ('AVG', [])
    >>> mass_regime_part('MONO(N15@87.5%)')
    ('MONO', [('N15', 0.875)])

    """
    ps = [ x.strip() for x in part_specification.partition('(') ]
    if ps[0] not in ('MONO', 'AVG'):
        raise ValueError("invalid mass regime list specification"
                         " (regime id must be 'MONO' or 'AVG')")
    if not ps[1]:
        return (ps[0], [])
    if ps[2][-1] != ')':
        raise ValueError("invalid mass regime list specification"
                         " (expected ')')")
    pps = [ x.strip() for x in ps[2][:-1].split(',') ]
    if len(pps) > 1:
        raise ValueError("invalid mass regime list specification"
                         " (multiple isotopes not yet implemented)")
    ppps = [ x.strip() for x in pps[0].partition('@') ]
    if not ppps[1]:
        raise ValueError("invalid mass regime list specification"
                         " (expected '@')")
    if ppps[0] not in ('N15',):
        raise ValueError("invalid mass regime list specification"
                         " (isotope id must currently be 'N15')")
    if ppps[2][-1] != '%':
        raise ValueError("invalid mass regime list specification"
                         " (expected '%')")
    prevalence = float(ppps[2][:-1]) / 100
    if not (0 <= prevalence <= 1):
        raise ValueError("invalid mass regime list specification"
                         " (prevalence must be in range 0-100%)")
    return (ps[0], [(ppps[0], prevalence)])

def mass_regime_list(mass_regime_list_specification):
    """Check and return a list of regime tuples (parent_regime,
    fragment_regime), where each regime is a tuple (id, [(isotope_id,
    prevalence), ...]). Multiple isotopes (when implemented) would be
    comma-separated.

    >>> pprint(mass_regime_list('AVG/MONO;MONO;MONO(N15@75%)'))
    [[('AVG', []), ('MONO', [])],
     [('MONO', []), ('MONO', [])],
     [('MONO', [('N15', 0.75)]), ('MONO', [('N15', 0.75)])]]

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

    # The first fragmentation regime should generally be MONO, so that
    # formulaic deltas with '!' do the expected thing.
    if result[0][1] != ('MONO', []):
        raise ValueError("first fragmentation regime was something other than"
                         " 'MONO' with no isotopes--this is almost certainly"
                         " not what was intended")
    return result


def parse_mod_term(s, is_potential=False):
    """Parse a modification term, returning a tuple (sign, mod, fixed_regime,
    residues, description).

    >>> parse_mod_term('-C2H3ON!@C')
    (-1, 'C2H3ON', True, 'C', None)
    >>> parse_mod_term('42@STY phosphorylation', is_potential=True)
    (1, 42.0, False, 'STY', 'phosphorylation')

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
    """Check and return a list of modification tuples.

    >>> fixed_mod_list('57@C')
    [(1, 57.0, False, 'C', None)]
    >>> fixed_mod_list('57@C,CH!@N desc')
    [(1, 57.0, False, 'C', None), (1, 'CH', True, 'N', 'desc')]
    >>> fixed_mod_list('57@C,58@C')
    Traceback (most recent call last):
        ...
    ValueError: invalid modification list specification '['C', 'C']' (duplicate residues prohibited)

    """
    if not specification:
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

    >>> potential_mod_list('')
    []
    >>> potential_mod_list('PO3H@STY; C2H2O@KST')
    [[(1, 'PO3H', False, 'STY', None)], [(1, 'C2H2O', False, 'KST', None)]]
    >>> potential_mod_list('PO3H@STY, C2H2O@KST')
    [[(1, 'PO3H', False, 'STY', None), (1, 'C2H2O', False, 'KST', None)]]
    >>> pprint(potential_mod_list('''(PO3H@STY phosphorylation;
    ...                               C2H2O@KST acetylation;
    ...                               CH2@AKST methylation),
    ...                              O@M oxidation'''))
    [[[[(1, 'PO3H', False, 'STY', 'phosphorylation')],
       [(1, 'C2H2O', False, 'KST', 'acetylation')],
       [(1, 'CH2', False, 'AKST', 'methylation')]],
      (1, 'O', False, 'M', 'oxidation')]]

    """
    if not specification:
        return []
    tree, remainder = parse_mod_disjunction(specification)
    if remainder:
        raise ValueError("invalid modification list specification"
                         " (unexpected '%s')" % remainder)
    debug("potential_mod_list:\n%s", pformat(tree))
    return tree


def p_positive(x): return x > 0
def p_negative(x): return x < 0
def p_nonnegative(x): return x >= 0
def p_proportion(x): return 0 <= x <= 1

# Configuration file parameter specification
#
# name -> (type, default, check_fn)
# type may be bool, int, float, str, modlist, or a tuple of values
# default, as a string value, or None if value must be explicitly specified
# check_fn is an optional function that returns True iff the value is valid

PARAMETER_INFO = {
    "databases" : (str, None),
    "decoy_locus_prefix" : (str, "SHUFFLED_"),
    "mass_regimes" : (mass_regime_list, "MONO"),
    "pervasive_mods" : (fixed_mod_list, ""),
    "potential_mods" : (potential_mod_list, ""),
    "potential_mod_limit" : (int, 2, p_nonnegative),
    "charge_limit" : (int, 3, p_positive),
    "min_peptide_length" : (int, 5, p_positive), # needed?
    "min_parent_spectrum_mass" : (float, 0, p_nonnegative),
    "max_parent_spectrum_mass" : (float, 10000, p_nonnegative),
    "TIC_cutoff_proportion" : (float, 0.98, p_proportion),
    "parent_mz_tolerance" : (float, 1.25, p_nonnegative),
    "fragment_mass_tolerance" : (float, 0.5, p_nonnegative),
    "intensity_class_count" : (int, 3, p_positive),
    "intensity_class_ratio" : (float, 2.0, p_positive), # really > 1.0?
    "best_result_count" : (int, 5, p_positive),
    }


def validate_parameters(parameters, parameter_info=PARAMETER_INFO):
    """Verify that parameters are valid, have valid values, and correspond to
    currently implemented functionality.  Values are converted, default values
    are filled in, and the resulting name/value dict returned.

    >>> sorted(validate_parameters({'spectrum, total peaks' : '40'},
    ...                            {'spectrum, total peaks'
    ...                             : (int, '50', p_positive),
    ...                             'output, spectra'
    ...                             : (bool, 'no')}).items())
    [('output, spectra', False), ('spectrum, total peaks', 40)]

    """

    pmap = {}
    for p_name, p_info in sorted(parameter_info.items()):
        type_ = p_info[0]
        default = p_info[1]
        check_fn = len(p_info) > 2 and p_info[2] or None

        v = parameters.get(p_name)
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

    unknown_parameters = set(parameters) - set(parameter_info)
    if unknown_parameters:
        warning("%s unknown parameter(s):\n %s"
                % (len(unknown_parameters),
                   pformat(sorted(list(unknown_parameters)))))
    return pmap


def generate_mass_bands(band_count, mass_list):
    """Yield (n, mass_lb, mass_ub) for each mass band, where n ranges from 1
    to band_count.  To generate the bands, the mass list (which is assumed
    already sorted by mass) is evenly partitioned into bands with masses in
    the range [mass_lb, mass_ub).

    >>> list(generate_mass_bands(1, [ float(x) for x in range(100) ]))
    [(1, 0.0, 100.0)]
    >>> list(generate_mass_bands(4, [ float(x) for x in range(100) ]))
    [(1, 0.0, 25.0), (2, 25.0, 50.0), (3, 50.0, 75.0), (4, 75.0, 100.0)]
    >>> list(generate_mass_bands(3, [ float(x) for x in range(100) ]))
    [(1, 0.0, 34.0), (2, 34.0, 68.0), (3, 68.0, 100.0)]

    """
    assert band_count > 0 and mass_list
    assert sorted(mass_list) == list(mass_list)
    band_size = int(math.ceil(float(len(mass_list)) / band_count))
    lb = mass_list[0]
    for bn in range(1, band_count):
        i = min(bn*band_size, len(mass_list)-1)
        ub = mass_list[i]
        yield bn, lb, ub
        lb = ub
    # Since the upper bound is exclusive, we add a little slop to make sure
    # the last spectrum is included.
    yield band_count, lb, round(mass_list[-1]+1)


def pythonize_swig_object(o, only_fields=None, skip_fields=[]):
    """Creates a pure Python copy of a SWIG object, so that it can be easily
    pickled, or printed (for debugging purposes).  Each SWIG object is
    pythonized as a dictionary, for convenience.  If provided, 'only_fields',
    limits the copy to the list of fields specified.  Otherwise,
    'skip_fields' if given is a list of methods not to include (this helps
    avoid infinite recursion).

    >>> pprint(pythonize_swig_object(cgreylag.score_stats(1, 1)))
    {'best_matches': [[{'mass_trace': [],
                        'missed_cleavage_count': -1,
                        'peptide_begin': [],
                        'peptide_sequence': '',
                        'predicted_parent_mass': 0.0,
                        'score': 0.0,
                        'sequence_name': [],
                        'spectrum_index': -1}]],
     'candidate_spectrum_count': 0,
     'combinations_searched': 0}

    """

    if isinstance(o, str):
        return o
    try:
        len(o)
    except TypeError:
        pass
    else:
        return list(pythonize_swig_object(x, only_fields, skip_fields)
                    for x in o)
    if hasattr(o, '__swig_getmethods__'):
        s = {}
        for a in o.__swig_getmethods__:
            if (only_fields != None and a in only_fields
                or only_fields == None and a not in skip_fields):
                s[a] = pythonize_swig_object(getattr(o, a), only_fields,
                                             skip_fields)
        return s
    return o

# A possible PCA (pyrrolidone carboxyl acid) modification accounts for
# circularization of the peptide N-terminal.  PCA mods are excluded if a
# static N-terminal mod has been specified.  Likewise, choosing a PCA mod will
# exclude choosing a potential N-terminal mod.  (The PCA mod for 'C' is
# dependent on a static mod of C+57 being in effect.)

def get_pca_table(mass_regimes):
    """Return a list of tuples (residues, parent delta, fragment delta)
    describing the PCA possibilities for each mass regime. Residues is a
    string specifying the residues this PCA delta applies to.
    """
    # FIX: According to Xtandem, C is only a candidate for PCA if
    # carboxyamidomethylated (C+57).  Currently we always search it.
    return [ [('', 0, 0),
              ('E', -1 * CP.parent_mass_regime[r].water_mass,
               -1 * CP.fragment_mass_regime[r].water_mass),
              ('QC', -1 * CP.parent_mass_regime[r].ammonia_mass,
               -1 * CP.fragment_mass_regime[r].ammonia_mass)]
             for r in range(len(mass_regimes)) ]


def enumerate_conjunction(mod_tree, limit, conjuncts=[]):
    if not mod_tree:
        if 0 < len(conjuncts) <= limit:
            yield conjuncts
        return
    first, rest = mod_tree[0], mod_tree[1:]
    if isinstance(first, list):
        for x in enumerate_disjunction(first, limit):
            for y in enumerate_conjunction(rest, limit, conjuncts + x):
                yield y
    else:
        for y in enumerate_conjunction(rest, limit, conjuncts):
            yield y
        for y in enumerate_conjunction(rest, limit, conjuncts + [first]):
            yield y

def enumerate_disjunction(mod_tree, limit=sys.maxint):
    """Generates the conjuncts for mod_tree that are no longer than limit.

    >>> list(enumerate_disjunction([['a'],['b'],['c']]))
    [[], ['a'], ['b'], ['c']]
    >>> list(enumerate_disjunction([[1,2,3]]))
    [[], [3], [2], [2, 3], [1], [1, 3], [1, 2], [1, 2, 3]]
    >>> list(enumerate_disjunction([[1,2,3],[4,5]]))
    [[], [3], [2], [2, 3], [1], [1, 3], [1, 2], [1, 2, 3], [5], [4], [4, 5]]
    >>> list(enumerate_disjunction([[1,2,3],[4,5]], limit=2))
    [[], [3], [2], [2, 3], [1], [1, 3], [1, 2], [5], [4], [4, 5]]
    >>> list(enumerate_disjunction([[1,2,3],[4,5]], limit=0))
    [[]]

    """
    assert isinstance(mod_tree, list)
    yield []
    for b in mod_tree:
        for s in enumerate_conjunction(b, limit):
            yield s

def get_mod_conjunct_triples(mod_tree, limit):
    """Return a triple (N, C, rest), where N (C) is a tuple of at most one N
    (C) conjunct, and rest is a tuple of conjuncts in a canonical order,
    ordered by increasing number of conjuncts, and with duplicates within and
    across removed.  A mass table is also appended to each conjunct.
    """
    # FIX: is there a more elegant way or place for all this?
    def enmass(t):
        def rmass(regime_index, par_frag, sign, delta, is_mono):
            global MASS_REGIME_ATOMIC_MASSES
            if isinstance(delta, str):
                if is_mono:
                    regime_index = 0
                return (sign * formula_mass(delta,
                           MASS_REGIME_ATOMIC_MASSES[regime_index][par_frag]))
            else:
                return sign * delta

        sign, delta, is_mono = t[:3]
        return t + (tuple((rmass(r, 0, sign, delta, is_mono),
                           rmass(r, 1, sign, delta, is_mono))
                          for r in range(len(XTP["mass_regimes"]))),)

    def triple(c):
        Ns = tuple(frozenset(enmass(x) for x in c if x[3] == '['))
        Cs = tuple(frozenset(enmass(x) for x in c if x[3] == ']'))
        assert len(Ns) <= 1 and len(Cs) <= 1
        rest = [ enmass(x) for x in c if x[3] not in '[]' ]
        rest = tuple(sorted(list(frozenset(rest))))
        return (Ns, Cs, rest)

    return sorted(list(frozenset(triple(conjunct)
                                 for conjunct
                                 in enumerate_disjunction(mod_tree, limit))),
                  key=lambda x: (sum(len(y) for y in x), x))


def gen_delta_bag_counts(i, remainder, bag):
    if i < 1:
        assert i == 0
        bag[0] = remainder
        yield tuple(bag)
        return
    for delta in range(1, remainder-i+1):
        bag[i] = delta
        for x in gen_delta_bag_counts(i-1, remainder-delta, bag):
            yield x

def generate_delta_bag_counts(mod_count, conjunct_length):
    """Generate all tuples of positive integers having length conjunct_length
    and sum mod_count.  As a special case, () is such a tuple having length 0
    and sum 0.

    >>> for i in range(6): print i, list(generate_delta_bag_counts(4, i))
    ...
    0 []
    1 [(4,)]
    2 [(3, 1), (2, 2), (1, 3)]
    3 [(2, 1, 1), (1, 2, 1), (1, 1, 2)]
    4 [(1, 1, 1, 1)]
    5 []
    >>> list(generate_delta_bag_counts(0, 0))
    [()]
    >>> list(generate_delta_bag_counts(0, 1))
    []

    """
    if conjunct_length == 0:
        return [()] if mod_count == 0 else []
    if mod_count < conjunct_length:
        return []
    return gen_delta_bag_counts(conjunct_length - 1, mod_count,
                                [0] * conjunct_length)


def set_context_conjuncts(context, mass_regime_index, N_cj, C_cj, R_cj):
    assert len(N_cj) < 1 and len(C_cj) < 1
    context.N_delta = 0
    if N_cj:
        context.N_delta = N_cj[0][5][mass_regime_index][1]
    context.C_delta = 0
    if C_cj:
        context.C_delta = C_cj[0][5][mass_regime_index][1]
    context.delta_bag_lookup.clear()
    context.delta_bag_lookup.resize(ord('Z')+1)
    context.delta_bag_delta.clear()
    for n, cj in enumerate(R_cj):
        context.delta_bag_delta.append(cj[5][mass_regime_index][1])
        for r in cj[3]:
            context.delta_bag_lookup[ord(r)] \
                = context.delta_bag_lookup[ord(r)] + (n,)


def search_all(options, context, score_statistics):
    """Search sequence database against searchable spectra."""

    mod_limit = XTP["potential_mod_limit"]

    mod_conjunct_triples = get_mod_conjunct_triples(XTP["potential_mods"],
                                                    mod_limit)
    info("%s unique potential modification conjuncts",
         len(mod_conjunct_triples))
    debug("mod_conjunct_triples (unique):\n%s", pformat(mod_conjunct_triples))

    mass_regimes = XTP["mass_regimes"]
    pca_table = get_pca_table(mass_regimes)
    debug("pca_table: %s", pca_table)

    total_combinations_searched = 0

    for mod_count in range(mod_limit + 1):
        context.mod_count = mod_count
        for mr_index, (mass_regime, pca_entry) in enumerate(zip(mass_regimes,
                                                                pca_table)):
            context.mass_regime_index = mr_index
            for pca_res, pca_parent_delta, pca_frag_delta in pca_entry:
                context.pca_residues = pca_res
                context.pca_delta = pca_frag_delta
                for cji, (N_cj, C_cj, R_cj) in enumerate(mod_conjunct_triples):
                    if pca_res and N_cj:
                        continue    # mutually exclusive, for now
                    set_context_conjuncts(context, mr_index, N_cj, C_cj, R_cj)
                    debug("mod_count: %s", mod_count)
                    debug("cj_triple: N=%s C=%s R=%s", N_cj, C_cj, R_cj)
                    for delta_bag in generate_delta_bag_counts(mod_count,
                                                               len(R_cj)):
                        debug("delta_bag: %s", delta_bag)

                        # this clear() avoids an SWIG/STL bug!?
                        context.delta_bag_count.clear()
                        context.delta_bag_count[:] = delta_bag

                        pmrf = CP.parent_mass_regime[mr_index].fixed_residue_mass
                        fmrf = CP.fragment_mass_regime[mr_index].fixed_residue_mass
                        p_fx = (pmrf[ord('[')]
                                + (N_cj and N_cj[0][5][mr_index][0] or 0)
                                + pmrf[ord(']')]
                                + (C_cj and C_cj[0][5][mr_index][0] or 0)
                                + pca_parent_delta + PROTON_MASS)
                        context.parent_fixed_mass = p_fx
                        f_N_fx = (fmrf[ord('[')]
                                  + (N_cj and N_cj[0][5][mr_index][1] or 0)
                                  + pca_frag_delta)
                        context.fragment_N_fixed_mass = f_N_fx
                        f_C_fx = (fmrf[ord(']')]
                                  + (C_cj and C_cj[0][5][mr_index][1] or 0)
                                  + CP.fragment_mass_regime[mr_index].water_mass)
                        context.fragment_C_fixed_mass = f_C_fx

                        info("MC=%s MR=%s PCA=%s CJ=%s DB=%s"
                             % (mod_count, mr_index, pca_res, cji, delta_bag))
                        debug("p_fx %s f_N_fx %s f_C_fx %s"
                              % (p_fx, f_N_fx, f_C_fx))

                        score_statistics.combinations_searched = 0
                        cgreylag.spectrum.search_runs(context,
                                                      score_statistics)
                        total_combinations_searched \
                            += score_statistics.combinations_searched
                        info("  %s candidate spectra examined, this bag",
                             score_statistics.combinations_searched)

    info('%s candidate spectra examined',
         score_statistics.candidate_spectrum_count)


def get_prefix_sequence(begin_pos, run_offset, sequence):
    """
    >>> get_prefix_sequence(0, 0, 'abcdefghi')
    '['
    >>> get_prefix_sequence(1, 0, 'abcdefghi')
    '[a'
    >>> get_prefix_sequence(3, 0, 'abcdefghi')
    '[abc'
    >>> get_prefix_sequence(4, 0, 'abcdefghi')
    'abcd'
    >>> get_prefix_sequence(4, 2, 'abcdefghi')
    'cdef'

    """
    prefix_start = run_offset + begin_pos - 4
    s = sequence[max(0, prefix_start):prefix_start+4]
    if len(s) < 4:
        s = '[' + s
    return s

def get_suffix_sequence(end_pos, run_offset, sequence):
    """
    >>> get_suffix_sequence(0, 0, 'abcdefghi')
    'abcd'
    >>> get_suffix_sequence(4, 0, 'abcdefghi')
    'efgh'
    >>> get_suffix_sequence(6, 0, 'abcdefghi')
    'ghi]'
    >>> get_suffix_sequence(8, 0, 'abcdefghi')
    'i]'
    >>> get_suffix_sequence(4, 2, 'abcdefghi')
    'ghi]'

    """
    suffix_start = run_offset + end_pos
    s = sequence[suffix_start:suffix_start+4]
    if len(s) < 4:
        s = s + ']'
    return s


def clean_string(v):
    """Strip and collapse internal whitespace.

    >>> clean_string(' one   two ')
    'one two'

    """
    if not v:
        return v
    return re.sub(r'[\s]+', ' ', v.strip())


def clean_defline(s):
    """Return the given string with tabs replaced by spaces and control
    and non-ASCII characters removed, then stripped.

    >>> tab=chr(9); clean_defline(' one' + tab + ' two three\001four ')
    'one  two threefour'

    """
    return re.sub(r'[^ -~]', '', s.replace('\t', ' ')).strip()


def abbrev_defline(s):
    """Return an abbreviated version of the defline--about 80 characters.

    >>> abbrev_defline('words ' * 10)
    'words words words words words words words words words words '
    >>> abbrev_defline('words ' * 20)
    'words words words words words words words words words words words words words words...'

    """
    ab = re.match(r'.{,80}\S{,170}', s).group(0)
    if len(ab) < len(s):
        ab += '...'
    return ab


def zopen(filename, mode='r', compresslevel=9):
    """Open a filename as with 'open', but using compression if indicated by
    the filename suffix."""
    if filename.endswith('.gz'):
        return gzip.GzipFile(filename, mode, compresslevel)
    else:
        return open(filename, mode)


def results_dump(score_statistics, searchable_spectra):
    """Return a result dict mapping spectrum names to (spectrum_metadata,
    best_matches) pairs.  (Unneeded fields are stripped.)
    """

    r = {}
    spectrum_metadata_fs = set(['name', 'file_id', 'mass', 'charge',
                                'total_ion_current', 'comparisons'])
    py_s_spectra = pythonize_swig_object(searchable_spectra,
                                         only_fields=spectrum_metadata_fs)
    py_matches = pythonize_swig_object(score_statistics.best_matches,
                                       skip_fields=['spectrum_index'])
    assert len(py_s_spectra) == len(py_matches)

    for sp_metadata, sp_matches in zip(py_s_spectra, py_matches):
        assert sp_metadata['name'] not in r, "duplicate spectrum name"
        r[sp_metadata['name']] = (sp_metadata, sp_matches)

    return r


def main(args=sys.argv[1:]):
    parser = optparse.OptionParser(usage=
                                   "usage: %prog [options] <job-id>"
                                   " <configuration-file> <ms2-file>...",
                                   description=__doc__, version=__version__)
    pa = parser.add_option
    pa("-P", "--parameter", nargs=2, dest="parameters", action="append",
       default=[],
       help="override a parameter in <parameter-file>, may be used multiple"
       " times", metavar="NAME VALUE")
    pa("-w", "--work-slice", nargs=2, type="float", dest="work_slice",
       help="search a subinterval [L:U) of the work space"
       " (where 0 <= L <= U <= 1) in standalone mode", metavar="L U")
    pa("-q", "--quiet", action="store_true", dest="quiet", help="no warnings")
    pa("-p", "--show-progress", action="store_true", dest="show_progress",
       help="show running progress")
    pa("--estimate", action="store_true", dest="estimate_only",
       help="just estimate the time required for the search")
    pa("-v", "--verbose", action="store_true", dest="verbose",
       help="be verbose")
    pa("--copyright", action="store_true", dest="copyright",
       help="print copyright and exit")
    pa("--debug", action="store_true", dest="debug",
       help="output debugging info")
    pa("--profile", action="store_true", dest="profile",
       help="dump Python profiling output to './greylag.prof.<pid>'")
    (options, args) = parser.parse_args(args=args)

    if options.copyright:
        print __copyright__
        sys.exit(0)

    if len(args) < 3:
        parser.print_help()
        sys.exit(1)

    job_id = args[0]
    configuration_fn = args[1]
    spectrum_fns = args[2:]

    if (any(True for f in spectrum_fns if not f.endswith('.ms2'))
        or (options.work_slice
            and not (0 <= options.work_slice[0]
                     <= options.work_slice[1] <= 1))):
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
    info("starting on %s", gethostname())

    # prevent format char problems
    if '%' in job_id:
        error("<job-id> may not contain '%'")

    # check spectrum basename uniqueness, as corresponding sqt files will be
    # in a single directory
    base_spectrum_fns = [ os.path.basename(fn) for fn in spectrum_fns ]
    if len(base_spectrum_fns) != len(set(base_spectrum_fns)):
        error("base spectrum filenames must be unique")

    # check -P names for validity
    bad_names = (set(n for n,v in options.parameters)
                 - set(PARAMETER_INFO.keys()))
    if bad_names:
        error("bad -P parameter names %s" % list(bad_names))

    # read params
    cp = ConfigParser.RawConfigParser()
    cp.optionxform = str                # be case-sensitive
    with open(configuration_fn) as configuration_file:
        cp.readfp(configuration_file)
    if not cp.has_section('greylag'):
        error("%s has no [greylag] section" % configuration_fn)
    parameters = dict(cp.items('greylag'))
    parameters.update(dict(options.parameters)) # command-line override
    global XTP
    XTP = validate_parameters(parameters)

    fixed_mod_map = dict((r[3], r) for r in XTP["pervasive_mods"])
    initialize_spectrum_parameters(options, XTP["mass_regimes"], fixed_mod_map)

    # read sequence dbs
    databases = XTP["databases"].split()
    # [(locusname, defline, seq, filename), ...]
    fasta_db = list(read_fasta_files(databases))
    # [(idno, offset, locusname, defline, seq, seq_filename), ...]
    db = []
    for idno, (locusname, defline, sequence, filename) in enumerate(fasta_db):
        db.extend(split_sequence_into_aa_runs(idno, locusname, defline,
                                              sequence, filename))
    db_residue_count = sum(len(dbi[3]) for dbi in db)

    info("read %s sequences (%s runs, %s residues)", len(fasta_db), len(db),
         db_residue_count)
    max_run_length = max(len(r[3]) for r in db)
    info("max run length is %s residues", max_run_length)
    if max_run_length > 2**31 - 1:
        error("runs longer than %s not yet supported", max_run_length)
    if not db:
        error("no database sequences")

    # read spectrum offset indices
    spectrum_offset_indices = []
    for spfn in spectrum_fns:
        idxfn = spfn + '.idx'
        with contextlib.closing(gzip.open(idxfn)) as idxf:
            idx = cPickle.load(idxf)
            # try to verify that index matches spectrum file
            if idx['file size'] != os.path.getsize(spfn):
                error("index '%s' needs rebuilding" % idxfn)
            spectrum_offset_indices.append(idx['offsets'])

    # FIX: assume standalone mode (for now)
    assert options.work_slice

    # read spectra per work slice
    spectra = read_spectra_slice(spectrum_fns, spectrum_offset_indices,
                                 options.work_slice)
    spectra.sort(key=lambda x: x.mass)

    if not spectra:
        warning("no input spectra")
    else:
        info("read %s spectra (mass range %s - %s)", len(spectra),
             spectra[0].mass, spectra[-1].mass)

    def peak_statistics(spectra):
        counts = [ len(sp.peaks) for sp in spectra ]
        counts.sort()
        n = len(counts)
        return (counts[0], counts[int(n*0.25)], counts[int(n*0.5)],
                counts[int(n*0.75)], counts[-1], sum(counts) / float(n))

    info("  peak stats: %s..%s..%s..%s..%s (mean=%.2f)"
         % peak_statistics(spectra))

    # filter and normalize spectra
    for sp in spectra:
        sp.filter_peaks(XTP["TIC_cutoff_proportion"],
                        CP.parent_mass_tolerance_max)
        sp.classify(XTP["intensity_class_count"], XTP["intensity_class_ratio"],
                    XTP["fragment_mass_tolerance"])

    min_psm = XTP["min_parent_spectrum_mass"]
    max_psm = XTP["max_parent_spectrum_mass"]
    # FIX: also filter by 1 + 2 + 4 rule?
    spectra = [ sp for sp in spectra
                if len(sp.peaks) >= 10 and min_psm <= sp.mass <= max_psm ]

    info("after filtering:")
    info("     %s spectra (mass range %s - %s)", len(spectra),
         spectra[0].mass, spectra[-1].mass)
    info("  peak stats: %s..%s..%s..%s..%s (mean=%.2f)"
         % peak_statistics(spectra))

    cgreylag.spectrum.set_searchable_spectra(spectra)
    score_statistics = cgreylag.score_stats(len(spectra),
                                            XTP["best_result_count"])

    if spectra:
        del spectra                     # release memory

        # FIX!!!
        # (cleavage_re, position of cleavage in cleavage_re)
        cleavage_motif = "[X]|[X]"
        cleavage_pattern, cleavage_pos = cleavage_motif_re(cleavage_motif)
        if cleavage_pos == None:
            error("cleavage site '%s' is missing '|'",
                  XTP["protein, cleavage site"])
        cleavage_pattern = re.compile(cleavage_pattern)

        context = cgreylag.search_context()
        for idno, offset, locusname, defline, seq, seq_filename in db:
            cp = []
            if cleavage_motif != "[X]|[X]":
                cp = list(generate_cleavage_points(cleavage_pattern,
                                                   cleavage_pos, seq))
            sr = cgreylag.sequence_run(idno, offset, seq, cp, locusname)
            context.sequence_runs.append(sr)
        context.maximum_missed_cleavage_sites = 1000000000 # FIX

        info("searching")
        search_all(options, context, score_statistics)
    else:
        warning("no spectra after filtering--search skipped")

    if options.estimate_only:
        print ("%.2f generic CPU hours"
               % (score_statistics.candidate_spectrum_count / 300.0e6))
        return

    # FIX!
    info("writing result file")
    result_fn = 'test.result.gz'
    with contextlib.closing(zopen(result_fn, 'w')) as result_file:
        d = { 'version' : __version__,
              'matches' : results_dump(score_statistics,
                                       cgreylag.cvar.spectrum_searchable_spectra),
              'total comparisons' : score_statistics.candidate_spectrum_count,
              'spectrum files' : base_spectrum_fns,
              'databases' : databases,
              'parameters' : XTP,
              'mass regime atomic masses' : MASS_REGIME_ATOMIC_MASSES,
              'proton mass' : PROTON_MASS,
              'argv' : sys.argv }
        cPickle.dump(d, result_file, cPickle.HIGHEST_PROTOCOL)
    info("finished, result file written to '%s'", result_fn)


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
            import cProfile
            import pstats
            report_fn = "greylag.prof.%s" % os.getpid()
            data_fn = report_fn + ".tmp"
            prof = cProfile.run('main()', data_fn)
            with open(report_fn, 'w') as report_f:
                try:
                    stats = pstats.Stats(data_fn, stream=report_f)
                    stats.strip_dirs()
                    stats.sort_stats('cumulative')
                    stats.print_stats(50)
                    stats.sort_stats('time')
                    stats.print_stats(50)
                    print "# profile report written to '%s'" % report_fn
                finally:
                    try:
                        os.remove(data_fn)
                    except:
                        pass
        else:
            main()
    except SystemExit:
        raise
    except:
        logging.exception("unhandled exception")
        sys.exit(1)
    finally:
        logging.shutdown()


# FIXES:
# - need to rigorously check for bad input in any file (fasta, spectrum, xml)
