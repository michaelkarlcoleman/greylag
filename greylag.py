# common greylag functions and constants

##  greylag, a collection of programs for MS/MS protein analysis
##  Copyright (C) 2006-2008  Stowers Institute for Medical Research
##
##  This program is free software: you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##  Contact: Mike Coleman
##           Stowers Institute for Medical Research
##           1000 East 50th Street
##           Kansas City, Missouri  64110
##           USA

from __future__ import with_statement


import logging; from logging import debug, info, warning
import math
import re
import sys

import cgreylag


VERSION = "0.1"


# handle to the singleton parameter object shared with the C++ module
CP = cgreylag.cvar.parameters_the


def set_logging(options):
    log_level = logging.WARNING
    if options.quiet:
        log_level = logging.ERROR
    if options.verbose:
        log_level = logging.INFO
    if options.debug:
        log_level = logging.DEBUG
    logfile = None
    if options.logfile:
        logfile = options.logfile
    logging.basicConfig(level=log_level, datefmt='%b %e %H:%M:%S',
                        format=('%(asctime)s [%(process)d]'
                                ' %(levelname)s: %(message)s'),
                        filename=logfile)


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


def initialize_spectrum_parameters(options, GLP, mass_regimes, fixed_mod_map):
    """Initialize parameters known to the spectrum module.
    fixed_mod_map maps, for example, 'M' to (1, 'O', False, 'M', 'oxidation').
    """

    # This function can be called multiple times to reinitialize for a new
    # search, though probably everything else (previous matches) ought to be
    # "forgotten" at that point, too.

    # This is the size of vectors that are indexed by residues (A-Z) or
    # special characters ('[]').
    RESIDUE_LIMIT = max(ord(c) for c in 'Z[]') + 1

    # These are currently monoisotopic.  (deuterium pointless?)
    CP.proton_mass = PROTON_MASS
    CP.hydrogen_mass = formula_mass("H")

    regime_manifest = []

    global MASS_REGIME_ATOMIC_MASSES
    # clear previous state
    MASS_REGIME_ATOMIC_MASSES = []
    CP.parent_mass_regime.clear()
    CP.fragment_mass_regime.clear()

    for rn, regime_pair in enumerate(mass_regimes):
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
                    m = (formula_mass(RESIDUE_FORMULA[r], atmass)
                         + GLP["mass_regime_debug_delta"])
                if n == 1:
                    regime_manifest.append((rn, r, m))
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
                CP.parent_mass_regime.append(creg)
            else:
                creg.fixed_residue_mass[ord('[')] -= CP.hydrogen_mass
                creg.fixed_residue_mass[ord(']')] -= creg.hydroxyl_mass
                CP.fragment_mass_regime.append(creg)
    for r in RESIDUES_W_BRACKETS:
        info('fixed mass %s: %s', r,
             [ ("%.6f" % CP.parent_mass_regime[rn].fixed_residue_mass[ord(r)],
                "%.6f" % CP.fragment_mass_regime[rn].fixed_residue_mass[ord(r)])
               for rn in range(len(mass_regimes)) ])
    for r in RESIDUES:
        for rn in range(len(mass_regimes)):
            # check for physically impossible/meaningless masses
            if CP.parent_mass_regime[rn].fixed_residue_mass[ord(r)] < 1.0:
                raise ValueError('bogus parent mass specification for %s' % r)
            if CP.fragment_mass_regime[rn].fixed_residue_mass[ord(r)] < 1.0:
                raise ValueError('bogus fragment mass specification for %s' % r)

    CP.parent_mass_tolerance_1 = GLP["parent_mz_tolerance"]
    CP.parent_mass_tolerance_max = (GLP["parent_mz_tolerance"]
                                    * GLP["charge_limit"])

    CP.fragment_mass_tolerance = GLP["fragment_mass_tolerance"]
    CP.intensity_class_count = GLP["intensity_class_count"]

    CP.minimum_peptide_length = GLP["min_peptide_length"]

    # CP.ln_factorial[n] == ln(n!)
    CP.ln_factorial.resize(int(GLP["max_parent_spectrum_mass"]
                               / GLP["fragment_mass_tolerance"] + 100), 0.0)
    for n in range(2, len(CP.ln_factorial)):
        CP.ln_factorial[n] = CP.ln_factorial[n-1] + math.log(n)

    # FIX or eliminate
    #CP.estimate_only = bool(options.estimate_only)
    #CP.show_progress = bool(options.show_progress)
    CP.estimate_only = False
    CP.show_progress = False

    return regime_manifest


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

def get_mod_conjunct_triples(mod_tree, limit, mass_regimes):
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
                          for r in range(len(mass_regimes))),)

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
        error('invalid cleavage motif pattern')
    i = 0
    for part in parts:
        if part == '|':
            if cleavage_pos != None:
                error ("invalid cleavage motif pattern" " (multiple '|'s)")
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


def file_sha1(filename):
    """Return the (binary) SHA1 digest of the given file."""
    try:
        import hashlib
    except:
        return "no checksum--libs missing"
    h = hashlib.sha1()
    h.update(open(filename).read())
    return h.digest()


def read_fasta_files(filenames):
    """Yield (locusname, defline, sequence, filename) tuples as read from
    FASTA files (uppercasing sequence)."""

    loci_seen = set()

    for filename in filenames:
        locusname, defline = None, None
        seqs = []
        with open(filename) as f:
            for line in f:
                line = line.strip()
                if line[:1] == '>':
                    if defline != None:
                        yield (locusname, defline, ''.join(seqs), filename)
                    elif seqs:
                        error("bad format: line precedes initial defline"
                              " in '%s'" % filename)
                    defline = line[1:]
                    locusname_rest = defline.split(None, 1)
                    if not locusname_rest:
                        error("empty locus name not allowed in '%s'" % filename)
                    locusname = locusname_rest[0]
                    if locusname in loci_seen:
                        error("locus name '%s' is not unique in the search"
                              " database(s) in '%s'" % (locusname, filename))
                    loci_seen.add(locusname)
                    seqs = []
                else:
                    seqs.append(line.upper())
            if defline:
                yield (locusname, defline, ''.join(seqs), filename)


