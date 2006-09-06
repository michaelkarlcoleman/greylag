#!/usr/bin/env python2.4

'''
Partial re-implementation of X!Tandem, with extra features.

'''

__version__ = "$Id$"


import bisect
import fileinput
import math
import operator
import optparse
import os
import os.path
from pprint import *
import re
import sys

import elementtree.ElementTree

import cxtpy


#import gc
#gc.disable()
#gc.set_debug(gc.DEBUG_STATS|gc.DEBUG_LEAK)


def warn(s):
    print >> sys.stderr, 'warning:', s
def error(s):
    sys.exit('error: ' + s)
def fileerror(s):
    error(s + (", at line %s of file '%s'"
               % (fileinput.filelineno(), fileinput.filename())))


# try:
#     import psyco
#     psyco.full()
#     warn('psyco enabled')
# except:
#     pass


# name -> value map of processed XML parameters
XTP = {}

# handle to the singleton parameter object shared with the C++ module
CP = cxtpy.cvar.parameters_the

PROTON_MASS = 1.007276

# atom -> (monoisotopic mass, average mass)
ATOMIC_MASS = {
    'H' : (  1.007825035, 1.00794  ),
    'O' : ( 15.99491463, 15.9994   ),
    'N' : ( 14.003074,   14.0067   ),
    'C' : ( 12.0,        12.0107   ),
    'S' : ( 31.9720707,  32.065    ),
    'P' : ( 30.973762,   30.973761 ),
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
    

def initialize_spectrum_parameters():
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
        # NB: swig currently only exposes inner vector as tuple
        # CP.potential_modification_mass[ord(residue)].append(modvalue)
        CP.potential_modification_mass[ord(residue)] \
            = CP.potential_modification_mass[ord(residue)] + (modvalue,)
    CP.potential_modification_mass_refine.resize(128)
    for residue, modvalue in XTP["refine, potential modification mass"]:
        # CP.potential_modification_mass_refine[ord(residue)].append(modvalue)
        CP.potential_modification_mass_refine[ord(residue)] \
            = CP.potential_modification_mass_refine[ord(residue)] + (modvalue,)

    CP.cleave_N_terminal_mass_change = XTP["protein, cleavage N-terminal mass change"]
    CP.cleave_C_terminal_mass_change = XTP["protein, cleavage C-terminal mass change"]

    CP.fragment_mass_error = XTP["spectrum, fragment mass error"]

    CP.proton_mass = PROTON_MASS
    CP.water_mass = formula_mass("H2O")

    # CP.factorial[n] == (double) n!
    CP.factorial.resize(100, 1.0)
    for n in range(2, len(CP.factorial)):
        CP.factorial[n] = CP.factorial[n-1] * n


def cleavage_motif_re(motif):
    """Return (regexp, pos), where regexp is a regular expression that will
    match a cleavage motif, and pos is the position of the cleavage with
    respect to the match (e.g., 1 for '[KR]|{P}', or None if absent).
    """
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


def read_spectra_from_ms2_file(fn):
    """Return a list of spectrum objects read from an ms2 file."""
    f = file(fn)
    spectra = []
    while True:
        s = cxtpy.spectrum()
        if not s.read(f):
            break
        spectra.append(s)
    return spectra


def generate_peptides(seq, cleavage_points, maximum_missed_cleavage_sites):
    """Yield (begin, end) for each apt peptide in sequence."""
    len_cp = len(cleavage_points)
    for begin_i in range(len_cp-1):
        for end_i in range(begin_i+1,
                           min(len_cp,
                               begin_i+2+maximum_missed_cleavage_sites)):
            begin, end = cleavage_points[begin_i], cleavage_points[end_i]
            if end - begin >= 4:        # make 4 a param!
                yield begin, end


def generate_cleavage_points(sequence):
    """Yields the offsets of the cleavages in sequence."""
    yield 0
    # if pattern length > 1, will miss overlapping matches!
    for m in re.finditer(r'[KR](?=[^P])', sequence):
        yield m.start() + 1
    yield len(sequence)

aa_sequence = re.compile(r'[ARNDCQEGHILKMFPSTWYV]+')

def split_sequence_into_aa_runs(idno, defline, sequence):
    """Returns a tuple (idno, start, defline, seq, is_N, is_C) for each
    contiguous run of residues in sequence, where 'start' is the position or
    'seq' in 'sequence', and 'is_N' and 'is_C' are flags indicating whether
    the the beginning/end of 'seq' is at the beginning/end of 'sequence'."""
    matches = list(aa_sequence.finditer(sequence))
    return [ (idno, m.start(), defline, m.group(),
              n==0 and m.start()==0,
              n==len(matches)-1 and m.end()==len(sequence))
             for n, m in enumerate(matches) ]

def read_fasta_files(filenames):
    """Yield (defline, sequence) pairs as read from FASTA files (uppercasing
    sequence)."""
    defline = None
    seqs = []
    for line in fileinput.input(filenames):
        line = line.strip()
        if line[:1] == '>':
            if defline:
                yield (defline, ''.join(seqs))
            elif seqs:
                fileerror("bad format: line precedes initial defline")
            defline = line[1:]
            seqs = []
        else:
            seqs.append(line.upper())
    if defline:
        yield (defline, ''.join(seqs))


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

# name -> (type, default, check_fn)
# type may be bool, int, float, str, modlist, or a tuple of values
# default, as a string value, or None if value must be explicitly specified
# check_fn is an optional function that returns True iff the value is valid
# *1: param required, as has differing defaults based on other params, in
# xtandem

def mod_list(modification_list_specification, unique_mods=True):
    """Check and return a list of (residue, modification) tuples.  If
    unique_mods is True, also check that each residue is specified only once."""
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
    return speclist

def potential_mod_list(modification_list_specification):
    """Check and return a list of (residue, modification) tuples."""
    return mod_list(modification_list_specification, unique_mods=False)
        

# "ni" means check verifies that "not implemented" functionality is not
# specified
def p_ni_true(x): return x
def p_ni_false(x): return not x
def p_ni_empty(x): return len(x) == 0
def p_ni_daltons(x): return x == "Daltons"
def p_ni_zero(x): return x == 0
def p_positive(x): return x > 0
def p_negative(x): return x < 0
def p_nonnegative(x): return x >= 0

XML_PARAMETER_INFO = {
    "list path, default parameters" : (str, ""),
    "list path, taxonomy information" : (str, None),
    "output, histogram column width" : (int, 30, p_positive),
    "output, histograms" : (bool, "no"),
    "output, log path" : (str, ""),     # ignored
    "output, maximum valid expectation value" : (float, None),
    "output, message" : (str, "."),     # ignored
    "output, one sequence copy" : (bool, "no"),
    "output, parameters" : (bool, "no"),
    "output, path hashing" : (bool, "no"),
    "output, path" : (str, ""),
    "output, performance" : (bool, "no"),
    "output, proteins" : (bool, "no"),
    "output, results" : (('all', 'valid', 'stochastic'), "all"),
    "output, sequence path" : (str, "", p_ni_empty),
    "output, sequences" : (bool, "no"),
    "output, sort results by" : (('protein', 'spectrum'), "spectrum"),
    "output, spectra" : (bool, "no"),
    "output, xsl path" : (str, ""),
    "protein, C-terminal residue modification mass" : (float, "0.0", p_ni_zero), # eliminate, obsolete
    "protein, N-terminal residue modification mass" : (float, "0.0", p_ni_zero), # eliminate, obsolete
    "protein, cleavage C-terminal mass change" : (float, formula_mass("OH")),
    "protein, cleavage N-terminal limit" : (int, "100000000", p_positive),
    "protein, cleavage N-terminal mass change" : (float, formula_mass("H")),
    "protein, cleavage semi" : (bool, "no", p_ni_false),
    "protein, cleavage site" : (("[RK]|{P}",), None),
    "protein, homolog management" : (bool, "no", p_ni_false),
    "protein, modified residue mass file" : (str, "", p_ni_empty),
    "protein, taxon" : (str, None),
    "refine" : (bool, "no", p_ni_false),
    "refine, cleavage semi" : (bool, "no", p_ni_false),
    "refine, maximum valid expectation value" : (float, None),
    "refine, modification mass" : (mod_list, "", p_ni_empty),
    "refine, point mutations" : (bool, "no", p_ni_false),
    "refine, potential C-terminus modifications": (potential_mod_list, "", p_ni_empty),
    "refine, potential N-terminus modifications": (potential_mod_list, "", p_ni_empty),
    "refine, potential N-terminus modification position limit" : (int, 50),
    "refine, potential modification mass" : (potential_mod_list, "", p_ni_empty),
    "refine, potential modification motif" : (str, "", p_ni_empty),
    "refine, sequence path" : (str, "", p_ni_empty),
    "refine, spectrum synthesis" : (bool, "no", p_ni_false),
    "refine, tic percent" : (float, 20.0, p_nonnegative),
    "refine, unanticipated cleavage" : (bool, p_ni_false),
    "refine, use potential modifications for full refinement" : (bool, "no", p_ni_false),
    "residue, modification mass" : (mod_list, ""),
    "residue, potential modification mass" : (potential_mod_list, ""),
    "residue, potential modification motif" : (str, "", p_ni_empty),
    "scoring, a ions" : (bool, "no", p_ni_false),
    "scoring, b ions" : (bool, "yes", p_ni_true),
    "scoring, c ions" : (bool, "no", p_ni_false),
    "scoring, cyclic permutation" : (bool, "no", p_ni_false),
    "scoring, include reverse" : (bool, "no", p_ni_false),
    "scoring, maximum missed cleavage sites" : (int, None, p_nonnegative),
    "scoring, minimum ion count" : (int, None, p_positive), # v-1 -> m_lIonCount!
    "scoring, pluggable scoring" : (bool, "no"), # ignored
    "scoring, x ions" : (bool, "no", p_ni_false),
    "scoring, y ions" : (bool, "yes", p_ni_true),
    "scoring, z ions" : (bool, "no", p_ni_false),
    "spectrum, check all charges" : (bool, "no", p_ni_false),
    "spectrum, dynamic range" : (float, "100", p_positive),
    "spectrum, fragment mass error units" : (("Daltons", "ppm"), "Daltons", p_ni_daltons),
    "spectrum, fragment mass error" : (float, "0.45", p_positive),
    "spectrum, fragment mass type" : (("average", "monoisotopic"), "monoisotopic"),
    "spectrum, homology error" : (float, "4.5", p_positive),
    "spectrum, maximum parent charge" : (int, "4", p_positive),
    "spectrum, minimum fragment mz" : (float, "200.0", p_positive),
    "spectrum, minimum parent m+h" : (float, "850.0", p_positive),
    "spectrum, minimum peaks" : (int, "5", p_positive),
    "spectrum, neutral loss mass" : (float, "0.0", p_nonnegative),
    "spectrum, neutral loss window" : (float, "0.0", p_nonnegative),
    "spectrum, parent monoisotopic mass error minus" : (float, None, p_negative), # *1
    "spectrum, parent monoisotopic mass error plus" : (float, None, p_positive), # *1
    "spectrum, parent monoisotopic mass error units" : (("Daltons", "ppm"), "Daltons", p_ni_daltons),
    "spectrum, parent monoisotopic mass isotope error" : (bool, "no", p_ni_false), # parent error should be <0.5Da
    "spectrum, path" : (str, None),
    "spectrum, sequence batch size" : (int, 1000), # ignored
    "spectrum, threads" : (int, 1),     # ignored
    "spectrum, total peaks" : (int, "50", p_positive),
    "spectrum, use conditioning" : (bool, "yes"),
    "spectrum, use contrast angle" : (bool, "no", p_ni_false),
    "spectrum, use neutral loss window" : (bool, "no", p_ni_false),
    "spectrum, use noise suppression" : (bool, "yes", p_ni_false),
}

def validate_parameters(parameters):
    """Verify that parameters are valid, have valid values, and correspond to
    currently implemented functionality.  Default values are filled in, and a
    name/value dict returned."""

    pmap = {}
    for name, info in sorted(XML_PARAMETER_INFO.items()):
        type_ = info[0]
        default = info[1]
        check_fn = len(info) > 2 and info[2] or None

        v = parameters.get(name)
        if isinstance(v, str):
            v = v.strip()
        if v == None:
            if default != None:
                warn("parameter '%s' defaulting to '%s'" % (name, default))
                v = default
            else:
                error("missing required parameter '%s'" % name)
        if isinstance(type_, tuple):
            if not v in type_:
                error("parameter '%s' value '%s' not in %s (feature not implemented?)"
                      % (name, v, type_))
        elif type_ == bool:
            v = { 'yes' : True, 'no' : False }.get(v)
            if v == None:
                error("parameter '%s' requires a value of 'yes' or 'no'")
        else:
            try:
                v = type_(v)
            except ValueError, e:
                error("parameter '%s' has value '%s' with invalid format [%s]"
                      % (name, v, e))
        if check_fn and not check_fn(v):
            error("parameter '%s' has invalid value '%s' (or feature not implemented)"
                  % (name, v))
        pmap[name] = v

    unknown_parameters = set(parameters) - set(XML_PARAMETER_INFO)
    if unknown_parameters:
        warn("%s unknown parameters:\n %s"
             % (len(unknown_parameters), pformat(sorted(list(unknown_parameters)))))

    # FIX: this would trip the "cyclic" param (assume Daltons)
    assert abs(pmap["spectrum, parent monoisotopic mass error plus"]) > 0.095
    assert abs(pmap["spectrum, parent monoisotopic mass error minus"]) > 0.095

    return pmap


# probably moves to C++
# just a stub, 0 == "no mod", first and last are "[" and "]"
def generate_mod_patterns(seq, begin, end):
    # just the no-mod pattern for now
    return ([0] * (end - begin + 2),)

# probably moves to C++, mods NYI
# parent mass assumed to be average (CHECK!)
def get_peptide_mod_mass(peptide_seq, potential_mod_pattern,
                     is_N, is_C):
    return 0


def get_mz(mass, charge):
    return mass/charge + CP.proton_mass


def synthetic_B_spectrum(peptide_mass, peptide_seq, peptide_mod_pattern,
                         peptide_is_N):
    m = 0 + (CP.cleave_N_terminal_mass_change - formula_mass("H"))
    if peptide_is_N:
        m += CP.modification_mass[ord('[')]
        if False:                       # [ is diff modded
            #m += CP.potential_modification_mass[ord('[')]...
            pass
    ladder = [None] * (len(peptide_seq)-1)
    for i in range(len(ladder)):
        print >> sys.stderr, "ladder step %s %s" % (peptide_seq[i], CP.monoisotopic_residue_mass[ord(peptide_seq[i])])
        m += CP.monoisotopic_residue_mass[ord(peptide_seq[i])]
        # m += the mod delta, too
        intensity = 1.0
        if i == 1:
            if peptide_seq[0] == 'P':
                intensity = 10.0
            else:
                intensity = 3.0
        ladder[i] = (m, intensity)
    return ladder


def synthetic_Y_spectrum(peptide_mass, peptide_seq, peptide_mod_pattern,
                         peptide_is_C):
    m = CP.water_mass + (CP.cleave_C_terminal_mass_change - formula_mass("OH"))
    if peptide_is_C:
        m += CP.modification_mass[ord(']')]
        if False:                       # ] is diff modded
            #m += CP.potential_modification_mass[ord(']')]...
            pass
    ladder = [None] * (len(peptide_seq)-1)
    for i in range(len(ladder)-1, -1, -1):
        print >> sys.stderr, "ladder step %s %s" % (peptide_seq[i+1], CP.monoisotopic_residue_mass[ord(peptide_seq[i+1])])
        m += CP.monoisotopic_residue_mass[ord(peptide_seq[i+1])]
        # m += the mod delta, too
        intensity = 1.0
        if i == 1:
            if peptide_seq[0] == 'P':
                intensity = 10.0
            else:
                intensity = 3.0
        ladder[len(ladder)-1-i] = (m, intensity)
    return ladder


# move intesity synthesis and get_mz down here?

def synthetic_spectrum(ion_type, peptide_mass, charge, peptide_seq,
                       peptide_mod_pattern, peptide_is_N, peptide_is_C):
    assert XTP["spectrum, fragment mass type"] == "monoisotopic", "frag avg NYI"
    if ion_type == 'B':
        ladder = synthetic_B_spectrum(peptide_mass, peptide_seq,
                                      peptide_mod_pattern, peptide_is_N)
    elif ion_type == 'Y':
        ladder = synthetic_Y_spectrum(peptide_mass, peptide_seq,
                                      peptide_mod_pattern, peptide_is_C)
    else:
        assert False, "unimplemented ion type"

    sp = cxtpy.spectrum(peptide_mass, charge)
    #sp.peaks[:] = [ cxtpy.peak(get_mz(mass, charge), intensity)
    #                for mass, intensity in ladder ]
    sp.set_peaks_from_matrix([ (get_mz(mass, charge), intensity)
                               for mass, intensity in ladder ])

    print >> sys.stderr, "# synth %s sp: mass = %s, charge = %s" % (ion_type, sp.mass, sp.charge)
    for p in sp.peaks:
        print >> sys.stderr, "#           %s %s" % (p.mz, p.intensity)
        
    
    return sp


# probably moves to C++, mods NYI
def synthetic_spectra(peptide_seq, peptide_mod_pattern, peptide_mass,
                      peptide_is_N, peptide_is_C, max_spectrum_charge):
    # This follows xtandem's assumption that precursors with charge z can
    # produce fragments with charge at most z-1, except that charge 1
    # precursors can produce charge 1 fragments.

    return [ (ion_type, [ (charge,
                           synthetic_spectrum(ion_type, peptide_mass, charge,
                                              peptide_seq, peptide_mod_pattern,
                                              peptide_is_N, peptide_is_C)) 
                          for charge in range(1, max(2, max_spectrum_charge)) ])
             for ion_type in ('B', 'Y') ]


# probably moves to C++, mods NYI
def score_spectrum(synth_spectra, spectrum, spectrum_charge):
    # IDEA: could we combine all these spectra and just do the correlation
    # once?

    hyper_score = 1.0
    convolution_score = 0
    ion_peaks = []
    ion_scores = []
    for ion_type, ion_spectra in synth_spectra:
        i_peaks = 0
        i_scores = 0
        for charge, ion_spectrum in ion_spectra:
            # FIX: move this test outward?
            if (spectrum_charge == 1 and charge > spectrum_charge
                or spectrum_charge > 1 and charge > spectrum_charge - 1):
                continue

            # convolution score is just sum over charges/ions
            # hyperscore is product of p! over charges/ions (where p is corr
            # peak count) times the convolution score (clipped to FLT_MAX)
            # > blurred!

            conv_score, common_peak_count \
                        = cxtpy.spectrum.score_similarity(ion_spectrum,
                                                             spectrum)
            i_peaks += common_peak_count
            i_scores += conv_score
            hyper_score *= CP.factorial[common_peak_count]
            convolution_score += conv_score
        ion_peaks.append((ion_type, i_peaks))
        ion_scores.append((ion_type, i_scores))

    # want to know:
    # - total peaks (over charges) for each ion_type
    # - total convolution score (over charges) for each ion_type
    # - maybe grand total convolution score
    # - hyper score: product of all factorial(peak count) times
    #                grand total convolution score
    #   [xtandem clips hyper score to FLTMAX]

    hyper_score *= convolution_score
    return hyper_score, convolution_score, ion_scores, ion_peaks


def scale_hyperscore(hyper_score):
    return int(round(4 * math.log10(hyper_score)))


def get_spectrum_expectation(hyper_score, hyperscore_histogram):
    """Return the expectation value for this hyperscore, based on the
    histogram."""
    scaled_hyper_score = scale_hyperscore(hyper_score)
    max_histogram_index = max(hyperscore_histogram)
    assert max_histogram_index < 1000
    # histogram as a list
    histogram = [0] * (max_histogram_index+1)
    for score, count in hyperscore_histogram.iteritems():
        histogram[score] = count

    counts = sum(histogram)
    if counts < 200:
        # use default line
        return 10.0 ** (3.5 + -0.18 * scaled_hyper_score)
    c = counts
    # survival is just a reversed cumulative distribution
    survival = [0] * len(histogram)
    for i in range(len(histogram)):
        survival[i] = c
        c -= histogram[i]
    # note: neither histogram nor survival have 0 as rightmost element

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
    return 10.0 ** (b + m * scaled_hyper_score)


def get_final_protein_expect(sp_count, valid_spectra, match_ratio, raw_expect,
                             db_length, candidate_spectra):
    # (From xtandem) Compensate for multiple peptides supporting a protein.
    # The expectation values for the peptides are combined with a simple
    # Bayesian model for the probability of having two peptides from the same
    # protein having the best score in different spectra.
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
    p = min(match_ratio / candidate_spectra, 0.9999999)
    r += (sp_count * math.log10(p)
          + (valid_spectra - sp_count) * math.log10(1 - p))
    return r
    

def main():
    parser = optparse.OptionParser(usage=
                                   "usage: %prog [options] <parameter-file>",
                                   description=__doc__)
    parser.add_option("-o", "--output", dest="output",
                      help="destination file [default as given in parameter"
                      " file, '-' for stdout]",
                      metavar="FILE")
    parser.add_option("-v", "--verbose", action="store_true",
                      dest="verbose", help="be verbose")
    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.print_help()
        sys.exit(1)


    # read params
    parameters = read_xml_parameters(args[0])
    default_parameter_fn = parameters.get("list path, default parameters")
    if default_parameter_fn:
        default_parameters = read_xml_parameters(default_parameter_fn)
        default_parameters.update(parameters)
        parameters = default_parameters
    global XTP
    XTP = validate_parameters(parameters)
    #pprint(XTP)

    if options.output:
        if options.output != '-':
            sys.stdout = open(options.output, 'w')
    else:
        output_path = XTP["output, path"]
        if output_path:
            sys.stdout = open(output_path, 'w')

    taxonomy = read_taxonomy(XTP["list path, taxonomy information"])

    initialize_spectrum_parameters()

    # read sequence dbs
    fasta_db = list(read_fasta_files(taxonomy[XTP["protein, taxon"]]))
    db = []
    for idno, (defline, sequence) in enumerate(fasta_db):
        db.extend(split_sequence_into_aa_runs(idno, defline, sequence))
    print "read %s sequences (%s runs)" % (len(fasta_db), len(db))
    if not db:
        error("no database sequences")

    db = [ (idno, offset, defline, seq, is_N, is_C,
            list(generate_cleavage_points(seq)))
           for idno, offset, defline, seq, is_N, is_C in db ]
    db_residue_count = sum(len(dbi[3]) for dbi in db)
    print "cleavage_points found"

    # read spectra
    spectrum_fn = XTP["spectrum, path"]
    spectra = read_spectra_from_ms2_file(spectrum_fn)
    if not spectra:
        error("no input spectra")
    print "read %s spectra" % len(spectra)

    spectra.sort(key=lambda x: x.mass)
    print "sorted spectra"

    # filter and normalize spectra
    # NI: optionally using "contrast angle" (mprocess::subtract)
    if XTP["spectrum, use conditioning"]:
        spectra = [ s for s in spectra
                    if s.filter_and_normalize(XTP["spectrum, minimum fragment mz"],
                                              XTP["spectrum, dynamic range"],
                                              XTP["spectrum, minimum peaks"],
                                              XTP["spectrum, total peaks"]) ]
    print "     %s spectra after filtering" % len(spectra)

    for s in spectra:
        print 'spectrum', s.mass, s.charge
        for p in s.peaks:
            print p.mz, p.intensity

    spectrum_mass_index = ([ (sp.mass, sp.charge, n) for n, sp in enumerate(spectra) ]
                           + [ (sp.secondary_mass, sp.secondary_charge, n)
                               for n, sp in enumerate(spectra)
                               if sp.secondary_charge ])
    spectrum_mass_index.sort()

    # spectrum index is wrt 'spectra' above
    # spectrum index -> [ <match info>, ... ]
    best_match = {}
    # spectrum index -> best_hyperscore
    best_score = {}
    # spectrum index -> (2nd-best hyperscore, 2nd-best convolution score)
    second_best_score = {}
    # (spectrum index, charge) -> (scaled, binned hyperscore -> count)
    hyperscore_histogram = dict(((n, charge), {})
                                for m, charge, n in spectrum_mass_index)

    peptide_count = 0
    candidate_spectrum_count = 0
    peptides_w_candidate_spectra = 0
    for idno, offset, defline, seq, is_N, is_C, cleavage_points in db:
        #print '#', idno, offset, defline, seq, is_N, is_C, cleavage_points
        sys.stderr.write('p')
        if idno > 1000:
            warn('stopping at 1000')
            return
        
        for begin, end in generate_peptides(seq, cleavage_points,
                                            XTP["scoring, maximum missed cleavage sites"]):
            peptide_count += 1
            peptide_seq = seq[begin:end]
            peptide_mass = cxtpy.get_peptide_mass(peptide_seq, is_N, is_C)
            # pyro?
            for potential_mod_pattern in generate_mod_patterns(peptide_seq,
                                                               begin, end):
                peptide_mod_mass = (peptide_mass +
                                    get_peptide_mod_mass(peptide_seq,
                                                         potential_mod_pattern,
                                                         is_N, is_C))
                sp_mass_lb = peptide_mod_mass - XTP["spectrum, parent monoisotopic mass error plus"]
                sp_mass_ub = peptide_mod_mass - XTP["spectrum, parent monoisotopic mass error minus"]
                candidate_spectra_info \
                    = spectrum_mass_index[bisect.bisect_left(spectrum_mass_index, (sp_mass_lb,))
                                          :bisect.bisect_right(spectrum_mass_index, (sp_mass_ub,))]
                if not candidate_spectra_info:
                    continue
                peptides_w_candidate_spectra += 1 # FIX: omit?
                candidate_spectrum_count += len(candidate_spectra_info)
                #print peptide_mod_mass, candidate_spectra_info, \
                #      [ spectra[c[2]] for c in candidate_spectra_info ]

                synth_sp = synthetic_spectra(peptide_seq,
                                             potential_mod_pattern, 
                                             peptide_mod_mass, is_N, is_C,
                                             max(charge for m, charge, n
                                                 in candidate_spectra_info))
                warn('peptide: %s' % peptide_seq)
                warn('ssp: %s' % synth_sp)

                for c_sp_mass, c_sp_charge, c_sp_n in candidate_spectra_info:
                    warn('candidate: %s' % spectra[c_sp_n])
                    hyper_score, convolution_score, ion_scores, ion_peaks \
                                 = score_spectrum(synth_sp, spectra[c_sp_n], c_sp_charge)
                    warn('score: %s %s %s %s' % (hyper_score, convolution_score, ion_scores, ion_peaks))
                    if convolution_score > 2:
                        sp_ion_count = sum(n for ion_type, n in ion_peaks)
                        # update spectrum histograms
                        hh = hyperscore_histogram[(c_sp_n, c_sp_charge)]
                        scaled_hyper_score = scale_hyperscore(hyper_score)
                        if scaled_hyper_score in hh:
                            hh[scaled_hyper_score] += 1
                        else:
                            hh[scaled_hyper_score] = 1
                        
                        # incr m_tPeptideScoredCount
                    # nyi: permute stuff
                    has_b_and_y = len([1 for ion_type, n in ion_peaks if n > 0]) >= 2
                    if not has_b_and_y:
                        continue
                    
                    # check that parent masses are within error range (isotope ni)
                    # already done above (why does xtandem do it here?)
                    # if check fails, only eligible for 2nd-best record

                    # Remember all of the highest-hyper-scoring matches
                    # against each spectrum.  These might be in multiple
                    # domains (pos, length, mods) in multiple proteins.
                    # (Note: this effectively chooses the best parent charge,
                    # too.)
                    if hyper_score >= best_score.get(c_sp_n, 0):
                        # something like this
                        # BEWARE: below code relies on this tuple's order
                        match_info = (hyper_score, convolution_score,
                                      ion_scores, ion_peaks, c_sp_mass, 
                                      c_sp_charge, c_sp_n, begin, peptide_seq,
                                      potential_mod_pattern, peptide_mod_mass,
                                      idno, offset, defline, seq)
                        if hyper_score > best_score.get(c_sp_n, 0):
                            best_score[c_sp_n] = hyper_score
                            best_match[c_sp_n] = [match_info]
                        else:
                            best_match[c_sp_n].append(match_info)
                    # Also remember just the scores of the 2nd-best
                    # hyper-scoring
                    elif hyper_score > second_best_score.get(c_sp_n, 0):
                        second_best_score[c_sp_n] = (hyper_score, convolution_score)
                    
    print 'generated %s peptides' % peptide_count
    print '          %s peptides have candidate spectra' % peptides_w_candidate_spectra
    sys.stdout.flush()

    # calculate exp for best match for each spectrum
    #  (CHECK: what about multiple matches?  charges?)
    # FIX: xtandem adds to a spectrum's histogram for each domain, but we keep
    # things separate

    # spectrum index -> expectation value
    expect = {}
    # spectrum index -> (scaled, binned hyperscore -> count)
    best_histogram = {}
    for sp_n, match_info_list in best_match.iteritems():
        sp_hyper_score, sp_charge = match_info_list[0][0], match_info_list[0][5]
        assert sp_n not in expect, "shouldn't see multiple charges here"
        hh = hyperscore_histogram[(sp_n, sp_charge)]
        best_histogram[sp_n] = hh
        expect[sp_n] = get_spectrum_expectation(sp_hyper_score, hh)

    warn("expect: %s" % expect)

    # NB: at this point we're using 'expect' to keep track of whether a
    # spectrum still passes or not
    
    # filter out spectra w/ low exp values??
    for sp_n in best_match:
        if expect[sp_n] > XTP["refine, maximum valid expectation value"]:
            del expect[sp_n]

    # protein id -> list of spectra
    best_protein_matches = {}
    for sp_n, match_info_list in best_match.iteritems():
        for match in match_info_list:
            best_protein_matches.setdefault(match[11], []).append(match)

    # calculate exp for each protein
    valid_spectra_count = sum(1 for sp_n, exp in expect.iteritems()
                              if exp <= XTP["output, maximum valid expectation value"])
    best_histogram_sum = sum(sum(h.itervalues())
                             for sp_n, h in best_histogram.iteritems()
                             if sp_n in expect)
    if not best_match:
        error("no matches found")
    match_ratio = float(best_histogram_sum) / len(best_match)

    if XTP["output, results"] == "valid":
        for sp_n in best_match:
            if sp_n in expect and expect[sp_n] > 0.95 * XTP["output, maximum valid expectation value"]:
                del expect[sp_n]

    # spectrum id's for spectra that are "repeats"; a repeat is a spectrum for
    # which a better corresponding spectrum exists having the same domain 0
    # match.
    repeats = set()
    for pid, matches in best_protein_matches.iteritems():
        domain_0_matches = {}           # spectrum id -> match
        for m in matches:
            if m[6] not in domain_0_matches:
                domain_0_matches[m[6]] = m
        for sid_x in domain_0_matches:
            if sid_x not in expect:
                continue
            for sid_y in domain_0_matches:
                if sid_y not in expect:
                    continue
                if sid_x < sid_y:
                    if (domain_0_matches[sid_x][7:9]
                        == domain_0_matches[sid_y][7:9]):
                        if expect[sid_x] > expect[sid_y]:
                            repeats.add(sid_y)
                        else:
                            repeats.add(sid_x)

    # protein id -> protein expectation value (log10 of expectation)
    raw_protein_expect = {}
    # protein id -> set of spectrum ids used in expectation value
    raw_protein_expect_spectra = {}
    for protein_id, matches in best_protein_matches.iteritems():
        for m in matches:
            spectrum_id = m[6]
            if spectrum_id in repeats or spectrum_id not in expect:
                continue
            log_expect = math.log10(expect[spectrum_id])
            if log_expect > -1.0:
                continue
            # this avoids adding matches for multiple domains (REWORK THIS?)
            if spectrum_id not in raw_protein_expect_spectra.get(protein_id, set()):
                raw_protein_expect[protein_id] = (raw_protein_expect.get(protein_id, 0)
                                                  + log_expect)
                raw_protein_expect_spectra.setdefault(protein_id,
                                                      set()).add(spectrum_id)
    
    bias = float(candidate_spectrum_count) / db_residue_count

    # protein id -> protein expectation value (log10 of expectation)
    protein_expect = {}
    for protein_id in raw_protein_expect:
        sp_count = len(raw_protein_expect_spectra[protein_id])
        protein_expect[protein_id] = get_final_protein_expect(sp_count,
                                                              len(expect),
                                                              match_ratio,
                                                              raw_protein_expect[protein_id],
                                                              len(fasta_db),
                                                              candidate_spectrum_count)
        if sp_count > 1:
            protein_length = len(fasta_db[protein_id][1])
            if protein_length * bias < 1.0:
                protein_expect[protein_id] += sp_count * math.log10(protein_length*bias)

    # protein id -> sum of peak intensities for supporting spectra
    intensity = {}
    # FIX: this is wrong because it counts all domains separately!
    for pid, matches in best_protein_matches.iteritems():
        for m in matches:
            sp_n = m[6]
            if sp_n in expect:
                intensity[pid] = (intensity.get(pid, 0.0)
                                  + sum(p.intensity for p in spectra[sp_n].peaks))

    # need ordering here?

    # temp output
    for sp_n, sp in enumerate(spectra):
        if sp_n in expect:
            print "%s %s" % (expect[sp_n], sp)
    


    # The "output," limit is the only one applied to proteins, but it is also
    # applied to peptides at one point.  The "refine," limit is only applied
    # to peptides.

    # limit to "output, maximum valid expectation value"


    # sort proteins by exp
    # within each protein, order spectra by start pos (then by what?)
    # handle domains!
    # output results



    # [strip control chars from deflines]
    


    # search
    #    [don't bother generating peptides that are two big (or too small?)
    #     for the given spectra, len must be >= 4 (make this a param)]

    # pyro_check: if no N-terminal mod otherwise specified, look for potential
    # mods for these N-terminal residues: Q or C+57 -> loss of ammonia,
    # E -> loss of water

    # don't score a spectrum unless it has a min # of peaks, and at least one
    # from each of ABC and XYZ





if __name__ == '__main__':
    main()
