#!/usr/bin/env python2.4

'''
Partial re-implementation of X!Tandem, with extra features.

'''

__version__ = "$Id$"


import bisect
import fileinput
import operator
import optparse
import os
import os.path
from pprint import *
import re
import sys

import elementtree.ElementTree

import spectrum


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
CP = spectrum.cvar.parameters_the

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
        # FIX: figure out why this doesn't work!
        # CP.potential_modification_mass[ord(residue)].append(modvalue)
        CP.potential_modification_mass[ord(residue)] \
            = list(CP.potential_modification_mass[ord(residue)]) + [modvalue]
    CP.potential_modification_mass_refine.resize(128)
    for residue, modvalue in XTP["refine, potential modification mass"]:
        # FIX: figure out why this doesn't work!
        # CP.potential_modification_mass_refine[ord(residue)].append(modvalue)
        CP.potential_modification_mass_refine[ord(residue)] \
            = list(CP.potential_modification_mass_refine[ord(residue)]) + [modvalue]

    CP.cleave_N_terminal_mass_change = XTP["protein, cleavage N-terminal mass change"]
    CP.cleave_C_terminal_mass_change = XTP["protein, cleavage C-terminal mass change"]

    CP.proton_mass = PROTON_MASS
    CP.water_mass = formula_mass("H2O")


def read_spectra_from_ms2_file(fn):
    """Return a list of spectrum objects read from an ms2 file."""
    f = file(fn)
    spectra = []
    while True:
        s = spectrum.spectrum()
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
    for m in re.finditer(r'[KR][^P]', sequence):
        yield(m.start() + 1)
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
    "refine, spectrum synthesis" : (bool, "no"),
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
def peptide_mass(peptide_seq, is_N, is_C):
    NC_mods = 0
    if is_N:
        NC_mods += CP.modification_mass[ord('[')]
    if is_C:
        NC_mods += CP.modification_mass[ord(']')]
    return (sum((CP.average_residue_mass[ord(r)]
                 + CP.modification_mass[ord(r)])
                for r in peptide_seq)
            + NC_mods
            + CP.cleave_N_terminal_mass_change
            + CP.cleave_C_terminal_mass_change
            + CP.proton_mass)
# probably moves to C++, mods NYI
# parent mass assumed to be average (CHECK!)
def peptide_mod_mass(peptide_seq, potential_mod_pattern,
                     is_N, is_C):
    return 0



def synthetic_B_spectrum(ion_type, peptide_mass, charge, peptide_seq,
                         peptide_mod_pattern, peptide_is_N):
    m = (CP.cleave_N_terminal_mass_change - formula_mass("H"))
    if peptide_is_N:
        m += CP.modification_mass[ord('[')]
        if False:                       # [ is diff modded
            #m += CP.potential_modification_mass[ord('[')]...
            pass
        
        
                       
                       



def synthetic_spectrum(ion_type, peptide_mass, charge, peptide_seq,
                       peptide_mod_pattern, peptide_is_N, peptide_is_C):
    if ion_type == 'B':
        return synthetic_B_spectrum(ion_type, peptide_mass, charge,
                                    peptide_seq, peptide_mod_pattern,
                                    peptide_is_N)
    elif ion_type == 'Y':
        return synthetic_Y_spectrum(ion_type, peptide_mass, charge,
                                    peptide_seq, peptide_mod_pattern,
                                    peptide_is_C)
    assert False, "unimplemented ion type"


# probably moves to C++, mods NYI
def score_spectrum(peptide_seq, peptide_mod_pattern, peptide_mass,
                   peptide_is_N, peptide_is_C, spectrum, spectrum_mass,
                   spectrum_charge):

    #search_charges = { 1 : [1], 2 : [1,2], 3 : [1,2,3], 4 : [1,2,3,4] }
    search_charges = { 1 : [1], 2 : [1], 3 : [1,2], 4 : [1,2,3] } # WRONG?
    
    for ion_type in ('B', 'Y'):
        for charge in search_charges[spectrum_charge]:
            synth_sp = synthetic_spectrum(ion_type, peptide_mass, charge,
                                          peptide_seq, peptide_mod_pattern,
                                          peptide_is_N, peptide_is_C)

            # convolution score is just sum over charges/ions
            # hyperscore is product of p! over charges/ions (where p is corr
            # peak count) times the convolution score (clipped to FLT_MAX)
            # > blurred!
            peaks, score = correlate_spectra(synth_sp, spectrum)

    # sum peaks, score over ?
    return hyperscore?


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
        error("sequence database is empty")

    db = [ (idno, offset, defline, seq, is_N, is_C,
            list(generate_cleavage_points(seq)))
           for idno, offset, defline, seq, is_N, is_C in db ]
    print "cleavage_points found"

    # read spectra
    spectrum_fn = XTP["spectrum, path"]
    spectra = read_spectra_from_ms2_file(spectrum_fn)
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

    spectrum_mass_index = ([ (sp.mass, sp.charge, n) for n, sp in enumerate(spectra) ]
                           + [ (sp.secondary_mass, sp.secondary_charge, n)
                               for n, sp in enumerate(spectra)
                               if sp.secondary_charge ])
    spectrum_mass_index.sort()

    peptide_count = 0
    peptides_w_candidate_spectra = 0
    for idno, offset, defline, seq, is_N, is_C, cleavage_points in db:
        #print '#', idno, offset, defline, seq, is_N, is_C, cleavage_points
        for begin, end in generate_peptides(seq, cleavage_points,
                                            XTP["scoring, maximum missed cleavage sites"]):
            peptide_count += 1
            peptide_seq = seq[begin:end]
            peptide_mass = peptide_mass(peptide_seq, is_N, is_C)
            # pyro?
            for potential_mod_pattern in generate_mod_patterns(peptide_seq,
                                                               begin, end):
                peptide_mod_mass = (peptide_mass +
                                    peptide_mod_mass(peptide_seq,
                                                     potential_mod_pattern,
                                                     is_N, is_C))
                sp_mass_lb = peptide_mod_mass - XTP["spectrum, parent monoisotopic mass error plus"]
                sp_mass_ub = peptide_mod_mass - XTP["spectrum, parent monoisotopic mass error minus"]
                candidate_spectra_info \
                    = spectrum_mass_index[bisect.bisect_left(spectrum_mass_index, (sp_mass_lb,))
                                          :bisect.bisect_right(spectrum_mass_index, (sp_mass_ub,))]
                if candidate_spectra_info:
                    peptides_w_candidate_spectra += 1
                #print peptide_mod_mass, candidate_spectra_info, \
                #      [ spectra[c[2]] for c in candidate_spectra_info ]

                for c_sp_mass, c_sp_charge, c_sp_n in candidate_spectra_info:
                    sc = score_spectrum(peptide_seq, potential_mod_pattern,
                                        peptide_mod_mass, is_N, is_C,
                                        spectra[c_sp_n], c_sp_mass,
                                        c_sp_charge) 

    print 'generated %s peptides' % peptide_count
    print '          %s peptides have candidate spectra' % peptides_w_candidate_spectra




    # - spectrum preparation tricks:
    # parent mass is a range, defined by errors
    # blurring of mz/I values when converting to integer (why?)

    # search
    # for each sequence in the database(s):
    #    find all position matches for all given (always-applied) motifs
    #    [skipping over '*'s (which are never missed cleaves)]
    #    for each possible cleavage peptide for the sequence:
    #    [subject to a limit on missed cleavages]
    #    [don't bother generating peptides that are two big (or too small?)
    #     for the given spectra, len must be >= 4 (make this a param)]



    # pyro_check: if no N-terminal mod otherwise specified, look for potential
    # mods for these N-terminal residues: Q or C+57 -> loss of ammonia,
    # E -> loss of water


    # don't score a spectrum unless it has a min # of peaks, and at least one
    # from each of ABC and XYZ


    # optional refinement step


    # output results
    # [strip control chars from deflines]




if __name__ == '__main__':
    main()
