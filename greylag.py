#!/usr/bin/env python2.4

'''
Partial re-implementation of X!Tandem, with extra features.

'''

__version__ = "$Id$"


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


def filter_spectrum(sp, XTP):
    """Returns a reason (a string) if spectrum to be dropped, else None.
    (This function alters the spectrum peaklist.)
    """

    # corresponds to mspectrumcondition::condition
    # FIX: currently this redundantly analyzes shared peaklists

    #max_parent_charge = parameters.get("spectrum, maximum parent charge", 4)
    #if sp.charge > max_parent_charge:
    #    return 'charge > %s' % max_parent_charge
    # noise suppression check (probably redundant)
    # if sp.mass < 500.0:     # m_fMinMass; "spectrum, minimum parent m+h"
    #    return 'parent mass < "spectrum, minimum parent m+h"'

    # remove_isotopes (XXX: redundant w/below?)
    # >>> removes multiple entries within 0.95 Da of each other, retaining
    # the highest value. this is necessary because of the behavior of some
    # peak finding routines in commercial software
    peakpairs = len(sp.peaks) - 1
    i = 0
    while i < peakpairs:
        p0, p1 = sp.peaks[i:i+2]
        if p1.mz - p0.mz < 0.95:
            peakpairs -= 1
            if p1.intensity > p0.intensity:
                del sp.peaks[i]
            else:
                del sp.peaks[i+1]
        else:
            i += 1

    # remove parent
    # >>> set up m/z regions to ignore: those immediately below the m/z of
    # the parent ion which will contain uninformative neutral loss ions,
    # and those immediately above the parent ion m/z, which will contain
    # the parent ion and its isotope pattern
    parentMZ = 1.00727 + (sp.mass - 1.00727) / sp.charge
    sp.peaks[:] = [ p for p in sp.peaks
                    if (parentMZ >= p.mz
                        and parentMZ - p.mz >= 50.0 / sp.charge
                        or parentMZ < p.mz
                        and p.mz - parentMZ >= 5.0 / sp.charge) ]

    # remove_low_masses
    LOWESTMASS = XTP["spectrum, minimum fragment mz"]
    sp.peaks[:] = [ p for p in sp.peaks if p.mz > LOWESTMASS ]

    # normalize spectrum
    # dynamic_range
    # * use the dynamic range parameter to set the maximum intensity value for the spectrum.
    # * then remove all peaks with a normalized intensity < 1
    DYNAMICRANGE = XTP["spectrum, dynamic range"]
    if not sp.peaks:
        return "no remaining peaks"
    factor = max(1.0, max(p.intensity for p in sp.peaks)) / DYNAMICRANGE
    sp.peaks[:] = [ (p.mz, p.intensity/factor) for p in sp.peaks
                    if p.intensity/factor >= 1.0 ]
    
    # reject the spectrum if there aren't enough peaks 
    MINIMUMPEAKS = XTP["spectrum, minimum peaks"]
    if len(sp.peaks) < MINIMUMPEAKS:
        return "fewer than %s remaining peaks" % MINIMUMPEAKS

    # check is_noise (NYI)
    # * is_noise attempts to determine if the spectrum is simply noise. if the spectrum
    # * does not have any peaks within a window near the parent ion mass, it is considered
    # * noise.
    #limit = sp.mass / sp.charge
    #if sp.charge < 3:
    #    limit = sp.mass - 600.0
    #if not [ True for mz, i in sp.peaks if mz > limit ]:
    #    return "looks like noise"

    # clean_isotopes removes peaks that are probably C13 isotopes
    # FIX: this seems almost identical to remove_isotopes above (kill one off?)
    peakpairs = len(sp.peaks) - 1
    i = 0
    while i < peakpairs:
        p0, p1 = sp.peaks[i:i+2]
        if p1.mz - p0.mz < 1.5:
            peakpairs -= 1
            if p1.intensity > p0.intensity:
                del sp.peaks[i]
            else:
                del sp.peaks[i+1]
        else:
            i += 1

    # keep the MAXPEAKS most intense peaks
    MAXPEAKS = XTP["spectrum, total peaks"]
    sp.peaks[:] = sorted(sp.peaks, key=lambda x: x.intensity, reverse=True)
    if MAXPEAKS < len(sp.peaks):
        del sp.peaks[MAXPEAKS:]
    sp.peaks[:] = sorted(sp.peaks)

    return None


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

aa_sequence = re.compile(r'[^ARNDCQEGHILKMFPSTWYV]+')

def split_sequence_into_aa_runs(sequence, defline):
    """Given an amino-acid sequence, uppercase it and split into runs of amino
    acids suitable for meaningful search (no intervening '*', etc.)."""
    return [ (m.group(), m.start(), defline)
             for m in aa_sequence.finditer(sequence) ]

def read_fasta_files(filenames):
    """Yield (defline, sequence) pairs as read from FASTA files."""
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
            seqs.append(line)
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
# type may be bool, int, float, str, or a tuple of values
# default, as a string value, or None if value must be explicitly specified
# check_fn is an optional function that returns True iff the value is valid
# *1: param required, as has differing defaults based on other params, in
# xtandem

# "ni" means check verifies that only implemented functionality is used
def p_ni_true(x): return x
def p_ni_false(x): return not x
def p_ni_emptystring(x): return x == ""
def p_ni_daltons(x): return x == "Daltons"
def p_positive(x): return x > 0
def p_negative(x): return x < 0
def p_nonnegative(x): return x >= 0

# FIX: add all params
XML_PARAMETER_INFO = {
    "list path, default parameters" : (str, ''),
    "list path, taxonomy information" : (str, None),
    "output, path" : (str, ''),
    "output, message" : (str, '.'),     # ignored
    "protein, C-terminal residue modification mass" : (float, "0.0"),
    "protein, N-terminal residue modification mass" : (float, "0.0"),
    "protein, cleavage C-terminal mass change" : (float, None), # actually calcMass("OH")
    "protein, cleavage N-terminal mass change" : (float, None), # actually calcMass("H")
    "protein, cleavage N-terminal limit" : (int, "100000000", p_positive),
    "protein, cleavage semi" : (bool, "no", p_ni_false),
    "protein, cleavage site" : (("[RK]|{P}",), None),
    "protein, homolog management" : (bool, "no", p_ni_false),
    "protein, modified residue mass file" : (str, "", p_ni_emptystring),
    "protein, taxon" : (str, None),
    "residue, modification mass" : (str, ""),
    "residue, potential modification mass" : (str, ""),
    "residue, potential modification motif" : (str, ""),
    "scoring, a ions" : (bool, "no", p_ni_false),
    "scoring, b ions" : (bool, "yes", p_ni_true),
    "scoring, c ions" : (bool, "no", p_ni_false),
    "scoring, cyclic permutation" : (bool, "no", p_ni_false),
    "scoring, include reverse" : (bool, "no", p_ni_false),
    "scoring, minimum ion count" : (int, None, p_positive), # v-1 -> m_lIonCount!
    "scoring, maximum missed cleavage sites" : (int, None, p_nonnegative),
    "scoring, pluggable scoring" : (bool, "no", p_ni_false),
    "scoring, x ions" : (bool, "no", p_ni_false),
    "scoring, y ions" : (bool, "yes", p_ni_true),
    "scoring, z ions" : (bool, "no", p_ni_false),
    "spectrum, dynamic range" : (float, "100", p_positive),
    "spectrum, fragment mass error" : (float, "0.45", p_positive),
    "spectrum, fragment mass type" : (("average", "monoisotopic"), "monoisotopic"),
    "spectrum, fragment monoisotopic mass error units" : (("Daltons", "ppm"), "Daltons", p_ni_daltons),
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
    "spectrum, use neutral loss window" : (bool, "no", p_ni_false),
    "spectrum, use noise suppression" : (bool, "yes", p_ni_false),
}

# name -> value map of processed XML parameters
XTP = {}


def validate_parameters(parameters):
    """Verify that parameters are valid, have valid values, and correspond to
    currently implemented functionality.  Default values are filled in, and a
    name/value dict returned."""

    pmap = {}
    for name, info in XML_PARAMETER_INFO.iteritems():
        type_ = info[0]
        default = info[1]
        check_fn = len(info) > 2 and info[2] or None

        v = parameters.get(name)
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

    # read sequence dbs
    fasta_db = list(read_fasta_files(taxonomy[XTP["protein, taxon"]]))
    db = [ split_sequence_into_aa_runs(sequence, defline)
           for sequence, defline in fasta_db ]
    print "read %s sequences (%s runs)" % (len(fasta_db), len(db))

    # read spectra
    spectrum_fn = XTP["spectrum, path"]
    spectra = read_spectra_from_ms2_file(spectrum_fn)
    print "read %s spectra" % len(spectra)

    # filter spectra
    # NI: optionally using "contrast angle" (mprocess::subtract)
    if XTP["spectrum, use conditioning"]:
        spectra = [ s for s in spectra if not filter_spectrum(s, XTP) ]
    pprint(spectra)
    print "     %s spectra after filtering" % len(spectra)


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
