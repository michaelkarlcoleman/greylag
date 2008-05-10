#!/usr/bin/env python

'''This program does the actual search of mass spectra against a sequence
database.  It listens for a connection from greylag-rally, which sends all of
the information necessary to do the search, and returns the search results on
the same connection.  It will continue to process further search requests on
the connection, timing out if none are received within the specified interval.
Then new connections will be accepted, one at a time.  If there are no new
connections within the specified interval, the program exits.

On startup, the program will attempt to listen on a sequence of ports,
depending on the values of --port and --port-count.  If these are 12345 and 4,
for example, it will try ports 12345, 12346, 12347, and 12348, using the first
available.  Since only one process is allowed to listen on each port, this
functions as a form of locking that will prevent more than --port-count
occurrences of this program from running on the host simultaneously.

At most, this program accesses the sequence databases through the filesystem.
'''

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


import cPickle as pickle
import errno
import logging; from logging import debug, info, warning
import optparse
import os
from pprint import pprint, pformat
import re
import socket
import sys

from greylag import *
import cgreylag


# gc possibly harms performance here, so disable it.  gc only matters for
# cycles, which we (hope we) don't create.  See the gc module docs.
import gc; gc.disable()


def error(s, *args):
    "fatal error --> exit with error"
    # if we're unit testing, just throw an exception
    if __name__ != "__main__":
        raise Exception((s + " (fatal error)") % args)
    logging.error(s, *args)
    sys.exit(1)

class ChaseException(Exception): pass

def chase_error(s, *args):
    "error --> disconnect from client"
    logging.error(s, *args)
    raise ChaseException((s + " (disconnecting)") % args)



# Try to drop dead immediately on SIGINT (control-C), instead of normal Python
# KeyboardInterrupt processing, since we may spend long periods of time in
# uninterruptible C++ calls.
try:
    import signal
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    #signal.signal(signal.SIGPIPE, signal.SIG_DFL)
except:
    pass


# name -> value map of processed greylag config parameters
GLP = {}


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
    if GLP["enable_pca_mods"]:
        return [ [('', 0, 0),
                  ('E', -1 * CP.parent_mass_regime[r].water_mass,
                   -1 * CP.fragment_mass_regime[r].water_mass),
                  ('QC', -1 * CP.parent_mass_regime[r].ammonia_mass,
                   -1 * CP.fragment_mass_regime[r].ammonia_mass)]
                 for r in range(len(mass_regimes)) ]
    return [ [('', 0, 0)] for r in range(len(mass_regimes)) ]


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


def set_context_conjuncts(context, mass_regime_index, R_cj):
    context.delta_bag_lookup.clear()
    context.delta_bag_lookup.resize(ord('Z')+1)
    context.delta_bag_delta.clear()
    for n, cj in enumerate(R_cj):
        context.delta_bag_delta.append(cj[6][mass_regime_index][1])
        for r in cj[3]:
            context.delta_bag_lookup[ord(r)] \
                = context.delta_bag_lookup[ord(r)] + (n,)


def search_all(context, mod_limit, mod_conjunct_triples, score_statistics):
    """Search sequence database against searchable spectra."""

    mass_regimes = GLP["mass_regimes"]
    pca_table = get_pca_table(mass_regimes)
    debug("pca_table: %s", pca_table)

    total_evaluation_count = 0

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
                    context.conjunct_index = cji
                    set_context_conjuncts(context, mr_index, R_cj)
                    debug("cj_triple: N=%s C=%s R=%s", N_cj, C_cj, R_cj)
                    for delta_bag in generate_delta_bag_counts(mod_count,
                                                               len(R_cj)):
                        # this clear() avoids a SWIG/STL bug!?
                        context.delta_bag_count.clear()
                        context.delta_bag_count[:] = delta_bag

                        parent_delta = sum(count*r_cj[6][mr_index][0]
                                           for count, r_cj
                                           in zip(delta_bag, R_cj))
                        debug("parent_delta: %s", parent_delta)

                        pmrf = CP.parent_mass_regime[mr_index].fixed_residue_mass
                        fmrf = CP.fragment_mass_regime[mr_index].fixed_residue_mass
                        p_fx = (pmrf[ord('[')]
                                + (N_cj[0][6][mr_index][0] if N_cj else 0)
                                + pmrf[ord(']')]
                                + (C_cj[0][6][mr_index][0] if C_cj else 0)
                                + parent_delta + pca_parent_delta + PROTON_MASS)
                        context.parent_fixed_mass = p_fx
                        f_N_fx = (fmrf[ord('[')]
                                  + (N_cj[0][6][mr_index][1] if N_cj else 0)
                                  + pca_frag_delta)
                        context.fragment_N_fixed_mass = f_N_fx
                        f_C_fx = (fmrf[ord(']')]
                                  + (C_cj[0][6][mr_index][1] if C_cj else 0)
                                  + CP.fragment_mass_regime[mr_index].water_mass)
                        context.fragment_C_fixed_mass = f_C_fx

                        info("MC=%s MR=%s PCA=%s CJ=%s DB=%s", mod_count,
                             mr_index, pca_res, cji, delta_bag)
                        debug("p_fx %s f_N_fx %s f_C_fx %s", p_fx, f_N_fx,
                             f_C_fx)

                        score_statistics.evaluation_count = 0
                        cgreylag.spectrum.search_runs(context,
                                                      score_statistics)
                        total_evaluation_count \
                            += score_statistics.evaluation_count
                        info("  %20s evaluations, this bag",
                             score_statistics.evaluation_count)

    info('%s candidate spectra examined',
         score_statistics.candidate_spectrum_count)
    info('%s total evaluations', score_statistics.evaluation_count)


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


def fix_up_flanking_residues(fasta_db, best_matches):
    """For each peptide match, choose one of the match loci and update the
    N_peptide_flank and C_peptide_flank residues according to the flanking
    residues at that locus.  If there are multiple loci, choose the most
    tryptic.  (Use '-' if there is no flanking residue or it is unknown.)
    """

    locus_sequence_by_name = dict((locusname, seq)
                                  for (locusname, defline, seq, filename)
                                  in fasta_db)

    def check_residue(r):
        return r if r in RESIDUES else '-'

    def count_tryptic_pairs(l):
        return sum(1 for N,C in l if N in ('K', 'R', '-') and C not in ('P',))

    for spectrum_best_matches in best_matches:
        for m in spectrum_best_matches:
            if not m['peptide_sequence']:
                continue
            N_residue = check_residue(m['peptide_sequence'][0])
            C_residue = check_residue(m['peptide_sequence'][-1])
            candidate_flanks = []
            for locus_index in range(len(m['sequence_name'])):
                l_sequence_name = m['sequence_name'][locus_index]
                l_peptide_begin = m['peptide_begin'][locus_index]

                locus_sequence = locus_sequence_by_name[l_sequence_name]
                N_flank = (locus_sequence[l_peptide_begin-1]
                           if l_peptide_begin > 0 else '-')
                peptide_end = l_peptide_begin+len(m['peptide_sequence'])
                C_flank = (locus_sequence[peptide_end]
                           if peptide_end < len(locus_sequence) else '-')
                candidate_flanks.append([(N_flank, N_residue),
                                         (C_residue, C_flank)])

            best_candidate = max(candidate_flanks, key=count_tryptic_pairs)
            m['N_peptide_flank'] = best_candidate[0][0]
            m['C_peptide_flank'] = best_candidate[-1][-1]


def pythonize_swig_object(o, only_fields=None, skip_fields=[]):
    """Creates a pure Python copy of a SWIG object, so that it can be easily
    pickled, or printed (for debugging purposes).  Each SWIG object is
    pythonized as a dictionary, for convenience.  If provided, 'only_fields',
    limits the copy to the list of fields specified.  Otherwise,
    'skip_fields' if given is a list of methods not to include (this helps
    avoid infinite recursion).  Callable sub-objects are skipped.

    >>> pprint(pythonize_swig_object(cgreylag.score_stats(1, 1)))
    {'best_matches': [[{'C_peptide_flank': '-',
                        'N_peptide_flank': '-',
                        'conjunct_index': -1,
                        'mass_regime_index': -1,
                        'mass_trace': [],
                        'pca_delta': 0.0,
                        'peptide_begin': [],
                        'peptide_sequence': '',
                        'predicted_parent_mass': 0.0,
                        'score': 0.0,
                        'sequence_name': [],
                        'spectrum_index': -1}]],
     'candidate_spectrum_count': 0,
     'evaluation_count': 0}

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
            if (not callable(getattr(o, a))
                and (only_fields != None and a in only_fields
                     or only_fields == None and a not in skip_fields)):
                s[a] = pythonize_swig_object(getattr(o, a), only_fields,
                                             skip_fields)
        return s
    return o


def results_dump(fasta_db, score_statistics, searchable_spectra):
    """Return a result dict mapping spectrum names to (spectrum_metadata,
    best_matches) pairs.

    Unneeded fields are stripped.  Some are always removed, while others are
    removed if they have a value matching their default (0 or []).
    Placeholders in the match list (which do not represent actual matches) are
    also removed.
    """

    # The goal of stripping is not so much to reduce the size of the glw
    # files--though it does do this--as to reduce the memory footprint of
    # loading a merged glw file representing an entire run.  It also probably
    # speeds things up a bit.

    r = {}
    spectrum_metadata_fs = set(['id', 'name', 'file_id', 'mass', 'charge',
                                'total_ion_current', 'comparisons'])
    py_s_spectra = pythonize_swig_object(searchable_spectra,
                                         only_fields=spectrum_metadata_fs)
    py_matches = pythonize_swig_object(score_statistics.best_matches,
                                       skip_fields=['spectrum_index'])
    assert len(py_s_spectra) == len(py_matches)

    # FIX: doing this here because it's much easier after pythonizing
    fix_up_flanking_residues(fasta_db, py_matches)

    # note that both strip functions modify their argument
    def strip_meta(meta):
        # these fields are redundant, could strip them
        #del meta['file_id']
        #del meta['name']
        #del meta['charge']
        return meta

    def strip_match(match):
        assert len(match['peptide_begin']) == len(match['sequence_name'])
        if len(match['peptide_begin']) == 0:
            return None
        if match['conjunct_index'] == 0:
            del match['conjunct_index']
        if match['mass_regime_index'] == 0:
            del match['mass_regime_index']
        if len(match['mass_trace']) == 0:
            del match['mass_trace']
        if match['pca_delta'] == 0:
            del match['pca_delta']
        return match

    for sp_metadata, sp_matches in zip(py_s_spectra, py_matches):
        sp_key = sp_metadata['id']
        assert sp_key not in r, "duplicate spectrum id"
        sp_matches_stripped = [ strip_match(m) for m in sp_matches ]
        sp_matches_stripped = [ m for m in sp_matches_stripped if m != None ]
        r[sp_key] = (strip_meta(sp_metadata), sp_matches_stripped)

    return r


def set_parameters(arg, options):
    global GLP
    GLP = arg
    fixed_mod_map = dict((r[3], r) for r in GLP["pervasive_mods"])
    regime_manifest = initialize_spectrum_parameters(options, GLP,
                                                     GLP["mass_regimes"],
                                                     fixed_mod_map)
    mod_conjunct_triples = get_mod_conjunct_triples(GLP["potential_mods"],
                                                    GLP["potential_mod_limit"],
                                                    GLP["mass_regimes"])
    if len(mod_conjunct_triples) >= sys.maxint-8:
        chase_error("too many conjunct triples")
    GLP[">mod_conjunct_triples"] = mod_conjunct_triples

    info("%s unique potential modification conjuncts",
         len(mod_conjunct_triples))
    debug("mod_conjunct_triples (unique):\n%s",
          pformat(mod_conjunct_triples))

    # (cleavage_re, position of cleavage in cleavage_re)
    cleavage_pattern, cleavage_position \
                      = cleavage_motif_re(GLP["cleavage_motif"])
    if cleavage_position == None:
        chase_error("cleavage site '%s' is missing '|'",
                    GLP["protein, cleavage site"])
    cleavage_pattern = re.compile(cleavage_pattern)
    GLP[">cleavage_pattern"] = cleavage_pattern
    GLP[">cleavage_position"] = cleavage_position


def set_sequences(arg):
    assert arg[0] in ('name', 'value')
    if arg[0] == 'name':
        checked = [ (db, file_sha1(db)) for db, cksum in arg[1] ]
        if checked != arg[1]:
            chase_error("database checksum does not match [%s]",
                        (checked, arg[1]))
        # [(locusname, defline, seq, filename), ...]
        fasta_db = list(read_fasta_files([db for db, cksum in arg[1]]))
    else:
        fasta_db = arg[1]
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
    if max_run_length > sys.maxint:
        chase_error("runs longer than %s not yet supported", max_run_length)
    if not db:
        chase_error("no database sequences")

    context = cgreylag.search_context()
    context.nonspecific_cleavage = (GLP["cleavage_motif"] == "[X]|[X]")
    for idno, offset, locusname, defline, seq, seq_filename in db:
        cp = []
        if not context.nonspecific_cleavage:
            cp = list(generate_cleavage_points(GLP[">cleavage_pattern"],
                                               GLP[">cleavage_position"], seq))
        sr = cgreylag.sequence_run(idno, offset, seq, cp, locusname)
        context.sequence_runs.append(sr)
    context.maximum_missed_cleavage_sites = GLP["maximum_missed_cleavage_sites"]
    return context, fasta_db


def make_swig_spectrum(py_spectrum):
    """Given a spectrum tuple (file_index, physical_index, index, name, mass,
    charge, peaks), return an equivalent SWIG spectrum object.
    """
    sp = cgreylag.spectrum(float(py_spectrum[4]), int(py_spectrum[5]))
    sp.file_id, sp.physical_id, sp.id, sp.name = py_spectrum[:4]

    peaks_as_string = py_spectrum[6]
    assert peaks_as_string[-1] == '\n', "malformed spectrum file"
    # ( (mz0, intensity0), (mz1, intensity1), ... )
    peaks = tuple(tuple(float(v) for v in pkline.split(None, 1))
                  for pkline in peaks_as_string.split('\n')[:-1])
    #debug("py_spectrum: %s" % (py_spectrum,))
    #debug("peaks: %s" % (peaks,))
    sp.set_peaks_from_matrix(peaks)
    return sp


def set_spectra(arg):
    # [ (file_index, physical_index, index, name, mass, charge, peaks), ... ]
    py_spectra = arg

    # FIX: how to handle input errors???

    ### create SWIG spectrum objects, for C++ access
    # FIX: are physical_index and index even needed?
    debug("parsing %s spectra", len(py_spectra))
    spectra = [ make_swig_spectrum(sp) for sp in py_spectra ]
    debug("parsed")
    spectra.sort(key=lambda x: x.mass)

    # FIX: send stats back to master?  move common code?
    def peak_statistics(spectra):
        counts = [ len(sp.peaks) for sp in spectra ]
        counts.sort()
        n = len(counts)
        return (counts[0], counts[int(n*0.25)], counts[int(n*0.5)],
                counts[int(n*0.75)], counts[-1], sum(counts) / float(n))

    def print_spectrum_statistics(spectra):
        info("  %s spectra (mass range %s - %s)", len(spectra),
             spectra[0].mass, spectra[-1].mass)
        info("  peak stats: %s..%s..%s..%s..%s (mean=%.2f)"
             % peak_statistics(spectra))

    if spectra:
        info("read spectra:")
        print_spectrum_statistics(spectra)
    else:
        warning("no input spectra")

    # filter and normalize spectra
    for sp in spectra:
        sp.filter_peaks(GLP["TIC_cutoff_proportion"],
                        CP.parent_mass_tolerance_max)
        sp.classify(GLP["intensity_class_count"], GLP["intensity_class_ratio"],
                    GLP["fragment_mass_tolerance"])

    min_psm = GLP["min_parent_spectrum_mass"]
    max_psm = GLP["max_parent_spectrum_mass"]

    # search only spectra that are of good enough quality
    # FIX: also filter by 1 + 2 + 4 rule?
    good_spectra = [ sp for sp in spectra
                     if (len(sp.peaks) >= 10
                         and min_psm <= sp.mass <= max_psm) ]
    noise_spectrums_ids = (set(sp.id for sp in spectra)
                           - set(sp.id for sp in good_spectra))

    if good_spectra:
        info("after filtering:")
        print_spectrum_statistics(good_spectra)
    else:
        info("no spectra pass filters")

    cgreylag.spectrum.set_searchable_spectra(good_spectra)
    return noise_spectrums_ids


def perform_search(state):

    results = { 'modification conjuncts' : GLP[">mod_conjunct_triples"],
                'argv' : sys.argv,
                'noise_spectrum_ids' : state['noise_spectrum_ids'],
                'matches' : {}, 'total comparisons' : 0 }

    spectra_count = len(cgreylag.cvar.spectrum_searchable_spectra)

    if not spectra_count:
        info("no spectra after filtering--search skipped")
        return results

    score_statistics = cgreylag.score_stats(spectra_count,
                                            GLP["best_result_count"])

    info("searching")
    search_all(state['context'], GLP["potential_mod_limit"],
               GLP[">mod_conjunct_triples"], score_statistics)

    info("returning result")
    results['matches'] = results_dump(state['fasta_db'], score_statistics,
                                      cgreylag.cvar.spectrum_searchable_spectra)
    results['total comparisons'] = score_statistics.candidate_spectrum_count

    #debug("sending search results:\n%s" % pformat(results))

    #     if options.estimate_only:
    #         # this factor is just an empirical guess
    #         print ("%.2f generic CPU hours"
    #                % (score_statistics.candidate_spectrum_count / 439.0e6))
    #         return

    return results



# inspired by Beazley's generator talk
def listen_for_connections(port, port_count):
    """Listen on a port from the range [port, port+port_count) and yield a
    series of connections."""
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)

    for p in range(port, port+port_count):
        try:
            s.bind(('', p))             # socket.INADDR_ANY?
            s.listen(5)
            info("listening on port %s", p)
            break
        except socket.error, e:
            if e[0] == errno.EADDRINUSE:
                continue
            raise
    else:
        error("could not listen, all specified ports in use")

    while True:
        client_socket, client_addr = s.accept()
        info("received connection from %s", client_addr)

        # try to better notice reboots, net failures, etc
        try:
            s.setsockopt(socket.SOL_SOCKET, socket.SO_KEEPALIVE,
                         s.getsockopt(socket.SOL_SOCKET,
                                      socket.SO_KEEPALIVE) | 1)
        except socket.error:
            pass

        yield client_socket


# FIX: give cleaner diagnostics on protocol failure?


def reset_state(state):
    state.clear()
    cgreylag.spectrum.set_searchable_spectra([])


def handle_command(options, state, command, arg):
    """Handle command/arg from client, returning (response, response_arg) or
    None.  State that persists across commands is stored in 'state', or on the
    C++ side with cgreylag.spectrum.set_searchable_spectra().
    """
    if command == 'parameters':
        reset_state(state)
        set_parameters(arg, options)
        return None
    elif command == 'sequences':
        context, fasta_db = set_sequences(arg)
        state['context'] = context
        state['fasta_db'] = fasta_db
        return None
    elif command == 'spectra':
        state['noise_spectrum_ids'] = set_spectra(arg)
        return None
    elif command == 'search':
        #currently ignore arg
        results = perform_search(state)
        return 'found', results
    else:
        assert False, "unknown command '%s'" % command


def main(args=sys.argv[1:]):
    parser = optparse.OptionParser(usage=
                                   "usage: %prog [options]",
                                   description=__doc__, version=VERSION)
    pa = parser.add_option
    DEFAULT_PORT = 10078
    pa("--port", dest="port", type="int", default=DEFAULT_PORT,
       help="first listener port [default=%s]" % DEFAULT_PORT)
    DEFAULT_PORT_COUNT = 4
    pa("--port-count", dest="port_count", type="int",
       default=DEFAULT_PORT_COUNT,
       help="number of ports to try [default=%s]" % DEFAULT_PORT_COUNT)
    #DEFAULT_TIMEOUT = 600
    #pa("--timeout", dest="timeout", type="int",
    #   default=DEFAULT_TIMEOUT,
    #   help="inactivity timeout [default=%ss]" % DEFAULT_TIMEOUT)
    pa("-q", "--quiet", action="store_true", dest="quiet", help="no warnings")
    pa("-p", "--show-progress", action="store_true", dest="show_progress",
       help="show running progress")
    pa("-v", "--verbose", action="store_true", dest="verbose",
       help="be verbose")
    pa("-l", "--logfile", dest="logfile",
       help="log to FILE instead of stderr", metavar="FILE")
    pa("--copyright", action="store_true", dest="copyright",
       help="print copyright and exit")
    pa("--debug", action="store_true", dest="debug",
       help="output debugging info")
    pa("--profile", action="store_true", dest="profile",
       help="dump Python profiling output to './greylag-chase.prof.<pid>'")
    (options, args) = parser.parse_args(args=args)

    if options.copyright:
        print __copyright__
        sys.exit(0)

    if len(args) > 0:
        parser.print_help()
        sys.exit(1)

    set_logging(options)

    info("starting on %s", socket.gethostname())

    state = {}
    for client in listen_for_connections(options.port, options.port_count):
        try:
            client_f_to = client.makefile('w')
            client_f_from = client.makefile('r')

            print >> client_f_to, "greylag %s ready" % VERSION
            client_f_to.flush()

            # new client--clear all prior state
            reset_state(state)

            while True:
                command_line = client_f_from.readline()
                if not command_line:
                    info("client closed connection")
                    break
                command, arg_length = command_line.split()
                arg = pickle.loads(client_f_from.read(int(arg_length)))
                try:
                    r = handle_command(options, state, command, arg)
                except ChaseException, e:
                    r = ('error', e)
                if not r:
                    continue
                response, response_arg = r
                if response == 'error':
                    print >> client_f_to, 'error', response_arg
                    break
                assert response == 'found'
                p_response_arg = pickle.dumps(response_arg)
                print >> client_f_to, 'found', len(p_response_arg)
                client_f_to.write(p_response_arg)
                client_f_to.flush()
        except socket.error, e:
            info("closing connection on error [%s]", e)
        except Exception, e:
            try: print >> client_f_to, ('error [%s "%s"]'
                                        % (sys.exc_info()[0], e))
            except: pass
            raise
        finally:
            try: client_f_to.close()
            except: pass
            try: client_f_from.close()
            except: pass
            try: client.close()
            except: pass

    info("exiting")


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
            report_fn = "greylag-chase.prof.%s" % os.getpid()
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
