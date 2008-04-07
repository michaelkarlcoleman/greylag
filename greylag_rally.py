#!/usr/bin/env python

'''This program coordinates the search of mass spectra by sending subsearches
to running instances of greylag-chase, which are expected to be listening on
the specified hosts/ports.

<job-id> is a unique identifier used as a prefix for work files.  The
<configuration-file> contains program options.  The spectra are in the
<ms2-file>s.

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


import asynchat
import asyncore
import ConfigParser
import contextlib
import cPickle as pickle
import errno
import exceptions
import heapq
import logging; from logging import debug, info, warning
import optparse
import os
from pprint import pprint, pformat
import re
import socket
import sys
import time

from greylag import *


# gc possibly harms performance here, so disable it.  gc only matters for
# cycles, which we (hope we) don't create.  See the gc module docs.
import gc; gc.disable()


def error(s, *args):
    "fatal error"
    # if we're unit testing, just throw an exception
    if __name__ != "__main__":
        raise Exception((s + " (fatal error)") % args)
    logging.error(s, *args)
    sys.exit(1)


# name -> value map of processed greylag config parameters
GLP = {}


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

    # The first fragmentation regime should generally be MONO, so that
    # formulaic deltas with '!' do the expected thing.
    if result[0][1] != ('MONO', []):
        raise ValueError("first fragmentation regime was something other than"
                         " 'MONO' with no isotopes--this is almost certainly"
                         " not what was intended")
    return result


def parse_mod_term(s, is_potential=False):
    """Parse a modification term, returning a tuple (sign, mod, fixed_regime,
    residues, description, marker).  The marker character must be printable
    and non-alphabetic.

    >>> parse_mod_term('-C2H3ON!@C')
    (-1, 'C2H3ON', True, 'C', None, None)
    >>> parse_mod_term("42@STY phosphorylation '*'", is_potential=True)
    (1, 42.0, False, 'STY', 'phosphorylation', '*')

    """

    m = re.match(r'^\s*(-|\+)?(([1-9][0-9.]*)|([A-Z][A-Z0-9]*))(!?)'
                 r'@([A-Z]+|\[|\])(\s+([A-Za-z0-9_]+))?'
                 r"(\s+'([[-`!-@{-~])')?\s*$", s)
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
    return (1 if mg[0] != '-' else -1, delta, mg[4] == '!', residues, mg[7],
            mg[9])


def fixed_mod_list(specification):
    """Check and return a list of modification tuples.

    >>> fixed_mod_list('57@C')
    [(1, 57.0, False, 'C', None, None)]
    >>> fixed_mod_list("57@C '*',CH!@N desc")
    [(1, 57.0, False, 'C', None, '*'), (1, 'CH', True, 'N', 'desc', None)]
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
    """Parse s looking for either a term or a parenthesized subexpression.
    Return the pair (tree of modification tuples, unparsed suffix of s).
    """
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
    """Parse s looking for a ','-separated sequence of subexpressions.  Return
    the pair (tree of modification tuples, unparsed suffix of s).
    """
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
    """Parse s looking for a ';'-separated sequence of subexpressions.  Return
    the pair (tree of modification tuples, unparsed suffix of s).
    """

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
    >>> pprint(potential_mod_list('PO3H@STY; C2H2O@KST'))
    [[(1, 'PO3H', False, 'STY', None, None)],
     [(1, 'C2H2O', False, 'KST', None, None)]]
    >>> pprint(potential_mod_list('PO3H@STY, C2H2O@KST'))
    [[(1, 'PO3H', False, 'STY', None, None),
      (1, 'C2H2O', False, 'KST', None, None)]]
    >>> pprint(potential_mod_list('''(PO3H@STY phosphorylation '*';
    ...                               C2H2O@KST acetylation '^';
    ...                               CH2@AKST methylation '#'),
    ...                              O@M oxidation '@'
    ...                           '''))
    [[[[(1, 'PO3H', False, 'STY', 'phosphorylation', '*')],
       [(1, 'C2H2O', False, 'KST', 'acetylation', '^')],
       [(1, 'CH2', False, 'AKST', 'methylation', '#')]],
      (1, 'O', False, 'M', 'oxidation', '@')]]

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
    "mass_regimes" : (mass_regime_list, "MONO"),
    "mass_regime_debug_delta" : (float, 0),
    "pervasive_mods" : (fixed_mod_list, ""),
    "potential_mods" : (potential_mod_list, ""),
    "potential_mod_limit" : (int, 2, p_nonnegative),
    "enable_pca_mods" : (bool, "yes"),
    "charge_limit" : (int, 3, p_positive),
    "min_peptide_length" : (int, 5, p_positive),
    "cleavage_motif" : (str, "[X]|[X]"),
    "maximum_missed_cleavage_sites" : (int, 1, p_nonnegative),
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
        check_fn = p_info[2] if len(p_info) > 2 else None

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


# about 100/s, probably fast enough (just transmitting takes longer?)
def read_spectra_from_ms2_files(spectrum_fns):
    """Read ms2 files, yielding tuples (file_index, physical_index, index,
    name, mass, charge, peaks).  The indexes are 0-based.  Mass and charge are
    returned as strings, and peaks is the newline-embedded peak block from the
    file.
    """
    header_re = re.compile(r'^:', re.MULTILINE)

    physical_index = 0
    index = 0
    pending_headers = []
    for file_index, spfn in enumerate(spectrum_fns):
        with open(spfn) as specfile:
            contents = specfile.read()
        for block in header_re.split(contents)[1:]:
            b = block.split('\n', 2)
            if len(b) != 3:
                error("malformed ms2 file '%s' (malformed header)" % spfn)
            name, h_line1, peaks = b
            h = h_line1.split(' ', 1)
            if len(h) != 2:
                error("malformed ms2 file '%s' (malformed mass/charge header)"
                      % spfn)
            mass, charge = h
            pending_headers.append((name, mass, charge))
            if peaks:
                for name, mass, charge in pending_headers:
                    yield (file_index, physical_index, index,
                           name, mass, charge, peaks)
                    index += 1
                pending_headers = []
                physical_index += 1
        if pending_headers:
            error("malformed ms2 file '%s' (headers at end?)" % spfn)


def check_ms2_files(spectrum_fns):
    ms2_re = re.compile(r'^((:.*\n[1-9]\d*\.\d+ [1-9]\d*\r?\n)+(\d+\.\d+ \d(\.\d+)?\r?\n)+)+$')

    for file_index, spfn in enumerate(spectrum_fns):
        with open(spfn) as specfile:
            contents = specfile.read()
        info("checking '%s' (%s bytes)" % (spfn, len(contents)))
        if not ms2_re.match(contents):
            error("invalid ms2 file '%s'" % spfn)


# def check_ms2_spectra(spectra):
#     mass_re = re.compile(r'[1-9]\d*\.\d+')
#     charge_re = re.compile(r'[1-9]\d*')
#     peaks_re = re.compile(r'([1-9]\d*\.\d+\ \d+(\.\d+)?\n)+')

#     for file_index, physical_index, index, name, mass, charge, peaks in spectra:
#         if not mass_re.match(mass):
#             error("invalid ms2 mass '%s'" % mass)
#         if not charge_re.match(charge):
#             error("invalid ms2 charge '%s'" % charge)
#         if not peaks_re.match(peaks):
#             error("invalid ms2 peaks '%s'" % peaks)


# FIX: optionally print peak statistics?  (slow?)
# def peak_statistics(spectra):
#     counts = [ len(sp.peaks) for sp in spectra ]
#     counts.sort()
#     n = len(counts)
#     return (counts[0], counts[int(n*0.25)], counts[int(n*0.5)],
#             counts[int(n*0.75)], counts[-1], sum(counts) / float(n))
#
# def print_spectrum_statistics(spectra):
#     info("  %s spectra (mass range %s - %s)", len(spectra),
#          spectra[0].mass, spectra[-1].mass)
#     info("  peak stats: %s..%s..%s..%s..%s (mean=%.2f)"
#          % peak_statistics(spectra))


class chase_client(asynchat.async_chat):

    # heapq of (deathtime, host, port) for connections that have died or
    # timed-out.  deathtime is the time of death.  It's up to external code to
    # recreate after a suitable wait, if desired.
    dead_clients = []

    # cached so that we don't have thousands of identical copies
    pickled_parameters = None
    pickled_sequences = None

    # heapq of (submitcount, lastsubmittime, spectrum) for all spectra not
    # yet searched.  submitcount is the number of time this spectrum has been
    # submitted to a chasers.  This can be greater than one because we plan
    # for any chaser to fail.  We don't decrement submitcount--it's just a
    # hint as to which spectra should be searched next.  lastsubmittime is
    # provided so that we can submit the least-recently submitted spectrum.
    spectrum_queue = None

    # number of spectra for which results not yet received
    spectra_to_go = None

    # The number of spectra submitted to clients for a search command.  This
    # can be adjusted upwards for efficiency.
    spectrum_batch_size = 1

    # Reply for each searched spectrum, indexed by spectrum.id (1-based).
    # None if no reply yet.
    results = []

    @classmethod
    def set_parameters(self, parameters):
        assert self.pickled_parameters == None, "no change once set"
        self.pickled_parameters = pickle.dumps(parameters)

    @classmethod
    def set_sequences(self, sequences):
        assert self.pickled_sequences == None, "no change once set"
        self.pickled_sequences = pickle.dumps(sequences)

    @classmethod
    def set_spectra(self, spectra):
        """Fix the list of spectra to be searched, which must be done before
        any clients are created (and cannot be changed later).
        """
        assert self.spectrum_queue == None, "no change once set"
        self.spectrum_queue = [ (0, 0, spectrum) for spectrum in spectra ]
        heapq.heapify(self.spectrum_queue)
        self.results = [None] * (len(spectra) + 1)
        self.spectra_to_go = len(spectra)


    def __init__(self, host, port):
        asynchat.async_chat.__init__(self)
        assert self.pickled_parameters and self.pickled_sequences
        assert self.spectrum_queue
        self.host = host
        self.port = port

        self.ibuffer = []
        self.set_terminator('\n')

        self.banner = None

        self.submit_count = 0           # number of submits to this client
        self.create_socket(socket.AF_INET, socket.SOCK_STREAM)
        self.set_keep_alive()
        self.connect((host, port))

    # override asynchat version to use a reasonable buffer size
    def push(self, data):
        self.producer_fifo.push(
            asynchat.simple_producer(data, buffer_size=self.ac_out_buffer_size))
        self.initiate_send()

    def __repr__(self):
        return "<%s connected to %s:%s>" % (self.__class__.__name__,
                                            self.host, self.port)


    def set_keep_alive(self):
        # try to better notice reboots, net failures, etc
        try:
            self.socket.setsockopt(
                socket.SOL_SOCKET, socket.SO_KEEPALIVE,
                self.socket.getsockopt(socket.SOL_SOCKET,
                                       socket.SO_KEEPALIVE) | 1
                )
        except socket.error:
            pass


    def _send_command(self, command, pickled_argument):
        self.push(command + ' ')
        self.push("%s\n" % len(pickled_argument))
        self.push(pickled_argument)

    def _submit_search(self):
        "send a search to client"
        self.submit_count += 1
        submit_spectra = []

        # start slow, to get chasers going immediately, then increase batch
        # size to increase master efficiency (and level load)
        #batch_size = min(self.spectrum_batch_size,
        #                 4**min(self.submit_count-1,10))
        batch_size = self.spectrum_batch_size

        if self.spectrum_queue:
            # on each retry round, divide batch size by four
            batch_size = max(1, (batch_size
                                 / 2**min(64, 2*self.spectrum_queue[0][0])))

        while self.spectrum_queue and len(submit_spectra) < batch_size:
            count, _submittime, spectrum = heapq.heappop(self.spectrum_queue)
            # spectrum = (file_index, physical_index, index,
            #             name, mass, charge, peaks)
            if self.results[spectrum[2]] != None:
                continue                # already have search result
            heapq.heappush(self.spectrum_queue,
                           (count+1, time.time(), spectrum))
            submit_spectra.append(spectrum)
        if submit_spectra:
            self._send_command('spectra', pickle.dumps(submit_spectra))
            self._send_command('search', pickle.dumps(None))
        else:
            self.close()
        debug("submitting %s spectra (%s not yet retired)"
              % (len(submit_spectra), len(self.spectrum_queue)))

    def _receive_response(self):
        r = ''.join(self.ibuffer)
        self.ibuffer = []
        return r


    def handle_connect(self):
        debug("connecting to %s:%s" % (self.host, self.port))

    def handle_expt(self):
        pass

    def handle_error(self):
        #self.handle_close()
        exc, why, _traceback = sys.exc_info()
        self.connected = False
        if exc == exceptions.KeyboardInterrupt: # FIX: works?
            error("received keyboard interrupt")
        if exc == socket.error:
            if why[0] == errno.ECONNREFUSED:
                debug("no chaser at %s:%s (connection refused)"
                      % (self.host, self.port))
            else:
                info("network error on connection %s:%s (%s %s)"
                     % (self.host, self.port, exc, why))
                debug("  traceback: %s" % asyncore.compact_traceback()[3])
        else:
            info("unexpected error on connection %s:%s (%s %s)"
                 % (self.host, self.port, exc, why))
            info("  traceback: %s" % asyncore.compact_traceback()[3])


    def handle_close(self):
        self.__class__.dead_clients.append((time.time(), self.host, self.port))
        asynchat.async_chat.handle_close(self)


    def collect_incoming_data(self, data):
        self.ibuffer.append(data)

    def found_terminator(self):
        if self.banner == None:
            self.banner = self._receive_response()
            banner_words = self.banner.split()
            if banner_words[0] != 'greylag':
                warning("%s:%s is not a greylag-chase" % (self.host, self.port))
                self.handle_close()
                return
            if banner_words[2:3] != ['ready']:
                debug("greylag-chase %s:%s is serving %s"
                      % (self.host, self.port, banner_words[3:4]))
                self.handle_close()
                return
            # now we can start talking
            self._send_command('parameters', self.pickled_parameters)
            self._send_command('sequences', self.pickled_sequences)
            self._submit_search()
            return

        if self.get_terminator() == '\n':
            reply = self._receive_response()
            reply_words = reply.split()
            if reply_words[0] == 'error':
                warning("%s:%s gave '%s'" % (self.host, self.port, reply))
                self.handle_close()
                return

            assert reply_words[0] == 'found' and len(reply_words) == 2
            self.set_terminator(int(reply_words[1]))
            return

        result = pickle.loads(self._receive_response())

        #debug("received: %s" % result)
        for noise_sp_id in result['noise_spectrum_ids']:
            if self.__class__.results[noise_sp_id] != None:
                warning("got differing filter results [id=%s]" % noise_sp_id)
                continue
            self.__class__.spectra_to_go -= 1
            self.__class__.results[noise_sp_id] = False

        for sp_id, sp_matches in result['matches'].iteritems():
            if self.__class__.results[sp_id] != None:
                if self.__class__.results[sp_id] != False:
                    pass
                    # FIX
                    #warning("got differing results [id=%s]" % sp_id)
                    #debug("r0:\n%s\nr1:\n%s\n"
                    #      % (self.__class__.results[sp_id],
                    #         sp_matches))
                else:
                    warning("got differing filter results [id=%s]" % sp_id)
            else:
                self.__class__.spectra_to_go -= 1
            self.__class__.results[sp_id] = sp_matches

        self.set_terminator('\n')
        self._submit_search()


def main(args=sys.argv[1:]):
    parser = optparse.OptionParser(usage=
                                   "usage: %prog [options]"
                                   " <configuration-file> <ms2-file>...",
                                   description=__doc__, version=VERSION)
    pa = parser.add_option
    DEFAULT_HOSTFILE = "/etc/greylag/hosts"
    pa("--hostfile", dest="hostfile", default=DEFAULT_HOSTFILE,
       help="file listing host:port locations where greylag-chase workers are"
       " listening.  [default '%s']" % DEFAULT_HOSTFILE)
    pa("--job-id", dest="job_id", default="unknown",
       help="used to generate unique output filenames [default 'unknown']")
    #pa("--estimate", action="store_true", dest="estimate_only",
    #   help="just estimate the time required for the search")
    pa("-v", "--verbose", action="store_true", dest="verbose",
       help="be verbose")
    pa("-q", "--quiet", action="store_true", dest="quiet", help="no warnings")
    pa("-l", "--logfile", dest="logfile",
       help="log to FILE instead of stderr", metavar="FILE")
    pa("-P", "--parameter", dest="parameters", action="append",
       default=[],
       help="override a parameter in <configuration-file>, may be used"
       " multiple times", metavar="NAME=VALUE")
    pa("--copyright", action="store_true", dest="copyright",
       help="print copyright and exit")
    pa("--debug", action="store_true", dest="debug",
       help="output debugging info")
    pa("--profile", action="store_true", dest="profile",
       help="dump Python profiling output to './greylag-rally.prof.<pid>'")
    (options, args) = parser.parse_args(args=args)

    if options.copyright:
        print __copyright__
        sys.exit(0)

    if len(args) < 2:
        parser.print_help()
        sys.exit(1)

    configuration_fn = args[0]
    spectrum_fns = args[1:]

    if (any(True for f in spectrum_fns if not f.endswith('.ms2'))
        or any(True for p in options.parameters if '=' not in p)):
        parser.print_help()
        sys.exit(1)

    set_logging(options)

    info("starting on %s", socket.gethostname())

    result_fn = 'chase_%s.glw' % options.job_id

    # read chaser (host, port) listener list
    try:
        with open(options.hostfile) as hostf:
            hosts = [ l.split(':', 1) for l in hostf ]
            hosts = [ (host, int(port)) for host, port in hosts ]
    except ValueError:
        error("invalid or empty host line in '%s'" % options.hostfile)

    if not hosts:
        error("no valid search hosts specified")

    # check spectrum basename uniqueness, as corresponding sqt files will be
    # in a single directory
    base_spectrum_fns = [ os.path.basename(fn) for fn in spectrum_fns ]
    if len(base_spectrum_fns) != len(set(base_spectrum_fns)):
        error("base spectrum filenames must be unique")

    options.parameters = [ p.split('=', 1) for p in options.parameters ]

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
    global GLP
    GLP = validate_parameters(parameters)

    # test parse cleavage motif here, so we notice errors early
    cleavage_motif_re(GLP["cleavage_motif"])

    fixed_mod_map = dict((r[3], r) for r in GLP["pervasive_mods"])
    # FIX: is this even necessary?
    regime_manifest = initialize_spectrum_parameters(options, GLP,
                                                     GLP["mass_regimes"],
                                                     fixed_mod_map)
    info("reading sequence databases")
    # read sequence dbs
    databases = GLP["databases"].split()
    # [(locusname, defline, seq, filename), ...]
    fasta_db = list(read_fasta_files(databases))

    # This is what we'll send to chasers.  Either (1) the already parsed
    # contents of the databases, which is slow to send, or (2) a list of
    # (database name, sha1 checksums), which requires that the files exist on
    # the chaser hosts.
    db_representation = ("value", fasta_db)
    if 1:                               # FIX: make this an option?
        db_representation = ("name",
                             [ (db, file_sha1(db)) for db in databases ])

    info("FIX checking spectra")
    #check_ms2_files(spectrum_fns)
    info("checked")

    info("reading spectra")
    spectra = list(read_spectra_from_ms2_files(spectrum_fns))

    if spectra:
        info("read %s spectra" % len(spectra))
        # print_spectrum_statistics(spectra)
    else:
        error("no input spectra")

    # FIX: gather these statistics from chasers?
    #if spectra:
    #    info("after filtering:")
    #    print_spectrum_statistics(spectra)

    score_statistics = [] * len(spectra)

    # now connect to chasers and loop until done (or timeout)
    chase_client.set_parameters(GLP)
    chase_client.set_sequences(db_representation)
    chase_client.set_spectra(spectra)

    # Trivial strategy: Hand each client 1/10th of their likely workload in
    # each batch.  If client crashes, we haven't lost too much work.
    # Making this value larger keeps polling efficient and lets us handle a
    # large number of clients.  If it's too large, though, we're need more RAM
    # for temporary storage of pythonized spectra.
    chase_client.spectrum_batch_size = min(100, max(1, int(0.1 * len(spectra)
                                                           / len(hosts))))

    # FIX: this is actually faster (?)
    #chase_client.spectrum_batch_size = 1

    # adjust for performance
    chase_client.ac_in_buffer_size = 62000
    chase_client.ac_out_buffer_size = 62000
    #chase_client.ac_in_buffer_size = 620000
    #chase_client.ac_out_buffer_size = 620000

    if spectra:
        for host, port in hosts:
            chase_client(host, port)
    else:
        warning("no spectra after filtering")

    # retry dead clients after 60s
    retryafter = 60
    #retryafter = 10

    start_time = time.time()
    laststatus = 0
    while True:
        asyncore.loop(count=1, timeout=retryafter, use_poll=True)
        if not chase_client.spectrum_queue:
            break

        now = time.time()

        if now - laststatus >= 10:
            eta_minutes = ((now - start_time)
                           / (len(spectra) - chase_client.spectra_to_go + 0.1)
                           * chase_client.spectra_to_go / 60.0)
            info("%s spectra to search, %s chasers, ETA %dm"
                 % (chase_client.spectra_to_go, len(asyncore.socket_map),
                    int(eta_minutes)))
            laststatus = now

        while (chase_client.dead_clients
               and chase_client.dead_clients[0][0] < now - retryafter):
            debug("died %s restart %s" % (chase_client.dead_clients[0][0], now))
            dt, host, port = heapq.heappop(chase_client.dead_clients)
            chase_client(host, port)
        if not asyncore.socket_map and chase_client.spectrum_queue:
            debug("sleeping")
            time.sleep(max(0, retryafter - (now
                                            - chase_client.dead_clients[0][0])))
        # failsafe, better performance?
        #time.sleep(0.1)


    # FIX: could close all client sockets here--worth doing?

    # FIX
    #total_comparisons = sum(ss.candidate_spectrum_count
    #                        for ss in score_statistics)
    #if options.estimate_only:
    #    # divisor is just a slightly informed guess
    #    print ("%.2f generic CPU hours" % (total_comparisons / 439.0e6))
    #    return

    info("writing result file '%s'", result_fn)

    # sp = (file_index, physical_index, index, name, mass, charge, peaks)
    result_matches = dict(((sp[0], sp[3]), chase_client.results[sp[2]])
                          for sp in spectra if chase_client.results[sp[2]])

    mod_conjunct_triples = get_mod_conjunct_triples(GLP["potential_mods"],
                                                    GLP["potential_mod_limit"],
                                                    GLP["mass_regimes"])

    d = { 'version' : VERSION,
          'matches' : result_matches,
          'total comparisons' : 0, # total_comparisons,
          'spectrum files' : base_spectrum_fns,
          'databases' : databases,
          'parameters' : GLP,
          'mass regime atomic masses' : MASS_REGIME_ATOMIC_MASSES,
          'mass regime manifest' : sorted(regime_manifest),
          'proton mass' : PROTON_MASS,
          'modification conjuncts' : mod_conjunct_triples,
          'argv' : sys.argv }
    with contextlib.closing(open(result_fn, 'w')) as result_file:
        pk = pickle.Pickler(result_file, pickle.HIGHEST_PROTOCOL)
        pk.fast = 1                     # stipulate no circular references
        pk.dump(d)

    info("finished")


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
            report_fn = "greylag-rally.prof.%s" % os.getpid()
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
