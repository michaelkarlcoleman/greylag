#!/usr/bin/env python

"""Create a censor file for the given FASTA sequence database.  This file
   notes subsequences that are redundant, and therefore need not be searched.
   The locations of the censored subsequences are noted so that correct search
   output can still be generated.  The censor file is a gzip'ed Python pickle
   file with suffix '.censor'.
"""

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


import contextlib
import cPickle
import gzip
import optparse
import os.path
import re
import sys

from greylag import VERSION

##FIX
import greylag_chase


def error(s):
    sys.exit('error: ' + s)


def main(args=sys.argv[1:]):
    parser = optparse.OptionParser(usage=
                                   "usage: %prog [options] <fasta-file>",
                                   description=__doc__)
    pa = parser.add_option
    pa("--copyright", action="store_true", dest="copyright",
       help="print copyright and exit")
    pa("--version", action="store_true", dest="version",
       help="print version and exit")
    (options, args) = parser.parse_args(args=args)

    if options.copyright:
        print __copyright__
        sys.exit(0)
    if options.version:
        print VERSION
        sys.exit(0)

    if len(args) != 1:
        parser.print_help()
        sys.exit(1)


    # [(locusname, defline, seq, filename), ...]
    fasta_db = list(greylag_chase.read_fasta_files([args[0]]))
    # [(idno, offset, locusname, defline, seq, seq_filename), ...]
    db = []
    for idno, (locusname, defline, sequence, filename) in enumerate(fasta_db):
        db.extend(greylag_chase.split_sequence_into_aa_runs(idno, locusname,
                                                            defline, sequence,
                                                            filename))
    db_residue_count = sum(len(dbi[3]) for dbi in db)

    print "# read %s sequences (%s runs, %s residues)" % (len(fasta_db), len(db),
                                                          db_residue_count)

    # each node within the tree is pair: the first element is a 26-ary list
    # holding child nodes; the second is a list of annotations
    LETTERS = 26
    censor_tree = ([None] * LETTERS, [])

    # minimum interesting peptide length (for statistics)
    ## FIX
    MIN_PEPTIDE_LENGTH = 0
    MAX_PEPTIDE_LENGTH = 10

    # statistics
    peptides_to_search = 0
    peptides_to_skip = 0
    #  residues in the above peptides
    residues_to_search = 0
    residues_to_skip = 0

    for idno, offset, locusname, defline, seq, seq_filename in db:
        sys.stderr.write('.')
        for start in xrange(len(seq)):
            suffix = seq[start:]
            sstart = 0
            ctp = censor_tree

            while suffix[sstart:]:
                # walk the previously seen part of suffix
                i = ord(suffix[sstart]) - ord('A')
                if not ctp[0][i]:
                    break       # prefix of suffix unseen
                ctp = ctp[0][i]
                sstart += 1
                #print 'in %s skipping %s' % (suffix, suffix[:sstart])
                if sstart >= MIN_PEPTIDE_LENGTH:
                    peptides_to_skip += 1
                    residues_to_skip += sstart

            while suffix[sstart:]:
                # add the rest of suffix (if any) to the censor tree
                i = ord(suffix[sstart]) - ord('A')
                assert not ctp[0][i]
                if sstart <= MAX_PEPTIDE_LENGTH:
                    ctp[0][i] = ([None] * LETTERS, [])
                    ctp = ctp[0][i]
                sstart += 1
                #print 'in %s adding %s' % (suffix, suffix[:sstart])
                if sstart >= MIN_PEPTIDE_LENGTH:
                    peptides_to_search += 1
                    residues_to_search += sstart

    #from pprint import pprint
    #pprint(censor_tree)

    print ('# search %s (%s residues), skip %s (%s residues)'
           % (peptides_to_search, residues_to_search, peptides_to_skip,
              residues_to_skip))
    print '# average search length %.2f' % (float(residues_to_search)
                                            / peptides_to_search)

    print '# residue search proportion %.3f' % (float(residues_to_search)
                                                / (residues_to_search
                                                   + residues_to_skip))

    sys.exit()


    for fn in args:
        with open(fn) as specfile:
            contents = specfile.read()
        specnames = set()
        with contextlib.closing(gzip.open(fn + '.idx', 'w')) as idx:
            ms = [ m for m in re.finditer(r'^:.*$', contents, re.MULTILINE) ]
            specnames = [ m.group() for m in ms ]
            if len(set(specnames)) < len(ms):
                error("duplicate spectrum names not allowed")
            if specnames != sorted(specnames):
                error("spectra must be ordered by name")
            offsets = [ m.start() for m in ms ]
            cPickle.dump({ 'offsets' : offsets,
                           'file size' : os.path.getsize(fn) },
                         idx, cPickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    main()
