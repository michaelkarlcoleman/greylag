#!/usr/bin/env python

"""

This script will accept FASTA formatted input and generate a shuffled sequence
for each input sequence.  (Thus, the output will have twice as many sequences
as the input.)  The prefix 'SHUFFLED_' will be added to the defline of each
shuffled sequence.  Also, all sequences will be re-wrapped.

The script will strip duplicate loci in the input, giving an error if the
duplicates, which have the same locus name, do not also have the same
sequence.

The randomness seed defaults to zero.  Use different values if you want
different shuffles for a given input.  Note that the random number generator
is reseeded after each sequence is generated.  This means that a particular
input sequence will always be mapped to the same shuffled sequence,
independent of what other sequences exist in the input, which is desirable for
some purposes.

"""


import hashlib
import optparse
import random
import re
import sys

import greylag


# change docstring above if this is changed!
SHUFFLE_PREFIX = 'SHUFFLED_'


def warn(message):
    print >> sys.stderr, "warning: %s" % message

def error(s):
    sys.exit('error: ' + s)


whitespace = re.compile('[ \t]+')

default_wrap=80


def main():
    parser = optparse.OptionParser(usage="usage: %prog [options] <file>",
                                   description=__doc__)
    parser.add_option("-r", "--reverse", action="store_true",
                      dest="reverse",
                      help="reverse the sequences, instead of shuffling")
    parser.add_option("-s", "--seed", type="int", dest="seed",
                      help="seed for randomness [default 0]",
                      default=0, metavar="N")
    parser.add_option("-n", "--no-original", action="store_true",
                      dest="no_original",
                      help="don't output original sequences")
    parser.add_option("-v", "--verbose", action="store_true",
                      dest="verbose", help="be verbose")
    parser.add_option("-w", "--wrap", dest="wrap", type="int",
                      default=default_wrap,
                      help="wrap sequence to specified width"
                      " [default %s, 0 means don't wrap at all]" % default_wrap,
                      metavar="COLUMNS")
    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.print_help()
        sys.exit(1)

    random.seed(options.seed)

    # locus id -> (defline, hash of sequence)
    seen = {}

    for locusname, defline, sequence, filename \
            in greylag.read_fasta_files([args[0]]):
        write_locus(options, seen, locusname, defline, sequence)


def write_locus(options, seen, locusname, defline, sequence):
    h = hashlib.sha1()
    h.update(sequence)
    sequence_hash = h.digest()

    if locusname in seen:
        seen_defline, seen_hash = seen[locusname]
        if seen_hash != sequence_hash:
            error("differing sequence for locus '%s'" % locusname)
        if options.verbose and seen_defline != defline:
            warn("differing deflines for locus '%s'" % locusname)
        return
    seen[locusname] = (defline, sequence_hash)

    s_list = list(sequence)
    random.seed(options.seed)
    if options.reverse:
        s_list.reverse()
    else:
        random.shuffle(s_list)
    shuffle_sequence = ''.join(s_list)
    shuffle_defline = SHUFFLE_PREFIX + locusname + ' FALSE POSITIVE'
    for d, s in [(defline, sequence), (shuffle_defline, shuffle_sequence)]:
        if options.no_original and d == defline:
            continue
        print '>' + d
        if options.wrap:
            for start in range(0, len(s), options.wrap):
                print s[start:start+options.wrap]
        else:
            print s


if __name__ == '__main__':
    main()
