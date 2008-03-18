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


import fileinput
import optparse
import random
import re
import sha
import sys


# change docstring above if this is changed!
SHUFFLE_PREFIX = 'SHUFFLED_'


def warn(message):
    print >> sys.stderr, "warning: %s [at %s:%s]" \
          % (message, fileinput.filename(), fileinput.filelineno())
def error(s):
    sys.exit('error: ' + s)


whitespace = re.compile('[ \t]+')

default_wrap=80


def main():
    parser = optparse.OptionParser(usage="usage: %prog [options] [<file>...]",
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

    random.seed(options.seed)

    defline = None
    seqs = []
    # locus id -> (defline, hash of sequence)
    seen = {}
    for line in fileinput.input(args):
        line = line.strip()
        if line[:1] == '>':
            if defline:
                out(defline, seqs, options, seen)
            elif seqs:
                warn("discarding sequence prior to initial defline")
            defline = line
            seqs = []
        else:
            seqs.append(re.sub(whitespace, '', line))
    if defline:
        out(defline, seqs, options, seen)


def out(defline, seqs, options, seen):
    sequence = ''.join(seqs)
    sequence_hash = sha.new(sequence).digest()
    locus_id = re.split(r'[ 	]+', defline, 1)[0]
    if locus_id in seen:
        seen_defline, seen_hash = seen[locus_id]
        if seen_hash != sequence_hash:
            error("differing sequence for locus '%s'" % locus_id)
        if options.verbose and seen_defline != defline:
            warn("differing deflines for locus '%s'" % locus_id)
        return
    seen[locus_id] = (defline, sequence_hash)

    s_list = list(sequence)
    random.seed(options.seed)
    if options.reverse:
        s_list.reverse()
    else:
        random.shuffle(s_list)
    shuffle_sequence = ''.join(s_list)
    shuffle_defline = '>' + SHUFFLE_PREFIX + locus_id[1:] + ' FALSE POSITIVE'
    for d, s in [(defline, sequence), (shuffle_defline, shuffle_sequence)]:
        if options.no_original and d == defline:
            continue
        print d
        if options.wrap:
            for start in range(0, len(s), options.wrap):
                print s[start:start+options.wrap]
        else:
            print s


if __name__ == '__main__':
    main()
