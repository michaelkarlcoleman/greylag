#!/usr/bin/env python

"""

This script will accept FASTA formatted input and generate a decoy sequence
for each input sequence.  (Thus, the output will have twice as many sequences
as the input.)  A prefix (e.g. 'SHUFFLED_') will be added to the defline of
each decoy sequence.  Also, all sequences will be re-wrapped.

The script will strip duplicate loci in the input, giving an error if the
duplicates, which have the same locus name, do not also have the same
sequence.

For shuffling, the randomness seed defaults to zero.  Use different values if
you want different shuffles for a given input.  Note that the random number
generator is reseeded after each sequence is generated.  This means that a
particular input sequence will always be mapped to the same shuffled sequence,
independent of what other sequences exist in the input, which is desirable for
some purposes.

"""


from collections import defaultdict
import hashlib
import optparse
import random
import re
import sys

import greylag


def warn(message):
    print >> sys.stderr, "warning: %s" % message

def error(s):
    sys.exit('error: ' + s)


whitespace = re.compile('[ \t]+')

default_wrap=80


class abstract_decoy_maker:
    def make(self, s):
        """Return a decoy sequence of the same length.  (Argument is a list of
        residues, and may be modified by this function.)
        """
        error("abstract method not implemented")

class reverse_decoy_maker(abstract_decoy_maker):
    def make(self, s):
        s.reverse()
        return s

class shuffle_decoy_maker(abstract_decoy_maker):
    def __init__(self, random_seed):
        self.random_seed = random_seed

    def make(self, s):
        random.seed(self.random_seed)
        random.shuffle(s)
        return s

class markov_decoy_maker(abstract_decoy_maker):
    def __init__(self, random_seed, length, original_sequence_file):
        random.seed(random_seed)
        self.length = length

        # for order 0 through self.length:
        # [ length-mer -> subsequent residues -> count, ... ]
        self.transition = [ defaultdict(lambda: defaultdict(int))
                            for i in range(length+1) ]
        for locusname, defline, sequence, filename \
                in greylag.read_fasta_files([original_sequence_file]):
            for order in range(length+1):
                seq = '-' * order + sequence
                for i in xrange(len(sequence)):
                    self.transition[order][seq[i:i+order]][seq[i+order]] += 1

        #import pprint
        #for n, t in enumerate(self.transition):
        #    print '### order', n
        #    for k in sorted(t.keys()):
        #        pprint.pprint((k, dict(t[k])))


    # FIX: this could be cached for better performance
    def _choose_next(self, hist):
        histpairs = hist.items()
        total = sum(x[1] for x in histpairs)
        choice = random.randrange(total)
        cum = 0
        for r, count in histpairs:
            cum += count
            if choice < cum:
                return r
        assert False

    def make(self, s):
        key = '-' * self.length
        result = []
        for i in xrange(len(s)):
            for order in range(self.length, -1, -1):
                k = key[-order:] if order else ''
                if k not in self.transition[order]:
                    continue
                r = self._choose_next(self.transition[order][k])
                break
            else:
                assert False

            result.append(r)
            if key:                     # if self.length==0, key is ''
                key = key[1:] + r
        return result


def add_sixmers(sixmer_set, sequence):
    for i in range(len(sequence) - 6 + 1):
        sixmer_set.add(sequence[i:i+6])


def write_locus(options, decoy_maker, seen, sixmers, locusname, defline,
                sequence):
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

    decoy_sequence = ''.join(decoy_maker.make(list(sequence)))
    decoy_defline = greylag.DEFAULT_DECOY_PREFIX + locusname + ' FALSE POSITIVE'

    add_sixmers(sixmers[0], sequence)
    add_sixmers(sixmers[1], decoy_sequence)

    for d, s in [(defline, sequence), (decoy_defline, decoy_sequence)]:
        if options.no_original and d == defline:
            continue
        print '>' + d
        if options.wrap:
            for start in range(0, len(s), options.wrap):
                print s[start:start+options.wrap]
        else:
            print s


def main():
    parser = optparse.OptionParser(usage="usage: %prog [options] <file>",
                                   description=__doc__)
    parser.add_option("-r", "--reverse", action="store_true",
                      dest="reverse",
                      help="reverse the sequences, instead of shuffling")
    parser.add_option("-m", "--markov", type="int", dest="markov_length",
                      help="generate Markov sequences, with memory of length"
                      " N, instead of shuffling", default=None, metavar="N")
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

    if (len(args) != 1
        or options.markov_length != None and options.reverse
        or options.markov_length != None and options.markov_length < 0):
        parser.print_help()
        sys.exit(1)

    if options.markov_length != None:
        decoy_maker = markov_decoy_maker(options.seed, options.markov_length,
                                         args[0])
    elif options.reverse:
        decoy_maker = reverse_decoy_maker()
    else:
        decoy_maker = shuffle_decoy_maker(options.seed)

    # locus id -> (defline, hash of sequence)
    seen = {}

    # real and decoy 6-mers seen
    sixmers = (set(), set())

    for locusname, defline, sequence, filename \
            in greylag.read_fasta_files([args[0]]):
        write_locus(options, decoy_maker, seen, sixmers,
                    locusname, defline, sequence)

    common_sixmers = sixmers[0] & sixmers[1]
    print >> sys.stderr, ("six-mers: %s real %s decoy %s both"
                          % (len(sixmers[0]) - len(common_sixmers),
                             len(sixmers[1]) - len(common_sixmers),
                             len(common_sixmers)))

if __name__ == '__main__':
    main()
