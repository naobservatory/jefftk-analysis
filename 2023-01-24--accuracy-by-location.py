#!/usr/bin/env python3

import sys
from Bio.SeqIO.FastaIO import SimpleFastaParser

seq_fname, = sys.argv[1:]
K=40
locs = {}

def rc(s):
    return "".join({'T':'A',
                    'G':'C',
                    'A':'T',
                    'C':'G',
                    'N':'N'}[x] for x in reversed(s))

with open(seq_fname) as inf:
    seq, = inf
    seq = seq.strip()

    n_kmers = len(seq) - K + 1
    
    for i in range(n_kmers):
        kmer_in = seq[i:i+K]
        for kmer in [kmer_in, rc(kmer_in)]:
            locs[kmer] = i

loc_observations_full_matches = [0] * n_kmers
loc_observations_any_matches = [0] * n_kmers
loc_transitions = [0] * (n_kmers + K)

def kmers(seq, offset=0):
    for i in range(len(seq) - K + 1):
        if i < offset: continue
        yield seq[i:i+K]

for title, seq in SimpleFastaParser(sys.stdin):
    total = 0
    matches = 0

    for prev_kmer, kmer, next_kmer in zip(
            kmers(seq, offset=0),
            kmers(seq, offset=1),
            kmers(seq, offset=2)):

        if kmer in locs and prev_kmer not in locs:
            loc_transitions[locs[kmer]] += 1
        if kmer in locs and next_kmer not in locs:
            loc_transitions[locs[kmer] + K] += 1
            
    for kmer in kmers(seq):
        total += 1
        if kmer in locs:
            matches += 1
            loc_observations_any_matches[locs[kmer]] += 1

    if total == matches:
        for kmer in kmers(seq):
            loc_observations_full_matches[locs[kmer]] += 1

for loc in range(n_kmers):
    print("%d\t%d\t%d\t%d" % (
        loc,
        loc_observations_full_matches[loc],
        loc_observations_any_matches[loc],
        loc_transitions[loc]))

    
                
    
