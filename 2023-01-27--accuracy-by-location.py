#!/usr/bin/env python3

import sys
from collections import Counter
from Bio.SeqIO.FastaIO import SimpleFastaParser

seq_fname, = sys.argv[1:]
K=40
locs = {}
fwd_locs = {}
rev_locs = {}

def rc(s):
    return "".join({'T':'A',
                    'G':'C',
                    'A':'T',
                    'C':'G',
                    'N':'N'}[x] for x in reversed(s))

with open(seq_fname) as inf:
    (_, seq), = SimpleFastaParser(inf)
    n_kmers = len(seq) - K + 1
    
    for i in range(n_kmers):
        kmer_in = seq[i:i+K]
        for kmer in [kmer_in, rc(kmer_in)]:
            locs[kmer] = i
        fwd_locs[kmer_in] = i
        rev_locs[rc(kmer_in)] = i

loc_observations_full_matches = [0] * n_kmers
loc_observations_any_matches = [0] * n_kmers
loc_read_coverage = [0] * n_kmers

def kmers(seq, offset=0):
    for i in range(len(seq) - K + 1):
        if i < offset: continue
        yield seq[i:i+K]

for title, seq in SimpleFastaParser(sys.stdin):
    total = 0
    matches = 0

    alignment_votes = Counter()
    for pos, kmer in enumerate(kmers(seq)):
        total += 1
        if kmer in locs:
            matches += 1
            loc_observations_any_matches[locs[kmer]] += 1
            if kmer in fwd_locs:
                alignment_votes[locs[kmer] - pos] += 1
            if kmer in rev_locs:
                alignment_votes[locs[kmer] + pos] += 1

    if total == matches:
        for kmer in kmers(seq):
            loc_observations_full_matches[locs[kmer]] += 1

    best_alignment, _ = alignment_votes.most_common(1)[0]
    for pos, _ in enumerate(seq):
        loc = best_alignment + pos
        if loc < len(loc_read_coverage):
            loc_read_coverage[loc] += 1
    
for loc in range(n_kmers):
    print("%d\t%d\t%d\t%d" % (
        loc,
        loc_observations_full_matches[loc],
        loc_observations_any_matches[loc],
        loc_read_coverage[loc]))

    
                
    
