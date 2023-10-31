#!/usr/bin/env python3

import sys
import regex
from Bio.SeqIO.QualityIO import FastqGeneralIterator

COLOR_RED_BOLD = '\x1b[1;31m'
COLOR_GREEN_BOLD = '\x1b[1;32m'
COLOR_BLUE_BOLD = '\x1b[1;34m'
COLOR_END = '\x1b[0m'

fastq1, fastq2 = sys.argv[1:]

reverse_reads = []

# These are the beginnings of the adapters AdapterRemoval2 identified
# automatically, and they're also the Illumina TruSeq DNA and RNA combinatorial
# dual indexing adapters.
fwd = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
rev = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

def match(haystack, needle):
    best_pos = None
    best_score = 0
    for pos in range(len(haystack) - len(needle) + 1): 
        score = 0
        for h, n in zip(haystack[pos:], needle):
            if h == n:
                score += 1
        if score > best_score:
            best_pos = pos
            best_score = score

    return best_pos, best_score / len(needle)

def color(seq, pos, adapter, include_brackets=True):
    bits = [seq[:pos]]
    if include_brackets:
        bits.append(COLOR_BLUE_BOLD)
        bits.append("[")
    for s, a in zip (seq[pos:], adapter):
        bits.append(COLOR_GREEN_BOLD if s == a else COLOR_RED_BOLD)
        bits.append(s)
    if include_brackets:
        bits.append(COLOR_BLUE_BOLD)
        bits.append("]")
    bits.append(COLOR_END)
    bits.append(seq[pos + len(adapter):])
    return "".join(bits)
                
mismatches = 0
read_pairs = 0
with open(fastq1) as inf1:
    with open(fastq2) as inf2:
        for (t1, s1, q1), (t2, s2, q2) in zip(
                FastqGeneralIterator(inf1),
                FastqGeneralIterator(inf2)):
            read_pairs += 1
            
            assert t1.split()[0] == t2.split()[0]

            pos1, score1 = match(s1, fwd)
            pos2, score2 = match(s2, rev)

            if score1 < .9 or score2 < .9:
                continue

            if pos1 != pos2:
                if True:
                    print(t1)
                    print(color(s1, pos1, fwd))
                    print(color(s2, pos2, rev))
                    print()
                mismatches += 1
                
print(mismatches, read_pairs)
