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

def rc(s):
    return "".join({'T':'A',
                    'G':'C',
                    'A':'T',
                    'C':'G',
                    'N':'N'}[x] for x in reversed(s))

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

mismatch_no_align_delta = []
mismatch_align_delta = []

max_n = 25
cur_n = 0
with open(fastq1) as inf1:
    with open(fastq2) as inf2:
        for (t1, s1, q1), (t2, s2, q2) in zip(
                FastqGeneralIterator(inf1),
                FastqGeneralIterator(inf2)):
            
            assert t1.split()[0] == t2.split()[0]

            pos1, score1 = match(s1, fwd)
            pos2, score2 = match(s2, rev)

            if score1 < .9 or score2 < .9:
                continue

            if pos1 != pos2 and pos1 > 0 and pos2 > 0:
                trimmed_s1 = s1[:pos1]
                trimmed_s2 = s2[:pos2]

                if len(trimmed_s2) > len(trimmed_s1):
                    trimmed_s1, trimmed_s2 = trimmed_s2, trimmed_s1
                
                rc_trimmed_s2 = rc(trimmed_s2)
                rc_alignment_pos, rc_alignment_score = match(
                    trimmed_s1, rc_trimmed_s2)

                if rc_alignment_score > 0.9:
                    if False:
                        print(t1)
                        print(color(s1, pos1, fwd))
                        print(color(s2, pos2, rev))
                        print(color(trimmed_s1, rc_alignment_pos, rc_trimmed_s2))
                        print()
                    mismatch_align_delta.append(abs(pos1 - pos2))
                else:
                    mismatch_no_align_delta.append(abs(pos1 - pos2))

                    cur_n += 1
                    if cur_n <= max_n:
                        print("@%s\n%s\n+\n%s\n@%s\n%s\n+\n%s" % (
                            t1, s1, q1,
                            t2, s2, q2))                        

#print(sorted(mismatch_no_align_delta))
#print(sorted(mismatch_align_delta))
