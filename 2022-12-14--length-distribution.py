#!/usr/bin/env python3

import glob
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from collections import Counter

all_lengths = {}
for fname in glob.glob("*_R1.fastq.gz"):
    basename = fname[:-len('_R1.fastq.gz')]
    lengths = Counter()
    
    for fastq_fname in [
            basename + '.collapsed',
            basename + '.collapsed.truncated',
            basename + '.paired.truncated',
            basename + '.singleton.truncated']:
        with open(fastq_fname) as inf:
            for (title, sequence, quality) in FastqGeneralIterator(inf):
                length = len(sequence)
                if fastq_fname.endswith(".paired.truncated"):
                    length = -1
                elif fastq_fname.endswith(".singleton.truncated"):
                    length = -2
                lengths[length] += 1

    all_lengths[basename] = lengths

max_len = 589
min_len = -2

rows = [["length"] + list(sorted(all_lengths))]
for length in range(min_len, max_len + 1):
    rows.append([length] +
                [all_lengths[basename][length]
                 for basename in sorted(all_lengths)])

with open('insert-lengths.tsv', 'w') as outf:
    for row in rows:
        outf.write("\t".join(str(x) for x in row) + "\n")
