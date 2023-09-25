#!/usr/bin/env python3
import sys
import glob
from Bio.SeqIO.QualityIO import FastqGeneralIterator

sample, read_id, = sys.argv[1:]

for fname in glob.glob("hvfastqs/%s.*.fastq" % sample):
    with open(fname) as inf:
        for (title, sequence, quality) in FastqGeneralIterator(inf):
            if title == read_id:
                print(sequence)
