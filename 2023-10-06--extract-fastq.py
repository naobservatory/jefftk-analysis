#!/usr/bin/env python3
import sys
import glob
from Bio.SeqIO.QualityIO import FastqGeneralIterator

read_ids = sys.argv[1:]

for fname in glob.glob("hvfastqs/*.fastq"):
    with open(fname) as inf:
        for (title, sequence, quality) in FastqGeneralIterator(inf):
            if title in read_ids:
                print("@%s\n%s\n+\n%s" % (
                    title, sequence, quality))
