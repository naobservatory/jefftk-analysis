#!/usr/bin/env python3

import os
import sys
import gzip
import json
import math
import random
import hashlib
import importlib
from collections import namedtuple, defaultdict, Counter
from Bio.SeqIO.QualityIO import FastqGeneralIterator

to_extract = []
target_read_ids = set()
samples = set()

def read_id_to_sample_id(read_id):
    return read_id.replace("M_", "").split(".")[0]

with open("reads-to-evaluate.tsv") as inf:
    for line in inf:
        read_id, length, score = line.rstrip("\n").split("\t")
        to_extract.append((read_id, length, score))
        target_read_ids.add(read_id)
        samples.add(read_id_to_sample_id(read_id))

FASTQ_DIR="../2023-11-20--evaluate-alignments"

reads = defaultdict(list)
for sample in samples:
    for suffix in [".pair1.truncated", ".pair2.truncated", ".collapsed"]:
        with gzip.open("%s/%s%s.gz" % (FASTQ_DIR, sample, suffix), "rt") as inf:
            for (title, sequence, quality) in FastqGeneralIterator(inf):
                read_id = title.split()[0]
                if read_id in target_read_ids:
                    reads[read_id].append(sequence)

for read_id, length, score in to_extract:
    print("%s (%s, %s)" % (read_id, length, score))
    for read in reads[read_id]:
        print(read)
