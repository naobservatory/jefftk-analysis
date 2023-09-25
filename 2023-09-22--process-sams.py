#!/usr/bin/env python3

import os
import glob
import json
import pysam
import Bio.SeqIO
from collections import Counter
import matplotlib.pyplot as plt

metadata = {}
for fname in ["papers", "bioprojects", "samples"]:
    with open("/Users/jeffkaufman/code/mgs-pipeline/dashboard/metadata_%s.json" %
              fname) as inf:
        metadata[fname] = json.load(inf)

genomes = {}
for genome in glob.glob("raw-genomes/*.fna"):
    with open(genome) as inf:
        for record in Bio.SeqIO.parse(inf, "fasta"):
            genomes[record.id] = record.seq
            
def process_sam(fname):
    sample = os.path.basename(fname)
    with pysam.AlignmentFile(fname, "r")  as sam:
        for record in sam:
            ref = sam.get_reference_name(record.reference_id)
            ref_seq = genomes[ref]
            qry_seq = record.query_sequence
            ref_pos = record.reference_start

            longest_alignment = max(
                [length for category, length in
                 record.cigartuples
                 if category == 0], default=0)
            longest_soft_clip = max(
                [length for category, length in
                 record.cigartuples
                 if category == 4], default=0)
            all_soft_clipped = sum(
                [length for category, length in
                 record.cigartuples
                 if category == 4])

            score = record.get_tag("AS") / (len(qry_seq) - all_soft_clipped)

            if longest_soft_clip < 30 or longest_alignment < 30: continue
            
            print("%.2f" % (score), longest_soft_clip, longest_alignment, record.query_name)

for paper in metadata["papers"]:
    if "Rothman" not in paper: continue
    for bioproject in metadata["papers"][paper]["projects"]:
        for sample in metadata["bioprojects"][bioproject]:
            if metadata["samples"][sample].get(
                    "enrichment", "") == "panel":
                continue
            
            fname = "hvsams/%s.sam" % sample
            if os.path.exists(fname):
                process_sam(fname)
        
