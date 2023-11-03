#!/usr/bin/env python3

# We want to know the fragment size distribution of each sample, broken down
# by kraken categorization.  Full species assignment is too much data, but we
# could at least do overall, bacteria, viruses, and human viruses.
#
# Fragment length data is currently only in the lengths of sequences in the
# cleaned/ files, while read-level species assignment is only in the processed/
# files.
#
# I think what we want is, for each sample:
#
# 1. Read through processed/ and take the first 1M read IDs for each of the
#    four categories.
# 2. Read through cleaned/ and get the lengths.
# 3. Make nice plots.
#
# I'll prototype this here, and then move it over to mgs-pipeline/run.py

import glob
import gzip
import subprocess
from collections import Counter

# We take the first N for each file, and then cut down to N max by weighting
# each file.
TARGET_LEN=10 # 1_000

human_viruses = set()
with open("/Users/jeffkaufman/code/mgs-pipeline/human-viruses.tsv") as inf:
    for line in inf:
        taxid, _ = line.split("\t")
        human_viruses.add(int(taxid))

parents = {}
with open("/Users/jeffkaufman/code/mgs-pipeline/dashboard/nodes.dmp") as inf:
    for line in inf:
        child_taxid, parent_taxid, *_ = \
            line.replace("\t|\n", "").split("\t|\t")
        parents[int(child_taxid)] = int(parent_taxid)
        
def taxid_under(clade, taxid):
    while taxid not in [0, 1]:
        if taxid == clade:
            return True
        taxid = parents[taxid]
    return False
        
def taxid_matches(taxid, category):
    if category == "all":
        return True
    if category == "humanviral":
        return taxid in human_viruses
    return taxid_under(
        {
            "bacterial": 2,
            "viral": 10239,
        }[category], taxid)

read_ids = {}
full_counts = Counter()
for fname in [
        "s3://nao-mgs/PRJNA729801/processed/SRR14530724.collapsed.truncated.kraken2.tsv.gz",
        "s3://nao-mgs/PRJNA729801/processed/SRR18341134.collapsed.truncated.kraken2.tsv.gz",
]:
    read_ids[fname] = {
        "all": [],
        "bacterial": [],
        "viral": [],
        "humanviral": [],
    }
    process = subprocess.Popen(["aws", "s3", "cp", fname, "-"],
                               stdout=subprocess.PIPE,
                               shell=False)
    with gzip.open(process.stdout, "rt") as inf:
        for line in inf:
            bits = line.removesuffix("\n").split("\t")
            read_id = bits[1]
            full_assignment = bits[2]

            taxid = int(full_assignment.split()[-1].rstrip(")"))

            for category in read_ids[fname]:
                if len(read_ids[fname][category]) < TARGET_LEN:
                    if taxid_matches(taxid, category):
                        read_ids[fname][category].append(read_id)
    for category in read_ids[fname]:
        full_counts[category] += len(read_ids[fname][category])

subsetted_ids = {}
for category, full_count in full_counts.items():
    subsetted_ids[category] = []
    for fname in read_ids:
        if full_count <= TARGET_LEN:
            subsetted_ids[category].extend(read_ids[fname][category])
        else:
            target = TARGET_LEN * len(read_ids[fname][category]) // full_count
            subsetted_ids[category].extend(read_ids[fname][category][:target])

import pprint
pprint.pprint(subsetted_ids)
        

    
