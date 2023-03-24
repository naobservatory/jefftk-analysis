#!/usr/bin/env python3

import os
import glob
from collections import Counter
from collections import defaultdict

# study -> accession -> n_reads
n_reads = defaultdict(dict)

for n_reads_fname in glob.glob(
        "/Users/jeffkaufman/code/mgs-pipeline/studies/*/metadata/*.n_reads"):
    study = n_reads_fname.split("/")[-3]
    accession = n_reads_fname.split("/")[-1].replace(".n_reads", "")
    with open(n_reads_fname) as inf:
        n_reads[study][accession] = int(inf.read())

studies = list(sorted(n_reads))

# taxid, name -> study -> relab
counts = defaultdict(Counter)

# taxid, name -> total
total_counts = Counter()

for study in studies:
    hv_totals = Counter()
    study_reads = 0

    for accession in n_reads[study]:
        hv = "%s.humanviruses.tsv" % accession

        if not os.path.exists(hv):
            continue
        
        with open(hv) as inf:
            for line in inf:
                taxid, count, name = line.strip().split("\t")
                hv_totals[taxid, name] += int(count)

        study_reads += n_reads[study][accession]

    for (taxid, name), count in hv_totals.items():
        counts[taxid, name][study] = count / study_reads
        total_counts[taxid, name] += count

print("name", "total", *studies, sep="\t")
        
for total_count, (taxid, name) in sorted([
        (count, key) for (key, count) in total_counts.items()], reverse=True):

    print("%s (%s)" % (name, taxid), end="\t")
    print("%s" % total_count, end="\t")
    
    cols = ["%.1e" % counts[taxid, name][study]
            for study in studies]
    print(*cols, sep="\t")
