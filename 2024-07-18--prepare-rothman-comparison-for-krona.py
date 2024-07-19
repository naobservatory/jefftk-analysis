#!/usr/bin/env python3

# aws s3 sync s3://nao-mgs/PRJNA729801/cladecounts-2024-06/ cladecounts/
import os
import sys
import gzip
import json
from collections import Counter, defaultdict

with open(os.path.expanduser(
        "~/code/mgs-pipeline/dashboard/metadata_samples.json")) as inf:
    metadata_samples = json.load(inf)

with open(os.path.expanduser(
        "~/code/mgs-pipeline/dashboard/metadata_bioprojects.json")) as inf:
    metadata_bioprojects = json.load(inf)
    
direct_assignments_by_category = defaultdict(Counter)
for fname in os.listdir("cladecounts"):
    sample = fname.removesuffix(".tsv.gz").removesuffix("-div0000")
    if sample not in metadata_bioprojects["PRJNA729801"]:
        continue
    
    if metadata_samples[sample].get("enrichment") == "panel":
        continue

    category = metadata_samples[sample]["fine_location"]

    with gzip.open(os.path.join("cladecounts", fname), "rt") as inf:
        for line in inf:
            (taxid,
             direct_assignments,
             direct_hits,
             clade_assignments,
             clade_hits) = [int (x) for x in line.rstrip("\n").split("\t")]

            direct_assignments_by_category[category][taxid] += direct_assignments

for category, counts in direct_assignments_by_category.items():
    with open("%s.tsv" % category, "w") as outf:
        for taxid, count in sorted(counts.items()):
            outf.write("%s\t%s\n" % (taxid, count))
