#!/usr/bin/env python3

import os
import glob
import json
from collections import Counter

MGS_PIPELINE_DIR="/Users/jeffkaufman/code/mgs-pipeline"
MGS_RESTRICTED_DIR="/Users/jeffkaufman/code/mgs-restricted"

with open(os.path.join(MGS_PIPELINE_DIR, "dashboard",
                       "metadata_papers.json")) as inf:
    metadata_papers = json.load(inf)

with open(os.path.join(MGS_PIPELINE_DIR, "dashboard",
                       "metadata_bioprojects.json")) as inf:
    metadata_bioprojects = json.load(inf)

with open(os.path.join(MGS_PIPELINE_DIR, "dashboard",
                       "metadata_samples.json")) as inf:
    metadata_samples = json.load(inf)

with open(os.path.join(MGS_RESTRICTED_DIR, "dashboard",
                       "metadata_papers.json")) as inf:
    metadata_papers.update(json.load(inf))

with open(os.path.join(MGS_RESTRICTED_DIR, "dashboard",
                       "metadata_bioprojects.json")) as inf:
    metadata_bioprojects.update(json.load(inf))

with open(os.path.join(MGS_RESTRICTED_DIR, "dashboard",
                       "metadata_samples.json")) as inf:
    metadata_samples.update(json.load(inf))

sample_to_paper = {}
for paper in metadata_papers:
    for bioproject in metadata_papers[paper]["projects"]:
        for sample in metadata_bioprojects[bioproject]:
            enrichment_suffix = ""
            if metadata_samples[sample].get("enrichment", "") == "panel":
                enrichment_suffix = " panel"

            sample_to_paper[sample] = "%s%s" % (
                paper, enrichment_suffix)

paper_viral = Counter()
paper_both = Counter()
unknown = 0
for fname in glob.glob("*.json"):
    with open(fname) as inf:
        sample, *_ = fname.split(".")
        paper = sample_to_paper.get(sample, None)
        if not paper:
            unknown += 1
            continue
        
        r = json.load(inf)
        paper_both[paper] += r["n_both"]
        paper_viral[paper] += r["n_viral"]

        if paper == "Rothman 2021":
            for seq_id, kraken_info in r["both"]:
                print(kraken_info)

print("%s unknown samples" % unknown)
        
for paper in sorted(paper_viral):
    print(paper, paper_both[paper], paper_viral[paper], sep="\t")

