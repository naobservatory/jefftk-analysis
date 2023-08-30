#!/usr/bin/env python3

import os
import glob
import json
import random
from collections import defaultdict


MGS_PIPELINE_DIR="/Users/jeffkaufman/code/mgs-pipeline"
MGS_RESTRICTED_DIR="/Users/jeffkaufman/code/mgs-restricted"

bioproject_to_s3_bucket = {}

with open(os.path.join(MGS_PIPELINE_DIR, "dashboard",
                       "metadata_papers.json")) as inf:
    metadata_papers = json.load(inf)

with open(os.path.join(MGS_PIPELINE_DIR, "dashboard",
                       "metadata_bioprojects.json")) as inf:
    metadata_bioprojects = json.load(inf)
    for bioproject in metadata_bioprojects:
        bioproject_to_s3_bucket[bioproject] = "nao-mgs"

with open(os.path.join(MGS_PIPELINE_DIR, "dashboard",
                       "metadata_samples.json")) as inf:
    metadata_samples = json.load(inf)

with open(os.path.join(MGS_RESTRICTED_DIR, "dashboard",
                       "metadata_papers.json")) as inf:
    metadata_papers.update(json.load(inf))

with open(os.path.join(MGS_RESTRICTED_DIR, "dashboard",
                       "metadata_bioprojects.json")) as inf:
    restricted_metadata_bioprojects = json.load(inf)
    metadata_bioprojects.update(restricted_metadata_bioprojects)
    for bioproject in restricted_metadata_bioprojects:
        bioproject_to_s3_bucket[bioproject] = "nao-restricted"

with open(os.path.join(MGS_RESTRICTED_DIR, "dashboard",
                       "metadata_samples.json")) as inf:
    metadata_samples.update(json.load(inf))

sample_to_paper = {}
sample_to_bioproject = {}
for paper in metadata_papers:
    for bioproject in metadata_papers[paper]["projects"]:
        for sample in metadata_bioprojects[bioproject]:
            enrichment_suffix = ""
            if metadata_samples[sample].get("enrichment", "") == "panel":
                enrichment_suffix = " panel"

            sample_to_bioproject[sample] = bioproject
            sample_to_paper[sample] = "%s%s" % (
                paper, enrichment_suffix)

combined_by_paper = defaultdict(list)
for fname in glob.glob("combined-vb/*.json"):
    with open(fname) as inf:
        sample, *_ = os.path.basename(fname).split(".")
        paper = sample_to_paper.get(sample, None)
        if not paper:
            continue
        
        for seq_id, kraken_info in json.load(inf)["both"]:
            combined_by_paper[paper].append((
                sample_to_bioproject[sample], sample, seq_id, kraken_info))

with open("targets.tsv", "w") as outf:
    for paper, combined in sorted(combined_by_paper.items()):
        random.shuffle(combined)
        for bioproject, sample, seq_id, kraken_info in combined[:10]:
            outf.write("%s\t%s\t%s\t%s\t%s\n" % (
                bioproject_to_s3_bucket[bioproject],
                bioproject,
                sample,
                seq_id,
                kraken_info))
