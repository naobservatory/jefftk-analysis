#!/usr/bin/env python3

import os
import json
import math
from collections import namedtuple

MGS_PIPELINE_DIR="/home/ec2-user/mgs-pipeline"
MGS_RESTRICTED_DIR="/home/ec2-user/mgs-restricted"

bioproject_to_s3_bucket = {}
bioproject_to_dashboard_dir = {}
with open(os.path.join(MGS_PIPELINE_DIR, "dashboard",
                       "metadata_papers.json")) as inf:
    metadata_papers = json.load(inf)

with open(os.path.join(MGS_PIPELINE_DIR, "dashboard",
                       "metadata_bioprojects.json")) as inf:
    metadata_bioprojects = json.load(inf)
    for bioproject in metadata_bioprojects:
        bioproject_to_s3_bucket[bioproject] = "nao-mgs"
        bioproject_to_dashboard_dir[bioproject] = os.path.join(
            MGS_PIPELINE_DIR, "dashboard")

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
        bioproject_to_dashboard_dir[bioproject] = os.path.join(
            MGS_RESTRICTED_DIR, "dashboard")

with open(os.path.join(MGS_RESTRICTED_DIR, "dashboard",
                       "metadata_samples.json")) as inf:
    metadata_samples.update(json.load(inf))

human_viruses = {}
with open(os.path.join(MGS_PIPELINE_DIR, "human-viruses.tsv")) as inf:
    for line in inf:
        taxid, name = line.rstrip("\n").split("\t")
        human_viruses[int(taxid)] = name

Alignment = namedtuple(
    "Alignment",
    ["read_id", "genome", "taxid", "cigar", "ref_start", "score", "length"])

def process_sample(bioproject, sample):
    dashboard_dir = bioproject_to_dashboard_dir[bioproject]
    hvr_fname = os.path.join(
        dashboard_dir, "hvreads", "%s.hvreads.json" % sample)
    if os.path.exists(hvr_fname):
        with open(hvr_fname) as inf:
            for value in json.load(inf).values():
                if type(value[0]) != type(1):
                    return True
                return False
for paper in sorted(metadata_papers):
    for bioproject in metadata_papers[paper]["projects"]:
        for sample in metadata_bioprojects[bioproject]:
            r = process_sample(bioproject, sample)
            if r in [True, False]:
                if r:
                    print(paper, bioproject, sample)
                break
