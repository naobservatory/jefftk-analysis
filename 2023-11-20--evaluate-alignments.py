#!/usr/bin/env python3

import os
import json
import math
from collections import namedtuple, defaultdict

MGS_PIPELINE_DIR="/Users/jeffkaufman/code/mgs-pipeline"
MGS_RESTRICTED_DIR="/Users/jeffkaufman/code/mgs-restricted"

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

def process_sample(paper, bioproject, sample):
    dashboard_dir = bioproject_to_dashboard_dir[bioproject]

    hvr = {}
    hvr_fname = os.path.join(
        dashboard_dir, "hvreads", "%s.hvreads.json" % sample)
    if os.path.exists(hvr_fname):
        with open(hvr_fname) as inf:
            hvr = json.load(inf)

    hvr = {}
    hvr_fname = os.path.join(
        dashboard_dir, "hvreads", "%s.hvreads.json" % sample)
    if os.path.exists(hvr_fname):
        with open(hvr_fname) as inf:
            hvr = json.load(inf)

    als = defaultdict(list)
    als_fname = os.path.join(
        dashboard_dir, "alignments", "%s.alignments.tsv" % sample)
    if os.path.exists(als_fname):
        with open(als_fname) as inf:
            for line in inf:
                (read_id, genome, taxid, cigar, ref_start,
                 score, length) = line.rstrip("\n").split("\t")
                als[read_id].append(Alignment(
                    read_id,
                    genome,
                    int(taxid),
                    cigar,
                    int(ref_start),
                    int(score),
                    int(length)))

    for al in als:
        assert al in hvr

    n_assigned_hv = 0
    n_aligned_hv = 0
    n_assigned_but_not_aligned = 0
    n_aligned_but_not_assigned = 0
        
    for read_id in sorted(hvr):
        assigned_taxid = hvr[read_id][0]

        assigned_hv = assigned_taxid in human_viruses
        aligned_hv = False
        if read_id in als:
            combined_score = sum(al.score for al in als[read_id])
            combined_length = sum(al.length for al in als[read_id])
            aligned_hv = combined_score / math.log(combined_length) > 22

        n_assigned_hv += int(assigned_hv)
        n_aligned_hv += int(aligned_hv)
        n_assigned_but_not_aligned += int(assigned_hv and not aligned_hv)
        n_aligned_but_not_assigned += int(aligned_hv and not assigned_hv)

    if False:
        print("%s\t%.0f%%\t%.0f%%\t%.0f%%\t%.0f%%" % (
            sample,
            100 * n_assigned_hv / len(hvr),
            100 * n_aligned_hv / len(hvr),
            100 * n_assigned_but_not_aligned / len(hvr),
            100 * n_aligned_but_not_assigned / len(hvr)))
        
    return (
        n_assigned_hv,
        n_aligned_hv,
        n_assigned_but_not_aligned,
        n_aligned_but_not_assigned,
        len(hvr))
        
print("paper",
      "assigned",
      "aligned",
      "as_not_al",
      "al_not_as",
      sep="\t")

for paper in sorted(metadata_papers):
    n_assigned_hv = n_aligned_hv = n_assigned_but_not_aligned = \
        n_aligned_but_not_assigned = n_hvr = 0
    
    for bioproject in metadata_papers[paper]["projects"]:
        for sample in metadata_bioprojects[bioproject]:
            if metadata_samples[sample].get("enrichment", None) == "panel":
                continue
            # skip non-influent samples
            if paper == "Bengtsson-Palme 2016":
                if not metadata_samples[sample]["fine_location"].startswith(
                        "Inlet"):
                    continue
            if paper == "Ng 2019":
                if metadata_samples[sample]['fine_location'] != "Influent":
                    continue
            if metadata_samples[sample].get(
                    "collection", "wastewater") != "wastewater":
                continue
            
            (s_n_assigned_hv,
             s_n_aligned_hv,
             s_n_assigned_but_not_aligned,
             s_n_aligned_but_not_assigned,
             s_n_hvr) = process_sample(paper, bioproject, sample)
            n_assigned_hv += s_n_assigned_hv
            n_aligned_hv += s_n_aligned_hv
            n_assigned_but_not_aligned += s_n_assigned_but_not_aligned
            n_aligned_but_not_assigned += s_n_aligned_but_not_assigned
            n_hvr += s_n_hvr

    if n_hvr > 0:
        print("%s\t%.0f%%\t%.0f%%\t%.0f%%\t%.0f%%" % (
            paper,
            100 * n_assigned_hv / n_hvr,
            100 * n_aligned_hv / n_hvr,
            100 * n_assigned_but_not_aligned / n_hvr,
            100 * n_aligned_but_not_assigned / n_hvr),
              flush=True)
