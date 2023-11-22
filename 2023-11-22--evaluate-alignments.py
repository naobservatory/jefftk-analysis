#!/usr/bin/env python3

import os
import json
import math
import random
import hashlib
from collections import namedtuple, defaultdict
from Bio.SeqIO.QualityIO import FastqGeneralIterator

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
    als_fname = os.path.join("%s.alignments.masked.tsv" % sample)
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

    n_hit_hv = 0
    n_aligned_hv = 0
    n_hit_but_not_aligned = 0
    n_aligned_but_not_hit = 0

    printed = 0

    # read_id -> info
    to_print = defaultdict(list)
    
    for read_id in sorted(als, key=lambda k: hashlib.sha256(k.encode('utf-8')).digest()):
        #assigned_taxid = hvr[read_id][0]
        #kraken_output = hvr[read_id][1]

        combined_score = sum(al.score for al in als[read_id])
        combined_length = sum(al.length for al in als[read_id])
        normalized_score = combined_score / math.log(combined_length)
        aligned_hv = normalized_score > 22

        if not aligned_hv:
            continue

        aligned_taxid, = set(al.taxid for al in als[read_id])
        
        if read_id not in hvr:
            to_print[read_id].append(read_id)
            to_print[read_id].append("bowtie: %.0f" % normalized_score)
            to_print[read_id].append("%s %s" % (
                aligned_taxid, human_viruses.get(aligned_taxid, "absent")))
            for al in als[read_id]:
                to_print[read_id].append(
                    "  %s %s %s %s" % (
                        al.taxid, al.cigar, al.genome, al.ref_start))
            printed += 1

            if printed == 10:
                break

    for suffix in [".pair1.truncated", ".pair1.truncated", ".collapsed"]:
        with open("%s%s" % (sample, suffix)) as inf:
            for (title, sequence, quality) in FastqGeneralIterator(inf):
                read_id = title.split()[0]
                if read_id in to_print:
                    to_print[read_id].append(sequence)


    for read_id, lines in sorted(to_print.items()):
        for line in lines:
            print(line)
        print()
        
for paper in sorted(metadata_papers):
    if "Rothman" not in paper:
        continue

    for bioproject in metadata_papers[paper]["projects"]:
        for sample in metadata_bioprojects[bioproject]:
            if sample != "SRR14530891":
                continue
            
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

            process_sample(paper, bioproject, sample)
