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

mode, = sys.argv[1:]

spec = importlib.util.spec_from_file_location(
    'tuning_scoring',
    os.path.join(os.path.dirname(__file__), "2023-12-07--tuning-scoring.py"))
tuning_scoring = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tuning_scoring)

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

plot_data_fname = "bowtie-scores-categorized"
plot_data = []
if os.path.exists(plot_data_fname):
    with open(plot_data_fname) as inf:
        plot_data = json.load(inf)

def load_alignments(fname):
    als = defaultdict(list)
    with gzip.open(fname, "rt") as inf:
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
    return als

def process_sample(paper, bioproject, sample):
    dashboard_dir = bioproject_to_dashboard_dir[bioproject]

    hvr = {}
    hvr_fname = os.path.join(
        dashboard_dir, "hvreads", "%s.hvreads.json" % sample)
    if not os.path.exists(hvr_fname):
        return

    with open(hvr_fname) as inf:
        hvr = json.load(inf)

    alignments_dir = os.path.join(dashboard_dir, "alignments")
    hv_alignments_fname = os.path.join(
        alignments_dir, "%s.hv.alignments.tsv.gz" % sample)
    human_alignments_fname = os.path.join(
        alignments_dir, "%s.human.alignments.tsv.gz" % sample)

    if not os.path.exists(hv_alignments_fname) or \
       not os.path.exists(human_alignments_fname):
        return

    hv_als = load_alignments(hv_alignments_fname)
    human_als = load_alignments(human_alignments_fname)

    #print("n_reads", metadata_samples[sample]["reads"])
    #print("n_hv_als", len(hv_als))
    #print("n_human_als", len(human_als))
    #print("n_hvreads", len(hvr))

    scores_to_print = []
    data = []
    seen = set()
    for read_id, hv_al in sorted(hv_als.items(),reverse=False):
        if read_id in human_als:
            continue

        score = sum(a.score for a in hv_al)
        length = sum(a.length for a in hv_al)

        if mode == "plot":
            if read_id in tuning_scoring.no_hv:
                is_hv = False
            elif read_id in tuning_scoring.yes_hv:
                is_hv = True
            else:
                continue
            seen.add(read_id)
            data.append([read_id, score, length, is_hv])
        elif mode == "sample":
            bend1_length = 40
            bend1_score = 60

            bend2_length = 80
            bend2_score = 95

            max_length = 200
            max_score = 100

            min_length_cutoff = 28
            uncollapsed_score = 269
            
            if length < min_length_cutoff:
                continue
            elif not read_id.startswith("M_"):
                border_score_for_length = uncollapsed_score
            elif length < bend1_length:
                border_score_for_length = bend1_score
            elif length < bend2_length:
                border_score_for_length = (
                    (length - bend1_length) / (bend2_length - bend1_length)
                    * (bend2_score - bend1_score)) + bend1_score
            else:
                border_score_for_length = (
                    (length - bend2_length) / (max_length - bend2_length)
                    * (max_score - bend2_score)) + bend2_score

            if read_id.startswith("M_"):
                #if abs(score - border_score_for_length) > 5:
                #    continue
                if score < border_score_for_length or \
                   score > border_score_for_length + 10:
                    continue
            else:
                #if abs(score - border_score_for_length) > 35:
                #    continue
                if score < border_score_for_length or \
                   score > border_score_for_length + 10:
                    continue
                
            should_skip = False
            for already_read_id, _, _, _ in plot_data:
                if read_id == already_read_id:
                    should_skip = True
                    break
            
            if read_id.startswith("M_") and False:
                for _, other_score, other_length, _ in plot_data:
                    if abs(other_length - length) < 2:
                        should_skip = True
                        break

            if should_skip:
                continue

            #if score < border_score_for_length:
            #    continue
            #if score > border_score_for_length + 10:
            #    continue

            #print("%s, %s; border is %.0f" % (
            #    length, score, border_score_for_length))
            
            scores_to_print.append((length, score, read_id))
            
    if mode == "sample":
        print("\n%s matching read pairs, %s collapsed and %s uncollapsed\n" % (
            len(scores_to_print),
            len([x for x in scores_to_print if x[-1].startswith("M_")]),
            len([x for x in scores_to_print if not x[-1].startswith("M_")])))
        input("press [enter] to continue, Ctrl+C to quit")
        
        random.shuffle(scores_to_print)
        targets = {}
        for length, score, read_id in scores_to_print[:20]:
            targets[read_id] = length, score
        
        reads = defaultdict(list)
        for suffix in [".pair1.truncated", ".pair2.truncated", ".collapsed"]:
            with gzip.open("%s%s.gz" % (sample, suffix), "rt") as inf:
                for (title, sequence, quality) in FastqGeneralIterator(inf):
                    read_id = title.split()[0]
                    if read_id in targets:
                        reads[read_id].append(sequence)
                        
        for read_id, (length, score) in sorted(targets.items()):
            print("%s (%s, %s)" % (read_id, length, score))
            for read in reads[read_id]:
                print(read)
            
    return data

scores = []
for paper in sorted(metadata_papers):
    if "Rothman" not in paper and "Brinch" not in paper:
        continue

    for bioproject in metadata_papers[paper]["projects"]:
        for sample in metadata_bioprojects[bioproject]:
            if mode == "plot":
                if sample not in [
                        "SRR14530882",
                        "SRR14530891",
                        "SRR14530884",
                        "SRR14530771",
                        "ERR3563070",
                ]:
                    continue
            elif mode == "sample":
                if sample != "ERR3563070":
                    continue

            if metadata_samples[sample].get("enrichment", None) == "panel":
                continue

            scores.extend(process_sample(paper, bioproject, sample))

if mode == "plot":
    with open(plot_data_fname, "w") as outf:
        json.dump(scores, outf, indent=2)
