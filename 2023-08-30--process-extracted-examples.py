#!/usr/bin/env python3

import os
import glob
import json
import random
from collections import defaultdict, Counter


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

examples_by_paper = defaultdict(list)
for fname in glob.glob("combined-vb-examples/*.txt"):
    with open(fname) as inf:
        s3_bucket, bioproject, sample, seq_id, kraken_info = next(
            inf).strip().split("\t")
        sequences = []
        provided_seq_id = None
        for line in inf:
            line = line.strip()
            if line == "--":
                continue
            elif line.startswith("@"):
                assert not provided_seq_id
                provided_seq_id = line.removeprefix("@").split()[0]
            else:
                assert provided_seq_id
                if seq_id == provided_seq_id:
                    sequences.append(line)
                provided_seq_id = None
        assert len(sequences) in [1, 2]
                    
        examples_by_paper[sample_to_paper[sample]].append((
            sample, seq_id, kraken_info, sequences))

ALL_COLORS = [
    '\x1b[1;31m',
    '\x1b[1;32m',
    '\x1b[1;33m',
    '\x1b[1;34m',
    '\x1b[1;35m',
    '\x1b[1;36m',
]
COLOR_END = '\x1b[0m'
K = 35
def colorize(kraken_info, sequences):
    hits = Counter()
    for token in kraken_info.split():
        hit, count = token.split(":")
        if hit in ["0", "1", "|", "A"]:
            continue
        hits[hit] += int(count)

    colors = {} # hit -> color
    available_colors = ALL_COLORS[:]
    for _, hit in sorted([
            (c, h) for (h, c) in hits.items()], reverse=True):
        if available_colors:
            colors[hit] = available_colors.pop(0)

    colorized_kraken_info = []
    for token in kraken_info.split():
        hit, count = token.split(":")
        color = colors.get(hit, None)
        if color:
            token = color + token + COLOR_END
        colorized_kraken_info.append(token)
    colorized_kraken_info = " ".join(colorized_kraken_info)

    color_indexes = {} # read_number, pos -> color

    pos = 0
    read_number = 0
    for segment in kraken_info.split():
        hit, count = segment.split(":")
        if hit == "|":
            read_number += 1
            pos = 0
            continue

        count = int(count)
        if hit in colors:
            for i in range(K):
                color_indexes[read_number, pos + i] = colors[hit]

        pos += count

    colorized_sequences = []
    for read_number, sequence in enumerate(sequences):
        colorized = []
        for pos, val in enumerate(sequence):
            color = color_indexes.get((read_number, pos), None)
            if read_number == 0:
                if not color:
                    colorized.append(val)
                else:
                    colorized.append(color)
                    colorized.append(val)
                    colorized.append(COLOR_END)
            else:
                c_val = {'A': 'T',
                         'T': 'A',
                         'C': 'G',
                         'G': 'C',
                         'N': 'N'}[val]
                if not color:
                    colorized.insert(0, c_val)
                else:
                    colorized.insert(0, COLOR_END)
                    colorized.insert(0, c_val)
                    colorized.insert(0, color)

        colorized_sequences.append("".join(colorized))
    
    return colorized_kraken_info, colorized_sequences
        
for paper in sorted(examples_by_paper):
    if "Johnson" not in paper or "panel" in paper: continue
    print(paper)
    for sample, seq_id, kraken_info, sequences in examples_by_paper[paper]:
        colorized_kraken_info, colorized_sequences = colorize(
            kraken_info, sequences)
        print("  %s" % seq_id)
        print("  %s" % colorized_kraken_info)
        for sequence in colorized_sequences:
            print("  %s" % (sequence))
