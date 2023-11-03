#!/usr/bin/env python3

import os
import gzip
import json
import math
from collections import Counter, defaultdict
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

bioproject_to_s3_bucket = {}
MGS_PIPELINE_DIR="/Users/jeffkaufman/code/mgs-pipeline"
MGS_RESTRICTED_DIR="/Users/jeffkaufman/code/mgs-restricted"

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

def populate_data(paper, samples, chosen_category):
    all_category_lengths = Counter()
    for sample in samples:
        with gzip.open("readlengths/%s.rl.json.gz" % sample, "rt") as inf:
            for category, category_lengths in json.load(inf).items():
                if category != chosen_category: continue
                for length, count in category_lengths.items():
                    if length != "NC":
                        length = int(length)
                        if length > 300:
                            # optimize graphs for comparision to 2x150
                            length = "NC"
                    all_category_lengths[length] += count
    total_counts = sum(all_category_lengths.values())
    nc_fraction = all_category_lengths["NC"]/total_counts
        
    if nc_fraction == 1:
        return

    max_len = max(length
                  for length in all_category_lengths
                  if length != "NC")
        
    xs = list(range(max_len + 1))
    ys = [100 * all_category_lengths[x] / total_counts for x in xs]

    label = paper + "\n(nc=%.0f%%)" % (100*nc_fraction)

    data.append((label, xs, ys))

for category in ["a", "v", "h", "b"]:
    # (plot label, xs, ys)
    data = []
    for paper in metadata_papers:
        try:
            samples = []
            for bioproject in metadata_papers[paper]["projects"]:
                for sample in metadata_bioprojects[bioproject]:
                    samples.append(sample)        
                
            populate_data(paper, samples, chosen_category=category)
        except FileNotFoundError:
            continue

    ncols = 4
    nrows = math.ceil(len(data) / ncols)
    fig, axs = plt.subplots(constrained_layout=True,
                            sharex=True,
                            figsize=(2*ncols, 2*nrows),
                            nrows=nrows,
                            ncols=ncols)
    fig.supxlabel("read length")
    fig.supylabel("percentage reads")
    fig.suptitle("Read lengths by paper (%s)" % (
        {
            "a" : "All Reads",
            "v": "Viral Reads",
            "b": "Bacterial Reads",
            "h": "Human Viral Reads"
        }[category]))

    
    for i, (label, xs, ys) in enumerate(sorted(data)):
        ax = axs[i // ncols][i % ncols] 
        ax.yaxis.set_major_formatter(mtick.PercentFormatter())
        ax.plot(xs, ys)
        ax.set_title(label)
    fig.savefig("read-lengths-by-paper-%s.png" % category, dpi=180)
    plt.clf()
