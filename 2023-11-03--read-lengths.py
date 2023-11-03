#!/usr/bin/env python3

import os
import json
import gzip
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

    
def graph_category_lengths(sample_lists, fname, chosen_category):
    fig, ax = plt.subplots(constrained_layout=True)
    ax.yaxis.set_major_formatter(mtick.PercentFormatter())
    fig.supxlabel("read length")
    fig.supylabel("percentage reads")
    for label, sample_list in sorted(sample_lists.items()):
        all_category_lengths = Counter()
        for sample in sample_list:
            with gzip.open("readlengths/%s.rl.json.gz" % sample, "rt") as inf:
                for category, category_lengths in json.load(inf).items():
                    if category != chosen_category: continue
                    for length, count in category_lengths.items():
                        if length != "NC":
                            length = int(length)
                        all_category_lengths[length] += count
        total_counts = sum(all_category_lengths.values())
        nc_fraction = all_category_lengths["NC"]/total_counts
        
        if nc_fraction == 1:
            continue

        min_len = min(length
                      for length in all_category_lengths
                      if length != "NC")
        max_len = max(length
                      for length in all_category_lengths
                      if length != "NC")
        
        xs = list(range(min_len, max_len + 1))
        ys = [100 * all_category_lengths[x] / total_counts for x in xs]

        label = label + " (nc=%.0f%%)" % (100*nc_fraction)
        
        ax.plot(xs, ys, label=label)
    ax.legend()
    fig.savefig(fname, dpi=180)
    plt.clf()

for paper in metadata_papers:
    try:
        sample_lists = defaultdict(list)
        for bioproject in metadata_papers[paper]["projects"]:
            for sample in metadata_bioprojects[bioproject]:
                if bioproject_to_s3_bucket[bioproject] == "nao-mgs":
                    label = bioproject
                else:
                    label = sample
                    if bioproject.startswith("NAO"):
                        label = label.removesuffix("-1").removesuffix("-2").removesuffix("_L1").removesuffix("_L2")
                sample_lists[label].append(sample)        
                
        graph_category_lengths(sample_lists,
                               fname=paper.replace(" ", "") + "-samples.png",
                               chosen_category="a")
    except FileNotFoundError:
        continue
