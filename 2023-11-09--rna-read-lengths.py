#!/usr/bin/env python3

import os
import gzip
import json
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from collections import Counter, defaultdict

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

# https://stackoverflow.com/questions/20601872/numpy-or-scipy-to-calculate-weighted-median
def weighted_quantiles_interpolate(values, weights, quantiles=0.5):
    i = np.argsort(values)
    weights = np.array(weights)
    c = np.cumsum(weights[i])
    q = np.searchsorted(c, quantiles * c[-1])
    return np.where(c[q]/c[-1] == quantiles, 0.5 *
                    (values[i[q]] + values[i[q+1]]), values[i[q]])
    
data = defaultdict(dict)
def populate_data(paper, samples, chosen_category):
    all_category_lengths = Counter()
    for sample in samples:
        with gzip.open("readlengths/%s.rl.json.gz" % sample, "rt") as inf:
            for category, category_lengths in json.load(inf).items():
                if category != chosen_category: continue

                #if max(int(x) for x in category_lengths if x != "NC") < 210:
                #    continue # skip 2x100 samples
                
                for length, count in category_lengths.items():
                    if length != "NC":
                        length = int(length)
                    all_category_lengths[length] += count
    total_counts = sum(all_category_lengths.values())
    if total_counts < 100:
        return
    
    nc_fraction = all_category_lengths["NC"]/total_counts

    if nc_fraction == 1:
        return

    max_len = max(length
                  for length in all_category_lengths
                  if length != "NC")
    xs = list(range(max_len + 1))
    ys = [100 * all_category_lengths[x] / total_counts for x in xs]

    weighted_median = weighted_quantiles_interpolate(xs, ys)
    print(paper, chosen_category, weighted_median)

    label = "nc=%.0f%%" % (100*nc_fraction)
    data[paper][chosen_category] = label, xs, ys

paper_samples = defaultdict(list)

for paper in metadata_papers:
    if metadata_papers[paper]["na_type"] in [
            "DNA", "DNA+RNA", "RNA+DNA"]:
        continue
    
    paper_subset_samples = defaultdict(list)
    missing = False
    for bioproject in metadata_papers[paper]["projects"]:
        for sample in metadata_bioprojects[bioproject]:
            if not os.path.exists("readlengths/%s.rl.json.gz" % sample):
                missing = True

            if metadata_samples[sample].get("enrichment", None) == "panel":
                continue

            if metadata_samples[sample].get(
                    "na_type", metadata_papers[paper]["na_type"]) != "RNA":
                continue

            if metadata_samples[sample].get(
                    "collection", "wastewater") != "wastewater":
                continue

            if "airport" in metadata_samples[sample]:
                continue

            label = paper
            if bioproject_to_s3_bucket[bioproject] == "nao-restricted":
                continue

            if "Rothman" in paper:
                with gzip.open(
                        "readlengths/%s.rl.json.gz" % sample, "rt") as inf:
                    rls = json.load(inf)
                    max_len = max(
                        max([int(x) for x in lengths if x != "NC"],
                             default=0)
                         for lengths in rls.values())
                    if max_len > 210:
                        label += " 2x150"
                    else:
                        label += " 2x100"
            
            paper_subset_samples[label].append(sample)
    if missing:
        # Respond to missing data by hiding the whole paper's chart
        continue

    for paper_subset in paper_subset_samples:
        paper_samples[paper_subset].extend(paper_subset_samples[paper_subset])

for paper, samples in paper_samples.items():
    for category in "av":
        populate_data(paper, samples, chosen_category=category)

ncols = 3
nrows = math.ceil(len(data) / ncols)
fig, axs = plt.subplots(constrained_layout=True,
                        sharex=True,
                        figsize=(4*ncols, 4*nrows),
                        nrows=nrows,
                        ncols=ncols)
fig.supxlabel("read length")
fig.supylabel("percentage reads")
fig.suptitle("RNA Read Lengths by Paper, All Reads vs Viral Reads")

for i, paper in enumerate(sorted(data)):
    ax = axs[i // ncols][i % ncols]
    ax.yaxis.set_major_formatter(mtick.PercentFormatter())
    title = [paper]
    for category, (label, xs, ys) in sorted(data[paper].items()):
        ax.plot(xs, ys, label=category)
        title.append("%s %s" % (category, label))
    ax.set_title("\n".join(title))
    ax.legend()
fig.savefig("rna-read-lengths.png", dpi=180)
plt.clf()
