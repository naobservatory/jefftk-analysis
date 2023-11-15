#!/usr/bin/env python3

import os
import gzip
import json
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from collections import Counter, defaultdict
from scipy.ndimage.filters import gaussian_filter1d

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

colors = {}

for paper in metadata_papers:
    if "Yang" not in paper:
        continue
    
    paper_subset_samples = defaultdict(list)
    missing = False
    location_counts = Counter()
    for bioproject in metadata_papers[paper]["projects"]:
        for sample in sorted(metadata_bioprojects[bioproject]):
            if not os.path.exists("readlengths/%s.rl.json.gz" % sample):
                missing = True

            location = metadata_samples[sample]["location"]
            location_counts[location] += 1
                
            label = "%s-%s" % (
                location, location_counts[location])
                
            paper_subset_samples[label].append(sample)
            colors[label] = {
                "Kashgar": "tab:blue",
                "Hotan": "tab:orange",
                "Mi-Dong": "tab:green",
            }[location]
    if missing:
        # Respond to missing data by hiding the whole paper's chart
        continue

    for paper_subset in paper_subset_samples:
        paper_samples[paper_subset].extend(paper_subset_samples[paper_subset])

for paper, samples in paper_samples.items():
    for category in "v":
        populate_data(paper, samples, chosen_category=category)

fig, ax = plt.subplots(constrained_layout=True,
                       figsize=(10,8))
fig.supxlabel("read length")
fig.supylabel("percentage of all viral reads")
fig.suptitle("Yang 2020 Viral RNA Read Lengths")

ax.yaxis.set_major_formatter(mtick.PercentFormatter())
for i, paper in enumerate(sorted(data)):
    (label, xs, ys), = data[paper].values() # only one category
    ax.plot(xs, ys, label=paper + " " + label, color=colors[paper])
ax.legend()
fig.savefig("rna-read-lengths-yang.png", dpi=180)
plt.clf()

fig, ax = plt.subplots(constrained_layout=True,
                       figsize=(8,5))
fig.supxlabel("read length")
fig.supylabel("percentage of all viral reads")
fig.suptitle("Yang 2020 Viral RNA Read Lengths")

ax.yaxis.set_major_formatter(mtick.PercentFormatter())
for i, paper in enumerate(sorted(data)):
    (label, xs, ys), = data[paper].values() # only one category

    ys_smooth = gaussian_filter1d(ys, sigma=5)
    ax.plot(xs, ys_smooth, label=paper + " " + label, color=colors[paper])
ax.legend()
fig.savefig("rna-read-lengths-yang-smooth.png", dpi=180)
plt.clf()
