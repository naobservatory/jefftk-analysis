#!/usr/bin/env python3

import os
import gzip
import json
import math
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

data = defaultdict(dict)
def populate_data(paper, samples, chosen_category):
    all_category_lengths = Counter()
    for sample in samples:
        with gzip.open("readlengths/%s.rl.json.gz" % sample, "rt") as inf:
            for category, category_lengths in json.load(inf).items():
                if category != chosen_category: continue
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
    max_len += 1
    all_category_lengths[max_len] = all_category_lengths["NC"]
    del all_category_lengths["NC"]

    xs = list(range(max_len + 1))
    ys = [all_category_lengths[x] / total_counts for x in xs]

    label = "nc=%.0f%%" % (100*nc_fraction)
    data[paper][chosen_category] = label, xs, ys

paper_samples = defaultdict(list)

for paper in metadata_papers:
    if not (paper.startswith("NAO-23") or paper.startswith("MJ")):
        continue
    
    paper_subset_samples = defaultdict(list)
    missing = False
    for bioproject in metadata_papers[paper]["projects"]:
        for sample in metadata_bioprojects[bioproject]:
            if not os.path.exists("readlengths/%s.rl.json.gz" % sample):
                missing = True

            paper_subset = sample
            if paper.startswith("NAO-23"):
                for suffix in ["-1", "-2", "_L1", "_L2"]:
                    paper_subset = paper_subset.removesuffix(suffix)            
                
            paper_subset_samples[paper_subset].append(sample)
    if missing:
        # Respond to missing data by hiding the whole paper's chart
        continue

    for paper_subset in paper_subset_samples:
        paper_samples[paper_subset].extend(paper_subset_samples[paper_subset])

for paper, samples in paper_samples.items():
    for category in "av":
        populate_data(paper, samples, chosen_category=category)

ncols = 4
nrows = math.ceil(len(data) / ncols)
fig, axs = plt.subplots(constrained_layout=True,
                        sharex=True,
                        figsize=(3*ncols, 3*nrows),
                        nrows=nrows,
                        ncols=ncols)
fig.supxlabel("read length")
fig.supylabel("percentage reads")
fig.suptitle("Read lengths by paper, All Reads vs Viral Reads")

for i, paper in enumerate(sorted(data)):
    ax = axs[i // ncols][i % ncols]
    ax.yaxis.set_major_formatter(mtick.PercentFormatter())
    title = [paper]
    for category, (label, xs, ys) in sorted(data[paper].items()):
        ax.plot(xs, ys, label=category)
        title.append("%s %s" % (category, label))
    ax.set_title("\n".join(title))
    ax.legend()
fig.savefig("read-lengths-by-sample-av.png", dpi=180)
plt.clf()
