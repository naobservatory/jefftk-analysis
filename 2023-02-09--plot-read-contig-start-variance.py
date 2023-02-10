#!/usr/bin/env python3

import sys
import glob
import json
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

samples = {} # accession_{1,2} -> vid -> f/r + loc -> count
for fname in sys.argv[1:]: 
    sample, _ = fname.split("_")
    with open(fname) as inf:
        samples[sample] = json.load(inf)

all_vids = set()
for vids in samples.values():
    all_vids.update(vids)
all_vids = list(sorted(all_vids))

def to_slug(vid):
    return vid.split(".")[0]



max_locs = {
    'f': Counter(), # vid -> max_loc
    'r': Counter(), # vid -> max_loc
}

for sample in samples:
    for vid in all_vids:
        for frloc, count in samples[sample][vid].items():
            fr, loc = frloc[0], int(frloc[1:])
            max_locs[fr][vid] = max(max_locs[fr][vid], loc)

dir_sample_sums = {
    'f': {}, # vid -> [sample1 fsum, sample2 fsum, ...]
    'r': {}, # vid -> [sample1 fsum, sample2 fsum, ...]
}

for vid in all_vids:
    for direction in ['f', 'r']:
        dir_sample_sums[direction][vid] = np.array([
            sum(
                count
                for frloc, count in samples[sample][vid].items()
                if frloc.startswith(direction)) + 1
            for sample in sorted(samples)])    
            
def get_variance(vid, direction):
    xs = []
    ys = []
    for loc in range(max_locs[direction][vid] + 1):
        frloc = "%s%s" % (direction, loc)
        vals = np.array([samples[sample][vid].get(frloc, 0)
                         for sample in sorted(samples)])
        normalized = vals / dir_sample_sums[direction][vid]
        mean = np.mean(normalized)
        errors = normalized - mean
        squared_errors = errors**2
        unbiased_variance = sum(squared_errors) / (len(vals)-1)
        
        xs.append(loc)
        ys.append(unbiased_variance)

    return xs, ys

fig, axs = plt.subplots(nrows=3,
                        ncols=3,
                        #sharex=True,
                        sharey=True,
                        figsize=(12, 12),
                        constrained_layout=True)

plt.suptitle("variance in read starts by location")
fig.supxlabel("position along genome")
fig.supylabel("variance")
for row, vid in enumerate(all_vids):
    slug = to_slug(vid)
    ax = axs[row%3][row//3]

    for direction in ['f', 'r']:
        xs, ys = get_variance(vid, direction)

        ax.set_yscale('log')
        ax.scatter(xs, ys, label=direction, marker=',', s=0.5)
        ax.legend()
        ax.set_title("%s" % slug)

fig.savefig("read-start-variance.png", dpi=180)
plt.clf()

