#!/usr/bin/env python3

import sys
import glob
import json
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

samples = {} # accession_{1,2} -> vid -> f/r + loc -> count
for fname in sys.argv[1:]: 
    accession, _ = fname.split("_")
    with open(fname) as inf:
        samples[accession] = json.load(inf)

all_vids = set()
for vids in samples.values():
    all_vids.update(vids)
all_vids = list(sorted(all_vids))

def to_slug(vid):
    return vid.split(".")[0]



fig, axs = plt.subplots(#nrows=9, ncols=1,
                        #figsize=(4*1, 4*9),
                        nrows=3, ncols=3,
                        figsize=(12, 12),
                        constrained_layout=True)

plt.suptitle("read starts by location")
fig.supxlabel("position along genome")
fig.supylabel("reads starting at this position")
for row, vid in enumerate(all_vids):
    slug = to_slug(vid)
    ax = axs[row%3][row//3]
    #ax = axs[row]

    for direction in ['f', 'r']:
        vals = Counter()
        for accession in samples:
            for frloc, count in samples[accession][vid].items():
                if frloc.startswith(direction):
                    loc = int(frloc[1:])
                    vals[loc] += count
        xs = []
        ys = []
        for loc in range(max(vals) + 1):
            xs.append(loc)
            ys.append(vals[loc])

        ax.plot(xs, ys, label=direction)
        ax.legend()
        ax.set_title("%s" % slug)

fig.savefig("read-start-positions.png", dpi=180)
plt.clf()

fig, axs = plt.subplots(#nrows=9, ncols=1,
                        #figsize=(4*1, 4*9),
                        nrows=3, ncols=3,
                        figsize=(12, 12),
                        constrained_layout=True)

plt.suptitle("fraction of reads observed starting at each location")
fig.supxlabel("read start locations, most to least frequent")
fig.supylabel("cumulative probility")
for row, vid in enumerate(all_vids):
    slug = to_slug(vid)
    ax = axs[row%3][row//3]
    #ax = axs[row]

    for direction in ['f', 'r']:
        vals = Counter()
        max_loc = 0
        for accession in samples:
            for frloc, count in samples[accession][vid].items():
                if frloc.startswith(direction):
                    loc = int(frloc[1:])
                    max_loc = max(max_loc, loc)
                    vals[loc] += count

        xs = []
        ys = []
        total = sum(vals.values())
        cumulative = 0
        percentiles = set([50, 80, 95])
        for i, val in enumerate(sorted(
                [vals[loc] for loc in range(max_loc + 1)], reverse=True)):
            cumulative += val
            ys.append(cumulative / total)
            xs.append(i)

            for percentile in list(percentiles):
                if cumulative / total * 100 > percentile:
                    print("%.1f%%\t%s\t%.0f" % (100 * i / max_loc, slug,
                                              cumulative / total * 100)) 
                    percentiles.remove(percentile)

        ax.plot(xs, ys, label=direction)
        ax.legend()
        ax.set_title("%s" % slug)

fig.savefig("read-start-positions-cdf.png", dpi=180)
plt.clf()
