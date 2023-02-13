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
                        nrows=6, ncols=3,
                        figsize=(12, 24),
                        constrained_layout=True)

plt.suptitle("read starts by location, per sample")
fig.supylabel("position along genome in this sample")
fig.supxlabel("reads starting at this position")
for direction in 'fr':
    for row, vid in enumerate(all_vids):

        slug = to_slug(vid)
        ax = axs[(row%3)*2 + (direction == 'f')][row//3]


        loc_vals = Counter()
        max_loc = 0
        for accession in samples:
            for frloc, count in samples[accession][vid].items():
                if frloc.startswith(direction):
                    loc = int(frloc[1:])
                    max_loc = max(max_loc, loc)
                    loc_vals[loc] += count

        loc_ranks = {} # loc -> rank
        for rank, (val, loc) in enumerate(sorted(
                [(loc_vals[loc], loc) for loc in range(max_loc + 1)],
                reverse=True)):
            loc_ranks[loc] = rank

        offset = 0
        for accession in samples:
            points = []
            vals = Counter()
            for frloc, count in samples[accession][vid].items():
                if frloc.startswith(direction):
                    loc = int(frloc[1:])
                    if loc < 0:
                        continue
                    points.append((loc_ranks[loc], count))

            if len(points) < max_loc / 4:
                continue

            points.sort()

            xs = [x for (x,y) in points]
            ys = [y for (x,y) in points]
            ys = [y / sum(ys) for y in ys]

            xs = [x + offset * (max(xs)/20) for x in xs]
            ys = [y + offset * (max(ys)/20) for y in ys]
            
            ax.plot(xs, ys)

            offset += 1

        ax.set_title("%s:%s" % (slug, direction))

fig.savefig("read-start-positions-by-sample.png", dpi=180)
plt.clf()
