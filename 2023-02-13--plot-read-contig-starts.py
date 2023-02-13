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
all_vids = [x for x in all_vids if "tbrv" in x]

def to_slug(vid):
    return vid.split(".")[0]

vid, = all_vids
slug = to_slug(vid)

totals = {'f': Counter(), 'r': Counter()}

for sample in samples:
    for direction in 'fr':
        for frloc, count in samples[sample][vid].items():
            if frloc.startswith(direction):
                totals[direction][sample] += count

sample_sizes = [(totals['f'][sample] + totals['r'][sample], sample)
                for sample in samples]
sample_sizes.sort(reverse=True)
    
                
N_ROWS=8
N_COLS=5
assert len(samples) == N_ROWS * N_COLS

fig, axs = plt.subplots(nrows=N_ROWS, ncols=N_COLS,
                        sharex=True, sharey=True,
                        figsize=(3*N_COLS, 3*N_ROWS),
                        constrained_layout=True)

plt.suptitle("predicting per-sample starts from other samples")
fig.supxlabel("relative frequency of this location in other samples")
fig.supylabel("relative frequency of this location in this sample")

for i, (sample_size, sample) in enumerate(sample_sizes):
    ax = axs[i//N_COLS][i%N_COLS]
    ax.set_xlim(xmax=2, xmin=0)
    ax.set_ylim(ymax=2, ymin=0)
    ax.xaxis.set_major_formatter(mtick.PercentFormatter())
    ax.yaxis.set_major_formatter(mtick.PercentFormatter())

    for direction in ['f', 'r']:
        # Determine relative frequencies of these locations in other samples
        other_loc_counts = Counter() # loc -> count
        this_loc_counts = Counter()  # loc -> count
        for other_sample in samples:
            if other_sample != sample:
                for frloc, count in samples[other_sample][vid].items():
                    if frloc.startswith(direction):
                        loc = int(frloc[1:])
                        other_loc_counts[loc] += (
                            count / totals[direction][other_sample]
                                  / (len(samples)-1))
                        
        for frloc, count in samples[sample][vid].items():
            if frloc.startswith(direction):
                loc = int(frloc[1:])
                this_loc_counts[loc] += (
                    count / totals[direction][sample])

        all_locs = set(other_loc_counts)
        all_locs.update(this_loc_counts)

        xs = []
        ys = []
        for loc in all_locs:
            xs.append(other_loc_counts[loc] * 100)
            ys.append(this_loc_counts[loc] * 100)

        ax.scatter(xs, ys, label=direction, marker=',', s=0.5)

    ax.legend()
    ax.set_title("%s n=%s" % (sample, sample_size))

fig.savefig("read-starts-tbrv.png", dpi=180)
plt.clf()
