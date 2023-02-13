#!/usr/bin/env python3

import sys
import json
import datetime
import itertools
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

with open("internal-counts.json") as inf:
    internal_counts = json.load(inf)

with open("prefix-counts.json") as inf:
    prefix_counts = json.load(inf)
    
fig, axs = plt.subplots(nrows=4,
                        ncols=2,
                        sharex=True, sharey=True,
                        figsize=(6,12),
                        constrained_layout=True)

plt.suptitle("patterns in read prefixes")
fig.supxlabel("kmers, most to least frequent")
fig.supylabel("proportion of data covered")

for i in range(8):
    k = i+1
    internal = internal_counts[i]
    prefixes = prefix_counts[i]

    ax = axs[i//2][i%2]

    for data, label in [
            [internal_counts[i], "internal"],
            [prefix_counts[i], "prefix"]]:
        xs = []
        ys = []

        ax.xaxis.set_major_formatter(mtick.PercentFormatter())
        ax.yaxis.set_major_formatter(mtick.PercentFormatter())

        raw = [(count, kmer)
               for (kmer, count) in data.items()
               if 'N' not in kmer]
        raw.sort(reverse=True)

        total = sum(data.values())
        cumulative = 0
        
        for n, (count, kmer) in enumerate(raw):            
            cumulative += count
            xs.append(100 * n / (4**k-1))
            ys.append(100 * cumulative / total)

        ax.plot(xs, ys, label=label)
    ax.plot([0,100], [0,100], color='black', label="x=y")
        
    ax.legend()
    ax.set_title("k=%s" % k)

fig.savefig("kmer-prefix-frequency.png", dpi=180)

