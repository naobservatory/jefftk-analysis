#!/usr/bin/env python3

# First:
#   scp 'REMOTE:kmer-egd/SEQUENCE.*' .

import sys
from collections import defaultdict, Counter

# ex: tomato.brown.rugose.nt.seq longest-timeseries.tsv
sequence, metadata = sys.argv[1:]
K=40

def rc(s):
    return "".join({'T':'A', 'G':'C', 'A':'T', 'C':'G', 'N':'N'}[x]
                   for x in reversed(s))

metadatas = defaultdict(list)
with open(metadata) as inf:
    for line in inf:
        accession, date, wtp = line.strip().split("\t")
        metadatas[wtp].append(date)

wtps = list(sorted(metadatas))
        
# Locations are relative to the forward sequence, and are counted from the end
# of the kmer closest to the start of the sequence.
loc = {}  # kmer -> loc

with open(sequence) as inf:
    line, = inf
    seq = line.strip()
    n_kmers = len(seq) - K + 1
    for i in range(n_kmers):
        kmer_in = seq[i:i+K]
        loc[kmer_in] = i
        loc[rc(kmer_in)] = i

# loc -> wtp -> [count1, count2, ...]
loc_vals = defaultdict(lambda: defaultdict(Counter))

for wtp in wtps:
    for prefix in "ACGT":
        with open("%s.%s.%s" % (sequence, wtp, prefix)) as inf:
            for line in inf:
                kmer, *counts = line.strip().split("\t")
                record = loc_vals[loc[kmer]][wtp]
                for i, count in enumerate(counts):
                    record[i] += int(count)

cols = {}
wtp_cols = {}

for wtp, dates in sorted(metadatas.items()):
    for date_index, date in enumerate(dates):
        wtp_cols[wtp, date] = [
            loc_vals[i][wtp][date_index]
            for i in range(max(loc_vals))]
        if wtp not in cols:
            cols[wtp] = [0]*max(loc_vals)
        for i in range(max(loc_vals)):
            cols[wtp][i] += loc_vals[i][wtp][date_index]

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mtick

colors = {'HTP': 'b', 'JWPCP': 'g', 'OC': 'r', 'SJ': 'c'}
style = {
    'HTP': 'solid',
    'JWPCP': 'dashed',
    'OC': 'dashdot',
    'SJ': 'dotted',
}


fig, ax = plt.subplots()

global_mean_col = [0]*n_kmers
for wtp, col in sorted(cols.items()):
    loc_percentages = [x / (n_kmers-1) * 100 for x in range(len(col))]
    adj_col = [x/sum(col)*len(col) for x in col]
    for i, x in enumerate(adj_col):
        global_mean_col[i] += x
    
    ax.plot(loc_percentages, adj_col,
            label=wtp,
            linewidth=0.75,
            color=mcolors.BASE_COLORS[colors[wtp]],
            linestyle=style[wtp])
global_mean_col = [x/len(wtps) for x in global_mean_col]

ax.set_ylabel("tbrv kmer relative abundance, average of each site's samples")
ax.set_xlabel("position along genome")
ax.set_title("per-wtp tbrv relative abundance by position")
ax.xaxis.set_major_formatter(mtick.PercentFormatter())
ax.legend()

fig.savefig("tbrv-wtp-averages.png", dpi=180)

plt.clf()
fig, ax = plt.subplots()
wtp_mean_cols = {}
for wtp, col in sorted(cols.items()):
    loc_percentages = [x / (n_kmers-1) * 100 for x in range(len(col))]
    adj_col = [x/sum(col)*len(col) for x in col]
    wtp_mean_cols[wtp] = adj_col

    adj_col = [x / global_mean_col[i] for i, x in enumerate(adj_col)]
    
    ax.plot(loc_percentages, adj_col,
            label=wtp,
            linewidth=0.75,
            color=mcolors.BASE_COLORS[colors[wtp]],
            linestyle=style[wtp])
ax.set_ylabel("tbrv kmer normalized abundance, site mean vs global mean")
ax.set_xlabel("position along genome")
ax.set_title("per-wtp tbrv normalized by position")
ax.xaxis.set_major_formatter(mtick.PercentFormatter())
ax.legend()
fig.savefig("tbrv-wtp-normalized.png", dpi=180)

for wtp in wtps:
    plt.clf()
    fig, ax = plt.subplots()

    ax.set_title("%s tbrv relative abundance by position" % wtp)

    for (wtp_2, date), col in sorted(wtp_cols.items()):
        if wtp_2 != wtp: continue
        
        loc_percentages = [x / (n_kmers-1) * 100 for x in range(len(col))]
        adj_col = [x/sum(col)*len(col) for x in col]
        
        ax.scatter(loc_percentages,
                   adj_col,
                   label=date,
                   marker=',',
                   #color=mcolors.BASE_COLORS[colors[wtp]],
                   s=[0.01]*len(adj_col),
                   )

    ax.set_xlabel("position along genome")
    ax.set_ylabel("tbrv kmer relative abundance")
    ax.xaxis.set_major_formatter(mtick.PercentFormatter())
    ax.legend(markerscale=25)

    fig.savefig("tbrv-%s-abundance.png" % wtp, dpi=180)


        
for wtp in wtps:
    plt.clf()
    fig, ax = plt.subplots()

    ax.set_title("%s tbrv normalized abundance by position" % wtp)

    for (wtp_2, date), col in sorted(wtp_cols.items()):
        if wtp_2 != wtp: continue
        
        loc_percentages = [x / (n_kmers-1) * 100 for x in range(len(col))]
        adj_col = [x/sum(col)*len(col) for x in col]
        adj_col = [x / wtp_mean_cols[wtp][i]
                   if wtp_mean_cols[wtp][i] > 0 else 0
                   for i, x in enumerate(adj_col)]
        
        ax.scatter(loc_percentages,
                   adj_col,
                   label=date,
                   marker=',',
                   s=[0.01]*len(adj_col),
                   )

    ax.set_xlabel("position along genome")
    ax.set_ylabel("tbrv kmer normalized (delta squared), sample vs site mean")
    ax.xaxis.set_major_formatter(mtick.PercentFormatter())
    ax.legend(markerscale=25, ncol=3, loc='lower center')
    ax.set_ylim(ymax=3, ymin=-1)

    fig.savefig("tbrv-%s-normalized.png" % wtp, dpi=180)


        
