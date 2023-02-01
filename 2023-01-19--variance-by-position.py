#!/usr/bin/env python3

# First:
#   scp 'REMOTE:kmer-egd/SEQUENCE.*' .

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as mtick
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

wtp_means = defaultdict(list)
wtp_variances = defaultdict(list)
wtp_locs = defaultdict(list)

wtp_norm_means = defaultdict(list)
wtp_norm_variances = defaultdict(list)
wtp_norm_locs = defaultdict(list)

for wtp in wtps:
    sample_sums = None
    for prefix in "ACGT":
        with open("%s.%s.%s" % (sequence, wtp, prefix)) as inf:
            for line in inf:
                kmer, *counts = line.strip().split("\t")

                if sample_sums is None:
                    sample_sums = [0]*len(counts)
                
                for i, x in enumerate(counts):
                    sample_sums[i] += int(x)

    for prefix in "ACGT":
        with open("%s.%s.%s" % (sequence, wtp, prefix)) as inf:
            for line in inf:
                kmer, *counts = line.strip().split("\t")
                vals = [int(x) for x in counts]
                n_vals = len(vals)

                mean = sum(vals) / n_vals
                squared_errors = [(x - mean)*(x-mean) for x in vals]
                unbiased_variance = sum(squared_errors) / (n_vals-1)

                wtp_means[wtp].append(mean)
                wtp_variances[wtp].append(unbiased_variance)
                wtp_locs[wtp].append(loc[kmer])

                vals = [x / sample_sums[i] for (i,x) in enumerate(vals)]

                mean = sum(vals) / n_vals
                squared_errors = [(x - mean)*(x-mean) for x in vals]
                unbiased_variance = sum(squared_errors) / (n_vals-1)

                wtp_norm_means[wtp].append(mean)
                wtp_norm_variances[wtp].append(unbiased_variance)
                wtp_norm_locs[wtp].append(loc[kmer])

 
for wtp in wtps:
    fig = plt.figure()
    ax = fig.add_subplot()
    
    ax.set_title("%s tbrv per-kmer variance as a function of abundance" % wtp)
    ax.set_xlabel("kmer abundance")
    ax.set_ylabel("unbiased variance")
    ax.scatter(wtp_means[wtp],
               wtp_variances[wtp],
               s=[0.1]*len(wtp_means[wtp]),
               marker=',')
    fig.savefig("tbrv-%s-variance.png" % wtp, dpi=180)

for wtp in wtps:
    fig = plt.figure()
    ax = fig.add_subplot()
    
    ax.set_title("%s tbrv per-kmer variance as a function of abundance" % wtp)
    ax.set_xlabel("kmer abundance")
    ax.set_ylabel("unbiased variance after adjusting for total tbrv abundance")
    ax.scatter(wtp_norm_means[wtp],
               wtp_norm_variances[wtp],
               s=[0.1]*len(wtp_norm_means[wtp]),
               marker=',')
    fig.savefig("tbrv-%s-norm-variance.png" % wtp, dpi=180)

#raise Exception("breakpoint")
                
for wtp in wtps:
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    
    ax.set_title("%s tbrv per-kmer variance as a function of abundance" % wtp)
    ax.set_xlabel("kmer abundance")
    ax.set_ylabel("unbiased variance")
    ax.set_zlabel("location in genome")
    ax.zaxis.set_major_formatter(mtick.PercentFormatter())
    ax.scatter(wtp_norm_means[wtp],
               wtp_norm_variances[wtp],
               [x/(n_kmers-1)*100 for x in wtp_norm_locs[wtp]],
               s=[0.3]*len(wtp_norm_means[wtp]),
               marker=',')
    #fig.savefig("tbrv-%s-variance.png" % wtp, dpi=180)

plt.show()

