#!/usr/bin/env python3

import glob
import json
import numpy as np
import matplotlib.pyplot as plt

def phred_to_q(phred_score):
    return ord(phred_score) - ord('!')

def average_quality(phred_counts):
    xs = []
    weights = []

    for phred_score, count in phred_counts.items():
        xs.append(phred_to_q(phred_score))
        weights.append(count)

    return np.average(xs, weights=weights)

data = {} # accession -> slug -> xs, ys
for fname in glob.glob("*.json"):
    accession, _ = fname.split(".")

    with open(fname) as inf:
        qc = json.load(inf)

        if "post_cleaning" not in qc: continue
        post_cleaning = qc["post_cleaning"]
        if "qualities" not in post_cleaning: continue
        qualities = post_cleaning["qualities"]

        total = sum(
            sum(counts.values())
            for slug in qualities
            for counts in qualities[slug])

        data[accession] = {}
        for slug in qualities:
            vals = qualities[slug]

            # ignore categories that cover very little data
            if sum(sum(val.values()) for val in vals) < total * 0.001: continue
            
            xs = [pos for (pos, _) in enumerate(vals)]
            ys = [average_quality(val) for val in vals]

            data[accession][slug] = xs, ys

for accession in sorted(data):
    fig, ax = plt.subplots(constrained_layout=True)
    for slug in sorted(data[accession]):
        xs, ys = data[accession][slug]
        ax.plot(xs, ys, label=slug)
    ax.legend()
    fig.savefig("%s.qc.q.png" % accession)
    plt.close()

slugs = sorted(set(
    slug
    for accession in data
    for slug in data[accession]))

for slug in slugs:
    fig, ax = plt.subplots(constrained_layout=True)
    for accession in sorted(data):
        xs, ys = data[accession][slug]
        ax.plot(xs, ys, label=accession)
    ax.legend()
    fig.savefig("multi.%s.qc.q.png" % slug)
    plt.close()
