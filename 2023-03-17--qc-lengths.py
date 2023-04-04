#!/usr/bin/env python3

import glob
import json
import numpy as np
import matplotlib.pyplot as plt

data = {} # accession -> xs, ys
for fname in glob.glob("*.json"):
    accession, _ = fname.split(".")

    with open(fname) as inf:
        qc = json.load(inf)

        if "post_cleaning" not in qc: continue
        post_cleaning = qc["post_cleaning"]
        if "lengths" not in post_cleaning: continue
        lengths = post_cleaning["lengths"]

        max_length = max(int(length) for length in lengths)
        
        xs = []
        ys = []
        for length in range(max_length+1):
            xs.append(length)
            ys.append(lengths.get(str(length), 0))

        data[accession] = xs, ys

for accession, (xs, ys) in sorted(data.items()):
    fig, ax = plt.subplots(constrained_layout=True)
    ax.plot(xs, ys)
    fig.savefig("%s.qc.png" % accession)
    plt.close()

fig, ax = plt.subplots(constrained_layout=True)
big_xs = np.array([])
big_ys = np.array([])
for accession, (xs, ys) in sorted(data.items()):
    xs = np.array(xs)
    ys = np.array(ys)

    print("%s: %s" % (accession, np.average(xs, weights=ys)))
    
    if len(xs) > len(big_xs):
        big_xs.resize(xs.shape)
        big_ys.resize(xs.shape)
    elif len(xs) < len(big_xs):
        xs.resize(big_xs.shape)
        ys.resize(big_xs.shape)

    big_xs += xs
    big_ys += ys

ax.plot(big_xs, big_ys)
fig.savefig("overall.qc.png")
plt.close()

fig, ax = plt.subplots(constrained_layout=True)
for accession, (xs, ys) in sorted(data.items()):
    scaled_ys = np.array(ys) / sum(ys)
    ax.plot(xs, scaled_ys, label=accession)
ax.legend()
fig.savefig("multi.qc.png")
plt.close()
