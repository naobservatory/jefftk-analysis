#!/usr/bin/env python3

import sys
import glob
import subprocess
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

fasta_pairs = sys.argv[1:]
pairs = [x.split(":") for x in fasta_pairs]
decorated = [] # average_alignment, pair, alignment

for pair in pairs:
    fasta1, fasta2 = pair
    moving_average_tsv = subprocess.check_output([
        "/Users/jeffkaufman/code/sequencing_tools/align.py",
        "--moving-average-tsv", "--min-score=-10000", fasta1, fasta2])
    x= []
    for line in moving_average_tsv.split(b"\n"):
        line = line.strip()
        if not line: continue
        position, match, average = line.split()
        if average == b'average': continue
        x.append(float(average) * 100)

    decorated.append((sum(x)/len(x), pair, x))

if len(decorated) > 1:
    fig, ax = plt.subplots(nrows=len(pairs), sharex=True,
                           figsize=(6,2*len(decorated)),
                           constrained_layout=True)
    plt.suptitle("Agreement along genome: assembled vs closest GenBank match")
    fig.supxlabel("pseudoposition")
    fig.supylabel("agreement")
else:
    fig, ax = plt.subplots()
    ax.set_xlabel("pseudoposition")
    ax.set_ylabel("agreement")

for i, (_, (fasta1, fasta2), x) in enumerate(decorated):
    vid1 = fasta1.replace(".fasta", "")
    vid2 = fasta2.replace(".fasta", "")

    label = vid1 + "-" + vid2

    loc_percentages = [v / len(x) * 100 for v in range(len(x))]

    axi = ax[i] if len(decorated) > 1 else ax

    axi.plot(loc_percentages, x)
    axi.xaxis.set_major_formatter(mtick.PercentFormatter())
    axi.yaxis.set_major_formatter(mtick.PercentFormatter())
    axi.set_ylim([-5, 105])
    axi.set_title(label)

fig.savefig("alignment-pairs-2.png", dpi=180)
