#!/usr/bin/env python3

from Bio import Align
from collections import Counter, defaultdict
from Bio.SeqIO.FastaIO import SimpleFastaParser

# format: label, start (one-indexed), end (one-indexed, inclusive)
references = {
    "pHR'-CMVLacZ": [
        [1, 697],
        [7168, 8025],
        [8557, 8629],
        [8630, 9181],
    ],
    "pCDH-EF1a-GaussiaSP-MCS-IRES-zeo": [
        [1, 688],
        [4327, 4444],
        [7168, 8025],
        [8561, 8666],
        [9059, 9181],
    ],
    "pLVX.TRE3G.eGFP": [
        [1, 688],
        [4328, 4478],
        [7168, 8025],
        [8444, 9181],
    ],
    "pLKO.1-ZsGreen1": [
        [1, 697],
        [4303, 4480],
        [7168, 8025],
        [8556, 8666],
        [9067, 9181],
    ],
    "pLenti6-V5-Puro": [
        [1, 697],
        [7168, 8025],
        [8557, 8614],
        [8615, 8666],
        [9067, 9181],
    ],
    "pHR-hSyn-EGFP": [
        [1, 688],
        [4204, 4481],
        [7168, 8025],
        [8557, 8666],
        [8632, 9181],
    ],
    "FUW-tetO-MCS": [
        [1, 697],
        [4298, 4480],
        [7168, 8025],
        [8443, 8944],
        [9057, 9181],
    ],
    "pRRL-sffv-eGFP-cmv-hsGDNF": [
        [1, 697],
        [4327, 4444],
        [7168, 8025],
        [8557, 8666],
        [9067, 9181],
    ],

}

reads = []
with open("p2ra_hiv_reads.txt") as inf:
    for line in inf:
        bits = line.rstrip("\n").split("\t")
        if bits[0] == "study_name":
            continue

        reads.append((bits[20], # fwd
                      bits[21], # rev
                      tuple(bits)))

with open("hiv.fasta") as inf:
    for title, reference_seq in SimpleFastaParser(inf):
        assert title == "AF033819.3 HIV-1, complete genome"

aligner = Align.PairwiseAligner()
aligner.end_gap_score = 0
aligner.match_score = 3
aligner.mismatch_score = -5
aligner.internal_open_gap_score = -10
aligner.internal_extend_gap_score = -2

def align(qry, ref):
    return aligner.align(qry, ref)[0]

studies = set()

ys_counter = defaultdict(Counter)
for fwd, rev, bits in reads:
    study = bits[0]
    observed_positions = set()
    for read in [fwd, rev]:
        for (qry_beg, qry_end), (ref_beg, ref_end) in zip(
                *align(read, reference_seq).aligned):
            for i in range(ref_end - ref_beg + 1):
                observed_positions.add(ref_beg + i)
    for observed_position in observed_positions:
        studies.add(study)
        ys_counter[study][observed_position] += 1

import matplotlib.pyplot as plt


study_segments = defaultdict(list)
for study in studies:
    xs = []
    ys = []
    study_segments[study].append((xs, ys))

    for x in range(len(reference_seq)):
        if ys_counter[study][x]:
            xs.append(x)
            y = ys_counter[study][x]
            ys.append(y)
        elif xs:
            xs = []
            ys = []
            study_segments[study].append((xs, ys))

fig, axs = plt.subplots(figsize=(12,(len(studies)+1)*2),
                       nrows=1+len(studies),
                       ncols=1,
                       sharex=True)
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
for color, study, ax in zip(colors, sorted(studies), fig.axes):
    ax.set_title(study)
    ax.set_ylabel("Observations")
    for xs, ys in study_segments[study]:
        ax.plot(xs, ys, color=color)

height = 1/len(references)
colors = ["pink",
          "lightgreen",
          "lightblue",
          "lightyellow",
          "oldlace",
          "lavender",
          "lightgrey",
          "azure",
          "honeydew"]

labels = []
ticks = []

for i, (color, (label, regions)) in enumerate(
        zip(colors,
            sorted(references.items()))):
    y_position = i / len(references)
    labels.append(label)
    ticks.append(y_position + height/2)
    for region in regions:
        if len(region) == 2:
            start, end_inclusive = region
            individual_label = None
        else:
            start, end_inclusive, individual_label = region

        length = end_inclusive - start + 1

        # convert to zero indexing
        start -= 1
        end_inclusive -= 1

        axs[-1].add_patch(plt.Rectangle(
            (start, y_position),
            length,
            height, facecolor=color,
            edgecolor='black'))

        if length > 100 and individual_label:
            axs[-1].text(start + length/2, y_position + (height/2),
                    individual_label,
                    ha='center', va='center')

axs[-1].set_yticks(ticks)
axs[-1].set_yticklabels(labels)

fig.supxlabel("Position along AF033819.3 HIV-1 genome")
fig.savefig("p2ra-hiv-coverage.png", dpi=180)

plt.clf()
