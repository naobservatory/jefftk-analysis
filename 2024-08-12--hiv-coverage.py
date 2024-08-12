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

ys_counter = Counter()
read_mapper = defaultdict(set)
for fwd, rev, bits in reads:
    observed_positions = set()
    for read in [fwd, rev]:
        for (qry_beg, qry_end), (ref_beg, ref_end) in zip(
                *align(read, reference_seq).aligned):
            for i in range(ref_end - ref_beg + 1):
                observed_positions.add(ref_beg + i)
    for observed_position in observed_positions:
        ys_counter[observed_position] += 1
        read_mapper[observed_position].add(bits)

import matplotlib.pyplot as plt

segments = []
xs = []
ys = []
segments.append((xs, ys))
max_y = 0

covered = set()

for reference in references.values():
    for start, end_inclusive in reference:
        # convert to zero indexing
        start -= 1
        end = end_inclusive
        for i in range(start, end):
            covered.add(i)
        
uncovered = set(range(len(reference_seq))) - covered
        
for x in range(len(reference_seq)):
    if ys_counter[x]:
        xs.append(x)
        y = ys_counter[x]
        max_y = max(max_y, y)
        ys.append(y)
    else:
        if xs:
            xs = []
            ys = []
            segments.append((xs, ys))

uncovered_reads = set()
for uncovered_pos in uncovered:
    uncovered_reads |= read_mapper.get(uncovered_pos, set())
    
for uncovered_read in sorted(uncovered_reads):
    print(*uncovered_read, sep="\t")
            
fig, ax = plt.subplots(figsize=(12,6),
                       nrows=2,
                       ncols=1,
                       sharex=True)

for xs, ys in segments:
    ax[0].plot(xs, ys, color="blue")

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

        ax[1].add_patch(plt.Rectangle(
            (start, y_position),
            length,
            height, facecolor=color,
            edgecolor='black'))

        if length > 100 and individual_label:
            ax[1].text(start + length/2, y_position + (height/2),
                    individual_label,
                    ha='center', va='center')

ax[1].set_yticks(ticks)
ax[1].set_yticklabels(labels)

fig.supxlabel("Position along AF033819.3 HIV-1 genome")
ax[0].set_ylabel("Observations")

fig.savefig("p2ra-hiv-coverage.png", dpi=180)

plt.clf()
