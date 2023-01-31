#!/usr/bin/env python3

import sys
import glob

slugs = sys.argv[1:]
viruses = {}
for slug in slugs:
    d = {
        "assembled": {},
        "genbank_top": {},
        "genbank_closest": {},
    }

    d["assembled"]["fasta"], = glob.glob("%s.*.assembled.fasta" % slug)
    d["genbank_top"]["fasta"] = d["assembled"]["fasta"].replace(
        ".assembled.", ".")
    d["genbank_closest"]["fasta"], = (
        set(glob.glob("%s.*.fasta" % slug)) -
        set([d["assembled"]["fasta"], d["genbank_top"]["fasta"]]))

    for k in d.values():
        k["loc_full"] = []
        k["loc_any"] = []
        k["loc_rc"] = []

        for fname in glob.glob(
                "reads-%s/*.locs.tsv" % k["fasta"].replace(".fasta", "")):
            with open(fname) as inf:
                for line in inf:
                    loc, obs_full, obs_any, obs_rc = [
                        int(x) for x in line.strip().split()]

                    if loc >= len(k["loc_full"]):
                        k["loc_full"].append(0)
                        k["loc_any"].append(0)
                        k["loc_rc"].append(0)

                    k["loc_full"][loc] += obs_full
                    k["loc_any"][loc] += obs_any
                    k["loc_rc"][loc] += obs_rc
    viruses[slug] = d

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

labels = [
    ("loc_full", "exact matches only"),
    ("loc_any", "any k-mer matches"),
    ("loc_rc", "whole reads where any k-mer matches"),
]

plt.close()
plt.clf()
fig, ax = plt.subplots()
tag="loc_full"
ax.set_ylabel("reads")
ax.set_xlabel("position along genome")
ax.set_title("reads by location (exact match)")
ax.xaxis.set_major_formatter(mtick.PercentFormatter())

for i, (slug, d) in enumerate(sorted(viruses.items())):
    k = d["assembled"]
    loc_percentages = [x / len(k[tag]) * 100
                       for x in range(len(k[tag]))]

    ax.plot(loc_percentages, k[tag], label=slug.upper())
    ax.legend()
    fig.savefig("overlay-%s.png" % "-".join(sorted(viruses)), dpi=180)
