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

for slug, d in viruses.items():
    for c, k in d.items():
        plt.close()
        plt.clf()
        fig, ax = plt.subplots()

        ax.set_ylabel("reads")
        ax.set_xlabel("position along genome")
        ax.set_title("%s: reads by location (%s)" % (slug.upper(), c))
        loc_percentages = [x / len(k["loc_full"]) * 100
                           for x in range(len(k["loc_full"]))]

        for tag, label in labels:
            ax.plot(loc_percentages, k[tag], label=label)

        ax.xaxis.set_major_formatter(mtick.PercentFormatter())
        ax.legend()
        fig.savefig("%s-%s-matches.png" % (slug, c), dpi=180)

    for tag, label in labels:
        plt.close()
        plt.clf()
        fig, ax = plt.subplots()

        ax.set_ylabel("reads")
        ax.set_xlabel("position along genome")
        ax.set_title("%s: reads by location (%s)" % (slug.upper(), label))

        for c, k in d.items():
            loc_percentages = [x / len(k[tag]) * 100
                               for x in range(len(k[tag]))]
            ax.plot(loc_percentages, k[tag], label=c)

        ax.xaxis.set_major_formatter(mtick.PercentFormatter())
        ax.legend()
        fig.savefig("%s-%s-matches.png" % (slug, tag), dpi=180)
        

plt.close()
plt.clf()
tag="loc_full"
tagname="whole reads where any k-mer matches"
fig, ax = plt.subplots(nrows=len(viruses), sharex=True, figsize=(6,12),
                       constrained_layout=True)
#fig.tight_layout()
plt.suptitle("reads by location (%s)" % tagname)
fig.supxlabel("scaled reads")
fig.supylabel("position along genome")
for i, (slug, d) in enumerate(sorted(viruses.items())):
    k = d["assembled"]
    loc_percentages = [x / len(k[tag]) * 100
                       for x in range(len(k[tag]))]

    #ax[i].grid(visible=True, axis='y')
    ax[i].plot(loc_percentages, k[tag], label=slug)
    ax[i].xaxis.set_major_formatter(mtick.PercentFormatter())
    ax[i].ticklabel_format(axis='y', style='sci')
    ax[i].set_title(slug)
    fig.savefig("all-matches.png", dpi=180)
