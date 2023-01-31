#!/usr/bin/env python3

import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from Bio.SeqIO.FastaIO import SimpleFastaParser

K=40

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

    for c, k in d.items():
        with open(k["fasta"]) as inf:
            (_, seq), = SimpleFastaParser(inf)
            k['seq'] = seq
    
    viruses[slug] = d

def kmers(seq, offset=0):
    for i in range(len(seq) - K + 1):
        if i < offset: continue
        yield seq[i:i+K]

def rc(s):
    return "".join({'T':'A', 'G':'C', 'A':'T', 'C':'G', 'N':'N'}[x]
                   for x in reversed(s))


def compare(slug1, slug2):
    fig = plt.figure()
    ax = fig.add_subplot()

    ax.set_xlabel(slug2)
    ax.set_ylabel(slug1)
    
    seq1 = viruses[slug1]["assembled"]["seq"]
    seq2 = viruses[slug2]["assembled"]["seq"]

    seq1_f_kmers = list(kmers(seq1))
    seq2_f_kmers = list(kmers(seq2))

    seq1_r_kmers = list(kmers(rc(seq1)))
    seq2_r_kmers = list(kmers(rc(seq2)))

    seq1_kmers = seq1_r_kmers + seq1_f_kmers
    seq2_kmers = seq2_r_kmers + seq2_f_kmers

   
    data = np.ndarray((len(seq1_kmers),len(seq2_kmers)), dtype=bool)

    for i in range(len(seq1_kmers)):
        for j in range(len(seq2_kmers)):
            data[i][j] = seq1_kmers[i] == seq2_kmers[j]

    axis_labels = ["-100%", "-50%", "0%", "50%", "100%"]
    axis_locations = np.arange(len(axis_labels))
    x_axis_locations = [
        x * len(seq2_kmers) / max(axis_locations) for x in axis_locations]
    y_axis_locations = [
        x * len(seq1_kmers) / max(axis_locations) for x in axis_locations]
            
    ax.matshow(data)
    ax.set_xticks(x_axis_locations)
    ax.set_xticklabels(axis_labels)
    ax.set_yticks(y_axis_locations)
    ax.set_yticklabels(axis_labels)

    fig.savefig("%s-%s-matrix.png" % (slug1, slug2), dpi=180)

compare('tbrv', 'tmv')


                            
