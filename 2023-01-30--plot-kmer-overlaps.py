#!/usr/bin/env python3

import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from Bio.SeqIO.FastaIO import SimpleFastaParser


def rc(s):
    return "".join({'T':'A', 'G':'C', 'A':'T', 'C':'G', 'N':'N'}[x]
                   for x in reversed(s))


def compare(fasta1, fasta2):
    vid1 = fasta1.replace(".fasta", "")
    vid2 = fasta2.replace(".fasta", "")

    with open(fasta1) as inf:
        (_, seq1), = SimpleFastaParser(inf)
    with open(fasta2) as inf:
        (_, seq2), = SimpleFastaParser(inf)

    fig = plt.figure()
    ax = fig.add_subplot()

    ax.set_xlabel(vid1)
    ax.set_ylabel(vid2)
        
    rc_seq1 = rc(seq1)
    rc_seq2 = rc(seq2)
   
    midpoint1 = len(seq1)
    midpoint2 = len(seq2)

    scalar = 32

    def to_pip(v):
        return int(v/scalar)
       
    data = np.ndarray((to_pip(midpoint1*2) + 1,
                       to_pip(midpoint2*2) + 1),
                      dtype=np.uint)

    for s1, s2, d1, d2 in [
            (seq1, seq2, 1, 1),
            (seq1, rc_seq2, 1, -1),
            (rc_seq1, rc_seq2, -1, -1),
            (rc_seq1, seq2, -1, 1),
            ]:
        s1_len = len(s1)
        s2_len = len(s2)
        for i1 in range(s1_len):
            if i1%100 == 0:
                print ("%s x %s: %s / %s" % (d1, d2, i1, s1_len))
            for i2 in range(s2_len):
                c1 = s1[i1]
                c2 = s2[i2]

                match = 0
                p1 = i1
                p2 = i2
            
                while c1 == c2:
                    match += 1
                    p1 += 1
                    p2 += 1
                    if p1 >= s1_len or p2 >= s2_len:
                        break
                    c1 = s1[p1]
                    c2 = s2[p2]

                if match > 0:
                    p1 = i1
                    p2 = i2
                    while True:
                        p1 -= 1
                        p2 -= 1
                        if p1 < 0 or p2 < 0:
                            break
                        c1 = s1[p1]
                        c2 = s2[p2]
                        if c1 != c2:
                            break
                        match += 1

                data[to_pip(midpoint1 + d1*i1)]\
                    [to_pip(midpoint2 + d2*i2)] += match

    axis_labels = ["-100%", "-50%", "0%", "50%", "100%"]
    axis_locations = np.arange(len(axis_labels))
    x_axis_locations = [
        v * len(data[0]) / max(axis_locations) for v in axis_locations]
    y_axis_locations = [
        v * len(data) / max(axis_locations) for v in axis_locations]

    ax.matshow(data)
    ax.set_xticks(x_axis_locations)
    ax.set_xticklabels(axis_labels)
    ax.set_yticks(y_axis_locations)
    ax.set_yticklabels(axis_labels)
 
    plt.gca().invert_yaxis()  # want 100% x 100% in upper right
    fig.savefig("matrix-%s-%s.png" % (vid1, vid2), dpi=180)

compare(*sys.argv[1:])


                            
