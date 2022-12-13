#!/usr/bin/env python3

import sys
import glob
import gzip
import json
from collections import defaultdict
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio import Align

fwds = defaultdict(set)
revs = defaultdict(set)

aligner = Align.PairwiseAligner()
# These are the scoring settings porechop uses by default.
# https://github.com/rrwick/Porechop/blob/master/porechop/porechop.py#L145
aligner.end_gap_score = 0
aligner.match_score = 3
aligner.mismatch_score = -6
aligner.internal_open_gap_score = -5
aligner.internal_extend_gap_score = -2

def rc(s):
    return "".join({'T':'A',
                    'G':'C',
                    'A':'T',
                    'C':'G',
                    'N':'N'}[x] for x in reversed(s))

match_table = []
with open('anjali-table.txt') as inf:
    for line in inf:
        assay_name, probe, fwd, rev = line.strip().split()
        match_table.append((assay_name, probe, fwd, rev))

fname, = sys.argv[1:]

with gzip.open(fname, mode='rt') as inf:
    with open(fname + ".matched-contents", "w") as outf:
        for (title, sequence, quality) in FastqGeneralIterator(inf):
            out = []
            for assay_name, probe, fwd, rev in match_table:
                for is_rc in [0, 1]:
                    probe_alignment = aligner.align(
                        rc(probe) if is_rc else probe, sequence)[0]
                    fwd_alignment = aligner.align(
                        rc(fwd) if is_rc else fwd, sequence)[0]
                    rev_alignment = aligner.align(
                        rc(rev) if is_rc else rev, sequence)[0]

                    out.extend((
                        "%.1f" % (probe_alignment.score / len(probe)),
                        "%.1f" % (fwd_alignment.score / len(fwd)),
                        "%.1f" % (rev_alignment.score / len(rev))))
            outf.write("\t".join(out) + "\n")
