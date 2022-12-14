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
    with open(fname + ".match.fasta", "w") as outf:
        for (title, sequence, quality) in FastqGeneralIterator(inf):
            sample_names = set()
            matches = set()
            highlights = []
            
            for assay_name, probe, fwd, rev in match_table:
                for is_rc in [0, 1]:
                    for needle, name in [
                            [probe, "probe"],
                            [fwd, "fwd"],
                            [rev, "rev"]]:
                        alignment = aligner.align(
                            rc(needle) if is_rc else needle, sequence)[0]
                        if alignment.score / len(needle) > 2.5:
                            highlights.append(alignment.aligned[-1])
                            sample_names.add(assay_name)
                            matches.add(name)

            if not highlights: continue
            
            sequence = list(sequence)
            for highlight in highlights:
                for begin, end in highlight:
                    for i in range(begin, end):
                        sequence[i] = sequence[i].lower()

            outf.write(">%s %s %s\n%s\n" % (
                title,
                ":".join(sorted(sample_names)),
                ":".join(sorted(matches)),
                "".join(sequence)))

                        

