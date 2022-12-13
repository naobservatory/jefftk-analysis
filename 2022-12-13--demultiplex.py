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

with open('barcodes.tsv') as inf:
    for line in inf:
        experimenter, sample_id, sample_name, fwd, rev = line.strip().split()
        fwd = fwd.upper()
        rev = rev.upper()

        fwds[fwd].add(sample_name)
        revs[rev].add(sample_name)

fname, = sys.argv[1:]

with gzip.open(fname, mode='rt') as inf:
    with open(fname + ".barcoding", "w") as outf:
        for (title, sequence, quality) in FastqGeneralIterator(inf):
            fwd_alignments = []
            for barcode in fwds:
                for is_rc in [0, 1]:
                    alignment = aligner.align(
                        rc(barcode) if is_rc else barcode, sequence)[0]
                    score = alignment.score / len(barcode)
                    fwd_alignments.append(
                        (score, alignment, barcode, is_rc))
            fwd_alignments.sort()
            fwd_best_score, fwd_best_alignment, fwd_best_barcode, fwd_is_rc = fwd_alignments[-1]

            rev_alignments = []
            for barcode in revs:
                for is_rc in [0, 1]:
                    alignment = aligner.align(
                        rc(barcode) if is_rc else barcode, sequence)[0]
                    score = alignment.score / len(barcode)
                    rev_alignments.append(
                        (score, alignment, barcode, is_rc))
            rev_alignments.sort()

            rev_best_score, rev_best_alignment, rev_best_barcode, rev_is_rc = rev_alignments[-1]

            outf.write("%.2f\t%.2f\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                fwd_best_score, rev_best_score,
                fwd_is_rc, rev_is_rc,
                fwd_best_barcode, rev_best_barcode,
                json.dumps(fwd_best_alignment.aligned.tolist()),
                json.dumps(rev_best_alignment.aligned.tolist())))
