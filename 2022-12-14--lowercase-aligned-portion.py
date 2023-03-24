#!/usr/bin/env python3

import json
import gzip
import glob
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def rc(s):
    return "".join({'T':'A',
                    'G':'C',
                    'A':'T',
                    'C':'G',
                    'N':'N'}[x] for x in reversed(s))

sample_barcodes = {}
with open("laura-sample-barcodes.txt") as inf:
    for line in inf:
        sample_name, fwd, rev = line.strip().split()
        sample_barcodes[sample_name] = fwd, rev

out_generous = []
out_strict = []
for fname in glob.glob("*.barcoding"):
    with open(fname) as inf:
        for lineno, line in enumerate(inf):
            (score_fwd, score_rev,
             fwd_is_rc, rev_is_rc,
             fwd_found, rev_found,
             fwd_alignment, rev_alignment) = line.strip().split("\t")

            sample_names = set()
            highlights = []
            
            for score, is_rc, found, alignment, fwdrev in [
                    (score_fwd, fwd_is_rc, fwd_found, fwd_alignment, 0),
                    (score_rev, rev_is_rc, rev_found, rev_alignment, 1)]:
                if float(score) > 2.5:
                    for sample_name in sample_barcodes:
                        if found == sample_barcodes[sample_name][fwdrev]:
                            break
                    else:
                        continue

                    sample_names.add(sample_name)
                    highlights.append(json.loads(alignment))

            if sample_names:
                original = fname.replace(".barcoding", "")
                with gzip.open(original, mode='rt') as inf2:
                    for seqno, (title, sequence, quality) in enumerate(
                            FastqGeneralIterator(inf2)):
                        if seqno == lineno:
                            sequence = list(sequence)
                            for _, highlight in highlights:
                                for begin, end in highlight:
                                    for i in range(begin, end):
                                        sequence[i] = sequence[i].lower()

                            row = ">%s %s\n%s" % (
                                title,
                                " ".join(sorted(sample_names)),
                                "".join(sequence))

                            # found both ends and they match the same sample
                            if len(highlights) == 2 and len(sample_name) == 1:
                                out_strict.append(row)
                            else:
                                out_generous.append(row)


with open("demultiplexed-strict.fasta", "w") as outf:
    for row in out_strict:
        outf.write(row)
        outf.write("\n")

with open("demultiplexed-generous.fasta", "w") as outf:
    for row in out_generous:
        outf.write(row)
        outf.write("\n")


