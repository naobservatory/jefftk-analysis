#!/usr/bin/env python3

import re
from Bio.SeqIO.FastaIO import SimpleFastaParser
from collections import defaultdict

def canonical_title(title):
    return re.sub(" basecall_model_version_id.*", "", title)

seqs = defaultdict(list)

with open("demultiplexed-strict.fasta") as inf:
    for title, sequence in SimpleFastaParser(inf):
        seqs[canonical_title(title)].append((title, sequence))
        
with open("demultiplexed-generous.fasta") as inf:
    for title, sequence in SimpleFastaParser(inf):
        seqs[canonical_title(title)].append((title, sequence))

with open("probes-and-primers.fasta") as inf:
    for title, sequence in SimpleFastaParser(inf):
        seqs[canonical_title(title)].append((title, sequence))

out = []
for title, records in seqs.items():
    if len(records) == 1: continue

    seq = None
    new_title = None
    for full_title, sequence in records:
        if not seq:
            seq = sequence
            new_title = full_title
        else:
            new_seq = []
            for c1, c2 in zip(seq, sequence):
                if c1.lower() != c2.lower():
                    raise Exception("Mismatch: %s %s %s" % (c1, c2, title))
                new_seq.append(max(c1, c2))  # "a" > "A"
            seq = "".join(new_seq)

            title_suffix, = re.findall("basecall_model_version_id=[^ ]*(.*)$",
                                       full_title)
            new_title += " " + title_suffix

    out.append((new_title, seq))
                
with open("combined.fasta", "w") as outf:
    for title, seq in out:
        outf.write(">%s\n%s\n" % (title, seq))
