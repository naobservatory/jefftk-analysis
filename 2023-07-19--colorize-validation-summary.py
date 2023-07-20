#!/usr/bin/env python3

import sys
import gzip
import json
from collections import defaultdict

COLOR_END = '\x1b[0m'
COLOR_RED = '\x1b[1;31m'
COLOR_BLUE = '\x1b[1;34m'
COLOR_MAGENTA = '\x1b[1;35m'

# For now just take the first reference genome for each

metadata_fname = "metadata.txt"
our_taxid_to_all_fname = "our_taxid_to_all.json"
taxid_to_fastas = defaultdict(list)
with open(metadata_fname) as inf:
  cols = None
  for line in inf:
    rows = line.strip().split("\t")
    if not cols:
      cols = rows
      continue

    taxid_to_fastas[int(rows[cols.index("taxid")])].append(
      rows[cols.index("local_filename")])

with open(our_taxid_to_all_fname) as inf:
  our_taxid_to_all = json.load(inf)

all_to_our_taxid = {}
for our_taxid, all_taxids in our_taxid_to_all.items():
  for all_taxid in all_taxids:
    all_to_our_taxid[all_taxid] = int(our_taxid)

def pick_genome(our_taxid):
  for any_taxid in sorted(our_taxid_to_all[str(our_taxid)]):
    if any_taxid in taxid_to_fastas:
      for fasta in sorted(taxid_to_fastas[any_taxid]):
        return fasta
  assert False

def load_genome(our_taxid):
  genome = []

  fasta_fname = pick_genome(our_taxid)
  o = gzip.open if fasta_fname.endswith(".gz") else open
  with o(fasta_fname) as inf:
    for line in inf:
      if type(line) == type(b""):
        line = line.decode("utf-8")
      line = line.strip()
      if line.startswith(">"): continue
      genome.append(line)
  return "".join(genome)

def rc(s):
  return "".join({'T':'A',
                  'G':'C',
                  'A':'T',
                  'C':'G',
                  'N':'N'}[x] for x in reversed(s))

validation_summary, = sys.argv[1:]
taxid = int(validation_summary.removesuffix(".validation-summary.tsv"))
genome = load_genome(taxid)

def color_alignment(read, canonical):
  out = []
  for c1, c2 in zip(read, canonical):
    if c1 != c2:
      out.append(COLOR_RED)
    out.append(c1)
    if c1 != c2:
      out.append(COLOR_END)
  return "".join(out)

with open(validation_summary) as inf:
  for line in inf:
    paper, read_id, read1_desc, read2_desc, read1, read2 = \
      line.removesuffix("\n").split("\t")

    pos1 = int(read1_desc.split("@")[-1])
    pos2 = int(read2_desc.split("@")[-1])
    
    print(
      COLOR_BLUE + paper + COLOR_END,
      COLOR_MAGENTA + read_id + COLOR_END,
      read1_desc,
      read2_desc,
      color_alignment(read1, genome[pos1:]),
      color_alignment(read2, genome[pos2:]))
