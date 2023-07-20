#!/usr/bin/env python3

# Assumes the following have already run:
#   mgs-pipeline/dashboard/prepare-dashboard-data.sh
#   jefftk-analysis/2023-07-14--verify-all-targets-in-refseq.py

import os
import re
import sys
import gzip
import json
import argparse
from collections import defaultdict, Counter


COLOR_END = '\x1b[0m'
COLOR_RED = '\x1b[1;31m'
COLOR_BLUE = '\x1b[1;34m'
COLOR_MAGENTA = '\x1b[1;35m'

parser = argparse.ArgumentParser()
parser.add_argument("--paper")
parser.add_argument("--sample")
parser.add_argument("--genome-fasta")
parser.add_argument("--taxid", type=int)
parser.add_argument("dashboard_dir")
args = parser.parse_args()

papers = [
  "Brinch 2020",
  "Rothman 2021",
  "Spurbeck 2023",
  "Crits-Christoph 2021",
]

with open(os.path.join(args.dashboard_dir, "metadata_papers.json")) as inf:
  metadata_papers = json.load(inf)

with open(os.path.join(args.dashboard_dir, "metadata_bioprojects.json")) as inf:
  metadata_bioprojects = json.load(inf)

with open(os.path.join(args.dashboard_dir, "metadata_samples.json")) as inf:
  metadata_samples = json.load(inf)

# output of p2ra/list_taxids.py
taxid_list="""
10298           hsv_1           HSV-1
10310           hsv_2           HSV-2
10359           cmv             CMV
10376           ebv             EBV
10407           hbv             HBV
10566           hpv             HPV
10632           jcv             JCV
10804           aav2            AAV2
11103           hcv             HCV
11320           influenza       Influenza A
11520           influenza       Influenza B
11676           hiv             HIV
68558           aav6            AAV6
82300           aav5            AAV5
122928          norovirus       Norovirus (GI)
122929          norovirus       Norovirus (GII)
493803          mcv             MCV
1891762         bkv             BKV
2697049         sars_cov_2      SARS-COV-2
"""

taxid_to_name = {}
for line in taxid_list.split("\n"):
  line = line.strip()
  if not line: continue

  taxid, human_readable = re.match(
    "^([0-9]+) +[a-z0-9_]+ +(.*)$", line).groups()
  taxid_to_name[int(taxid)] = human_readable

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

def pick_genomes(our_taxid):
  if args.genome_fasta:
    yield args.genome_fasta
  for any_taxid in sorted(our_taxid_to_all[str(our_taxid)]):
    if any_taxid in taxid_to_fastas:
      for fasta in sorted(taxid_to_fastas[any_taxid]):
        yield fasta

def load_genomes(our_taxid):
  genomes = []

  for fasta_fname in pick_genomes(our_taxid):
    genome = []
    o = gzip.open if fasta_fname.endswith(".gz") else open
    with o(fasta_fname) as inf:
      for line in inf:
        if type(line) == type(b""):
          line = line.decode("utf-8")
        line = line.strip()
        if line.startswith(">"): continue
        genome.append(line)
    genomes.append("".join(genome))

  genomes.sort()
  return genomes

def rc(s):
  return "".join({'T':'A',
                  'G':'C',
                  'A':'T',
                  'C':'G',
                  'N':'N'}[x] for x in reversed(s))

DUP_LEN=25
def count_dups(hvr_fname):
  if os.path.exists(hvr_fname):
    with open(hvr_fname) as inf:
      hvr = json.load(inf)
  else:
    hvr = {}

  by_start_end = defaultdict(list) # start, end -> read id
  for read_id, (kraken_info, *reads) in sorted(hvr.items()):
    assert reads
    if len(reads) == 1:
      read, =reads
      if len(read) < DUP_LEN:
        continue
      start = read[DUP_LEN:]
      end = read[:-DUP_LEN]
    else:
      fwd, rev = reads
      if len(fwd) < DUP_LEN or len(rev) < DUP_LEN:
        continue
      start = fwd[DUP_LEN:]
      end = rev[DUP_LEN:]

    start_rc = rc(start)
    if start_rc < start:
      start = rc(end)
      end = start_rc

    by_start_end[start, end].append(read_id)

  read_ids_by_kraken_info = defaultdict(
    list) # kraken_info -> non-duplicate read ids
  for (start, end), read_ids in sorted(by_start_end.items()):
    first_kraken_info = hvr[read_ids[0]][0]
    read_ids_by_kraken_info[first_kraken_info].append(read_ids[0])
  return read_ids_by_kraken_info

def load_reads(read_id_records):
  loaded = {}
  papers_by_read_id = {}
  for paper, bioproject, sample, read_id in read_id_records:
    papers_by_read_id[read_id] = paper
    hvr_fname = os.path.join(
      args.dashboard_dir, "hvreads/%s.hvreads.json" % sample)
    with open(hvr_fname) as inf:
      loaded[read_id] = json.load(inf)[read_id]

  return loaded, papers_by_read_id

def compute_score(genome, pos, candidate):
  matches = 0
  for genome_base, candidate_base in zip(genome[pos:], candidate):
    if genome_base == candidate_base:
      matches += 1
  return matches

def score_candidate(genome, candidate):
  best_pos = 0
  best_score = 0
  for pos in range(len(genome) - len(candidate)):
    score = compute_score(genome, pos, candidate)
    if score > best_score:
      best_pos = pos
      best_score = score
  return best_score, best_pos

def align_read_to_one_genome(genome, read_info):
  if len(read_info) == 2:
    kraken_info, fwd = read_info
    rev = ""
  else:
    kraken_info, fwd, rev = read_info

  # For now, ignore reverse reads

  best_pos = 0
  best_candidate = None
  best_score = 0
  best_label = None
  for candidate, label in [
      (fwd, 'fwd'),
      (rc(fwd), 'rc_fwd'),
      (rev, 'rev'),
      (rc(rev), 'rc_rev')
  ]:
    score, pos = score_candidate(genome, candidate)
    if score > best_score:
      best_candidate = candidate
      best_pos = pos
      best_score = score
      best_label = label

  other, other_label = {
    'fwd': (rc(rev), "rc_rev"),
    'rc_rev': (fwd, "fwd"),
    'rc_fwd': (rev, "rev"),
    'rev': (rc(fwd), "rc_fwd"),
  }[best_label]

  if not other:
    assert fwd
    assert not rev

    return best_score, 0, best_pos, best_candidate, 0, "", best_label, None

  other_score, other_pos = score_candidate(genome, other)
  return best_score, other_score, best_pos, best_candidate, \
    other_pos, other, best_label, other_label

def align_read(genomes, read_info):
  best_score = 0
  best_result = None
  best_genome_index = None
  for genome_index, genome in enumerate(genomes):
    result = align_read_to_one_genome(genome, read_info)
    score1, score2, *_ = result
    score = score1 + score2
    if score1 + score2 > best_score:
      best_result = result
      best_score = score
      best_genome_index = genome_index
      
  return best_result, genome_index

def offset(pos, fwd):
  return " "*pos + fwd

def color_alignment(read, canonical):
  out = []
  for c1, c2 in zip(read, canonical):
    if c1 != c2:
      out.append(COLOR_RED)
    out.append(c1)
    if c1 != c2:
      out.append(COLOR_END)
  return "".join(out)

def validate(our_taxid):
  genomes = load_genomes(our_taxid)
  reads, papers_by_read_id = load_reads(
    determine_non_duplicate_read_ids(our_taxid))

  if not papers_by_read_id:
    return

  with open("%s.validation-summary.tsv" % our_taxid, "w") as outf:
    for read_id, read_info in sorted(reads.items()):
      if args.sample and not (
          read_id.startswith(args.sample) or
          read_id.removeprefix("M_").startswith(args.sample)):
        continue

      if args.paper and args.paper != papers_by_read_id[read_id]:
        continue

      (score1, score2,
       pos1, read1,
       pos2, read2,
       label1, label2), genome_index = align_read(genomes, read_info)

      outf.write("\t".join((
        papers_by_read_id[read_id],
        read_id,
        str(genome_index),
        "%s:%s/%s @%s" % (label1, score1, len(read1), pos1),
        "%s:%s/%s @%s" % (label2, score2, len(read2), pos2),
        read1, read2)) + "\n")
      outf.flush()

def determine_non_duplicate_read_ids(our_taxid):
  ndrids_fname = "%s.ndrids.tsv" % our_taxid
  if not os.path.exists(ndrids_fname):
    non_duplicate_read_ids = []
    for paper in sorted(papers):
      n = 0
      for bioproject in sorted(metadata_papers[paper]["projects"]):
        for sample in sorted(metadata_bioprojects[bioproject]):
          sample_details = metadata_samples[sample]
          if sample_details.get("enrichment", None) == "panel":
            continue

          allmatches_fname = os.path.join(
            args.dashboard_dir, "allmatches/%s.allmatches.tsv" % sample)
          hvr_fname = os.path.join(
            args.dashboard_dir, "hvreads/%s.hvreads.json" % sample)
          read_ids_by_kraken_info = count_dups(hvr_fname)

          with open(allmatches_fname) as inf:
            for line in inf:
              line = line.strip()
              if not line: continue

              _, _, name_and_taxid, _, kraken_info = line.split("\t")
              taxid = int(name_and_taxid.split()[-1].replace(")", ""))

              if taxid not in all_to_our_taxid:
                continue
              if all_to_our_taxid[taxid] != our_taxid:
                continue

              if not read_ids_by_kraken_info[kraken_info]:
                print("  skipping duplicate")
                continue
              non_duplicate_read_ids.append((
                paper, bioproject, sample,
                read_ids_by_kraken_info[kraken_info].pop()))
              n += 1
      print(" ", paper, n)

    with open(ndrids_fname, "w") as outf:
      for paper, bioproject, sample, read_id in sorted(non_duplicate_read_ids):
        outf.write("%s\t%s\t%s\t%s\n" % (
          paper, bioproject, sample, read_id))

  with open(ndrids_fname) as inf:
    return [
      line.strip().split("\t")
      for line in inf]

for our_taxid in taxid_to_name:
  if args.taxid and args.taxid != our_taxid: continue
  if taxid_to_name[our_taxid] == "AAV6":
    continue # Not in RefSeq  validate(our_taxid)

  validate(our_taxid)
