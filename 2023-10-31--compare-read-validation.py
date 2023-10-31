#!/usr/bin/env python3

# Assumes the following have already run:
#   mgs-pipeline/dashboard/prepare-dashboard-data.sh

import os
import re
import sys
import gzip
import json
from collections import defaultdict, Counter

papers = [
  "Brinch 2020",
  "Rothman 2021",
  "Spurbeck 2023",
  "Crits-Christoph 2021",
]

dashboard_dir="/Users/jeffkaufman/code/mgs-pipeline/dashboard/"

with open(os.path.join(dashboard_dir, "metadata_papers.json")) as inf:
  metadata_papers = json.load(inf)

with open(os.path.join(dashboard_dir, "metadata_bioprojects.json")) as inf:
  metadata_bioprojects = json.load(inf)

with open(os.path.join(dashboard_dir, "metadata_samples.json")) as inf:
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

children = defaultdict(set)
with open(os.path.join(dashboard_dir, "nodes.dmp")) as inf:
  for line in inf:
    child_taxid, parent_taxid, rank, *_ = \
      line.replace("\t|\n", "").split("\t|\t")
    child_taxid = int(child_taxid)
    parent_taxid = int(parent_taxid)
    children[parent_taxid].add(child_taxid)

detailed_taxids = defaultdict(set) # target taxid -> all taxids in clade

def populate_detailed_taxids(root_taxid, taxid):
  detailed_taxids[root_taxid].add(taxid)
  for child in children[taxid]:
    populate_detailed_taxids(root_taxid, child)

broad_taxid_from_detailed = {}
for taxid in taxid_to_name:
  populate_detailed_taxids(taxid, taxid)
  for broad_taxid in detailed_taxids[taxid]:
    broad_taxid_from_detailed[broad_taxid] = taxid

simon_details = {}
with open("simon-excluded-read-ids.txt") as inf:
  for line in inf:
    read_id = line.strip()
    sample = line.split(".")[0].removeprefix("M_")

    with open(os.path.join(
        dashboard_dir, "hvreads", "%s.hvreads.json" % sample)) as hvreads:
      hvrec = json.load(hvreads)[read_id]
      detailed_taxid = hvrec[0]
      if detailed_taxid not in broad_taxid_from_detailed:
        continue

      simon_details[read_id] = hvrec

jeff_details = {}
with open("jeff-excluded-read-ids.txt") as inf:
  for line in inf:
    read_id = line.strip()
    sample = line.split(".")[0].removeprefix("M_")

    with open(os.path.join(
        dashboard_dir, "hvreads", "%s.hvreads.json" % sample)) as hvreads:
      hvrec = json.load(hvreads)[read_id]
      detailed_taxid = hvrec[0]
      if detailed_taxid not in broad_taxid_from_detailed:
        assert False

      jeff_details[read_id] = hvrec

for read_id in jeff_details:
  if read_id in simon_details: continue
  assert False
      
for broad_taxid in taxid_to_name:
  print(broad_taxid, taxid_to_name[broad_taxid])

  for read_id in simon_details:
    if read_id in jeff_details: continue

    hvrec = simon_details[read_id]
    detailed_taxid = hvrec[0]
    if broad_taxid_from_detailed[detailed_taxid] != broad_taxid:
      continue

    reads = hvrec[2:]
    
    print("  S", read_id)
    for read, quality in reads:
      print("    ", read)


    
      
  
      
      

