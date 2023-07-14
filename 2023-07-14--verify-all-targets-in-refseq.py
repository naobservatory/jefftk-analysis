#!/usr/bin/env python3
import os
import re
import json
import subprocess
from collections import defaultdict

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

taxids_fname = "taxids.txt"
our_taxid_to_all_fname = "our_taxid_to_all.json"
our_taxid_to_all = defaultdict(list)
if not os.path.exists(taxids_fname):
    fetch = []
    for taxid in taxid_to_name:
        fetch.append(taxid)
        our_taxid_to_all[taxid].append(taxid)

        output = subprocess.check_output(["gimme_taxa.py", str(taxid)])
        output = output.decode("utf-8")

        for line in output.split("\n"):
            line = line.strip()
            if line.startswith("parent_taxid") or not line:
                continue
            _, descendent_taxid, descendent_name = line.split("\t")
            descendent_taxid = int(descendent_taxid)
            fetch.append(descendent_taxid)
            our_taxid_to_all[taxid].append(descendent_taxid)
    with open(taxids_fname, "w") as outf:
        for taxid in fetch:
            outf.write("%s\n" % taxid)

    with open(our_taxid_to_all_fname, "w") as outf:
        json.dump(our_taxid_to_all, outf)

metadata_fname = "metadata.txt"
if not os.path.exists(metadata_fname):
    subprocess.check_call([
        "ncbi-genome-download",
        "--taxids", taxids_fname,
        "--formats", "fasta",
        "--metadata-table", metadata_fname,
        "viral",
    ])

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
        
for our_taxid in taxid_to_name:
    seen = False
    for any_taxid in our_taxid_to_all[str(our_taxid)]:
        if any_taxid in taxid_to_fastas:
            seen = True
    if not seen:
        print("Missing %s (%s)" % (our_taxid, taxid_to_name[our_taxid]))

        
        
    
