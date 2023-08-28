#!/usr/bin/env python3

import sys
import json

MGS_PIPELINE_DIR="../mgs-pipeline"

# taxid -> parent
parents = {}
# taxid -> rank
ranks = {}
with open("%s/dashboard/nodes.dmp" % MGS_PIPELINE_DIR) as inf:
    for line in inf:
        child_taxid, parent_taxid, rank, *_ = \
            line.replace("\t|\n", "").split("\t|\t")
        parents[int(child_taxid)] = int(parent_taxid)
        ranks[int(child_taxid)] = rank

TAXID_ROOT=1
TAXID_BACTERIA=2
TAXID_VIRUSES=10239

bacteria = set()
viruses = set()

for taxid in parents:
    tmp = taxid
    while tmp != TAXID_ROOT:
        tmp = parents[tmp]
        if tmp == TAXID_BACTERIA:
            bacteria.add(taxid)
            break
        elif tmp == TAXID_VIRUSES:
            viruses.add(taxid)
            break

n_viral = 0
both = []
for line in sys.stdin:
    try:
        _, seq_id, _, _, kraken_info = line.strip().split("\t")
    except Exception:
        print(line)
        raise

    viral_count = 0
    bacterial_count = 0
    for segment in kraken_info.split():
        hit, count = segment.split(":")
        if hit == "A":
            continue # ambiguous nucleotide
        
        if hit == ":":
            continue # fwd/rev division
    
        hit = int(hit)
        count = int(count)

        if hit in viruses:
            viral_count += count
        elif hit in bacteria:
            bacterial_count += count

    if viral_count >= 10:
        n_viral += 1

        if bacterial_count >= 10:
            both.append((seq_id, kraken_info))

json.dump({
    "n_both": len(both),
    "n_viral": n_viral,
    "both": both,
    }, sys.stdout)

    
    
