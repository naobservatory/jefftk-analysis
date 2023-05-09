#!/usr/bin/env python3

import sys

parents = {}  # child_taxid -> parent_taxid
with open("nodes.dmp") as inf:
    for line in inf:
        child_taxid, parent_taxid, rank, *_ = \
            line.replace("\t|\n", "").split("\t|\t")
        child_taxid = int(child_taxid)
        parent_taxid = int(parent_taxid)
        parents[child_taxid] = parent_taxid

for line in sys.stdin:
    line = line.strip()
    if line:
        taxid = int(line)

        taxids = []
        while parents[taxid] != taxid:
            taxids.append(taxid)
            taxid = parents[taxid]

        print("\t".join(str(x) for x in taxids))
            
        
    
