#!/usr/bin/env python3

import sys
from collections import defaultdict

target_taxid, = sys.argv[1:]

children = defaultdict(list)  # parent_taxid -> [child_taxid]
with open("/Users/jeffkaufman/code/mgs-pipeline/dashboard/nodes.dmp") as inf:
    for line in inf:
        child_taxid, parent_taxid, rank, *_ = \
            line.replace("\t|\n", "").split("\t|\t")
        child_taxid = child_taxid
        parent_taxid = parent_taxid
        children[parent_taxid].append(child_taxid)

all_taxids = set()
def collect(taxid, depth=0):
    #print(" "* depth, taxid)
    all_taxids.add(taxid)
    for child in children[taxid]:
        if child not in all_taxids:
            collect(child, depth=depth+1)
collect(target_taxid)
print(",".join(sorted(all_taxids)))
