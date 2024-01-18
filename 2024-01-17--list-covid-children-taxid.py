#!/usr/bin/env python3

from collections import defaultdict

COVID=2697049
children = defaultdict(set)
with open("/home/ec2-user/mgs-pipeline/dashboard/nodes.dmp") as inf:
    for line in inf:
        child_taxid, parent_taxid, rank, *_ = \
            line.replace("\t|\n", "").split("\t|\t")
        child_taxid = int(child_taxid)
        parent_taxid = int(parent_taxid)
        children[parent_taxid].add(child_taxid)

# prints empty set, telling us the taxonomy doesn't
# currently recognize any childen
print(children[COVID])

        
