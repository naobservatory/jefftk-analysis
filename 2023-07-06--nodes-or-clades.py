#!/usr/bin/env python3

# Answer the question: do any of the taxids listed in human-viruses.tsv have
# children that are not listed in human-viruses.tsv?  If so we need to be
# looking for these children too, which we currently don't do.

hv = {}
with open("/Users/jeffkaufman/code/mgs-pipeline/human-viruses.tsv") as inf:
    for line in inf:
        taxid, name = line.strip().split("\t")
        hv[int(taxid)] = name

children = {}
with open("/Users/jeffkaufman/code/mgs-pipeline/dashboard/nodes.dmp") as inf:
    for line in inf:
        child_taxid, parent_taxid, rank, *_ = \
            line.replace("\t|\n", "").split("\t|\t")
        child_taxid = int(child_taxid)
        parent_taxid = int(parent_taxid)

        if parent_taxid not in children:
            children[parent_taxid] = []
        children[parent_taxid].append(child_taxid)

def validate(top_taxid, top_name, taxid):
    for child in children.get(taxid, []):
        if child not in hv:
            print("For %s (%s) we're missing %s under %s (%s)" % (
                top_name, top_taxid, child, hv[taxid], taxid))
        else:
            validate(top_taxid, top_name, child)
        
for taxid, name in sorted(hv.items()):
    validate(taxid, name, taxid)
