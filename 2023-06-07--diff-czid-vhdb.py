#!/usr/bin/env python3

# Viruses from https://czid.org/pathogen_list
# "It serves to highlight organisms with known pathogenicity in humans, but may
# not be fully comprehensive."
with open("czid.txt") as inf:
    czid_taxids = set()
    czid_names = {}
    for n, line in enumerate(inf):
        line = line.strip()
        if n % 2 == 0:
            name = line
        else:
            assert line.startswith("Tax ID: ")
            taxid = int(line.replace("Tax ID: ", ""))
            
            czid_taxids.add(taxid)
            czid_names[taxid] = name

# From Virus-Host DB: https://www.genome.jp/virushostdb/
# Taken from https://github.com/naobservatory/mgs-pipeline/commit/df7588fcad49c99b69b27589489229e3d890da38
with open("vhdb.txt") as inf:
    vhdb_taxids = set(int(x) for x in inf)

# taxid -> name
names = {}
with open("/Users/jeffkaufman/code/mgs-pipeline/dashboard/names.dmp") as inf:
    for line in inf:
        taxid, name, unique_name, name_class = line.replace(
            "\t|\n", "").split("\t|\t")
        taxid = int(taxid)
        if taxid not in names or name_class == "scientific name":
            names[taxid] = name

if False:
    only_czid = czid_taxids - vhdb_taxids
    if only_czid:
        print("Only in CZID:")
        for taxid in sorted(only_czid):
            print("  %s\t%s" % (taxid, names.get(taxid, czid_names[taxid])))

    only_vhdb = vhdb_taxids - czid_taxids
    if only_vhdb:
        print("Only in Virus-Host DB:")
        for taxid in sorted(only_vhdb):
            print("  %s\t%s" % (taxid, names[taxid]))
            
# The problem is, sometimes one of them picks a higher clade than the other.
# We don't want to say it's a real disagreement when CZID picks
# 28876 (Rotavirus B) and VHDB picks 10942 (Human rotavirus B).
#
# Let's call something a disagreement if one names a taxid and the other
# doesn't name that taxid or anything in its clade (and the reverse)

parents = {}  # child_taxid -> parent_taxid
with open("/Users/jeffkaufman/code/mgs-pipeline/dashboard/nodes.dmp") as inf:
    for line in inf:
        child_taxid, parent_taxid, rank, *_ = \
            line.replace("\t|\n", "").split("\t|\t")
        child_taxid = int(child_taxid)
        parent_taxid = int(parent_taxid)
        parents[child_taxid] = parent_taxid

def is_ancestor(taxid_child, taxid_parent):
    while taxid_child != 1:
        if taxid_child == taxid_parent:
            return True
        taxid_child = parents[taxid_child]
    return False

def has_any_ancestors(target, other):
    for other_taxid in other:
        if is_ancestor(target, other_taxid):
            return True
    return False

def has_any_descendants(target, other):
    for other_taxid in other:
        if is_ancestor(other_taxid, target):
            return True
    return False

def disagreements(target, other):
    for taxid in target:
        try:
            if has_any_ancestors(taxid, other): continue
            if has_any_descendants(taxid, other): continue
        except KeyError:
            yield taxid

        yield taxid

only_czid = disagreements(czid_taxids, vhdb_taxids)
only_vhdb = disagreements(vhdb_taxids, czid_taxids)

if only_czid:
    print("Only in CZID:")
    for taxid in sorted(only_czid):
        print("  %s\t%s" % (taxid, names.get(taxid, czid_names[taxid])))

if only_vhdb:
    print("Only in Virus-Host DB:")
    for taxid in sorted(only_vhdb):
        print("  %s\t%s" % (taxid, names[taxid]))
