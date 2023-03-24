#!/usr/bin/env python3

# wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
# unzip taxdmp.zip

# taxid -> Node
nodes = {}

class Node:
    def __init__(self, taxid, name):
        self.taxid = taxid
        self.name = name
        self.parent = None
        self.rank = None
        self.children = []

with open("names.dmp") as inf:
    for line in inf:
        taxid, name, unique_name, name_class = line.replace(
            "\t|\n", "").split("\t|\t")
        taxid = int(taxid)
        nodes[taxid] = Node(taxid, name)
        
with open("nodes.dmp") as inf:
    for line in inf:
        taxid, parent_taxid, rank, *_ = \
            line.replace("\t|\n", "").split("\t|\t")
        taxid = int(taxid)
        parent_taxid = int(parent_taxid)
        nodes[taxid].parent = nodes[parent_taxid]
        nodes[taxid].rank = rank
        nodes[taxid].parent.children.append(nodes[taxid])

root = nodes[1]
        
def print_tree(node, indent=0):
    print("  "*indent + node.name)
    for child in node.children:
        if child is node: continue
        print_tree(child, indent = indent+1)

print_tree(root)
