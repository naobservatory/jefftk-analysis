#!/usr/bin/env python3

import sys
import math
import matplotlib as mpl
from collections import defaultdict, Counter

VIRUS=10239
BACTERIA=2

COLOR_RED = '\x1b[1;31m'
COLOR_YELLOW = '\x1b[1;33m'
COLOR_CYAN = '\x1b[1;36m'
COLOR_END = '\x1b[0m'

in1, = sys.argv[1:]

human_viruses = set()
with open("human-viruses.tsv") as inf:
    for line in inf:
        human_viruses.add(int(line.strip().split()[0]))

children = defaultdict(list)  # parent_taxid -> [children]
with open("dashboard/nodes.dmp") as inf:
    for line in inf:
        child_taxid, parent_taxid, rank, *_ = \
            line.replace("\t|\n", "").split("\t|\t")
        child_taxid = int(child_taxid)
        parent_taxid = int(parent_taxid)
        if child_taxid != parent_taxid:
            children[parent_taxid].append(child_taxid)

def populate_set(taxid, s):
    s.add(taxid)
    for child in children[taxid]:
        populate_set(child, s)

viruses = set()
populate_set(VIRUS, viruses)

bacteria = set()
populate_set(BACTERIA, bacteria)

viruses = set()
def populate_set(taxid):
    viruses.add(taxid)
    for child in children[taxid]:
        populate_set(child)
populate_set(VIRUS)

names = defaultdict(str)
with open("dashboard/names.dmp") as inf:
    for line in inf:
        taxid, name, _, category, *_ = \
            line.replace("\t|\n", "").split("\t|\t")
        taxid = int(taxid)
        if not names[taxid] or category == "scientific name":
            names[taxid] = name

nreads = 0
nviruses = None
nbacteria = None
to_print = [] # key, val
with open(in1) as inf:
    for line in inf:
        taxid, direct_assignments, direct_hits, \
            clade_assignments, clade_hits = line.strip().split("\t")
        taxid = int(taxid)
        direct_assignments = int(direct_assignments)
        clade_assignments = int(clade_assignments)
        if taxid in [0, 1]:
            nreads += clade_assignments
        elif taxid == VIRUS:
            nviruses = clade_assignments
        elif taxid == BACTERIA:
            nbacteria = clade_assignments

        if taxid in human_viruses:
            color = COLOR_RED
        elif taxid in viruses:
            color = COLOR_YELLOW
        elif taxid not in bacteria:
            color = COLOR_CYAN
        else:
            color = ""

        to_print.append((direct_assignments,
                         "%s %s%s%s (%s)" % (
                             str(direct_assignments).rjust(8),
                             color,
                             names[taxid],
                             COLOR_END,
                             taxid)))

print("Viruses:  %s (%.2f%%)" % (nviruses, 100*nviruses/nreads))
print("Bacteria: %s (%.2f%%)" % (nbacteria, 100*nbacteria/nreads))
print()

for _, line in sorted(to_print, reverse=True):
    print(line)
