#!/usr/bin/env python3

import os
import glob
import json
from collections import Counter
from collections import defaultdict

"""
Generate data for interactive HTML for exploring the taxonomy of human viruses
identified in the samples.

Columns are bioprojects, rows are taxonomic nodes, cells are relative abundances.

Rows can be expanded or collapsed to show child nodes.
"""

MGS_PIPELINE_DIR="/Users/jeffkaufman/code/mgs-pipeline"

human_viruses = set()   # {taxid}
with open("%s/human-viruses.tsv" % MGS_PIPELINE_DIR) as inf:
    for line in inf:
        taxid, name = line.strip().split("\t")
        taxid = int(taxid)
        human_viruses.add(taxid)

parents = {}  # child_taxid -> parent_taxid
with open("nodes.dmp") as inf:
    for line in inf:
        child_taxid, parent_taxid, rank, *_ = \
            line.replace("\t|\n", "").split("\t|\t")
        child_taxid = int(child_taxid)
        parent_taxid = int(parent_taxid)
        parents[child_taxid] = parent_taxid
        
# taxid -> node
nodes = {}

mentioned_taxids = set()
for virus_taxid in human_viruses:
    taxid = virus_taxid
    mentioned_taxids.add(taxid)
    while taxid != 1:
        taxid = parents[taxid]
        mentioned_taxids.add(taxid)

for taxid in mentioned_taxids:
    nodes[taxid] = [taxid]

for virus_taxid in human_viruses:
    taxid = virus_taxid
    node = nodes[taxid]
    while taxid != 1:
        taxid = parents[taxid]
        if node in nodes[taxid]:
            break
        nodes[taxid].append(node)
        node = nodes[taxid]

# Format: [taxid, children]
# Ex: [1, [10239, [10472,
#       [10473, [46014],
#               [10474, [10475, ...],
#                       [693626, ...],
#                       [1299307, ...]]], ...], ...], ...]
tree = nodes[1]
        
names = {}  # taxid -> name
with open("names.dmp") as inf:
    for line in inf:
        taxid, name, unique_name, name_class = line.replace(
            "\t|\n", "").split("\t|\t")
        taxid = int(taxid)
        if taxid in mentioned_taxids:
            names[taxid] = name
        
# project -> accession -> n_reads
project_sample_reads = defaultdict(dict)
for n_reads_fname in glob.glob(
        "%s/bioprojects/*/metadata/*.n_reads" % MGS_PIPELINE_DIR):
    project = n_reads_fname.split("/")[-3]
    accession = n_reads_fname.split("/")[-1].replace(".n_reads", "")
    if os.path.exists("%s.humanviruses.tsv" % accession):
        with open(n_reads_fname) as inf:
            project_sample_reads[project][accession] = int(inf.read())
        
projects = list(sorted(project_sample_reads))

bioproject_names = {} # project -> bioproject name
for project in projects:
    with open("%s/bioprojects/%s/metadata/name.txt" % (
            MGS_PIPELINE_DIR, project)) as inf:
        bioproject_names[project] = inf.read().strip()

# virus -> project -> [count, relab]
virus_project_counts = {}
for virus_taxid in human_viruses:
    virus_project_counts[virus_taxid] = {}
    for project in projects:
        project_count = 0
        project_total = 0
        for accession in project_sample_reads[project]:
            project_total += project_sample_reads[project][accession]
            hv = "%s.humanviruses.tsv" % accession
            with open(hv) as inf:
                for line in inf:
                    taxid, count, name = line.strip().split("\t")
                    taxid = int(taxid)
                    count = int(count)
                    if taxid == virus_taxid:
                        project_count += count
        if project_count > 0:
            virus_project_counts[virus_taxid][project] = (
                project_count, project_count / project_total)

with open("data.js", "w") as outf:
    for name, val in [
            ("virus_project_counts", virus_project_counts),
            ("projects", projects),
            ("names", names),
            ("bioproject_names", bioproject_names),
            ("tree", tree)]:
        outf.write("%s=%s;\n" % (name, json.dumps(
            val, indent=None if name == "tree" else 2)))

               
