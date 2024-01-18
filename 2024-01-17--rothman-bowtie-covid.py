#!/usr/bin/env python3
import os
import re
import sys
import json
from collections import defaultdict, Counter

sample, = sys.argv[1:]

COVID_VIRAL_ANCESTORY = [
    2697049, # Severe acute respiratory syndrome coronavirus 2
    694009,  # Severe acute respiratory syndrome-related coronavirus

    # No hits on any of the below, but keep them anyway

    2509511, # Sarbecovirus
    694002,  # Betacoronavirus
    2501931, # Orthocoronavirinae
    11118,   # Coronaviridae
    76804,   # Nidovirales
    2732506, # Pisoniviricetes
    2732408, # Pisuviricota
    2732396, # Orthornavirae
    2559587, # Riboviria
    10239,   # Viruses
]

cdata = {
    "covid-bowtie": defaultdict(list),
    "covid-kraken": set(),
    "non-covid-kraken": set(),
    "human": set(),
}
cdata["processed"] = {}

with open("%s.covid.alignments.tsv" % sample) as inf:
    for line in inf:
        read_id, genome, taxid, cigar, start_pos, part_score, part_length = \
            line.rstrip("\n").split("\t")
        cdata["covid-bowtie"][read_id].append((
            int(part_score),
            int(part_length)))

human_alignments_fname = "%s.covid-human.alignments.tsv" % sample
if os.path.exists(human_alignments_fname):
    with open(human_alignments_fname) as inf:
        for line in inf:
            read_id, *_ = line.split("\t")
            cdata["human"].add(read_id)

kraken_assignments_fname = "%s.kraken.tsv" % sample
if os.path.exists(kraken_assignments_fname):
    with open(kraken_assignments_fname) as inf:
        for line in inf:
            _, read_id, assignment, *_ = line.rstrip("\n").split("\t")
            taxid, = re.findall(".*taxid ([0-9]+).*", assignment)
            taxid = int(taxid)

            if taxid in COVID_VIRAL_ANCESTORY:
                if taxid == 2697049:
                    cdata["covid-kraken"].add(read_id)
            elif taxid != 0:
                cdata["non-covid-kraken"].add(read_id)

def interpret_score(score, length):
    bend1_length = 40
    bend1_score = 60

    bend2_length = 80
    bend2_score = 95

    max_length = 200
    max_score = 100

    min_length_cutoff = 28
    uncollapsed_score = 269

    if length < min_length_cutoff:
        return False
    elif not collapsed:
        return score >= uncollapsed_score
    elif length < bend1_length:
        return score >= bend1_score
    elif length < bend2_length:
        return score >= (
            (length - bend1_length) / (bend2_length - bend1_length)
            * (bend2_score - bend1_score)) + bend1_score
    else:
        return score >= (
            (length - bend2_length) / (max_length - bend2_length)
            * (max_score - bend2_score)) + bend2_score

r = Counter()

for read_id in sorted(cdata["covid-bowtie"]):
    if read_id in cdata["covid-kraken"]:
        r["n_kraken_decided_covid"] += 1
    
    collapsed = read_id.startswith("M_")

    score = 0
    length = 0
    for part_score, part_length in cdata["covid-bowtie"][read_id]:
        score += part_score
        length += part_length

    score_decision = interpret_score(score, length)

    decision = score_decision
    if score_decision:
        if read_id in cdata["human"]:
            r["n_rejected_for_human_bowtie"] += 1
            decision = False

        if read_id in cdata["non-covid-kraken"]:
            r["n_rejected_for_non_covid_kraken"] += 1
            decision = False
    
    if decision:
        r["n_decided_covid"] += 1
        if read_id not in cdata["covid-kraken"]:
            r["n_covid_missed_by_kraken"] += 1
    else:
        if read_id in cdata["covid-kraken"]:    
            r["n_covid_rejected_by_bowtie"] += 1

with open("%s.stats.json" % sample, "w") as outf:
    json.dump(r, outf)
