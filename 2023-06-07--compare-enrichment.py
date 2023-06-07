#!/usr/bin/env python3

import json
import glob
import gzip

# https://www.illumina.com/products/by-type/sequencing-kits/library-prep-kits/respiratory-virus-oligo-panel.html
listed_targets = [
    "Human coronavirus 229E",
    "Human coronavirus NL63",
    "Human coronavirus OC43",
    "Human coronavirus HKU1",
    "SARS-CoV-2",
    "Human adenovirus B1",
    "Human adenovirus C2",
    "Human adenovirus E4",
    "Human bocavirus 1 (Primate bocaparvovirus 1 isolate st2)",
    "Human bocavirus 2c PK isolate PK-5510",
    "Human bocavirus 3",
    "Human parainfluenza virus 1",
    "Human parainfluenza virus 2",
    "Human parainfluenza virus 3",
    "Human parainfluenza virus 4a",
    "Human metapneumovirus (CAN97-83)",
    "Respiratory syncytial virus (type A)",
    "Human Respiratory syncytial virus 9320 (type B)",
    "Influenza A virus (A/Puerto Rico/8/1934(H1N1))",
    "Influenza A virus (A/Korea/426/1968(H2N2))",
    "Influenza A virus (A/New York/392/2004(H3N2))",
    "Influenza A virus (A/goose/Guangdong/1/1996(H5N1))",
    "Human bocavirus 4 NI strain HBoV4-NI-385",
    "KI polyomavirus Stockholm 60",
    "WU Polyomavirus",
    "Human parechovirus type 1 PicoBank/HPeV1/a",
    "Human parechovirus 6",
    "Human rhinovirus A89",
    "Human rhinovirus C (strain 024)",
    "Human rhinovirus B14",
    "Human enterovirus C104 strain: AK11",
    "Human enterovirus C109 isolate NICA08-4327",
    "Influenza A virus (A/Zhejiang/DTID-ZJU01/2013(H7N9))",
    "Influenza A virus (A/Hong Kong/1073/99(H9N2))",
    "Influenza A virus (A/Texas/50/2012(H3N2))",
    "Influenza A virus (A/Michigan/45/2015(H1N1))",
    "Influenza B virus (B/Lee/1940)",
    "Influenza B virus (B/Wisconsin/01/2010)",
    "Influenza B virus (B/Brisbane/60/2008)",
    "Influenza B virus (B/Colorado/06/2017)",
    "Influenza B virus (B/Washington/02/2019)",
]

other_targets = [
    "Pepper mild mottle virus",
    "Influenza A virus",
    "Influenza B virus",
    "Bacteria",
    "Viruses",
    "Norovirus",    
]

targets = listed_targets + other_targets

with open("/Users/jeffkaufman/code/mgs-pipeline/dashboard/"
          "taxonomic_names.json") as inf:
    taxonomic_names = json.load(inf)

target_taxids = {}
    
for target_name in targets:
    for taxid, names in taxonomic_names.items():
        if target_name in names:
            target_taxids[target_name] = taxid

with open("/Users/jeffkaufman/code/mgs-pipeline/dashboard/"
          "metadata_samples.json") as inf:
    sample_metadata = json.load(inf)

# taxid -> {"viral": unenriched clade assignments,
#           "panel": enriched clade assignments}
taxid_clade_assignments = {}
for taxid in target_taxids.values():
    taxid_clade_assignments[taxid] = {"viral": 0, "panel": 0}

# date, wtp -> {"viral": sample, "panel": sample}
sample_pairs = {}
for cladecounts in glob.glob("PRJNA729801-cladecounts/*.tsv.gz"):
    sample = cladecounts.split("/")[-1].split(".")[0]

    key = (sample_metadata[sample]["date"],
           sample_metadata[sample]["fine_location"])
    if key not in sample_pairs:
        sample_pairs[key] = {"viral": None, "panel": None}

    sample_pairs[key][sample_metadata[sample]["enrichment"]] = sample

reads = {"viral": 0,
         "panel": 0}
for sample_pair in sample_pairs.values():
    if not sample_pair["viral"] or not sample_pair["panel"]: continue
    for enrichment, sample in sample_pair.items():
        reads[enrichment] += sample_metadata[sample]["reads"]

        with gzip.open("PRJNA729801-cladecounts/%s.tsv.gz" % sample) as inf:
            for line in inf:
                taxid, _, _, clade_assigments, _ = line.decode(
                    "utf-8").strip().split("\t")
                if taxid in taxid_clade_assignments:
                    taxid_clade_assignments[
                        taxid][enrichment] += int(clade_assigments)

print("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
    "target",
    "taxid",
    "panel",
    "viral",
    "ratio",
    "panel ra",
    "viral ra",
))
for target_name, target_taxid in sorted(target_taxids.items()):
    viral = taxid_clade_assignments[target_taxid]["viral"]
    panel = taxid_clade_assignments[target_taxid]["panel"]
    
    if viral and panel:
        raw_ratio = panel / viral * reads["viral"] / reads["panel"]
        if raw_ratio > 10:
            ratio = "%.0f" % raw_ratio
        if raw_ratio > 1:
            ratio = "%.1f" % raw_ratio
        else:
            ratio = "%.2f" % raw_ratio
    elif panel:
        ratio = "<%s" % round(panel * reads["viral"] / reads["panel"])
    elif viral:
        ratio = "0"
    else:
        ratio = "n/a"

    print("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (
        target_name,
        target_taxid,
        panel,
        viral,
        ratio,
        panel / reads["panel"],
        viral / reads["viral"],
    ))
                
