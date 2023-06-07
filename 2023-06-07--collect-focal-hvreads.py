#!/usr/bin/env python3
import json
import sys
import glob
import subprocess

focal_pathogens = {
    "sars_cov_2": 2697049,
    "hiv": 11676,
    "influenza_a": 11320,
    "influenza_b": 11520,
    "norovirus_all": 142786,
    "norovirus_g1": 122928,
    "norovirus_g2": 122929,
    "hcv": 11103,
}

bioprojects = {
    "crits-christoph": "PRJNA661613",
    "rothman": "PRJNA729801",
    "spurbeck": "PRJNA924011",
}

with open("/Users/jeffkaufman/code/mgs-pipeline/dashboard/"
          "metadata_samples.json") as inf:
    metadata_samples = json.load(inf)

target_taxids = set()
for taxid in focal_pathogens.values():
    target_taxids.add(taxid)
    for child_taxid in subprocess.check_output([
        "/Users/jeffkaufman/code/jefftk-analysis/"
        "2023-06-07--taxonomic-children.py",
            str(taxid)]).decode('utf-8').split(","):
        target_taxids.add(int(child_taxid))

for bioproject in bioprojects.values():
    matches = {}
    for hvreads in glob.glob("%s-hvreads/*.hvreads.json" % bioproject):
        sample = hvreads.split("/")[-1].split(".")[0]
        if metadata_samples[sample].get("enrichment", "") != "viral":
            continue
        print(hvreads)
        with open(hvreads) as inf:
            for read_name, read_info in json.load(inf).items():
                kraken = read_info[0]
                for hit in kraken.split():
                    taxid, length = hit.split(":")
                    if taxid in ["|", "A"]: continue
                    taxid = int(taxid)
                    if taxid in target_taxids:
                        matches[read_name] = read_info
                        break
    if matches:
        with open("%s-focal2.json" % bioproject, "w") as outf:
            json.dump(matches, outf)
            
