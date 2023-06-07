#!/usr/bin/env python3

import os
import json
import subprocess

focal_pathogens = {
    # 115 in c-c, 408 in rothman, 18 in spurbeck
    "sars_cov_2": 2697049,
    # one in rothman, one in c-c
    "hiv": 11676,
    # none
    "influenza_a": 11320,
    # none
    "influenza_b": 11520,
    # 6 in c-c, 1600 in rothman, 94 in spurbeck    
    "norovirus_all": 142786,
    # 1 in c-c, 91 in rothman, 15 in spurbeck
    "norovirus_g1": 122928,
    # 4 in c-c, 1300 in rothman, 72 in spurbeck
    "norovirus_g2": 122929,
    # none
    "hcv": 11103,
}

bioprojects = {
    "crits-christoph": "PRJNA661613",
    "rothman": "PRJNA729801",
    "spurbeck": "PRJNA924011",
}

focal_pathogens_full_fname = "focal_pathogens_full.json"

if os.path.exists(focal_pathogens_full_fname):
    with open(focal_pathogens_full_fname) as inf:
        focal_pathogens_full = json.load(inf)
else:
    focal_pathogens_full = {}
    for focal_pathogen_name, focal_pathogen_taxid in focal_pathogens.items():
        focal_pathogens_full[focal_pathogen_name] = [
            int(taxid)
            for taxid in subprocess.check_output([
                    "/Users/jeffkaufman/code/jefftk-analysis/"
                    "2023-06-07--taxonomic-children.py",
                    str(focal_pathogen_taxid)]).decode('utf-8').split(",")]
    with open(focal_pathogens_full_fname, "w") as outf:
        json.dump(focal_pathogens_full, outf)

K=35 # kraken default
for focal_pathogen_name in focal_pathogens:
    target_taxid = focal_pathogens[focal_pathogen_name]
    all_taxids = set(focal_pathogens_full[focal_pathogen_name])
    
    genomes = []

    with open("metadata-%s.tsv" % focal_pathogens[focal_pathogen_name]) as inf:
        cols = None
        for line in inf:
            row = line.strip().split("\t")
            if not cols:
                cols = row
                continue

            genomes.append(row[cols.index("local_filename")])

    for bioproject_name, bioproject_id in bioprojects.items():
        matches = []
        with open("%s-focal2.json" % bioproject_id) as inf:
            for read_name, read_info in json.load(inf).items():
                hit_details =  []
                kraken = read_info[0]
                pos = 0
                fwd_rev = 1
                for hit in kraken.split():
                    if hit == "|:|":
                        pos = 0
                        fwd_rev = 2
                        continue
                    taxid, length = hit.split(":")

                    length = int(length)
                    if taxid == "A": continue
                    
                    taxid = int(taxid)
                    if taxid in all_taxids:
                        match_bases = read_info[fwd_rev][pos:pos+length+K]
                        hit_details.append((taxid, fwd_rev, pos, match_bases))

                    pos += length

                if hit_details:
                    matches.append((
                        read_name,
                        read_info,
                        hit_details,
                    ))

        print(focal_pathogen_name, bioproject_name, len(matches))
        
