#!/usr/bin/env python3
import os
import re
import gzip
import json
import subprocess
from collections import defaultdict

# output of p2ra/list_taxids.py
taxid_list="""
10298           hsv_1           HSV-1
10310           hsv_2           HSV-2
10359           cmv             CMV
10376           ebv             EBV
10407           hbv             HBV
10566           hpv             HPV
10632           jcv             JCV
10804           aav2            AAV2
11103           hcv             HCV
11320           influenza       Influenza A
11520           influenza       Influenza B
11676           hiv             HIV
68558           aav6            AAV6
82300           aav5            AAV5
122928          norovirus       Norovirus (GI)
122929          norovirus       Norovirus (GII)
493803          mcv             MCV
1891762         bkv             BKV
2697049         sars_cov_2      SARS-COV-2
"""

taxid_to_name = {}
for line in taxid_list.split("\n"):
    line = line.strip()
    if not line: continue

    taxid, human_readable = re.match(
        "^([0-9]+) +[a-z0-9_]+ +(.*)$", line).groups()
    taxid_to_name[int(taxid)] = human_readable

metadata_fname = "metadata.txt"
our_taxid_to_all_fname = "our_taxid_to_all.json"
taxid_to_fastas = defaultdict(list)
with open(metadata_fname) as inf:
    cols = None
    for line in inf:
        rows = line.strip().split("\t")
        if not cols:
            cols = rows
            continue

        taxid_to_fastas[int(rows[cols.index("taxid")])].append(
            rows[cols.index("local_filename")])

with open(our_taxid_to_all_fname) as inf:
    our_taxid_to_all = json.load(inf)

READ_LENGTH=100
def to_reads(s):
    for i in range(len(s) - READ_LENGTH + 1):
        yield s[i:i+READ_LENGTH]

def fake_reads(fake_reads_fname, our_taxid):
    with gzip.open(fake_reads_fname, "w") as outf:
        for any_taxid in our_taxid_to_all[str(our_taxid)]:
            if any_taxid in taxid_to_fastas:
                for fasta in sorted(taxid_to_fastas[any_taxid]):
                    with gzip.open(fasta) as inf:
                        raw_genome = []
                        for line in inf.read().decode("utf-8").split("\n"):
                            if line.startswith(">"):
                                continue
                            raw_genome.append(line.strip())
                        for n, simulated_read in enumerate(
                                to_reads("".join(raw_genome))):
                            outf.write(("@%s:%s\n%s\n+\n%s\n" % (
                                os.path.basename(
                                    fasta.replace("_genomic.fna.gz", "")),
                                n,
                                simulated_read,
                                "F"*(len(simulated_read)))
                                        ).encode("utf-8"))

def classify_reads(fake_reads_fname, classified_fname):
    subprocess.check_call([
        "/home/ec2-user/kraken2/kraken2",
        "--db", "/run/kraken-db/",
        "--output", classified_fname,
        "--use-names",
        fake_reads_fname
    ])

for our_taxid in taxid_to_name:
    if our_taxid == 68558:
        continue # no genomes

    alias_taxids = our_taxid_to_all[str(our_taxid)]

    fake_reads_fname = "%s.fake.reads.fastq.gz" % our_taxid
    if not os.path.exists(fake_reads_fname):
        fake_reads(fake_reads_fname, our_taxid)

    classified_fname = fake_reads_fname + ".classified"
    if not os.path.exists(classified_fname):
        classify_reads(fake_reads_fname, classified_fname)

    total = 0
    classification_counts = defaultdict(int)
    with open(classified_fname) as inf:
        for line in inf:
            total += 1
            desc = line.split("\t")[2]
            classified_taxid = int(desc.split("(taxid ")[-1].replace(")", ""))
            if classified_taxid in alias_taxids:
                desc = "%s (taxid %s)" % (taxid_to_name[our_taxid], our_taxid)
            classification_counts[desc] += 1

    count_classification = [
        (count, classification)
        for (classification, count) in classification_counts.items()]
    count_classification.sort(reverse=True)

    print("%s (%s)" % (taxid_to_name[our_taxid], our_taxid))
    for count, classification in count_classification:
        print ("  %s (%.1f%%)\t%s" % (
            count, 100*count/total, classification))
