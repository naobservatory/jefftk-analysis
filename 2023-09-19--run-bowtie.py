#!/usr/bin/env python3

import glob
import os.path
import subprocess
from collections import defaultdict

samples = defaultdict(dict)
for fastq in glob.glob("hvfastqs/*.fastq"):
    sample, category, _ = os.path.basename(fastq).rsplit(".", 2)
    samples[sample][category] = fastq

for sample, fastqs in sorted(samples.items()):
    out = "hvsams/%s.sam" % sample
    if os.path.exists(out):
        continue
    
    cmd = [
        "/home/ec2-user/bowtie2-2.5.1-linux-x86_64/bowtie2",
        "--local",
        "-x", "human-viruses",
        "--very-sensitive-local",
        "--score-min", "G,1,0",
        "--mp", "2,0",
        "--no-unal",
        "--threads", "24",
        "-S", out,
    ]
    assert ("pair1" in fastqs) == ("pair2" in fastqs)
    if "pair1" in fastqs:
        cmd.extend(["-1", fastqs["pair1"],
                    "-2", fastqs["pair2"]])
    if "combined" in fastqs:
        cmd.extend(["-U", fastqs["combined"]])

    try:
        subprocess.check_call(cmd)
    except Exception:
        print(" ".join(cmd))
        raise
