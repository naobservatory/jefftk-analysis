#!/usr/bin/env python3

import time
import subprocess

results_fname = "/home/ec2-user/null-n-threads.txt"

while True:
    start = time.time()

    print("running...")
    subprocess.check_call(
        "paste"
        "  <(aws s3 cp"
        "      's3://nao-mgs/PRJNA661613/raw/SRR23998356_1.fastq.gz' -"
        "      | gunzip)"
        "  <(aws s3 cp"
        "      's3://nao-mgs/PRJNA661613/raw/SRR23998356_2.fastq.gz' -"
        "      | gunzip) > /dev/null",
            shell=True,
            executable='/bin/bash')

    end = time.time()

    with open(results_fname, "a") as outf:
        outf.write("%.0f\n" % (end - start))
