#!/usr/bin/env python3

import time
import subprocess

results_fname = "/home/ec2-user/bowtie-n-threads.txt"

while True:
    for n_threads in range(1, 29):
        start = time.time()

        print("%s..." % n_threads)
        subprocess.check_call(
            "/home/ec2-user/bowtie2-2.5.2-linux-x86_64/bowtie2"
            " --threads %s"
            " --no-unal"
            " --no-sq"
            " -S /dev/null"
            " -x /home/ec2-user/mgs-pipeline/bowtie/chm13.draft_v1.0_plusY"
            " -1 <(aws s3 cp"
            "      's3://nao-mgs/PRJNA661613/raw/SRR23998356_1.fastq.gz' -"
            "      | gunzip)"
            " -2 <(aws s3 cp"
            "      's3://nao-mgs/PRJNA661613/raw/SRR23998356_2.fastq.gz' -"
            "      | gunzip)" % n_threads,
            shell=True,
            executable='/bin/bash')

        end = time.time()

        with open(results_fname, "a") as outf:
            outf.write("%s %.0f\n" % (n_threads, end - start))
