#!/usr/bin/env python3

import re
import datetime
from Bio.SeqIO.QualityIO import FastqGeneralIterator

SAMPLE_RATE_HZ=5000

def ts(yyyy, mt, dd, hr, mn, sec, ms):
    return datetime.datetime(
        int(yyyy), int(mt), int(dd),
        int(hr), int(mn), int(sec),
        int(ms)).timestamp()

def process(inf):
    total_reads = 0
    total_bases = 0
    total_duration = 0

    for line in inf:
        if not line.startswith("["):
            continue

        yyyy, mt, dd, hr, mn, sec, ms, level, log = re.match(
            "\[(\d\d\d\d)-(\d\d)-(\d\d) (\d\d?):(\d\d):(\d\d)[.](\d\d\d)\]"
            " \[([a-z]*)\] (.*)", line).groups()

        if "downloading" in log:
            assert "_400bps_" in log
            speed_bps = 400
        elif "Running" in log:
            fname = log.split('"')[5]
            sample =  fname.replace("pod5_samples/", "").replace(".pod5", "")
            fastq = sample + ".fastq"
            n_bases = 0
            n_reads = 0
            with open(fastq) as inf2:
                for (title, sequence, quality) in FastqGeneralIterator(inf2):
                    n_bases += len(sequence)
                    n_reads += 1

        elif "batch size" in log:
            start_ts = ts(yyyy, mt, dd, hr, mn, sec, ms)

        elif "Basecalled @ Samples/s" in log:
            samples_s = float(log.replace(
                "> Basecalled @ Samples/s:", "").strip())
            end_ts = ts(yyyy, mt, dd, hr, mn, sec, ms)

        elif "Finished" in log:
            duration = end_ts - start_ts
            if False:
                print("%s:" % sample)
                print("  %s reads, %s bases, %.0f mean len, %.0fs, %.0f bps" % (
                    n_reads, n_bases, n_bases / n_reads, duration,
                    n_bases / duration))

            total_reads += n_reads
            total_bases += n_bases
            total_duration += duration

    print("%.0f\t%.0f" % (
        total_bases / total_duration,
        total_duration / total_bases * 48_000_000_000 / 60 / 60))

print("system", "model", "bps", "hr/48Gb", sep="\t")
for system in ["g5.xlarge", "mac"]:
    for speed in ["fast", "hac", "sup"]:
        print(system, speed, sep="\t", end="\t")
        with open("results.%s/dorado-output-%s.txt" % (system, speed)) as inf:
            process(inf)
