#!/usr/bin/env bash

REMOTE="$1"

for wtp in OC SJ HTP JWPCP ; do
    fname="clean-TS-ss-$wtp-mvl-uniq-c.tsv"
    if [ ! -e $fname ] ; then
        aws s3 cp s3://prjna729801/clean-TS-ss-$wtp-mvl-uniq-c.tsv .
    fi
    ~/code/jefftk-analysis/2023-01-17--kmer-plotting.py \
        $wtp $fname $wtp-1000.png 1000 500000
    ~/code/jefftk-analysis/2023-01-17--kmer-plotting.py \
        $wtp $fname $wtp-100.png 100 10000
    ~/code/jefftk-analysis/2023-01-17--kmer-plotting.py \
        $wtp $fname $wtp-40.png 40 400

    fname="tomato.brown.rugose.nt.seq.$wtp.uniq-c.mvl"
    if [ ! -e $fname ] ; then
        scp "$REMOTE:kmer-egd/$fname" .
    fi
    ~/code/jefftk-analysis/2023-01-17--kmer-plotting.py \
        $wtp $fname $wtp-tbrv-1000.png 1000 500000
    
done
