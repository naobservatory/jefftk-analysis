#!/usr/bin/env bash

out=~/combined-vb/$(basename $1 | sed s/.kraken2.tsv.gz/.json/)
if [ ! -s $out ]; then
    aws s3 cp $1 - \
        | gunzip \
        | ./2023-08-28--combined-viral-bacterial.py \
              > $out
fi
