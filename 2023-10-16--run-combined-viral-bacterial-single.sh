#!/usr/bin/env bash

out=combined-vb/$(basename "$1" | sed s/.kraken2.tsv.gz/.json/)
echo $out
if [ ! -s $out ]; then
    ~/code/jefftk-analysis/2023-10-16--combined-viral-bacterial.py $1 $out
fi
