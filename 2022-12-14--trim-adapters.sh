#!/bin/bash

adapter_trimming=../wastewater_viromics_sarscov2/adapter_trimming
determine_adapter=${adapter_trimming}/determine_adapter.py

for fname in *_R1.fastq.gz; do
    in1=$fname
    in2=${fname/R1/R2}
    basename=${fname/_R1.fastq.gz}
    AdapterRemoval \
         --file1 $in1 \
         --file2 $in2 \
         --basename $basename \
         --trimns \
         --trimqualities \
         --collapse \
         --interleaved-output \
         --adapter1 $(cat $in1 | gunzip | $determine_adapter - fwd) \
         --adapter2 $(cat $in2 | gunzip | $determine_adapter - rev)
done

