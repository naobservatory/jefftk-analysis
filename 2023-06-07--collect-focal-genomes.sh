#!/usr/bin/env bash

for taxid in 2697049 11676 11320 11520 142786 122928 122929 11103; do
    all_taxids=$(~/code/jefftk-analysis/2023-06-07--taxonomic-children.py $taxid)
    ncbi-genome-download \
        --taxids $all_taxids \
        --metadata-table metadata-$taxid.tsv \
        --formats fasta \
        viral \
        --flat-output \
        -v \
        --output-folder refseq
done
