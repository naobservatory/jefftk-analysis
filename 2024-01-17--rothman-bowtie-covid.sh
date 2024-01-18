#!/usr/bin/env bash

# Rothman 2021
METADATA=~/mgs-pipeline/bioprojects/PRJNA729801/metadata/metadata.tsv

for fname in $(aws s3 ls s3://nao-mgs/PRJNA729801/alignments/ \
                   | awk '{print $NF}' \
                   | grep hv.alignments.tsv.gz); do
    sample=${fname/.hv.alignments.tsv.gz}
    enrichment=$(cat $METADATA  | awk -F '\t' '($1=="'$sample'"){print $NF}')
    if [[ $enrichment != 0 ]]; then
       continue
    fi
    hv_out=$sample.covid.alignments.tsv
    if [[ ! -f $hv_out ]]; then
        # Covid is 2697049 and the taxonomy doesn't currently have any children
        # for it so it's enough to look only for that taxid.
        aws s3 cp s3://nao-mgs/PRJNA729801/alignments/$fname - \
            | gunzip \
            | awk -F'\t' '($3==2697049){print}' > $hv_out
    fi

    human_out=$sample.covid-human.alignments.tsv
    if [[ ! -f $human_out ]]; then
        aws s3 cp s3://nao-mgs/PRJNA729801/alignments/${fname/.hv./.human.} - \
            | gunzip \
            | ~/jefftk-analysis/2024-01-17--get-aligment-lines-matching-ids.py \
                  $hv_out 0 \
                  > $human_out
        touch $human_out
    fi

    kraken_out=$sample.kraken.tsv
    if [[ ! -f $kraken_out ]]; then
        for kraken_fname in $(aws s3 ls \
                                  s3://nao-mgs/PRJNA729801/processed/$sample \
                                  | awk '{print $NF}'); do
            aws s3 cp s3://nao-mgs/PRJNA729801/processed/$kraken_fname - \
                | gunzip \
                | ~/jefftk-analysis/2024-01-17--get-aligment-lines-matching-ids.py \
                      $hv_out 1 \
                      >> $kraken_out
        done
        touch $kraken_out
    fi
done
    
for fname in *.covid.alignments.tsv; do
    sample=${fname/.covid.alignments.tsv}
    if [[ ! -f $sample.stats.json ]]; then
        ~/jefftk-analysis/2024-01-17--rothman-bowtie-covid.py $sample
    fi
done

~/jefftk-analysis/2024-01-17--rothman-bowtie-covid-summary.py
