#!/usr/bin/env bash

if [ ! -d hvreads ]; then
    mkdir hvreads

    for bucket in nao-mgs nao-restricted; do
        for bioproject in $(aws s3 ls s3://$bucket/ | awk '{print $NF}'); do
            for hvreads in $(aws s3 ls s3://$bucket/${bioproject}hvreads/ \
                                 | awk '{print $NF}'); do
                echo s3://$bucket/${bioproject}hvreads/$hvreads
            done
        done
    done | xargs -P 32 -I {} aws s3 cp {} hvreads/
fi

if [ ! -d hvfastqs ]; then
    mkdir hvfastqs
    ls hvreads | \
        xargs -P 32 -I {} \
              ~/jefftk-analysis/2023-09-07--json-to-fasta.py hvreads hvfastqs {}
fi

if [ ! -e observed-human-virus-taxids.txt ]; then
    ~/jefftk-analysis/2023-09-07--determine-hv-taxids.py \
        hvreads/ \
        ~/mgs-pipeline/human-viruses.tsv \
        observed-human-virus-taxids.txt
fi

~/jefftk-analysis/2023-09-07--get-genomes.py 

if [ ! -d raw-genomes ]; then
    mkdir raw-genomes
    for x in $(find refseq/ | grep gz$); do
        zcat "$x" > raw-genomes/$(basename ${x/.fna.gz/.fna})
    done
fi

if [ ! -e human-viruses.1.bt2 ]; then
    ~/bowtie2-2.5.1-linux-x86_64/bowtie2-build \
        -f \
        --threads 32 \
        $(find raw-genomes/ | grep .fna$ | tr '\n' ',') \
        human-viruses
fi

if [ ! -e hvsams ]; then
    mkdir hvsams
    ~/jefftk-analysis/2023-09-19--run-bowtie.py
fi
