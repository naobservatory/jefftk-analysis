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

if [ ! -d hvfastas ]; then
    mkdir hvfastas
    ls hvreads | \
        xargs -P 32 -I {} \
              ~/jefftk-analysis/2023-09-07--json-to-fasta.py hvreads hvfastas {}
fi

if [ ! -e observed-human-virus-taxids.txt ]; then
    ~/jefftk-analysis/2023-09-07--determine-hv-taxids.py \
        hvreads/ \
        ~/mgs-pipeline/human-viruses.tsv \
        observed-human-virus-taxids.txt
fi


~/jefftk-analysis/2023-09-07--get-genomes.py 
