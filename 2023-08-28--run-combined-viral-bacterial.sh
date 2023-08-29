#!/usr/bin/env bash

mkdir -p ~/combined-vb/
for bucket in nao-mgs nao-restricted; do
    for bioproject in $(aws s3 ls s3://$bucket/ | awk '{print $NF}'); do
        for kraken in $(aws s3 ls s3://$bucket/${bioproject}processed/ \
                            | awk '{print $NF}'); do
            echo s3://$bucket/${bioproject}processed/$kraken
        done
    done
done \
    | xargs -P 32 -I {} \
            ./2023-08-28--run-combined-viral-bacterial-single.sh {}
