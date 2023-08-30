#!/usr/bin/env bash
S3_BUCKET=$(echo "$1" | awk -F "\t" '{print $1}')
BIOPROJECT=$(echo "$1" | awk -F "\t" '{print $2}')
SAMPLE=$(echo "$1" | awk -F "\t" '{print $3}')
SEQ_ID=$(echo "$1" | awk -F "\t" '{print $4}')
KRAKEN_INFO=$(echo "$1" | awk -F "\t" '{print $5}')

OUT=combined-vb-examples/$SAMPLE.$RANDOM$RANDOM.txt
echo "$1" > $OUT
for fasta_gz in $(aws s3 ls s3://$S3_BUCKET/$BIOPROJECT/cleaned/ \
                      | awk '{print $NF}' \
                      | grep -v settings \
                      | grep ^$SAMPLE); do
    aws s3 cp s3://$S3_BUCKET/$BIOPROJECT/cleaned/$fasta_gz - \
        | gunzip \
        | grep ^@"$SEQ_ID" -A 1 \
               >> $OUT
done
