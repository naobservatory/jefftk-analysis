#!/usr/bin/env bash
BIOPROJECT="$1"
mkdir -p $BIOPROJECT
cd $BIOPROJECT
aws s3 ls s3://nao-mgs/$BIOPROJECT/cleaned/ | \
    awk '{print $NF}' | \
    grep settings | \
    xargs -I {} -P 32 aws s3 cp s3://nao-mgs/$BIOPROJECT/cleaned/{} .
gunzip *.gz
