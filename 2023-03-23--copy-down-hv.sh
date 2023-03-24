#!/usr/bin/env bash

for study in $(aws s3 ls s3://nao-mgs/ | awk '{print $NF}'); do
    for hv in $(aws s3 ls s3://nao-mgs/${study}humanviruses/ | \
                    awk '{print $NF}'); do
        if [ ! -e $hv ]; then
            echo s3://nao-mgs/${study}humanviruses/$hv
        fi
    done

    for vc in $(aws s3 ls s3://nao-mgs/${study}viruscounts/ | \
                    awk '{print $NF}'); do
        if [ ! -e $vc ]; then
            echo s3://nao-mgs/${study}viruscounts/$vc
        fi
    done
done | xargs -I {} -P 16 aws s3 cp {} .


