#!/usr/bin/env bash
mkdir -p combined-vb-examples
cat targets.tsv \
    | xargs -P 16 -I {} \
            ~/jefftk-analysis/2023-08-30--extract-examples-single.sh {}
