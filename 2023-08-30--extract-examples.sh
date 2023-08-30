#!/usr/bin/env bash
mkdir -p combined-vb-examples
cat targets.tsv \
    | head -n 8 \
    | xargs -P 16 -I {} \
            ~/code/jefftk-analysis/2023-08-30--extract-examples-single.sh {}
