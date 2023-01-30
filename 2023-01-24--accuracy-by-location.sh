#!/usr/bin/env bash

#SLUG=tbrv
SLUG=pmmov

ls ~/data/prjna729801/$SLUG/*.fasta | \
    xargs -P8 -I {} bash -c "cat {} | \
      python3 ~/code/jefftk-analysis/2023-01-24--accuracy-by-location.py \
              ~/data/prjna729801/$SLUG/seq > \
             {}.locs.tsv"

