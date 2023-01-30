#!/usr/bin/env bash

ls ~/data/prjna729801/ | grep fasta$ | sed s/.fasta// | while read slug; do
    ls ~/data/prjna729801/reads-$slug/*.fasta | \
        xargs -P8 -I {} bash -c "cat {} | \
          python3 ~/code/jefftk-analysis/2023-01-27--accuracy-by-location.py \
                ~/data/prjna729801/$slug.fasta > \
               {}.locs.tsv"
done
