#!/usr/bin/env python3

import os
import glob
import json
from collections import Counter
import matplotlib.pyplot as plt

metadata = {}
for fname in ["papers", "bioprojects", "samples"]:
    with open("/Users/jeffkaufman/code/mgs-pipeline/dashboard/metadata_%s.json" %
              fname) as inf:
        metadata[fname] = json.load(inf)

for panel_enrichment in [True, False]:
    for paper in metadata["papers"]:
        n_reads = 0
        all_count = 0
        good_count = 0
        good_clip_count = 0
        for bioproject in metadata["papers"][paper]["projects"]:
            for sample in metadata["bioprojects"][bioproject]:
                is_panel_enriched = metadata["samples"][sample].get(
                    "enrichment", "") == "panel"
                if panel_enrichment != is_panel_enriched:
                    continue
                fname = "samcounts/%s.sam.json" % sample
                if os.path.exists(fname):
                    n_reads += metadata["samples"][sample]["reads"]
                    with open(fname) as inf:
                        for key, count in json.load(inf).items():
                            alignment, clip = [int(x) for x in key.split("-")]
                            all_count += count
                            if alignment > 45:
                                good_count += count
                                if clip > 45:
                                    good_clip_count += count
        if n_reads > 0:
            print(paper +
                  (" panel" if panel_enrichment else ""),
                  n_reads, all_count, good_count, good_clip_count,
                  sep="\t")
