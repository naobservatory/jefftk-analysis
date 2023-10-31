#!/usr/bin/env python3

import os
import glob
import json
from collections import Counter

MGS_PIPELINE_DIR="/Users/jeffkaufman/code/mgs-pipeline"
MGS_RESTRICTED_DIR="/Users/jeffkaufman/code/mgs-restricted"

sample_viral = Counter()
sample_both = Counter()
for fname in glob.glob("combined-vb/*.json"):
    with open(fname) as inf:
        sample, *_ = os.path.basename(fname).split(".")
        
        r = json.load(inf)
        sample_both[sample] += r["n_both"]
        sample_viral[sample] += r["n_viral"]
        
for sample in sorted(sample_viral):
    print(sample, sample_both[sample], sample_viral[sample],
          sample_viral[sample]/ sample_both[sample], sep="\t")

