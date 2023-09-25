#!/usr/bin/env python3

import os
import glob
import json
import pysam
from collections import Counter

for samfile in glob.glob("hvsams/*.sam"):
    sample = os.path.basename(samfile)

    counts = Counter()
    
    with pysam.AlignmentFile(samfile, "r")  as sam:
        for record in sam:
            total_length = sum(length for _, length in record.cigartuples)
            longest_alignment = max(
                [length for category, length in record.cigartuples
                 if category == 0], default=0)
            longest_soft_clip = max(
                [length for category, length in record.cigartuples
                 if category == 4], default=0)

            counts["%s-%s" % (longest_alignment, longest_soft_clip)] += 1

            #print(total_length, longest_alignment, longest_soft_clip,
            #      record.cigarstring)
        
    with open("samcounts/%s.json" % sample, "w") as outf:
        json.dump(counts, outf)
