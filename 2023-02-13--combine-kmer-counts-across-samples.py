import sys
import json
from collections import Counter

MAX_KMER_LEN=8

kmer_counts = [Counter() for x in range(MAX_KMER_LEN)]

in_fnames = sys.argv[1:]
for fname in in_fnames:
    with open(fname) as inf:
        for k, counts in enumerate(json.load(inf)):
            for kmer, count in counts.items():
                kmer_counts[k][kmer] += count

print(json.dumps(kmer_counts))
