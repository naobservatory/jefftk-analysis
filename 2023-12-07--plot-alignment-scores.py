#!/usr/bin/env python3

import json
from collections import Counter

with open("tmp.Rothman2021.bowtie-scores-after-human-removal") as inf:
    raw_scores = json.load(inf)

scores = {}
adjusted_scores = Counter()    
for rs, count in raw_scores.items():
    score, length = [int(x) for x in rs.split("-")]
    scores[score, length] = count
    
    adjusted_scores[score] += count

import matplotlib.pyplot as plt
fig, ax = plt.subplots()

xs = sorted(x for x in adjusted_scores if x >= 40*2)
ys = [adjusted_scores[x] for x in xs]
print(sum(ys))

ax.plot(xs, ys)
plt.show()
