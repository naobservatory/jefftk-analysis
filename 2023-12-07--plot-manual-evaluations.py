#!/usr/bin/env python3

import json
import math

with open("bowtie-scores-categorized") as inf:
    evaluations = json.load(inf)

evaluations = [x for x in evaluations if "ERR3563070" in x[0]]
#evaluations = [x for x in evaluations
#               if "SRR14530891" in x[0]
#               or "SRR14530882" in x[0]
#               or "SRR14530771" in x[0]]
    
import matplotlib.pyplot as plt
fig, ax = plt.subplots()

for is_hv in [True, False]:
    xs = []
    ys = []
    for read_id, score, length, is_hv_val in evaluations:
        if is_hv_val != is_hv:
            continue
        if not read_id.startswith("M_"):
            continue
        xs.append(length)
        ys.append(score)
    
    ax.scatter(xs, ys,
               label="hv" if is_hv else "no",
               marker="x" if is_hv else ",")

xs = [25, 40, 80, 200]
ys = [60, 60, 95, 100]
ax.plot(xs, ys, label="div")
ax.legend()
ax.set_xlabel("length")
ax.set_ylabel("score")

plt.show()
plt.clf()

uncollapsed_scores_by_is_hv = {
    True: [],
    False: [],
}

for read_id, score, length, is_hv_val in evaluations:
    if read_id.startswith("M_"):
        continue
    uncollapsed_scores_by_is_hv[is_hv_val].append(score)

uncollapsed_scores_by_is_hv[True].sort()
uncollapsed_scores_by_is_hv[False].sort()

best_cutoff = None
best_overall_score = None
for cutoff in \
        set(uncollapsed_scores_by_is_hv[True]) | \
        set(uncollapsed_scores_by_is_hv[False]):
    overall_score = 0
    for score in uncollapsed_scores_by_is_hv[True]:
        if score < cutoff:
            overall_score -= 1
        else:
            overall_score += 1
    for score in uncollapsed_scores_by_is_hv[False]:
        if score < cutoff:
            overall_score += 1
        else:
            overall_score -= 1
    if best_overall_score is None or overall_score > best_overall_score:
        best_overall_score = overall_score
        best_cutoff = cutoff

print("Full length cutoff: %s (overall score=%s)" % (
    best_cutoff, best_overall_score))

fig, ax = plt.subplots()
for is_hv in [True, False]:
    xs = []
    ys = []
    for i, score in enumerate(uncollapsed_scores_by_is_hv[is_hv]):
        xs.append(i)
        ys.append(score)
    ax.scatter(xs, ys,
               label="hv" if is_hv else "no",
               marker="x" if is_hv else ",")
    
xs = [0, max(len(uncollapsed_scores_by_is_hv[True]),
             len(uncollapsed_scores_by_is_hv[False]))]
ys = [best_cutoff, best_cutoff]
ax.plot(xs, ys, label="div")
ax.set_xlabel("length")
ax.set_ylabel("y-axis is not meaningful")
ax.legend()   
plt.show()
plt.clf()
