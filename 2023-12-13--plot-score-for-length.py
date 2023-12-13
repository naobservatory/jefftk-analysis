#!/usr/bin/env python3

import json
import math
import matplotlib.pyplot as plt

with open("evaluation_notes.json") as inf:
    evaluations = json.load(inf)

fig, ax = plt.subplots()

for is_hv in [True, False]:
    xs = []
    ys = []

    for sample in evaluations:
        for read_id in evaluations[sample]:
            r = evaluations[sample][read_id]
            if "hv" not in r:
                continue

            is_hv_val = r["hv"]
            notes = r["bowtie-on-kraken-unclassified"]
            if "non-aligned" in notes:
                continue
            #if "kraken_assigned_nonhv" in notes:
            #    continue
            if "human" in notes:
                continue


            score = notes["hv_score"]
            length = notes["hv_length"]

            if is_hv_val != is_hv:
                continue
            if not notes.get("collapsed"):
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

for sample in evaluations:
    for read_id in evaluations[sample]:
        r = evaluations[sample][read_id]
        if "hv" not in r:
            continue

        is_hv_val = r["hv"]
        notes = r["bowtie-on-kraken-unclassified"]
        if "non-aligned" in notes:
            continue
        score = notes["hv_score"]
        length = notes["hv_length"]

        if notes.get("collapsed"):
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
