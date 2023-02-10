#!/usr/bin/env python3

import sys
import glob
import json
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

samples = {} # accession_{1,2} -> vid -> f/r + loc -> count
for fname in sys.argv[1:]: 
    sample, _ = fname.split("_")
    with open(fname) as inf:
        samples[sample] = json.load(inf)

all_vids = set()
for vids in samples.values():
    all_vids.update(vids)
all_vids = list(sorted(all_vids))

def to_slug(vid):
    return vid.split(".")[0]
    
fig, axs = plt.subplots(nrows=3,
                        ncols=3,
                        #sharex=True, sharey=True,
                        figsize=(12, 12),
                        constrained_layout=True)

def predict_helper(vid, sample, direction):
    other_vals = Counter()
    this_vals = Counter()
    for other_sample in samples:
        if sample == other_sample: continue
        for frloc, count in samples[other_sample][vid].items():
            if frloc.startswith(direction):
                loc = int(frloc[1:])
                other_vals[loc] += count
    for frloc, count in samples[sample][vid].items():
        if frloc.startswith(direction):
            loc = int(frloc[1:])
            this_vals[loc] += count

    predictions = []
    errors = []

    n_locs = max(*this_vals, *other_vals) + 1
    for loc in range(n_locs):
        # fake_prediction = total_this * other_vals[loc] / total_other
        
        near_total_this = sum(v for (l, v) in this_vals.items() if l != loc)
        near_total_other = sum(v for (l, v) in other_vals.items() if l != loc)
        real_prediction = near_total_this * other_vals[loc] / near_total_other
        
        predictions.append(real_prediction)
        errors.append(real_prediction - this_vals[loc])

    errors = np.array(errors)
    msq_errors = errors**2

    return np.array(predictions), msq_errors

def predict(vid, direction):
    global_predictions = None
    global_sq_errors = None
    for sample in samples:
        predictions, sq_errors = predict_helper(vid, sample, direction)

        if global_predictions is None:
            global_predictions = predictions
        else:
            global_predictions += predictions

        if global_sq_errors is None:
            global_sq_errors = sq_errors
        else:
            global_sq_errors += sq_errors

    global_predictions *= 1/len(samples)
    global_rms_errors = np.sqrt(global_sq_errors * 1/len(samples))

    return global_predictions, global_rms_errors        

plt.suptitle("read starts by location")
fig.supxlabel("position along genome")
fig.supylabel("prediction error for whether reads start at this position")
for row, vid in enumerate(all_vids):
    slug = to_slug(vid)
    ax = axs[row%3][row//3]

    for direction in ['f', 'r']:
        _, global_rms_errors = predict(vid, direction)

        xs = []
        ys = []
        for loc in range(len(global_rms_errors)):
            xs.append(loc)
            ys.append(global_rms_errors[loc])

        ax.scatter(xs, ys, label=direction, marker=',', s=0.5)
        ax.legend()
        ax.set_title("%s" % slug)

fig.savefig("read-start-prediction-errors.png", dpi=180)
plt.clf()
