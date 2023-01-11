#!/usr/bin/env python3

import sys
import math
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np
from collections import defaultdict, Counter
from matplotlib import rcParams

rcParams.update({'figure.autolayout': True})

metadata, fname_out_base = sys.argv[1:]

wtps = []
wtp_dates = defaultdict(list) # wtp -> [date1, ...\]
labels = []
with open(metadata) as inf:
    for line in inf:
        accession, date, wtp = line.strip().split()
        labels.append((wtp, date))
        wtp_dates[wtp].append(date)
        if wtp not in wtps:
            wtps.append(wtp)

global_corr_total_unweighted = Counter()
global_corr_total_mean_weighted = Counter()
global_corr_total_mean2_weighted = Counter()
global_corr_count = Counter()

wtp_corr_total_unweighted = Counter()
wtp_corr_total_mean_weighted = Counter()
wtp_corr_total_mean2_weighted = Counter()
wtp_corr_count = Counter()

for line in sys.stdin:
    line = line.strip()
    if line.isdigit():
        continue # skip timestamps
    (wtp_idx_1,
     wtp_idx_2,
     sample_idx_1,
     sample_idx_2,
     global_corr_total_unweighted_v,
     global_corr_total_mean_weighted_v,
     global_corr_total_mean2_weighted_v,
     global_corr_count_v,
     wtp_corr_total_unweighted_v,
     wtp_corr_total_mean_weighted_v,
     wtp_corr_total_mean2_weighted_v,
     wtp_corr_count_v,
    ) = line.split()

    wtp_idx_1, wtp_idx_2, sample_idx_1, sample_idx_2 = [
        int(x) for x in [wtp_idx_1, wtp_idx_2, sample_idx_1, sample_idx_2]]

    key = wtp_idx_1, wtp_idx_2, sample_idx_1, sample_idx_2
    global_corr_total_unweighted[key] += float(global_corr_total_unweighted_v)
    global_corr_total_mean_weighted[key] += float(global_corr_total_mean_weighted_v)
    global_corr_total_mean2_weighted[key] += float(global_corr_total_mean2_weighted_v)
    global_corr_count[key] += float(global_corr_count_v)
    wtp_corr_total_unweighted[key] += float(wtp_corr_total_unweighted_v)
    wtp_corr_total_mean_weighted[key] += float(wtp_corr_total_mean_weighted_v)
    wtp_corr_total_mean2_weighted[key] += float(wtp_corr_total_mean2_weighted_v)
    wtp_corr_count[key] += float(wtp_corr_count_v)

def div_or_zero(num, denom):
    if denom == 0:
        return 0
    return num/denom

unweighted_data = []
mean_weighted_data = []
mean2_weighted_data = []
wtp_datas = defaultdict(list) # wtp -> data for wtp

for y_pos, (wtp1, date1) in enumerate(labels):
    unweighted_row = []
    mean_weighted_row = []
    mean2_weighted_row = []
    wtp_row = []

    for x_pos, (wtp2, date2) in enumerate(labels):
        key = (
            wtps.index(wtp1),
            wtps.index(wtp2),
            wtp_dates[wtp1].index(date1),
            wtp_dates[wtp2].index(date2))

        if x_pos > y_pos and False:
            # make it triangular
            unweighted_row.append(math.nan)
            mean_weighted_row.append(math.nan)
            mean2_weighted_row.append(math.nan)
            continue

        count = global_corr_count[key]
        unweighted_row.append(
            div_or_zero(global_corr_total_unweighted[key], count))
        mean_weighted_row.append(
            div_or_zero(global_corr_total_mean_weighted[key], count))
        mean2_weighted_row.append(
            div_or_zero(wtp_corr_total_mean2_weighted[key], count))

        if wtp1 == wtp2:
            # unweighted only
            wtp_row.append(div_or_zero(
                wtp_corr_total_unweighted[key],
                wtp_corr_count[key]))

    unweighted_data.append(unweighted_row)
    mean_weighted_data.append(mean_weighted_row)
    mean2_weighted_data.append(mean2_weighted_row)
    wtp_datas[wtp1].append(wtp_row)

text_labels = ["%s %s" % label for label in labels]

datasets = [
    [unweighted_data, "unweighted", text_labels],
    [mean_weighted_data, "mean weighted", text_labels],
    [mean2_weighted_data, "mean2 weighted", text_labels],
]

for wtp, wtp_data in wtp_datas.items():
    datasets.append([wtp_data, "%s unweighted" % wtp,
                     [label_date
                      for (label_wtp, label_date)
                      in labels
                      if label_wtp == wtp]])

for data, name, use_labels in datasets[:]:
    log_data = []
    for row in data:
        log_row = []
        for col_idx in range(len(row)):
            abs_v = abs(row[col_idx])
            sign_v = 1 if row[col_idx] > 0 else -1
            if abs_v < 1:
                val = 0
            else:
                val = sign_v * math.log(abs_v, 10)
            log_row.append(val)
        log_data.append(log_row)
    datasets.append([log_data, "%s log" % name, use_labels])
    
for data, name, use_labels in datasets:
    df = pd.DataFrame(data, columns=use_labels, index=use_labels)
    sns.set(font_scale=0.5)
    heatmap = sns.heatmap(df, cmap=sns.color_palette("vlag", as_cmap=True),
                          center=0, xticklabels=True, yticklabels=True)

    if len(use_labels) == len(labels):
        prev_wtp, _ = labels[0]
        for i, (wtp, date) in enumerate(labels):
            if wtp != prev_wtp:
                heatmap.axhline(i, color="black")
                heatmap.axvline(i, color="black")
            prev_wtp = wtp
    heatmap.set_title("Sample correlations, %s" % name)
    heatmap.figure.savefig(fname_out_base + "_%s.png" %
                           name.replace(" ", "_"), dpi=180)
    plt.clf()
