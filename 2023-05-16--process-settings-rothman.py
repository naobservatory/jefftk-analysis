#!/usr/bin/env python3

import glob
import os.path
import json
from collections import defaultdict
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

with open(os.path.expanduser(
        "~/code/mgs-pipeline/dashboard/metadata_samples.json")) as inf:
    sample_metadata = json.load(inf)

overall = {}
by_site = defaultdict(dict)

# site -> date -> discard_fraction
discard_fraction_by_date_and_site = defaultdict(dict)

for settings in glob.glob("*.settings"):
    sample, _ = os.path.splitext(settings)

    if sample_metadata[sample]["enrichment"] == "panel": continue
    
    in_length_distribution = False
    cols = None

    data = []
    with open(settings) as inf:
        for line in inf:
            if line.startswith("[Length distribution]"):
                in_length_distribution = True
                continue
            if not in_length_distribution:
                continue

            row = line.strip().split("\t")
            if not cols:
                cols = row
                continue

            data.append([int(x) for x in row])

    site = sample_metadata[sample]["fine_location"]
    date = sample_metadata[sample]["date"]

    fig, ax = plt.subplots(constrained_layout=True)
    xs = [r[0] for r in data]
    ys = [r[-4] for r in data]
    ax.plot(xs, ys)
    ax.set_title("Lengths: %s %s" % (site, date))
    fig.savefig("%s.collapsed.lengths.png" % sample, dpi=180)
    plt.clf()
    plt.close()

    fig, ax = plt.subplots(constrained_layout=True)
    xs = [r[0] for r in data]
    ys = [r[1] for r in data]
    ax.plot(xs, ys)
    ax.set_title("Lengths: %s %s" % (site, date))
    fig.savefig("%s.mate1.lengths.png" % sample, dpi=180)
    plt.clf()
    plt.close()

    total_all = 0
    total_discarded = 0
    for length, *rest in data:
        if length not in overall:
            overall[length] = rest
        else:
            for i in range(len(rest)):
                overall[length][i] += rest[i]

        if length not in by_site[site]:
            by_site[site][length] = rest
        else:
            for i in range(len(rest)):
                by_site[site][length][i] += rest[i]

        total_all += rest[-1]
        total_discarded += rest[-2]
        
    discard_fraction_by_date_and_site[site][date] = total_discarded / total_all

fig, ax = plt.subplots(constrained_layout=True)
ax.set_title("Discard fraction by site")
ax.yaxis.set_major_formatter(mtick.PercentFormatter())
for site in sorted(discard_fraction_by_date_and_site):
    xs = []
    ys = []
    for date in sorted(discard_fraction_by_date_and_site[site]):
        xs.append(date)
        ys.append(discard_fraction_by_date_and_site[site][date] * 100)
    ax.scatter(xs, ys, label="%s" % (site))
ax.legend()
fig.savefig("discard-fraction-by-site-over-time.png", dpi=180)
plt.clf()
plt.close()

    
fig, ax = plt.subplots(constrained_layout=True)
ax.set_title("Lengths by Site")
ax.set_yscale("log")
for site in sorted(by_site):
    xs = []
    ys = []
    for length, record in sorted(by_site[site].items()):
        xs.append(length)
        ys.append(record[-1])

    total_y = sum(ys)
    ys = [v/total_y for v in ys]
        
    ax.plot(xs, ys, label="%s" % (site))
ax.legend()
fig.savefig("lengths-by-site.png", dpi=180)
plt.clf()
plt.close()
