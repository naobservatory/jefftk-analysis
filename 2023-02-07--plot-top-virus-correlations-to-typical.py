#!/usr/bin/env python3

import sys
import json
import datetime
import itertools
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

# accession -> date, site
metadata = {}
with open("longest-timeseries.tsv") as inf:
    for line in inf:
        accession, date, site = line.strip().split()
        metadata[accession] = (
            datetime.datetime.fromisoformat(date), site)

accession_counts = {} # accession -> vid -> count        
for accession in metadata:
    with open("2023-02-07--top-counts/%s.json" % accession) as inf:
        accession_counts[accession] = json.load(inf)

all_vids = set()
for vids in accession_counts.values():
    all_vids.update(vids)

all_vids = list(sorted(all_vids))

def to_slug(vid):
    return vid.split(".")[0]

first_date, *_, last_date = sorted(
    date for (date, _) in metadata.values())

all_sites =  sorted(set(
        site for (_, site) in metadata.values()))
first_site, *_, last_site = all_sites

def to_date_color(accession):
    date, _ = metadata[accession]
    return (last_date - date) / (last_date - first_date)


site_colors=cm.rainbow([i/(len(all_sites)-1) for i in range(len(all_sites))])
                       
def to_site_color(accession):
    _, site = metadata[accession]
    return site_colors[all_sites.index(site)]

# What are normal ratios?



fig, ax = plt.subplots(nrows=3,
                       ncols=3,
                       figsize=(9,9),
                       constrained_layout=True)

plt.suptitle("per-sample abundance vs average abundance")
fig.supxlabel("mean abundance of other viruses")
fig.supylabel("abundance of this virus")
for i, vid1 in enumerate(all_vids):
    slug = to_slug(vid1)

    xs = []
    ys = []
    cs = []

    for accession, counts in accession_counts.items():
        ys.append(counts[vid1])

        total = 0
        count = 0
        for vid2 in all_vids:
            if vid1 == vid2:
                continue
            total += counts[vid2]
            count += 1

        xs.append(total / count)
        cs.append(to_site_color(accession))
        
    a = ax[i//3][i%3]

    a.scatter(xs, ys, c=cs)
    a.set_title("%s" % (slug))

fig.savefig("corr-avg.png", dpi=180)
plt.clf()
