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

fig, ax = plt.subplots(nrows=len(all_vids),
                       ncols=len(all_vids),
                       #sharex=True, sharey=True,
                       figsize=(2.5*len(all_vids), 2*len(all_vids)),
                       constrained_layout=True)

plt.suptitle("per-sample abundance correlations, date-colored")
fig.supxlabel("abundance")
fig.supylabel("abundance")
for row, vid1 in enumerate(all_vids):
    slug1 = to_slug(vid1)
    for col, vid2 in enumerate(all_vids):
        slug2 = to_slug(vid2)
        a = ax[row, col]

        xs = []
        ys = []
        cs = []
        for accession, counts in accession_counts.items():
            xs.append(counts[vid1])
            ys.append(counts[vid2])
            cs.append(to_date_color(accession))

        a.scatter(xs, ys, c=cs)
        a.set_title("%s-%s" % (slug1, slug2))

fig.savefig("corr-date.png", dpi=180)

plt.clf()

fig, ax = plt.subplots(nrows=len(all_vids),
                       ncols=len(all_vids),
                       #sharex=True, sharey=True,
                       figsize=(2.5*len(all_vids), 2*len(all_vids)),
                       constrained_layout=True)

plt.suptitle("per-sample abundance correlations, site-colored")
fig.supxlabel("abundance")
fig.supylabel("abundance")
for row, vid1 in enumerate(all_vids):
    slug1 = to_slug(vid1)
    for col, vid2 in enumerate(all_vids):
        slug2 = to_slug(vid2)
        a = ax[row, col]

        xs = []
        ys = []
        cs = []
        for accession, counts in accession_counts.items():
            xs.append(counts[vid1])
            ys.append(counts[vid2])
            cs.append(to_site_color(accession))

        a.scatter(xs, ys, c=cs)
        a.set_title("%s-%s" % (slug1, slug2))

fig.savefig("corr-site.png", dpi=180)


