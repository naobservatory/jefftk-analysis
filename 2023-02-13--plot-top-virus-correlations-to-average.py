#!/usr/bin/env python3

import sys
import json
import datetime
import itertools
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

# accession -> date, site, nreads
metadata = {}
with open("longest-timeseries-nreads.tsv") as inf:
    for line in inf:
        accession, date, site, nreads = line.strip().split()
        metadata[accession] = (
            datetime.datetime.fromisoformat(date), site, int(nreads))

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
    date for (date, _, _) in metadata.values())

all_sites =  sorted(set(
        site for (_, site, _) in metadata.values()))
first_site, *_, last_site = all_sites

def to_date_color(accession):
    date, _, _ = metadata[accession]
    return (last_date - date) / (last_date - first_date)


site_colors=cm.rainbow([i/(len(all_sites)-1) for i in range(len(all_sites))])
                       
def to_site_color(accession):
    _, site, _ = metadata[accession]
    return site_colors[all_sites.index(site)]

fig, ax = plt.subplots(nrows=3,
                       ncols=3,
                       figsize=(9,9),
                       constrained_layout=True)

def relative_abundance(vid, accession):
    absolute_abundance = accession_counts[accession][vid]
    _, _, total_reads = metadata[accession]
    return absolute_abundance/total_reads

def relative_relative_abundance(vid, accession):
    return relative_abundance(vid, accession) / np.mean([
        relative_abundance(vid, other_accession)
        for other_accession in metadata
        if other_accession != accession])

def mean_relative_relative_abundance(accession):
    return np.mean([relative_relative_abundance(vid, accession)
                    for vid in all_vids])

plt.suptitle("per-sample abundance vs average abundance")
fig.supxlabel("mean relative relative abundance")
fig.supylabel("relative abundance of this virus")
for i, vid1 in enumerate(all_vids):
    slug = to_slug(vid1)

    xs = []
    ys = []
    cs = []

    for accession, counts in accession_counts.items():
        ys.append(relative_relative_abundance(vid1, accession))
        xs.append(mean_relative_relative_abundance(accession))
        cs.append(to_date_color(accession))
        
    a = ax[i//3][i%3]

    a.scatter(xs, ys, c=cs)
    a.set_title("%s" % (slug))

fig.savefig("corr-date-avg.png", dpi=180)

plt.clf()

fig, ax = plt.subplots(nrows=3,
                       ncols=3,
                       figsize=(9,9),
                       constrained_layout=True)

plt.suptitle("per-sample abundance vs average abundance")
fig.supxlabel("mean relative relative abundance")
fig.supylabel("relative abundance of this virus")
for i, vid1 in enumerate(all_vids):
    slug = to_slug(vid1)

    xs = []
    ys = []
    cs = []

    for accession, counts in accession_counts.items():
        ys.append(relative_relative_abundance(vid1, accession))
        xs.append(mean_relative_relative_abundance(accession))
        cs.append(to_site_color(accession))
        
    a = ax[i//3][i%3]

    a.scatter(xs, ys, c=cs)
    a.set_title("%s" % (slug))

fig.savefig("corr-site-avg.png", dpi=180)
