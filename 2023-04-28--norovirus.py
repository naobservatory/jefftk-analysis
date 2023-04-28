#!/usr/bin/env python3

import math
import matplotlib.pyplot as plt
from collections import defaultdict, Counter


state_year_counts = defaultdict(Counter)
state_month_counts = defaultdict(Counter)

# Downloaded from https://wwwn.cdc.gov/norsdashboard/ 2023-04-28.
with open("cdc-nors-outbreak-data.tsv") as inf:
    cols = None
    for line in inf:
        row = line.strip().split("\t")
        if not cols:
            cols = row
            continue
            
        year = row[cols.index("Year")]
        month = row[cols.index("Month")]
        date = "%s-%s" % (year, month.zfill(2))
        state = row[cols.index("State")]
        etiology = row[cols.index("Etiology")]
        illnesses = row[cols.index("Illnesses")]
        if "Norovirus" not in etiology: continue
        genotype = row[cols.index("Serotype or Genotype")]

        state_year_counts[state][year] += 1
        state_month_counts[state][date] += 1

# Ignore any states that haven't reported at least ten times ever.
MIN_COUNT=10
low_states = set()
for state in state_year_counts:
    if sum(state_year_counts[state].values()) < MIN_COUNT:
        low_states.add(state)
for state in low_states:
    del state_year_counts[state]
    del state_month_counts[state]

ncols=5
nrows=math.ceil(len(state_year_counts)/ncols)
        
fig, ax = plt.subplots(nrows=nrows,
                       ncols=ncols,
                       figsize=(6*ncols, 6*nrows),
                       constrained_layout=True)
plt.suptitle("CDC NORS Reports")
fig.supxlabel("time")
fig.supylabel("reports")
for n, state in enumerate(sorted(state_year_counts)):
    a = ax[n // ncols][n % ncols]
    a.set_title(state)
    
    xs = []
    ys = []
    for year in range(2009, 2022):
        xs.append(year)
        ys.append(state_year_counts[state][str(year)])
    a.plot(xs, ys)
fig.savefig("nors-counts-over-time.png")
plt.clf()

fig, ax = plt.subplots(nrows=1,
                       ncols=1,
                       figsize=(12,6),
                       constrained_layout=True)

plt.suptitle("CDC NORS Reports: CA")
fig.supxlabel("time")
fig.supylabel("reports")
state = "California"
xs = []
ys = []
for date, count in sorted(state_month_counts[state].items()):
    xs.append(date)
    ys.append(count)
ax.plot(xs, ys)
fig.savefig("nors-counts-over-time-ca.png")
plt.clf()

fig, ax = plt.subplots(nrows=1,
                       ncols=1,
                       figsize=(12,6),
                       constrained_layout=True)

plt.suptitle("CDC NORS Reports: CA monthly")
fig.supxlabel("month")
fig.supylabel("reports")
state = "California"
for year in range(2009, 2022):
    xs = []
    ys = []
    for month in range(1,13):
        date = "%s-%s" % (year, str(month).zfill(2))
        count = state_month_counts[state][date]
        xs.append(month)
        ys.append(count)
    ax.plot(xs, ys, label=year)
ax.legend()
fig.savefig("nors-counts-monthly-ca.png")
plt.clf()

fig, ax = plt.subplots(nrows=1,
                       ncols=1,
                       figsize=(12,6),
                       constrained_layout=True)

plt.suptitle("CDC NORS Reports: national monthly")
fig.supxlabel("month")
fig.supylabel("reports")
for year in range(2009, 2022):
    xs = []
    ys = []
    for month in range(1,13):
        date = "%s-%s" % (year, str(month).zfill(2))
        count = 0
        for state in state_month_counts:
            count += state_month_counts[state][date]
        xs.append(month)
        ys.append(count)
    ax.plot(xs, ys, label=year)
ax.legend()
fig.savefig("nors-counts-monthly-national.png")
plt.clf()
