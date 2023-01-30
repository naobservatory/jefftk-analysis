#!/usr/bin/env python3

import sys

plot_prefix, *fnames = sys.argv[1:]

loc_observations_full_matches = []
loc_observations_any_matches = []
loc_transitions = []

for fname in fnames:
    with open(fname) as inf:
        for line in inf:
            loc, obs_full, obs_any, obs_trans = [
                int(x) for x in line.strip().split()]

            if loc >= len(loc_observations_full_matches):
                loc_observations_full_matches.append(0)
                loc_observations_any_matches.append(0)
                loc_transitions.append(0)
            
            loc_observations_full_matches[loc] += obs_full
            loc_observations_any_matches[loc] += obs_any
            loc_transitions[loc] += obs_trans

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

fig, ax = plt.subplots()

ax.set_ylabel("reads")
ax.set_xlabel("position along genome")
ax.set_title("reads by location")            
loc_percentages = [x / len(loc_observations_full_matches) * 100
                   for x in range(len(loc_observations_full_matches))]

ax.plot(loc_percentages,
        loc_observations_full_matches,
        label="exact matches only")
ax.plot(loc_percentages,
        loc_observations_any_matches,
        label="any k-mer matches")
ax.xaxis.set_major_formatter(mtick.PercentFormatter())
ax.legend()
fig.savefig("%s-matches.png" % plot_prefix, dpi=180)

plt.clf()
fig, ax = plt.subplots()

ax.set_ylabel("reads with a disagreement at this position")
ax.set_xlabel("position along genome")
ax.set_title("genome disagreements by location")

loc_percentages = [x / len(loc_transitions) * 100
                   for x in range(len(loc_transitions))]
ax.bar(loc_percentages,
       loc_transitions,
       width=0.1)
ax.xaxis.set_major_formatter(mtick.PercentFormatter())
fig.savefig("%s-disagreements.png" % plot_prefix, dpi=180)

