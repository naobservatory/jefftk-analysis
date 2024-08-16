#!/usr/bin/env python3

import random
import numpy as np
from collections import Counter
from statistics import geometric_mean

rng = np.random.default_rng()

concentrations = ["71k/mL", "7.1k/mL", "710/mL"]

# concentration -> number of flagged reads at each of the four junctions
flags =  {
    "710/mL": [9, 1, 5, 7],
    "7.1k/mL": [327, 75, 151, 184],
    "71k/mL": [1097, 238, 508, 700],
}

depths = {
    "710/mL": 2.51e8,
    "7.1k/mL": 1.1e9,
    "71k/mL": 2.39e8,
}

spike_in_relative_abundance = {
    "710/mL": 4e-5,
    "7.1k/mL": 1e-4,
    "71k/mL": 3e-3,
}

# Grimm et al. 2023 estimate for SARS-CoV-2 RAi(1%) with Rothman-style
# sequencing.
target_ra = 1e-7

# How much would we actually spike in relative to what we did spike in, if we
# wanted to hit the target relative abundance?
adjustments = {}
for concentration in concentrations:
    adjustments[concentration] = \
        target_ra / spike_in_relative_abundance[concentration]

ra_flags = {}
for concentration in concentrations:
    ra_flags[concentration] = [
        x/depths[concentration] for x in flags[concentration]]

def simulate_one(concentration, read_depth):
    for ra in ra_flags[concentration]:
        abundance = rng.poisson(ra * read_depth * adjustments[concentration])
        if abundance >= 2:
            return True
    return False

n_simulations_per_scenario = 100000

print("Read Depth", *concentrations, "combined", sep="\t")
for read_depth in np.logspace(np.log10(1e11), np.log10(1e7), 100):
    row = [int(read_depth)]
    for concentration in concentrations + ["combined"]:
        t = 0
        for _ in range(n_simulations_per_scenario):
            chosen_concentration = concentration
            if concentration == "combined":
                chosen_concentration = random.choice(concentrations)
            
            if simulate_one(chosen_concentration, read_depth):
                t += 1
        row.append(t / n_simulations_per_scenario)
    print(*row, sep="\t")
