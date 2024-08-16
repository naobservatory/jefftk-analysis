#!/usr/bin/env python3

import numpy as np
from collections import Counter
from statistics import geometric_mean

rng = np.random.default_rng()

# concentration -> number of flagged reads at each of the four junctions
flags =  {
    "710/mL": [9, 1, 5, 7],
    "7100/mL": [327, 75, 151, 184],
    "71000/mL": [1097, 238, 508, 700],
}

depths = {
    "710/mL": 2.51e8,
    "7100/mL": 1.1e9,
    "71000/mL": 2.39e8,
}

ra_flags = {}
for concentration in flags:
    ra_flags[concentration] = [
        x/depths[concentration] for x in flags[concentration]]

def simulate_one(concentration, read_depth):
    for ra in ra_flags[concentration]:
        abundance = rng.poisson(ra * read_depth)
        if abundance >= 2:
            return True
    return False

n_simulations_per_scenario = 100000

print("Read Depth", *sorted(flags), sep="\t")
for read_depth in np.logspace(np.log10(1e9), np.log10(1e3), 100):
    row = [int(read_depth)]
    for concentration in sorted(flags):
        t = 0
        for _ in range(n_simulations_per_scenario):
            if simulate_one(concentration, read_depth):
                t += 1
        row.append(t / n_simulations_per_scenario)
    print(*row, sep="\t")
