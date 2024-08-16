#!/usr/bin/env python3

import numpy as np
from collections import Counter
from statistics import geometric_mean

rng = np.random.default_rng()

# redo all of this with
#  - poisson
#  - relative abundances
#  - target concentrations

# concentration -> number of flagged reads at each of the four junctions
flags =  {
    "710/mL": [9, 1, 5, 7],
    "7100/mL": [327, 75, 151, 184],
    "71000/mL": [1097, 238, 508, 700],
}

def simulate_one(concentration, subsampling):
    real_flags = flags[concentration]
    n_real = sum(flags[concentration])
    n_simulated = round(n_real * subsampling)

    if n_simulated < 1:
        return False

    probabilities = [n / n_real for n in real_flags]

    observations = Counter()
    for observation in rng.choice(len(real_flags), size=n_simulated, replace=True,
                                  p=probabilities):
        observations[observation] += 1

    (junction, count), = observations.most_common(1)

    return count > 1


n_simulations_per_scenario = 10000

subsamples = [
    0.237,
    0.154,
    0.1,
    0.0750,
    0.0562,
    0.0421,
    0.0316,
    0.0237,
    0.0154,
    0.01,
    0.00750,
    0.00562,
]

for s1, s2 in zip(subsamples, subsamples[1:]):
    subsamples.append(round(geometric_mean([s1, s2]), ndigits=5))
subsamples.sort(reverse=True)

for s1, s2 in zip(subsamples, subsamples[1:]):
    subsamples.append(round(geometric_mean([s1, s2]), ndigits=5))
subsamples.sort(reverse=True)

for s1, s2 in zip(subsamples, subsamples[1:]):
    subsamples.append(round(geometric_mean([s1, s2]), ndigits=5))
subsamples.sort(reverse=True)




print("subsample", *sorted(flags), sep="\t")
for subsample in subsamples:
    row = [subsample]
    # get all the concentrations into the same range for a faster simuation
    effective_subsample = subsample
    for concentration in sorted(flags):
        t = 0
        for _ in range(n_simulations_per_scenario):
            if simulate_one(concentration, effective_subsample):
                t += 1
        row.append(t / n_simulations_per_scenario)
        effective_subsample = subsample / 10
    print(*row, sep="\t")
