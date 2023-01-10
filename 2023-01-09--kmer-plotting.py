#!/usr/bin/env python3

import sys
import math
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np

wtp, fname_in, fname_out_log, fname_out_lin = sys.argv[1:]

MIN_X = -1.05
MAX_X = 4.02

MIN_Y = 0
MAX_Y = 3.17

def to_xpip(val):
    return round(val*100 - MIN_X*100)
def to_ypip(val):
    return round(val*100 - MIN_Y*100)

def xpip_to_real(xpip):
    xval = (xpip + MIN_X*100) / 100
    return 10**xval  # log(mean, 10) -> mean

def ypip_to_real(ypip):
    yval = (ypip + MIN_Y*100) / 100
    return yval

data = []
for y_val in range(to_ypip(MIN_Y), to_ypip(MAX_Y)):
    row = []
    for x_val in range(to_xpip(MIN_X), to_xpip(MAX_X)):
        row.append(0)
    data.append(row)

with open(fname_in) as inf:
    for line in inf:
        count, log_mean, normalized_variance = line.strip().split()
        count = int(count)
        log_mean = float(log_mean)
        normalized_variance = float(normalized_variance)

        data[to_ypip(normalized_variance)][to_xpip(log_mean)] += count

for y_pip in range(len(data)):
    for x_pip in range(len(data[y_pip])):
        val = data[y_pip][x_pip]
        if val > 0:
            val = math.log(val, 10) + 1
        data[y_pip][x_pip] = val


x = [xpip_to_real(xpip) for xpip in range(len(data[0]))]
y = [ypip_to_real(ypip) for ypip in range(len(data))]
fig, ax = plt.subplots()
ax.pcolormesh(x, y, data)
ax.set_ylabel("normalized variance (rms error / mean abundance)")
ax.set_xlabel("mean abundance")
ax.set_title("%s per-kmer normalized variance as a function of abundance" % wtp)
ax.set_xscale("log")
ax.set_xlim(xmin=5)
fig.savefig(fname_out_log, dpi=600)

ax.set_xscale("linear")
ax.set_xlim(xmin=MIN_X,xmax=20)
fig.savefig(fname_out_lin, dpi=600)
