#!/usr/bin/env python3

import sys
import math
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np

wtp, fname_in, fname_prefix = sys.argv[1:]

X_SCALAR=100
Y_SCALAR=1

min_x = max_x = min_y = max_y = None

with open(fname_in) as inf:
    for line in inf:
        _, log_mean, normalized_variance = line.strip().split()
        log_mean = float(log_mean)
        normalized_variance = float(normalized_variance)
        if min_x is None or log_mean < min_x:
            min_x = log_mean
        if max_x is None or log_mean > max_x:
            max_x = log_mean
        if min_y is None or normalized_variance < min_y:
            min_y = normalized_variance
        if max_y is None or normalized_variance > max_y:
            max_y = normalized_variance
                        
def to_xpip(val):
    return round((val- min_x)*X_SCALAR)
def to_ypip(val):
    return round((val - min_y)*Y_SCALAR)

def xpip_to_real(xpip):
    xval = (xpip + min_x*X_SCALAR) / X_SCALAR
    return 10**xval  # log(mean, 10) -> mean

def ypip_to_real(ypip):
    yval = (ypip + min_y*Y_SCALAR) / Y_SCALAR
    return yval

data = []
for y_val in range(to_ypip(min_y), to_ypip(max_y)+1):
    row = []
    for x_val in range(to_xpip(min_x), to_xpip(max_x)+1):
        row.append(0)
    data.append(row)

with open(fname_in) as inf:
    for line in inf:
        count, log_mean, normalized_variance = line.strip().split()
        count = int(count)
        log_mean = float(log_mean)
        normalized_variance = float(normalized_variance)

        try:
            data[to_ypip(normalized_variance)][to_xpip(log_mean)] += count
        except IndexError:
            print(log_mean, to_xpip(log_mean), min_x,max_x)
            print(normalized_variance, to_ypip(normalized_variance), min_y,max_y)
            print(len(data))
            print(len(data[0]))
            raise

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
ax.set_ylabel("variance / mean abundance")
ax.set_xlabel("mean abundance")
ax.set_title("%s per-kmer normalized variance as a function of abundance" % wtp)
ax.set_xscale("log")
fig.savefig(fname_prefix + "-log.png", dpi=180)

ax.set_ylim(ymax=200)
fig.savefig(fname_prefix + "-log-trunc.png", dpi=180)

ax.set_xscale("linear")
ax.set_xlim(xmax=600)
ax.set_ylim(ymax=600)
fig.savefig(fname_prefix + "-lin-trunc.png", dpi=180)
