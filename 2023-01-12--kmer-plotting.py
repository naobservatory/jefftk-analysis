#!/usr/bin/env python3

import sys
import math
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np

wtp, fname_in, fname_prefix = sys.argv[1:]

X_SCALAR=1
Y_SCALAR=0.001

min_x = max_x = min_y = max_y = None

with open(fname_in) as inf:
    for line in inf:
        _, mean, variance = line.strip().split()
        mean = float(mean)
        variance = float(variance)
        if min_x is None or mean < min_x:
            min_x = mean
        if max_x is None or mean > max_x:
            max_x = mean
        if min_y is None or variance < min_y:
            min_y = variance
        if max_y is None or variance > max_y:
            max_y = variance

def to_xpip(val):
    return round((val - min_x)*X_SCALAR)
def to_ypip(val):
    return round((val - min_y)*Y_SCALAR)

def xpip_to_real(xpip):
    return (xpip + min_x*X_SCALAR) / X_SCALAR

def ypip_to_real(ypip):
    return (ypip + min_y*Y_SCALAR) / Y_SCALAR

data = []
for y_val in range(to_ypip(min_y), to_ypip(max_y)+1):
    row = []
    for x_val in range(to_xpip(min_x), to_xpip(max_x)+1):
        row.append(0)
    data.append(row)

with open(fname_in) as inf:
    for line in inf:
        count, mean, variance = line.strip().split()
        count = int(count)
        mean = float(mean)
        variance = float(variance)

        try:
            data[to_ypip(variance)][to_xpip(mean)] += count
        except IndexError:
            print(mean, to_xpip(mean), min_x,max_x)
            print(variance, to_ypip(variance), min_y,max_y)
            print(len(data))
            print(len(data[0]))
            raise

if True:
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
ax.set_ylabel("variance")
ax.set_xlabel("mean abundance")
ax.set_title("%s per-kmer normalized variance as a function of abundance" % wtp)
fig.savefig(fname_prefix + "-lin-full.png", dpi=180)

ax.set_xlim(xmax=250)
ax.set_ylim(ymax=1e4)
fig.savefig(fname_prefix + "-lin-trunc.png", dpi=180)

#    ax.set_xscale("linear")
#    ax.set_xlim(xmax=600)
#    ax.set_ylim(ymax=600)
#    fig.savefig(fname_prefix + "-lin-trunc.png", dpi=180)
