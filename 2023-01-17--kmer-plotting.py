#!/usr/bin/env python3

import sys
import math
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np

wtp, fname_in, fname_out, xmax, ymax = sys.argv[1:]

xmax = int(xmax)
ymax = int(ymax)

min_x = max_x = min_y = max_y = None

all_x = set()
all_y = set()

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

        all_x.add(mean)
        all_y.add(variance)

min_x = min_y = 0
        
all_x = list(sorted(all_x))
all_y = list(sorted(all_y))

incr_x = min((a-b) for (a,b) in zip(all_x[1:], all_x))
incr_y = min((a-b) for (a,b) in zip(all_y[1:], all_y))
            
if max_x > xmax:
    max_x = xmax
if max_y > ymax:
    max_y = ymax

xpips = round((max_x - min_x + 1) / incr_x)
ypips = round((max_y - min_y + 1) / incr_y)

while xpips > 1000:
    xpips //= 2

while ypips > 1000:
    ypips //= 2

X_SCALAR = xpips / (max_x - min_x + 1)
Y_SCALAR = ypips / (max_y - min_y + 1)

def to_xpip(val):
    return int((val - min_x)*X_SCALAR)
def to_ypip(val):
    return int((val - min_y)*Y_SCALAR)

def xpip_to_real(xpip):
    return (xpip + min_x*X_SCALAR) / X_SCALAR

def ypip_to_real(ypip):
    return (ypip + min_y*Y_SCALAR) / Y_SCALAR

data = []
for y_val in range(ypips):
    row = []
    for x_val in range(xpips):
        row.append(0)
    data.append(row)

print(len(data[0]), len(data))
    
with open(fname_in) as inf:
    for line in inf:
        count, mean, variance = line.strip().split()
        count = int(count)
        mean = float(mean)
        variance = float(variance)

        xpip = to_xpip(mean)
        ypip = to_ypip(variance)
        
        if xpip >= xpips or ypip >= ypips:
            continue
        
        try:
            data[ypip][xpip] += count
        except IndexError:
            print(mean, to_xpip(mean), min_x,max_x)
            print(variance, to_ypip(variance), min_y,max_y)
            print(len(data))
            print(len(data[0]))
            raise

if True:
    for y_pip in range(ypips):
        for x_pip in range(xpips):
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
ax.set_xlim(xmax=xmax)
ax.set_ylim(ymax=ymax)
fig.savefig(fname_out, dpi=180)
