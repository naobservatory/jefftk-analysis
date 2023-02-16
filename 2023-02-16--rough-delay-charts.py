#!/usr/bin/env python3

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

plt.xkcd()
plt.rcParams["font.family"] = "serif"

fig, ax = plt.subplots(constrained_layout=True)

ax.set_title("Rough Benefit of Catching a Pandemic Sooner\n"
             "For illustrative purposes only")
ax.set_xlabel("time")

def b(t):
    return 100-sum(d(i) for i in range(t+1))

def d(t):
    if t == 0:
        return 0
    if t == 1:
        return 0.1
    if t == 2:
        return 0.3
    if t == 3:
        return 0.5
    if t == 4:
        return 1
    if t == 5:
        return 2
    if t == 6:
        return 4
    if t == 7:
        return 8
    if t == 8:
        return 16
    if t == 9:
        return 32
    if t == 10:
        return 48
    if t == 11:
        return 52
    if t == 12:
        return 40
    if t == 13:
        return 30
    if t == 14:
        return 25
    if t == 15:
        return 23
    return d(t-1)*.95

    #return b(t-1) - b(t)

xs = []
ys_b = []
ys_d = []

for t in range(100):
    xs.append(t)
    ys_b.append(b(t))
    ys_d.append(d(t))

ax1 = ax
ax2 = ax1.twinx()

ax1.set_xticks([])
ax1.set_yticks([])
ax2.set_yticks([])

ax1.plot(xs, ys_b, label="b(t): benefit of catching at this stage", color="blue")
ax2.plot(xs, ys_d, label="d(t): cost of delaying one day, b'(t)", color="orange")

ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

fig.legend(loc="center right")
    
fig.savefig("delay.png", dpi=180)
