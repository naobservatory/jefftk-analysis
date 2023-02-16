#!/usr/bin/env python3

import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

plt.xkcd()
plt.rcParams["font.family"] = "sans"

fig, axs = plt.subplots(constrained_layout=True,
                        ncols=4,
                        nrows=1,
                        figsize=(8, 3))

plt.suptitle("Rough Benefit of Flagging a Pandemic Sooner\n"
             "(illustrative purposes only)")
fig.supxlabel("time")

def b(t, d):
    return 100-sum(d(i) for i in range(t+1))

def d1(t):
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
        return 46
    if t == 11:
        return 52
    if t == 12:
        return 54
    if t == 13:
        return 52
    if t == 14:
        return 48
    if t == 15:
        return 42
    if t == 16:
        return 28
    if t == 17:
        return 26
    if t == 18:
        return 25
    return d(t-1)*.98

    #return b(t-1) - b(t)

def d2(t):
    if t == 0:
        return 1

    if t < 50:
        return d2(t-1)*1.5
    return 0

def d3(t):
    if t == 0:
        return 1

    if t < 30:
        return d3(t-1)*1.5
    return d3(t-1)*.95

def d4(t):
    if t < 90:
        return t
    return 0

twins = [(ax_n, ax_n.twinx())
         for ax_n in axs]

labeled = False
for d, (ax1, ax2) in zip(
        [d4, d2, d3, d1],
        twins):
    xs = []
    ys_b = []
    ys_d = []

    for t in range(100):
        xs.append(t)
        ys_b.append(b(t, d))
        ys_d.append(d(t))

    ax1.set_xticks([])
    ax1.set_yticks([])
    ax2.set_yticks([])

    if labeled:
        ax2.plot(xs, ys_b, color="blue")
        ax1.plot(xs, ys_d, color="orange")
    else:
        ax2.plot(xs, ys_b, label="b(t)", color="blue")
        ax1.plot(xs, ys_d, label="d(t)", color="orange")
        labeled = True

    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

fig.legend()

fig.savefig("delay.png", dpi=180)
