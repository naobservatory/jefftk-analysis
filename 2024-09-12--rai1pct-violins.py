#!/usr/bin/env python3

import numpy as np

rothman_sars_cov_2 = []
uci_flu_a = []
uci_flu_b = []
mu_flu_a = []
mu_flu_b = []

SARS_COV_2="2697049"
FLU_A="11320"
FLU_B="11520"

with open("/Users/jeffkaufman/code/p2ra-restricted/"
          "model_output/fits.tsv") as inf:
    col = None
    for line in inf:
        row = line.rstrip("\n").split("\t")
        if not col:
            col = row
            continue

        def v(x):
            return row[col.index(x)]
        
        if v("location") != "Overall":
            continue

        rai100 = float(v("ra_at_1in100"))
        
        if v("study") == "rothman" and v("taxids") == SARS_COV_2:
            rothman_sars_cov_2.append(rai100)
        elif v("study") == "rothman-nao" and v("taxids") == FLU_A:
            uci_flu_a.append(rai100)
        elif v("study") == "rothman-nao" and v("taxids") == FLU_B:
            uci_flu_b.append(rai100)
        elif v("study") == "johnson-nao" and v("taxids") == FLU_A:
            mu_flu_a.append(rai100)
        elif v("study") == "johnson-nao" and v("taxids") == FLU_B:
            mu_flu_b.append(rai100)

data = [uci_flu_b, uci_flu_a, mu_flu_b, mu_flu_a, rothman_sars_cov_2]
for x in data:
    x.sort()
            
def pct(p, d):
    return d[round(p / 100 * len(d))]

print("percentile",
      "ScV2",
      "MU Flu A",
      "MU Flu B",
      "UCI Flu A",
      "UCI Flu B",
      sep="\t")
for percentile in [5, 25, 50, 75, 95]:
    print("%s%%\t%.1e\t%.1e\t%.1e\t%.1e\t%.1e" % (
        percentile,
        pct(percentile, rothman_sars_cov_2),
        pct(percentile, mu_flu_a),
        pct(percentile, mu_flu_b),
        pct(percentile, uci_flu_a),
        pct(percentile, uci_flu_b)))
            
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(8,4),
                       tight_layout=True)

labels = [
    'Unpublished UCI Data\nFlu B\n  6 read pairs',
    'Unpublished UCI Data\nFlu A\n  81 read pairs',
    'Unpublished MU Data\nFlu B\n  162 read pairs',
    'Unpublished MU Data\nFlu A\n  288 read pairs',
    'Rothman et al. (2021)\nSARS-CoV-2\n  384 read pairs',
]

log_data = [np.log10(d) for d in data]
parts = ax.violinplot(log_data, positions=range(len(data)), vert=False,
                      showextrema=False, showmedians=True)

parts["cmedians"].set_edgecolor('black')

for part, label in zip(parts['bodies'], labels):
    if label.startswith("Rothman"):
        part.set_facecolor('red')
    elif label.startswith("Unpublished MU"):
        part.set_facecolor('blue')
    elif label.startswith("Unpublished UCI"):
        part.set_facecolor('green')
        

# Customize the plot
ax.set_yticks(range(len(labels) + 1))
ax.set_yticklabels(labels + [""])
ax.set_title('RAi(1%): expected relative abundance at 1% weekly incidence')

log_ticks = [-11, -10, -9, -8, -7, -6]
ax.set_xlim([min(log_ticks), max(log_ticks)])
ax.set_xticks(log_ticks)
ax.set_xticklabels([f'$10^{{{x}}}$' for x in log_ticks])

PRELIMINARY_ESTIMATE_FLU_A=3e-8
ax.axvline(x=np.log10(PRELIMINARY_ESTIMATE_FLU_A), color='lightgrey',
           linestyle='--', linewidth=1,
           ymin=0,
           ymax=0.93,
           )
ax.text(np.log10(PRELIMINARY_ESTIMATE_FLU_A) + 0.03,
        ax.get_ylim()[1] - 0.2,
        'Previous Flu A Estimate',
        ha='center', va='top', fontsize=8)

plt.savefig("rai1pct-violins.png", dpi=300)
plt.clf()
