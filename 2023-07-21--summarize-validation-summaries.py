#!/usr/bin/env python3

import re
import sys
import glob
import gzip
import json
import math
from collections import defaultdict, Counter

COLOR_END = '\x1b[0m'
COLOR_RED = '\x1b[1;31m'
COLOR_BLUE = '\x1b[1;34m'
COLOR_MAGENTA = '\x1b[1;35m'

# output of p2ra/list_taxids.py
taxid_list="""
10298           hsv_1           HSV-1
10310           hsv_2           HSV-2
10359           cmv             CMV
10376           ebv             EBV
10407           hbv             HBV
10566           hpv             HPV
10632           jcv             JCV
10804           aav2            AAV2
11103           hcv             HCV
11320           influenza       Influenza A
11520           influenza       Influenza B
11676           hiv             HIV
68558           aav6            AAV6
82300           aav5            AAV5
122928          norovirus       Norovirus (GI)
122929          norovirus       Norovirus (GII)
493803          mcv             MCV
1891762         bkv             BKV
2697049         sars_cov_2      SARS-COV-2
"""

taxid_to_name = {}
for line in taxid_list.split("\n"):
  line = line.strip()
  if not line: continue

  taxid, human_readable = re.match(
    "^([0-9]+) +[a-z0-9_]+ +(.*)$", line).groups()
  taxid_to_name[int(taxid)] = human_readable

# ourtaxid, percent_identical
scores = []
taxid_counts = Counter()
for validation_summary in glob.glob("*.validation-summary.tsv"):
  taxid = int(validation_summary.removesuffix(".validation-summary.tsv"))
  with open(validation_summary) as inf:
    for line in inf:
      line = line.removesuffix("\n")
      paper, read_id, genome_index, read1_desc, read2_desc, read1, read2 = \
        line.split("\t")

      matches1 = int(read1_desc.split(":")[-1].split("/")[0])
      matches2 = int(read2_desc.split(":")[-1].split("/")[0])

      total1 = int(read1_desc.split(" ")[0].split("/")[-1])
      total2 = int(read2_desc.split(" ")[0].split("/")[-1])

      percent_identical = (matches1 + matches2) / (total1 + total2)
      #print("%0.5f" % percent_identical)

      scores.append((percent_identical, taxid))
      taxid_counts[taxid] += 1

all_taxids = [taxid for taxid in taxid_counts if taxid_counts[taxid] > 0]

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

fig, ax = plt.subplots(constrained_layout=True)
ax.xaxis.set_major_formatter(mtick.PercentFormatter())
ax.yaxis.set_major_formatter(mtick.PercentFormatter())

xs = []
ys = []
for n, (percent_identical, _) in enumerate(sorted(scores)):
  xs.append(100*percent_identical)
  ys.append(100*n/len(scores))
ax.plot(xs, ys)

fig.savefig("percent-identical-cdf.png", dpi=180)
plt.clf()

ncols = 3
nrows = math.ceil(len(all_taxids)/3)
fig, axs = plt.subplots(constrained_layout=True,
                        ncols=ncols,
                        nrows=nrows,
                        figsize=(2.5*ncols, 2.5*nrows),
                        sharex=True,
                        sharey=True)
plt.suptitle("Quality of Initial Matches")
fig.supxlabel("percentage of bases matching reference genome")
fig.supylabel("percentage of matches with at least this quality")

named_taxids = [
  (taxid_to_name[taxid], taxid) for taxid in all_taxids]

def collapse(xs_in, ys_in):
  xs_out = []
  ys_out = []
  last_x = None
  for x, y in zip(reversed(xs_in), reversed(ys_in)):
    if last_x == x:
      continue
    last_x = x
    xs_out.append(x)
    ys_out.append(y)

  return list(reversed(xs_out)), list(reversed(ys_out))

for i, (name, taxid) in enumerate(sorted(named_taxids)):
  ax = axs[i//ncols][i%ncols]

  ax.xaxis.set_major_formatter(mtick.PercentFormatter())
  ax.yaxis.set_major_formatter(mtick.PercentFormatter())

  xs = []
  ys = []
  matching = [percent_identical
              for (percent_identical, record_taxid)
              in scores
              if record_taxid == taxid]
  assert matching

  for n, percent_identical in enumerate(sorted(matching)):
    xs.append(100*percent_identical)
    ys.append(100*(n+1)/(len(matching)))
  ax.set_title("%s\n(%s initial match%s)" % (
    name.replace("(", "").replace(")", ""),
    len(matching),
    "es" if len(matching) > 1 else ""))

  xs, ys = collapse(xs, ys)
  
  if len(xs) > 1:
    ax.plot(xs, ys)
  else:
    ax.scatter(xs, ys)

fig.savefig("percent-identical-multi-cdf.png", dpi=180)
plt.clf()
