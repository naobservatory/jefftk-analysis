import re
import math
import glob
import scipy
import datetime
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from collections import defaultdict, Counter

# Ignore the 16 out of 94 samples that are really very small.
MIN_VIRAL_COUNT = 30000

# In calibrating, require seeing at least this many reads of a species in every
# sample.
MIN_CALIBRATION = 20

metadata = {}
wtps = set()
with open("/Users/jeffkaufman/code/"
          "mgs-pipeline/studies/PRJNA729801/metadata/metadata.tsv") as inf:
    for line in inf:
        accession, date, wtp = line.strip().split("\t")
        metadata[accession] = date, wtp
        wtps.add(wtp)

wtps = list(sorted(wtps))

lengths = {}
with open("read-lengths.txt") as inf:
    for line in inf:
        accession, length = line.strip().split()
        lengths[accession] = int(length) - 2

# Raw counts from the input files.
virus_to_accession = defaultdict(Counter)
accession_to_virus = defaultdict(Counter)

# Viral relative abundance: raw counts normalized by the total number of
# recognized viral reads.
#
#                                       reads(virus, sample)
#  relative_abundance(virus, sample) =  ----------------------
#                                       reads(virus=*, sample)
#
rel_virus_to_accession = defaultdict(Counter)
rel_accession_to_virus = defaultdict(Counter)

for fname in glob.glob("*.viruscounts.tsv"):
    accession = fname.split(".")[0]

    if accession in ["SRR21452137", "SRR21452135"] and True:
        # These two samples look very different from the rest, and also have
        # accession numbers in a different range. SRR21452137 is especially
        # strange, but both are weird.  Contamination or very different
        # processing?
        continue

    total = 0
    with open(fname) as inf:
        for line in inf:
            count, virus = line.strip().split("\t")
            total += int(count)

    if total < MIN_VIRAL_COUNT and False:
        continue

    with open(fname) as inf:
        for line in inf:
            count, virus = line.strip().split("\t")
            count = int(count)

            if (accession == "SRR21452135" and
                "Cucumber green mottle mosaic virus" in virus and
                False):
                count = 12000

            virus_to_accession[virus][accession] = \
                accession_to_virus[accession][virus] = count

    for virus in accession_to_virus[accession]:
        rel_virus_to_accession[virus][accession] = \
            rel_accession_to_virus[accession][virus] = \
                accession_to_virus[accession][virus] / total

viruses = list(sorted(virus_to_accession))
accessions = list(sorted(accession_to_virus))

z_virus_to_accession = defaultdict(dict)
z_accession_to_virus = defaultdict(dict)
for virus in viruses:
    zscores = scipy.stats.zscore([virus_to_accession[virus][accession]
                                  for accession in accessions])
    for accession, zscore in zip(accessions, zscores):
        z_virus_to_accession[virus][accession] = \
            z_accession_to_virus[accession][virus] = zscore

# For now, calibration viruses are ones that appear 20+ times in every sample
# and we weigh them all equally.  The idea is, if we have 20 observations then
# +/- one would just have us be off by 5%, which isn't that bad.  It feels like
# it should be possible to use every species for calibration, where more
# abundant ones are weighted more heavily to represent us having more
# confidence that we've measured their means correctly, but I'm not all that
# sure how to do it.
#
# Some half-baked thoughts on this.  Consider a vector of read counts:
#
#     a   b    c     d      e
#    [1, 10, 100, 1000, 10000]
#
# We are very confident that the poisson noise introduced by sampling is
# minimal by the time you get to 100 reads, and is already small at 10, but
# pretty significant at 1.  In calibrating we probably want to weigh d and e at
# ~100%, c and very nearly 100%, b at something like 50%, and a at almost
# nothing.  That is, if we were expecting one read and we got 0 or 4 that's not
# very meaningful, but if we were expecting 10k reads and got 1k or 40k that's
# very meaningful.

calibration_viruses = [
    virus
    for virus in viruses
    if all(accession_to_virus[accession][virus] >= MIN_CALIBRATION
           for accession in accessions)]

calibration_averages = {}
for virus in calibration_viruses:
    calibration_averages[virus] = sum(
        virus_to_accession[virus].values()) / len(accessions)

accession_scalars = {} # accession -> calibration scalar
for accession in accessions:
    accession_scalars[accession] = sum(
        virus_to_accession[virus][accession] / calibration_averages[virus]
        for virus in calibration_viruses) / len(calibration_viruses)

if False:
    for accession in accessions:
        print("%.2f%% %s %s" % (
            accession_scalars[accession],
            accession,
            sum(accession_to_virus[accession].values())))

date_wtp_accession = []
wtp_date_accession = []
len_wtp_date_accession = []
for accession in accessions:
    date, wtp = metadata[accession]
    date_wtp_accession.append((date, wtp, accession))
    wtp_date_accession.append((wtp, date, accession))
    len_wtp_date_accession.append((lengths[accession], wtp, date, accession))

date_wtp_accession.sort()
wtp_date_accession.sort()
len_wtp_date_accession.sort()
accessions_by_date = [accession for date, wtp, accession in date_wtp_accession]
accessions_by_wtp = [accession for wtp, date, accession in wtp_date_accession]
accessions_by_len_wtp = [accession for length, wtp, date, accession
                         in len_wtp_date_accession]

# Viral relative relative abundance:
#
#  relative_relative_abundance(virus, sample) =
#
#             relative_abundance(virus, sample)
#        --------------------------------------------
#        average(relative_abundance(virus, sample=*))
#
rel_rel_virus_to_accession = defaultdict(Counter)
rel_rel_accession_to_virus = defaultdict(Counter)

for virus in viruses:
    total_relative_abundance = 0
    for accession in accessions:
        total_relative_abundance += rel_virus_to_accession[virus][accession]
    average_relative_abundance = total_relative_abundance / len(accessions)

    for accession in accessions:
        rel_rel_virus_to_accession[virus][accession] = \
            rel_rel_accession_to_virus[accession][virus] = (
                rel_virus_to_accession[virus][accession] /
                average_relative_abundance)

if False:
    for accession in accession_to_virus:
        print("%s: %s %.2f%% %.2f%%" % (
            accession,
            accession_to_virus[accession][
                "Tomato brown rugose fruit virus (taxid 1761477)"],
            100*rel_accession_to_virus[accession][
                "Tomato brown rugose fruit virus (taxid 1761477)"],
            100*rel_rel_accession_to_virus[accession][
                "Tomato brown rugose fruit virus (taxid 1761477)"]))

# Alternatively, instead of normalizing we could plot cosine similarity (dot
# product).  This doesn't care about the lengths of the vectors (total
# abundances) just how they compare in terms of where they point.

accession_to_vector = {}
for accession in accessions:
    accession_to_vector[accession] = np.array([
        accession_to_virus[accession][virus]
        for virus in viruses])
total_vector = np.add.reduce(list(accession_to_vector.values()))

def cosine_similarity(a, b):
    return np.dot(a, b) / (np.linalg.norm(a)*np.linalg.norm(b))

# [(average similarity, accession), ...]
typicalities = [
    (np.mean([cosine_similarity(accession_to_vector[accession1],
                                accession_to_vector[accession2])
              for accession2 in accessions]),
     accession1)
    for accession1 in accessions]
typicalities.sort()
accessions_by_similarity = [
    accession for similarity, accession in typicalities]

# Start with the most typical accession
greedy_similarity_ordering = [(typicalities[-1][1])]
while len(greedy_similarity_ordering) < len(accessions):
    # add remaining accessions in order of how similar they are to the previous
    best_similarity = None
    best_accession = None
    for accession in accessions:
        if accession not in greedy_similarity_ordering:
            similarity = cosine_similarity(
                accession_to_vector[accession],
                accession_to_vector[greedy_similarity_ordering[-1]])
            if best_similarity is None or similarity > best_similarity:
                best_similarity = similarity
                best_accession = accession
    greedy_similarity_ordering.append(best_accession)


fig, ax = plt.subplots(constrained_layout=True, figsize=(6,12))
labels = [
    "%s %s" % (metadata[accession][1], metadata[accession][0])
    for similarity, accession in typicalities]
widths = [similarity for similarity, accession in typicalities]
ys = np.arange(len(widths))
ax.set_ylabel("all samples, least to most similar")
ax.set_xlabel("mean cosine similarity")
ax.set_title("mean cosine similarity by sample")
ax.barh(y=ys, width=widths, tick_label=labels)
ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize = 8)
fig.savefig("similarities.png", dpi=180)
plt.clf()


for label, ordered_accessions in [
        ("date", accessions_by_date),
        ("wtp", accessions_by_wtp),
        ("len_wtp", accessions_by_len_wtp),
        ("similarity", accessions_by_similarity),
        ("greedy", greedy_similarity_ordering),
        ("accession", accessions)]:
    labels = []
    for accession in ordered_accessions:
        date, wtp = metadata[accession]
        labels.append("%s %s" % (wtp, date))

    data = []
    for accession1 in ordered_accessions:
        row = []
        for accession2 in ordered_accessions:
            row.append(cosine_similarity(
                accession_to_vector[accession1],
                accession_to_vector[accession2]))

        data.append(row)

    fig, ax = plt.subplots(figsize=(12,10))
    df = pd.DataFrame(data, columns=labels, index=labels)

    heatmap = sns.heatmap(df, xticklabels=True, yticklabels=True, ax=ax)
    heatmap.set_xticklabels(heatmap.get_xmajorticklabels(), fontsize = 4)
    heatmap.set_yticklabels(heatmap.get_ymajorticklabels(), fontsize = 4)

    heatmap.set_title("Sample correlations by %s, cosine similarity" % label)
    heatmap.figure.savefig("cs-%s.png" % label, dpi=180)
    plt.clf()

fig, ax = plt.subplots(figsize=(10,10))
for wtp in wtps:
    xs = []
    ys = []
    for accession in accessions:
        date, a_wtp = metadata[accession]
        if a_wtp != wtp: continue

        xs.append(datetime.datetime.fromisoformat(date))
        ys.append(cosine_similarity(
            accession_to_vector[accession], total_vector))

    ax.scatter(xs, ys, label=wtp)
ax.legend()
ax.set_xlabel("sample date")
ax.set_ylabel("sample cosine similarity to average across all samples")
fig.savefig("scatter.png", dpi=180)
plt.clf()

exit(0)

def short_name(virus):
    return re.sub("[(][^)]*[)]", "", virus)

n_cols = 2
n_rows = math.ceil(len(calibration_viruses) / n_cols)

fig, axs = plt.subplots(figsize=(4*n_cols, 3*n_rows),
                        nrows=n_rows,
                        ncols=n_cols,
                        constrained_layout=True)
plt.suptitle("calibration")
fig.supxlabel("absolute abundance")
fig.supylabel("expected abundance")
for i, virus in enumerate(calibration_viruses):
    ax = axs[i % n_rows][i // n_rows]
    ax.set_xscale("log")
    ax.set_yscale("log")
    xs = []
    ys = []
    pmmov_ys = []

    for accession in accessions:
        # Now make a prediction from all the other viruses and samples.
        xs.append(virus_to_accession[virus][accession])

        # Using each other calibration virus in turn, how much would we
        # predict?  We predict as:
        #
        #   prediction = reads(other_virus, this_sample) *
        #       average(reads(this_virus, other_sample) /
        #               reads(other_virus, other_sample)
        #
        # Then we compute an average of those predictions.

        predictions = []
        for other_virus in calibration_viruses:
            if virus == other_virus: continue

            prediction = (
                virus_to_accession[other_virus][accession] *
                np.mean([
                    virus_to_accession[virus][other_accession] /
                    virus_to_accession[other_virus][other_accession]
                    for other_accession in accessions
                    if other_accession != accession]))
            predictions.append(prediction)
            if "Pepper" in other_virus:
                pmmov_ys.append(prediction)

        ys.append(np.mean(predictions))

    ax.scatter(xs, ys, label="calibrated")

    if "Pepper" not in virus:
        ax.scatter(xs, pmmov_ys, label="pmmov")

    ys = xs
    ax.scatter(xs, ys, label="perfect")
    ax.set_title(short_name(virus))
    ax.legend()

fig.savefig("cv.png", dpi=180)
plt.clf()
