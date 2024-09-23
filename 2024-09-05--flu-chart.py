#!/usr/bin/env python3

import re
import os
import csv
import gzip
import json
import datetime
from collections import defaultdict, Counter

TARGET_STATES = {"Missouri": "MU",
                 "California": "UCI"}
FLU_A=11320
FLU_B=11520
H5N1=102793
TARGET_TAXIDS = {FLU_A: "Influenza A",
                 FLU_B: "Influenza B"}

POPULATION = {
    "Missouri": 6_178_000,
    "California": 39_030_000,
}

# https://ndc.services.cdc.gov/wp-content/uploads/2021/02/MMWR_Week_overview.pdf
# "The first day of any MMWR week is Sunday. MMWR week numbering is
#  sequential beginning with 1 and incrementing with each week to a
#  maximum of 52 or 53. MMWR week #1 of an MMWR year is the first week
#  of the year that has at least four days in the calendar year. For
#  example, if January 1 occurs on a Sunday, Monday, Tuesday or
#  Wednesday, the calendar week that includes January 1 would be MMWR
#  week #1. If January 1 occurs on a Thursday, Friday, or Saturday, the
#  calendar week that includes January 1 would be the last MMWR week of
#  the previous year (#52 or #53). Because of this rule, December 29,
#  30, and 31 could potentially fall into MMWR week #1 of the following
#  MMWR year."
# Another way of saying this is that in years where Jan 1 falls on Thursday,
# Friday, Saturday, or Sunday then the first MMWR week of the year starts on
# the first Sunday of the year, while if it falls on Monday, Tuesday, or
# Wednesday then the first week starts with the last Sunday of the previous
# year.
def parse_mmwr_week(year: int, week: int) -> datetime.date:
    first_day_of_year = datetime.date(year, 1, 1)
    first_sunday_of_year = datetime.date(
        year, 1, 7 - first_day_of_year.weekday()
    )
    assert first_sunday_of_year.weekday() == 6

    if first_day_of_year.weekday() in [3, 4, 5, 6]:
        # First sunday of the year is included in MMWR week 1.
        return first_sunday_of_year + datetime.timedelta(weeks=week - 1)
    else:
        # First sunday of the year is not included in MMWR week 1.
        return first_sunday_of_year + datetime.timedelta(weeks=week - 2)

def load_weekly_data():
    # state -> taxid -> date -> num positive
    output = defaultdict(lambda: defaultdict(Counter))

    # Downloaded 2023-05-08 from
    # https://gis.cdc.gov/grasp/fluview/fluportaldashboard.html with options:
    #   Select Data Source:
    #      [x] WHO/NREVSS   [ ] ILINet
    #   [x] State
    #   Select Regions: [x] Select All
    #   Select Seasons:
    #     [x] 2023-2024
    #     [x] 2022-2023
    #     [x] 2021-2022
    #     [x] 2020-2021
    #     [x] 2019-2020
    #
    # This is CDC data on the number of positive tests by week, for Flu A vs
    # Flu B.
    #
    # This comes from clinical labs; I can't find week-level data for the
    # public health labs, but there's already a lot here.
    with open("/Users/jeffkaufman/code/p2ra-restricted/prevalence-data/"
              "CDC_WHO_NREVSS_Clinical_Labs.csv") as inf:
        cols = None
        for row in csv.reader(inf):
            if row[0].startswith("*"):
                continue  # initial comment

            if not cols:
                cols = row
                continue

            state = row[cols.index("REGION")]
            year = int(row[cols.index("YEAR")])
            mmwr_week = int(row[cols.index("WEEK")])
            parsed_start = parse_mmwr_week(year, mmwr_week)
            total_tests = row[cols.index("TOTAL SPECIMENS")]
            positive_a = row[cols.index("TOTAL A")]
            positive_b = row[cols.index("TOTAL B")]

            if state not in TARGET_STATES:
                continue

            if total_tests == positive_a == positive_b == "X":
                continue

            output[state][11320][parsed_start] = int(positive_a)
            output[state][11520][parsed_start] = int(positive_b)

    return output

# state -> taxid -> date -> num positive
weekly_data = load_weekly_data()

with open("/Users/jeffkaufman/code/mgs-restricted/dashboard/"
          "metadata_papers.json") as inf:
    metadata_papers = json.load(inf)

with open("/Users/jeffkaufman/code/mgs-restricted/dashboard/"
          "metadata_bioprojects.json") as inf:
    metadata_bioprojects = json.load(inf)

with open("/Users/jeffkaufman/code/mgs-restricted/dashboard/"
          "metadata_samples.json") as inf:
    metadata_samples = json.load(inf)

parents = {}
with open("/Users/jeffkaufman/code/mgs-pipeline/dashboard/nodes.dmp") as inf:
    for line in inf:
        child_taxid, parent_taxid, *_ = \
            line.replace("\t|\n", "").split("\t|\t")
        parents[int(child_taxid)] = int(parent_taxid)

# state -> date -> reads
sample_depth = defaultdict(Counter)

# state -> date -> taxid -> observations
observations = defaultdict(lambda: defaultdict(Counter))

# state -> date -> taxid -> observations
observations_dedup = defaultdict(lambda: defaultdict(Counter))

sample_states = {} # sample -> state
sample_dates = {} # sample -> date
state_dates = defaultdict(set) # state -> dates

seen = {} # min(25bp), max(25bp) -> first read id, first sample
DUP_SUBSET_START = 1
DUP_SUBSET_END = DUP_SUBSET_START + 25
def rc(s):
    return "".join(
        {"T": "A", "G": "C", "A": "T", "C": "G", "N": "N"}[x]
        for x in reversed(s))

DATA_DIR="/Users/jeffkaufman/work/2024-09-07--flu-segment-coverage"
for paper in metadata_papers:
    flu_fname = "%s/%s.flu.tsv" % (DATA_DIR, paper)
    if not os.path.exists(flu_fname):
        continue

    for bioproject in metadata_papers[paper]["projects"]:
        for sample in metadata_bioprojects[bioproject]:
            state = metadata_samples[sample].get("state")
            date = metadata_samples[sample].get("date")

            if (state and date and state in TARGET_STATES and
                len(date) == len("2020-01-01")):

                date = datetime.date.fromisoformat(date)

                sample_depth[state][date] += metadata_samples[
                    sample]["reads"]
                sample_states[sample] = state
                sample_dates[sample] = date
                state_dates[state].add(date)

        with gzip.open(os.path.join(
                "/Users/jeffkaufman/code/p2ra-restricted/bioprojects",
                bioproject, "hv_clade_counts.tsv.gz"), "rt") as inf:
            col = None
            for line in inf:
                row = line.rstrip("\n").split("\t")
                if not col:
                    col = row
                    continue

                sample = row[col.index("sample")]
                taxid = int(row[col.index("taxid")])
                count = int(row[col.index("n_reads_clade")])
                if taxid in TARGET_TAXIDS and sample in sample_states:
                    observations[sample_states[sample]][
                        sample_dates[sample]][taxid] += count

        with open(flu_fname) as inf:
            for line in inf:
                read_id, sample, taxid, genome_id, fwd, rev = \
                    line.rstrip().split("\t")
                assert sample in metadata_samples
                if sample not in sample_states:
                    continue                
                
                taxid = int(taxid)

                key1 = (fwd[DUP_SUBSET_START:DUP_SUBSET_END] + "-" +
                        rc(rev[DUP_SUBSET_START:DUP_SUBSET_END]))
                key2 = (rev[DUP_SUBSET_START:DUP_SUBSET_END] + "-" +
                        rc(fwd[DUP_SUBSET_START:DUP_SUBSET_END]))
                key = min(key1, key2)
                if key in seen:
                    continue # only count unique reads
                else:
                    seen[key] = read_id, sample, fwd, rev

                while taxid not in [0, 1, FLU_A, FLU_B]:
                    taxid = parents[taxid]
                assert taxid in [FLU_A, FLU_B]

                observations_dedup[sample_states[sample]][
                    sample_dates[sample]][taxid] += count
                
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

for taxid in TARGET_TAXIDS:
    fig, axs = plt.subplots(figsize=(8, 4*len(TARGET_STATES)),
                            nrows=len(TARGET_STATES),
                            sharex=True,
                            tight_layout=True)
    for ax, (state, dataset) in zip(fig.axes, sorted(TARGET_STATES.items())):
        color = {
            "California": "#508cf3",
            "Missouri": "#ec5143",
        }[state]

        linestyle = {
            FLU_A: "-",
            FLU_B: "--",
        }[taxid]
        
        ax2 = ax.twinx()

        ax.set_ylabel("Relative abundance")
        ax2.set_ylabel("State-level weekly positive tests per 100k")

        xs = list(sorted(state_dates[state]))
        if not xs:
            continue

        ys = [observations[state][date][taxid] /
              sample_depth[state][date]
              for date in xs]

        ax.scatter(xs, ys, label="Relative\nAbundance", color="black")

        min_date = datetime.date.fromisoformat("2023-10-01")
        max_date = max(xs) + datetime.timedelta(days=6)

        xs = [date for date in weekly_data[state][taxid]
              if min_date <= date <= max_date]
        ys = [weekly_data[state][taxid][date] *
              100_000 / POPULATION[state]
              for date in xs]

        ax2.plot(xs, ys, label="Positive Tests", color=color, linestyle=linestyle)

        study_start = datetime.date.fromisoformat("2023-12-01")
        study_end = datetime.date.fromisoformat("2024-05-01")
        ax.axvspan(study_start, study_end, alpha=0.1, color='gray',
                   label='Analysis\nPeriod')

        ax.legend(loc="upper left")
        ax2.legend(loc='upper right')

        ax.set_title({
            "UCI": "University of California, Irvine",
            "MU": "University of Missouri",
        }[dataset])

        ax2.set_ylim(ymin=0)
        ax.set_ylim(ymin=0)

        ax.xaxis.set_major_locator(mdates.MonthLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=20, ha='right')

    fig.suptitle(TARGET_TAXIDS[taxid])
    fig.savefig("relative-abundance-positive-tests-%s.png" % taxid, dpi=300)
    plt.clf()
