#!/usr/bin/env python3

import datetime

costs = []

months = [
    "Ignore",
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec",
]

with open("Sequencing Cost - Data Table.tsv") as inf:
    for line in inf:
        line = line.removesuffix("\n")

        bits = line.split("\t")
        date, cost_per_mb, cost_per_genome = bits

        if date == "Date":
            continue

        month_full, yy = date.split("-")

        month_numeric = months.index(month_full)

        yyyy = 2000 + int(yy)

        costs.append((
            datetime.date(year=yyyy,
                          month=month_numeric,
                          day=1),
            float(cost_per_mb.removeprefix("$").replace(",", ""))))

def find_cost(year, month):
    target_date = datetime.date(year=year, month=month, day=1)

    best_cost = None
    best_date = None
    for candidate_date, cost in costs:
        if (best_date is None or
            abs((target_date - candidate_date).days) <
            abs((target_date - best_date).days)):

            best_cost = cost
            best_date = candidate_date

    return best_cost, best_date
        
with open("results-20230720-131657 - results-20230720-131657.tsv") as inf:
    for line in inf:
        line = line.removesuffix("\n")
        bits = line.split("\t")

        year, month, mbases, *_ = bits
        if year == "Year":
            continue
        
        year = int(year)
        month = int(month)
        mbases = int(mbases)

        cost, reference_date = find_cost(year, month)

        print(year, month, mbases, cost, reference_date, sep="\t")
