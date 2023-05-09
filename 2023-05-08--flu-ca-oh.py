#!/usr/bin/env python3

import csv
from collections import Counter

weekly_n = Counter()
weekly_ca = Counter()
weekly_oh = Counter()

with open("CDC_WHO_NREVSS_Clinical_Labs.csv") as inf:
    cols = None
    for row in csv.reader(inf):
        if row[0].startswith("*"):
            continue  # initial comment

        if not cols:
            cols = row
            continue

        region = row[cols.index("REGION")]
        year = int(row[cols.index("YEAR")])
        mmwr_week = int(row[cols.index("WEEK")])
        total_tests = row[cols.index("TOTAL SPECIMENS")]
        positive_a = row[cols.index("TOTAL A")]
        positive_b = row[cols.index("TOTAL B")]

        date = "%s-%s" % (year, str(mmwr_week).zfill(2))
        
        if total_tests == positive_a == positive_b == "X":
            continue

        weekly_n[date] += int(positive_a) + int(positive_b)
        if region == "California":
            weekly_ca[date] += int(positive_a) + int(positive_b)
        elif region == "Ohio":
            weekly_oh[date] += int(positive_a) + int(positive_b)


for date in weekly_n:
    print("%s\t%s\t%s\t%s" % (
        date, weekly_n[date], weekly_ca[date], weekly_oh[date]))
        
