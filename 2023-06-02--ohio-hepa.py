#!/usr/bin/env python3
import sys

sys.path.append("/Users/jeffkaufman/code/p2ra/")
import populations

with open("ohio-hepa.tsv") as inf:
    for line in inf:
        line = line.strip()
        if not line:
            continue
        county_short, cases = line.split("\t")
        cases = int(cases)

        county_long = "%s County" % county_short

        pop = populations.us_population(
            state="Ohio", county=county_long, year=2020)

        cases_per_100k_y = cases / pop.people * 100_000 / 4.5

        print ("%.2f per 100k/y (n=%s) %s" % (
            cases_per_100k_y, cases, county_short))
