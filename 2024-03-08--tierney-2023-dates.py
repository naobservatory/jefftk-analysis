#!/usr/bin/env python3
import json
from collections import Counter

with open("dashboard/metadata_samples.json") as inf:
    metadata_samples = json.load(inf)

with open("dashboard/metadata_bioprojects.json") as inf:
    metadata_bioprojects = json.load(inf)

with open("dashboard/metadata_papers.json") as inf:
    metadata_papers = json.load(inf)

bioproject, = metadata_papers["Tierney 2023"]["projects"]
print(bioproject)
reads_by_month = Counter()

for sample in metadata_bioprojects[bioproject]:
    date = metadata_samples[sample]["date"]
    if len(date) == 4:
        continue

    month = date[:7]

    reads_by_month[month] += metadata_samples[sample]["reads"]

for month, count in sorted(reads_by_month.items()):
    print(month, count, sep="\t")




