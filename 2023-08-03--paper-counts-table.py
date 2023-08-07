#!/usr/bin/env python3

import json
with open("metadata_papers.json") as inf:
    papers = json.load(inf)

with open("metadata_bioprojects.json") as inf:
    bioprojects = json.load(inf)

with open("metadata_samples.json") as inf:
    samples = json.load(inf)

with open("human_virus_sample_counts.json") as inf:
    human_virus_sample_counts = json.load(inf)

with open("comparison_sample_counts.json") as inf:
    comparison_sample_counts = json.load(inf)

def determine_enrichment(sample_attrs, paper_name):
    return samples[sample].get("enrichment", "")

VIRUSES="10239"
print("paper\tcountry\treads\thuman viral relative abundance")
for paper_name, paper_attrs in sorted(papers.items()):
    if paper_attrs["link"] == "personal communication": continue

    def include(sample):
        if samples[sample].get("enrichment", None) == "panel":
            return False
        
        if paper_name == "Bengtsson-Palme 2016":
            return samples[sample]["fine_location"].startswith("Inlet")

        if paper_name == "Ng 2019":
            return samples[sample]['fine_location'] == "Influent"

        return True

    reads = 0
    viral_reads = 0
    human_viral_reads = 0
    countries = set()

    for bioproject in paper_attrs["projects"]:
        for sample in bioprojects[bioproject]:
            if not include(sample): continue

            countries.add(samples[sample]["country"])
            
            reads += samples[sample]["reads"]
            for human_virus, sample_counts in human_virus_sample_counts.items():
                human_viral_reads += sample_counts.get(sample, 0)

            viral_reads += comparison_sample_counts[VIRUSES].get(sample, 0)

    if not countries:
        print(paper_name)

    name = paper_name
    print("%s\t%s\t%s\t%.10f" % (
        name,
        "Multiple" if len(countries) > 1 else next(iter(countries)),
        reads,
        human_viral_reads / reads,
    ))
