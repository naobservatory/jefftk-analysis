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
    enrichment = samples[sample].get("enrichment", "")
    method = samples[sample].get("method", "")

    if paper_name in ["Ng 2019", "Wang 2022"]:
        return "retentate"
    elif paper_name in ["Bengtsson-Palme 2016",
                        "Brinch 2020",
                        "Munk 2022",
                        "Petersen 2015",
                        "Fierer 2022",
                        "Hendriksen 2019",
                        ]:
        return "pellet"
    elif paper_name in ["Cui 2023"]:
        return "supernatant"

    elif paper_name in ["Rothman 2021"]:
        if enrichment == "panel":
            return "panel"
        else:
            return "filtrate"

    elif paper_name in ["Crits-Christoph 2021"]:
        prep = {
            "amicon": "amicon",
            "COL": "col",
            "MOS": "mos"
        }[method]
        if enrichment == "panel":
            prep = "panel-" + prep
        return prep

    elif paper_name in ["McCall 2023"]:
        if enrichment == "panel":
            return "panel"
        else:
            return "filtrate"

    elif paper_name in ["Yang 2020"]:
        return "viral"

    elif paper_name in ["Spurbeck 2023"]:
        return method

    elif paper_name in ["Johnson 2023"]:
        return "filtrate"
        
    return enrichment

VIRUSES="10239"
print("paper\tenrichment\treads\tviral relative abundance\thuman viral relative abundance")
for paper_name, paper_attrs in sorted(papers.items()):
    #if paper_attrs["link"] == "personal communication": continue

    def include(sample):
        if paper_name == "Bengtsson-Palme 2016":
            return samples[sample]["fine_location"].startswith("Inlet")

        if paper_name == "Ng 2019":
            return samples[sample]['fine_location'] == "Influent"
        
        return True

    enrichments = {}
    for bioproject in paper_attrs["projects"]:
        for sample in bioprojects[bioproject]:
            if not include(sample): continue

            enrichment = determine_enrichment(samples[sample], paper_name)
            enrichments[enrichment] = {
                "reads": 0,
                "viral_reads": 0,
                "human_viral_reads": 0,
            }

    for bioproject in paper_attrs["projects"]:
        for sample in bioprojects[bioproject]:
            if not include(sample): continue

            enrichment = determine_enrichment(samples[sample], paper_name)
            enrichments[enrichment]["reads"] += samples[sample]["reads"]
            for human_virus, sample_counts in human_virus_sample_counts.items():
                enrichments[enrichment][
                    "human_viral_reads"] += sample_counts.get(sample, 0)

            enrichments[enrichment][
                "viral_reads"] += comparison_sample_counts[VIRUSES].get(sample, 0)

    for enrichment in sorted(enrichments):
        reads = enrichments[enrichment]["reads"]
        viral_reads = enrichments[enrichment]["viral_reads"]
        human_viral_reads = enrichments[enrichment]["human_viral_reads"]

        print("%s\t%s\t%s\t%.10f\t%.10f" % (
            paper_name,
            enrichment,
            reads,
            viral_reads / reads,
            human_viral_reads / reads,
        ))
