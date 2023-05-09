import os
import gzip
import json
import subprocess

dashboard = os.path.expanduser("~/code/mgs-pipeline/dashboard/")

with open(os.path.join(dashboard, "human_virus_sample_counts.json")) as inf:
    human_virus_sample_counts = json.load(inf)

with open(os.path.join(dashboard, "metadata_samples.json")) as inf:
    metadata_samples = json.load(inf)

with open(os.path.join(dashboard, "metadata_bioprojects.json")) as inf:
    metadata_bioprojects = json.load(inf)

with open(os.path.join(dashboard, "metadata_papers.json")) as inf:
    metadata_papers = json.load(inf)

spurbeck_bioproject, = metadata_papers["Spurbeck 2023"]["projects"]
samples = metadata_bioprojects[spurbeck_bioproject]

target_cc_taxids = {
    694003: "betacoronavirus 1",
    95341: "sapovirus",
    39733: "astroviridae",
    694009: "sars-cov-2",
    142786: "norovirus",
    12234: "tobamovirus",
    2559587: "other riboviria",
    10239:   "other viruses",
}

target_hv_taxids = {
    2559587: "riboviria",
    10239:   "viruses"
}

# sample -> taxid -> count
ccs = {}

for sample in samples:
    cladecounts = "%s.tsv.gz" % sample
    if not os.path.exists(cladecounts):
        subprocess.check_call(
            ["aws", "s3", "cp", "s3://nao-mgs/%s/cladecounts/%s" % (
                spurbeck_bioproject, cladecounts), "."])

    counts = {}
    with gzip.open(cladecounts) as inf:
        for line in inf:
            taxid, _, _, clade_assignments, _ = line.strip().split()
            taxid = int(taxid)
            clade_assignments = int(clade_assignments)

            if taxid in target_cc_taxids:
                counts[taxid] = clade_assignments
    ccs[sample] = counts

import matplotlib.pyplot as plt

sites = set(metadata_samples[sample]["fine_location"]
            for sample in samples)
for site in list(sites):
    if site not in "EFGH":
        sites.remove(site)
markers=("d", "v", "s", "*", "^")

for i, taxid in enumerate(sorted(target_cc_taxids)):
    fig, ax = plt.subplots(constrained_layout=True)
    plt.suptitle("relative abundances by sampling site")
    fig.supylabel("relative abundance of clade")
    fig.supxlabel("date")
    plt.xticks(rotation=45, ha='right')

    ax.set_title(target_cc_taxids[taxid])
    ax.set_yscale("log")
    
    for site_index, site in enumerate(sorted(sites)):
        xs = []
        ys = []
        for sample in samples:
            if metadata_samples[sample]["fine_location"] == site:
                xs.append(metadata_samples[sample]["date"])

                val = ccs[sample].get(taxid, 0)
                if target_cc_taxids[taxid] == "other riboviria":
                    for other_taxid, other_taxid_name in target_cc_taxids.items():
                        if other_taxid_name not in [
                                "other riboviria",
                                "other viruses"]:
                            val -= ccs[sample].get(other_taxid, 0)
                                
                if target_cc_taxids[taxid] == "other viruses":
                    val -= ccs[sample].get(2559587, 0) # riboviria
                
                ys.append(val /
                          metadata_samples[sample]["reads"])
        ax.scatter(xs, ys, label=site, marker=markers[site_index])
    ax.legend()
    
    fig.savefig("efgh-%s.png" % ( target_cc_taxids[taxid].replace(" ", "-")),
                dpi=180)
    plt.clf()
    
