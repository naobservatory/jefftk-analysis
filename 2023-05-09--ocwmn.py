#!/usr/bin/env python3

import json

site_populations = {
    "A": 365_000,
    "B":  54_000,
    "C": 323_000,
    "D":  11_000,
    "E":  25_000,
    "F": 650_000,
    "G":  45_000,
    "H": 655_000,
    "I":  46_000,
    "J":  65_000,
}

# site -> delta -> [name, county]
distances = {}
for site in site_populations:
    distances[site] = {}

# From https://storymaps.arcgis.com/stories/0fc051c015dc4e8b9a80b76cdf406b24 with
# https://www.arcgis.com/sharing/rest/content/items/5967d1c921424060a47241fa98a286a0/data?f=json -O ocwmn.json
with open("ocwmn.json") as inf:
    for feature in json.load(inf)["operationalLayers"][2][
            "featureCollection"]["layers"][0]["featureSet"]["features"]:
        pop = feature["attributes"]["Population_Served"]
        name = feature["attributes"]["Facility_Name"]
        county = feature["attributes"]["County"]

        for site in site_populations:
            distances[site][abs(site_populations[site] - pop)] = (
                name, county)

for site in site_populations:
    print("%s:" % site)
    for n, (distance, (name, county)) in enumerate(
            sorted(distances[site].items())):
        if n > 5: break
        print("  %s %s; %s" % (str(distance).rjust(10), name, county))
        
