#!/usr/bin/env python3

import os
import glob
import json
import pprint
from collections import defaultdict

hiv_hep = defaultdict(list)
for hvr in glob.glob(
        "/Users/jeffkaufman/code/mgs-pipeline/dashboard/hvreads/*.hvreads.json"
) + glob.glob(
    "/Users/jeffkaufman/code/mgs-restricted/dashboard/hvreads/*.hvreads.json"):
    sample = os.path.basename(hvr).removesuffix(".hvreads.json")
    with open(hvr) as inf:
        for read_id, raw_details in json.load(inf).items():
            details = raw_details[:]
            if type(details[0]) == type(1):
                assignment = details.pop(0)
            kraken_info = details.pop(0)
            taxids = set()
            for token in kraken_info.split():
                taxid, count = token.split(":")

                if taxid.isdigit():
                    taxids.add(int(taxid))

            if 11676 in taxids and 35269 in taxids:
                raw_details.insert(0, read_id)
                hiv_hep[sample].append(raw_details)

with open("hiv-hepb.json", "w") as outf:
    json.dump(hiv_hep, outf)
