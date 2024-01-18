#!/usr/bin/env python3
import glob
import json
from collections import Counter

result = Counter()
for fname in glob.glob("*.stats.json"):
    with open(fname) as inf:
        for stat, value in json.load(inf).items():
            result[stat] += value

import pprint
pprint.pprint(result)
