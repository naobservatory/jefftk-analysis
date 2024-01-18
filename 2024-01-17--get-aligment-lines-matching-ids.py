#!/usr/bin/env python3
import sys

target_ids_file, match_column = sys.argv[1:]
match_column = int(match_column)

targets = set()
with open(target_ids_file) as inf:
    for line in inf:
        targets.add(line.split()[0])

for line in sys.stdin:
    if line.split()[match_column] in targets:
        sys.stdout.write(line)
