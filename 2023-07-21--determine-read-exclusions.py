#!/usr/bin/env python3

import re
import sys
import glob
import gzip
import json
import math
from collections import defaultdict, Counter

# output of p2ra/list_taxids.py
taxid_list="""
10298           hsv_1           HSV-1
10310           hsv_2           HSV-2
10359           cmv             CMV
10376           ebv             EBV
10407           hbv             HBV
10566           hpv             HPV
10632           jcv             JCV
10804           aav2            AAV2
11103           hcv             HCV
11320           influenza       Influenza A
11520           influenza       Influenza B
11676           hiv             HIV
68558           aav6            AAV6
82300           aav5            AAV5
122928          norovirus       Norovirus (GI)
122929          norovirus       Norovirus (GII)
493803          mcv             MCV
1891762         bkv             BKV
2697049         sars_cov_2      SARS-COV-2
"""

taxid_to_name = {}
for line in taxid_list.split("\n"):
  line = line.strip()
  if not line: continue

  taxid, human_readable = re.match(
    "^([0-9]+) +[a-z0-9_]+ +(.*)$", line).groups()
  taxid_to_name[int(taxid)] = human_readable

def accept(taxid, percent_identical):
  if taxid == 10804: # AAV2
    #  Checked 9 reads, including the worst ones.  All look real, and
    #  include sizeable chunks that are great matches for AAV2 and not for
    #  anything else.  The low-scoring AAV2 reads often have long
    #  suspicious-looking runs of G with some C mixed in.
    assert percent_identical >= 37
    return True
  elif taxid == 82300: # AAV2
    # Only one read
    # No good matches; all blastn matches are bacterial
    assert percent_identical == 40
    return False
  elif taxid in [1891762, 10632, 493803]: # BKV, JCV, MCV
    # Checked all the worst-ranked ones with blast, and they still look
    # good.  Probably lots of variety that isn't captured by the RefSeq
    # genome?
    return True
  elif taxid == 10359: # CMV
    if percent_identical >= 86:
      return True
    elif percent_identical <= 47:
      # Checking blast these all look bacterial.
      return False
    else:
      assert False
  elif taxid == 10376: # EBV
    if percent_identical >= 69:
      return True
    elif percent_identical == 47:
      # Only one read didn't look like EBV; probably bacterial
      return False
    else:
      assert False
  elif taxid == 11103: # HCV
    # Only one read, looks bacterial
    assert percent_identical == 40
    return False
  elif taxid == 11676: # HIV
    assert percent_identical >= 75
    # All three reads looked good on manual inspection
    # The only reason one got 75% was that adapter trimming didn't work
    # on this read.
    return True
  elif taxid == 10566: # HPV
    assert percent_identical == 100
    return True
  elif taxid == 10298: # HSV-1
    if percent_identical >= 69:
      return True
    elif percent_identical <= 51:
      # Two of these, look bacterial
      return False
    else:
      assert False
  elif taxid == 11320: # Flu A
    # Only one read
    # Some kind of bacteria per megablast
    assert percent_identical == 37
    return False
  elif taxid in [122928, 122929]: # Norovirus
    # Manually checked and all low-ranked ones look good
    return True
  elif taxid == 2697049: # SARS-CoV-2
    if percent_identical >= 64:
      return True
    elif percent_identical == 57:
      # Better match for bacterial than SARS-CoV-2, though could be chimeric
      return False
    else:
      assert False

  assert False          

exclusions = []
for validation_summary in glob.glob("*.validation-summary.tsv"):
  taxid = int(validation_summary.removesuffix(".validation-summary.tsv"))
  with open(validation_summary) as inf:
    for line in inf:
      line = line.removesuffix("\n")
      paper, read_id, genome_index, read1_desc, read2_desc, read1, read2 = \
        line.split("\t")

      matches1 = int(read1_desc.split(":")[-1].split("/")[0])
      matches2 = int(read2_desc.split(":")[-1].split("/")[0])

      total1 = int(read1_desc.split(" ")[0].split("/")[-1])
      total2 = int(read2_desc.split(" ")[0].split("/")[-1])

      percent_identical = round(100*(matches1 + matches2) / (total1 + total2))

      if not accept(taxid, percent_identical):
        exclusions.append((read_id, taxid, taxid_to_name[taxid]))

with open("excluded-read-ids.tsv", "w") as outf:
  for exclusion in sorted(exclusions):
    outf.write("%s\t%s\t%s\n" % exclusion)
        
