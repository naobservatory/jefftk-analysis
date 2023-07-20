#!/usr/bin/env python3

import argparse

translation_table = """
    AAs  = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
  Starts = ---M------**--*----M---------------M----------------------------
  Base1  = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
  Base2  = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
  Base3  = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
"""

line_hash = {}
for line in translation_table.split("\n"):
    line = line.strip()
    if not line:
        continue
    name, details = line.split("=")
    line_hash[name.strip()] = details.strip()

codon_to_amino_acid = {}

for amino_acid, base1, base2, base3 in zip(
        line_hash["AAs"],
        line_hash["Base1"],
        line_hash["Base2"],
        line_hash["Base3"]):
    codon_to_amino_acid[base1 + base2 + base3] = amino_acid

def rc(s):
  return "".join({'T':'A',
                  'G':'C',
                  'A':'T',
                  'C':'G',
                  'N':'N'}[x] for x in reversed(s))

parser = argparse.ArgumentParser()
parser.add_argument("genome")
args = parser.parse_args()

for genome in [
        args.genome[0:],
        args.genome[1:],
        args.genome[2:],
        rc(args.genome)[0:],
        rc(args.genome)[1:],
        rc(args.genome)[2:]
]:
    protein = []
    for codon_start in range(0, len(genome), 3):
        codon = genome[codon_start : codon_start + 3]
        if len(codon) == 3:
            protein.append(codon_to_amino_acid[codon])
    print("".join(protein))
            
