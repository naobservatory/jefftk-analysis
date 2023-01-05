import sys
from Bio.SeqIO.FastaIO import SimpleFastaParser

kmer, = sys.argv[1:]

for title, sequence in SimpleFastaParser(sys.stdin):
    if kmer in sequence:
        print(title)
