import sys
from Bio.SeqIO.FastaIO import SimpleFastaParser

f1, f2, K = sys.argv[1:]
K = int(K)

def to_kmers(fname):
    with open(fname) as inf:
        for title, sequence in SimpleFastaParser(inf): pass

    return set(sequence[i:i+K]
               for i in range(len(sequence) - K + 1))

f1_kmers = to_kmers(f1)
f2_kmers = to_kmers(f2)

for kmer in sorted(list(f1_kmers.intersection(f2_kmers))):
    print(kmer)
    
