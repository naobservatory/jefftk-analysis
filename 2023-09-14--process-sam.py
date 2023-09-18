#!/usr/bin/env python3

from Bio import SeqIO
from Bio import AlignIO
import pysam
import matplotlib.pyplot as plt

genomes = {}
with open("GCF_002374975.1_ASM237497v1_genomic.fna") as inf:
    for record in SeqIO.parse(inf, "fasta"):
        genomes[record.id] = record.seq
        
with pysam.AlignmentFile("SRR12204734.pair.sam", "r") as alignments:
    fig, axs = plt.subplots(nrows=20, figsize=(10, 20),
                            constrained_layout=True)
    n = 0

    for alignment in alignments:
        ref = alignments.get_reference_name(alignment.reference_id)
        if ref not in genomes:
            continue

        if n != 8:
            n += 1
            continue
        
        reference_sequence = genomes[ref]
        aligned_sequence = alignment.query_sequence
        
        alignment_start = alignment.reference_start
        alignment_end = alignment.reference_end

        # Plot the reference sequence
        axs[n].plot(range(alignment_start, alignment_end + 1),
                    reference_sequence[alignment_start:alignment_end + 1],
                    label="Reference", color="blue")

        print(reference_sequence[alignment_start:alignment_end + 1])
        print(aligned_sequence)
            
        
        # Plot the aligned sequence
        alignment_position = alignment.reference_start
        for cigar_operation in alignment.cigartuples:
            op_type, op_length = cigar_operation
            if (op_type == 0 or  # Match or Mismatch
                op_type == 7 or  # Match
                op_type == 8):   # Mismatch
                axs[n].plot(range(alignment_position,
                                  alignment_position + op_length),
                            list(aligned_sequence[:op_length]),
                            label="Aligned Sequence",
                            color="red")
                aligned_sequence = aligned_sequence[op_length:]
            elif (op_type == 1 or  # Insertion
                  op_type == 4):   # Soft clipping
                axs[n].plot([alignment_position - 0.5,
                             alignment_position + 0.5], ['*', '*'],
                            label="Insertion", color="green")
                aligned_sequence = aligned_sequence[op_length:]
            elif (op_type == 2 or  # Deletion
                  op_type == 3):   # Skipped region from reference
                axs[n].plot([alignment_position - 0.5,
                             alignment_position + 0.5], ['*', '*'],
                            label="Deletion", color="purple")
                alignment_position += op_length
            elif (op_type == 5 or  # Hard clipping
                  op_type == 6):   # Padding
                pass
            else:
                print(op_type)
                assert(False)

        n += 1

        if n >= 20:
            break
        
    plt.show()
