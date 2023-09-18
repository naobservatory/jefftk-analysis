#!/usr/bin/env python3

import sys
import glob
import pysam
import Bio.SeqIO
import matplotlib.pyplot as plt

sam_fname, *read_ids = sys.argv[1:]

genomes = {}
for genome in glob.glob("*.fna"):
    with open(genome) as inf:
        for record in Bio.SeqIO.parse(inf, "fasta"):
            genomes[record.id] = record.seq
        
with pysam.AlignmentFile(sam_fname, "r") as alignments:
    for alignment in alignments:
        if alignment.query_name not in read_ids:
            continue
        
        ref = alignments.get_reference_name(alignment.reference_id)
        if ref not in genomes:
            raise Exception("Unknown reference genome %s" % ref)

        fig, ax = plt.subplots(figsize=(10, 2), constrained_layout=True)

        
        reference_sequence = genomes[ref]
        aligned_sequence = alignment.query_sequence
        
        alignment_start = alignment.reference_start
        alignment_end = alignment.reference_end

        print(aligned_sequence)
        
        if alignment.reference_end is None:
            continue

        # Plot the reference sequence
        ax.plot(range(alignment_start, alignment_end + 1),
                reference_sequence[alignment_start:alignment_end + 1],
                label="Reference", color="blue")

        #print(reference_sequence[alignment_start:alignment_end + 1])
        #print(aligned_sequence)

        # print(alignment.query_name)
        # Plot the aligned sequence
        alignment_position = alignment.reference_start
        for cigar_operation in alignment.cigartuples:
            op_type, op_length = cigar_operation
            print(op_type, op_length)
            if (op_type == 0 or  # Match or Mismatch
                op_type == 7 or  # Match
                op_type == 8):   # Mismatch
                ax.plot(range(alignment_position,
                              alignment_position + op_length),
                        list(aligned_sequence[:op_length]),
                        label="Aligned Sequence",
                        color="red")
                aligned_sequence = aligned_sequence[op_length:]
                alignment_position += op_length
            elif (op_type == 1 or  # Insertion
                  op_type == 4):   # Soft clipping
                ax.plot([alignment_position - 0.5,
                         alignment_position + 0.5], ['*', '*'],
                        label="Insertion", color="green")
                aligned_sequence = aligned_sequence[op_length:]
            elif (op_type == 2 or  # Deletion
                  op_type == 3):   # Skipped region from reference
                ax.plot([alignment_position - 0.5,
                         alignment_position + 0.5], ['*', '*'],
                        label="Deletion", color="purple")
                alignment_position += op_length
            elif (op_type == 5 or  # Hard clipping
                  op_type == 6):   # Padding
                pass
            else:
                print(op_type)
                assert(False)

            print(alignment_position)
        
        plt.show()
        plt.close()
