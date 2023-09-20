#!/usr/bin/env python3

import sys
import glob
import pysam
import Bio.SeqIO

sam_fname, *read_ids = sys.argv[1:]

genomes = {}
for genome in glob.glob("raw-genomes/*.fna"):
    with open(genome) as inf:
        for record in Bio.SeqIO.parse(inf, "fasta"):
            genomes[record.id] = record.seq

# copied from icdiff
def get_columns():
    def ioctl_GWINSZ(fd):
        try:
            import fcntl
            import termios
            import struct

            cr = struct.unpack(
                "hh", fcntl.ioctl(fd, termios.TIOCGWINSZ, "1234")
            )
        except Exception:
            return None
        return cr

    cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
    if cr and cr[1] > 0:
        return cr[1]
    return 80

COLUMNS = get_columns()

with pysam.AlignmentFile(sam_fname, "r") as alignments:
    for alignment in alignments:
        if alignment.query_name not in read_ids:
            continue

        ref = alignments.get_reference_name(alignment.reference_id)
        if ref not in genomes:
            raise Exception("Unknown reference genome %s" % ref)

        if alignment.reference_end is None:
            continue

        ref_seq = genomes[ref]
        qry_seq = alignment.query_sequence

        ref_pos = alignment.reference_start
        qry_pos = 0
        scr_pos = 0

        ref_row = []
        mid_row = []
        qry_row = []

        print(alignment.query_name)
        def maybe_print(force=False):
            global scr_pos
            if scr_pos >= COLUMNS or (force and mid_row):
                ref_line = "".join(ref_row)
                mid_line = "".join(mid_row)
                qry_line = "".join(qry_row)
                
                print(ref_pos - len(ref_line.replace("-", "")),
                      "...",
                      ref_pos)
                print(ref_line)
                print(mid_line)
                print(qry_line)
                print(qry_pos - len(qry_line.replace("-", "")),
                      "...",
                      qry_pos)
                print()

                ref_row.clear()
                mid_row.clear()
                qry_row.clear()
                scr_pos = 0

        for cigar_index, cigar_operation in enumerate(alignment.cigartuples):
            op_type, op_length = cigar_operation
            if (op_type == 0 or  # Match or Mismatch
                op_type == 7 or  # Match
                op_type == 8):   # Mismatch
                for _ in range(op_length):
                    base_ref = ref_seq[ref_pos]
                    base_qry = qry_seq[qry_pos]

                    ref_row.append(base_ref)
                    qry_row.append(base_qry)

                    if base_ref == base_qry:
                        mid_row.append("|")
                    else:
                        mid_row.append(" ")

                    ref_pos += 1
                    qry_pos += 1
                    scr_pos += 1

                    maybe_print()
            elif op_type == 4 and cigar_index == 0: # Soft clipping at beginning
                ref_pos -= op_length
                for _ in range(op_length):
                    base_ref = ref_seq[ref_pos]
                    base_qry = qry_seq[qry_pos]

                    ref_row.append(base_ref)
                    qry_row.append(base_qry)

                    if base_ref == base_qry:
                        mid_row.append(":")
                    else:
                        mid_row.append(" ")

                    ref_pos += 1
                    qry_pos += 1
                    scr_pos += 1

                    maybe_print()

            elif op_type == 4 and cigar_index == len(alignment.cigartuples) -1:
                # Soft clipping at end
                for _ in range(op_length):
                    base_ref = ref_seq[ref_pos]
                    base_qry = qry_seq[qry_pos]

                    ref_row.append(base_ref)
                    qry_row.append(base_qry)

                    if base_ref == base_qry:
                        mid_row.append(":")
                    else:
                        mid_row.append(" ")

                    ref_pos += 1
                    qry_pos += 1
                    scr_pos += 1

                    maybe_print()

            elif (op_type == 1 or  # Insertion
                  op_type == 4):   # Soft clipping
                for _ in range(op_length):
                    ref_row.append("-")
                    mid_row.append(" ")
                    base_qry = qry_seq[qry_pos]
                    qry_row.append(base_qry)
                    qry_pos += 1
                    scr_pos += 1
                    maybe_print()

            elif (op_type == 2 or  # Deletion
                  op_type == 3):   # Skipped region from reference
                for _ in range(op_length):
                    base_ref = ref_seq[ref_pos]
                    ref_row.append(base_ref)
                    qry_pos += 1
                    scr_pos += 1
                    mid_row.append(" ")
                    qry_row.append("-")
                    maybe_print()

            elif (op_type == 5 or  # Hard clipping
                  op_type == 6):   # Padding
                pass
            else:
                print(op_type)
                assert(False)

        maybe_print(force=True)
