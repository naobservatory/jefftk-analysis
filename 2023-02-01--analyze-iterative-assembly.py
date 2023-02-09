final_contig, intermediate_contigs = sys.argv[1:]

max_intermediate = 0
for contig in intermediate_contigs:
    n = int(contig.split("."))
    max_intermediate = max(n, max_intermediate)
