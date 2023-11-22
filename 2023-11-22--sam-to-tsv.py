#!/usr/bin/env python3

import sys
import json
import pysam

sam_in, tsv_out = sys.argv[1:]

with open("/home/ec2-user/mgs-pipeline/bowtie/genomeid-to-taxid.json") as inf:
    genomeid_to_taxid = json.load(inf)
    
with pysam.AlignmentFile(sam_in, "r")  as sam:
    with open(tsv_out, "w") as outf:
        for record in sam:
            genomeid = sam.get_reference_name(record.reference_id)
            taxid, genome_name = genomeid_to_taxid[genomeid]
            outf.write(
                "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (
                    record.query_name,
                    genomeid,
                    taxid,
                    record.cigarstring,
                    record.reference_start,
                    record.get_tag("AS"),
                    len(record.query_sequence)))
