#!/usr/bin/env python3

import random

GENOME_LEN=10000
N_CRITICAL_SECTIONS=10
CRITICAL_SECTION_LEN=100
READ_LEN=200
K=30

def simulate():
    critical_sections = [
        (start := random.randrange(GENOME_LEN - CRITICAL_SECTION_LEN),
         start + CRITICAL_SECTION_LEN)
        for _ in range(N_CRITICAL_SECTIONS)
    ]

    read_start = random.randrange(GENOME_LEN - READ_LEN)
    read_end = read_start + READ_LEN
    for critical_start, critical_end in critical_starts:
        if read_start < 

        

simulate()

        
