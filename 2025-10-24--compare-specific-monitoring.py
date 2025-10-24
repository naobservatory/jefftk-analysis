#!/usr/bin/env python3

import random

READ_LENGTH = 170
TARGET_COVERAGE = 2
TARGET_OBSERVATIONS = 2
SIMULATIONS = 10_000

def simulate_single_position(genome_length):
    possible_read_starts = genome_length - READ_LENGTH
    target_position = random.randrange(genome_length)
    target_observations = 0
    n = 0
    while True:
        n += 1
        read_start = random.randrange(possible_read_starts)
        if read_start <= possible_read_starts < read_start + READ_LENGTH:
            target_observations += 1

            if target_observations == TARGET_OBSERVATIONS:
                return n

def calculate_coverage(genome_length):
    # Each read gives us READ_LENGTH bp, and to reach TARGET_COVERAGE we need
    # to see a number of bp from the genome that adds up to TARGET_COVERAGE.
    return TARGET_COVERAGE * genome_length / READ_LENGTH

def start():
    print("genome_length",
          "single_position_2x",
          "2x_coverage",
          sep="\t")
    for genome_length in range(5000, 100_000, 100):
        if genome_length > 10_000 and genome_length % 1000 != 0:
            continue

        single_position_results = [
            simulate_single_position(genome_length)
            for _ in range(SIMULATIONS)]
        single_position_results.sort()
        single_position_result = single_position_results[
            len(single_position_results)//2]

        print(genome_length,
              single_position_result,
              round(calculate_coverage(genome_length)),
              sep="\t", flush=True)

if __name__ == "__main__":
    start()
