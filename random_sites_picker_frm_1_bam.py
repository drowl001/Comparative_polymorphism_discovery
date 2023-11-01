#!/usr/bin/env python

import random
import sys
import pysam

def random_positions_from_bam(bam_file, output_file, num_positions):
    positions = []
    bam = pysam.AlignmentFile(bam_file, "rb")
    for read in bam:
        if not read.is_unmapped:  # Skip unmapped reads
            positions.append(f"{read.reference_name}\t{read.reference_start}")

    unique_positions = list(set(positions))  # Remove duplicates
    if len(unique_positions) < num_positions:
        print("Error: Insufficient unique positions in the input BAM file.")
        sys.exit(1)

    random.shuffle(unique_positions)
    selected_positions = unique_positions[:num_positions]

    with open(output_file, 'w') as outfile:
        outfile.write('\n'.join(selected_positions))

# Example usage
if len(sys.argv) != 3:
    print("Usage: python script.py input.bam output.txt")
    sys.exit(1)

input_bam = sys.argv[1]
output_file = sys.argv[2]
num_positions = 262

random_positions_from_bam(input_bam, output_file, num_positions)
print(f"Randomly selected positions written to {output_file}.")
