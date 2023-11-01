#!/usr/bin/env python

import random
import sys
import pysam

def random_positions_from_bam(bam_files, output_file):
    num_files = len(bam_files)
    num_positions = int(481 / num_files)  # Adjusted number of positions per file

    positions = []
    for bam_file in bam_files:
        bam = pysam.AlignmentFile(bam_file, "rb")
        for read in bam:
            if not read.is_unmapped:  # Skip unmapped reads
                positions.append(f"{read.reference_name}\t{read.reference_start}")
        bam.close()

    unique_positions = list(set(positions))  # Remove duplicates
    if len(unique_positions) < num_files * num_positions:
        print("Error: Insufficient unique positions in the input BAM files.")
        sys.exit(1)

    random.shuffle(unique_positions)
    selected_positions = unique_positions[:num_files * num_positions]

    with open(output_file, 'w') as outfile:
        outfile.write('\n'.join(selected_positions))

    print(f"Randomly selected positions written to {output_file}.")

# Example usage
if len(sys.argv) < 3:
    print("Usage: python script.py input1.bam input2.bam ... output.txt")
    sys.exit(1)

input_bam_files = sys.argv[1:-1]
output_file = sys.argv[-1]

random_positions_from_bam(input_bam_files, output_file)
