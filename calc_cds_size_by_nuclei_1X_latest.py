#!/usr/bin/env python

import sys
import pysam

def calculate_cds_coverage(bam_filename, bed_filename, fasta_filename, output_filename):
    # Load the reference genome
    fasta_file = pysam.FastaFile(fasta_filename)

    # Initialize a dictionary to store coverage per coding region position
    coding_region_coverage = {}

    # Open the BAM file for reading
    bamfile = pysam.AlignmentFile(bam_filename, "rb")

    # Read the BED file and calculate coverage
    with open(bed_filename, 'r') as bed_file:
        for line in bed_file:
            chrom, start, end = line.strip().split('\t')
            start = int(start)
            end = int(end)

            # Initialize coverage counter
            coverage = 0

            # Calculate coverage within the region
            for pileupcolumn in bamfile.pileup(chrom, start, end+1):
                coverage = pileupcolumn.nsegments
                if coverage >= 1:
                    for pos in range(max(start, pileupcolumn.pos), min(end + 1, pileupcolumn.pos + 1)):
                        coding_region_coverage[(chrom, pos)] = 1

    # Calculate total bases covered at 1X coverage
    total_bases_1x_coverage = len(coding_region_coverage)

    # Write the result to the output file
    with open(output_filename, 'w') as output_file:
        output_file.write(f"Total bases in coding region covered at 1X: {total_bases_1x_coverage}\n")

    # Close the BAM file
    bamfile.close()

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python script.py BAM_FILE BED_FILE FASTA_FILE OUTPUT_FILE")
    else:
        bam_filename = sys.argv[1]
        bed_filename = sys.argv[2]
        fasta_filename = sys.argv[3]
        output_filename = sys.argv[4]
        calculate_cds_coverage(bam_filename, bed_filename, fasta_filename, output_filename)
