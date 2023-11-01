#!/usr/bin/env python

#usage: ./subset_vcf_file <positions.txt> <vcf file> <output vcf file>


import argparse

# Function to read positions from the text file into a list of tuples
def read_positions(file_path):
    positions = []
    with open(file_path, 'r') as f_pos:
        for line in f_pos:
            chrom, pos = line.strip().split('\t')
            pos = int(pos)
            positions.append((chrom, pos))
    return positions

# Function to process the input VCF file and filter lines based on positions
def filter_vcf(input_file_path, output_file_path, positions):
    with open(input_file_path, 'r') as f_in, open(output_file_path, 'w') as f_out:
        # Loop over input file
        for line in f_in:
            # Write header lines to output file
            if line.startswith('#'):
                f_out.write(line)
            else:
                # Parse chromosome name and position from input line
                chrom, pos, *_ = line.split('\t')
                pos = int(pos)
                # Check if position is in the sorted positions list
                if (chrom, pos) in positions:
                    # Write matching line to output file
                    f_out.write(line)

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Filter VCF file based on positions')
    parser.add_argument('positions_file', help='Path to the text file with positions')
    parser.add_argument('input_vcf', help='Path to the input VCF file')
    parser.add_argument('output_vcf', help='Path to the output VCF file')
    args = parser.parse_args()

    # Read positions from the text file into a list of tuples
    positions = read_positions(args.positions_file)

    # Sort positions by chromosome and position
    positions.sort()

    # Filter the input VCF file and write to the output VCF file
    filter_vcf(args.input_vcf, args.output_vcf, positions)

if __name__ == '__main__':
    main()
