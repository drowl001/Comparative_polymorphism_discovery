#!/usr/bin/env python

#usage: ./vcf_compare_2files.py  <vcf file1> <vcf file2>


import argparse

def parse_vcf(file):
    """
    Parse a VCF file and return a set of (chromosome, position) tuples
    """
    sites = set()
    with open(file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom, pos = fields[0], int(fields[1])
            sites.add((chrom, pos))
    return sites

def write_vcf(sites, input_file, output_file):
    """
    Write a VCF file containing only the sites in the given set
    """
    with open(input_file) as f, open(output_file, 'w') as out:
        for line in f:
            if line.startswith('#'):
                out.write(line)
                continue
            fields = line.strip().split('\t')
            chrom, pos = fields[0], int(fields[1])
            if (chrom, pos) in sites:
                out.write(line)

def main(vcf1, vcf2):
    # Parse the VCF files
    vcf1_sites = parse_vcf(vcf1)
    vcf2_sites = parse_vcf(vcf2)

    # Find the common sites between the two VCF files
    common_sites = vcf1_sites.intersection(vcf2_sites)

    # Find the sites unique to each VCF file
    vcf1_unique = vcf1_sites - vcf2_sites
    vcf2_unique = vcf2_sites - vcf1_sites

    # Find the intersection sites between the two VCF files
    intersection_sites = vcf1_sites.intersection(vcf2_sites)

    # Write the results to output VCF files
    write_vcf(common_sites, vcf1, 'common_sites.vcf')
    write_vcf(vcf1_unique, vcf1, 'vcf1_unique.vcf')
    write_vcf(vcf2_unique, vcf2, 'vcf2_unique.vcf')
    write_vcf(intersection_sites, vcf1, 'intersection_sites.vcf')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compare two VCF files')
    parser.add_argument('vcf1', help='Path to first VCF file')
    parser.add_argument('vcf2', help='Path to second VCF file')
    args = parser.parse_args()
    main(args.vcf1, args.vcf2)
