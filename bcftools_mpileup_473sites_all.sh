#!/bin/bash -l


bcftools mpileup --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  -f DAOM_197198_chromosome_assembly.fasta -T  random_sites_from13_bams_473.txt1 -b bam_file_list_CDS.txt | bcftools query --format '%CHROM\t%POS\t%REF\t%DP\t%ALT[\t%AD]\n' > random_sites_from13_bams_473_bcftulsmpileup.txt

bcftools mpileup --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  -f DAOM_197198_chromosome_assembly.fasta -T 	random_sites_from13_bams_473.txt2 -b bam_file_list_CDS.txt | bcftools query --format '%CHROM\t%POS\t%REF\t%DP\t%ALT[\t%AD]\n' > random_sites_from13_bams_473_bcftulsmpileup.txt2

bcftools mpileup --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  -f DAOM_197198_chromosome_assembly.fasta -T 	random_sites_from13_bams_473.txt3 -b bam_file_list_CDS.txt | bcftools query --format '%CHROM\t%POS\t%REF\t%DP\t%ALT[\t%AD]\n' > random_sites_from13_bams_473_bcftulsmpileup.txt3

bcftools mpileup --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  -f DAOM_197198_chromosome_assembly.fasta -T 	random_sites_from13_bams_473.txt4 -b bam_file_list_CDS.txt | bcftools query --format '%CHROM\t%POS\t%REF\t%DP\t%ALT[\t%AD]\n' > random_sites_from13_bams_473_bcftulsmpileup.txt4

bcftools mpileup --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  -f DAOM_197198_chromosome_assembly.fasta -T 	random_sites_from13_bams_473.txt5 -b bam_file_list_CDS.txt | bcftools query --format '%CHROM\t%POS\t%REF\t%DP\t%ALT[\t%AD]\n' > random_sites_from13_bams_473_bcftulsmpileup.txt5

bcftools mpileup --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  -f DAOM_197198_chromosome_assembly.fasta -T 	random_sites_from13_bams_473.txt6 -b bam_file_list_CDS.txt | bcftools query --format '%CHROM\t%POS\t%REF\t%DP\t%ALT[\t%AD]\n' > random_sites_from13_bams_473_bcftulsmpileup.txt6

bcftools mpileup --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  -f DAOM_197198_chromosome_assembly.fasta -T 	random_sites_from13_bams_473.txt7 -b bam_file_list_CDS.txt | bcftools query --format '%CHROM\t%POS\t%REF\t%DP\t%ALT[\t%AD]\n' > random_sites_from13_bams_473_bcftulsmpileup.txt7

bcftools mpileup --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  -f DAOM_197198_chromosome_assembly.fasta -T 	random_sites_from13_bams_473.txt8 -b bam_file_list_CDS.txt | bcftools query --format '%CHROM\t%POS\t%REF\t%DP\t%ALT[\t%AD]\n' > random_sites_from13_bams_473_bcftulsmpileup.txt8

bcftools mpileup --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  -f DAOM_197198_chromosome_assembly.fasta -T 	random_sites_from13_bams_473.txt9 -b bam_file_list_CDS.txt | bcftools query --format '%CHROM\t%POS\t%REF\t%DP\t%ALT[\t%AD]\n' > random_sites_from13_bams_473_bcftulsmpileup.txt9

bcftools mpileup --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  -f DAOM_197198_chromosome_assembly.fasta -T 	random_sites_from13_bams_473.txt10 -b bam_file_list_CDS.txt | bcftools query --format '%CHROM\t%POS\t%REF\t%DP\t%ALT[\t%AD]\n' > random_sites_from13_bams_473_bcftulsmpileup.txt10
