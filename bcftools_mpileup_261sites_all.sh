#!/bin/bash -l


bcftools mpileup --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  -f DAOM_197198_chromosome_assembly.fasta -T 261_random_RiPE.bam.txt1 coding_reads_frm_RiPE_bwa_mdsRG.bam | bcftools query --format '%CHROM\t%POS\t%REF\t%DP\t%ALT[\t%AD]\n' > 261_random_RiPE.bam_sites1_bcftulsmpileup.txt

bcftools mpileup --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  -f DAOM_197198_chromosome_assembly.fasta -T 261_random_RiPE.bam.txt2 coding_reads_frm_RiPE_bwa_mdsRG.bam | bcftools query --format '%CHROM\t%POS\t%REF\t%DP\t%ALT[\t%AD]\n' > 261_random_RiPE.bam_sites2_bcftulsmpileup.txt


bcftools mpileup --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  -f DAOM_197198_chromosome_assembly.fasta -T 261_random_RiPE.bam.txt3 coding_reads_frm_RiPE_bwa_mdsRG.bam | bcftools query --format '%CHROM\t%POS\t%REF\t%DP\t%ALT[\t%AD]\n' > 261_random_RiPE.bam_sites3_bcftulsmpileup.txt


bcftools mpileup --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  -f DAOM_197198_chromosome_assembly.fasta -T 261_random_RiPE.bam.txt4 coding_reads_frm_RiPE_bwa_mdsRG.bam | bcftools query --format '%CHROM\t%POS\t%REF\t%DP\t%ALT[\t%AD]\n' > 261_random_RiPE.bam_sites4_bcftulsmpileup.txt


bcftools mpileup --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  -f DAOM_197198_chromosome_assembly.fasta -T 261_random_RiPE.bam.txt5 coding_reads_frm_RiPE_bwa_mdsRG.bam | bcftools query --format '%CHROM\t%POS\t%REF\t%DP\t%ALT[\t%AD]\n' > 261_random_RiPE.bam_sites5_bcftulsmpileup.txt


bcftools mpileup --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  -f DAOM_197198_chromosome_assembly.fasta -T 261_random_RiPE.bam.txt6 coding_reads_frm_RiPE_bwa_mdsRG.bam | bcftools query --format '%CHROM\t%POS\t%REF\t%DP\t%ALT[\t%AD]\n' > 261_random_RiPE.bam_sites6_bcftulsmpileup.txt


bcftools mpileup --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  -f DAOM_197198_chromosome_assembly.fasta -T 261_random_RiPE.bam.txt7 coding_reads_frm_RiPE_bwa_mdsRG.bam | bcftools query --format '%CHROM\t%POS\t%REF\t%DP\t%ALT[\t%AD]\n' > 261_random_RiPE.bam_sites7_bcftulsmpileup.txt


bcftools mpileup --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  -f DAOM_197198_chromosome_assembly.fasta -T 261_random_RiPE.bam.txt8 coding_reads_frm_RiPE_bwa_mdsRG.bam | bcftools query --format '%CHROM\t%POS\t%REF\t%DP\t%ALT[\t%AD]\n' > 261_random_RiPE.bam_sites8_bcftulsmpileup.txt

bcftools mpileup --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  -f DAOM_197198_chromosome_assembly.fasta -T 261_random_RiPE.bam.txt9 coding_reads_frm_RiPE_bwa_mdsRG.bam | bcftools query --format '%CHROM\t%POS\t%REF\t%DP\t%ALT[\t%AD]\n' > 261_random_RiPE.bam_sites9_bcftulsmpileup.txt

bcftools mpileup --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR  -f DAOM_197198_chromosome_assembly.fasta -T 261_random_RiPE.bam.txt10 coding_reads_frm_RiPE_bwa_mdsRG.bam | bcftools query --format '%CHROM\t%POS\t%REF\t%DP\t%ALT[\t%AD]\n' > 261_random_RiPE.bam_sites10_bcftulsmpileup.txt
