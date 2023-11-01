#!/bin/bash -l


ref_genome= DAOM_197198_chromosome_assembly.fasta

bwa index $ref_genome


# Filenames for PE raw reads
r1= ADNg-Ri-PE_ACTGAT_L001_R1.fastq.gz
r2= ADNg-Ri-PE_ACTGAT_L001_R2.fastq.gz

# align, convert SAM to BAM, and sort the BAM file
bwa mem $ref_genome $r1 $r2 | samtools view -bSu | samtools sort -o RiPE_bwa.bam

# Index the sorted BAM file
samtools index RiPE_bwa.bam

# Mark duplicates
java -jar picard.jar MarkDuplicates I=RiPE_bwa.bam O=RiPE_bwa_md.bam M=RiPE_md_metrics.txt

# Sort order coordinates
java -jar picard.jar SortSam I=RiPE_bwa_md.bam O=RiPE_bwa_mds.bam SORT_ORDER=coordinate

# Add read groups
java -jar picard.jar AddOrReplaceReadGroups I=RiPE_bwa_mds.bam O=RiPE_bwa_mdsRG.bam RGID=RiPE RGLB=RiPE RGPL=illumina RGPU=RiPE RGSM=RiPE

# Index the final BAM file
samtools index RiPE_bwa_mdsRG.bam
