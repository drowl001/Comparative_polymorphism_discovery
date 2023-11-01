#!/bin/bash -l


ref_genome= DAOM_197198_chromosome_assembly.fasta

bwa index $ref_genome

# Loop over each pair of PE reads
for i in {1..24}; do

  # Filenames for PE raw reads
  r1= rirreg${i}_R1.fastq.gz
  r2= rirreg${i}_R2.fastq.gz

  # align, convert SAM to BAM, and sort the BAM file
  bwa mem $ref_genome $r1 $r2 | samtools view -bSu | samtools sort -o rirreg${i}_bwa.bam

  # Index the sorted BAM file
  samtools index rirreg${i}_bwa.bam

  # Mark duplicates
  java -jar picard.jar MarkDuplicates I=rirreg${i}_bwa.bam O=rirreg${i}_bwa_md.bam M=rirreg${i}_md_metrics.txt

  # Sort order coordinates
  java -jar picard.jar SortSam I=rirreg${i}_bwa_md.bam O=rirreg${i}_bwa_mds.bam SORT_ORDER=coordinate

  # Add read groups
  java -jar picard.jar AddOrReplaceReadGroups I=rirreg${i}_bwa_mds.bam O=rirreg${i}_bwa_mdsRG.bam RGID=rirreg${i} RGLB=rirreg${i} RGPL=illumina RGPU=rirreg${i} RGSM=rirreg${i}

  # Index the final BAM file
  samtools index rirreg${i}_bwa_mdsRG.bam 

done
