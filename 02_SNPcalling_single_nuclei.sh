#!/bin/bash -l

freebayes -f DAOM_197198_chromosome_assembly.fasta  -p 10 -F 0.1 -0 -u  -b rirreg10_bwa_mdsRG.bam rirreg11_bwa_mdsRG.bam  rirreg13_bwa_mdsRG.bam rirreg14_bwa_mdsRG.bam rirreg15_bwa_mdsRG.bam  rirreg19_bwa_mdsRG.bam rirreg1_bwa_mdsRG.bam  rirreg3_bwa_mdsRG.bam rirreg4_bwa_mdsR$

vcffilter -f "QUAL > 30 & TYPE = snp" rirreg.w3.fb.p10.vcf > rirreg.w3.fb.p10_filtered.vcf
