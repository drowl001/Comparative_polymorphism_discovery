#!/bin/bash -l

freebayes -f DAOM_197198_chromosome_assembly.fasta  -p 10 -J -K -F 0.1 -0 -u  -b RiPE_bwa_mdsRG.bam > rirreg.chen.fb.RiPE.vcf

vcffilter -f "QUAL > 30 & TYPE = snp" rirreg.chen.fb.RiPE.vcf > rirreg.chen.fb.RiPE_filtered.vcf
