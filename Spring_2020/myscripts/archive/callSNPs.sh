#!/bin/bash

########################################################################
# Script for calling variants using bcftools mpileup and call commands
# Feb 2018
# S.R. Keller
########################################################################


#!/bin/bash

# convert to bam

samtools view -b -@ 4 WA_PP1_F1.sam -o WA_PP1_F1.bam

# check stats on bam alignment file

samtools flagstat WA_PP1_F1.bam

# fix mate pairs

samtools fixmate WA_PP1_F1.bam WA_PP1_F1_fixmate.bam

# sort

samtools  sort -@ 4 WA_PP1_F1_fixmate.bam -o WA_PP1_F1_fixmate.sorted.bam 

# mark and remove duplicates

# /data/popgen/sambamba_v0.6.0 markdup -t 4 -r -p WA_PP1_F1_fixmate.sorted.bam WA_PP1_F1_fixmate.sorted.rmdup.bam

samtools rmdup WA_PP1_F1_fixmate.sorted.bam WA_PP1_F1_fixmate.sorted.rmdup.bam

# index

samtools index WA_PP1_F1_fixmate.sorted.rmdup.bam

# call SNPs with bcftools

bcftools mpileup \
  -Ou -f /data/project_data/beetles/reference/OTAU.fna /data/project_data/beetles/sam/WA_PP1_F1_fixmate.sorted.rmdup.bam \
  -C 50 --min-MQ 20 --min-BQ 20 --threads 4 --skip-indels --annotate AD,DP | \
  bcftools call -Ov -mv --format-fields GQ >WA_PP1_F1.vcf
 
