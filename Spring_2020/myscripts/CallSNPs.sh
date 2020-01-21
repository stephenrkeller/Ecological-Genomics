#!/bin/bash

#!/bin/bash

# convert to bam

samtools view -b -@ 4 WA_PP1_F1.sam -o WA_PP1_F1.bam

# check stats on bam alignment file

samtools flagstat WA_PP1_F1.bam >WA_PP1_F1.bamstats
samtools depth WA_PP1_F1.bam >WA_PP1_F1.depth

# fix mate pairs

samtools fixmate WA_PP1_F1.bam WA_PP1_F1_fixmate.bam

# sort

samtools  sort -@ 4 WA_PP1_F1_fixmate.bam -o WA_PP1_F1_fixmate.sorted.bam 

# mark and remove duplicates

samtools rmdup WA_PP1_F1_fixmate.sorted.bam WA_PP1_F1_fixmate.sorted.rmdup.bam

# index

samtools index WA_PP1_F1_fixmate.sorted.rmdup.bam

# call SNPs with samtools/bcftools


