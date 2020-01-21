#!/bin/bash

# Convert sam file to binary bam format

samtools view -b -@ 4 WA_PP1_F1_bwamem.sam -o WA_PP1_F1_bwamem.bam

# Check stats on our alignment

samtools flagstat WA_PP1_F1_bwamem.bam


# each position, how many reads
samtools depth WA_PP1_F1_bwamem.bam

# Fix mate pairs between paired end reads that have been orphaned
samtools fixmate WA_PP1_F1_bwamem.bam WA_PP1_F1_bwamem.fm.bam

# Sort alignment
samtools sort -@ 4 WA_PP1_F1_bwamem.fm.bam WA_PP1_F1_bwamem.fm.sort.bam

# Remove PCR duplicates
# samtools rmdup WA_PP1_F1_bwamem.fm.sort.bam WA_PP1_F1_bwamem.fm.sort.rmdup.bam

# Index final alignment for fast processing, makes has table
samtools index WA_PP1_F1_bwamem.fm.sort.bam

# Pause here; look at result before proceeding to SNP calling




