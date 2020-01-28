#!/bin/bash

# Pipeline to process and map RS exome sequences

# Set your repo address here -- double check carefully!
myrepo="/users/s/r/srkeller/Ecological_Genomics/Spring_2020/myresults"

# Each student gets assigned a population to work with:
mypop="AB"

#Directory with demultiplexed fastq files
input="/data/project_data/RS_ExomeSeq/fastq/edge_fastq/pairedcleanreads/${mypop}"


# Output dir to store mapping files (bam)
output="/data/project_data/RS_ExomeSeq/mapping"


#  Trim reads and count numbers of cleaned and paired reads

source ./trimmedReadCounts.sh


#  Map reads to ref genome using BWA

source /users/s/r/srkeller/Ecological_Genomics/Spring_2020/myscripts/mapping.sh


# Take sequence alignment  (sam) files and convert to bam>sort>remove PCR dups>sort again>index
# Calculate alignment stats for each individual and create table for my population

source /users/s/r/srkeller/Ecological_Genomics/Spring_2020/myscripts/process_bam.sh


