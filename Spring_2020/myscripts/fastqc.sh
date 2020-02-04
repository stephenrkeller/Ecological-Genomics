#!/bin/bash

cd ~/Ecological_Genomics/Spring_2020/myresults/

# I'm creating a new dir to store my results

mkdir fastqc

for file in /data/project_data/RS_ExomeSeq/fastq/edge_fastq/AB*fastq.gz

do

fastqc ${file} -o fastqc/

done


