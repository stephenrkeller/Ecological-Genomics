#!/bin/bash 
cd /data/project_data/RS_ExomeSeq/fastq/edge_fastq/

ls -l AB*fastq.gz

for file in AB*fastq.gz

do

 fastqc "$file" -o ~/myresults/fastqc
 echo -e "\n\n"

done




