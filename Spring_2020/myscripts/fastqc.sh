#!/bin/bash 

# Pipeline to assess sequencing quality of RS exome sequences

# Set your repo address here -- double check carefully!
myrepo="/users/s/r/srkeller/Ecological_Genomics/Spring_2020/myresults"

# Make a new folder within 'myresults' to hold your fastqc outputs
mkdir fastqc

# Each student gets assigned a population to work with:
mypop="AB" 

#Directory with demultiplexed fastq files
input="/data/project_data/RS_ExomeSeq/fastq/edge_fastq/${mypop}"

#  Run fastqc program in a loop, processing all files that contain your population code
for file in ${input}*fastq.gz

do

 fastqc ${file} -o ${myrepo}/fastqc
 echo -e "\n\n   Results saved to ${myrepo}/fastqc...on to the next one!   \n\n"

done




