#!/bin/bash   
 
cd /data/project_data/RS_ExomeSeq/fastq/edge_fastq  
 
mkdir pairedcleanreads
mkdir unpairedcleanreads

for R1 in AB*R1_fastq.gz  

do 
 
 R2=${R1/_R1_fastq.gz/_R2_fastq.gz}
 short=`echo $R1 | cut -c1-5`
 echo $short 
java -classpath /data/popgen/Trimmomatic-0.33/trimmomatic-0.33.jar org.usadellab.trimmomatic.TrimmomaticPE \
        -threads 10 \
        -phred33 \
         "$R1" \
         "$R2" \
         /data/project_data/RS_ExomeSeq/fastq/edge_fastq/pairedcleanreads/"$short"_R1.cl.pd.fq \
         /data/project_data/RS_ExomeSeq/fastq/edge_fastq/unpairedcleanreads/"$short"_R1.cl.un.fq \
         /data/project_data/RS_ExomeSeq/fastq/edge_fastq/pairedcleanreads/"$short"_R2.cl.pd.fq \
         /data/project_data/RS_ExomeSeq/fastq/edge_fastq/unpairedcleanreads/"$short"_R2.cl.un.fq \
        ILLUMINACLIP:/data/popgen/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 \
        LEADING:20 \
        TRAILING:20 \
        SLIDINGWINDOW:6:20 \
        HEADCROP:12 \
        MINLEN:35 
 
done 
