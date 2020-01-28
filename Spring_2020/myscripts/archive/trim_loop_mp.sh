#!/bin/bash   
 
cd /data/project_data/fastq/  
 
for f1 in *_R1.fq.gz  
 
do 
 
 f2=${f1%%_R1.fq.gz}"_R2.fq.gz"
 short=`echo $f1 | cut -c1-11`
 echo $short 
java -classpath /data/popgen/Trimmomatic-0.33/trimmomatic-0.33.jar org.usadellab.trimmomatic.TrimmomaticPE \
        -threads 10 \
        -phred33 \
         "$f1" \
         "$f2" \
         /data/project_data/fastq/cleanreads2/"$short"_R1.cl.pd.fq \
         /data/project_data/fastq/cleanreads2/"$short"_R1.cl.un.fq \
         /data/project_data/fastq/cleanreads2/"$short"_R2.cl.pd.fq \
         /data/project_data/fastq/cleanreads2/"$short"_R2.cl.un.fq \
        ILLUMINACLIP:/data/popgen/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 \
        LEADING:20 \
        TRAILING:20 \
        SLIDINGWINDOW:6:20 \
        HEADCROP:12 \
        MINLEN:35 
>> log.txt
 
done 
