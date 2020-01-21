#!/bin/bash
      java -classpath /data/popgen/Trimmomatic-0.33/trimmomatic-0.33.jar org.usadellab.trimmomatic.TrimmomaticPE \
     		-threads 1 \
     		-phred33 \
    		 /data/project_data/beetles/rawdata/WA_PP1_YOURSAMPLE_R1.fastq.gz \
    		 /data/project_data/beetles/rawdata/WA_PP1_YOURSAMPLE_R2.fastq.gz \
    		 ~/cleanreads/"WA_PP1_YOURSAMPLE_R1_clean_paired.fa" \
    		 ~/cleanreads/"WA_PP1_YOURSAMPLE_R1_clean_unpaired.fa" \
    		 ~/cleanreads/"WA_PP1_YOURSAMPLE_R2_clean_paired.fa" \
    		 ~/cleanreads/"WA_PP1_YOURSAMPLE_R2_clean_unpaired.fa" \
    		 ILLUMINACLIP:/data/popgen/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 \
       		 LEADING:28 \
             TRAILING:28 \
             SLIDINGWINDOW:6:28 \
             HEADCROP:12 \
             MINLEN:35 \
