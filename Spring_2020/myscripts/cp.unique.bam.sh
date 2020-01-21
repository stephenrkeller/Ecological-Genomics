#!/bin/bash

# Copy 1 bam file per individual over to reads2snps directory for mapping and variant calling


INDS=`cat 93.sam.bam.inds.txt`

cd /data/project_data/sam/

for IND in $INDS

do

cp $IND /data/project_data/snps/reads2snps/

done

 
