#!/bin/bash

# PBIO/BIO 381: Apring 2017
# Update, March 22, 2018

# This script runs the program 'Bayescan' designed to identify markers under diversifying or balancing selection by analyzing SNP divergence among populations 

for thin in 10 20 30
do
bayescan /data/project_data/beetles/snps/OTAU_2018_reads2snps_DP10GP95_biallelic_MAF01_Miss0.8.vcf.recode.vcf.bayescan \
 -threads 2 \
 -thin $thin \
 -pr_odds 10 \
 -od ~/myresults \
 -o bayescantest${thin}
done


