#!/bin/bash

java -Xmx2G -jar /data/popgen/PGDSpider_2.0.9.0/PGDSpider2-cli.jar \
-inputfile ~/myresults/OTAU_2018_reads2snps_DP10GP95_biallelic_MAF01_Miss0.8.vcf.recode.vcf -inputformat VCF \
-outputfile ~/myresults/test.bayescan \
-outputformat GESTE_BAYE_SCAN -spid ~/myscripts/vcf2bayescan.spid
