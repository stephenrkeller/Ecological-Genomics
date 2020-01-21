#!/bin/bash

java -Xmx512M -jar /data/popgen/PGDSpider_2.0.9.0/PGDSpider2-cli.jar -inputfile ./out.recode.vcf -inputformat VCF -outputfile ./Centaurea_C1_SNP_genotyping_matrix_50pct.vcf.vcf.geno -outputformat EIGENSOFT -spid ./vcf2admixture_Cent.spid
