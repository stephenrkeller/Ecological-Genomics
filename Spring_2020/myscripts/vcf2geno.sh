#!/bin/bash

java -Xmx2G -jar /data/popgen/PGDSpider_2.0.9.0/PGDSpider2-cli.jar -inputfile ~/myresults/out.recode.vcf -inputformat VCF -outputfile ~/myresults/out.recode.vcf.geno -outputformat EIGENSOFT -spid ~/myscripts/beetle.spid
