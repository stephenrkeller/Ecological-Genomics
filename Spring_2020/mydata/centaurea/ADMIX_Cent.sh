#!/bin/bash

# Run ADMIXTURE to determine the number of genetic clusters in the SNP data, 
# and the ancestry proportions of each individual

# Remember the utility of 'for loops'?

for K in 1 2 3 4 5 6 7 8 9 10

do 

admixture -C 0.000001 --cv -j20 ~/Centaurea_C1_SNP_genotyping_matrix_50pct.vcf.vcf.geno $K \
| tee log${K}.out

done

grep -h CV log*.out >chooseK_Cent.txt



