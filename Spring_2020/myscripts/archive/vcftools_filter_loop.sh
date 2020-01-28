#!/bin/bash

for DP in {5..6}
do
vcftools --vcf OTAU_2018_bcftools.vcf --min-alleles 2 --max-alleles 2 --mac 2 \
  --minDP $DP --site-pi --out $DP
done


