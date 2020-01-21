#!/bin/bash

for GQ in 10 15 20 25

do
vcftools --vcf OTAU_2018_samtools.vcf \
  --min-alleles 2 \
  --max-alleles 2 \
  --minGQ $GQ \
  --max-missing 0.5 \
  --mac 2 \
  --het \
  --out $GQ
done


