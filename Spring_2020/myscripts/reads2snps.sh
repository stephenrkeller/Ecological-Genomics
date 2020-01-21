#!/bin/bash

########################################################################
# Calling variants using reads2snps
# Feb 2018
# S.R. Keller
########################################################################


./reads2snps/reads2snp_2.0.64.bin  \
  -min 5 \
  -th1 0.90 \
  -nbth 15\
  -bamlist ~/bam/bamlist.txt  \
  -bamref /data/otau/reference/OTAU.fna \
  -out OTAU_2018_reads2snps \
