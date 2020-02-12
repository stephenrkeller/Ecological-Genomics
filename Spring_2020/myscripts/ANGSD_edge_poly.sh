#!/bin/bash

myrepo="/users/s/r/srkeller/Ecological_Genomics/Spring_2020"

REF="/data/project_data/RS_ExomeSeq/ReferenceGenomes/Pabies1.0-genome_reduced.fa"

# List of bam files for this region
ls /data/project_data/RS_ExomeSeq/mapping/BWA/*sorted.rm*.bam > ${myrepo}/myresults/EDGE_bam.list


ANGSD -b ${myrepo}/myresults/EDGE_bam.list \
-ref ${REF} -anc ${REF} -out ${myrepo}/myresults/EDGE_poly \
-nThreads 10 \
-remove_bads 1 \
-C 50 \
-baq 1 \
-minMapQ 20 \
-minQ 20 \
-setMinDepth 3 \
-minInd 2 \
-setMinDepthInd 1 \
-setMaxDepthInd 17 \
-skipTriallelic 1 \
-GL 1 \
-doHWE 1 \
-doGeno 32 \
-doPost 1 \
-doCounts 1 \
-doMajorMinor 1 \
-doMaf 1 \
-doGlf 2 \
-dumpCounts 2 \
-doDepth 1 \
-SNP_pval 1e-6

# PCA - Estimating the covariance matrix with pcangsd
python /data/popgen/pcangsd/pcangsd.py -beagle ${myrepo}/myresults/EDGE_poly.beagle.gz -o ${myrepo}/myresults/_covmatrix -threads 15

python /data/popgen/pcangsd/pcangsd.py -beagle ${myrepo}/myresults/EDGE_poly.beagle.gz -o EDGE_poly_covmatrix -threads 15

