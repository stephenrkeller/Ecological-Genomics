#!/bin/bash

myrepo="/users/s/r/srkeller/Ecological_Genomics/Spring_2020"

OUT="/data/project_data/RS_ExomeSeq/ANGSD/div"

REF="/data/project_data/RS_ExomeSeq/ReferenceGenomes/Pabies1.0-genome_reduced.fa"


# List of bam files for this region
ls /data/project_data/RS_ExomeSeq/mapping/BWA/*sorted.rm*.bam > ${OUT}/EDGE_bam.list


ANGSD -b ${OUT}/EDGE_bam.list \
-ref ${REF} -anc ${REF} \
-out ${OUT}/EDGE_all_Folded \
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
-doCounts 1 \
-doMajorMinor 1 \
-doMaf 1 \
-doSaf 1 \
-fold 1

#Estimation of the SFS for all sites using the FOLDED SFS
realSFS ${OUT}/EDGE_all_Folded.saf.idx -maxIter 1000 -tole 1e-6 -P 1 > ${OUT}/EDGE_all_Folded.sfs

# Estimate thetas and stats using the SFS as a prior

ANGSD -b ${OUT}/EDGE_bam.list \
-ref ${REF} -anc ${REF} \
-out ${OUT}/EDGE_all_Folded \
-nThreads 10 \
-remove_bads 1 \
-C 50 \
-baq 1 \
-minMapQ 20 \
-minQ 20 \
-minInd 2 \
-setMinDepthInd 1 \
-setMaxDepthInd 17 \
-setMinDepth 3 \
-skipTriallelic 1 \
-GL 1 \
-doCounts 1 \
-doMajorMinor 1 \
-pest ${OUT}/EDGE_all_Folded.sfs \
-doSaf 1 \
-doThetas 1 \
-fold 1

thetaStat do_stat ${OUT}/EDGE_all_Folded.thetas.idx
thetaStat print ${OUT}/EDGE_all_Folded.thetas.idx

