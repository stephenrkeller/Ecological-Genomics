#!/bin/bash

myrepo="/users/s/r/srkeller/Ecological_Genomics/Spring_2020"

OUT="/data/project_data/RS_ExomeSeq/ANGSD"

REF="/data/project_data/RS_ExomeSeq/ReferenceGenomes/Pabies1.0-genome_reduced.fa"


# List of bam files for this region
#ls /data/project_data/RS_ExomeSeq/mapping/BWA/*sorted.rm*.bam > ${OUT}/EDGE_bam.list


ANGSD -b ${myrepo}/myresults/EDGE_bam.list \
-ref ${REF} -anc ${REF} -out ${OUT}/PCA/EDGE_poly \
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
-doGeno 32 \
-doPost 1 \
-doCounts 1 \
-doMajorMinor 1 \
-doMaf 1 \
-doGlf 2 \
-SNP_pval 1e-6 \


##################################
#  PCA using ngsTools/ngsCovar
##################################


gunzip ${OUT}/PCA/EDGE_poly.geno.gz

NSITES=`cat ${OUT}/PCA/EDGE_poly.geno | tail -n+2 | wc -l`

/data/popgen/ngsTools/ngsPopGen/ngsCovar -probfile ${OUT}/PCA/EDGE_poly.geno -outfile ${OUT}/PCA/EDGE_poly.covar -nind 108 -nsites $NSITES -call 0 -norm 0


# PCA - Estimating the covariance matrix with pcangsd -- Not working!
#python /data/popgen/pcangsd/pcangsd.py -beagle ${myrepo}/myresults/EDGE_poly.beagle.gz -o ${myrepo}/myresults/_covmatrix -threads 15

#python /data/popgen/pcangsd/pcangsd.py -beagle ${myrepo}/myresults/EDGE_poly.beagle.gz -o EDGE_poly_covmatrix -threads 15



############################################

# Do association testing for plant height:

############################################

blups="/data/project_data/RS_ExomeSeq/ANGSD/blups"

OUT="/data/project_data/RS_ExomeSeq/ANGSD"

contigs="/data/project_data/RS_ExomeSeq/ANGSD/contig_splits"

mycontigs="xaa"


ANGSD -b ${OUT}/EDGE_bam.list \
-ref ${REF} \
-out ${OUT}/GWAS/Ht_EDGE_MD_PC12_${mycontigs} \
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
-doPost 1 \
-doCounts 1 \
-doMajorMinor 1 \
-doMaf 1 \
-SNP_pval 1e-6 \
-yQuant ${blups}/Ht_EDGE_MD.blups \
-doAsso 5 \
-rf ${contigs}/${mycontigs} \
-cov ${OUT}/PCA/EDGE_PC12.txt








