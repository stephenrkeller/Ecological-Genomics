#!/bin/bash

myrepo="/users/s/r/srkeller/Ecological_Genomics/Spring_2020"

OUT="/data/project_data/RS_ExomeSeq/ANGSD/PCA"

REF="/data/project_data/RS_ExomeSeq/ReferenceGenomes/Pabies1.0-genome_reduced.fa"


# List of bam files for this region
ls /data/project_data/RS_ExomeSeq/mapping/BWA/*sorted.rm*.bam > ${OUT}/EDGE_bam.list


ANGSD -b ${myrepo}/myresults/EDGE_bam.list \
-ref ${REF} -anc ${REF} -out ${OUT}/EDGE_poly \
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


gunzip ${OUT}/EDGE_poly.geno.gz

NSITES=`cat ${OUT}/EDGE_poly.geno | tail -n+2 | wc -l`

/data/popgen/ngsTools/ngsPopGen/ngsCovar -probfile ${OUT}/EDGE_poly.geno -outfile ${OUT}/EDGE_poly.covar -nind 108 -nsites $NSITES -call 0 -norm 0


# PCA - Estimating the covariance matrix with pcangsd -- Not working!
#python /data/popgen/pcangsd/pcangsd.py -beagle ${myrepo}/myresults/EDGE_poly.beagle.gz -o ${myrepo}/myresults/_covmatrix -threads 15

#python /data/popgen/pcangsd/pcangsd.py -beagle ${myrepo}/myresults/EDGE_poly.beagle.gz -o EDGE_poly_covmatrix -threads 15



############################################

# Do association testing for plant height:

############################################

blups="/data/project_data/RS_ExomeSeq/ANGSD/blups"
contigs="/data/project_data/RS_ExomeSeq/ANGSD/contig_splits"
mycontigs="xaa"
PCA="/data/project_data/RS_ExomeSeq/ANGSD/PCA"


ANGSD -b ${myrepo}/myresults/EDGE_bam.list \
-ref ${REF} \
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
-yQuant ${blups}/Ht_EDGE_VTMDNC.blups \
-doAsso 5 \
-rf ${contigs}/${mycontigs} \
-cov ${PCA}/EDGE_PC1.txt \
-out ${myrepo}/Ht_EDGE_all_pc1_${mycontigs}







