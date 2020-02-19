#!/bin/bash

myrepo="/users/s/r/srkeller/Ecological_Genomics/Spring_2020"

#mkdir ${myrepo}/myresults/ANGSD

output="${myrepo}/myresults/ANGSD"

mypop="AB"


# List of bam files for this pop
ls /data/project_data/RS_ExomeSeq/mapping/BWA/${mypop}_*sorted.rm*.bam >${output}/${mypop}_bam.list

REF="/data/project_data/RS_ExomeSeq/ReferenceGenomes/Pabies1.0-genome_reduced.fa"

# Estimating GL's and allele frequencies for all sites with ANGSD
ANGSD -b ${output}/${mypop}_bam.list \
-ref ${REF} -anc ${REF} \
-out ${output}/${mypop}_allsites \
-nThreads 1 \
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
# -SNP_pval 1e-6

#Estimation of the SFS for all sites
realSFS ${output}/${mypop}_allsites.saf.idx -maxIter 1000 -tole 1e-6 -P 1 > ${output}/${mypop}_allsites.sfs

# Estimating thetas for all sites, using the SFS from above as a prior to estimate the GL's
ANGSD -b ${output}/${mypop}_bam.list \
-ref ${REF} -anc ${REF} \
-out ${output}/${mypop}_allsites \
-nThreads 1 \
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
-pest ${output}/${mypop}_allsites.sfs \
-doSaf 1 \
-doThetas 1

thetaStat do_stat ${output}/${mypop}_allsites.thetas.idx

#############################################################
#### Estimating thetas for all sites, using the FOLDED SFS
#### Not strictly necessary, but I was curious...
#############################################################

ANGSD -b ${output}/${mypop}_bam.list \
-ref ${REF} \
-out ${output}/${mypop}_outFold \
-nThreads 1 \
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
# -SNP_pval 1e-6

#Estimation of the SFS for all sites using the FOLDED SFS
realSFS ${output}/${mypop}_outFold.saf.idx -maxIter 1000 -tole 1e-6 -P 1 > ${output}/${mypop}_outFold.sfs

ANGSD -b ${output}/${mypop}_bam.list \
-ref ${REF} \
-out ${output}/${mypop}_outFold \
-nThreads 1 \
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
-pest ${output}/${mypop}_outFold.sfs \
-doSaf 1 \
-doThetas 1 \
-fold 1

thetaStat do_stat ${output}/${mypop}_outFold.thetas.idx
thetaStat print ${output}/${mypop}_outFold.thetas.idx

#################################################################################


### Estimating Fst estimates between your focal population and all the others ###
for POP2 in `cat ${myrepo}/mydata/edgepops_no${mypop}.txt`
do
	realSFS ${output}/${mypop}_allsites.saf.idx ${output}/${POP2}_allsites.saf.idx -P 1 > ${output}/${mypop}_${POP2}_allsites.sfs
	
	realSFS fst index ${output}/${mypop}_allsites.saf.idx ${output}/${POP2}_allsites.saf.idx -sfs ${output}/${mypop}_${POP2}_allsites.sfs -fstout ${output}/${mypop}_${POP2}_allsites -whichFst 1
	
	realSFS fst print ${output}/${mypop}_${POP2}_allsites.fst.idx > ${output}/${mypop}_${POP2}_allsites.fst
done







