

############################################
#####     For the whole Edge region    #####
############################################

# Region name
REGION="EDGE"

# List of bam files for this region
ls /data/project_data/RS_ExomeSeq/mapping/BWA/*sorted.rmdup.bam > ${REGION}_bam.list

# Reference
REF="/data/project_data/RS_ExomeSeq/ReferenceGenomes/Pabies1.0-genome_reduced.fa"

# Estimating genotype likelihoods and other outputs with ANGSD
ANGSD -b ./${REGION}_bam.list -GL 1 -out ./${REGIONS}all -ref ${REF} -anc ${REF} -nThreads 15 -uniqueOnly 1 -remove_bads 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -skipTriallelic 1 -setMinDepthInd 1 -setMaxDepthInd 17 -setMinDepth 3 -minInd 2 -doHWE 1 -doGeno 32 -doPost 1 -doCounts 1 -doMaf 1 -doMajorMinor 1 -doGlf 1,2 -doSaf 1

# Estimating genotype likelihoods and other outputs with ANGSD
ANGSD -b ./${REGION}_bam.list -GL 1 -out ./${REGIONS}poly -ref ${REF} -anc ${REF} -nThreads 15 -uniqueOnly 1 -remove_bads 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -skipTriallelic 1 -setMinDepthInd 1 -setMaxDepthInd 17 -setMinDepth 3 -minInd 2 -doHWE 1 -doGeno 32 -doPost 1 -doCounts 1 -doMaf 1 -doMajorMinor 1 -doGlf 1,2 -doSaf 1 -SNP_pval 1e-6

#EM optimization of the sfs
realSFS ./${REGION}poly.saf.idx -maxIter 50000 -tole 1e-6 -P 4 > ./${REGIONS}poly.sfs

# Estimating thetas
ANGSD -b ./${REGIONS}_bam.list -ref ${REF} -anc ${REF} -out ./${REGIONS} -pest ./${REGIONS}.sfs -GL 1 -doSaf 1 -doThetas 1 -nThreads 1 -uniqueOnly 1 -remove_bads 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 2 -setMinDepthInd 1 -setMaxDepthInd 17 -setMinDepth 3 -skipTriallelic 1 
thetaStat do_stat ./${REGION}.thetas.idx

# PCA - Estimating the covariance matrix with pcangsd
python /data/popgen/pcangsd/pcangsd.py -beagle ${REGION}poly.beagle.gz -o ./${REGION}_covmatrix -threads 15




#######################################
#####     For each populations    #####
#######################################

# Population name
POP="AB"

# List of bam files for this population
ls /data/project_data/RS_ExomeSeq/mapping/BWA/${POP}_*sorted.rmdup.bam > ${POP}_bam.list

# Reference
REF="/data/project_data/RS_ExomeSeq/ReferenceGenomes/Pabies1.0-genome_reduced.fa"

# Estimating genotype likelihoods and other outputs with ANGSD
ANGSD -b ./${POP}_bam.list -ref ${REF} -anc ${REF} -out ./${POP} -nThreads 1 -uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 2 -setMinDepthInd 1 -setMaxDepthInd 17 -setMinDepth 3 -skipTriallelic 1 -GL 1 -doHWE 1 -doGeno 32 -doPost 1 -doCounts 1 -doMaf 1 -doMajorMinor 1 -doGlf 1 -dumpCounts 2 -doQsDist 1 -doDepth 1 -doSaf 1 #-SNP_pval 1e-6 
# if you want to restrict the estimation of the genotype likelihoods to some selected sites
# add the options "-sites ./selected_sites.txt" (tab delimited file with the position of the site in column 1 and the chromosome in column 2)  and "-rf ./selected_chromosome.chrs" (just the list of unique chromosome)

# Optimization and estimation of the SFS
realSFS ./${POP}.saf.idx -maxIter 50000 -tole 1e-6 -P 2 > ./${POP}.sfs

# Estimating thetas, need to use the sfs as prior to estimate the genotype likelihoods
ANGSD -b ./${POP}_bam.list -ref ${REF} -anc ${REF} -out ./${POP} -pest ./${POP}.sfs -GL 1 -doSaf 1 -doThetas 1 -nThreads 1 -uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 2 -setMinDepthInd 1 -setMaxDepthInd 17 -setMinDepth 3 -skipTriallelic 1 -doCounts 1 -doDepth 1 -dumpCounts 2 -doQsDist 1 -doMajorMinor 1
thetaStat do_stat ./${POP}.thetas.idx

### Estimating Fst estimates between one population and all the others ###
for POP2 in BFA # or from a file `cat ./list_pop`
do
	realSFS ./${POP}.saf.idx ./${POP2}.saf.idx -P 4 > ./${POP}_${POP2}.sfs
	realSFS fst index ./${POP}.saf.idx ./${POP2}.saf.idx -sfs ./${POP}_${POP2}.sfs -fstout ./${POP}_${POP2} -whichFst 1
	realSFS fst print ./${POP}_${POP2}.fst.idx > ./${POP}_${POP2}.fst
done
