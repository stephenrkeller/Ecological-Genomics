#!/bin/bash

# Calling SNP variants using the reads2snp pipeline developed by N. Galtier and implemented in Gayral et al. 2014

# Create a list of input bam files per individual for merging with samtools

#cd /data/project_data/snps/reads2snps

#NAMES=`cat unique_inds.txt`

#for NAME in $NAMES
#do
#  cd /data/project_data/sam
#  ls $NAME*.sam.bam >/data/project_data/snps/reads2snps/"$NAME.bamlist.txt"
#done

###################################
##  Step 1: Sort bam files ##
###################################

#cd /data/project_data/snps/reads2snps

#for file in /data/project_data/sam/*bwaaln.sam.bam
#do
#  sambamba_v0.6.0 sort -m 3G -t 8 -p --tmpdir=/data/project_data/snps/reads2snps/tmp/ "$file" "$file.sorted.bam"
#done

##################################################
##  Step 2: Merge sorted bam files for each ind ##
##################################################

#cd /data/project_data/sam

#chmod 775 *.sorted.bam
#chmod 775 *.sorted.bam.bai

#for NAME in $NAMES
#do
#  ls $NAME*.sorted.bam >"$NAME.sorted.bamlist.txt"
#  INPUT=`cat $NAME*.sorted.bamlist.txt`
#  sambamba_v0.6.0 merge -t 16 -p "$NAME.merged.sorted.bam" $INPUT
#done

#############################
## Step 3: Remove PCR dups ##
#############################

cd /data/project_data/sam

#for file in *.merged.sorted.bam
#do 
#sambamba_v0.6.0 markdup -t 16 -r 34.merged.sorted.bam 34.merged.sorted.bam.merged.sorted.mrkdup.bam

sambamba_v0.6.0 markdup -t 20 -r -p --tmpdir=/data-tmp/ --hash-table-size 1000000 --overflow-list-size 1000000 10.merged.sorted.bam 10.merged.sorted.bam.merged.sorted.mrkdup.bam
#done

chmod 775 10.merged.sorted.bam.merged.sorted.mrkdup.bam

# This is now ready to pass to reads2snps for mapping

###################################################################
## Step 4: Create a list of the bam files for input to reads2snp ##
###################################################################

# The bam list file must be a two-colon file including one line per analyzed individual, such as:

# file1.bam   indiv1_name
# file2.bam   indiv2_name

# To make this list from the directory where the sorted bam files are located:

#cd /data/project_data/sam

#ls *merged.sorted.mrkdup.bam >mergebamlist.txt
#grep -o -P '\d\d\_.*bam' bamlist.txt >bamfiles.txt
#grep -o -P '\d\d' mergebamlist.txt >mergebamnames.txt
#paste mergebamlist.txt mergebamnames.txt >SSW_byind.txt

#############################################################
## Step 5: Call SNPs from bam files and reference assembly ##
#############################################################

#cd /data/project_data/snps/reads2snps/

/data/popgen/reads2snp_2.0/reads2snp_2.0.64.bin \
   -bamlist SSW_by24inds.txt  \
   -bamref /data/project_data/assembly/08-11-35-36_cl20_longest_orfs_gene.cds

