###########################################################################
#####  Script used for the production of the results described in     #####
#####       Capblancq et al. 2020 Evolutionary Applications           #####
###########################################################################

###########################
#### MAPPING THE READS ####
###########################


#Directory with demultiplexed fastq.gz files
input="../../datashare/Spruce/exome_capture/fastq/cleaned_reads/"

#Directory for outputs
output="./Mapping_WES/ref_Pabies"


############ STEP 1 ##############
### Sorting and filtering data ###

# numer of reads in individuals
rm ${output}/NBreads_1P.txt
rm ${output}/NBreads_2P.txt
rm ${output}/NBreads_names.txt
for forward in ${input}*_fastq_1P.gz
do
	reverse=${forward/_fastq_1P.gz/_fastq_2P.gz}
	f=${forward/_fastq_1P.gz/}
	name=`basename ${f}`
	echo ${name} >> ${output}/NBreads_names.txt
	echo $(zcat ${forward} | wc -l)/4|bc >> ${output}/NBreads_1P.txt
	echo $(zcat ${reverse} | wc -l)/4|bc >> ${output}/NBreads_2P.txt
done

R --vanilla <<EOF

	reads_R1 = read.table("${output}/NBreads_1P.txt",header = FALSE)
	reads_R2 = read.table("${output}/NBreads_2P.txt",header = FALSE)
	ind_names = read.table("${output}/NBreads_names.txt", header = FALSE)
	bwa_res = cbind (ind_names, reads_R1, reads_R2)
	colnames(bwa_res) = c("Individuals", "Total_reads", "Total_mapped", "Total_mapq10", "Total_mapq30")
	write.table(bwa_res, "${output}/bwa_results.txt")
EOF


########### STEP 2 ###############
### Mapping the reads with bwa ###

#Reference for mapping
ref="./ReferenceGenomes/Pabies1.0-genome_reduced.fa"

#number of CPU used
t=10

# Indexing the genome
#bwa index ${ref}

# Aligning individual sequences to the reference
for forward in ${input}*_fastq_1P.gz
do
	f=${forward/_fastq_1P.gz/}
	name=`basename ${f}`
	echo "@ Aligning $name..."
	reverse=${forward/_fastq_1P.gz/_fastq_2P.gz}
    bwa mem -t ${t} -M -a ${ref} ${forward} ${reverse} > ${output}/BWA/${name}.sam
done

### Sorting SAM files and converting to BAM files
for f in ${output}/BWA/*.sam
do
	out=${f/.sam/}
	samtools view -bS ${f} -o ${out}.bam
	samtools sort ${out}.bam -o ${out}.sorted.bam
	rm ${out}.bam
done

# Stats on bwa alignments
rm ${output}/BWA/res.aln.reads.out
rm ${output}/BWA/names.txt
for INDIV in ${output}/BWA/*.sorted.bam
do
	f=${INDIV/.bam/}
	name=`basename ${f}`
	echo ${name} >> ${output}/BWA/names.txt
	samtools flagstat ${INDIV} | awk 'NR>=6&&NR<=13 {print $1}' | column -x
done >> ${output}/BWA/res.aln.reads.out

R --vanilla <<EOF

	align_res = read.table("${output}/BWA/res.aln.reads.out",header = FALSE)
	ind_names = read.table("${output}/BWA/names.txt", header = FALSE)
	bwa_res = cbind (ind_names, align_res)
	colnames(bwa_res) = c("Individuals", "Total_reads", "R1_reads", "R2_reads", "Properly_paired", "ReadMapped_MateUnmapped", "ReadUnmapped_MateMapped", "MateDifferentContig", "MateDifferentContig_Q5")
	Percent_properlypaired = as.numeric(bwa_res[,5]) / as.numeric(bwa_res[,2])
	bwa_res = data.frame(bwa_res, Percent_properlypaired=Percent_properlypaired)
	write.table(bwa_res, "${output}/BWA/bwa_results.txt")
EOF

# Reads mapping quality scores after bwa alignment
rm ${output}/BWA/reads_mapping_Qscores.txt
for file in ${output}/BWA/*.sorted.bam
do
	samtools view ${file} | awk '$5>0{c1++}; $5>29{c29++}; $5>=0{c0++}; $5>19{c19++}; $5>9{c9++} END {print c0 " " c1 " " c9 " " c19 " " c29}'
done >> ${output}/BWA/reads_mapping_Qscores.txt

# Nucleotide coverage on bwa .bam files
rm ${output}/BWA/mean_coverage.txt
for file in ${output}/BWA/*.sorted.bam
do
	samtools depth ${file} | awk '{sum+=$3} END {print sum/NR}'
done >> ${output}/BWA/mean_coverage.txt

R --vanilla <<EOF

	reads_Q = read.table("${output}/BWA/reads_mapping_Qscores.txt",header = FALSE)
	mean_cov = read.table("${output}/BWA/mean_coverage.txt",header = FALSE)
	ind_names = read.table("${output}/BWA/names.txt", header = FALSE)
	qual_res = cbind (ind_names, reads_Q, mean_cov)
	colnames(qual_res) = c("Individuals", "Total_reads", "Total_mapped", "Total_mapq10", "Total_mapq20", "Total_mapq30", "Mean_cov")
	write.table(qual_res, "${output}/BWA/bwa_quality_results.txt")
EOF


### Removing PCR duplicates
for file in ${output}/BWA/*.sorted.bam
do
	sambamba-0.6.8 markdup -r -t 2 ${file} ${file}.rmdup.bam
	samtools sort ${file}.rmdup.bam -o ${file}.rmdup.bam.sorted.bam
done

for file in ${output}/BWA/*.sorted.bam.rmdup.bam.sorted.bam
do
	f=${file/.sorted.bam.rmdup.bam.sorted.bam/}
	mv ${file} ${f}.final.bam
done

for file in ${output}/BWA/*.final.bam
do
	samtools index ${file}
done

# Deleting all the sam files and other temporary files to save space
rm ${output}/BWA/*.sam
rm ${output}/BWA/*.sorted.bam
rm ${output}/BWA/*.rmdup.bam
rm ${output}/BWA/*.rmdup.bam.bai

# Stats after PCR duplicates removal
rm ${output}/BWA/rmdup.res.aln.reads.out
rm ${output}/BWA/rmdup.names.txt
for INDIV in ${output}/BWA/*.final.bam
do
	f=${INDIV/final.bam/}
	name=`basename ${f}`
	echo ${name} >> ${output}/BWA/rmdup.names.txt
	samtools flagstat ${INDIV} | awk 'NR>=6&&NR<=13 {print $1}' | column -x
done >> ${output}/BWA/rmdup.res.aln.reads.out

R --vanilla <<EOF

	align_res = read.table("${output}/BWA/rmdup.res.aln.reads.out",header = FALSE)
	ind_names = read.table("${output}/BWA/rmdup.names.txt", header = FALSE)
	bwa_res = cbind (ind_names, align_res)
	colnames(bwa_res) = c("Individuals", "Total_reads", "R1_reads", "R2_reads", "Properly_paired", "ReadMapped_MateUnmapped", "ReadUnmapped_MateMapped", "MateDifferentContig", "MateDifferentContig_Q5")
	Percent_properlypaired = as.numeric(bwa_res[,5]) / as.numeric(bwa_res[,2])
	bwa_res = data.frame(bwa_res, Percent_properlypaired=Percent_properlypaired)
	write.table(bwa_res, "${output}/BWA/bwa_rmdup_results.txt")
EOF

# Reads mapping quality scores
rm ${output}/BWA/rmdup_reads_mapping_Qscores.txt
for file in ${output}/BWA/*.final.bam
do
	samtools view ${file} | awk '$5>0{c1++}; $5>29{c29++}; $5>=0{c0++}; $5>19{c19++}; $5>9{c9++} END {print c0 " " c1 " " c9 " " c19 " " c29}'
done >> ${output}/BWA/rmdup_reads_mapping_Qscores.txt

# Nucleotide coverage
rm ${output}/BWA/rmdup_mean_coverage.txt
for file in ${output}/BWA/*.final.bam
do
	samtools depth ${file} | awk '{sum+=$3} END {print sum/NR}'
done >> ${output}/BWA/rmdup_mean_coverage.txt

R --vanilla <<EOF

	reads_Q = read.table("${output}/BWA/rmdup_reads_mapping_Qscores.txt",header = FALSE)
	mean_cov = read.table("${output}/BWA/rmdup_mean_coverage.txt",header = FALSE)
	ind_names = read.table("${output}/BWA/rmdup.names.txt", header = FALSE)
	qual_res = cbind (ind_names, reads_Q, mean_cov)
	colnames(qual_res) = c("Individuals", "Total_reads", "Total_mapped", "Total_mapq10", "Total_mapq20", "Total_mapq30", "Mean_cov")
	write.table(qual_res, "${output}/BWA/bwa_rmdup_quality_results.txt")
EOF



##################################
####### DOWNSTREAM ANALYSES ######
##################################


##########################################
##### Genotype likelihood per REGION #####

# List of .bam files for each region
cat ./REGIONS/CORE_bam.list
cat ./REGIONS/MARGIN_bam.list
cat ./REGIONS/EDGE_bam.list

# Estimating the genotype likelihoods and SFS for each regional population separately with all the sites (both mono- and poly-morphic sites by removing the SNP_val option)
N=3
for POP in CORE MARGIN EDGE
do
    (
	count=`cat ./REGIONS/${POP}_bam.list | wc -l`
	MAX=$(($count+$count+$count+$count+$count))
	~/TOOLS/angsd/angsd -b ./REGIONS/${POP}_bam.list -ref ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/Pabies1.0-genome_reduced.fa -anc ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/Pabies1.0-genome_reduced.fa -out ./REGIONS/${POP}/${POP} -nThreads 4 -uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -minInd 2 -setMinDepthInd 2 -setMaxDepthInd 17 -setMinDepth 15 -setMaxDepth ${MAX} -skipTriallelic 1  -minMaf 0.003 -GL 1 -doHWE 1 -doGeno 32 -doPost 1 -doCounts 1 -doMaf 1 -doMajorMinor 1 -doGlf 1 -dumpCounts 2 -doQsDist 1 -doDepth 1 -doSaf 1 -fold 0 -hetbias_pval 0.001 -dosnpstat 1
    echo "starting task $POP.."
    sleep $(( (RANDOM % 3) + 1))
    ) &
    # allow only to execute $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
        # wait only for first job
        wait -n
    fi
done

# Finding loci intersection among regions
zcat ./REGIONS/CORE/CORE.mafs.gz | awk 'BEGIN { OFS = ":" }{ print $1,$2 }' | sed '1d' | sort > ./REGIONS/CORE/CORE_sites.txt
zcat ./REGIONS/EDGE/EDGE.mafs.gz | awk 'BEGIN { OFS = ":" }{ print $1,$2 }' | sed '1d' | sort > ./REGIONS/EDGE/EDGE_sites.txt
zcat ./REGIONS/MARGIN/MARGIN.mafs.gz | awk 'BEGIN { OFS = ":" }{ print $1,$2 }' | sed '1d' | sort > ./REGIONS/MARGIN/MARGIN_sites.txt

comm -12 ./REGIONS/CORE/CORE_sites.txt ./REGIONS/EDGE/EDGE_sites.txt > ./REGIONS/CORE_EDGE_sites.txt
comm -12 ./REGIONS/CORE_EDGE_sites.txt ./REGIONS/MARGIN/MARGIN_sites.txt > ./REGIONS/CORE_EDGE_MARGIN_sites.txt

sed 's/:/\t/' ./REGIONS/CORE_EDGE_MARGIN_sites.txt | sort -b -k1,1 > ./REGIONS/intersect.txt
cut -f1 ./REGIONS/intersect.txt | uniq | sort > ./REGIONS/intersect.chrs
~/TOOLS/angsd/angsd sites index ./REGIONS/intersect.txt

# Number of total sites with monomorphic and plolymorphic sites
#cat ./REGIONS/intersect.txt | wc -l # 30 823 071 sites

# Re-launching the genotype likelihoods estimation by keeping only intersecting sites (monomorphic and polymorphic) and making the optimization step on the SFS
N=3
for POP in CORE EDGE MARGIN 
do
    (
	~/TOOLS/angsd/angsd -b ./REGIONS/${POP}_bam.list -GL 1 -out ./REGIONS/${POP}/${POP}_intersect -ref ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/Pabies1.0-genome_reduced.fa -anc ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/Pabies1.0-genome_reduced.fa -sites ./REGIONS/intersect.txt -rf ./REGIONS/intersect.chrs -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -skipTriallelic 1 -doHWE 1 -doGeno 32 -doPost 1 -doCounts 1 -doMaf 1 -doMajorMinor 1 -doGlf 1 -doSaf 1 -fold 0 
    #EM optimization of the sfs
	~/TOOLS/angsd/misc/realSFS ./REGIONS/${POP}/${POP}_intersect.saf.idx -maxIter 50000 -tole 1e-6 -P 4 > ./REGIONS/${POP}/${POP}_intersect.sfs
    echo "starting task $POP.."
    sleep $(( (RANDOM % 3) + 1))
    ) &
    # allow only to execute $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
        # wait only for first job
        wait -n
    fi
done


##############################
##### Thetas, Tajima's D #####

# estimating thetas from bam files using the optimized SFS as a prior and keeping only the intersecting sites
N=3
for POP in EDGE MARGIN CORE 
do
    (
	~/TOOLS/angsd/angsd -bam ./REGIONS/${POP}_bam.list -out ./REGIONS/${POP}/${POP}_intersect -doThetas 1 -doSaf 1 -pest ./REGIONS/${POP}/${POP}_intersect.sfs -anc ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/Pabies1.0-genome_reduced.fa -ref ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/Pabies1.0-genome_reduced.fa -GL 1 -nThreads 4 -sites ./REGIONS/intersect.txt -rf ./REGIONS/intersect.chrs -uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20
	~/TOOLS/angsd/misc/thetaStat do_stat ./REGIONS/${POP}/${POP}_intersect.thetas.idx
    sleep $(( (RANDOM % 3) + 1))
    ) &
    # allow only to execute $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
        # wait only for first job
        wait -n
    fi
done

~/TOOLS/angsd/misc/thetaStat do_stat ./REGIONS/CORE/CORE_intersect.thetas.idx
##################################
##### Linkage disequilibrium #####

### LD per region

# Genotype likelihoods for intersected and polymorphic sites (beagle) 
N=3
for POP in MARGIN CORE EDGE 
do
    (
    #~/TOOLS/angsd/angsd -b ./REGIONS/${POP}_bam.list -GL 1 -out ./REGIONS/${POP}/${POP}_intersect_poly -ref ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/Pabies1.0-genome_reduced.fa -anc ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/Pabies1.0-genome_reduced.fa -sites ./REGIONS/intersect.txt -rf ./REGIONS/intersect.chrs -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -skipTriallelic 1 -SNP_pval 1e-6 -doMaf 1 -doMajorMinor 1 -doGlf 2
	# Estimating LD
	zcat REGIONS/${POP}/${POP}_intersect_poly.mafs.gz | awk '{print $1, $2}' > ./REGIONS/${POP}/pos_${POP}_intersect_poly.txt
	count=`cat ./REGIONS/${POP}_bam.list | wc -l`
    N_SITES=$((`zcat ./REGIONS/${POP}/${POP}_intersect_poly.mafs.gz | wc -l`-1))
	~/TOOLS/ngsLD/ngsLD --geno ./REGIONS/${POP}/${POP}_intersect_poly.beagle.gz --probs --n_ind ${count} --n_sites ${N_SITES} --pos ./REGIONS/${POP}/pos_${POP}_intersect_poly.txt --out ./REGIONS/${POP}/${POP}_intersect_poly_maf0.5.LD --min_maf 0.05
    sleep $(( (RANDOM % 3) + 1))
    ) &
    # allow only to execute $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
        # wait only for first job
        wait -n
    fi
done


R --vanilla <<EOF

	# Margin data
	data_margin <- read.table("./REGIONS/MARGIN/MARGIN_maf0.5.LD", nrows = 1000000)
	data_margin <- data_margin[-which(as.character(data_margin$V1)=="chr"),]

	# Edge data
	data_edge <- read.table("./REGIONS/EDGE/EDGE_maf0.5.LD", nrows = 1000000)
	data_edge <- data_edge[-which(as.character(data_edge$V1)=="chr"),]

	# Core data
	data_core <- read.table("./REGIONS/CORE/CORE_maf0.5.LD", nrows = 1000000)
	data_core <- data_core[-which(data_core$V1=="chr"),]

	####################################
	### Fitting an exponential decay ###
	
	# Function of LD decay (from https://www.pnas.org/content/98/20/11479)

	### CORE ###
	# Finding the initial parameters a and b
	Cstart <- c(C=0.1)
	n = 178
	modelC <- nls(LD ~ ((10+C*dist)/((2+C*dist)*(11+C*dist)))*(1+((3+C*dist)*(12+12*C*dist+(C*dist)^2))/(n*(2+C*dist)*(11+C*dist))), data = data.frame(dist=data_core[,5],LD=data_core[,6]), start = Cstart, control=nls.control(maxiter=100)) 
	# extract rho, the recombination parameter, 4Nr
	new.rho <- summary(modelC)$parameters[1]
	# feed in the new value of rho to obtain LD values adjusted for their distances along the chromosome/genome
	fpoints<-((10+new.rho*data_core[,5])/((2+new.rho*data_core[,5])*(11+new.rho*data_core[,5])))*(1+((3+new.rho*data_core[,5])*(12+12*new.rho*data_core[,5]+(new.rho*data_core[,5])^2))/(n*(2+new.rho*data_core[,5])*(11+new.rho*data_core[,5])))
	# final table
	ld.df <- data.frame(distance=data_core[,5],fpoints,rep("Core",length(data_core[,5])))
	colnames(ld.df) <- c("distance","fpoints","Region")
	ld.df_core<-ld.df[order(ld.df$distance),]

	### MARGIN ###
	Cstart <- c(C=0.1)
	n = 51
	modelC <- nls(LD ~ ((10+C*dist)/((2+C*dist)*(11+C*dist)))*(1+((3+C*dist)*(12+12*C*dist+(C*dist)^2))/(n*(2+C*dist)*(11+C*dist))), data = data.frame(dist=data_margin[,5],LD=data_margin[,6]), start = Cstart, control=nls.control(maxiter=100)) 
	new.rho <- summary(modelC)$parameters[1]
	fpoints<-((10+new.rho*data_margin[,5])/((2+new.rho*data_margin[,5])*(11+new.rho*data_margin[,5])))*(1+((3+new.rho*data_margin[,5])*(12+12*new.rho*data_margin[,5]+(new.rho*data_margin[,5])^2))/(n*(2+new.rho*data_margin[,5])*(11+new.rho*data_margin[,5])))
	ld.df <- data.frame(distance=data_margin[,5],fpoints,rep("Margin",length(data_margin[,5])))
	colnames(ld.df) <- c("distance","fpoints","Region")
	ld.df_margin<-ld.df[order(ld.df$distance),]

	### EDGE ###
	Cstart <- c(C=0.1)
	n = 110
	modelC <- nls(LD ~ ((10+C*dist)/((2+C*dist)*(11+C*dist)))*(1+((3+C*dist)*(12+12*C*dist+(C*dist)^2))/(n*(2+C*dist)*(11+C*dist))), data = data.frame(dist=data_edge[,5],LD=data_edge[,6]), start = Cstart, control=nls.control(maxiter=100)) 
	new.rho <- summary(modelC)$parameters[1]
	fpoints<-((10+new.rho*data_edge[,5])/((2+new.rho*data_edge[,5])*(11+new.rho*data_edge[,5])))*(1+((3+new.rho*data_edge[,5])*(12+12*new.rho*data_edge[,5]+(new.rho*data_edge[,5])^2))/(n*(2+new.rho*data_edge[,5])*(11+new.rho*data_edge[,5])))
	ld.df <- data.frame(distance=data_edge[,5],fpoints,rep("Edge",length(data_edge[,5])))
	colnames(ld.df) <- c("distance","fpoints","Region")
	ld.df_edge<-ld.df[order(ld.df$distance),]

	##############
	#### Plot ####

	cols <- c('Core'="#FEDF00", 'Margin'="#377EB8", 'Edge'="#4DAF4A")

	pdf("LD_decay_maf0.5.perRegion.pdf")
	plot(ld.df_core$distance,ld.df_core$fpoints, xlim=c(0,500), type = "n", xlab = "Distance (bp)", ylab = "LD (r2)")
	lines(ld.df_core$distance,ld.df_core$fpoints, lwd=2,col=cols[1], xlim=c(0,500))
	lines(ld.df_margin$distance,ld.df_margin$fpoints, lwd=2,col=cols[2], xlim=c(0,500))
	lines(ld.df_edge$distance,ld.df_edge$fpoints, lwd=2,col=cols[3], xlim=c(0,500))
	legend('topright','groups',c("Core","Margin","Edge"), lwd = 2, col=cols) #, ncol=5,nrow=2,bty ="n")
	dev.off()

EOF


####################################
##### PCA on the full sampling #####

### Pruning the full sampling dataset

# List of .bam files used
ls ~/mydata/Thibaut/WES_mapping/Mapping_WES/ref_Pabies/BWA/*.final.bam > ./all_bam.list

# Estimating genotype likelihood in output .beagle format for the polymorphic sites common in the three regional populations
count=`cat ./all_bam.list | wc -l`
~/TOOLS/angsd/angsd -b ./all_bam.list -ref ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/Pabies1.0-genome_reduced.fa -anc ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/Pabies1.0-genome_reduced.fa -fai ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/Pabies1.0-genome_reduced.fa.fai -out ./Full_Sampling_intersect_poly -sites ./REGIONS/intersect.txt -rf ./REGIONS/intersect.chrs -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -skipTriallelic 1 -SNP_pval 1e-6 -minMaf 0.003 -GL 1 -doHWE 1 -doGeno 32 -doPost 1 -doCounts 1 -doMaf 1 -doMajorMinor 1 -doGlf 2 -dovcf 1 -nInd ${count} 

# 926107 polymorphic sites

# LD estimation
zcat ./Full_Sampling_intersect.mafs.gz | awk '{print $1, $2}' > ./pos_Full_Sampling_intersect.txt
count=`cat ./all_bam.list | wc -l`
N_SITES=$((`zcat ./Full_Sampling_intersect.mafs.gz | wc -l`-1))
~/TOOLS/ngsLD/ngsLD --geno ./Full_Sampling_intersect.beagle.gz --probs --n_ind ${count} --n_sites ${N_SITES} --pos ./pos_Full_Sampling_intersect.txt --out ./Full_Sampling_intersect_maf0.5.LD --min_maf 0.05

# Selecting the sites to keep after pruning the dataset 
perl ~/TOOLS/ngsLD/scripts/prune_graph.pl --in_file ./Full_Sampling_intersect_maf0.5.LD --max_kb_dist 5 --min_weight 0.5 --out ./Full_SamplingLD_unlinked.id

##### ICI ne marche pas a cause d'un package manquant... #######


### PCA on the full sampling pruned dataset

# Estimating genotype likelihood in output .beagle format for the prunned dataset
count=`cat ./all_bam.list | wc -l`
~/TOOLS/angsd/angsd -b ./all_bam.list -ref ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/Pabies1.0-genome_reduced.fa -anc ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/Pabies1.0-genome_reduced.fa -fai ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/Pabies1.0-genome_reduced.fa.fai -out ./Full_Sampling_intersect_prunned -sites Full_SamplingLD_unlinked.id -rf Full_SamplingLD_unlinked.id -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -skipTriallelic 1 -SNP_pval 1e-6 -minMaf 0.003 -GL 1 -doHWE 1 -doGeno 32 -doPost 1 -doCounts 1 -doMaf 1 -doMajorMinor 1 -doGlf 2 -dovcf 1 -nInd ${count} 

# Estimating the covariance matrix with pcangsd ()
mkdir ./PCA
python ~/TOOLS/pcangsd/pcangsd.py -beagle ./Full_Sampling_intersect.beagle.gz -o ./PCA/Full_Sampling_intersect.covmatrix -threads 8 

# Anylysing the covariance matrix and plotting the results in R

R --vanilla <<EOF

	library(ade4)

	## PCA from covariance matrix
	TAB <- read.table("./PCA/Full_Sampling_intersect_prunned.covmatrix")
	pca <- eigen(TAB)

	## Info about individuals
	names <- read.table("names_bam.list")
	info_inds<-read.table("../../Info_samples_revised.txt", header=T)
	info_inds<-info_inds[match(as.vector(names[,1]), as.character(info_inds$Family)),]

	## Explain variance
	var <- pca$values/sum(pca$values)
	barplot(pca$values[1:100])
	barplot(var[1:100])
	
	## Figure
	pdf("./PCA_RedSpruce.pdf")
	plot(pca$vectors[,c(1,2)], type = "n", xlim = c(-0.1,0.1), xlab = "PC1 (15.8%)", ylab = "PC2 (0.5%)", main = "Complete sampling (3.5 millions SNP)")
	par(xpd=FALSE)
	abline(h = 0, lty=5, col = "black")
	abline(v = 0, lty=5, col = "black")
	s.class(pca$vectors[which(info_inds$Region=="M"),c(1,2)], fac=as.factor(info_inds$Site[which(info_inds$Region=="M")]), cellipse = 0, col = rep("#4DAF4A", length(levels(as.factor(info_inds$Site[which(info_inds$Region=="M")])))), add.plot = T, clabel = 0.8)
	s.class(pca$vectors[which(info_inds$Region=="C"),c(1,2)], fac=as.factor(info_inds$Site[which(info_inds$Region=="C")]), cellipse = 0, col = rep("#F9D017", length(levels(as.factor(info_inds$Site[which(info_inds$Region=="C")])))), add.plot = T, clabel = 0.8)
	s.class(pca$vectors[which(info_inds$Region=="E"),c(1,2)], fac=as.factor(info_inds$Site[which(info_inds$Region=="E")]), cellipse = 0, col = rep("#377EB8", length(levels(as.factor(info_inds$Site[which(info_inds$Region=="E")])))), add.plot = T, clabel = 0.8)
	legend(0.07,0.26, legend = c("Core", "Margin", "Edge"), border = "black", pch = 21, col = "black", pt.bg = c("#F9D017", "#4DAF4A", "#377EB8"), title = "Regions", text.font = 4, title.adj = 0.2)
	dev.off()

EOF


######################################
### Fst among regional populations ###

### Pruning the dataset

# Re-launching the genotype likelihoods estimation by keeping only intersecting and prunned sites
N=3
for POP in CORE EDGE MARGIN 
do
    (
	~/TOOLS/angsd/angsd -b ./REGIONS/${POP}_bam.list -GL 1 -out ./REGIONS/${POP}/${POP}_intersect_prunned -ref ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/Pabies1.0-genome_reduced.fa -anc ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/Pabies1.0-genome_reduced.fa -sites ./Full_SamplingLD_unlinked.id -rf ./Full_SamplingLD_unlinked.id -nThreads 10 -uniqueOnly 1 -remove_bads 1 -trim 1 -C 50 -baq 1 -minMapQ 20 -minQ 20 -skipTriallelic 1 -doHWE 1 -doGeno 32 -doPost 1 -doCounts 1 -doMaf 1 -doMajorMinor 1 -doGlf 1 -doSaf 1 -fold 0 
    #EM optimization of the sfs
	~/TOOLS/angsd/misc/realSFS ./REGIONS/${POP}/${POP}_intersect.saf.idx -maxIter 50000 -tole 1e-6 -P 4 > ./REGIONS/${POP}/${POP}_intersect.sfs
    echo "starting task $POP.."
    sleep $(( (RANDOM % 3) + 1))
    ) &
    # allow only to execute $N jobs in parallel
    if [[ $(jobs -r -p | wc -l) -gt $N ]]; then
        # wait only for first job
        wait -n
    fi
done


# Estimating the 2d SFS as a priori, using only 1000000 sites for memory limitations
~/TOOLS/angsd/misc/realSFS ./REGIONS/CORE/CORE_intersect_prunned.saf.idx ./REGIONS/EDGE/EDGE_intersect_prunned.saf.idx -P 4 -nSites 1000000 > ./REGIONS/CORE.EDGE.ml
~/TOOLS/angsd/misc/realSFS ./REGIONS/CORE/CORE_intersect_prunned.saf.idx ./REGIONS/MARGIN/MARGIN_intersect_prunned.saf.idx -P 4 -nSites 1000000 > ./REGIONS/CORE.MARGIN.ml
~/TOOLS/angsd/misc/realSFS ./REGIONS/MARGIN/MARGIN_intersect_prunned.saf.idx ./REGIONS/EDGE/EDGE_intersect_prunned.saf.idx -P 4 -nSites 1000000 > ./REGIONS/MARGIN.EDGE.ml

# Prepare the fst 
~/TOOLS/angsd/misc/realSFS fst index ./REGIONS/CORE/CORE_intersect_prunned.saf.idx ./REGIONS/EDGE/EDGE_intersect_prunned.saf.idx ./REGIONS/MARGIN/MARGIN_intersect_prunned.saf.idx -sfs ./REGIONS/CORE.EDGE.ml -sfs ./REGIONS/CORE.MARGIN.ml -sfs ./REGIONS/MARGIN.EDGE.ml -fstout Regions_Fst

# Get the global estimate
~/TOOLS/angsd/misc/realSFS fst stats Regions_Fst.fst.idx 



##################################
##### Stairway-plot analysis #####

# Using the SFS for each regional population obtained from the estimation and optimization of the SFS (see section Thetas, Tajima's D)
./REGIONS/CORE_intersect.sfs
./REGIONS/EDGE_intersect.sfs
./REGIONS/MARGIN_intersect.sfs

# Only not exomic sites ?????????

# Then see the following scripts: ?? / ?? / ?? 


#########################
##### Dadi analysis #####

## 3DSFS estimate for dadi analysis
~/TOOLS/angsd/misc/realSFS dadi -P 8 ./REGIONS/CORE/CORE_intersect.saf.idx ./REGIONS/MARGIN/MARGIN_intersect.saf.idx ./REGIONS/EDGE/EDGE_intersect.saf.idx -sfs ./REGIONS/CORE/CORE_intersect.sfs -sfs ./REGIONS/MARGIN/MARGIN_intersect.sfs -sfs ./REGIONS/EDGE/EDGE_intersect.sfs -ref ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/Pabies1.0-genome_reduced.fa -anc ~/mydata/Thibaut/WES_mapping/ReferenceGenomes/Pabies1.0-genome_reduced.fa -nSites 10000000 > ./REGIONS/CORE_MARGIN_EDGE_intersect.sfs

## To convert the sfs file into dadi like sfs 
~/mydata/Thibaut/RedSpruce_demography/DADI/ref_Pabies/realsfs2dadi.pl ./REGIONS/CORE_MARGIN_EDGE_intersect.sfs 178 52 110 > ~/mydata/Thibaut/RedSpruce_demography/DADI/ref_Pabies/RedSpruce.sfs

# Then see the following scripts: ?? / ?? / ?? 








count=`cat ./REGIONS/CORE_bam.list | wc -l`
N_SITES=$((`zcat ./REGIONS/CORE/CORE_intersect.mafs.gz | wc -l`-1))
~/TOOLS/ngsLD/ngsLD --geno ./REGIONS/CORE/CORE_intersect.beagle.gz --probs --n_ind ${count} --n_sites ${N_SITES} --pos ./REGIONS/CORE/pos_CORE_intersect.txt --out ./REGIONS/CORE/CORE_intersect_maf0.5.LD --min_maf 0.05
    
