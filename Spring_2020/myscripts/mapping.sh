#!/bin/bash

# Each student gets assigned a population to work with:
mypop="AB_05"

#Directory with demultiplexed fastq files
input="/data/project_data/RS_ExomeSeq/fastq/edge_fastq/pairedcleanreads/${mypop}"


# Output dir to store mapping files (bam)
output="/data/project_data/RS_ExomeSeq/mapping"


# Reference genome for aligning our reads
# Note -- this is a reduced version of the full Picea abies genome (>20 Gb!), containing just scaffolds with probes for our exome seqs
ref="/data/project_data/RS_ExomeSeq/ReferenceGenomes/Pabies1.0-genome_reduced.fa"


# Indexing the genome -- already done.  In the future, you'll need this step if working on a new project/genome
#bwa index ${ref}

# Aligning individual sequences to the reference

for forward in ${input}*_R1.cl.pd.fq
do
	reverse=${forward/_R1.cl.pd.fq/_R2.cl.pd.fq}
	f=${forward/_R1.cl.pd.fq/}
	name=`basename ${f}`
	echo "@ Aligning $name..."
	bwa mem -t 1 -M ${ref} ${forward} ${reverse} > ${output}/BWA/${name}.sam
done



