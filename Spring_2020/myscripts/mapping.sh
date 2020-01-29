#!/bin/bash

# Reference genome for aligning our reads

# Note -- this is a reduced version of the full Picea abies genome (>20 Gb!), containing just scaffolds with probes for our exome seqs
ref="/data/project_data/RS_ExomeSeq/ReferenceGenomes/Pabies1.0-genome_reduced.fa"

#number of CPU used -- set conservatively
t=1

# Indexing the genome -- already done.  In the future, you'll need this step if working on a new project/genome
#bwa index ${ref}

# Aligning individual sequences to the reference

for forward in ${input}*_R1.cl.pd.fq
do
	reverse=${forward/_R1.cl.pd.fq/_R2.cl.pd.fq}
	f=${forward/_R1.cl.pd.fq/}
	name=`basename ${f}`
	echo "@ Aligning $name..."
	bwa mem -t ${t} -M -a ${ref} ${forward} ${reverse} > ${output}/BWA/${name}.sam
done


