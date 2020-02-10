#!/bin/bash

# Indexing the genome -- already done.  In the future, you'll need this step if working on a new project/genome
#bwa index ${ref}

ref="/data/project_data/RS_ExomeSeq/ReferenceGenomes/Pabies1.0-genome_reduced.fa"

# Aligning individual sequences to the reference

for forward in ${input}*_R1.cl.pd.fq
do
	reverse=${forward/_R1.cl.pd.fq/_R2.cl.pd.fq}
	f=${forward/_R1.cl.pd.fq/}
	name=`basename ${f}`
	echo "@ Aligning $name..."
	bwa mem -t 1 -M ${ref} ${forward} ${reverse} > ${output}/BWA/${name}.sam
done



