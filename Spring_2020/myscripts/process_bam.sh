#!/bin/bash

### Sorting SAM files and converting to BAM files
for f in ${output}/BWA/*.sam
do
	out=${f/.sam/}
	sambamba-0.7.1-linux-static view -S --format=bam ${f} -o ${out}.bam
	samtools sort ${out}.bam -o ${out}.sorted.bam
	rm ${out}.bam
done

### Removing PCR duplicates
for file in ${output}/BWA/${mypop}*.sorted.bam
do
	f=${file/.sorted.bam/}
	sambamba-0.7.1-linux-static markdup -r -t 10 ${file} ${f}.rmdup.bam
	samtools sort ${f}.rmdup.bam -o ${f}.sorted.rmdup.bam
done


### Stats on bwa alignments
for file in ${output}/BWA/${mypop}*.sorted.rmdup.bam
do
	f=${file/.sorted.rmdup.bam/}
	name=`basename ${f}`
	echo ${name} >> ${myrepo}/myresults/${mypop}.names.txt
	samtools flagstat ${file} | awk 'NR>=5&&NR<=13 {print $1}' | column -x
done >> ${myrepo}/myresults/${mypop}.flagstats.txt


### Nucleotide coverage
for file in ${output}/BWA/${mypop}*.sorted.rmdup.bam
do
	samtools depth ${file} | awk '{sum+=$3} END {print sum/NR}'
done >> ${myrepo}/myresults/${mypop}.coverage.txt



