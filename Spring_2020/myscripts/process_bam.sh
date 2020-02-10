#!/bin/bash

### Sorting SAM files and converting to BAM files
for f in ${output}/BWA/${mypop}*.sam
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
	sambamba-0.7.1-linux-static markdup -r -t 1 ${file} ${f}.sorted.rmdup.bam
done

### Indexing for fast lookup

for file in ${output}/BWA/${mypop}*.sorted.rmdup.bam
do
	samtools index ${file}
done




