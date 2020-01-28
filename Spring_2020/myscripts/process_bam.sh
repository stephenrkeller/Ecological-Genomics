#!/bin/bash

### Removing PCR duplicates
for file in ${output}/BWA/${mypop}*.sorted.bam
do
	f=${file/.sorted.bam/}
	sambamba-0.7.1-linux-static markdup -r -t 2 ${file} ${f}.rmdup.bam
	samtools sort ${f}.rmdup.bam -o ${f}.rmdup.sorted.bam
done


# Stats on bwa alignments
for INDIV in ${output}/BWA/${mypop}*.rmdup.sorted.bam
do
	f=${INDIV/.rmdup.sorted.bam/}
	name=`basename ${f}`
	echo ${name} >> ${myrepo}/${mypop}.names.txt
#	samtools flagstat ${INDIV} | awk 'NR>=5&&NR<=13 {print $1}' | column -x
done >> ${myrepo}/${mypop}.flagstats.txt


# Reads mapping quality scores
for file in ${output}/BWA/${mypop}*.rmdup.sorted.bam
do
	samtools view ${file} | awk '$5>=0{c0++}; $5>0{c1++}; $5>9{c9++}; $5>19{c19++}; $5>29{c29++};  END {print c0 " " c1 " " c9 " " c19 " " c29}'
done >> ${myrepo}/${mypop}.Qscores.txt

# Nucleotide coverage
for file in ${output}/BWA/${mypop}*.rmdup.sorted.bam
do
	samtools depth ${file} | awk '{sum+=$3} END {print sum/NR}'
done >> ${myrepo}/${mypop}.coverage.txt

R --vanilla <<EOF

	reads_Q = read.table("${myrepo}/${mypop}.Qscores.txt",header = FALSE)
	mean_cov = read.table("${myrepo}/${mypop}.coverage.txt",header = FALSE)
	ind_names = read.table("${myrepo}/${mypop}.names.txt", header = FALSE)
	qual_res = cbind (ind_names, reads_Q, mean_cov)
	colnames(qual_res) = c("Individuals", "Total_reads", "Total_mapped", "Total_mapq10", "Total_mapq20", "Total_mapq30", "Mean_cov")
	write.table(qual_res, "${myrepo}/${mypop}_MappingResults.txt")
EOF

for file in ${output}/BWA/${mypop}*.rmdup.sorted.bam
do
	samtools index ${file}
done

# Deleting temporary files to save space
#rm ${output}/BWA/${mypop}*.sam
#rm ${output}/BWA/${mypop}*.sam
#rm ${output}/BWA/${mypop}*.rmdup.bam
#rm ${output}/BWA/${mypop}*.rmdup.bam.bai
