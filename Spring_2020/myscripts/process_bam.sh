#!/bin/bash

### Removing PCR duplicates
for file in ${output}/BWA/${mypop}*.sorted.bam
do
	f=${file/.sorted.bam/}
	sambamba-0.7.1-linux-static markdup -r -t 10 ${file} ${f}.rmdup.bam
	samtools sort ${f}.rmdup.bam -o ${f}.sorted.rmdup.bam
done


# Stats on bwa alignments
for file in ${output}/BWA/${mypop}*.sorted.rmdup.bam
do
	f=${file/.sorted.rmdup.bam/}
	name=`basename ${f}`
	echo ${name} >> ${myrepo}/myresults/${mypop}.names.txt
	samtools flagstat ${file} | awk 'NR>=5&&NR<=13 {print $1}' | column -x
done >> ${myrepo}/myresults/${mypop}.flagstats.txt


# Reads mapping quality scores
for file in ${output}/BWA/${mypop}*.sorted.rmdup.bam
do
	samtools view ${file} | awk '$5>=0{c0++}; $5>0{c1++}; $5>9{c9++}; $5>19{c19++}; $5>29{c29++};  END {print c0 " " c1 " " c9 " " c19 " " c29}'
done >> ${myrepo}/myresults/${mypop}.Qscores.txt

# Nucleotide coverage
for file in ${output}/BWA/${mypop}*.sorted.rmdup.bam
do
	samtools depth ${file} | awk '{sum+=$3} END {print sum/NR}'
done >> ${myrepo}/myresults/${mypop}.coverage.txt

#R --vanilla <<EOF
#
#	reads_Q = read.table("${myrepo}/myresults/${mypop}.Qscores.txt",header = FALSE)
#	mean_cov = read.table("${myrepo}/myresults/${mypop}.coverage.txt",header = FALSE)
#	ind_names = read.table("${myrepo}/myresults/${mypop}.names.txt", header = FALSE)
#	qual_res = cbind (ind_names, reads_Q, mean_cov)
#	colnames(qual_res) = c("Individuals", "Total_reads", "Total_mapped", "Total_mapq10", "Total_mapq20", "Total_mapq30", "Mean_cov")
#	write.table(qual_res, "${myrepo}/myresults/${mypop}_MappingResults.txt")
#EOF

for file in ${output}/BWA/${mypop}*.sorted.rmdup.bam
do
	samtools index ${file}
done

# Deleting temporary files to save space
#rm ${output}/BWA/${mypop}*.sam
#rm ${output}/BWA/${mypop}*.sam
#rm ${output}/BWA/${mypop}*.rmdup.bam
#rm ${output}/BWA/${mypop}*.rmdup.bam.bai
