#!/bin/bash

# Count reads per individual and read-pair
for forward in ${input}*_R1.cl.pd.fq
do
	reverse=${forward/_R1.cl.pd.fq/_R2.cl.pd.fq}
	f=${forward/_R1.cl.pd.fq/}
	name=`basename ${f}`
	echo ${name} >> ${myrepo}/myresults/TrimmedReadsInds.txt
	echo $(cat ${forward} | wc -l)/4 | bc >> ${myrepo}/myresults/TrimmedReadsCounts_R1.txt
	echo $(cat ${reverse} | wc -l)/4 | bc >> ${myrepo}/myresults/TrimmedReadsCounts_R2.txt
done


R --vanilla <<EOF

	reads_R1 = read.table("${myrepo}/myresults/TrimmedReadsCounts_R1.txt",header = FALSE)
	reads_R2 = read.table("${myrepo}/myresults/TrimmedReadsCounts_R2.txt",header = FALSE)
	ind_names = read.table("${myrepo}/myresults/TrimmedReadsInds.txt", header = FALSE)
	reads = cbind(ind_names, reads_R1, reads_R2)
	colnames(reads) = c("Individuals", "R1_counts", "R2_counts")
	write.table(reads, "${myrepo}/myresults/TrimmedReadCounts_table.txt", quote=F, row.names=F)
EOF
