#!/bin/bash

# Set your repo address here -- double check carefully!
myrepo="/users/s/r/srkeller/Ecological_Genomics/Spring_2020"

# Each student gets assigned a population to work with:
mypop="AB"


# Output dir to store mapping files (bam)
output="/data/project_data/RS_ExomeSeq/mapping"


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


### Stats on bwa alignments

echo -e "Num.reads R1 R2 Paired MateMapped Singletons MateMappedDiffChr" \
  >${myrepo}/myresults/${mypop}.flagstats.txt

for file in ${output}/BWA/${mypop}*.sorted.rmdup.bam
do
	f=${file/.sorted.rmdup.bam/}
	name=`basename ${f}`
	echo ${name} >> ${myrepo}/myresults/${mypop}.names.txt
	samtools flagstat ${file} | awk 'NR>=6&&NR<=12 {print $1}' | column -x
done >> ${myrepo}/myresults/${mypop}.flagstats.txt


### Nucleotide coverage
for file in ${output}/BWA/${mypop}*.sorted.rmdup.bam
do
	samtools depth ${file} | awk '{sum+=$3} END {print sum/NR}'
done >> ${myrepo}/myresults/${mypop}.coverage.txt




