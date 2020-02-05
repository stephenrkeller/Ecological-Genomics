#!/bin/bash

# Set your repo address here -- double check carefully!
myrepo="/users/s/r/srkeller/Ecological_Genomics/Spring_2020"

# Each student gets assigned a population to work with:
mypop="AB_05"


# Output dir to store mapping files (bam)
output="/data/project_data/RS_ExomeSeq/mapping"


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




