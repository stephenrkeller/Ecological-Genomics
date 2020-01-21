#!/bin/bash

cd /data/project_data/snps/reads2snps

NAMES=`cat unique_inds.txt`

for NAME in $NAMES
do
  cd /data/project_data/sam/
  ls $NAME*.sam.bam >/data/project_data/snps/reads2snps/"$NAME.bamlist.txt"
done





#for NAME in $NAMES
#do
#  echo "$NAME.bam"
#done

