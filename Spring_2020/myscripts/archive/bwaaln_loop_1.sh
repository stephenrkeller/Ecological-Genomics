#!/bin/bash

# To run it from the present directory: ./bwaaln.sh > output.bwaaln.txt
cd /data/project_data/fastq/cleanreads

MyLeftReads=`find . -maxdepth 1 -name "1*_R1*.fq.gz"`
 
for myLeft in $MyLeftReads
do
    #The mate file will have a _2 instead of a _1. This changes the first occurence of _1 with _2.
    myRight=${myLeft/_R1/_R2}
    #Is there a mate file?
    if [ -f $myRight ]
    then
        echo $myLeft $myRight
        myShort=`echo $myLeft | cut -c3-13`
        echo $myShort
        bwa aln -t 4 /data/project_data/assembly/08-11-35-36_cl20_longest_orfs.cds $myLeft > /data/project_data/sam/$myLeft".sai" 
        bwa aln -t 4 /data/project_data/assembly/08-11-35-36_cl20_longest_orfs.cds $myRight > /data/project_data/sam/$myRight".sai" 
        bwa sampe -r '@RG\tID:'"$myShort"'\tSM:'"$myShort"'\tPL:Illumina' \
                 -P /data/project_data/assembly/08-11-35-36_cl20_longest_orfs.cds /data/project_data/sam/$myLeft".sai" /data/project_data/sam/$myRight".sai" \
                 $myLeft \
                 $myRight > /data/project_data/sam/$myShort"_bwaaln.sam"
        rm /data/project_data/sam/$myLeft".sai"
        rm /data/project_data/sam/$myRight".sai"
    fi
done
