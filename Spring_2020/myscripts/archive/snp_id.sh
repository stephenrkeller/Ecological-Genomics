#!/bin/bash

cd /data/project_data/sam

#for file in *.sam

#do

#samtools view -@ 16 -bS "$file" >"$file.bam"

#done

#Merge all the ban files into one
#samtools merge merged.bam *.bam

#Make sure reads are in the correct order, probably unnecessary given the cleaning step
#samtools fixmate merged.bam merged.fixmate.bam

#Sort the reads with a different, faster program... 
#Also need to specify a different tmp directory or it will fill the rhel-root...
#sambamba_v0.6.0 sort -m 8G -t 8 -p --tmpdir=/data/project_data/sam/tmp/ merged.fixmate.bam merged.fixmate.sorted.sambada

# Remove PCR duplicates (same as Picard tools)
#sambamba_v0.6.0 markdup -t 8 --tmpdir=/data/project_data/sam/tmp/ merged.fixmate.sorted.bam merged.fixmate.sorted.mrkdup.bam

#samtools_019 rmdup merged.fixmate.sorted.bam merged.fixmate.sorted.rmdup.bam

#samtools index merged.fixmate.sorted.rmdup.bam

# Now for the SNP detection
samtools mpileup -uf /data/project_data/assembly/08-11-35-36_cl20_longest_orfs_gene.cds merged.fixmate.sorted.rmdup.bam | bcftools call -mv -> Poch_merged.raw.vcf
