## Author: Stephen R. Keller  
## Ecological Genomics:   

### Overall Description of notebook      

My notes on coding and data analysis for the Spring 2020 Ecological Genomics course @ UVM


### Date started: 2020-Jan-22
### Date end:   (year-month-day)    

### Philosophy   
Science should be reproducible and one of the best ways to achieve this is by logging research activities in a notebook. Because science/biology has increasingly become computational, it is easier to document computational projects in an electronic form, which can be shared online through Github.    

### Helpful features of the notebook     

**It is absolutely critical for your future self and others to follow your work.**     

* The notebook is set up with a series of internal links from the table of contents.    
* All notebooks should have a table of contents which has the "Page", date, and title (information that allows the reader to understand your work).     
* Also, one of the perks of keeping all activities in a single document is that you can **search and find elements quickly**.     
* Lastly, you can share specific entries because of the three "#" automatically creates a link when the notebook renders on github.      


<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.  


### Table of contents for 60 entries (Format is *Page: Date(with year-month-day). Title*)        
* [Page 1: 2020-01-22](#id-section1). Intro to Github, Markdown, and UNIX command-line
* [Page 2: 2020-01-29](#id-section2). FastQC, read trimming, and mapping to reference genome
* [Page 3:](#id-section3).
* [Page 4:](#id-section4).
* [Page 5:](#id-section5).
* [Page 6:](#id-section6).
* [Page 7:](#id-section7).
* [Page 8:](#id-section8).
* [Page 9:](#id-section9).
* [Page 10:](#id-section10).
* [Page 11:](#id-section11).
* [Page 12:](#id-section12).
* [Page 13:](#id-section13).
* [Page 14:](#id-section14).
* [Page 15:](#id-section15).
* [Page 16:](#id-section16).
* [Page 17:](#id-section17).
* [Page 18:](#id-section18).
* [Page 19:](#id-section19).
* [Page 20:](#id-section20).
* [Page 21:](#id-section21).
* [Page 22:](#id-section22).
* [Page 23:](#id-section23).
* [Page 24:](#id-section24).
* [Page 25:](#id-section25).
* [Page 26:](#id-section26).
* [Page 27:](#id-section27).
* [Page 28:](#id-section28).
* [Page 29:](#id-section29).
* [Page 30:](#id-section30).
* [Page 31:](#id-section31).
* [Page 32:](#id-section32).
* [Page 33:](#id-section33).
* [Page 34:](#id-section34).
* [Page 35:](#id-section35).
* [Page 36:](#id-section36).
* [Page 37:](#id-section37).
* [Page 38:](#id-section38).
* [Page 39:](#id-section39).
* [Page 40:](#id-section40).
* [Page 41:](#id-section41).
* [Page 42:](#id-section42).
* [Page 43:](#id-section43).
* [Page 44:](#id-section44).
* [Page 45:](#id-section45).
* [Page 46:](#id-section46).
* [Page 47:](#id-section47).
* [Page 48:](#id-section48).
* [Page 49:](#id-section49).
* [Page 50:](#id-section50).
* [Page 51:](#id-section51).
* [Page 52:](#id-section52).
* [Page 53:](#id-section53).
* [Page 54:](#id-section54).
* [Page 55:](#id-section55).
* [Page 56:](#id-section56).
* [Page 57:](#id-section57).
* [Page 58:](#id-section58).
* [Page 59:](#id-section59).
* [Page 60:](#id-section60).

------
<div id='id-section1'/>

### Page 1: 2018-01-24. Notes on using Github, markdown, and the UNIX command-line

* Ask Lauren to give a brief intro on using git to create a repo and document your work in a lab notebook; push to origin on github

* I'll lead a tutorial on logging into the class unix server, doing some basic unix navigation and file manipulation, and writing simple for loops to run Fastqc.  Notes for all the 2020 tutorials are posted here: [2020 tutorial page](https://pespenilab.github.io/Ecological-Genomics/Tutorials.html)



* Went through  a bit of unix code showing how to log-in to the server

* Showed where the shared project data files live

* Went through unix commands:  `pwd, cd, ls, grep, cp, wc, pipe (|), redirect to file (>)`

* Showed how to use git to add, commit, and push to server

* Posted to Slack afterwards about using `git pull` prior to merge any changes from github first prior to pushing commits to the cloud

* On next Wednesday, pick up with using vim to edit files (using .bashrc as example) and then working with fastq files

* Plan should be to walk-through fastqs, assign population samples to students, have them run fastqc, have them design a trimomatic script to clean the reads up, then set up a bash script in `screen` to map reads with bwa and calculate genotype likelihoods 


------
<div id='id-section2'/>
### Page 2: 2018-01-29 

Here's the code we used to run FastQC:

```
#!/bin/bash 

# Pipeline to assess sequencing quality of RS exome sequences

# Set your repo address here -- double check carefully!
myrepo="/users/s/r/srkeller/Ecological_Genomics/Spring_2020"

# Make a new folder within 'myresults' to hold your fastqc outputs
mkdir ${myrepo}/myresults/fastqc

# Each student gets assigned a population to work with:
mypop="AB" 

#Directory with demultiplexed fastq files
input="/data/project_data/RS_ExomeSeq/fastq/edge_fastq/${mypop}"

#  Run fastqc program in a loop, processing all files that contain your population code
for file in ${input}*fastq.gz

do

 fastqc ${file} -o ${myrepo}/myresults/fastqc
 echo -e "\n\n   Results saved to ${myrepo}/myresults/fastqc...on to the next one!   \n\n"

done
```


Example output from FastQC -- copied to my 'docs' folder on GitHub for webpage display:

* [AB_05_R1_fastq](https://stephenrkeller.github.io/Ecological_Genomics/AB_05_R1_fastq_fastqc.html)
* [AB_05_R2_fastq](https://stephenrkeller.github.io/Ecological_Genomics/AB_05_R2_fastq_fastqc.html)

We then moved onto read trimming using Trimmomatic:

```
#!/bin/bash   
 
cd /data/project_data/RS_ExomeSeq/fastq/edge_fastq  
 
mkdir pairedcleanreads
mkdir unpairedcleanreads

for R1 in AB*R1_fastq.gz  

do 
 
 R2=${R1/_R1_fastq.gz/_R2_fastq.gz}
 short=`echo $R1 | cut -c1-5`
 echo $short 
java -classpath /data/popgen/Trimmomatic-0.33/trimmomatic-0.33.jar org.usadellab.trimmomatic.TrimmomaticPE \
        -threads 10 \
        -phred33 \
         "$R1" \
         "$R2" \
         /data/project_data/RS_ExomeSeq/fastq/edge_fastq/pairedcleanreads/"$short"_R1.cl.pd.fq \
         /data/project_data/RS_ExomeSeq/fastq/edge_fastq/unpairedcleanreads/"$short"_R1.cl.un.fq \
         /data/project_data/RS_ExomeSeq/fastq/edge_fastq/pairedcleanreads/"$short"_R2.cl.pd.fq \
         /data/project_data/RS_ExomeSeq/fastq/edge_fastq/unpairedcleanreads/"$short"_R2.cl.un.fq \
        ILLUMINACLIP:/data/popgen/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 \
        LEADING:20 \
        TRAILING:20 \
        SLIDINGWINDOW:6:20 \
        HEADCROP:12 \
        MINLEN:35 
 
done 

```

Lastly, we'll start mapping to the reference genome using BWA and use a bash script to process the resulting bam files to:
* sort
* mark and remove PCR duplicates
* sort again, and then index for quick processing later

Here's the mapping code:

```
#!/bin/bash

# Reference genome for aligning our reads

# Note -- this is a reduced version of the full Picea abies genome (>20 Gb!), containing just scaffolds with probes for our exome seqs
ref="/data/project_data/RS_ExomeSeq/ReferenceGenomes/Pabies1.0-genome_reduced.fa"

#number of CPU used -- set conservatively
t=1

# Indexing the genome -- already done.  In the future, you'll need this step if working on a new project/genome
#bwa index ${ref}

# Aligning individual sequences to the reference

for forward in ${input}*_R1.cl.pd.fq
do
	reverse=${forward/_R1.cl.pd.fq/_R2.cl.pd.fq}
	f=${forward/_R1.cl.pd.fq/}
	name=`basename ${f}`
	echo "@ Aligning $name..."
	bwa mem -t ${t} -M -a ${ref} ${forward} ${reverse} > ${output}/BWA/${name}.sam
done

### Sorting SAM files and converting to BAM files
###  Note, a similar program to samtools that is faster for 'view', 'flagstat', and 'markdup' is sambamba.  Use it when possible.

for f in ${output}/BWA/*.sam
do
	out=${f/.sam/}
	sambamba-0.7.1-linux-static view -S --format=bam ${f} -o ${out}.bam
	samtools sort ${out}.bam -o ${out}.sorted.bam
	rm ${out}.bam
done

```

And lastly the post-processing of the bam files:

```
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

R --vanilla <<EOF

	reads_Q = read.table("${myrepo}/myresults/${mypop}.Qscores.txt",header = FALSE)
	mean_cov = read.table("${myrepo}/myresults/${mypop}.coverage.txt",header = FALSE)
	ind_names = read.table("${myrepo}/myresults/${mypop}.names.txt", header = FALSE)
	qual_res = cbind (ind_names, reads_Q, mean_cov)
	colnames(qual_res) = c("Individuals", "Total_reads", "Total_mapped", "Total_mapq10", "Total_mapq20", "Total_mapq30", "Mean_cov")
	write.table(qual_res, "${myrepo}/myresults/${mypop}_MappingResults.txt")
EOF

for file in ${output}/BWA/${mypop}*.sorted.rmdup.bam
do
	samtools index ${file}
done

```

To make our scripts not too liong and cumbersome, and to enable more focused troubleshooting if bugs arise, we separated each of these major functions into the own scripts, and then ran them using a shell-script wrapper:

```
#!/bin/bash

# Pipeline to process and map RS exome sequences

# Set your repo address here -- double check carefully!
myrepo="/users/s/r/srkeller/Ecological_Genomics/Spring_2020"

# Each student gets assigned a population to work with:
mypop="AB"

#Directory with demultiplexed fastq files
input="/data/project_data/RS_ExomeSeq/fastq/edge_fastq/pairedcleanreads/${mypop}"


# Output dir to store mapping files (bam)
output="/data/project_data/RS_ExomeSeq/mapping"


#  Trim reads and count numbers of cleaned and paired reads

cd ${myrepo}/myscripts

source ./trimmedReadCounts.sh


#  Map reads to ref genome using BWA

source ./mapping.sh


# Take sequence alignment  (sam) files and convert to bam>sort>remove PCR dups>sort again>index
# Calculate alignment stats for each individual and create table for my population

source ./process_bam.sh

```


------
<div id='id-section3'/>
### Page 3:
------
<div id='id-section4'/>
### Page 4:
------
<div id='id-section5'/>
### Page 5:
------
<div id='id-section6'/>
### Page 6:
------
<div id='id-section7'/>
### Page 7:
------
<div id='id-section8'/>
### Page 8:
------
<div id='id-section9'/>
### Page 9:
------
<div id='id-section10'/>
### Page 10:
------
<div id='id-section11'/>
### Page 11:
------
<div id='id-section12'/>
### Page 12:
------
<div id='id-section13'/>
### Page 13:
------
<div id='id-section14'/>
### Page 14:
------
<div id='id-section15'/>
### Page 15:
------
<div id='id-section16'/>
### Page 16:
------
<div id='id-section17'/>
### Page 17:
------
<div id='id-section18'/>
### Page 18:
------
<div id='id-section19'/>
### Page 19:
------
<div id='id-section20'/>
### Page 20:
------
<div id='id-section21'/>
### Page 21:
------
<div id='id-section22'/>
### Page 22:
------
<div id='id-section23'/>
### Page 23:
------
<div id='id-section24'/>
### Page 24:
------
<div id='id-section25'/>
### Page 25:
------
<div id='id-section26'/>
### Page 26:
------
<div id='id-section27'/>
### Page 27:
------
<div id='id-section28'/>
### Page 28:
------
<div id='id-section29'/>
### Page 29:
------
<div id='id-section30'/>
### Page 30:
------
<div id='id-section31'/>
### Page 31:
------
<div id='id-section32'/>
### Page 32:
------
<div id='id-section33'/>
### Page 33:
------
<div id='id-section34'/>
### Page 34:
------
<div id='id-section35'/>
### Page 35:
------
<div id='id-section36'/>
### Page 36:
------
<div id='id-section37'/>
### Page 37:
------
<div id='id-section38'/>
### Page 38:
------
<div id='id-section39'/>
### Page 39:
------
<div id='id-section40'/>
### Page 40:
------
<div id='id-section41'/>
### Page 41:
------
<div id='id-section42'/>
### Page 42:
------
<div id='id-section43'/>
### Page 43:
------
<div id='id-section44'/>
### Page 44:
------
<div id='id-section45'/>
### Page 45:
------
<div id='id-section46'/>
### Page 46:
------
<div id='id-section47'/>
### Page 47:
------
<div id='id-section48'/>
### Page 48:
------
<div id='id-section49'/>
### Page 49:
------
<div id='id-section50'/>
### Page 50:
------
<div id='id-section51'/>
### Page 51:
------
<div id='id-section52'/>
### Page 52:
------
<div id='id-section53'/>
### Page 53:
------
<div id='id-section54'/>
### Page 54:
------
<div id='id-section55'/>
### Page 55:
------
<div id='id-section56'/>
### Page 56:
------
<div id='id-section57'/>
### Page 57:
------
<div id='id-section58'/>
### Page 58:
------
<div id='id-section59'/>
### Page 59:
------
<div id='id-section60'/>
### Page 60:

------
