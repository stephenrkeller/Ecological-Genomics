#!/bin/bash

# All bash scripts start with the head line above. It tells the unix interpreter to read the following script as "bash" commands, and points it to the directory where the bash program lives

# Comments can be preceeded with '#'

# In this script, we want to automate running FASTQC to check for run quality over many different runs. This will save us time since we don't have to use a point and clock GUI interface, and we can automatically feed it as many files as we want within a single script.

# First step is to cd to where the data live. This wouldn't be necessary if the directory where you execute the script is the same as where the data are stored, but that's not always the case. To avoid confusion, it's a good habitat to always give your scripts explicit paths to data and programs.

cd /data/example_data/fastq

# Let's first see a list of what fastq files are here (recall the use of the '*' wildcard)

ls -l *fastq*

# We can print a few lines of data to screen to see what they look like. We'll first use the 'echo' tool, which simply prints to screen whatever your write.

echo -e "\n\n Here's what the fastq format looks like: \n\n"

# We can open up the compressed ata in tar.gz format without having to uncompress it. This is very handy, and saves lots of space by keeping the data compressed all the time. The tool 'zcat' is just like 'cat', but works on gz-ipped files

zcat srkeller_GBS_20150617_GBS9L1A_R1.fastq.gz | head -n 40

# Notice what we did in the above command? First, we used 'zcat' to open the file, then we sent the contents of the file on to the 'head' command using the vertical line chatacter '|', also known as the pipe. Using the pipe let's you directly send the output of one command as input to the next. You will use this all the time!

echo -e "\n\n Now, let's run 'FASTQC' to look at the quality of the run: \n\n"

# We can do an operation on many files at a time using 'for' loops. In this case, we're going to run FASTQC on each file in our directory that has a name starting with "poplar" and ending in .gz (note the use of the * wildcard). 

mkdir ~/fastqc_output

for file in ./*/poplar*.gz

do

  fastqc "$file"  -o ~/fastqc_output
  echo -e "\n\n"

done

# Notice the use of the '-o' flag, which directs the output to a use-defined directory. Here, we choose our home directory. A better choice might be to make a subdir within home specifically for FASTQC output, and direct all output there.

# Lastly, we can use 'scp' to copy the output files from our home directoty to our website so we can see them anywhere.

mkdir ~/public_html/fastqc

scp ~/fastqc_output/*.html srkeller@zoo.uvm.edu:~/public_html/fastqc/
