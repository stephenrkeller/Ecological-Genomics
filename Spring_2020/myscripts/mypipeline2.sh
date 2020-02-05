#!/bin/bash

# We'll use this as a wrapper to run our different mapping scripts

# Path to my repo:
myrepo="/users/s/r/srkeller/Ecological_Genomics/Spring_2020"

# My population:

mypop="AB"

# Directory to our cleaned and paired reads:

input="/data/project_data/RS_ExomeSeq/fastq/edge_fastq/pairedcleanreads/${mypop}"

# Directory to store the outputs of our mapping

output="/data/project_data/RS_ExomeSeq/mapping"

# Run mapping.sh

source ./mapping.sh

# Run the post-processing steps

source ./process_bam.sh


