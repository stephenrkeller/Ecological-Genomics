#!/bin/bash

# The command below uses blastn to search the contigs in the Centaurea Velvet de-novo reference GBS assembly fasta file
# (-query) to the already formatted Genbank nt database (-db) downloaded on Aug 11, 2017.
# You can enter 'blastn --help' for a list of the parameters.
# We choose the tab-delimited output format (6) and to only help the top hit (-max_target_seqs)
# and only if it has a minimum evalue of 0.00001.

# Set the path to the taxdb for taxonomy search
export BLASTDB=$HOME/Centaurea/

blastn -query ~/Centaurea/contigs.fa \
       -db /data/archive/databases/nt/nt \
       -out ~/Centaurea/BLAST_Velvet_nt.outfmt6 \
       -outfmt '6 qseqid sseqid evalue pident length mismatch gapopen qstart qend sstart send stitle sskingdoms staxids sscinames' \
       -num_threads 20 \
       -max_target_seqs 1
