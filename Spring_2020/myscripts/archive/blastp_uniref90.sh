#!/bin/bash

# This single line using the blastp command below will compare your transcript fasta file
# (-query) to the already formatted uniref90 database (-db).
# You can enter 'blastp --help' for a list of the parameters.
# We choose the tab-delimited output format (6) and to only help the top hit (-max_target_seqs)
# and only if it has a minimum evalue of 0.00001.

blastp -query /data/project_data/assembly/08-11-35-36_cl20_longest_orfs_gene.pep \
       -db /data/archive/databases/uniref90/uniprot_uniref90.trinotate.pep \
       -out /data/project_data/assembly/blast/blastp_vs_uniref90.outfmt6 \
       -outfmt 6 \
       -num_threads 2 \
       -max_target_seqs 1
