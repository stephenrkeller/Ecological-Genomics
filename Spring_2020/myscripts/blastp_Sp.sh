#!/bin/bash

# This single line using the blastp command below will compare your transcript fasta file
# (-query) to the already formatted purple sea urchin, S. purpuratus protein database (-db).
# You can enter 'blastp --help' for a list of the parameters.
# We choose the tab-delimited output format (6) and to only help the top hit (-max_target_seqs)
# and only if it has a minimum evalue of 0.00001.

blastp -query /data/project_data/assembly/08-11-35-36_cl20_longest_orfs_gene.pep \
       -db /data/project_data/assembly/database/SPU_peptide.fasta \
       -out /data/project_data/assembly/blast/blastp_vs_Spurp.outfmt6 \
       -outfmt 6 \
       -num_threads 2 \
       -max_target_seqs 1
