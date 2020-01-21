#!/bin/bash

echo "Entering ~/bayenv2/pops"
cd ~/bayenv2/pops/

for f in *.pop
do
  echo "######################################################"
  echo "  "
  echo "Computing Allele Counts for Population $f"
  echo "  "
  echo "######################################################"
  vcftools --vcf /data/project_data/bayenv2/Intergenic/balsam_336inds_42pops_1353snp_Intergenic.vcf --keep $f --counts --out ../counts/$f
done

echo "Bulk Renaming Allele Count Files"
cd ~/bayenv2/counts/
rename 'pop.frq.' '' *.count
echo "Done!  Allele counts are in pop.count files!!"
