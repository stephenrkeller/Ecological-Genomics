#!/bin/bash

cd ~/bayenv2/counts/

mkdir cut 
mkdir cut/ALT
mkdir cut/REF

sed -i -e "1d" *.count

for f in *.count
do
  cut -f5 $f > cut/REF/$f.REF.cut
  cut -f6 $f > cut/ALT/$f.ALT.cut
done


paste REF/*.cut > REF/REF_ALL.cut
paste ALT/*.cut > ALT/ALT_ALL.cut

cd ~/bayenv2/counts/
ls *.count > popnames
sed -i -- 's/.count//g' popnames
tr "\n" "\t" < popnames > popnames2
     
cat popnames2 REF/REF_ALL.cut > REF/REF_ALL2.cut
cat popnames2 ALT/ALT_ALL.cut > ALT/ALT_ALL2.cut


