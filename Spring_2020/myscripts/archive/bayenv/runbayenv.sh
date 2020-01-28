#!/bin/bash

for s in SNP*.txt
do
  bayenv2 -i $s -m MYNAME.covmat -p 42 -r RANDNUM -t -e ENVAR_PC3 -k 100000 -n 3 -f -c -X -o $s
done
