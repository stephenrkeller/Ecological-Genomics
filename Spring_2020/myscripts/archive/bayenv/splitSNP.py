#!/usr/bin/env python

# Input file
input = raw_input("Name of the allele counts file:  ")
ac = open(input)

count = 1

while True:
	line1 = ac.readline()
	if line1 == '':
	  break
	line2 = ac.readline()
	acout = open("SNP_%d.txt" %count, 'w')
	acout.write(line1)
	acout.write(line2)
	acout.close()
	count = count+1
	
	
