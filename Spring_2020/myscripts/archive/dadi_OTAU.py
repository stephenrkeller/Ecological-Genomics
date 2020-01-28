#!/bin/python

import dadi
import numpy

# The following only need to be done once, after converting from vcf to SNP format.
# Once the frequency spectrum format conversion has been written to file and saved, you can just pick up the script by importing the fs within python

#dd = dadi.Misc.make_data_dict('../myresults/out.recode.vcf.data')

#fs = dadi.Spectrum.from_data_dict(dd, pop_ids=['IT','NC','WA'], projections=[35,35,35], polarized=False)

#fs.to_file('OTAU.fs')

# Import the fs file
fs = dadi.Spectrum.from_file('OTAU.fs')

print "Wright's FST = ", fs.Fst() # Get Wright's Fst, assuming random mating within pops
print
print "Num. segregating sites = ", fs.S() # Get estimate of number of segregating sites
print
#ns = fs.samples_sizes
#print "Sample sizes = ", ns

import demography_OTAU

pts_l = [40,50,60]

func = demography_OTAU.OutOfItaly

paramlist = ["Ital PopSize","Initial NC PopSize","NC PopSize","Initial WA PopSize","WA PopSize","Bottleneck PopSize","Bottleneck timing", "Time of NC/WA split"]

upper_bound = [100,100,100,100,100,100,1,1]
lower_bound = [1e-2,1e-2,1e-2,1e-2,1e-2,1e-2,0,0]

p0 = [2,0.1,2,0.1,2,0.1,0.2,0.1]

func_ex = dadi.Numerics.make_extrap_log_func(func)

print('Beginning optimization ************************************************')
popt = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l, 
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(p0), maxiter=5)
# The verbose argument controls how often progress of the optimizer should be
# printed. It's useful to keep track of optimization process.
print('Finshed optimization **************************************************')



