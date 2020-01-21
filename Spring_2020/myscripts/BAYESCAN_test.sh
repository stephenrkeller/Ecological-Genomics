#!/bin/bash

# PBIO/BIO 381: Apring 2017

# This script runs the program 'Bayescan' designed to identify markers under diversifying or balancing selection by analyzing SNP divergence among populations 


bayescan SSW_all_biallelic.MAF0.02.Miss0.8.recode.bayescan \
 -threads 20 \
# The number of threads you want to use for processing
 -n 5000 \
# The number of samples of the mcmc chain to output for calculating the posterior probabilities
 -thin 10 \
# The thinning interval between recorded samples for the mcmc chain. The total number of mcmc iterations is the number of samples * the thinning interval
 -nbp 5 \
# The number of pilot runs to make for tuning the initial model parameters
 -pilot 5000 \
# The length of the pilot runs
 -burn 5000 \
# The length of the burn-in period to discard before recording samples from the posterior
 -pr_odds 10 \
# The prior odds for the probability that a locus is neutral and not under any selection
 -od ./
# Sets the directory for putting results
