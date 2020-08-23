# Supplementary R code for "Bayesian adaptive decision-theoretic designs
# for multi-arm multi stage clinical trials" by
# Bassi, Berkhof, de Jong, van de Ven
#
# correspondence to: p.vandeven@amsterdamumc.nl

# This script requires that the R script 'Functions for selecting the best from
# three, four or five experimental arms.R' was previously run and packages snow, 
# doSNOW and doParallel for parallel computing were installed.
# The required R script can be downloaded from 
# https://github.com/PeterMVanDeVen/Decision-Theoretic-MAMS-trials

# This example illustrates the use of the following functions for the 
# setting with four experimental arms: 
#
# create.d  : Simulates trial data;
# prdrop    : Evaluates characteristics for trials when arms can be 
#             dropped;
# prnodrop  : Evaluates characteristics for trials without dropping of 
#             arms;
# pcalc     : Summarizes trial characteristics
# reevaluate: Reevaluates characteristics of trials obtained by prdrop
#             and prnodrop for higher thresholds C/Q  
#
# Required functions provided in separate R script

# Packages snow, doSNOW an doParallel are required for parallel 
# computing 

library(snow)
library(doSNOW)
library(doParallel)

# Specify the number of trials to be simulated (simulation size)

Ntrials <- 1000

# Specify the response rate vector under which the frequentist 
# properties of the trial are to be evaluated

resp <- c(0.5,0.8,0.9,0.5)

# Specify the cap (= maximum number of subjects included in the trial)
# In this example the cap is set at a high number that is not reached 
# for the response rate vector and threshold (C/Q) used for stopping

maxpt <- 1000

# Set margin delta at 0 for symmetric setting with only experimental
# arms

delta <- 0

# Set dropctrl at TRUE to allow dropping of all experimental arms 

dropctrl <- TRUE

# Specify batch size for stage 1 
# 'burn' should be specified as the total batch size divided by
# the number of arms at start.
# 'burn <- 4' here corresponds to 16 subjects in the first stage

burn <- 4

# Specify the batch size for stage 2, 3, etc. 
# 'batch' should be specified as the total batch size divided 
# the number of arms at start
# 'batch <- 4' here corresponds to 16 subjects in each of the stages 
# 2, 3, ...

batch <- 4

# Specify the minimally required increase in probability of a correct
# decision that is required for continuation
# Note that in the manuscript this threshold gamma is referred to as 
# C/Q

gamma <- 1/1000

# Set the parameters a and b for the Beta(a,b) prior 
# The prior is common for all arms
# We use a uniform prior and set a = b = 1

prior <- c(1,1)

# Generate data for (Ntrials) trials using the create.d function

Ydata <- create.d(resp, maxpt, Ntrials)

# Evaluate all trials simulated using the prdrop function for the
# setting where dropping of arms is allowed
# This is the most computationally intensive part. Under settings considered and using 8 cores on a i5
# processor laptop a call to prdrop required approximely 2 hours computation time (1000 trials)
# Computation increases with smaller stages sizes and smaller gamma

trial.out <- prdrop(Ydata, delta, gamma, burn, Ntrials, batch, prior, maxpt, dropctrl)

# Evaluate all trials simulated using the prnodrop function for the
# setting where dropping of arms is not allowed
# Under settings considered and using 8 cores on a i5
# processor laptop a call to prnodrop required approximely 1 hour 
# computation time (1000 trials)


trialnodrop.out <- prnodrop(Ydata, delta, gamma, burn, Ntrials, batch, prior, maxpt)

# Summarize the results using the pcalc function 
# Function returns the probability of each final decision and a 
# vector containing total sample sizes for all trials

res_drop   <- pcalc(trial.out, prior, delta, nr_dec = 4)
res_nodrop <- pcalc(trialnodrop.out, prior, delta, nr_dec = 4)

# Look at the results

res_drop
res_nodrop

# Look at proportion of each final decision
# res_drop$dec and res_nodrop$dec are vectors (p1,p2,p3,p4) where
# p1: proportion of trials in which experimental arm 1 is selected
# p2: proportion of trials in which experimental arm 2 is selected
# p3: proportion of trials in which experimental arm 3 is selected 
# p4: proportion of trials in which experimental arm 4 is selected 

res_drop$dec
res_nodrop$dec

# Specifically look at proportion of a trials with a correct decision
# (p3)

res_drop$dec[3]
res_nodrop$dec[3]

# res_drop$N and res_nodrop$N contain the trial sizes
# Look at average trial size

mean(res_drop$N)
mean(res_nodrop$N)

# We now show how to reevaluate trials using higher threshold 
# C/Q=1/500 yielding shorter trials

new_gamma <- 1/500

# We reevaluate the simulate trials with and without dropping with 
# the new threshold C/Q = 1/500

trial.out.new.gamma <- reevaluate(trial.out, Ydata,  new_gamma)
trialnodrop.out.new.gamma <- reevaluate(trialnodrop.out, Ydata,  new_gamma)


# We summarize the results for the higher threshold using the 
# pcalc function 

res_drop.new.gamma   <- pcalc(trial.out.new.gamma, prior, delta, nr_dec = 4)
res_nodrop.new.gamma <- pcalc(trialnodrop.out.new.gamma, prior, delta, nr_dec = 4)

# Look at the results for the higher threshold

res_drop.new.gamma 
res_nodrop.new.gamma

# Specifically look at proportion of a trials with a correct decision 
# (p3)

res_drop.new.gamma$dec[3]
res_nodrop.new.gamma$dec[3]



