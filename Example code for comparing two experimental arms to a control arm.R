# Supplementary R code for "Bayesian adaptive decision-theoretic designs
# for multi-arm multi stage clinical trials" by
# Bassi, Berkhof, de Jong, van de Ven
#
# Correspondence/enquiries: p.vandeven@amsterdamumc.nl


# The script requires that the R script 'Functions for comparing two experimental
# arms to a control arm.R' was previously run and packages snow, doSNOW and
# doParallel for parallel computing were installed.
# The required R script can be downloaded from 
# https://github.com/PeterMVanDeVen/Decision-Theoretic-MAMS-trials


# The script illustrates the use of the following functions for the 
# setting with two experimental arms and a control arm: 
#
# create.d      : simulates trial data;
# prsinglestage : evaluation characterics for single-stage trials with
#               : final decision based on minimization of expected loss
# prdrop        : evaluates characteristics for trials when arms can be 
#                 dropped;
# prnodrop      : evaluates characteristics for trials without dropping of 
#                 arms;
# pcalc         : summarizes trial characteristics
# reevaluate    : Reevaluates characteristics of trials obtained by prdrop
#                 and prnodrop for higher thresholds C/Q  


# Packages snow, doSNOW an doParallel are required for parallel 
# computing 

library(snow)
library(doSNOW)
library(doParallel)

# Set random seed

set.seed(123456)

# Specify the number of trials to be simulated (simulation size)

Ntrials <- 1000

# Specify the response rate vector under which the frequentist 
# properties of the trial are to be evaluated
# First entry corresponds to control arm 
# We consider both a null scenario (resp0) and an alternative
# scenario (respA)

resp0 <- c(0.20,0.20,0.20)
resp1 <- c(0.20,0.20,0.35)

# Specify the cap (= maximum number of subjects included in the trial)
# In this example the cap is set at a high number that is not reached 
# for the response rate vector and threshold (C/Q) used for stopping

maxpt <- 400

# Specify the margin delta above which experimental treatments are
# considered superior to control

delta <- 0.1

# Specify whether early dropping of the control arm is allowed

dropctrl <- FALSE

# Specify batch size for stage 1 
# 'burn' should be specified as the total batch size divided by
# the number of arms at start.
# 'burn <- 48' here corresponds to 3*48=144 subjects in the first stage

burn <- 48

# Specify the batch size for stage 2, 3, etc. 
# 'batch' should be specified as the total batch size divided 
# the number of arms at start
# 'batch <- 12' here corresponds to 36 subjects in each of the stages 
# 2, 3, ...

batch <- 12

# Specify the minimally required increase in probability of a correct
# decision that is required for continuation
# Note that in the manuscript this threshold gamma is referred to as 
# C/Q

gamma <- 0.0015

# Set the parameters a and b for the Beta(a,b) prior 
# The prior is common for all arms
# We use a uniform prior and set a = b = 1

prior <- c(1,1)

# Generate data for (Ntrials) trials using the create.d function
# for null scenario (Ydata0) and alternative scenario (Ydata1)

Ydata0 <- create.d(resp0, maxpt, Ntrials)
Ydata1 <- create.d(resp1, maxpt, Ntrials)

#  Look at type I and type II error of single stage design with 396 patients

results_single_stage0 <- prsinglestage(Ydata0, delta, 396/3, Ntrials, prior, nr_dec = 4)
results_single_stage1 <- prsinglestage(Ydata1, delta, 396/3, Ntrials, prior, nr_dec = 4)

# Look at proportion of each final decision
# results_single_stage0$dec and results_single_stage1$dec are vectors (p1,p2,p3,p4) where
# p1: proportion of trials in which control arm is selected
# p2: proportion of trials in which only first experimental arm is
#     declared superior to control
# p3: proportion of trials in which only second experimental arm is 
#     declared superior to control
# p4: proportion of trials in which both experimental arms are 
#     declared superior to control

results_single_stage0$dec
results_single_stage1$dec

# Calculate type I and type II error probabilities for the single-stage design with maximum sample size

TypeIerror_singlestage = 1-results_single_stage0$dec[1]; TypeIerror_singlestage
TypeIIerror_singlestage = 1-sum(results_single_stage1$dec[3:4]);TypeIIerror_singlestage

# Evaluate all trials simulated using the prdrop function for the
# setting where dropping of arms is allowed (separately for both
# scenarios)
#
# This is the most computationally intensive part. Under settings considered and using 8 cores on a i5
# processor laptop each call to prdrop required 20-30 minutes computation time (1000 trials)
# Computation increases with smaller stages sizes and smaller gamma

trial.out0 <- prdrop(Ydata0, delta, gamma, burn, Ntrials, batch, prior, maxpt, dropctrl)
trial.out1 <- prdrop(Ydata1, delta, gamma, burn, Ntrials, batch, prior, maxpt, dropctrl)

# Summarize the results using the pcalc function 
# Function returns the probability of each final decision and a 
# vector containing total sample sizes for all trials

res0 <- pcalc(trial.out0, prior, delta, nr_dec = 4)
res1 <- pcalc(trial.out1, prior, delta, nr_dec = 4)

# Look at proportion of each final decision
# res0$dec and res1$dec are vectors (p1,p2,p3,p4) where
# p1: proportion of trials in which control arm is selected
# p2: proportion of trials in which only first experimental arm is
#     declared superior to control
# p3: proportion of trials in which only second experimental arm is 
#     declared superior to control
# p4: proportion of trials in which both experimental arms are 
#     declared superior to control

res0$dec
res1$dec

# Obtain familywise type I error probability

FWER = 1-res0$dec[1]; FWER

# Obtain proportion of trials declaring only second experimental arm 
# superior to control (correct decision as in simulation II)

p.correct = res1$dec[3]; p.correct

# Obtain proportion of trials in which second experimental arm is   
# declared superior to control (power in practical example)

power = res1$dec[3]+res1$dec[4]; power
typeIIerror = 1-power; typeIIerror

# res0$N and res1$N contains trial sizes
# Look at average trial size

mean(res0$N)
mean(res1$N)

# We now show how to reevaluate trials using higher threshold 
# C/Q=0.01 

new_gamma <- 0.01

# We reevaluate the simulated trials both under the null and alternative 
# scenario using the new threshold C/Q = 0.01

trial.out0.new.gamma <- reevaluate(trial.out0, Ydata0, new_gamma)
trial.out1.new.gamma <- reevaluate(trial.out1, Ydata1, new_gamma)

# We summarize the results for the higher threshold using the 
# pcalc function 

res0.new.gamma <- pcalc(trial.out0.new.gamma, prior, delta, nr_dec = 4)
res1.new.gamma <- pcalc(trial.out1.new.gamma, prior, delta, nr_dec = 4)

# Look at the results for the higher threshold
# Proportion of final decisions

res0.new.gamma$dec 
res1.new.gamma$dec

# Mean trial sizes

mean(res0.new.gamma$N)
mean(res1.new.gamma$N)

# Obtain familywise type I error probability

FWER.new.gamma = 1-res0.new.gamma$dec[1]; FWER.new.gamma

# Obtain proportion of trials declaring only second experimental arm 
# superior to control (correct decision as in simulation II)

p.correct.new.gamma = res1.new.gamma$dec[3]; p.correct.new.gamma

# Obtain proportion of trials in which second experimental arm is   
# declared superior to control (power in practical example)

power.new.gamma = res1.new.gamma$dec[3]+res1.new.gamma$dec[4]; power.new.gamma
typeIIerror.new.gamma = 1-power.new.gamma; typeIIerror.new.gamma 




## Next we illustrate the setting without considering dropping of arms

# The function prnodrop can be used in the setting where dropping of arms is not considered
#
# The prnodrop function requires less computational time than the prdrop as at each stage
# only a single option for continuation (continue with all arms) needs to be considered 
# Under settings considered and using 8 cores on a i5 processor laptop 
# each call to prdrop required 10-20 minutes computation time (1000 trials)
# Computation increases with smaller stages sizes and smaller gamma

trial_nodrop.out0 <- prnodrop(Ydata0, delta, gamma, burn, Ntrials, batch, prior, maxpt)
trial_nodrop.out1 <- prnodrop(Ydata1, delta, gamma, burn, Ntrials, batch, prior, maxpt)

# Again we summarize the results using the pcalc function 
# This function returns the probability of each final decision and a 
# vector containing total sample sizes for all trials

res_nodrop0 <- pcalc(trial_nodrop.out0, prior, delta, nr_dec = 4)
res_nodrop1 <- pcalc(trial_nodrop.out1, prior, delta, nr_dec = 4)

# Look at proportion of each final decision

res_nodrop0$dec
res_nodrop1$dec

# Look at mean trial size 

mean(res_nodrop0$N)
mean(res_nodrop1$N)

# Obtain familywise type I error probability

FWER = 1-res_nodrop0$dec[1]; FWER

# Obtain proportion of trials declaring only second experimental arm 
# superior to control (correct decision as in simulation II)

p.correct = res_nodrop1$dec[3]; p.correct

# Obtain proportion of trials in which second experimental arm is   
# declared superior to control (power in practical example)

power = res_nodrop1$dec[3]+res_nodrop1$dec[4]; power
typeIIerror = 1-power; typeIIerror

# As type I error was found to be below the desired level, we can increase the threshold
# We use the reevaluate function to do this

new_gamma <- 0.014

trial_nodrop.out0.new.gamma <- reevaluate(trial_nodrop.out0, Ydata0,  new_gamma)
trial_nodrop.out1.new.gamma <- reevaluate(trial_nodrop.out1, Ydata1,  new_gamma)

# We again summarize the results for the higher threshold using the 
# pcalc function 

res_nodrop0.new.gamma <- pcalc(trial_nodrop.out0.new.gamma, prior, delta, nr_dec = 4)
res_nodrop1.new.gamma <- pcalc(trial_nodrop.out1.new.gamma, prior, delta, nr_dec = 4)

# Look at proportion of decisions for the higher threshold

res_nodrop0.new.gamma$dec 
res_nodrop1.new.gamma$dec

# Compute mean trial sizes

mean(res_nodrop0.new.gamma$N)
mean(res_nodrop1.new.gamma$N)

# Obtain proportion of trials declaring only second experimental arm 
# superior to control (correct decision as in simulation II)

p.correct = res_nodrop1.new.gamma$dec[3]; p.correct

# Obtain proportion of trials in which second experimental arm is   
# declared superior to control (power in practical example)

power = res_nodrop1.new.gamma$dec[3]+res_nodrop1.new.gamma$dec[4]; power
typeIIerror = 1-power; typeIIerror


