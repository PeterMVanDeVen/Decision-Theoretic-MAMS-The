# Supplementary R code for "Bayesian adaptive decision-theoretic designs
# for multi-arm multi stage clinical trials" by
# Bassi, Berkhof, de Jong, van de Ven
#
# Correspondence/enquiries: p.vandeven@amsterdamumc.nl

# This script requires that the R script 'Functions for selecting the best from
# three, four or five experimental arms.R' was previously run and packages snow,
# doSNOW and doParallel for parallel computing were installed.
# The required R script can be downloaded from
# https://github.com/PeterMVanDeVen/Decision-Theoretic-MAMS-trials

# This example illustrates the use of the following functions for the
# setting with three experimental arms:
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

resp <- c(0.2,0.8,0.9,0.2,0.8)

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
# 'burn <- 4' here corresponds to 12 subjects in the first stage

burn <- 4

# Specify the batch size for stage 2, 3, etc.
# 'batch' should be specified as the total batch size divided
# the number of arms at start
# 'batch <- 3' here corresponds to 12 subjects in each of the stages
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

Ydata2 <- create.d(resp, maxpt, Ntrials)


# Evaluate all trials simulated using the prnodrop function for the
# setting where dropping of arms is not allowed
# Under settings considered and using 8 cores on a i5
# processor laptop a call to prnodrop required approximely 30 minutes
# computation time (1000 trials)

trial.out <- prnodrop(
    Y = Ydata2,
    gamma = gamma,
    burn = burn,
    batch = batch,
    prior = prior,
    maxpt = maxpt,
    n_cores = 4,
    AR = FALSE
)

# trialdrop.out <- prdrop(Ydata2, delta, gamma, burn, Ntrials, batch, prior, maxpt)
# trialnodrop.out <- prnodrop(Ydata2, delta, gamma, burn, Ntrials, batch, prior, maxpt)
# trialnodropAR.out2 <- prnodrop(Ydata2, delta, gamma, burn, Ntrials, batch, prior, maxpt, AR=TRUE)

# Summarize the results using the pcalc function
# Function returns the probability of each final decision and a
# vector containing total sample sizes for all trials

res_drop <- pcalc(trialdrop.out, prior, delta, nr_dec = 5)
res_nodrop <- pcalc(trialnodrop.out, prior, delta, nr_dec = 5)
res_nodropAR2 <- pcalc(trialnodropAR.out2, prior, delta, nr_dec = 5)

# Look at the results
res_drop
res_drop$dec[3]
mean(res_drop$N)

res_nodrop
res_nodrop$dec[3]
mean(res_nodrop$N)

res_nodropAR2
res_nodropAR2$dec[3]
mean(res_nodropAR2$N)

psucces<-c()
psuccesAR<-c()
psuccesdrop<-c()
pcorrect<-c()
pcorrectAR<-c()
pcorrectdrop<-c()
meanN<-c()
meanNAR<-c()
meanNdrop<-c()
prop_allocated_best<-c()
prop_allocated_bestAR<-c()
prop_allocated_bestdrop<-c()


for(gamma in (100:1)/1000){
out<-reevaluate(trialnodrop.out, Ydata2, gamma)
res<-pcalc(out, prior, delta, nr_dec=5)

psucces <- c(psucces, res$prop_succes)
pcorrect <- c(pcorrect, res$dec[3])
prop_allocated_best<-c(prop_allocated_best, res$prop_allocated[3])
meanN <- c(meanN, mean(res$N))

out<-reevaluate(trialnodropAR.out2, Ydata2, gamma)
res<-pcalc(out, prior, delta, nr_dec=5)

psuccesAR <- c(psuccesAR, res$prop_succes)
pcorrectAR <- c(pcorrectAR, res$dec[3])
prop_allocated_bestAR<-c(prop_allocated_bestAR, res$prop_allocated[3])
meanNAR <- c(meanNAR, mean(res$N))

out<-reevaluate(trialdrop.out, Ydata2, gamma)
res<-pcalc(out, prior, delta, nr_dec=5)

psuccesdrop <- c(psuccesdrop, res$prop_succes)
pcorrectdrop <- c(pcorrectdrop, res$dec[3])
prop_allocated_bestdrop<-c(prop_allocated_bestdrop, res$prop_allocated[3])
meanNdrop <- c(meanNdrop, mean(res$N))




}

par(mfrow=c(1,2))

plot(meanNAR, pcorrectAR, type='l', col='blue', xlim=c(0, 250), ylim=c(0,1),ylab='proportion of correct decisions', xlab='average trial size')
lines(meanN, pcorrect, col='red')
lines(meanNdrop, pcorrectdrop, col='green')

plot(meanNAR, prop_allocated_bestAR, type='l', col='blue', xlim=c(0, 250), ylim=c(0,1),ylab='proportion allocated to best arm', xlab='average trial size')
lines(meanN, prop_allocated_best, col='red')
lines(meanNdrop, prop_allocated_bestdrop, col='green')




resp <- c(0.2,0.6,0.8,0.2,0.6)



########

trialnodrop.out <- prnodrop(Ydata2, delta, gamma, burn, Ntrials, batch, prior, maxpt)

res_nodrop   <- pcalc(trialnodrop.out, prior, delta, nr_dec = 5)

res_nodrop$dec[3]
mean(res_nodrop$N)

res<-apply(trialnodropAR.out2[1,]$n,2, sum)
for(i in 2:5000){
res<-res+apply(trialnodropAR.out2[i,]$n,2, sum)}
res<-res/5000
round(res,1)
sum(round(res,1))

N<-c();
for(i in 30:10){
trialAR_reev<-reevaluate(trialnodropAR.out2, Ydata2, i/10000)
res_nodropAR2_reev <- pcalc(trialAR_reev, prior, delta, nr_dec = 5)

N<-c(N,mean(res_nodropAR2_reev$N))}

plot(30:10/10000, N)

i=22
trialAR_reev<-reevaluate(trialnodropAR.out2, Ydata2, i/10000)
res_nodropAR2_reev <- pcalc(trialAR_reev, prior, delta, nr_dec = 5)

mean(res_nodropAR2_reev$N)

res_nodropAR2_reev$dec
