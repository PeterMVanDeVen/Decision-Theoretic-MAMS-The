# Supplementary R code for "Bayesian adaptive decision-theoretic designs
# for multi-arm multi stage clinical trials" 



# This script contains all the functions required for the setting where
# two experimental arms are compared to control.
#
# The following functions are included:
#
# create.d  : Simulates trial data;
# prdrop    : Evaluates characteristics for trials when arms can be 
#             dropped;
# prnodrop  : Evaluates characteristics for trials without dropping of 
#             arms;
# pcalc     : Summarizes trial characteristics;
# reevaluate: Reevaluates characteristics of trials obtained by prdrop
#             and prnodrop for higher thresholds C/Q;  
# evaldrop  : Calculates proportion of correct decisions in case of 
#		: continuation (when arms can be dropped)
# evalnodrop: Calculates proportion of correct decisions in case of 
#		: continuation (in trials without dropping of arms)
# prsinglestage: Summarizes trial characteristics of a single-stage
#		: trials with final decision based on minimization of 
#		: loss
# Tdata     : Determines the correct decision given a true parameter 
#           : vector theta.
# Ttheta    : Determines the the loss-minimizing decision based on the 
#           : data given a true parameter vector theta.



# create.d
#
# Function create.d simulates data for [Ntrials] trials with [Nperarm] patients per arm and response vector [resp]
#
# Required input:
# resp    : Vector containing the true response/success probabilities
#         : in each of the arms
# Nperarm : Number of patients per arm
# Ntrials : Number of trials to be simulated (simulation size)

create.d <- function(resp, Nperarm, Ntrials){
  K <- length(resp) #no. arms
  Y <- array(dim=c(Nperarm, K, Ntrials))
  for(l in 1:Ntrials){
    for(k in 1:K){
      Y[,k,l] <- rbinom(n = Nperarm, size = 1, prob = resp[k])
    }
  }
  return(Y)
}



# prdrop
#
# Function prdrop evaluates simulated trials using the decision-theoretic criteria for continuation and dropping of arms.
#
# Required input:
# Y       : Trial data created using the create.d function
# delta   : Margin for superiority of experimental arms over the 
#         : control arm 
#         : Specify as 0 for symmetric case with experimental 
#         : treatments only
# gamma   : Minimally required increase (C/Q) in probability of 
#         : correct decision required for continuation
# burn    : Batch size per arm for stage 1 
#         : If specified as 0 then batch size (batch) specified for 
#         : subsequent stages is also used for stages 1. 
# batch   : Batch size per arm for stage 2,3,4... 
#         : Total batch size is determined as number of arms at the  
#				  : start multiplied by value of batch 
#			    : In case arms have been dropped patients are divided 
#         : equally over the remaining arms
# Ntrials : Number of simulated trials in Y
# prior   : Parameters c(a,b) for Beta(a,b) prior 
# maxpt   : Maximum for total number of patients that can be 
#         : included in the trial (cap)
# dropctrl: Specify as TRUE if dropping of the control arm is allowed
#         : or if only experimental arms are considered
#
# Uses functions: evaldrop, Ttheta, Tdata
#
# Returns: list containing for each simulated trial the number of 
#          patients allocated to each arm and the number of
#          responders per arm

prdrop <- function(Y, delta = 0, gamma, burn = 0, Ntrials, batch, prior, maxpt, dropctrl = TRUE){
  cl<-makeCluster(8) #change the 8 to your number of CPU cores
  registerDoParallel(cl)
  LL <- array(dim=c(Ntrials,(2*(dim(Y)[2]) + 3)))
  print(Sys.time())
  LL <- foreach(i=1:Ntrials, .combine = rbind, .export=c('evaldrop', 'Ttheta', 'Tdata')) %dopar% {
    n <- batch  # batch size at start (PER ARM)
    arms <- dim(Y)[2]
    y <- Y[,,i]   
    TR <- 1
    dat <- rep(0, arms)
    len <- rep(max(burn-n, 0),arms)
    active=rep(TRUE, arms)
    drop <- 0
    ben <- NULL
    isdr <- NULL
    ns<-c();
    first=TRUE;
    cdt <- (sum(len)+arms*n <= maxpt)
    while((gamma==0 | TR > gamma) & cdt){
      if(drop[1] == 0){
        len <- len+n
          dat <- NULL
        for(j in 1:length(len)){
          dat <- c(dat, sum(y[1:len[j], j]))
        }
        pr <- evaldrop(y = dat, n1 = len, n2 = rep(n, arms), prior, delta, gamma, dropctrl)
        TR <- max(pr$ben.cont - pr$ben.stop, na.rm = T)
        drop <- which.max(pr$ben.cont - pr$ben.stop)
        drop <- drop - 1
        ben <- c(ben, TR)
        isdr <- c(isdr, drop)
        if(first) {nstage <- len; first=FALSE} else {nstage <- rep(n, arms)};
        ns<-rbind(ns, nstage)
        active[drop]=FALSE
      } else {
        old_len=len;
        n_active=sum(active)
        active_arms = (1:arms)[active]
        len[active_arms] <- len[active_arms]+ floor((n*arms)/(arms-length(drop)))
        n_remain = n*arms - floor((n*arms)/(arms-length(drop)))*n_active
        if(n_remain > 0) {arms_with_min = active_arms[which(min(len[active_arms])==len[active_arms])];
        n_active_arms_with_minimum = length(arms_with_min)
          if(n_remain == n_active_arms_with_minimum) {len[arms_with_min]=len[arms_with_min]+1} else 
            if(n_remain < n_active_arms_with_minimum) {randsamp = sample(arms_with_min, n_remain); len[randsamp]=len[randsamp]+1} else
              {len[arms_with_min]=len[arms_with_min]+1; 
               n_remain = n_remain - n_active_arms_with_minimum; 
               randsamp = sample(active_arms, n_remain); len[randsamp]=len[randsamp]+1}}
        dat <- NULL
        for(j in 1:length(len)){
          dat <- c(dat, sum(y[1:len[j], j]))
        }
        nf <- len-old_len
        pr <- evaldrop(dat, len, nf, prior, delta, gamma, dropctrl)
        TR <- max(pr$ben.cont - pr$ben.stop, na.rm = T)
        dp <- which.max(pr$ben.cont - pr$ben.stop)
        nstage <- nf;
        ns <- rbind(ns, nstage)
        if(dp != 1){
          drop <- c(drop,dp-1)
          isdr <- c(isdr,dp-1);
          active[drop]=FALSE
        } else(isdr <- c(isdr, 0))
        ben <- c(ben, TR)
      }
      cdt <- (sum(len)+arms*n <= maxpt)
    }
    c(len,dat, list(benefit = ben, drop=isdr, n = ns))
  }
  print(Sys.time())
  stopCluster(cl)
  rm(cl)
  return(LL)
 }

# prnodrop
# 
# Function prnodrop evaluates simulated trials using the decision-theoretic criteria 
# for continuation when dropping of arms is not allowed.
#
# Required input:
# y       : Vector with data observed so far in a single trial
#         : (number of responders/successes in each arm)
# n1      : Vector with total number of patients for which outcomes
#         : have been observed in each arm of the trial
# n2      : Vector with planned number of patients to be included in 
#         : each arm in the additional stage of the trial
# prior   : Parameters c(a,b) for Beta(a,b) prior 
# delta   : Margin for superiority of experimental arms over the 
#         : control arm 
#         : Specify as 0 for symmetric case with experimental 
#         : treatments only
# gamma   : Minimally required increase (C/Q) in probability of 
#         : correct decision required for continuation
#
# Uses functions: evalnodrop, Ttheta, Tdata
#
# Returns: list containing for each simulated trial the number of 
#          patients allocated to each arm and the number of
#          responders per arm

prnodrop <- function(Y, delta = 0, gamma, burn = 0, Ntrials, batch, prior, maxpt){
  cl<-makeCluster(8) #change the 8 to your number of CPU cores
  registerDoParallel(cl)
  LL <- array(dim=c(Ntrials,(2*(dim(Y)[2]) + 3)))
  print(Sys.time())
  LL <- foreach(i=1:Ntrials, .combine = rbind, .export=c('evalnodrop', 'Ttheta', 'Tdata')) %dopar% {
    n <- batch  # batch size at start (PER ARM)
    arms <- dim(Y)[2]
    y <- Y[,,i]
    TR <- 1
    dat <- rep(0, arms)
    len <- rep(max(burn-n, 0),arms)
    ben <- NULL
    isdr <- NULL
    ns<-c();
    first=TRUE;
    cdt <- (sum(len)+arms*n <= maxpt)
    while((gamma==0 | TR > gamma) & (sum(len) + arms*n <= maxpt) != 0){
      len <- len+n
      dat <- NULL
      for(j in 1:length(len)){
        dat <- c(dat, sum(y[1:len[j], j]))
      }
      pr <- evalnodrop(y = dat, n1 = len, n2 = rep(n, arms), prior, delta, gamma)
      TR <- max(pr$ben.cont - pr$ben.stop)
      ben <- c(ben, TR)
      isdr <- c(isdr, 0)      
      if(first) {nstage <- len; first=FALSE} else {nstage=rep(n, arms)};ns <- rbind(ns,nstage);
    }  
    c(len,dat, list(benefit = ben, drop=isdr, n = ns))
  }
  print(Sys.time())
  stopCluster(cl)
  rm(cl)
  return(LL)}


# pcalc
#
# Function pcalc computes the probability of each final decision being the correct decision and 
# the total trial sizes using the output of the prdrop and prnodrop functions.
#
# Required input:
# LL      : List generated as output of prdrop, prnodrop or reevaluate 
#         : functions
# prior   : Parameters c(a,b) for Beta(a,b) prior 
# delta   : Margin for superiority of experimental arms over the 
#         : control arm as used in prdrop and prnodrop functions
# nr_dec  : Number of possible final decisions
#
# Uses function: Tdata
#
# Returns: list with vector 'dec' containing proportion of trials in 
#          which each of the final decisions is made and vector 'N' 
#          with total sample sizes of the trials

pcalc <- function(LL, prior, delta = 0, nr_dec){
	if(!is.numeric(LL)){
		Yrec <- matrix(as.numeric(LL[,1:(dim(LL)[2]-3)]), ncol = (dim(LL)[2]-3)) 
	} else{
		Yrec <- LL
	}   
	L <- dim(Yrec)[1]
	K <- dim(Yrec)[2]
	y <- Yrec[,((K/2)+1):K]
	n <- Yrec[,1:(K/2)]
	Pr <- array()
	for(j in 1:L){
		Pr[j] <- Tdata(y = y[j,], n = n[j,], prior, delta)
	}
	tbl <- table(factor(Pr, levels = 1:nr_dec))/L
	return (list(dec = tbl, N = apply(n,1,sum)))} 



# reevaluate
#
# Function reevaluate reevaluates final decisions and total trial sizes for simulated trials 
# for higher treshold gamma (C/Q) using output of prdrop or prnodrop function.
#
# Required input:
# trialres: List generated as output of prdrop of prnodrop functions
# Y       : Simulated trial data used to obtain trialres
# burn    : Batch size per arm for stage 1 
#         : If specified as 0 then batch size (batch) specified for 
#         : subsequent stages is also used for stages 1. 
# batch   : Batch size per arm for stage 2,3,4... 
#         : Total batch size is determined as number of arms at the  
#				  : start multiplied by value of batch 
#			    : In case arms have been dropped patients are divided 
#         : equally over the remaining arms
# newgamma: New minimally required increase (C/Q) in probability  
#         : of correct decision required for continuation
#
# Returns : List containing for each simulated trial the number of 
#           patients allocated to each arm and the number of
#           responders per arm

reevaluate <- function(trialres, Y, newgamma){
 Q <- dim(trialres)[2]
 L <- dim(trialres)[1]
   arms <- (Q-3)/2
   Yrec <- array(dim = c(L,Q-3))
   for(i in 1:L){
      ben <- unlist(trialres[i,Q-2])
      drp <- unlist(trialres[i,Q-1])
      nstage <- matrix(unlist(trialres[i,Q]), ncol=arms, byrow=FALSE)
      if(length(which(ben <= newgamma))==0) {st<-length(ben)} else {st<-min(which(ben <= newgamma))}
      # with new gamma stopped after st stages
      if(st>1) {npt<-apply(nstage[1:st,], 2, sum)} else {npt<-nstage[1,]}  
      dat<-c()
      for(j in 1:arms) {
       dat[j]=sum(Y[1:npt[j],j,i])}  
       Yrec[i,] <- c(npt, dat)
    } #endfor
return(Yrec)}



# evaldrop
#
# Function evaldrop computes the proportion of a correct decision in case the trial is continued for
# an additional stage in settings where dropping of arms is allowed.
#
# Required input:
# y       : Vector with data observed so far in a single trial
#         : (number of responders/successes in each arm)
# n1      : Vector with total number of patients for which outcomes
#         : have been observed in each arm of the trial
# n2      : Vector with planned number of patients to be included in 
#         : each arm in the additional stage of the trial
# prior   : Parameters c(a,b) for Beta(a,b) prior 
# delta   : Margin for superiority of experimental arms over the 
#         : control arm 
#         : Specify as 0 for symmetric case with experimental 
#         : treatments only
# gamma   : Minimally required increase (C/Q) in probability of 
#         : correct decision required for continuation
# dropctrl: Specify as TRUE if dropping of the control arm is allowed
#         : or if only experimental arms are considered
#
# Uses functions: Ttheta, Tdata

evaldrop <- function(y, n1, n2, prior, delta = 0, gamma, dropctrl = TRUE){
   if(delta < 0){ 
     cat('ERROR: delta must be a non-negative number')
     return(NULL)}
   ndr <- sum(n2 == 0)
   dr <- which(n2 == 0)
   J <- length(y)          # total number of arms
   active <- rep(TRUE, J)
   active[dr]=FALSE
   active_arms = (1:J)[active]
   Jnew <- length(y) - ndr # number of active arms
   if(gamma>0) {G <- ceiling(2.5*(1/gamma))} else {G<-10000} # number of draws per arm
   a <- prior[1]
   b <- prior[2]
   benefit <- NA
   theta <- matrix(rbeta(G*J, shape1 = a + y, shape2 = b + n1 - y), ncol = J, byrow = T)
   Tt <- apply(theta, 1, Ttheta, delta)
   Ty <- Tdata(y = y, n = n1, prior, delta = delta)
   n2new=rep(0,J)
   n2new[active_arms] = floor(sum(n2)/Jnew)
   n_remain = sum(n2)-sum(n2new)
   if(n_remain > 0) {len=n1+n2new;arms_with_min = active_arms[which(min(len[active_arms])==len[active_arms])];
        n_active_arms_with_minimum = length(arms_with_min)
          if(n_remain == n_active_arms_with_minimum) {len[arms_with_min]=len[arms_with_min]+1; n2new=len-n1} else 
            if(n_remain < n_active_arms_with_minimum) {randsamp = sample(arms_with_min, n_remain); len[randsamp]=len[randsamp]+1; n2new=len-n1} else
              {len[arms_with_min]=len[arms_with_min]+1; 
               n_remain = n_remain - n_active_arms_with_minimum; 
               randsamp = sample(active_arms, n_remain); len[randsamp]=len[randsamp]+1}; n2new=len-n1}
   Ynew <- matrix(rbinom(n = G*J, size = n2new, prob = t(theta)), ncol = J, byrow = T)
   Tynew <- apply( t(t(Ynew) + y), 1, Tdata, n1+n2new, prior, delta)
   b0 <- sum(Ty == Tt)/G
   benefit[1] <- sum(Tynew == Tt)/G
   if(Jnew > 2){
       if(ndr == 0){
         ifelse(dropctrl, idx <- (1:J), idx <- (2:J))  
            } else{ifelse(dropctrl, idx <- (1:J)[-(dr)], idx <- (2:J)[-(dr-1)])}
      for(k in idx){
         n2drop=rep(0,J)
         n2drop[active_arms] = floor(sum(n2)/(Jnew-1))
         n2drop[k]=0
         n_remain = sum(n2)-sum(n2drop) 
         active_arms[which(min(len[active_arms])==len[active_arms])]
          if(n_remain > 0) {len=n1+n2drop;active_after_drop=active; active_after_drop[k]=FALSE; active_arms_after_drop = (1:J)[active_after_drop];
                arms_not_dropped_with_min = active_arms_after_drop[which(min(len[active_arms_after_drop])==len[active_arms_after_drop])];
                n_active_arms_with_minimum = length(arms_not_dropped_with_min)
                if(n_remain == n_active_arms_with_minimum) {len[arms_not_dropped_with_min]=len[arms_not_dropped_with_min]+1; n2drop=len-n1} else 
                if(n_remain < n_active_arms_with_minimum) {randsamp = sample(arms_not_dropped_with_min, n_remain); len[randsamp]=len[randsamp]+1; n2drop=len-n1} else
              {len[arms_not_dropped_with_min]=len[arms_not_dropped_with_min]+1; 
               n_remain = n_remain - n_active_arms_with_minimum; 
               randsamp = sample(active_arms_after_drop, n_remain); len[randsamp]=len[randsamp]+1}; n2drop=len-n1}
        Ynew1d <- matrix(rbinom(n = G*J, size = n2drop, prob = t(theta)), ncol = J, byrow = T)
        Tynew <- apply( t(t(Ynew1d) + y), 1, Tdata, n1 + n2drop, prior,  delta)
        benefit[k+1] <- sum(Tynew == Tt)/G
      }
   }
return(list(ben.stop = b0, ben.cont = benefit))}




# evalnodrop
#
# Function evalnodrop computes the proportion of a correct decision in case the trial is continued for an 
# additional stage in settings where dropping of arms is not allowed.
#
# Required input:
# y       : Vector with data observed so far in a single trial
#         : (number of responders/successes in each arm)
# n1      : Vector with total number of patients for which outcomes
#         : have been observed in each arm of the trial
# n2      : Vector with planned number of patients to be included in 
#         : each arm in the additional stage of the trial
# prior   : Parameters c(a,b) for Beta(a,b) prior 
# delta   : Margin for superiority of experimental arms over the 
#         : control arm 
#         : Specify as 0 for symmetric case with experimental 
#         : treatments only
# gamma   : Minimally required increase (C/Q) in probability of 
#         : correct decision required for continuation
#
# Uses functions: Ttheta, Tdata

evalnodrop <- function(y, n1, n2, prior, delta = 0, gamma){
  if(delta < 0){ 
    cat("ERROR: delta must be a non-negative number \n")
    return(NULL)
  }
  J <- length(y)
  if(gamma>0) {G <- ceiling(2.5*(1/gamma))} else {G<-10000} # number of draws per arm
  a <- prior[1]
  b <- prior[2]
  benefit <- NA
  theta <- matrix(rbeta(G*J, shape1 = a + y, shape2 = b + n1 - y), ncol = J, byrow = T)
  Tt <- apply(theta, 1, Ttheta, delta)
  Ty <- Tdata(y = y, n = n1, prior, delta = delta)
  Ynew <- matrix(rbinom(n = G*J, size = n2, prob = t(theta)), ncol = J, byrow = T)
  Tynew <- apply( t(t(Ynew) + y), 1, Tdata, n1+n2, prior, delta)
  b0 <- sum(Ty == Tt)/G
  benefit[1] <- sum(Tynew == Tt)/G
  return(list(ben.stop = b0, ben.cont = benefit))}

# prsinglestage
#
# Function prsinglestage determines the correct decision given a true parameter vector theta. The function
# given here is for the loss function that appears as equation (3) in the manuscript and can be used 
# when comparing two experimental arms to a control arm.
#
# Y       : Trial data created using the create.d function
# delta   : Margin for superiority of experimental arms over the 
#         : control arm 
#         : Specify as 0 for symmetric case with experimental 
#         : treatments only
# burn    : Sample size per arm 
# Ntrials : Number of simulated trials in Y
# prior   : Parameters c(a,b) for Beta(a,b) prior 
# nr_dec  : Number of possible final decisions
#
# Uses function: Tdata
#
# Returns: list with vector 'dec' containing proportion of trials in 
#          which each of the final decisions is made and 'N' 
#          which is the total sample sizes of the trials

prsinglestage <- function(Y, delta = 0, burn, Ntrials, prior, nr_dec){
  cl<-makeCluster(2) #change the 2 to your number of CPU cores
  registerDoParallel(cl)
  LL <- array(dim=c(Ntrials,1))
  print(Sys.time())
  LL <- foreach(i=1:Ntrials, .combine = rbind, .export=c('Tdata')) %dopar% {
    arms <- dim(Y)[2]
    n <- rep(burn, arms) # batch size at start (PER ARM)
    y <- apply(Y[,,i][1:burn,],2,sum)
    Ty = Tdata(y, n, prior, delta)
    c(decisions=Ty)
    }
   print(Sys.time())
   stopCluster(cl)
   rm(cl)
   tbl <- table(factor(LL[,1], levels = 1:nr_dec))/Ntrials
   return(list(dec = tbl, N = arms*burn))
}

# Tdata
#
# Function Tdata determines the loss-minimizing decision based on the data. The function given here is for the loss function
# that appears as equation (3) in the manuscript and can be used when comparing two experimental arms to a control arm.
#
# Required input:
# y      : Vector with data observed so far in a single trial
#        : (number of responders/successes in each arm)
# n      : Vector with total number of patients for which outcomes
#        : have been observed in each arm of the trial
# prior  : Parameters c(a,b) for Beta(a,b) prior 
# delta  : Margin for superiority of experimental arms over the 
#        : control arm 
# Returns: 1 if the expected loss minimizing decision is declaring no 
#            experimental superior to control 
#          2 if declaring only the first experimental arm (arm 2) 
#            superior to control minimizes expected loss
#          3 if declaring only the second experimental arm (arm 3)  
#            superior to control minimizes expected loss
#          4 if declaring both experimental arms superior 
#            to control minimizes expected loss

Tdata <- function(y, n, prior, delta = 0){
  if(length(y) != length(n)){
    cat("ERROR: dimension of y and n must coincide")
    return(NULL)
  }
  a <- prior[1]
  b <- prior[2]
  G <- 1000
  theta <- matrix(rbeta(G*3, shape1 = a + y, shape2 = b + n - y), ncol = 3, byrow = T)
  prob <- c(mean(theta[,2]-theta[,1] <= delta & theta[,3]-theta[,1] <= delta), 
            mean(theta[,2]-theta[,1] > delta & theta[,3]-theta[,1] <= delta),
            mean(theta[,2]-theta[,1] <= delta & theta[,3]-theta[,1] > delta),
            mean(theta[,2]-theta[,1] > delta & theta[,3]-theta[,1] > delta))
  best <- max(prob); indices_max = which(prob==best)
  if (length(indices_max)>1) {tm=sample(indices_max,1)} else {tm = indices_max[1]} 
  return(tm)}



# Ttheta
#
# Function Ttheta determines the correct decision given a true parameter vector theta. The function
# given here is for the loss function that appears as equation (3) in the manuscript and can be used 
# when comparing two experimental arms to a control arm.
#
# Required input:
# theta  : Vector containing response/success probabilities for the  
#        : three arms 
# delta  : Margin for superiority of experimental arms over the 
#        : control arm 
# Returns: 1 if no experimental arm is superior to control 
#          2 if only first experimental (arm 2) is superior to control
#          3 if only second experimental (arm 3) is superior to 
#            control
#          4 if both experimental arms are superior to control

Ttheta <- function(theta, delta = 0){
  if(delta < 0){
    cat("ERROR: delta must be a non-negative number \n")
    return(NULL)
  }
  tm <- theta - theta[1]
  pos <- 1 + sum(which(tm[-1] > delta)) 
  return(pos) 
}



