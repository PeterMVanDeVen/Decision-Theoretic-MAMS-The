# Supplementary R code for "Bayesian adaptive decision-theoretic designs
# for multi-arm multi stage clinical trials"


# This example contains all functions required for the setting with
# multiple experimental arms.
#
# Following functions are included:
#
# create.d  : Simulates trial data;
# prnodrop  : Evaluates characteristics for trials without dropping of
#             arms;
# pcalc     : Summarizes trial characteristics;
# reevaluate: Reevaluates characteristics of trials obtained by prdrop
#             and prnodrop for higher thresholds C/Q;
# evalnodrop: Calculates proportion of correct decisions in case of
#		      continuation (in trials without dropping of arms)
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
    Y <- array(dim = c(Nperarm, K, Ntrials))
    for (l in 1:Ntrials) {
        for (k in 1:K) {
            Y[,k,l] <- rbinom(n = Nperarm, size = 1, prob = resp[k])
        }
    }
    return(Y)
}


# prnodrop
#
# Function prnodrop evaluates simulated trials using the decision-theoretic criteria
# for continuation when dropping of arms is not allowed.
#
# Required input:
# Y       : Vector with data observed so far in a single trial
#         : (number of responders/successes in each arm)
# n1      : Vector with total number of patients for which outcomes
#         : have been observed in each arm of the trial
# n2      : Vector with planned number of patients to be included in
#         : each arm in the additional stage of the trial
# prior   : Parameters c(a,b) for Beta(a,b) prior
# gamma   : Minimally required increase (C/Q) in probability of
#         : correct decision required for continuation
#
# Uses functions: evalnodrop, Ttheta, Tdata
#
# Returns: list containing for each simulated trial the number of
#          patients allocated to each arm and the number of
#          responders per arm


prnodrop <- function(
    Y,
    gamma,
    burn = 0,
    batch,
    prior,
    maxpt,
    n_cores,
    AR = FALSE
){
    cl <- makeCluster(min(detectCores(), n_cores)) #change the 8 to your number of CPU cores
    registerDoParallel(cl)
    Ntrials <- dim(Y)[3]
    LL <- array(dim = c(Ntrials,(2*(dim(Y)[2]) + 3)))
    print(Sys.time())

    LL <- foreach(i = 1:Ntrials, .combine = rbind, .export = c('evalnodrop', 'Ttheta', 'Tdata')) %dopar% {
        G <- 1000
        n <- batch  # batch size at start (PER ARM)
        arms <- dim(Y)[2]
        y <- Y[,,i]
        TR <- 1
        dat <- rep(0, arms)
        len <- rep(max(burn - n, 0),arms)
        len2 = rep(n, arms);
        ben <- NULL
        a <- prior[1]
        b <- prior[2]
        ns <- matrix(ncol = arms, nrow = 0)
        ys <- matrix(ncol = arms, nrow = 0)
        first <- TRUE
        while ((gamma == 0 | TR > gamma) & (sum(len) + arms*n <= maxpt) != 0) {
            len <- len + len2
            old_len <- len2
            dat <- NULL
            for (j in 1:length(len)) {
                dat <- c(dat, sum(y[1:len[j], j]))
            }
            if (AR == TRUE) {
                theta <- matrix(rbeta(G*arms, shape1 = a + dat, shape2 = b + len - dat), ncol = arms, byrow = T)
                Tt <- Ttheta(theta)
                ntotal = n*arms
                theta_ThompsonSampling <- matrix(rbeta(ntotal*arms, shape1 = a + dat, shape2 = b + len - dat), ncol = arms, byrow = T)
                allocationAR <- apply(theta_ThompsonSampling, 1, which.max)
                len2 <- table(c(allocationAR,(1:arms))) - rep(1,arms)
            }
            #ntotal = n*arms
            #sqrtprob=sqrt((table(c(Tt,1:arms))-1)/length(Tt))   # this is used specifically for selecting the best
            #probs = sqrtprob/sum(sqrtprob)
            #len2 = t(rmultinom(1, ntotal, probs))}
            pr <- evalnodrop(y = dat, n1 = len, n2 = len2, prior, gamma, AR)
            TR <- max(pr$ben.cont - pr$ben.stop)
            ben <- c(ben, TR)
            if (first) {
                nstage <- as.vector(len)
                ystage <- dat
                first <- FALSE
            } else {
                nstage <- as.vector(old_len)
                ystage <- dat - apply(ys,2,sum)
            }
            ns <- rbind(ns,nstage)
            ys <- rbind(ys, ystage)
        }
        c(
            as.vector(len),
            as.vector(dat),
            list(benefit = ben, n = ns, y = ys)
        )
    }
    print(Sys.time())
    stopCluster(cl)
    rm(cl)
    return(LL)
}


prAr <- function(
    Y,
    gamma,
    burn = 0,
    batch,
    prior,
    maxpt,
    n_cores
){
    cl <- makeCluster(min(detectCores(), n_cores)) #change the 8 to your number of CPU cores
    registerDoParallel(cl)
    Ntrials <- dim(Y)[3]
    LL <- array(dim = c(Ntrials,(2*(dim(Y)[2]) + 3)))
    print(Sys.time())

    LL <- foreach(i = 1:Ntrials, .combine = rbind, .export = c('evalFullAr', 'Ttheta')) %dopar% {
        G <- 1000
        n <- batch  # batch size at start (PER ARM)
        arms <- dim(Y)[2]
        y <- Y[,,i]
        TR <- 0
        dat <- rep(0, arms)
        len <- rep(max(burn - n, 0),arms)
        len2 = rep(n, arms);
        ben <- NULL
        a <- prior[1]
        b <- prior[2]
        ns <- matrix(ncol = arms, nrow = 0)
        ys <- matrix(ncol = arms, nrow = 0)
        first <- TRUE
        while ((gamma == 0 | TR < gamma) & (sum(len) + arms*n <= maxpt) != 0) {
            len <- len + len2
            old_len <- len2
            dat <- NULL
            for (j in 1:length(len)) {
                dat <- c(dat, sum(y[1:len[j], j]))
            }
            theta <- matrix(rbeta(G*arms, shape1 = a + dat, shape2 = b + len - dat), ncol = arms, byrow = T)
            Tt <- Ttheta(theta)
            ntotal = n*arms
            theta_ThompsonSampling <- matrix(rbeta(ntotal*arms, shape1 = a + dat, shape2 = b + len - dat), ncol = arms, byrow = T)
            allocationAR <- apply(theta_ThompsonSampling, 1, which.max)
            len2 <- table(c(allocationAR,(1:arms))) - rep(1,arms)

            TR <- evalFullAr(y = dat, n = len, prior)
            ben <- c(ben, TR)
            if (first) {
                nstage <- as.vector(len)
                ystage <- dat
                first <- FALSE
            } else {
                nstage <- as.vector(old_len)
                ystage <- dat - apply(ys,2,sum)
            }
            ns <- rbind(ns,nstage)
            ys <- rbind(ys, ystage)
        }
        c(
            as.vector(len),
            as.vector(dat),
            list(benefit = ben, n = ns, y = ys)
        )
    }
    print(Sys.time())
    stopCluster(cl)
    rm(cl)
    return(LL)
}


# pcalc
#
# Function pcalc computes the probability of each final decision being the correct decision and
# the total trial sizes using the output of the prdrop and prnodrop functions.
#
# Required input:
# LL      : List generated as output of prdrop, prnodrop or reevaluate
#         : functions
# prior   : Parameters c(a,b) for Beta(a,b) prior
# nr_dec  : Number of possible final decisions
#
# Uses function: Tdata
#
# Returns: list with vector 'dec' containing proportion of trials in
#          which each of the final decisions is made and vector 'N'
#          with total sample sizes of the trials
## dimension of Yrec determined assuming that $y part of output (fixed for prnodrop, not for prdrop)


pcalc <- function(LL, prior, nr_dec){
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
        Pr[j] <- Tdata(y = y[j,], n = n[j,], prior)
    }
    tbl <- table(factor(Pr, levels = 1:nr_dec))/L
    prop_allocated = apply(n,2,sum)/sum(n)
    prop_succes = sum(y)/sum(n)
    return (list(dec = tbl, prop_allocated = round(prop_allocated,3), prop_succes= round(prop_succes,3), N = apply(n,1,sum)))}



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
        nstage <- matrix(unlist(trialres[i,Q-1]), ncol=arms, byrow=FALSE)
        if(length(which(ben <= newgamma))==0) {st<-length(ben)} else {st<-min(which(ben <= newgamma))}
        # with new gamma stopped after st stages
        if(st>1) {npt<-apply(nstage[1:st,], 2, sum)} else {npt<-nstage[1,]}
        dat<-c()
        for(j in 1:arms) {
            dat[j]=sum(Y[1:npt[j],j,i])}
        Yrec[i,] <- c(npt, dat)
    } #endfor
    return(Yrec)}

reevaluateAr <- function(trialres, Y, newgamma){
    Q <- dim(trialres)[2]
    L <- dim(trialres)[1]
    arms <- (Q-3)/2
    Yrec <- array(dim = c(L,Q-3))
    for(i in 1:L){
        ben <- unlist(trialres[i,Q-2])
        nstage <- matrix(unlist(trialres[i,Q-1]), ncol=arms, byrow=FALSE)
        if(length(which(ben > newgamma))==0) {st<-length(ben)} else {st<-min(which(ben > newgamma))}
        # with new gamma stopped after st stages
        if(st>1) {npt<-apply(nstage[1:st,], 2, sum)} else {npt<-nstage[1,]}
        dat<-c()
        for(j in 1:arms) {
            dat[j]=sum(Y[1:npt[j],j,i])}
        Yrec[i,] <- c(npt, dat)
    } #endfor
    return(Yrec)}

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
# gamma   : Minimally required increase (C/Q) in probability of
#         : correct decision required for continuation
#
# Uses functions: Ttheta, Tdata

evalnodrop <- function(y, n1, n2, prior, gamma, AR = FALSE) {

    K <- length(y)
    # Number of draws per arm
    if (gamma > 0) {
        G <- ceiling(2.5*(1/gamma))
    } else {
        G <- 10000
    }

    a <- prior[1]
    b <- prior[2]

    theta <-
        matrix(
            rbeta(
                G*K,
                shape1 = a + y,
                shape2 = b + n1 - y
            ),
            ncol = K,
            byrow = TRUE
        )

    Tt <- Ttheta(theta)

    Ty <- Tdata(y = y, n = n1, prior)

    if (AR == FALSE) {
        Ynew <- matrix(
            rbinom(
                n = G*K,
                size = n2,
                prob = t(theta)
            ),
            ncol = K,
            byrow = TRUE
        )
        Tynew <- apply(t(t(Ynew) + y), 1, Tdata, n1 + n2, prior)
    } else {
        n2snew <- matrix(ncol = G, nrow = K)
        for (j in 1:G) {
            theta_ThompsonSampling <- matrix(rbeta(sum(n2)*K, shape1 = a + y, shape2 = b + n - y), ncol = K, byrow = TRUE)
            allocationAR <- apply(theta_ThompsonSampling, 1, which.max)
            len2 = table(c(allocationAR, (1:K))) - rep(1,K)
            n2snew[, j] <- len2
        }
        #ntotal = sum(n2)
        #sqrtprob=sqrt((table(c(Tt,1:K))-1)/length(Tt))   # this is used specifically for selecting the best
        #probs = sqrtprob/sum(sqrtprob)
        #n2snew = rmultinom(G, ntotal, probs)
        Ynew <- matrix(rbinom(n = G*K, size = as.vector(n2snew), prob = t(theta)), ncol = K, byrow = T)
        ynew = t(t(Ynew) + y)
        nnew = matrix(rep(n1,G), byrow = TRUE, ncol = K) + t(n2snew)
        Tynew = c()
        for (i in 1:G) {
            Tynew[i] <- Tdata(ynew[i,], nnew[i,], prior)
        }
    }

    b0 <- sum(Ty == Tt)/G
    b1 <- sum(Tynew == Tt)/G

    return(list(ben.stop = b0, ben.cont = 1))

}

evalFullAr <- function(y, n, prior) {

    if (length(y) != length(n)) {
        cat("ERROR: dimension of y and n must coincide")
        return(NULL)
    }
    a <- prior[1]
    b <- prior[2]
    G <- 10000
    K <- length(y)

    theta <-
        matrix(
            rbeta(
                G*K,
                shape1 = a + y,
                shape2 = b + n - y
            ),
            ncol = K,
            byrow = TRUE
        )

    best_arm <-
        factor(
            apply(theta, 1, function(x) which.max(x)),
            levels = 1:K
        )
    probs_arm <- as.vector(table(best_arm) / G)

    best <- max(probs_arm)
    return(best)

}

# Tdata
#
# Function Tdata determines the loss-minimizing decision based on the data. The function given here is for the
# loss function that appears as equation (2) in the manuscript and can be used when selecting from a set
# of multiple experimental arms the one that has highest posterior probability of maximizing the response rate.
#
# Required input:
# y      : Vector with data observed so far in a single trial
#        : (number of responders/successes in each arm)
# n      : Vector with total number of patients for which outcomes
#        : have been observed in each arm of the trial
# prior  : Parameters c(a,b) for Beta(a,b) prior
#
# Returns: 1 if experimental arm 1 has highest posterior probability
#            of having the highest response rate
#          2 if experimental arm 2 has highest posterior probability
#            of having the highest response rate
#          3 if experimental arm 3 has highest posterior probability
#            of having the highest response rate
#          etc

Tdata <- function(y, n, prior) {
    if (length(y) != length(n)) {
        cat("ERROR: dimension of y and n must coincide")
        return(NULL)
    }
    a <- prior[1]
    b <- prior[2]
    G <- 1000
    K <- length(y)

    theta <-
        matrix(
            rbeta(
                G*K,
                shape1 = a + y,
                shape2 = b + n - y
            ),
            ncol = K,
            byrow = TRUE
        )

    best_arm <-
        factor(
            apply(theta, 1, function(x) which.max(x)),
            levels = 1:K
        )
    probs_arm <- as.vector(table(best_arm) / G)

    best <- max(probs_arm)
    indices_max <- which(probs_arm == best)

    if (length(indices_max) > 1) {
        tm <- sample(indices_max,1)
    } else {
        tm <- indices_max[1]
    }

    return(tm)

}

# Ttheta
#
# Function Ttheta determines the correct decision given a true parameter vector theta. The function given here
# is for loss function that appears as equation (2) in the manuscript and can be used when selecting among a set
# of multiple experimental arms the one that maximized the response rate.
#
# Required input:
# theta  : Vector containing response/success probabilities for the
#        : three arms
#
# Returns: 1 if experimental arm 1 has highest response rate
#          2 if experimental arm 2 has highest response rate
#          3 if experimental arm 3 has highest response rate
#          etc

Ttheta <- function(theta) {
    if (is.null(dim(theta))) {
        pos <- which.max(theta)
    } else {
        pos <- apply(theta, 1, which.max)
    }
    return(pos)
}


