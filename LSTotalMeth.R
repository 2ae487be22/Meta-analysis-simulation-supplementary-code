### Remove previous variables
rm(list = ls())

StartTimeTotal <- proc.time()

#### Libraries, set seed, set cores ----
library(data.table)
library(doParallel)
library(foreach)
library(doRNG)
library(copula)
library(compiler)
library(metafor)
enableJIT(3)

set.seed(123)

# Number of cores for parallel
num.Cores <- 16
c1 <- makeCluster(num.Cores)

#### Declare variables ----

# Reps = number of repetitions of experiment
Reps = 10000

# k = number of studies in series
Studies = c(3,5,10,30,50,100)

# subj = number of subjects in study
Subj <- list(as.integer(c(100,100)), as.integer(c(20,100)), as.integer(c(250, 1000)), as.numeric(c(4.7, 1.2)))

# sd = study level standard deviation
True.sd = sqrt(1)

# theta = population level log(OR) - this should be considered more purely on the log scale
theta = c(log(0.25), log(0.8), log(1), log(1.25), log(4))

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw)
tau.sq = c(0, 0.008, 0.04, 3.04)

# Frequency of event averaged across 2 arms (before applying change due to theta) = EvFreq
EvFreq = c(0.1, 0.3, 0.5)

# controlProp = proportion of total sample in control arm
controlProp = 0.5

## Boundary of step function on p value, causing severity of publication bias
Severity.boundary <- c(0.05, 0.2)

# Set up strength of publication bias selection
Begg_a <- 0.5
Begg_b <- 3
Begg_c <- -0.3
Begg_sided <- 1

# Set up within study reporting bias
Tested.outcomes <- 5
Sd.split <- 0.8

# Size of per unit bias increase
Bias.multiple <- 0.85

#### Functions ----

findvaluesLogSim <- function(ID_num){
  dummy1 = ID_num %% Reps
  repetitions = ifelse(dummy1 == 0, Reps, dummy1)
  
  intermed <- Studies * Reps
  intermed <- c(0,cumsum(intermed))
  dummy2 = ID_num %% (Reps * sum(Studies))
  numstudies = Studies[(min(which(intermed >= dummy2))-1)]
  
  dummy3 = ID_num %% (Reps * sum(Studies) * length(EvFreq))
  dummy4 = ifelse(dummy3 == 0, length(EvFreq), (dummy3 %/% (Reps * sum(Studies))) + 1)
  eventfreq = EvFreq[dummy4]
  
  dummy5 = ID_num %% (Reps * sum(Studies) * length(EvFreq) * length(tau.sq))
  dummy6 = ifelse(dummy5 == 0, length(tau.sq), (dummy5 %/% (Reps * sum(Studies) * length(EvFreq))) + 1)
  hetero = tau.sq[dummy6]
  
  dummy7 = ID_num %% (Reps * sum(Studies) * length(EvFreq) * length(tau.sq) * length(theta))
  dummy8 = ifelse(dummy7 == 0, length(theta), (dummy7 %/% (Reps * sum(Studies) * length(EvFreq)* length(tau.sq))) + 1)
  truevalue = theta[dummy8]
  
  dummy9 = ID_num %% (Reps * sum(Studies) * length(EvFreq) * length(tau.sq) * length(theta) * length(Subj))
  dummy10 = ifelse(dummy9 == 0, length(Subj), (dummy9 %/% (Reps * sum(Studies) * length(EvFreq)* length(tau.sq) * length(theta))) + 1)
  subjects <- tryCatch(Subj[[dummy10]][1], error = function(e) {try(Subj[dummy10], silent=TRUE)})
  
  return(list(reps = repetitions, subj = subjects, theta = truevalue, tau2 = hetero, numstud = numstudies, evfreq = eventfreq))
}

findIDLogSim <- function(repetitions, subjects, truevalue, hetero, numstudies, eventfreq){
  IDnumber <- integer()
  for (o in 1:numstudies){
    counter.dummy <- as.integer((match(subjects, Subj)-1) * length(EvFreq) * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                  (match(truevalue, theta)-1) * length(tau.sq) * length(EvFreq) * sum(Studies) * Reps +
                                  (match(hetero, tau.sq)-1) * length(EvFreq) * sum(Studies) * Reps +
                                  (match(eventfreq, EvFreq)-1) * sum(Studies) * Reps +
                                  (sum(Studies[0:(match(numstudies, Studies)-1)]) + o -1) * Reps +
                                  repetitions
    )
    IDnumber <- append(IDnumber, counter.dummy)
  }
  return(IDnumber)
}


### Log_Odds_Ratio function
Log_Odds_Ratio <- function(StudySize, Log_O_R, Heterogeneity, Control_Prop, mu){
  StudyLogOR <- rnorm(1, Log_O_R, sqrt(Heterogeneity))
  Group1Size <- as.integer(Control_Prop*StudySize)
  Group2Size <- Group1Size
  EventFreq <- log(mu / (1 - mu))
  Pic <- exp(EventFreq - 0.5*StudyLogOR) / (1 + exp(EventFreq - 0.5*StudyLogOR))
  Group1Out1 <- as.integer(rbinom(1, Group1Size, Pic))
  Group1Out2 <- as.integer(Group1Size - Group1Out1)
  Pit <- exp(EventFreq + 0.5*StudyLogOR)  / (1 + exp(EventFreq + 0.5*StudyLogOR))
  Group2Out1 <- as.integer(rbinom(1, Group2Size, Pit))
  Group2Out2 <- as.integer(Group2Size - Group2Out1)
  if (Group1Out1 == 0 | Group2Out1 == 0 | Group1Out2 == 0 | Group2Out2 == 0){
    Group1Out1 <- Group1Out1 + 0.5
    Group2Out1 <- Group2Out1 + 0.5
    Group1Out2 <- Group1Out2 + 0.5
    Group2Out2 <- Group2Out2 + 0.5
  }
  return(c(Group1Out1, Group1Out2, Group2Out1, Group2Out2, Group1Size, Group2Size))
}

### Function for multiple outcomes

LogOR_mult_out <- function(StudySize, Theta, Heterogeneity, Control_Prop, mu, frac, num.times){
  StudyLogOR <- rnorm(1, Theta, sqrt(Heterogeneity))
  Group1Size <- as.integer(Control_Prop*StudySize)
  Group2Size <- Group1Size
  EventFreq <- log(mu / (1 - mu))
  Pic <- exp(EventFreq - 0.5*StudyLogOR) / (1 + exp(EventFreq - 0.5*StudyLogOR))
  Pit <- exp(EventFreq + 0.5*StudyLogOR)  / (1 + exp(EventFreq + 0.5*StudyLogOR))
  z <- normalCopula(param = frac, dim = num.times)
  Z <- rCopula(Group1Size, z)
  ControlGroup <- qbinom(Z, size=1, prob=Pic)
  CGO1 <- apply(ControlGroup, 2, sum)
  CGO2 <- Group1Size - CGO1
  y <- normalCopula(param = frac, dim = num.times)
  Y <- rCopula(Group1Size, y)
  TreatmentGroup <- qbinom(Y,size=1,prob=Pit)
  TGO1 <- apply(TreatmentGroup, 2, sum)
  TGO2 <- Group2Size - TGO1
  o <- data.table(CGO1, CGO2, TGO1, TGO2)
  o[CGO1 == 0] <- o[CGO1 == 0] + 0.5
  o[CGO2 == 0] <- o[CGO2 == 0] + 0.5
  o[TGO1 == 0] <- o[TGO1 == 0] + 0.5
  o[TGO2 == 0] <- o[TGO2 == 0] + 0.5
  Study.est <- log((o$TGO1/o$TGO2) / (o$CGO1/o$CGO2))
  Study.se <- sqrt(1/o$TGO1 + 1/o$TGO2 + 1/o$CGO1 + 1/o$CGO2)
  Study.p.val <- pnorm(Study.est/Study.se)
  return(c(o[order(Study.p.val)], Group1Size))
}

anyNA <- function(x) {
  i <- 1
  repeat {
    if (is.na(x[i])) return(TRUE)
    i <- i + 1
    if (i > length(x)) return(FALSE)
  }
}

.psort <- function(x,y) {
  
  ### t(apply(xy, 1, sort)) would be okay, but problematic if there are NAs;
  ### either they are removed completely (na.last=NA) or they are always put
  ### first/last (na.last=FALSE/TRUE); but we just want to leave the NAs in
  ### their position!
  
  if (is.null(x) || length(x) == 0) ### need to catch this
    return(NULL)
  
  if (missing(y)) {
    if (is.matrix(x)) {
      xy <- x
    } else {
      xy <- rbind(x) ### in case x is just a vector
    }
  } else {
    xy <- cbind(x,y)
  }
  
  n <- nrow(xy)
  
  for (i in seq_len(n)) {
    if (anyNA(xy[i,]))
      next
    xy[i,] <- sort(xy[i,])
  }
  
  colnames(xy) <- NULL
  
  return(xy)
  
}

mod.hc <- function(object, digits, transf, targs, control, tau2est, ...) {
  
  if (!inherits(object, "rma.uni"))
    stop("Argument 'object' must be an object of class \"rma.uni\".")
  
  if (inherits(object, "rma.ls"))
    stop("Method not yet implemented for objects of class \"rma.ls\". Sorry!")
  
  x <- object
  
  if (!x$int.only)
    stop("Method only applicable for models without moderators.")
  
  if (missing(digits))
    digits <- x$digits
  
  if (missing(transf))
    transf <- FALSE
  
  if (missing(targs))
    targs <- NULL
  
  yi <- x$yi
  vi <- x$vi
  k  <- length(yi)
  
  if (k == 1)
    stop("Stopped because k = 1.")
  
  if (!x$allvipos)
    stop("Cannot use method when one or more sampling variances are non-positive.")
  
  level <- ifelse(x$level > 1, (100-x$level)/100, ifelse(x$level > .5, 1-x$level, x$level))
  
  if (missing(control))
    control <- list()
  
  ###
  
  ### set control parameters for uniroot() and possibly replace with user-defined values
  con <- list(tol=.Machine$double.eps^0.25, maxiter=1000, verbose=FALSE)
  con[pmatch(names(control), names(con))] <- control
  
  ###
  
  ### original code by Henmi & Copas (2012), modified by Michael Dewey, small adjustments
  ### for consistency with other functions in the metafor package by Wolfgang Viechtbauer
  
  wi <- 1/vi ### fixed effects weights
  
  W1 <- sum(wi)
  W2 <- sum(wi^2) / W1
  W3 <- sum(wi^3) / W1
  W4 <- sum(wi^4) / W1
  
  ### fixed-effects estimate of theta
  beta <- sum(wi*yi) / W1
  
  ### Q statistic
  Q <- sum(wi * ((yi - beta)^2))
  
  ### DL estimate of tau^2
  ###### Modified here to take REML tau2
  tau2 <- max(0, tau2est)
  
  vb  <- (tau2 * W2 + 1) / W1 ### estimated Var of b
  se  <- sqrt(vb)             ### estimated SE of b
  VR  <- 1 + tau2 * W2        ### estimated Var of R
  SDR <- sqrt(VR)             ### estimated SD of R
  
  ### conditional mean of Q given R=r
  EQ <- function(r)
    (k - 1) + tau2 * (W1 - W2) + (tau2^2)*((1/VR^2) * (r^2) - 1/VR) * (W3 - W2^2)
  
  ### conditional variance of Q given R=r
  VQ <- function(r) {
    rsq <- r^2
    recipvr2 <- 1 / VR^2
    2 * (k - 1) + 4 * tau2 * (W1 - W2) +
      2 * tau2^2 * (W1*W2 - 2*W3 + W2^2) +
      4 * tau2^2 * (recipvr2 * rsq - 1/VR) * (W3 - W2^2) +
      4 * tau2^3 * (recipvr2 * rsq - 1/VR) * (W4 - 2*W2*W3 + W2^3) +
      2 * tau2^4 * (recipvr2 - 2 * (1/VR^3) * rsq) * (W3 - W2^2)^2
  }
  
  scale <- function(r){VQ(r)/EQ(r)}   ### scale parameter of the gamma distribution
  shape <- function(r){EQ(r)^2/VQ(r)} ### shape parameter of the gamma distribution
  
  ### inverse of f
  finv <- function(f)
    (W1/W2 - 1) * ((f^2) - 1) + (k - 1)
  
  ### equation to be solved
  eqn <- function(x) {
    integrand <- function(r) {
      pgamma(finv(r/x), scale=scale(SDR*r), shape=shape(SDR*r))*dnorm(r)
    }
    integral <- integrate(integrand, lower=x, upper=Inf)$value
    val <- integral - level / 2
    #cat(val, "\n")
    val
  }
  
  t0 <- try(uniroot(eqn, lower=0, upper=2, tol=con$tol, maxiter=con$maxiter))
  
  if (inherits(t0, "try-error"))
    stop("Error in uniroot().")
  
  t0 <- t0$root
  u0 <- SDR * t0 ### (approximate) percentage point for the distribution of U
  
  ###
  
  ci.lb <- beta - u0 * se ### lower CI bound
  ci.ub <- beta + u0 * se ### upper CI bound
  
  beta.rma  <- x$beta
  se.rma    <- x$se
  ci.lb.rma <- x$ci.lb
  ci.ub.rma <- x$ci.ub
  
  ### if requested, apply transformation to yi's and CI bounds
  
  if (is.function(transf)) {
    if (is.null(targs)) {
      beta      <- sapply(beta, transf)
      beta.rma  <- sapply(beta.rma, transf)
      se        <- NA
      se.rma    <- NA
      ci.lb     <- sapply(ci.lb, transf)
      ci.ub     <- sapply(ci.ub, transf)
      ci.lb.rma <- sapply(ci.lb.rma, transf)
      ci.ub.rma <- sapply(ci.ub.rma, transf)
    } else {
      beta      <- sapply(beta, transf, targs)
      beta.rma  <- sapply(beta.rma, transf, targs)
      se        <- NA
      se.rma    <- NA
      ci.lb     <- sapply(ci.lb, transf, targs)
      ci.ub     <- sapply(ci.ub, transf, targs)
      ci.lb.rma <- sapply(ci.lb.rma, transf, targs)
      ci.ub.rma <- sapply(ci.ub.rma, transf, targs)
    }
  }
  
  ### make sure order of intervals is always increasing
  
  tmp <- .psort(ci.lb, ci.ub)
  ci.lb <- tmp[,1]
  ci.ub <- tmp[,2]
  
  tmp <- .psort(ci.lb.rma, ci.ub.rma)
  ci.lb.rma <- tmp[,1]
  ci.ub.rma <- tmp[,2]
  
  ###
  
  res <- list(beta=beta, se=se, ci.lb=ci.lb, ci.ub=ci.ub,
              beta.rma=beta.rma, se.rma=se.rma, ci.lb.rma=ci.lb.rma, ci.ub.rma=ci.ub.rma,
              method="DL", method.rma=x$method, tau2=tau2, tau2.rma=x$tau2, digits=digits)
  
  class(res) <- "hc.rma.uni"
  return(res)
  
}

##### Parallel Loop Simulation ----

StartTime <- proc.time()

registerDoParallel(c1)
set.seed(9120)
LogOR.Sim.Results <- foreach (m = 1:Reps, .combine=rbind, .packages = c("data.table", "copula", "metafor"), 
                              .export = c("Studies", "Subj", "True.sd",
                                          "theta", "tau.sq", "controlProp", "Severity.boundary", "Begg_a", 
                                          "Begg_b", "Begg_sided", "Tested.outcomes", "Sd.split",
                                          "Bias.multiple", "Log_Odds_Ratio", "EvFreq", "LogOR_mult_out" , "Begg_c", "anyNA", ".psort", "mod.hc")
) %dorng% {
  
  # ID different for analysis
  ID = length(tau.sq) * length(Subj) * length(theta) * length(Studies) * length(EvFreq)
  
  LogOR.Sim <- data.table(
    Unique_ID = integer(length = ID),
    FE_Estimate = numeric(length = ID),
    FE_se = numeric(length = ID),
    REML_Estimate = numeric(length = ID),
    REML_se = numeric(length = ID),
    REML_tau2 = numeric(length = ID),
    DL_Estimate = numeric(length = ID),
    DL_se = numeric(length = ID),
    DL_tau2 = numeric(length = ID),
    DL_I2 = numeric(length = ID),
    HC_DL_se = numeric(length = ID),
    HC_DL_CIlb = numeric(length = ID),
    HC_DL_CIub = numeric(length = ID),
    HC_REML_se = numeric(length = ID),
    HC_REML_CIlb = numeric(length = ID),
    HC_REML_CIub = numeric(length = ID),
    KH_REML_CIlb = numeric(length = ID),
    KH_REML_CIub = numeric(length = ID),
    KH_REML_se = numeric(length = ID),
    KH_DL_CIlb = numeric(length = ID),
    KH_DL_CIub = numeric(length = ID),
    KH_DL_se = numeric(length = ID),
    Moreno_Estimate = numeric(length = ID),
    Moreno_se = numeric(length = ID),
    Mult_se = numeric(length = ID),
    Num_exc = integer(length = ID)
  )
  dummy.counter <- 1
  
  for (l in tau.sq){
    
    for (j in EvFreq){
      
      for (i in Subj){
        
        for(k in theta){
          
          for (n in Studies){
          
          Study_estimate <- numeric(length = n)
          Study_sd <- numeric(length = n)
          ni <- integer(length = n)
          excluded <- integer(length = n)
          
          for (o in 1:n){
            
            #Select sample size
            if (is.integer(i[1]) == TRUE){
              Study_patientnumber <- as.integer( (runif(1, sqrt(i[1]), sqrt(i[2])))^2 )
            } else {
              Study_patientnumber <- round(rlnorm(1, meanlog = i[1], sdlog = i[2]) + 2)
            }
            
            Number.of.biases <- rbinom(1, 2, 1/(Study_patientnumber^0.3))
            
            x <- Log_Odds_Ratio(Study_patientnumber, log(exp(k) * (Bias.multiple^Number.of.biases)), l, controlProp, j)
            x.N <- 2*x[5]
            
            study.mean <- log((x[3]/x[4]) / (x[1]/x[2]))
            study.se <- sqrt(1/x[1] + 1/x[2] + 1/x[3] + 1/x[4])
            
            Study_estimate[o] <- study.mean
            Study_sd[o] <- study.se
            excluded[o] <- ifelse(((x[1] == 0.5 & x[3] == 0.5) | (x[2] == 0.5 & x[4] == 0.5)), 1, 0 )
            ni[o] <- x.N
            
          }
          
          counter <- as.integer((apply(sapply(Subj, function(vec) {i %in% vec}), 1, which.max)[1]-1) * length(EvFreq) * length(theta) * length(tau.sq) * length(Studies) * Reps + 
                                  (match(k, theta)-1) * length(tau.sq) * length(EvFreq) * length(Studies) * Reps +
                                  (match(l, tau.sq)-1) * length(EvFreq) * length(Studies) * Reps +
                                  (match(j, EvFreq)-1) * length(Studies) * Reps +
                                  (match(n, Studies)-1) * Reps +
                                  m)
          
          
          #find.temp.values <- findIDLogSim(m, i, k, l, n, j)
          
          ### Temporary data.table
          #temp.data <- LogOR.Simulation[J(find.temp.values)]
          temp.data <- data.table(Study_estimate, Study_sd, ni, excluded)
          
          #           temp.data <- temp.data[ !((Study_G1O1 == 0.5 & Study_G2O1 == 0.5) | ( ((Study_n/2 + 1) - Study_G1O1) == 0.5 & ((Study_n/2 + 1) - Study_G2O1) == 0.5 ))]
          #           Excluded.studies <- n - dim(temp.data)[1]
          
          temp.data <- temp.data[excluded != 1]
          Excluded.studies <- n - dim(temp.data)[1]
          
          if (dim(temp.data)[1] < 2){
            LogOR.Sim[dummy.counter, `:=` (Unique_ID = counter,
                                           FE_Estimate = NA,
                                           FE_se = NA,
                                           REML_Estimate = NA,
                                           REML_se = NA,
                                           REML_tau2 = NA,
                                           DL_Estimate = NA,
                                           DL_se = NA,
                                           DL_tau2 = NA,
                                           DL_I2 = NA,
                                           HC_DL_se = NA,
                                           HC_DL_CIlb = NA,
                                           HC_DL_CIub = NA,
                                           HC_REML_se = NA,
                                           HC_REML_CIlb = NA,
                                           HC_REML_CIub = NA,
                                           KH_REML_CIlb = NA,
                                           KH_REML_CIub = NA,
                                           KH_REML_se = NA,
                                           KH_DL_CIlb = NA,
                                           KH_DL_CIub = NA,
                                           KH_DL_se = NA,
                                           Moreno_Estimate = NA,
                                           Moreno_se = NA,
                                           Mult_se = NA,
                                           Num_exc = Excluded.studies
            )]
            dummy.counter <- dummy.counter + 1
          } else {
            #Fixed and random effects
            ma.fe <- tryCatch({
              rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2 , method = "FE")
            },
            error = function(e){
              return(list(b = NA,  se = NA))
            },
            warning = function(w){
              return(list(list(b = NA,  se = NA)))
            }
            )
            
            ma.reml <- tryCatch({
              rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2  , method = "REML", control = list(stepadj = 0.5))
            },
            error = function(e){
              return(list(b = NA, tau2 = NA, se = NA))
            },
            warning = function(w){
              return(list(list(b = NA, tau2 = NA, se = NA)))
            }
            )
            
            ma.DL <- tryCatch({
              rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2  , method = "DL")
            },error = function(e){
              return(list(b = NA, tau2 = NA, se = NA, I2 = NA))
            },warning = function(w){
              return(list(list(b = NA, tau2 = NA, se = NA, I2 = NA)))
            })
            
            # Henmi & Copas
            
            ma.hc.DL <- tryCatch({
              hc(ma.DL)
            },error = function(e){
              return(list(se = NA, ci.lb = NA, ci.ub = NA))
            },warning = function(w){
              return(list(list(se = NA, ci.lb = NA, ci.ub = NA)))
            })
            
            ma.hc.REML <- tryCatch({
              mod.hc(ma.reml, tau2est = ma.reml$tau2)
            },error = function(e){
              return(list(se = NA, ci.lb = NA, ci.ub = NA))
            },warning = function(w){
              return(list(list(se = NA, ci.lb = NA, ci.ub = NA)))
            })
            
            # Knapp Hartung
            
            ma.reml.kh <- tryCatch({
              rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2  , method = "REML", knha = TRUE, control = list(stepadj = 0.5))
            },
            error = function(e){
              return(list(b = NA, tau2 = NA, se = NA, ci.lb = NA, ci.ub = NA))
            },
            warning = function(w){
              return(list(b = NA, tau2 = NA, se = NA, ci.lb = NA, ci.ub = NA))
            }
            )
            
            ma.DL.kh <- tryCatch({rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2  , method = "DL", knha = TRUE)}
                                 , error = function(e){return(list(se = NA, ci.lb = NA, ci.ub = NA))
                                 })
            
            ## Doi
            # estimate is equal to fixed effect, as are weights
            doi.var.DL <- tryCatch(sum( ( as.vector(weights(ma.fe, type = "diagonal")/100)^2 ) * (temp.data$Study_sd^2 + ma.DL$tau2) ), error=function(err) NA)
            doi.var.REML <- tryCatch(sum( ( as.vector(weights(ma.fe, type = "diagonal")/100)^2 ) * (temp.data$Study_sd^2 + ma.reml$tau2) ), error=function(err) NA)
            
            ## Moreno (?D-var) - not exactly clear which implementation is being used is likely equation 2a
            moreno.est <- tryCatch({ma.moren <- regtest(ma.fe , predictor = "vi", model = "lm")
            c(ma.moren$fit[[5]][1],ma.moren$fit[[5]][3])
            }, error=function(err) c(NA,NA))
            
            ## Mawdesley
            
            ma.mult <- tryCatch({
              mawd.lm <- lm(temp.data$Study_estimate ~ 1, weights = 1/(temp.data$Study_sd^2))
              sm.mawd.lm <- summary(mawd.lm)
              ifelse(mean(sm.mawd.lm$residuals^2) < 1, phi.est  <- 1, phi.est <- mean(sm.mawd.lm$residuals^2))
              rma.uni(temp.data$Study_estimate, temp.data$Study_sd^2 * phi.est , method = "FE")}
              , error = function(e){return(list(se = NA))
              })
            
            
            LogOR.Sim[dummy.counter, `:=` (Unique_ID = counter,
                                           FE_Estimate = ma.fe[[1]][1],
                                           FE_se = ma.fe$se,
                                           REML_Estimate = ma.reml$b[1],
                                           REML_se = ma.reml$se,
                                           REML_tau2 = ma.reml$tau2,
                                           DL_Estimate = ma.DL[[1]][1],
                                           DL_se = ma.DL$se,
                                           DL_tau2 = ma.DL$tau2,
                                           DL_I2 = ma.DL$I2,
                                           HC_DL_se = ma.hc.DL$se,
                                           HC_DL_CIlb = ma.hc.DL$ci.lb,
                                           HC_DL_CIub = ma.hc.DL$ci.ub,
                                           HC_REML_se = ma.hc.REML$se,
                                           HC_REML_CIlb = ma.hc.REML$ci.lb,
                                           HC_REML_CIub = ma.hc.REML$ci.ub,
                                           KH_REML_CIlb = ma.reml.kh$ci.lb,
                                           KH_REML_CIub = ma.reml.kh$ci.ub,
                                           KH_REML_se = ma.reml.kh$se,
                                           KH_DL_CIlb = ma.DL.kh$ci.lb,
                                           KH_DL_CIub = ma.DL.kh$ci.ub,
                                           KH_DL_se = ma.DL.kh$se,
                                           IVHet_DL_var = doi.var.DL,
                                           IVHet_REML_var = doi.var.REML,
                                           Moreno_Estimate = moreno.est[1],
                                           Moreno_se = moreno.est[2],
                                           Mult_se = ma.mult$se,
                                           Num_exc = Excluded.studies
            )]
            
            dummy.counter <- dummy.counter + 1
          }
          }
        }
      }
    }
  }
  LogOR.Sim
}

stopCluster(c1)

#### Need to then sort final table and add values for rep number, rep subj, rep theta, rep tau2, rep numstudies
LogOR.Sim.Results <- LogOR.Sim.Results[order(Unique_ID)]


##### Need to re append values - specific to analysis
ID =  length(Subj) * length(theta) * length(tau.sq) * length(EvFreq) * Reps * length(Studies)

Subj <- c(100, 20, 250, 4.7)

LogOR.Sim.Results$Rep_Number =  rep(1:Reps, times = ID/Reps)
LogOR.Sim.Results$Rep_NumStudies = rep(rep(Studies, each = Reps), times = ID/(Reps*length(Studies)))
LogOR.Sim.Results$Rep_ev_freq = rep(rep(EvFreq, each = Reps * length(Studies)), times = ID/(Reps*length(Studies)*length(EvFreq)))
LogOR.Sim.Results$Rep_tau.sq = rep(rep(tau.sq, each = Reps * length(Studies)*length(EvFreq)), times = ID/(Reps*length(Studies)*length(tau.sq)*length(EvFreq)))
LogOR.Sim.Results$Rep_theta = rep( rep(theta, each = Reps * length(Studies) * length(tau.sq)*length(EvFreq)), times = length(Subj))
LogOR.Sim.Results$Rep_Subj = rep(Subj, each = ID / length(Subj))

TimeTakenSim <- proc.time() - StartTime

#write.csv(LogOR.Sim.Results, file = "LSTotalV2.csv")
saveRDS(LogOR.Sim.Results, file = "LSTotalMethRDS")

#### Timings ----

TimeTakenTotal <- proc.time() - StartTimeTotal

TimeTakenSim
TimeTakenTotal

#### Checking values ----

sum(is.na(LogOR.Sim.Results))
LogOR.Sim.Results[is.na(LogOR.Sim.Results$REML_Est)]
summary(LogOR.Sim.Results)
