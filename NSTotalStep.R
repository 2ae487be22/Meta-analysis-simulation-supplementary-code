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

# subj = number of subjects in study, likely to be distributed
Subj <- list(as.integer(c(60,60)), as.integer(c(20,100)), as.integer(c(250, 1000)), as.numeric(c(4.2, 1.1)))

# sd = study level standard deviation
True.sd = sqrt(1)

# theta = population level mean - need good sense of range for SMD
theta = c( -0.76,  -0.12,  0, 0.12, 0.76)

# tau.sq = between studies variance (can be squared due to sqrt() in normal draw), ?to be distributed
tau.sq = c(0, 0.005, 0.022, 1.676)

# controlProp = proportion of total sample in control arm
controlProp = 0.5

## Boundary of step function on p value, causing severity of publication bias
Severity.boundary <- c(0.05, 0.2)

# Set up strength of publication bias selection IF STILL USING
Begg_a <- 0.5
Begg_b <- 3
Begg_c <- -0.3
Begg_sided <- 1

# Set up within study reporting bias - this is now one sided
Tested.outcomes <- 5
Sd.split <- 0.6

# Size of per unit bias increase
Bias.multiple <- c(0, log(0.85)/(-1.81), log(0.7225)/(-1.81))

#### Functions ----

findvaluesNSim <- function(ID_num){
  dummy1 = ID_num %% Reps
  repetitions = ifelse(dummy1 == 0, Reps, dummy1)
  
  intermed <- Studies * Reps
  intermed <- c(0,cumsum(intermed))
  dummy2 = ID_num %% (Reps * sum(Studies))
  numstudies = Studies[(min(which(intermed >= dummy2))-1)]
  
  dummy5 = ID_num %% (Reps * sum(Studies) * length(tau.sq))
  dummy6 = ifelse(dummy5 == 0, length(tau.sq), (dummy5 %/% (Reps * sum(Studies))) + 1)
  hetero = tau.sq[dummy6]
  
  dummy7 = ID_num %% (Reps * sum(Studies) * length(tau.sq) * length(theta))
  dummy8 = ifelse(dummy7 == 0, length(theta), (dummy7 %/% (Reps * sum(Studies) * length(tau.sq))) + 1)
  truevalue = theta[dummy8]
  
  dummy9 = ID_num %% (Reps * sum(Studies) * length(tau.sq) * length(theta) * length(Subj))
  dummy10 = ifelse(dummy9 == 0, length(Subj), (dummy9 %/% (Reps * sum(Studies) * length(tau.sq) * length(theta))) + 1)
  subjects <- tryCatch(Subj[[dummy10]][1], error = function(e) {try(Subj[dummy10], silent=TRUE)})
  
  return(list(reps = repetitions, subj = subjects, theta = truevalue, tau2 = hetero, numstud = numstudies))
}

findIDNSim <- function(repetitions, subjects, truevalue, hetero, numstudies){
  IDnumber <- integer()
  for (o in 1:numstudies){
    counter.dummy <- as.integer((match(subjects, Subj)-1)  * length(theta) * length(tau.sq) * sum(Studies) * Reps + 
                                  (match(truevalue, theta)-1) * length(tau.sq) * sum(Studies) * Reps +
                                  (match(hetero, tau.sq)-1) * sum(Studies) * Reps +
                                  (sum(Studies[0:(match(numstudies, Studies)-1)]) + o -1) * Reps +
                                  repetitions
    )
    IDnumber <- append(IDnumber, counter.dummy)
  }
  return(IDnumber)
}


### Unstandardised mean differrence function

UMD <- function(StudySize, Theta, Heterogeneity, Control_Prop, sd){
  StudyUMD <- rnorm(1, Theta, sqrt(Heterogeneity))
  Group1Size <- as.integer(Control_Prop*StudySize)
  Group2Size <- Group1Size
  ControlGroup <- rnorm(Group1Size, -StudyUMD/2, sd)
  TreatmentGroup <- rnorm(Group2Size, StudyUMD/2, sd)
  Studymean <- mean(TreatmentGroup) - mean(ControlGroup)
  Studysd <- sqrt( (var(ControlGroup) * (Group1Size - 1) + var(TreatmentGroup) * (Group2Size-1))/ (Group1Size + Group2Size -2) * (1/Group1Size + 1/Group2Size))
  return(c(Studymean, Studysd))
}

### UMD function with multiple outcome bias with frac being sd in first level, num.times = number of outcomes simulated
# outputs vectors ordered by p-val

UMD.mult.out <- function(StudySize, Theta, Heterogeneity, Control_Prop, total.sd, frac, num.times){
  StudyUMD <- rnorm(1, Theta, sqrt(Heterogeneity))
  Group1Size <- as.integer(Control_Prop*StudySize)
  Group2Size <- Group1Size
  z <- normalCopula(param = frac, dim = num.times)
  Z <- rCopula(Group1Size, z)
  ControlGroup <- qnorm(Z, mean = -StudyUMD/2, sd = total.sd)
  y <- normalCopula(param = frac, dim = num.times)
  Y <- rCopula(Group1Size, y)
  TreatmentGroup <- qnorm(Y, mean = StudyUMD/2, sd = total.sd)
  Studymean <- apply(TreatmentGroup,2,mean) - apply(ControlGroup, 2, mean)
  Studysd <- sqrt( (apply(TreatmentGroup, 2, var) * (Group1Size - 1) + apply(TreatmentGroup, 2, var) * (Group2Size-1))/ (Group1Size + Group2Size -2) * (1/Group1Size + 1/Group2Size))
  Begg_p <- pnorm(Studymean/Studysd)
  return(list(Studymean[order(Begg_p)], Studysd[order(Begg_p)]))
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

#### Parallel Sim Loop ----

StartTime <- proc.time()

registerDoParallel(c1)
set.seed(3465)
Normal.Sim.Results <- foreach (m = 1:Reps, .combine=rbind, .packages = c("data.table", "copula", "metafor"), 
                               .export = c("Studies", "Subj", "True.sd",
                                           "theta", "tau.sq", "controlProp", "UMD", "Severity.boundary", "Begg_a", 
                                           "Begg_b", "Begg_sided", "Tested.outcomes", "Sd.split",
                                           "Bias.multiple", "UMD.mult.out", "Begg_c", "anyNA", ".psort", "mod.hc")
) %dorng% {
  
  # ID different for analysis
  ID = length(tau.sq) * length(Subj) * length(theta) * length(Studies)
  
  Normal.Sim <- data.table(
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
    Mult_se = numeric(length = ID)
  )
  dummy.counter <- 1
  
  for (l in tau.sq){
    
    for (n in Studies){
      
      for (i in Subj){
        
        for(k in theta){
        
        Study_estimate <- numeric(length = n)
        Study_sd <- numeric(length = n)
        ni <- integer(length = n)
        
        for (o in 1:n){
          
          #Select sample size
          if (is.integer(i[1]) == TRUE){
            #Study_patientnumber <- i[1]
            Study_patientnumber <- as.integer( (runif(1, sqrt(i[1]), sqrt(i[2])))^2 )
          } else {
            Study_patientnumber <- round(rlnorm(1, meanlog = i[1], sdlog = i[2]) + 4)
          }
          
          repeat{
            
            Study_summary <- UMD(Study_patientnumber, k, l, controlProp, True.sd)
            Study_mean <- Study_summary[1]
            Study_StanDev <- Study_summary[2]
            
            Begg_p <- pnorm(Study_mean/(Study_StanDev))
            
            Step_weight <- ifelse(Begg_p < Severity.boundary[1], 1, ifelse(Begg_p < Severity.boundary[2], 0.75, 0.25))
            
            if(rbinom(1,1, Step_weight) == 1 ){break}
            
          }
          
          Study.n <- as.integer(0.5*Study_patientnumber) * 2
          
          Study_estimate[o] <- Study_mean
          Study_sd[o] <- Study_StanDev
          ni[o] <- Study.n
        }
        
        ## Counter without number of studies
        counter <- as.integer((apply(sapply(Subj, function(vec) {i %in% vec}), 1, which.max)[1]-1) * length(theta) * length(tau.sq) * length(Studies) * Reps + 
                                (match(k, theta)-1) * length(tau.sq) * length(Studies) * Reps +
                                (match(l, tau.sq)-1) * length(Studies) * Reps +
                                (match(n, Studies)-1) * Reps + 
                                m
        )
        
        
        ### Temporary data.table
        temp.data <- data.table(Study_estimate, Study_sd, ni)
        
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
        
        
        Normal.Sim[dummy.counter, `:=` (Unique_ID = counter,
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
                                        Moreno_Estimate = moreno.est[1],
                                        Moreno_se = moreno.est[2],
                                        Mult_se = ma.mult$se
        )]
        
        dummy.counter <- dummy.counter + 1
        }
      }
    }
  }
  Normal.Sim
}

stopCluster(c1)

Normal.Sim.Results <- Normal.Sim.Results[order(Unique_ID)]

##### Need to re append values - specific to analysis
ID =  length(Subj) * length(controlProp) * length(theta) * length(tau.sq) * Reps * length(Studies)

Subj <- c(60, 20, 250, 4.2)

Normal.Sim.Results$Rep_Number =  rep(1:Reps, times = ID/Reps)
Normal.Sim.Results$Rep_NumStudies = rep(rep(Studies, each = Reps), times = ID/(Reps*length(Studies)))
Normal.Sim.Results$Rep_tau.sq = rep(rep(tau.sq, each = Reps * length(Studies)), times = ID/(Reps*length(Studies)*length(tau.sq)))
Normal.Sim.Results$Rep_theta = rep( rep(theta, each = Reps * length(Studies) * length(tau.sq)), times = length(Subj))
Normal.Sim.Results$Rep_Subj = rep(Subj, each = ID / length(Subj))

TimeTakenAn <- proc.time() - StartTime


#write.csv(Normal.Sim.Results, file = "NSB0V1An.csv")
saveRDS(Normal.Sim.Results, file = "NSStepRDS")


### Checking values

sum(is.na(Normal.Sim.Results))
Normal.Sim.Results[is.na(Normal.Sim.Results$REML_Est)]

#### Timings ----

TimeTakenTotal <- proc.time() - StartTimeTotal

TimeTakenAn
TimeTakenTotal