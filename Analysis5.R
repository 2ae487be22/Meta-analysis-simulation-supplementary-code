#### Looped Analysis UMD

### Remove previous variables
rm(list = ls())
#rm(list=setdiff(ls(), "UMD.Complete"))

#### Libraries, set seed, set cores ----
library(data.table)
library(doParallel)
library(foreach)
library(doRNG)
library(copula)
library(compiler)
library(ggplot2)
library(stargazer)
library(reshape2)

#### Functions -----

CI.betw <- function(a, b, c){
  output <- ifelse(a >= b & a <= c, 1, 0)
  return(output)
}

CI.updown <- function(a, b, c){
  output <- ifelse(a >= b & a <= c, 0, ifelse(a < b, -1, 1))
  return(output)
}

#### UMD Declare variables ----

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
Tested.outcomes <- 2
Sd.split <- 0.6

# Size of per unit bias increase
Bias.multiple <- c(0, log(0.85)/(-1.81), log(0.7225)/(-1.81))

Global_list_directory <- c("D:/NormalAltOutcome1",
                           "D:/NormalMethods5",
                           "D:/NormalNoBias6",
                           "D:/NormalStep5")

Global_list_filenames <- c("NSTotalOutRDS",
                           "NSTotalMethRDS",
                           "NSTotalV2RDS",
                           "NSStepRDS")

Global_list_label <- c("AltOut",
                       "Method",
                       "None",
                       "Step")

### Declare empty data

UMD.Complete <- data.table()

for (i in 1:4){
  
  rm(list = c("Normal.Sim.Results", "Total.results"))
  
  setwd(paste(Global_list_directory[i]))  


#### UMD Import data here ----

system.time(Normal.Sim.Results <- readRDS(file = paste(Global_list_filenames[i]) ))

Normal.Sim.Results <- data.table(Normal.Sim.Results)

##### Confidence intervals ----

sig.level <- (1 - 0.05/2)


Normal.Sim.Results[, FE_CIlb := FE_Estimate - qnorm(sig.level) * FE_se]
Normal.Sim.Results[, FE_CIub := FE_Estimate + qnorm(sig.level) * FE_se]

Normal.Sim.Results[, REML_CIlb := REML_Estimate - qnorm(sig.level) * REML_se]
Normal.Sim.Results[, REML_CIub := REML_Estimate + qnorm(sig.level) * REML_se]

Normal.Sim.Results[, Doi_CIlb := FE_Estimate - qnorm(sig.level) * HC_DL_se]
Normal.Sim.Results[, Doi_CIub := FE_Estimate + qnorm(sig.level) * HC_DL_se]

### Does Moreno use z-score?

Normal.Sim.Results[, Moreno_CIlb := Moreno_Estimate - qnorm(sig.level) * Moreno_se]
Normal.Sim.Results[, Moreno_CIub := Moreno_Estimate + qnorm(sig.level) * Moreno_se]

Normal.Sim.Results[, Mult_CIlb := FE_Estimate - qnorm(sig.level) * Mult_se]
Normal.Sim.Results[, Mult_CIub := FE_Estimate + qnorm(sig.level) * Mult_se]


Normal.Sim.Results[, DL_CIlb := DL_Estimate - qnorm(sig.level) * DL_se]
Normal.Sim.Results[, DL_CIub := DL_Estimate + qnorm(sig.level) * DL_se]

#### Testing large scale summaries ----

Total.results <- Normal.Sim.Results[, .(Bias_FE = mean(FE_Estimate, na.rm = TRUE) - Rep_theta, Bias_REML = mean(REML_Estimate, na.rm = TRUE) - Rep_theta,
                                        Bias_DL = mean(DL_Estimate, na.rm = TRUE) - Rep_theta, Bias_Moreno = mean(Moreno_Estimate, na.rm = TRUE) - Rep_theta,
                                        MSE_FE = mean((FE_Estimate - Rep_theta)^2, na.rm = TRUE), MSE_REML = mean((REML_Estimate - Rep_theta)^2, na.rm = TRUE),
                                        MSE_DL = mean((DL_Estimate - Rep_theta)^2, na.rm = TRUE), "MSE Moreno" = mean((Moreno_Estimate- Rep_theta)^2, na.rm = TRUE) ,
                                        Coverage_FE = mean(CI.betw(Rep_theta, FE_CIlb, FE_CIub), na.rm = TRUE), Coverage_REML = mean(CI.betw(Rep_theta, REML_CIlb, REML_CIub), na.rm = TRUE),
                                        Coverage_DL = mean(CI.betw(Rep_theta, DL_CIlb, DL_CIub), na.rm = TRUE), Coverage_HC_DL = mean(CI.betw(Rep_theta, HC_DL_CIlb, HC_DL_CIub), na.rm = TRUE),
                                        Coverage_HC_REML = mean(CI.betw(Rep_theta, HC_REML_CIlb, HC_REML_CIub), na.rm = TRUE),
                                        Coverage_KH_DL = mean(CI.betw(Rep_theta, KH_DL_CIlb, KH_DL_CIub), na.rm = TRUE),
                                        Coverage_KH_REML = mean(CI.betw(Rep_theta, KH_REML_CIlb, KH_REML_CIub), na.rm = TRUE),
                                        Coverage_IVHet = mean(CI.betw(Rep_theta, Doi_CIlb, Doi_CIub), na.rm = TRUE),
                                        Coverage_Moreno = mean(CI.betw(Rep_theta, Moreno_CIlb, Moreno_CIub), na.rm = TRUE),
                                        Coverage_Mult = mean(CI.betw(Rep_theta, Mult_CIlb, Mult_CIub), na.rm = TRUE),
                                        Outside_FE = mean(CI.updown(Rep_theta, FE_CIlb, FE_CIub), na.rm = TRUE), Outside_REML = mean(CI.updown(Rep_theta, REML_CIlb, REML_CIub), na.rm = TRUE),
                                        Outside_DL = mean(CI.updown(Rep_theta, DL_CIlb, DL_CIub), na.rm = TRUE), Outside_HC_DL = mean(CI.updown(Rep_theta, HC_DL_CIlb, HC_DL_CIub), na.rm = TRUE),
                                        Outside_HC_REML = mean(CI.updown(Rep_theta, HC_REML_CIlb, HC_REML_CIub), na.rm = TRUE),
                                        Outside_KH_DL = mean(CI.updown(Rep_theta, KH_DL_CIlb, KH_DL_CIub), na.rm = TRUE),
                                        Outside_KH_REML = mean(CI.updown(Rep_theta, KH_REML_CIlb, KH_REML_CIub), na.rm = TRUE),
                                        Outside_IVHet = mean(CI.updown(Rep_theta, Doi_CIlb, Doi_CIub), na.rm = TRUE),
                                        Outside_Moreno = mean(CI.updown(Rep_theta, Moreno_CIlb, Moreno_CIub), na.rm = TRUE),
                                        Outside_Mult = mean(CI.updown(Rep_theta, Mult_CIlb, Mult_CIub), na.rm = TRUE),
                                        DL_I2 = mean(DL_I2, na.rm = TRUE), 
                                        EmpSE_Bias_FE = sd(FE_Estimate, na.rm = TRUE), EmpSE_Bias_DL = sd(DL_Estimate, na.rm = TRUE),
                                        EmpSE_Bias_REML = sd(REML_Estimate, na.rm = TRUE), EmpSE_Bias_Moreno = sd(Moreno_Estimate, na.rm = TRUE)
),
by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq, Rep_Subj)]

##### Add extra variable

Total.results$Bias_type <- as.factor(rep(paste(Global_list_label[i]), dim(Total.results)[1]))
#UMD.Complete <- Total.results
UMD.Complete <- rbindlist(list(UMD.Complete, Total.results), use.names = TRUE, fill = FALSE)

}

#### Write out ----

setwd("D:/UMD")

write.csv(UMD.Complete, "UMDResults_Looped.csv")
