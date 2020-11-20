#### Looped Analysis UMD

### Remove previous variables
rm(list = ls())
#rm(list=setdiff(ls(), "LogOR.Complete"))

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

#### LOR Declare variables ----

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

Global_list_directory <- c("D:/LogAltOutcome1",
                           "D:/LogMethods4",
                           "D:/LogNoBias5",
                           "D:/LogStep4")

Global_list_filenames <- c("LSTotalOutRDS",
                           "LSTotalMethRDS",
                           "LSTotalV2RDS",
                           "LSTotalStepRDS")

Global_list_label <- c("AltOut",
                       "Method",
                       "None",
                       "Step")

### Declare empty data

LogOR.Complete <- data.table()

for (i in 1:6){

rm(list = c("LogOR.Sim.Results", "Total.results", "Extra.values"))

setwd(paste(Global_list_directory[i]))  

#### LOR Import data here ----

sig.level <- (1 - 0.05/2)

system.time(LogOR.Sim.Results <- readRDS(file = paste(Global_list_filenames[i]) ) )

LogOR.Sim.Results <- data.table(LogOR.Sim.Results)

Extra.values <- LogOR.Sim.Results[, .(DL_I2 = mean(DL_I2, na.rm = TRUE),
                                      Av_exc = mean(Num_exc, na.rm = TRUE)),
                                  by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq, Rep_Subj, Rep_ev_freq)]

LogOR.Sim.Results[, c("Unique_ID", "Num_exc", "DL_I2") := NULL]

##### Confidence intervals ----


## Data.table attempt
LogOR.Sim.Results[, FE_CIlb := FE_Estimate - qnorm(sig.level) * LogOR.Sim.Results$FE_se]
LogOR.Sim.Results[, FE_CIub := FE_Estimate + qnorm(sig.level) * LogOR.Sim.Results$FE_se]
LogOR.Sim.Results[, REML_CIlb := REML_Estimate - qnorm(sig.level) * LogOR.Sim.Results$REML_se]
LogOR.Sim.Results[, REML_CIub := REML_Estimate + qnorm(sig.level) * LogOR.Sim.Results$REML_se]
LogOR.Sim.Results[, DL_CIlb := DL_Estimate - qnorm(sig.level) * LogOR.Sim.Results$DL_se]
LogOR.Sim.Results[, DL_CIub := DL_Estimate + qnorm(sig.level) * LogOR.Sim.Results$DL_se]
LogOR.Sim.Results[, Doi_CIlb := FE_Estimate - qnorm(sig.level) * LogOR.Sim.Results$HC_DL_se]
LogOR.Sim.Results[, Doi_CIub := FE_Estimate + qnorm(sig.level) * LogOR.Sim.Results$HC_DL_se]
LogOR.Sim.Results[, Moreno_CIlb := Moreno_Estimate - qnorm(sig.level) * LogOR.Sim.Results$Moreno_se]
LogOR.Sim.Results[, Moreno_CIub := Moreno_Estimate + qnorm(sig.level) * LogOR.Sim.Results$Moreno_se]
LogOR.Sim.Results[, Mult_CIlb := FE_Estimate - qnorm(sig.level) * LogOR.Sim.Results$Mult_se]
LogOR.Sim.Results[, Mult_CIub := FE_Estimate + qnorm(sig.level) * LogOR.Sim.Results$Mult_se]

#### Testing large scale summaries ----

Total.results <- LogOR.Sim.Results[, .(Bias_FE = mean(FE_Estimate, na.rm = TRUE) - Rep_theta, Bias_REML = mean(REML_Estimate, na.rm = TRUE) - Rep_theta,
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
                                       EmpSE_Bias_FE = sd(FE_Estimate, na.rm = TRUE), EmpSE_Bias_DL = sd(DL_Estimate, na.rm = TRUE),
                                       EmpSE_Bias_REML = sd(REML_Estimate, na.rm = TRUE), EmpSE_Bias_Moreno = sd(Moreno_Estimate, na.rm = TRUE)
),
by = .(Rep_theta, Rep_NumStudies, Rep_tau.sq, Rep_Subj, Rep_ev_freq)]

Total.results <- merge(Total.results, Extra.values, by = c("Rep_theta", "Rep_NumStudies", "Rep_tau.sq", "Rep_Subj", "Rep_ev_freq"))


##### Add extra variable

Total.results$Bias_type <- as.factor(rep(paste(Global_list_label[i]), dim(Total.results)[1]))
#LogOR.Complete <- Total.results
LogOR.Complete <- rbindlist(list(LogOR.Complete, Total.results), use.names = TRUE, fill = FALSE)

}

#### Write out ----

setwd("D:/LOR")

write.csv(LogOR.Complete, "LogORResults_Looped.csv")
