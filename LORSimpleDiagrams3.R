#### Analysis

### Remove previous variables
rm(list = ls())

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


Global_list_formalname <- c("Adjusted Outcome Bias",
                            "Methodological Bias",
                            "No Bias",
                            "Step Bias")

###  Loop ----

for (xi in 1:4){
  
  rm(list=setdiff(ls(), c("Global_list_directory", "Global_list_filenames", "Global_list_formalname", "xi")))
  
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
  
  sig.level <- (1 - 0.05/2)
  
  FileLoc <- Global_list_directory[xi]
  setwd(paste(FileLoc))
  
  
  #### LOR Import data here ----
  
  system.time(LogOR.Sim.Results <- readRDS(file = paste(Global_list_filenames[xi]) ))
  
  LogOR.Sim.Results <- data.table(LogOR.Sim.Results)
  
  #### Renaming variables ----
  
  # # Theta in log conditions
  # theta = c("log(0.25)", "log(0.8)", "log(1)", "log(1.25)", "log(4)")
  # LogOR.Sim.Results[, Rep_theta := factor(Rep_theta, labels = c("log(0.25)", "log(0.8)", "log(1)", "log(1.25)", "log(4)"))]
  # 
  
  LogOR.Sim.Results[, Rep_theta_word := factor(Rep_theta, labels = c("log(0.25)", "log(0.8)", "log(1)", "log(1.25)", "log(4)"))]
  
  
  # Subject in all conditions
  Subj <- c("Empirical", "Small", "Fixed", "Large")
  LogOR.Sim.Results[, Rep_Subj := factor(Rep_Subj, labels = c("Empirical", "Small", "Fixed", "Large"))]
  
  
  #### Specify type of bias and directory ----
  getwd()
  Type.of.bias <- Global_list_formalname[xi]
  mainDir <- file.path(getwd(), Type.of.bias)
  
  dir.create(file.path(mainDir), recursive = TRUE, showWarnings = FALSE)
  
  
  
  #### New code -----------------------------------
  
  ## Set held conditions
  
  c <- Subj[1]
  a <- theta[3]
  e <- Studies[2]
  b <- tau.sq[3]
  d <- EvFreq[2]
  
  subDir <- file.path("AllHeld")
  
  dir.create(file.path(mainDir, subDir), recursive = TRUE, showWarnings = FALSE)
  
  setwd(file.path(mainDir, subDir))
  
  #### Vary by Subj ----
  
  
  #------------------- Bias
  
  An.Cond <- LogOR.Sim.Results[Rep_NumStudies == e & Rep_tau.sq == b & Rep_theta == a & Rep_ev_freq == d]
  
  Bias.values.total <- An.Cond[, .(FE = FE_Estimate - Rep_theta,
                                   DL = DL_Estimate - Rep_theta, "Egger-Var" = Moreno_Estimate - Rep_theta),
                               by = .(Rep_Subj)]
  
  Bias.values.total2 <- melt(Bias.values.total, id = c("Rep_Subj"))
  
  p1 <- ggplot(Bias.values.total2, aes(x = as.factor(Rep_Subj), y = value, fill = variable ) ) +
    #geom_boxplot(alpha=1, outlier.shape = NA, position=position_dodge(1)) +
    geom_boxplot(alpha=1, outlier.shape = NA, coef = 0) +
    #ggtitle("Bias") + 
    xlab("Size of Studies") + ylab("Bias") + 
    theme_bw() + labs(fill = "Method") +
    coord_cartesian(ylim = c(-0.75,0.75)) #+
  #facet_grid(Rep_theta~ Rep_tau.sq) + 
  #ggtitle(paste("Bias in", c, "studies:", Type.of.bias, sep = " "))
  
  ggsave(paste("LOR", Type.of.bias, "BiasBoxPlot", "SubjVarying", ".pdf", sep = ""), p1, dpi = 800, device = "pdf", width = 8.01, 
         height = 5.67, units = "in")
  
  
  Av.Bias.values.total <- An.Cond[, .(FE = mean(FE_Estimate - Rep_theta),
                                      DL = mean(DL_Estimate - Rep_theta), "Egger-Var" = mean(Moreno_Estimate - Rep_theta)),
                                  by = .(Rep_Subj)]
  
  Av.Bias.values.total2 <- melt(Av.Bias.values.total, id = c("Rep_Subj"))
  
  p2 <- ggplot(Av.Bias.values.total2, aes(x = as.factor(Rep_Subj), y = value) ) +
    geom_point(aes(shape = variable, colour = variable), position=position_dodge(0.3), size = 3) +
    #geom_boxplot(alpha=1, outlier.shape = NA, position=position_dodge(1)) +
    #geom_boxplot(alpha=1, outlier.shape = NA, coef = 0) +
    #ggtitle("Bias") + 
    xlab("Size of Studies") + ylab("Bias") + 
    theme_bw() +
    coord_cartesian(ylim = c(-0.200,0.025)) + labs(colour = "Method", shape = "Method") #+
  #facet_grid(Rep_theta~ Rep_tau.sq) + 
  #ggtitle(paste("Bias in", c, "studies:", Type.of.bias, sep = " "))
  
  #p16v2
  
  ggsave(paste("LOR", Type.of.bias, "BiasPointPlot", "SubjVarying", ".pdf", sep = ""), p2, dpi = 800, device = "pdf", width = 8.01, 
         height = 5.67, units = "in")
  
  
  #------------------- MSE
  
  MSE1.values <- An.Cond[, .(FE = mean((FE_Estimate - Rep_theta)^2, na.rm = TRUE),
                             DL = mean((DL_Estimate - Rep_theta)^2, na.rm = TRUE), "Egger-Var" = mean((Moreno_Estimate- Rep_theta)^2, na.rm = TRUE) ),
                         by = .(Rep_Subj)]
  
  MSE1.values2 <- melt(MSE1.values, id = c("Rep_Subj"))
  
  p3 <- ggplot(MSE1.values2, aes(x = as.factor(Rep_Subj), y = value, colour = variable)) +
    geom_point(aes(shape = variable, colour = variable), position=position_dodge(0.3), size = 3) + theme_bw() +
    xlab("Size of Studies") + ylab("MSE (log scale)")  + scale_y_log10() + coord_cartesian(ylim = c(0.001, 1)) + # scale_colour_grey()+
    labs(colour = "Method", shape = "Method")
  #facet_grid(Rep_Subj~ Rep_tau.sq) + ggtitle(paste("MSE with", e, "studies:", Type.of.bias, sep = " "))#+ coord_cartesian(ylim = c(0, 1))
  
  #p15
  
  ggsave(paste("LOR", Type.of.bias, "MSEPointPlot", "SubjVarying", ".pdf", sep = ""), p3, dpi = 800, device = "pdf", width = 8.01, 
         height = 5.67, units = "in")
  
  
  #------------------- Coverage
  
  An.Cond[, FE_CIlb := FE_Estimate - qnorm(sig.level) * FE_se]
  An.Cond[, FE_CIub := FE_Estimate + qnorm(sig.level) * FE_se]
  
  An.Cond[, Doi_CIlb := FE_Estimate - qnorm(sig.level) * HC_DL_se]
  An.Cond[, Doi_CIub := FE_Estimate + qnorm(sig.level) * HC_DL_se]
  
  An.Cond[, Mult_CIlb := FE_Estimate - qnorm(sig.level) * Mult_se]
  An.Cond[, Mult_CIub := FE_Estimate + qnorm(sig.level) * Mult_se]
  
  
  An.Cond[, DL_CIlb := DL_Estimate - qnorm(sig.level) * DL_se]
  An.Cond[, DL_CIub := DL_Estimate + qnorm(sig.level) * DL_se]
  
  Coverage.values <- An.Cond[, .(FE = mean(CI.betw(Rep_theta, FE_CIlb, FE_CIub), na.rm = TRUE),
                                 DL = mean(CI.betw(Rep_theta, DL_CIlb, DL_CIub), na.rm = TRUE), "HC DL" = mean(CI.betw(Rep_theta, HC_DL_CIlb, HC_DL_CIub), na.rm = TRUE),
                                 "KH DL" = mean(CI.betw(Rep_theta, KH_DL_CIlb, KH_DL_CIub), na.rm = TRUE),
                                 IVHet = mean(CI.betw(Rep_theta, Doi_CIlb, Doi_CIub), na.rm = TRUE),
                                 Mult = mean(CI.betw(Rep_theta, Mult_CIlb, Mult_CIub), na.rm = TRUE)
  ),
  by = .(Rep_Subj)]
  
  Coverage.values2<- melt(Coverage.values, id = c("Rep_Subj"))
  
  p4 <-  ggplot(Coverage.values2, aes(x = as.factor(Rep_Subj), y = value, colour = variable, shape = variable)) +
    geom_point(alpha = 1, position=position_dodge(width=0.2), size = 3)  + labs(colour = "Method", shape = "Method") +
    xlab("Size of Studies")  + ylab("Coverage")  + geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") + 
    coord_cartesian(ylim = c(0.5, 1)) + scale_shape_manual(values=c(16, 17, 15, 4, 7, 8)) + theme_bw() #+  scale_colour_grey()
  
  #p19
  
  ggsave(paste("LOR", Type.of.bias, "CovPointPlot", "SubjVarying", ".pdf", sep = ""), p4, dpi = 800, device = "pdf", width = 8.01, height = 5.67, units = "in")
  
  
  
  #### Vary by number of studies ----
  
  An.Cond <- LogOR.Sim.Results[Rep_Subj == c & Rep_tau.sq == b & Rep_theta == a & Rep_ev_freq == d]
  
  Bias.values.total <- An.Cond[, .(FE = FE_Estimate - Rep_theta,
                                   DL = DL_Estimate - Rep_theta, "Egger-Var" = Moreno_Estimate - Rep_theta),
                               by = .(Rep_NumStudies)]
  
  Bias.values.total2 <- melt(Bias.values.total, id = c("Rep_NumStudies"))
  
  p5 <- ggplot(Bias.values.total2, aes(x = as.factor(Rep_NumStudies), y = value, fill = variable ) ) +
    geom_boxplot(alpha=1, outlier.shape = NA, coef = 0) +
    xlab("Number of Studies") + ylab("Bias") + 
    theme_bw() +
    coord_cartesian(ylim = c(-0.75,0.75)) + labs(fill = "Method")
  
  ggsave(paste("LOR", Type.of.bias, "BiasBoxPlot", "NumStudiesVarying", ".pdf", sep = ""), p5, dpi = 800, device = "pdf", width = 8.01, 
         height = 5.67, units = "in")
  
  
  Av.Bias.values.total <- An.Cond[, .(FE = mean(FE_Estimate - Rep_theta),
                                      DL = mean(DL_Estimate - Rep_theta), "Egger-Var" = mean(Moreno_Estimate - Rep_theta)),
                                  by = .(Rep_NumStudies)]
  
  Av.Bias.values.total2 <- melt(Av.Bias.values.total, id = c("Rep_NumStudies"))
  
  p6 <- ggplot(Av.Bias.values.total2, aes(x = Rep_NumStudies, y = value) ) +
    geom_point(aes(shape = variable, colour = variable), position=position_dodge(0.3), size = 3) +
    #scale_linetype_manual(values=c("solid", "twodash", "dotted")) +
    geom_line(aes(colour = variable)) + 
    xlab("Number of Studies") + ylab("Bias") + 
    theme_bw() +
    coord_cartesian(ylim = c(-0.200,0.025)) + labs(colour = "Method", shape = "Method")
  
  #p16v2
  
  ggsave(paste("LOR", Type.of.bias, "BiasPointPlot", "NumStudiesVarying", ".pdf", sep = ""), p6, dpi = 800, device = "pdf", width = 8.01, 
         height = 5.67, units = "in")
  
  
  #------------------- MSE
  
  MSE1.values <- An.Cond[, .(FE = mean((FE_Estimate - Rep_theta)^2, na.rm = TRUE),
                             DL = mean((DL_Estimate - Rep_theta)^2, na.rm = TRUE), "Egger-Var" = mean((Moreno_Estimate- Rep_theta)^2, na.rm = TRUE) ),
                         by = .(Rep_NumStudies)]
  
  MSE1.values2 <- melt(MSE1.values, id = c("Rep_NumStudies"))
  
  p7 <- ggplot(MSE1.values2, aes(x = Rep_NumStudies, y = value, colour = variable)) +
    geom_point(aes(shape = variable, colour = variable), position=position_dodge(0.3), size = 3) +
    geom_line(aes(colour = variable)) + 
    theme_bw() +
    xlab("Number of Studies") + ylab("MSE (log scale)")  + scale_y_log10() + coord_cartesian(ylim = c(0.001, 1)) + # scale_colour_grey()+
    labs(colour = "Method", shape = "Method")
  #facet_grid(Rep_Subj~ Rep_tau.sq) + ggtitle(paste("MSE with", e, "studies:", Type.of.bias, sep = " "))#+ coord_cartesian(ylim = c(0, 1))
  
  #p15
  
  ggsave(paste("LOR", Type.of.bias, "MSEPointPlot", "NumStudiesVarying", ".pdf", sep = ""), p7, dpi = 800, device = "pdf", width = 8.01, 
         height = 5.67, units = "in")
  
  #------------------- Coverage
  
  An.Cond[, FE_CIlb := FE_Estimate - qnorm(sig.level) * FE_se]
  An.Cond[, FE_CIub := FE_Estimate + qnorm(sig.level) * FE_se]
  
  An.Cond[, Doi_CIlb := FE_Estimate - qnorm(sig.level) * HC_DL_se]
  An.Cond[, Doi_CIub := FE_Estimate + qnorm(sig.level) * HC_DL_se]
  
  An.Cond[, Mult_CIlb := FE_Estimate - qnorm(sig.level) * Mult_se]
  An.Cond[, Mult_CIub := FE_Estimate + qnorm(sig.level) * Mult_se]
  
  
  An.Cond[, DL_CIlb := DL_Estimate - qnorm(sig.level) * DL_se]
  An.Cond[, DL_CIub := DL_Estimate + qnorm(sig.level) * DL_se]
  
  Coverage.values <- An.Cond[, .(FE = mean(CI.betw(Rep_theta, FE_CIlb, FE_CIub), na.rm = TRUE),
                                 DL = mean(CI.betw(Rep_theta, DL_CIlb, DL_CIub), na.rm = TRUE), "HC DL" = mean(CI.betw(Rep_theta, HC_DL_CIlb, HC_DL_CIub), na.rm = TRUE),
                                 "KH DL" = mean(CI.betw(Rep_theta, KH_DL_CIlb, KH_DL_CIub), na.rm = TRUE),
                                 IVHet = mean(CI.betw(Rep_theta, Doi_CIlb, Doi_CIub), na.rm = TRUE),
                                 Mult = mean(CI.betw(Rep_theta, Mult_CIlb, Mult_CIub), na.rm = TRUE)
  ),
  by = .(Rep_NumStudies)]
  
  Coverage.values2<- melt(Coverage.values, id = c("Rep_NumStudies"))
  
  p8 <-  ggplot(Coverage.values2, aes(x = Rep_NumStudies, y = value, colour = variable, shape = variable)) +
    geom_point(alpha = 1, position=position_dodge(width=0.2), size = 3)  +
    geom_line(aes(colour = variable)) + 
    labs(colour = "Method", shape = "Method") +
    xlab("Number of Studies")  + ylab("Coverage")  + geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") + 
    coord_cartesian(ylim = c(0.5, 1)) + scale_shape_manual(values=c(16, 17, 15, 4, 7, 8)) + theme_bw() #+  scale_colour_grey()
  
  #p19
  
  ggsave(paste("LOR", Type.of.bias, "CovPointPlot", "NumStudiesVarying", ".pdf", sep = ""), p8, dpi = 800, device = "pdf", width = 8.01, height = 5.67, units = "in")
  
  
  
  #### Vary by Theta ----
  
  An.Cond <- LogOR.Sim.Results[Rep_NumStudies == e & Rep_Subj == c & Rep_tau.sq == b & Rep_ev_freq == d]
  
  Bias.values.total <- An.Cond[, .(FE = FE_Estimate - Rep_theta,
                                   DL = DL_Estimate - Rep_theta, "Egger-Var" = Moreno_Estimate - Rep_theta),
                               by = .(Rep_theta_word)]
  
  Bias.values.total2 <- melt(Bias.values.total, id = c("Rep_theta_word"))
  
  p9 <- ggplot(Bias.values.total2, aes(x = as.factor(Rep_theta_word), y = value, fill = variable ) ) +
    geom_boxplot(alpha=1, outlier.shape = NA, coef = 0) +
    xlab("Theta") + ylab("Bias") + 
    theme_bw() +
    coord_cartesian(ylim = c(-0.75,0.75)) + labs(fill = "Method")
  
  #p16
  
  ggsave(paste("LOR", Type.of.bias, "BiasBoxPlot", "ThetaVarying", ".pdf", sep = ""), p9, dpi = 800, device = "pdf", width = 8.01, 
         height = 5.67, units = "in")
  
  Av.Bias.values.total <- An.Cond[, .(FE = mean(FE_Estimate - Rep_theta),
                                      DL = mean(DL_Estimate - Rep_theta), "Egger-Var" = mean(Moreno_Estimate - Rep_theta)),
                                  by = .(Rep_theta, Rep_theta_word)]
  
  Av.Bias.values.total2 <- melt(Av.Bias.values.total, id = c("Rep_theta", "Rep_theta_word"))
  
  p10 <- ggplot(Av.Bias.values.total2, aes(x = Rep_theta, y = value) ) +
    geom_point(aes(shape = variable, colour = variable), position=position_dodge(0.05), size = 3) +
    xlab("Theta") + ylab("Bias") + 
    geom_line(aes(colour = variable)) + 
    theme_bw() + #scale_x_continuous(breaks = c(log(0.25), log(0.8), log(1), log(1.25), log(4)), labels = c("ln(0.25)", "ln(0.8)", "ln(1)", "ln(1.25)", "ln(4)")) +
    coord_cartesian(ylim = c(-0.200,0.025)) + labs(colour = "Method", shape = "Method") #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  #p16v2
  
  ggsave(paste("LOR", Type.of.bias, "BiasPointPlot", "ThetaVarying", ".pdf", sep = ""), p10, dpi = 800, device = "pdf", width = 8.01, 
         height = 5.67, units = "in")
  
  #------------------- MSE
  
  MSE1.values <- An.Cond[, .(FE = mean((FE_Estimate - Rep_theta)^2, na.rm = TRUE),
                             DL = mean((DL_Estimate - Rep_theta)^2, na.rm = TRUE), "Egger-Var" = mean((Moreno_Estimate- Rep_theta)^2, na.rm = TRUE) ),
                         by = .(Rep_theta)]
  
  MSE1.values2 <- melt(MSE1.values, id = c("Rep_theta"))
  
  p11 <- ggplot(MSE1.values2, aes(x = Rep_theta, y = value, colour = variable)) +
    geom_point(aes(shape = variable, colour = variable), position=position_dodge(0.05), size = 3) +  
    geom_line(aes(colour = variable)) + 
    theme_bw() + #scale_x_continuous(breaks = c(log(0.25), log(0.8), log(1), log(1.25), log(4)), labels = c("ln(0.25)", "ln(0.8)", "ln(1)", "ln(1.25)", "ln(4)")) +
    xlab("Theta") + ylab("MSE (log scale)")  + scale_y_log10() + coord_cartesian(ylim = c(0.001, 1)) + # scale_colour_grey()+
    labs(colour = "Method", shape = "Method")
  #facet_grid(Rep_Subj~ Rep_tau.sq) + ggtitle(paste("MSE with", e, "studies:", Type.of.bias, sep = " "))#+ coord_cartesian(ylim = c(0, 1))
  
  #p15
  
  ggsave(paste("LOR", Type.of.bias, "MSEPointPlot", "ThetaVarying", ".pdf", sep = ""), p11, dpi = 800, device = "pdf", width = 8.01, 
         height = 5.67, units = "in")
  
  #------------------- Coverage
  
  An.Cond[, FE_CIlb := FE_Estimate - qnorm(sig.level) * FE_se]
  An.Cond[, FE_CIub := FE_Estimate + qnorm(sig.level) * FE_se]
  
  An.Cond[, Doi_CIlb := FE_Estimate - qnorm(sig.level) * HC_DL_se]
  An.Cond[, Doi_CIub := FE_Estimate + qnorm(sig.level) * HC_DL_se]
  
  An.Cond[, Mult_CIlb := FE_Estimate - qnorm(sig.level) * Mult_se]
  An.Cond[, Mult_CIub := FE_Estimate + qnorm(sig.level) * Mult_se]
  
  
  An.Cond[, DL_CIlb := DL_Estimate - qnorm(sig.level) * DL_se]
  An.Cond[, DL_CIub := DL_Estimate + qnorm(sig.level) * DL_se]
  
  Coverage.values <- An.Cond[, .(FE = mean(CI.betw(Rep_theta, FE_CIlb, FE_CIub), na.rm = TRUE),
                                 DL = mean(CI.betw(Rep_theta, DL_CIlb, DL_CIub), na.rm = TRUE), "HC DL" = mean(CI.betw(Rep_theta, HC_DL_CIlb, HC_DL_CIub), na.rm = TRUE),
                                 "KH DL" = mean(CI.betw(Rep_theta, KH_DL_CIlb, KH_DL_CIub), na.rm = TRUE),
                                 IVHet = mean(CI.betw(Rep_theta, Doi_CIlb, Doi_CIub), na.rm = TRUE),
                                 Mult = mean(CI.betw(Rep_theta, Mult_CIlb, Mult_CIub), na.rm = TRUE)
  ),
  by = .(Rep_theta)]
  
  Coverage.values2<- melt(Coverage.values, id = c("Rep_theta"))
  
  p12 <-  ggplot(Coverage.values2, aes(x = Rep_theta, y = value, colour = variable, shape = variable)) +
    geom_point(alpha = 1, position=position_dodge(width=0.05), size = 3)  +
    geom_line(aes(colour = variable)) + 
    labs(colour = "Method", shape = "Method") +
    xlab("Theta")  + ylab("Coverage")  + geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") + 
    coord_cartesian(ylim = c(0.5, 1)) + scale_shape_manual(values=c(16, 17, 15, 4, 7, 8)) + theme_bw() #+ scale_x_continuous(breaks = c(log(0.25), log(0.8), log(1), log(1.25), log(4)), labels = c("ln(0.25)", "ln(0.8)", "ln(1)", "ln(1.25)", "ln(4)")) #+  scale_colour_grey()
  
  #p19
  
  ggsave(paste("LOR", Type.of.bias, "CovPointPlot", "ThetaVarying", ".pdf", sep = ""), p12, dpi = 800, device = "pdf", width = 8.01, height = 5.67, units = "in")
  
  
  
  
  #### Vary by tau.sq ----
  
  An.Cond <- LogOR.Sim.Results[Rep_NumStudies == e & Rep_Subj == c & Rep_theta == a & Rep_ev_freq == d]
  
  Bias.values.total <- An.Cond[, .(FE = FE_Estimate - Rep_theta,
                                   DL = DL_Estimate - Rep_theta, "Egger-Var" = Moreno_Estimate - Rep_theta),
                               by = .(Rep_tau.sq)]
  
  Bias.values.total2 <- melt(Bias.values.total, id = c("Rep_tau.sq"))
  
  p13 <- ggplot(Bias.values.total2, aes(x = as.factor(Rep_tau.sq), y = value, fill = variable ) ) +
    geom_boxplot(alpha=1, outlier.shape = NA, coef = 0) +
    xlab("Tau squared") + ylab("Bias") + 
    theme_bw() +
    coord_cartesian(ylim = c(-0.75,0.75)) + labs(fill = "Method")
  
  #p16
  
  ggsave(paste("LOR", Type.of.bias, "BiasBoxPlot", "Tau2Varying", ".pdf", sep = ""), p13, dpi = 800, device = "pdf", width = 8.01, 
         height = 5.67, units = "in")
  
  Av.Bias.values.total <- An.Cond[, .(FE = mean(FE_Estimate - Rep_theta),
                                      DL = mean(DL_Estimate - Rep_theta), "Egger-Var" = mean(Moreno_Estimate - Rep_theta)),
                                  by = .(Rep_tau.sq)]
  
  Av.Bias.values.total2 <- melt(Av.Bias.values.total, id = c("Rep_tau.sq"))
  
  p14 <- ggplot(Av.Bias.values.total2, aes(x = as.factor(Rep_tau.sq), y = value) ) +
    geom_point(aes(shape = variable, colour = variable), position=position_dodge(0.3), size = 3) +
    xlab("Tau squared") + ylab("Bias") + 
    theme_bw() +
    coord_cartesian(ylim = c(-0.200,0.025)) + labs(colour = "Method", shape = "Method")
  
  #p16v2
  
  ggsave(paste("LOR", Type.of.bias, "BiasPointPlot", "Tau2Varying", ".pdf", sep = ""), p14, dpi = 800, device = "pdf", width = 8.01, 
         height = 5.67, units = "in")
  
  #------------------- MSE
  
  MSE1.values <- An.Cond[, .(FE = mean((FE_Estimate - Rep_theta)^2, na.rm = TRUE),
                             DL = mean((DL_Estimate - Rep_theta)^2, na.rm = TRUE), "Egger-Var" = mean((Moreno_Estimate- Rep_theta)^2, na.rm = TRUE) ),
                         by = .(Rep_tau.sq)]
  
  MSE1.values2 <- melt(MSE1.values, id = c("Rep_tau.sq"))
  
  p15 <- ggplot(MSE1.values2, aes(x = as.factor(Rep_tau.sq), y = value, colour = variable)) +
    geom_point(aes(shape = variable, colour = variable), position=position_dodge(0.1), size = 3) + theme_bw() +
    xlab("Tau squared") + ylab("MSE (log scale)")  + scale_y_log10() + coord_cartesian(ylim = c(0.001, 1)) + # scale_colour_grey()+
    labs(colour = "Method", shape = "Method")
  #facet_grid(Rep_Subj~ Rep_tau.sq) + ggtitle(paste("MSE with", e, "studies:", Type.of.bias, sep = " "))#+ coord_cartesian(ylim = c(0, 1))
  
  #p15
  
  ggsave(paste("LOR", Type.of.bias, "MSEPointPlot", "Tau2Varying", ".pdf", sep = ""), p15, dpi = 800, device = "pdf", width = 8.01, 
         height = 5.67, units = "in")
  
  #------------------- Coverage
  
  An.Cond[, FE_CIlb := FE_Estimate - qnorm(sig.level) * FE_se]
  An.Cond[, FE_CIub := FE_Estimate + qnorm(sig.level) * FE_se]
  
  An.Cond[, Doi_CIlb := FE_Estimate - qnorm(sig.level) * HC_DL_se]
  An.Cond[, Doi_CIub := FE_Estimate + qnorm(sig.level) * HC_DL_se]
  
  An.Cond[, Mult_CIlb := FE_Estimate - qnorm(sig.level) * Mult_se]
  An.Cond[, Mult_CIub := FE_Estimate + qnorm(sig.level) * Mult_se]
  
  
  An.Cond[, DL_CIlb := DL_Estimate - qnorm(sig.level) * DL_se]
  An.Cond[, DL_CIub := DL_Estimate + qnorm(sig.level) * DL_se]
  
  Coverage.values <- An.Cond[, .(FE = mean(CI.betw(Rep_theta, FE_CIlb, FE_CIub), na.rm = TRUE),
                                 DL = mean(CI.betw(Rep_theta, DL_CIlb, DL_CIub), na.rm = TRUE), "HC DL" = mean(CI.betw(Rep_theta, HC_DL_CIlb, HC_DL_CIub), na.rm = TRUE),
                                 "KH DL" = mean(CI.betw(Rep_theta, KH_DL_CIlb, KH_DL_CIub), na.rm = TRUE),
                                 IVHet = mean(CI.betw(Rep_theta, Doi_CIlb, Doi_CIub), na.rm = TRUE),
                                 Mult = mean(CI.betw(Rep_theta, Mult_CIlb, Mult_CIub), na.rm = TRUE)
  ),
  by = .(Rep_tau.sq)]
  
  Coverage.values2<- melt(Coverage.values, id = c("Rep_tau.sq"))
  
  p16 <-  ggplot(Coverage.values2, aes(x = as.factor(Rep_tau.sq), y = value, colour = variable, shape = variable)) +
    geom_point(alpha = 1, position=position_dodge(width=0.2), size = 3)  + labs(colour = "Method", shape = "Method") +
    xlab("Tau squared")  + ylab("Coverage")  + geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") + 
    coord_cartesian(ylim = c(0.5, 1)) + scale_shape_manual(values=c(16, 17, 15, 4, 7, 8)) + theme_bw() #+  scale_colour_grey()
  
  #p19
  
  ggsave(paste("LOR", Type.of.bias, "CovPointPlot", "Tau2Varying", ".pdf", sep = ""), p16, dpi = 800, device = "pdf", width = 8.01, height = 5.67, units = "in")
  
  
  #### Vary by event rate ----
  
  An.Cond <- LogOR.Sim.Results[Rep_NumStudies == e & Rep_Subj == c & Rep_theta == a & Rep_tau.sq == b]
  
  Bias.values.total <- An.Cond[, .(FE = FE_Estimate - Rep_theta,
                                   DL = DL_Estimate - Rep_theta, "Egger-Var" = Moreno_Estimate - Rep_theta),
                               by = .(Rep_ev_freq)]
  
  Bias.values.total2 <- melt(Bias.values.total, id = c("Rep_ev_freq"))
  
  p17 <- ggplot(Bias.values.total2, aes(x = as.factor(Rep_ev_freq), y = value, fill = variable ) ) +
    geom_boxplot(alpha=1, outlier.shape = NA, coef = 0) +
    xlab("Event rate") + ylab("Bias") + 
    theme_bw() +
    coord_cartesian(ylim = c(-0.75,0.75)) + labs(fill = "Method")
  
  #p16
  
  ggsave(paste("LOR", Type.of.bias, "BiasBoxPlot", "EventRateVarying", ".pdf", sep = ""), p17, dpi = 800, device = "pdf", width = 8.01, 
         height = 5.67, units = "in")
  
  Av.Bias.values.total <- An.Cond[, .(FE = mean(FE_Estimate - Rep_theta),
                                      DL = mean(DL_Estimate - Rep_theta), "Egger-Var" = mean(Moreno_Estimate - Rep_theta)),
                                  by = .(Rep_ev_freq)]
  
  Av.Bias.values.total2 <- melt(Av.Bias.values.total, id = c("Rep_ev_freq"))
  
  p18 <- ggplot(Av.Bias.values.total2, aes(x = as.factor(Rep_ev_freq), y = value, group = variable) ) +
    geom_point(aes(shape = variable, colour = variable), position=position_dodge(0.3), size = 3) +
    xlab("Event rate") + ylab("Bias") + 
    geom_line(aes(colour = variable), position=position_dodge(0.3)) + 
    theme_bw() +
    coord_cartesian(ylim = c(-0.200,0.025)) + labs(colour = "Method", shape = "Method")
  
  # p10 <- ggplot(Av.Bias.values.total2, aes(x = Rep_theta, y = value) ) +
  #   geom_point(aes(shape = variable, colour = variable), position=position_dodge(0.05), size = 3) +
  #   xlab("Theta") + ylab("Bias") + 
  #   geom_line(aes(colour = variable)) + 
  #   theme_bw() + #scale_x_continuous(breaks = c(log(0.25), log(0.8), log(1), log(1.25), log(4)), labels = c("ln(0.25)", "ln(0.8)", "ln(1)", "ln(1.25)", "ln(4)")) +
  #   coord_cartesian(ylim = c(-0.200,0.025)) + labs(colour = "Method", shape = "Method") #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
  # 
  
  #p16v2
  
  ggsave(paste("LOR", Type.of.bias, "BiasPointPlot", "EventRateVarying", ".pdf", sep = ""), p18, dpi = 800, device = "pdf", width = 8.01, 
         height = 5.67, units = "in")
  
  #------------------- MSE
  
  MSE1.values <- An.Cond[, .(FE = mean((FE_Estimate - Rep_theta)^2, na.rm = TRUE),
                             DL = mean((DL_Estimate - Rep_theta)^2, na.rm = TRUE), "Egger-Var" = mean((Moreno_Estimate- Rep_theta)^2, na.rm = TRUE) ),
                         by = .(Rep_ev_freq)]
  
  MSE1.values2 <- melt(MSE1.values, id = c("Rep_ev_freq"))
  
  p19 <- ggplot(MSE1.values2, aes(x = as.factor(Rep_ev_freq), y = value, colour = variable, group = variable)) +
    geom_point(aes(shape = variable, colour = variable), position=position_dodge(0.1), size = 3) + theme_bw() +
    xlab("Event rate") + ylab("MSE (log scale)")  + scale_y_log10() + coord_cartesian(ylim = c(0.001, 1)) + # scale_colour_grey()+
    geom_line(aes(colour = variable), position=position_dodge(0.1)) +
    labs(colour = "Method", shape = "Method")
  #facet_grid(Rep_Subj~ Rep_tau.sq) + ggtitle(paste("MSE with", e, "studies:", Type.of.bias, sep = " "))#+ coord_cartesian(ylim = c(0, 1))
  
  #p15
  
  ggsave(paste("LOR", Type.of.bias, "MSEPointPlot", "EventRateVarying", ".pdf", sep = ""), p19, dpi = 800, device = "pdf", width = 8.01, 
         height = 5.67, units = "in")
  
  #------------------- Coverage
  
  An.Cond[, FE_CIlb := FE_Estimate - qnorm(sig.level) * FE_se]
  An.Cond[, FE_CIub := FE_Estimate + qnorm(sig.level) * FE_se]
  
  An.Cond[, Doi_CIlb := FE_Estimate - qnorm(sig.level) * HC_DL_se]
  An.Cond[, Doi_CIub := FE_Estimate + qnorm(sig.level) * HC_DL_se]
  
  An.Cond[, Mult_CIlb := FE_Estimate - qnorm(sig.level) * Mult_se]
  An.Cond[, Mult_CIub := FE_Estimate + qnorm(sig.level) * Mult_se]
  
  
  An.Cond[, DL_CIlb := DL_Estimate - qnorm(sig.level) * DL_se]
  An.Cond[, DL_CIub := DL_Estimate + qnorm(sig.level) * DL_se]
  
  Coverage.values <- An.Cond[, .(FE = mean(CI.betw(Rep_theta, FE_CIlb, FE_CIub), na.rm = TRUE),
                                 DL = mean(CI.betw(Rep_theta, DL_CIlb, DL_CIub), na.rm = TRUE), "HC DL" = mean(CI.betw(Rep_theta, HC_DL_CIlb, HC_DL_CIub), na.rm = TRUE),
                                 "KH DL" = mean(CI.betw(Rep_theta, KH_DL_CIlb, KH_DL_CIub), na.rm = TRUE),
                                 IVHet = mean(CI.betw(Rep_theta, Doi_CIlb, Doi_CIub), na.rm = TRUE),
                                 Mult = mean(CI.betw(Rep_theta, Mult_CIlb, Mult_CIub), na.rm = TRUE)
  ),
  by = .(Rep_ev_freq)]
  
  Coverage.values2<- melt(Coverage.values, id = c("Rep_ev_freq"))
  
  p20 <-  ggplot(Coverage.values2, aes(x = as.factor(Rep_ev_freq), y = value, colour = variable, shape = variable, group = variable)) +
    geom_point(alpha = 1, position=position_dodge(width=0.2), size = 3)  + labs(colour = "Method", shape = "Method") +
    xlab("Event rate")  + ylab("Coverage")  + geom_hline(aes(yintercept = 0.95), alpha = 0.3, linetype = "longdash") + 
    geom_line(aes(colour = variable), position=position_dodge(width=0.2)) +
    coord_cartesian(ylim = c(0.5, 1)) + scale_shape_manual(values=c(16, 17, 15, 4, 7, 8)) + theme_bw() #+  scale_colour_grey()
  
  #p19
  
  ggsave(paste("LOR", Type.of.bias, "CovPointPlot", "EventRateVarying", ".pdf", sep = ""), p20, dpi = 800, device = "pdf", width = 8.01, height = 5.67, units = "in")
  
  
  #### Combination plot ----
  
  library(gridExtra)
  library(gtable)
  library(grid)
  
  
  legend2 <- gtable_filter(ggplotGrob(p2), "guide-box")
  
  # Bias
  
  a <- grid.arrange(p2 + ylab(NULL) + theme(legend.position="none"),
                    p6 + ylab(NULL) + theme(legend.position="none"),
                    p10 + ylab(NULL) + theme(legend.position="none"),
                    p14 + ylab(NULL) + theme(legend.position="none"),
                    p18 + ylab(NULL) + theme(legend.position="none"),
                    legend2,
                    nrow = 3, ncol = 2,
                    left = textGrob("Bias", rot = 90, vjust = 1)
  )
  
  ggsave(paste("LOR", Type.of.bias, "BiasPlot", ".pdf", sep = ""), a, dpi = 800, device = "pdf", width = 1.5*8.01, height = 2.25*5.67, units = "in")
  ggsave(paste("LOR", Type.of.bias, "BiasPlot", ".eps", sep = ""), a, dpi = 800, device = "eps", width = 1.5*8.01, height = 2.25*5.67, units = "in")
  
  
  # MSE
  
  legend2 <- gtable_filter(ggplotGrob(p3), "guide-box")
  
  
  a <- grid.arrange(p3 + ylab(NULL) + theme(legend.position="none"),
                    p7 + ylab(NULL) + theme(legend.position="none"),
                    p11 + ylab(NULL) + theme(legend.position="none"),
                    p15 + ylab(NULL) + theme(legend.position="none"),
                    p19 + ylab(NULL) + theme(legend.position="none"),
                    legend2,
                    nrow = 3, ncol = 2,
                    left = textGrob("MSE (log scale)", rot = 90, vjust = 1)
  )
  
  ggsave(paste("LOR", Type.of.bias, "MSEPlot", ".pdf", sep = ""), a, dpi = 800, device = "pdf", width = 1.5*8.01, height = 2.25*5.67, units = "in")
  ggsave(paste("LOR", Type.of.bias, "MSEPlot", ".eps", sep = ""), a, dpi = 800, device = "eps", width = 1.5*8.01, height = 2.25*5.67, units = "in")
  
  
  # Coverage
  
  legend2 <- gtable_filter(ggplotGrob(p4), "guide-box")
  
  
  a <- grid.arrange(p4 + ylab(NULL) + theme(legend.position="none"),
                    p8 + ylab(NULL) + theme(legend.position="none"),
                    p12 + ylab(NULL) + theme(legend.position="none"),
                    p16 + ylab(NULL) + theme(legend.position="none"),
                    p20 + ylab(NULL) + theme(legend.position="none"),
                    legend2,
                    nrow = 3, ncol = 2,
                    left = textGrob("Coverage", rot = 90, vjust = 1)
  )
  
  ggsave(paste("LOR", Type.of.bias, "CoveragePlot", ".pdf", sep = ""), a, dpi = 800, device = "pdf", width = 1.5*8.01, height = 2.25*5.67, units = "in")
  ggsave(paste("LOR", Type.of.bias, "CoveragePlot", ".eps", sep = ""), a, dpi = 800, device = "eps", width = 1.5*8.01, height = 2.25*5.67, units = "in")
  
}
