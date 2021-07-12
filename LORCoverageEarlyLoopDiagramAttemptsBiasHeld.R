### Remove previous variables
rm(list = ls())

nestedloop <- function(x,
                       varnames, sign=rep(1, length(varnames)),
                       varlabels=NULL){
  ##
  if (!inherits(x, "data.frame"))
    stop("Argument 'x' must be a data.frame.")
  ##
  mo <- matrix(sign,
               nrow=dim(x)[[1]], ncol=length(varnames),
               byrow=TRUE)
  xo <- x[,varnames]
  ##
  ## Re-ordering:
  res <- x[do.call(order, mo*xo),]
  ##
  attr(res, "varnames")  <- varnames
  attr(res, "varlabels") <- varlabels
  attr(res, "sign") <- sign
  ##
  class(res) <- c("nestedloop", class(res))
  ##
  res
}



lines.nestedloop <- function(x,
                             varnames=attributes(x)$varnames,
                             varlabels=attributes(x)$varlabels,
                             which="v",
                             col=if (which=="r") "#999999" else "black",
                             ymin.refline, ymax.refline,
                             cex.ref=0.9,
                             log=TRUE,
                             ...){
  ##
  nvar <- length(varnames)
  ##
  if (length(col)==1)
    col <- rep(col, nvar)
  ##
  if (which=="v"){
    ##
    ## Vertical lines
    ##
    nlen <- rep(NA, nvar)
    ##
    for (i in 1:nvar)
      nlen[i] <- length(unique(x[,varnames[i]]))
    ##
    cnlen <- cumprod(nlen)
    ##
    for (i in (nvar-1):1)
      abline(v=cnlen[nvar]*(0:cnlen[i])/cnlen[i]+1,
             col=col[i])
  }
  else if (which=="r"){
    ##
    ## Reference lines
    ##
    if (is.null(varlabels))
      varlabels <- varnames
    ##
    labels.varnames <- rep("", nvar)
    ##
    for (i in 1:length(varnames)){
      if (is.factor(x[,varnames[i]]))
        varvals <- unique(x[,varnames[i]])
      else{
        varvals <- format(unique(x[,varnames[i]]))
        varvals <- sub("^[[:space:]]*(.*?)[[:space:]]*$",
                       "\\1",
                       varvals,
                       perl=TRUE)
      }
      ##
      labels.varnames[i] <- paste(varlabels[i],
                                  " (",
                                  paste(varvals,
                                        collapse=", "),
                                  ")", sep="")
    }
    ##
    if (log){
      ymax <- log(ymax.refline)
      ymin <- log(ymin.refline)
    }
    else{
      ymax <- ymax.refline
      ymin <- ymin.refline
    }
    ##
    distance <- (ymax-ymin)/nvar
    ##
    ypos <- ymax-0.2*distance-(1/nvar)*(0:(nvar-1))*(ymax-ymin)
    ypos.ref.max <- ypos-0.20*distance
    ypos.ref.min <- ypos-0.75*distance
    ##
    if (log){
      ypos <- exp(ypos)
      ypos.ref.max <- exp(ypos.ref.max)
      ypos.ref.min <- exp(ypos.ref.min)
    }
    ##
    for (i in 1:nvar){
      ##print(c(ypos.ref.max[i], ypos.ref.min[i]))
      ##
      ##abline(h=ypos[i], lwd=1, col="red")
      ##abline(h=ypos.ref.max[i], lwd=1, col="blue")
      ##abline(h=ypos.ref.min[i], lwd=1, col="green")
      text(1, ypos[i], labels.varnames[i], adj=0, cex=cex.ref)
      ##
      xvar <- x[,varnames[i]]
      if (is.factor(xvar))
        xvar <- as.numeric(xvar)
      xvar <- ypos.ref.min[i] +
        (xvar-min(xvar))/(max(xvar)-min(xvar))*
        (ypos.ref.max[i]-ypos.ref.min[i])
      lines(xvar, col=col[i], type="s", lwd=1)
    }
  }
  ##
  invisible(NULL)
}


panel.nestedloop <- function(x,
                             varnames=attributes(x)$varnames,
                             varlabels=attributes(x)$varlabels,
                             which="v",
                             col=if (which=="r") "#999999" else "black",
                             ymin.refline, ymax.refline,
                             cex.ref=0.9,
                             log=TRUE,
                             ...){
  ##
  nvar <- length(varnames)
  ##
  if (length(col)==1)
    col <- rep(col, nvar)
  ##
  if (which=="v"){
    ##
    ## Vertical lines
    ##
    nlen <- rep(NA, nvar)
    ##
    for (i in 1:nvar)
      nlen[i] <- length(unique(x[,varnames[i]]))
    ##
    cnlen <- cumprod(nlen)
    ##
    for (i in (nvar-1):1)
      abline(v=cnlen[nvar]*(0:cnlen[i])/cnlen[i]+1,
             col=col[i])
  }
  else if (which=="r"){
    ##
    ## Reference lines
    ##
    if (is.null(varlabels))
      varlabels <- varnames
    ##
    labels.varnames <- rep("", nvar)
    ##
    for (i in 1:length(varnames)){
      if (is.factor(x[,varnames[i]]))
        varvals <- unique(x[,varnames[i]])
      else{
        varvals <- format(unique(x[,varnames[i]]))
        varvals <- sub("^[[:space:]]*(.*?)[[:space:]]*$",
                       "\\1",
                       varvals,
                       perl=TRUE)
      }
      ##
      labels.varnames[i] <- paste(varlabels[i],
                                  " (",
                                  paste(varvals,
                                        collapse=", "),
                                  ")", sep="")
    }
    ##
    if (log){
      ymax <- log(ymax.refline)
      ymin <- log(ymin.refline)
    }
    else{
      ymax <- ymax.refline
      ymin <- ymin.refline
    }
    ##
    distance <- (ymax-ymin)/nvar
    ##
    ypos <- ymax-0.2*distance-(1/nvar)*(0:(nvar-1))*(ymax-ymin)
    ypos.ref.max <- ypos-0.20*distance
    ypos.ref.min <- ypos-0.75*distance
    ##
    if (log){
      ypos <- exp(ypos)
      ypos.ref.max <- exp(ypos.ref.max)
      ypos.ref.min <- exp(ypos.ref.min)
    }
    ##
    for (i in 1:nvar){
      ##print(c(ypos.ref.max[i], ypos.ref.min[i]))
      ##
      ##abline(h=ypos[i], lwd=1, col="red")
      ##abline(h=ypos.ref.max[i], lwd=1, col="blue")
      ##abline(h=ypos.ref.min[i], lwd=1, col="green")
      ltext(1, ypos[i], labels.varnames[i], adj=c(0, 0.5), cex=cex.ref)
      ##
      xvar <- x[,varnames[i]]
      if (is.factor(xvar))
        xvar <- as.numeric(xvar)
      xvar <- ypos.ref.min[i] +
        (xvar-min(xvar))/(max(xvar)-min(xvar))*
        (ypos.ref.max[i]-ypos.ref.min[i])
      llines(xvar, col=col[i], type="s", lwd=1)
    }
  }
  ##
  invisible(NULL)
}


### Attempting nested loop ----

### Set up ----

library(data.table)


#### Easy access to values ----
# UMD

setwd("D:/LOR")

Subj <- c("None", "Step", "Method", "AltOut")
Ref <- c(0.1, 0.3, 0.5)

for(a in Subj){
  for (b in Ref){
  Sum_results2 <- data.table(read.csv("LogORResults_Looped.csv"))
  # Sum_results <- Sum_results2[, `:=` (Bias_FE = Bias_FE + Rep_theta,
  #                                     Bias_REML = Bias_REML + Rep_theta,
  #                                     Bias_DL = Bias_DL + Rep_theta,
  #                                     Bias_Moreno = Bias_Moreno + Rep_theta)]
  # 
  # Sum_results <- data.frame(Sum_results[Bias_type == "None"])
  
  #Sum_results$Bias_type <- ordered(Sum_results$Bias_type, levels = c("None", "Step", "ModBegg", "Method", "AltOut", "Out"))
  #Sum_results$Bias_typenum <- as.numeric(Sum_results$Bias_type)
  
  
  Sum_results <- Sum_results2[, `:=` (RMSE_FE = sqrt(MSE_FE),
                                      RMSE_REML = sqrt(MSE_REML),
                                      RMSE_DL = sqrt(MSE_DL),
                                      RMSE_Moreno = sqrt(MSE.Moreno))]
  
  Sum_results <- data.frame(Sum_results[Bias_type == a & Rep_ev_freq == b])
  
  
  # nldata <- nestedloop(Sum_results,
  #                      varnames=c("Rep_NumStudies", "Rep_tau.sq","Rep_Subj", "Rep_theta", "Rep_ev_freq"),
  #                      varlabels=
  #                        c("Number of studies" , 
  #                          "Heterogeneity",
  #                          "Study size", "LOR", 
  #                          "Event Rate"
  #                        ),
  #                      sign=c(1, 1, 1, 1))
  
  nldata <- nestedloop(Sum_results,
                       varnames=c("Rep_theta", "Rep_NumStudies", "Rep_tau.sq","Rep_Subj"),
                       varlabels=
                         c("LOR", "Number of studies" , 
                           "Heterogeneity",
                           "Study size"
                         ),
                       sign=c(1, 1, 1, 1))
  
  
  ##
  ## R object pd2 is used in nested-loop plot (Figure 2)
  ##
  pd2 <- nldata
  ##
  ## Use labels instead of numeric values for rho2 and p.c
  ##
  # pd2$rho2 <- factor(pd2$Bias_typenum,
  #                    levels=c(1,2,3,4,5,6),
  #                    labels=c("None", "Step", "ModBegg", "Method", "AltOut", "Out"))
  pd2$Rep_Subj <- factor(pd2$Rep_Subj,
                         levels=c(4.7, 20, 100, 250),
                         labels=c("Empirical", "Small", "Fixed at median", "Large"))
  
  
  pdf(paste0("EvFreq ", b, " LOR_Coverage_Plot_Bias_",a,".pdf"), paper="a4r", width=18, height=15)
  ##
  par(pty="m")
  ##
  ## Create skeleton of nested-loop plot using standard R plot function
  ##
  plot(pd2$Rep_theta,
       #log=c("y"), 
       type="n",
       ylim=c(0.5, 1.05), bty="n",
       xlab="5 x 4 x 4 x 6 = 480 ordered scenarios",
       ylab="Coverage",
       las=1, xaxt="n")
  ##
  ## Add vertical lines (using R function lines.nestedloop)
  ##
  lines(pd2, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
  ##
  ## Add reference lines (using R function lines.nestedloop)
  ##
  lines(pd2, which="r",
        ymin.refline=1, ymax.refline=1.07,
        cex.ref=0.7)
  
  
  ### Adding alpha
  library(RColorBrewer)
  
  colour1 <- brewer.pal(5, "Dark2")
  
  add.alpha <- function(col, alpha=1){
    if(missing(col))
      stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, 
          function(x) 
            rgb(x[1], x[2], x[3], alpha=alpha))  
  }
  
  colour2 <- add.alpha(colour1, 1)
  
  ##
  ## Estimates and legend (using standard R functions lines and legend)
  ##
  abline(h = 0.95, untf = FALSE, lty = "dotted")
  #lines(y=0, col="black", lwd=2, type="s")   # True theta
  #lines(pd2$Coverage_FE, col=colour2[1], type="s", lwd = 0.65)    # Peto
  lines(pd2$Coverage_REML, col=colour2[1], type="s", lwd = 0.75)   # Trimfill
  lines(pd2$Coverage_IVHet, col=colour2[2], type="s", lwd = 0.75)      # Peters
  lines(pd2$Coverage_HC_DL, col=colour2[3], type="s", lwd = 0.75)           # Limit meta-analysis, method 1
  lines(pd2$Coverage_KH_REML, col=colour2[4], type="s", lwd = 0.75)      # Limit meta-analysis, method 2
  lines(pd2$Coverage_Mult, col=colour2[5], type="s", lwd = 0.75) 
  
  #lines(pd2$Coverage_Moreno, col=colour2[6], type="s", lwd = 0.65)     # Limit meta-analysis, method 3
  
  
  ##
  legend(1, 0.63,
         lwd=c(2, rep(1, 6)),
         col=c(colour1),
         cex=0.9,
         bty="n",
         c("REML", "IVHet", "HC DL", "KH REML", "Mult"))
  ##
  dev.off()
}
}

#### Simpler diagram

for(a in Subj){
  for (b in Ref){
  Sum_results2 <- data.table(read.csv("LogORResults_Looped.csv"))
  # Sum_results <- Sum_results2[, `:=` (Bias_FE = Bias_FE + Rep_theta,
  #                                     Bias_REML = Bias_REML + Rep_theta,
  #                                     Bias_DL = Bias_DL + Rep_theta,
  #                                     Bias_Moreno = Bias_Moreno + Rep_theta)]
  # 
  # Sum_results <- data.frame(Sum_results[Bias_type == "None"])
  
  #Sum_results$Bias_type <- ordered(Sum_results$Bias_type, levels = c("None", "Step", "ModBegg", "Method", "AltOut", "Out"))
  #Sum_results$Bias_typenum <- as.numeric(Sum_results$Bias_type)
  
  
  Sum_results <- Sum_results2[, `:=` (RMSE_FE = sqrt(MSE_FE),
                                      RMSE_REML = sqrt(MSE_REML),
                                      RMSE_DL = sqrt(MSE_DL),
                                      RMSE_Moreno = sqrt(MSE.Moreno))]
  
  Sum_results <- data.frame(Sum_results[Bias_type == a & Rep_ev_freq == b])
  
  
  # nldata <- nestedloop(Sum_results,
  #                      varnames=c("Rep_NumStudies", "Rep_tau.sq","Rep_Subj", "Rep_theta", "Rep_ev_freq"),
  #                      varlabels=
  #                        c("Number of studies" , 
  #                          "Heterogeneity",
  #                          "Study size", "LOR", 
  #                          "Event Rate"
  #                        ),
  #                      sign=c(1, 1, 1, 1))
  
  nldata <- nestedloop(Sum_results,
                       varnames=c("Rep_NumStudies", "Rep_theta", "Rep_tau.sq","Rep_Subj"),
                       varlabels=
                         c("Number of studies",  "LOR", 
                           "Heterogeneity",
                           "Study size"
                         ),
                       sign=c(1, 1, 1, 1))
  
  
  ##
  ## R object pd2 is used in nested-loop plot (Figure 2)
  ##
  pd2 <- nldata
  ##
  ## Use labels instead of numeric values for rho2 and p.c
  ##
  # pd2$rho2 <- factor(pd2$Bias_typenum,
  #                    levels=c(1,2,3,4,5,6),
  #                    labels=c("None", "Step", "ModBegg", "Method", "AltOut", "Out"))
  pd2$Rep_Subj <- factor(pd2$Rep_Subj,
                         levels=c(4.7, 20, 100, 250),
                         labels=c("Empirical", "Small", "Fixed at median", "Large"))
  
  
  pdf(paste0("Alt EvFreq ", b, " LOR_Coverage_Plot_Bias_", a,".pdf"), paper="a4r", width=18, height=15)
  ##
  par(pty="m")
  ##
  ## Create skeleton of nested-loop plot using standard R plot function
  ##
  plot(pd2$Rep_theta,
       #log=c("y"), 
       type="n",
       ylim=c(0.5, 1.05), bty="n",
       xlab="5 x 4 x 4 x 6 = 480 ordered scenarios",
       ylab="Coverage",
       las=1, xaxt="n")
  ##
  ## Add vertical lines (using R function lines.nestedloop)
  ##
  lines(pd2, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
  ##
  ## Add reference lines (using R function lines.nestedloop)
  ##
  lines(pd2, which="r",
        ymin.refline=0.98, ymax.refline=1.07,
        cex.ref=0.7)
  
  
  ### Adding alpha
  library(RColorBrewer)
  
  colour1 <- brewer.pal(5, "Dark2")
  
  add.alpha <- function(col, alpha=1){
    if(missing(col))
      stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, 
          function(x) 
            rgb(x[1], x[2], x[3], alpha=alpha))  
  }
  
  colour2 <- add.alpha(colour1, 1)
  
  ##
  ## Estimates and legend (using standard R functions lines and legend)
  ##
  abline(h = 0.95, untf = FALSE, lty = "dotted")
  #lines(y=0, col="black", lwd=2, type="s")   # True theta
  #lines(pd2$Coverage_FE, col=colour2[1], type="s", lwd = 0.65)    # Peto
  lines(pd2$Coverage_REML, col=colour2[1], type="s", lwd = 0.75)   # Trimfill
  lines(pd2$Coverage_IVHet, col=colour2[2], type="s", lwd = 0.75)      # Peters
  lines(pd2$Coverage_HC_DL, col=colour2[3], type="s", lwd = 0.75)           # Limit meta-analysis, method 1
  lines(pd2$Coverage_KH_REML, col=colour2[4], type="s", lwd = 0.75)      # Limit meta-analysis, method 2
  lines(pd2$Coverage_Mult, col=colour2[5], type="s", lwd = 0.75) 
  
  #lines(pd2$Coverage_Moreno, col=colour2[6], type="s", lwd = 0.65)     # Limit meta-analysis, method 3
  
  
  ##
  legend(1, 0.63,
         lwd=c(2, rep(1, 6)),
         col=c(colour1),
         cex=0.9,
         bty="n",
         c("REML", "IVHet", "HC DL", "KH REML", "Mult"))
  ##
  dev.off()
  }
}


## Reorder dataset with simulation results according to
## simulation parameters theta, rho2, p.c, tau2, k.
##
## Simulation parameters theta, p.c and k are sorted in
## decreasing order (see argument sign).
##
# nldata <- nestedloop(res,
#                      varnames=c("theta", "rho2", "p.c", "tau2", "k"),
#                      varlabels=
#                        c("Odds ratio", "Selection",
#                          "Control event proportion", "Heterogeneity",
#                          "Number of studies"),
#                      sign=c(-1, 1, -1, 1, -1))