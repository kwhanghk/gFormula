####################################################################################
#   A Set of Functions for
#    "G-formula with EAGeR data" 
####################################################################################

# install required packages
list.of.packages <- c("ranger", "glmnet", "gbm")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(ranger)
library(glmnet)
library(gbm)

####################################################################################

# Function to create (\bar{X}_t, \bar{A}_t)
buildCovMat <- 
  function(t, id.t, dat, x.cov.bs, x.cov.td, td.cov.names = NULL, append = FALSE, no.last.tdv = FALSE, maxT)
  {
    # DESCRYPTION:
    #    Returns cbind((\bar{X}_t, \bar{A}_t)
    #    which is the historical values of covariates and exposures
    #    of given id set up to time t
    
    # ARGUMENTS:
    #    id.t = good set of "id" available at time t
    
    # FUNCTION:  
    if (t<1) return(NULL)
    
    x.bs.t <- x.cov.bs[x.cov.bs$id %in% id.t, ]
    x.bs.t <- x.bs.t[order(x.bs.t$id),]; x.bs.t <- x.bs.t[,-1];
    if (is.null(td.cov.names)) {
      td.cov.names <- colnames(x.cov.td)[3:length(colnames(x.cov.td))]
    }
    x.td.t.pool <- x.cov.td[x.cov.td$id %in% id.t,]
    x.td.t.pool <- x.td.t.pool[order(x.td.t.pool$study_month, x.td.t.pool$id),]
    
    # define a local function that returns td cov values at time t_ (<t)
    td.time.t <- function(t_, td.cov.names) { 
      x.td.t = x.td.t.pool[x.td.t.pool$study_month == t_, td.cov.names]
      td.names = NULL
      for(nm in td.cov.names) td.names = c(td.names, paste(c(nm, t_), collapse = "_"))
      colnames(x.td.t) <- td.names
      return(x.td.t)
    }
    
    x.td.t.bar <- td.time.t(1, td.cov.names)
    if (t>1) {
      for (k in 2:t) {
        k.temp <- td.time.t(k, td.cov.names)
        x.td.t.bar <- cbind(x.td.t.bar, k.temp)
      }
    }
    
    if (append == T & t < maxT) { # force to repeat the final value until maxT (not used here)
      for (k in (t+1):maxT) {
        k.temp <- td.time.t(k, td.cov.names)
        x.td.t.bar <- cbind(x.td.t.bar, k.temp)
      }
    }
    
    if (no.last.tdv) x.td.t.bar <- x.td.t.bar[,-ncol(x.td.t.bar)]
    
    return(cbind(x.bs.t, x.td.t.bar))
    
  }

####################################################################################
####################################################################################

# Function to identify good set of 'id' at time t
goodIdSet <- 
  function(t, dat, lag = 0)
  {
    # DESCRYPTION:
    #    Returns the vector of id's that have never withdrawn
    #    until time t
    
    # ARGUMENTS:
    #    t & dat
    
    # FUNCTION:   
    if (lag == 0) {
      id.t <- dat[dat$study_month == t & dat$withdraw == 0,"id"] # never withdrawed until time t
    } else {
      id.t <- dat[dat$study_month == (t-lag) & dat$withdraw == 0,"id"] # never withdrawed until time t-lag
      if ((t-lag) <  1) {
        id.t <- union(dat[dat$study_month == t & dat$withdraw == 1,"id"]
                      , goodIdSet(t, dat))
      }
    }
                    
    # id.t <- setdiff(dat[dat$study_month == t & dat$withdraw == 0,"id"], # never withdrawed
    #                 dat[dat$study_month == t & dat$all.y == 0 & dat$conceived == 0,"id"] # not conceived before eof
    #                )
    
    return(unique(id.t))
    
  }

####################################################################################
####################################################################################

resample.list <-
  function(eff.ids) 
  {
    # DESCRYPTION: resample the original set of id
    
    sampled.rows <- sample(eff.ids, size=length(eff.ids),replace=TRUE)
    return(sampled.rows)
  }

####################################################################################
####################################################################################

resample.datMat <- 
  function(boot.list, datMat) 
  {
    # DESCRYPTION: return the dataset for bootstrapping corresponding to the id set of 
    # boot.list
    
    resampled.mat = NULL
    for (id.i in boot.list) {
      resampled.mat <- rbind(resampled.mat, datMat[datMat$id == id.i,])
    }
    return(resampled.mat)
  }

####################################################################################
####################################################################################

g.formula.outcome <- 
  function(times.vec, dat, x.cov.bs, x.cov.td, td.cov.names = NULL, no.last.tdv = FALSE, lag = 0, maxT, outcome,
           aBar, method = "RF")
  {
    # DESCRYPTION:
    #    Returns the pointwise g-formula estimate for all timepoints in the study
    #    
    
    # ARGUMENTS:
    #    times.vec {1,2, ... , T}
    #    dat
    #    x.cov.bs
    #    x.cov.td
    #    td.cov.names: names of the time dependent covariates to include; NULL if all
    #    no.last.tdv: do not use A_t if ture (i.e. only use (\bar{X_t}, \bar{A_{t-1}}))
    #    lag: time lag from the function goodIdSet
    #    maxT
    #    outcome: outcome of interest (Y)
    #    aBar = treatment sequence (vector for deterministic static intervention)
    #    method = regression method 
    
    # OUTPUT:
    #    two 1xT vectors (in the list) whose elements are the pointwise 
    #    g-formula estimates at each timepoint
    
    # FUNCTION:  
    
    ntimes <- max(times.vec)
    result.E.Ybar.a <- rep(0,ntimes)
    
    for (mT in times.vec) {
      
      outmod <- vector("list", mT)
      cat("\n", "fitting sequential regressions for T =", mT, "\n"); 
      sub.T.vec <- times.vec[times.vec==mT]:tail(times.vec, n=1)
      
      for (t in sub.T.vec){
        # good id set at time t:
        id.t <- goodIdSet(t, dat, lag) 
        cat("good id set (for fitting) at time", t, "=", length(id.t), "\n"); flush.console()
        if(length(id.t) == 0) {
          break;
        }
        
        # first we need to fit regression with this id set
        rspv <- dat[(dat$id %in% id.t) & dat$study_month==t, outcome]
        cat("sum of", outcome ,"at time", t, "=", sum(rspv), "\n"); flush.console()
        cov.dat.t <- buildCovMat(t, id.t, dat, x.cov.bs, x.cov.td, 
                                 td.cov.names = td.cov.names, append = FALSE, no.last.tdv = no.last.tdv, maxT = maxT)
        if (t == mT) {
          if (method == "RF") {
            outmod[[t]] <- ranger(rspv ~ ., data = cbind(cov.dat.t, rspv), write.forest = TRUE)
          }
          if (method == "GBM") {
            outmod[[t]] <-  gbm(rspv ~ ., data=cbind(cov.dat.t, rspv), shrinkage = 0.0005,
                                interaction.depth = 2, n.trees=2*10^3,  distribution = "bernoulli", verbose = F)
          }
          if (method == "lasso") { # TODO: convert all factor variables to numeric variables
            cov.dat.t$eligibility <- as.numeric(factor(cov.dat.t$eligibility)) - 1
            cov.dat.t <- as.matrix(cov.dat.t)
            outmod[[t]] <- tryCatch(
              cv.glmnet(cov.dat.t, y=rspv, alpha=1, family='binomial', type.measure="auc"),
              error = function(e) {
                NA
              }
            )
          }
        } else {
          if (method == "RF") {
            outmod[[t]] <- ranger(rspv.a ~ ., data = cbind(cov.dat.t, rspv.a), write.forest = TRUE)
          }
          if (method == "GBM") {
            outmod[[t]] <-  gbm(rspv.a ~ ., data=cbind(cov.dat.t, rspv.a), shrinkage = 0.0005,
                                interaction.depth = 2, n.trees=2*10^3,  distribution = "gaussian", verbose = F)
          }
          if (method == "lasso") {
            cov.dat.t$eligibility <- as.numeric(factor(cov.dat.t$eligibility)) - 1
            cov.dat.t <- as.matrix(cov.dat.t)
            outmod[[t]] <- tryCatch(
              cv.glmnet(cov.dat.t, y=rspv.a, alpha=1),
              error = function(e) {
                NA
              }
            )
          }
        }

        # Next, make a prediction with \bar{C}_{t-1} = 1 or 0
        id.t.pred <- if (t > 1) goodIdSet(t-1, dat, lag) else id.t
        cat("good id set (for predicting) at time", t, "=", length(id.t.pred), "\n"); flush.console()
        cov.dat.t.pred <- buildCovMat(t, id.t.pred, dat, x.cov.bs, x.cov.td, 
                                      td.cov.names = td.cov.names, append = FALSE, no.last.tdv = no.last.tdv, maxT = maxT)
        cov.dat.t.new <- cov.dat.t.pred  
        
        compl.j.idx = paste(c("compliance", t), collapse = "_")
        if (t > 1) {
          for (j in (t-1):1) {
            compl.j.idx <-  c(compl.j.idx, paste(c("compliance", j), collapse = "_"))
          }
        }
        cov.dat.t.new[,compl.j.idx] <- data.frame(t(aBar[length(compl.j.idx):1]))
          
        if (method == "RF") {
          Eya <- predict(outmod[[t]], data = cov.dat.t.new)$predictions
        }
        if (method == "GBM") {
          Eya <- predict(outmod[[t]], newdata=cov.dat.t.new, type="response", n.trees=2*10^3)
        }
        if (method == "lasso") {
          if(is.na(outmod[[t]][1])) {
            Eya <- ifelse(t==mT, rspv, rspv.a)
          } else {
            cov.dat.t.new$eligibility <- as.numeric(factor(cov.dat.t.new$eligibility)) - 1
            cov.dat.t.new <- as.matrix(cov.dat.t.new)
            Eya <- predict(outmod[[t]], newx=cov.dat.t.new, s = "lambda.min", type="response")
          }
        }
      
        rspv.a <- Eya
        
      }
      
      result.E.Ybar.a[mT] <- mean(rspv.a)

    }
    
    return(result.E.Ybar.a)
    
  }

####################################################################################
####################################################################################

g.formula.outcome.c.nc <- 
  function(times.vec, dat, x.cov.bs, x.cov.td, td.cov.names = NULL, no.last.tdv = FALSE, lag = 0, maxT, outcome,
           method = "RF")
  {
    # DESCRYPTION:
    #    Returns the pointwise g-formula estimate for all timepoints in the study
    #    
    
    # ARGUMENTS:
    #    times.vec {1,2, ... , T}
    #    dat
    #    x.cov.bs
    #    x.cov.td
    #    td.cov.names: names of the time dependent covariates to include; NULL if all
    #    no.last.tdv: do not use A_t if ture (i.e. only use (\bar{X_t}, \bar{A_{t-1}}))
    #    lag: time lag from the function goodIdSet
    #    maxT
    #    outcome: outcome of interest (Y)
    
    # OUTPUT:
    #    two 1xT vectors (in the list) whose elements are the pointwise 
    #    g-formula estimates at each timepoint
    
    # FUNCTION:  
    
    ntimes <- max(times.vec)
    result.E.Ybar.a0 <- result.E.Ybar.a1 <- rep(0,ntimes)
    
    for (mT in times.vec) {
      
      outmod.0 <- vector("list", mT); outmod.1 <- vector("list", mT);
      cat("\n", "fitting sequential regressions for T =", mT, "\n"); 
      sub.T.vec <- times.vec[times.vec==mT]:tail(times.vec, n=1)
      
      for (t in sub.T.vec){
        # good id set at time t:
        id.t <- goodIdSet(t, dat, lag) 
        cat("good id set (for fitting) at time", t, "=", length(id.t), "\n"); flush.console()
        if(length(id.t) == 0) {
          break;
        }
        
        # first we need to fit regression with this id set
        rspv <- dat[(dat$id %in% id.t) & dat$study_month==t, outcome]
        cat("sum of", outcome ,"at time", t, "=", sum(rspv), "\n"); flush.console()
        cov.dat.t <- buildCovMat(t, id.t, dat, x.cov.bs, x.cov.td, 
                                 td.cov.names = td.cov.names, append = FALSE, no.last.tdv = no.last.tdv, maxT = maxT)
        
        if (t == mT) {
          if (method == "RF") {
            outmod.1[[t]] <- outmod.0[[t]] <- ranger(rspv ~ ., data = cbind(cov.dat.t, rspv), 
                                                     write.forest = TRUE)
          }
          if (method == "GBM") {
            outmod.1[[t]] <- outmod.0[[t]] <-  gbm(rspv ~ ., data=cbind(cov.dat.t, rspv), shrinkage = 0.0005,
                                                   interaction.depth = 2, n.trees=5*10^3, cv.folds = 5, 
                                                   distribution = "bernoulli", verbose = F)
          }
          if (method == "lasso") {
            cov.dat.t$eligibility <- as.numeric(factor(cov.dat.t$eligibility)) - 1
            cov.dat.t <- as.matrix(cov.dat.t)
            outmod.1[[t]] <- outmod.0[[t]] <- tryCatch(
              cv.glmnet(cov.dat.t, y=rspv, alpha=0.5, family='binomial', type.measure="deviance", nfolds = 5),
              error = function(e) {
                NA
              }
            )
          }
        } else {
          if (method == "RF") {
            outmod.1[[t]]  <- ranger(rspv.1 ~ ., data = cbind(cov.dat.t, rspv.1), write.forest = TRUE)
            outmod.0[[t]]  <- ranger(rspv.0 ~ ., data = cbind(cov.dat.t, rspv.0), write.forest = TRUE)
          }
          if (method == "GBM") {
            outmod.1[[t]] <-  gbm(rspv.1 ~ ., data=cbind(cov.dat.t, rspv.1), shrinkage = 0.0005, cv.folds = 5,
                                  interaction.depth = 2, n.trees=5*10^3,  distribution = "gaussian", verbose = F)
            outmod.0[[t]] <-  gbm(rspv.0 ~ ., data=cbind(cov.dat.t, rspv.0), shrinkage = 0.0005, cv.folds = 5,
                                  interaction.depth = 2, n.trees=5*10^3,  distribution = "gaussian", verbose = F)
          }
          if (method == "lasso") {
            cov.dat.t$eligibility <- as.numeric(factor(cov.dat.t$eligibility)) - 1
            cov.dat.t <- as.matrix(cov.dat.t)
            outmod.1[[t]] <- tryCatch(
              cv.glmnet(cov.dat.t, y=rspv.1, alpha=0.5, type.measure="deviance", nfolds = 5),
              error = function(e) {
                NA
              }
            )
            outmod.0[[t]] <- tryCatch(
              cv.glmnet(cov.dat.t, y=rspv.0, alpha=0.5, type.measure="deviance", nfolds = 5),
              error = function(e) {
                NA
              }
            )
          }
        }
 
        # Next, make a prediction with \bar{C}_{t-1} = 1 or 0
        id.t.pred <- if (t > 1) goodIdSet(t-1, dat, lag) else id.t
        cat("good id set (for predicting) at time", t, "=", length(id.t.pred), "\n"); flush.console()
        cov.dat.t.pred <- buildCovMat(t, id.t.pred, dat, x.cov.bs, x.cov.td, 
                                      td.cov.names = td.cov.names, append = FALSE, no.last.tdv = no.last.tdv, maxT = maxT)
        cov.dat.t.new1 <- cov.dat.t.new0 <- cov.dat.t.pred  
        
        compl.j.idx = paste(c("compliance", t), collapse = "_")
        if (t > 1) {
          for (j in (t-1):1) {
            compl.j.idx <-  c(compl.j.idx, paste(c("compliance", j), collapse = "_"))
          }
        }
        cov.dat.t.new1[,compl.j.idx] <- 1 # complete compliance
        cov.dat.t.new0[,compl.j.idx] <- 0 # complete non-compliance
        
        if (method == "RF") {
          Eya1 <- predict(outmod.1[[t]], data = cov.dat.t.new1)$predictions
          Eya0 <- predict(outmod.0[[t]], data = cov.dat.t.new0)$predictions
        }
        if (method == "GBM") {
          best.tree.1 <- gbm.perf(outmod.1[[t]]); best.tree.0 <- gbm.perf(outmod.0[[t]]);
          Eya1 <- predict(outmod.1[[t]], newdata=cov.dat.t.new1, type="response", n.trees=best.tree.1)
          Eya0 <- predict(outmod.0[[t]], newdata=cov.dat.t.new0, type="response", n.trees=best.tree.0)
        }
        if (method == "lasso") {
          if(is.na(outmod.1[[t]][1])) {
            Eya1 <- ifelse(t==mT, rspv, rspv.1)
            #Eya1 <- rspv
          } else {
            cov.dat.t.new1$eligibility <- as.numeric(factor(cov.dat.t.new1$eligibility)) - 1
            cov.dat.t.new1 <- as.matrix(cov.dat.t.new1)
            Eya1 <- predict(outmod.1[[t]], newx=cov.dat.t.new1, s = "lambda.min", type="response")
          }
          if(is.na(outmod.0[[t]][1])) {
            Eya0 <- ifelse(t==mT, rspv, rspv.0)
            #Eya0 <- rspv
          } else {
            cov.dat.t.new0$eligibility <- as.numeric(factor(cov.dat.t.new0$eligibility)) - 1
            cov.dat.t.new0 <- as.matrix(cov.dat.t.new0)
            Eya0 <- predict(outmod.0[[t]], newx=cov.dat.t.new0, s = "lambda.min", type="response")
          }
        }
        
        rspv.1 <- Eya1
        rspv.0 <- Eya0

      }
      
      result.E.Ybar.a1[mT] <- mean(rspv.1)
      result.E.Ybar.a0[mT] <- mean(rspv.0)
    }
    
    result <- list(result.E.Ybar.a1, result.E.Ybar.a0)
    names(result) <- c("E.Ybar.a1", "E.Ybar.a0")
    return(result)

  }

####################################################################################
####################################################################################

bootstrap.plot <- 
  function(times.vec, result.E.Ybar.a.B, ylim.up = 1, outcome, alpha, method = "pivotal", lgd.position = "topleft")
  {
    # DESCRYPTION:
    #    Returns plots with bootstrapped trajectories 
    #     - use pivotal CI
    
    
    # ARGUMENTS:
    #    ylim.up: upper limit of Y
    #    alpha: we want to compute (1 - alpha)x100% CI
    
    # FUNCTION:    
    
    N = dim(result.E.Ybar.a.B)[1]
    quartz()
    for (k in 1:N) {
      thick <- if(k==1) 3 else 0.1
      col <- if(k==1) "Green" else rgb(0, 1, 0, 0.3)
      plot(result.E.Ybar.a.B[k,], type = "l", lwd = thick, col = col, ylim = c(0,ylim.up),  xlab = "T", 
           ylab=paste("Probability of ", outcome))
      par(new=T)
    }
    par(new=F)
    
    if (method == "pivotal") {
      sd.a <- apply(result.E.Ybar.a.B[-1,], 2, quantile, probs = c(alpha/2, 1-alpha/2))
      E.Ybar.a.u <- 2*result.E.Ybar.a.B[1,] - sd.a[1,]
      E.Ybar.a.l <- 2*result.E.Ybar.a.B[1,] - sd.a[2,]
    } 
    if (method == "percentile") {
      sd.a <- apply(result.E.Ybar.a.B, 2, quantile, probs = c(alpha/2, 1-alpha/2))
      E.Ybar.a.u <- sd.a[2,]
      E.Ybar.a.l <- sd.a[1,]
    }
    if (method == "studentized") { # TODO
      E.Ybar.a.u <- sd.a[2,]
      E.Ybar.a.l <- sd.a[1,]
    }
   
    quartz()
    plot(E.Ybar.a.u, type = "n", lty=2, lwd = 2, col = "Green", ylim = c(0,ylim.up), xlab = "T",
         ylab=paste("Probability of ", outcome)); par(new=T)
    plot(E.Ybar.a.l, type = "n", lty=2, lwd = 2, col = "Green", ylim = c(0,ylim.up), xlab = "T",
         ylab=paste("Probability of ", outcome)); par(new=T)
    polygon(c(rev(times.vec), times.vec), c(E.Ybar.a.u, rev(E.Ybar.a.l)), col=rgb(0, 1, 0, 0.4), border = NA)
    par(new=T)
    plot(result.E.Ybar.a.B[1,], type = "l", lwd = 3, col = "Green", ylim = c(0,ylim.up),  xlab = "T", 
         ylab=paste("Probability of ", outcome))
    
    mat.ci <- rbind(E.Ybar.a.u,E.Ybar.a.l)
    rownames(mat.ci) <- c("CI upper", "CI lower")
    print(mat.ci)
    return(mat.ci)
    
  }

####################################################################################
####################################################################################

bootstrap.plot.c.nc <- 
  function(times.vec, result.E.Ybar.a1.B, result.E.Ybar.a0.B, ylim.up = 1, outcome, alpha, lgd.position = "topleft")
  {
    # DESCRYPTION:
    #    Returns plots with bootstrapped trajectories 
    #     - use pivotal CI
    
    
    # ARGUMENTS:
    #    ylim.up: upper limit of Y
    #    alpha: we want to compute (1 - alpha)x100% CI
    
    # FUNCTION:    
    
    N = dim(result.E.Ybar.a1.B)[1]
    quartz()
    for (k in 1:N) {
      thick <- if(k==1) 3 else 0.1
      col1 <- if(k==1) "Blue" else rgb(0, 0, 1, 0.3)
      col0 <- if(k==1) "Red" else rgb(1, 0, 0, 0.3)
      plot(result.E.Ybar.a1.B[k,], type = "l", lwd = thick, col = col1, ylim = c(0,ylim.up), xlab = "T", 
           ylab=paste("Probability of ", outcome))
      par(new=T)
      plot(result.E.Ybar.a0.B[k,], type = "l", lwd = thick, col = col0, ylim = c(0,ylim.up),  xlab = "T", 
           ylab=paste("Probability of ", outcome))
      legend(lgd.position, c("Compliance","Non-compliance"), # puts text in the legend
             lty=c(1,1), # gives the legend appropriate symbols (lines)
             lwd=c(3,3),col=c("blue","red")) # gives the legend lines the correct color and width
      par(new=T)
    }
    par(new=F)
    
    sd.a1 <- apply(result.E.Ybar.a1.B, 2, quantile, probs = c(alpha/2, 1-alpha/2))
    sd.a0 <- apply(result.E.Ybar.a0.B, 2, quantile, probs = c(alpha/2, 1-alpha/2))
    
    E.Ybar.a1.u <- 2*result.E.Ybar.a1.B[1,] - sd.a1[1,]
    E.Ybar.a1.d <- 2*result.E.Ybar.a1.B[1,] - sd.a1[2,]
    E.Ybar.a0.u <- 2*result.E.Ybar.a0.B[1,] - sd.a0[1,]
    E.Ybar.a0.d <- 2*result.E.Ybar.a0.B[1,] - sd.a0[2,]
    
    quartz()
    plot(E.Ybar.a1.u, type = "n", lty=2, lwd = 2, col = "Blue", ylim = c(0,ylim.up), xlab = "T",
         ylab=paste("Probability of ", outcome)); par(new=T)
    plot(E.Ybar.a1.d, type = "n", lty=2, lwd = 2, col = "Blue", ylim = c(0,ylim.up), xlab = "T",
         ylab=paste("Probability of ", outcome)); par(new=T)
    polygon(c(rev(times.vec), times.vec), c(E.Ybar.a1.u, rev(E.Ybar.a1.d)), col=rgb(0, 0, 1, 0.4), border = NA)
    par(new=T)
    plot(E.Ybar.a0.u, type = "n", lty=2, lwd = 2, col = "Red", ylim = c(0,ylim.up), xlab = "T",
         ylab=paste("Probability of ", outcome)); par(new=T)
    plot(E.Ybar.a0.d, type = "n", lty=2, lwd = 2, col = "Red", ylim = c(0,ylim.up), xlab = "T",
         ylab=paste("Probability of ", outcome)); par(new=T)
    polygon(c(rev(times.vec), times.vec), c(E.Ybar.a0.u, rev(E.Ybar.a0.d)), col=rgb(1, 0, 0, 0.4), border = NA)
    par(new=T)
    plot(result.E.Ybar.a1.B[1,], type = "l", lwd = 3, col = "Blue", ylim = c(0,ylim.up), xlab = "T", 
         ylab=paste("Probability of ", outcome))
    par(new=T)
    plot(result.E.Ybar.a0.B[1,], type = "l", lwd = 3, col = "Red", ylim = c(0,ylim.up),  xlab = "T", 
         ylab=paste("Probability of ", outcome))
    legend(lgd.position, c("Compliance","Non-compliance"), 
           lty=c(1,1), 
           lwd=c(2.5,2.5),col=c("blue","red")) 
    
    mat.ci.c <- rbind(E.Ybar.a1.u, E.Ybar.a1.d)
    mat.ci.nc <- rbind(E.Ybar.a0.u, E.Ybar.a0.d)
    rownames(mat.ci.c) <- c("CI upper for compliance", "CI lower for compliance")
    rownames(mat.ci.nc) <- c("CI upper for non-compliance", "CI lower for non-compliance")
    print(mat.ci.c); print(mat.ci.nc);
    
    return(list(mat.ci.c, mat.ci.nc))
    
  }