#####################################################################################  
#                   R Code for "G-formula with longitudinal data"                   #                
#                   Kwhangho Kim (kwhanghk@stat.cmu.edu)                            #
#                                                                                   #
#                                                      Final Edit : June/20/2017    #
#####################################################################################  


# -----------------------------------------------------------------------------------
# data structure: assume our data storage follows 
# (id, time, last, covariates, W, Y) long form longitudinal right censored structure
# where W: censorship variable ; Y: outcome variables of interest ; 
#       last : flag for last observation
#
# * No missing data assumed
# example for illustration: EAGeR data modified by Ashley
# -----------------------------------------------------------------------------------


# setwd("C:/Users/J Kim/Dropbox/R_Working/2017/ADA")
# setwd("~/Documents/R_Working/ADA")
# source("function_gFormula_0620.R")



#####################################################################################  
# 1) Point estimate of mean function of Y when we assign a sequence of treatment 
#    \bar{a}_T (deterministic static intervention)
##################################################################################### 


# -----------------------------------------------------------------------------------
# setup: manipulate data structure so that it becomes balanced across id
# intput: aspirin (raw data)
# output = aspirin (elongated), exclude list
# -----------------------------------------------------------------------------------

# 
dataStr <- 
  function(aspirin)
  {
    
    # TODO: time <- study_month
    #       last <- create last if there's no last   
    
    uniq_ids <- unique(aspirin$id)
    # length(uniq_ids)  #[1] 1228
    
    n0 <- length(unique(aspirin$id)) 
    all.T <- table(aspirin[aspirin$last==1, "study_month"])
    # barplot(cumsum(all.T))
    all.T.mat <- t(as.matrix(all.T))
    unique.times <- as.numeric(colnames(all.T.mat))
    ntimes <- max(unique.times) # max {T}
    times.vec <- rev(unique.times) # {T}
    
    # elongate each data up to the end of obsevation in the way that the appended data become
    # a constant equal to its final value after the time point at which no more 
    # data is collected on the subject (van der Laan & Robins, 2003 p.314)
    
    aspirin_ext <- aspirin
    for (i in uniq_ids) {
      last_dt <- aspirin[aspirin$id == i & aspirin$last == 1, ]
      if (last_dt$study_month == ntimes) {
        next
      } else {
        last_dt$last <- 0
        for (t in (last_dt$study_month+1):ntimes) {
          last_dt$study_month <- t
          aspirin_ext <- rbind(aspirin_ext, last_dt)
        }
      }
    }
    
    aspirin <- aspirin_ext[order(aspirin_ext$id, aspirin_ext$study_month),]
    
    
    if (dim(aspirin)[1] != n0 * ntimes) {
      # if inconsistency occurs due to missing data, remove them
      
      # type 1 irregularities: identify subject with missing values  
      # i.e. identify {id} with missing values before its final value
      exclude.list.1 = NULL
      for (m in times.vec) {
        id.m.vec <- aspirin[aspirin$last == 1 & aspirin$study_month == m, "id"]
        for (id.i in id.m.vec) {
          if (dim(aspirin[aspirin$id == id.i,])[1] != ntimes) {
            exclude.list.1 <- c(exclude.list.1, id.i)
          }
        }
      }
      exclude.list.1 <- unique(exclude.list.1)
      
      ## identify {id} whose outcome do not fall into one of the four categories until the end of the study
      # exclude.list.2 = NULL
      # for (m in times.vec) {
      #   id.m.vec <- aspirin[dat$last == 1 & dat$all.y != 1, "id"]
      #   for (id.i in id.m.vec) {
      #     if (dim(aspirin[aspirin$id == id.i,])[1] != m) {
      #       exclude.list.2 <- c(exclude.list.2, id.i)
      #     }
      #   }
      # }
      # exclude.list.2 <- unique(exclude.list.2)
      
      # exclude.list <- c(exclude.list.1, exclude.list.2)
      exclude.list <- c(exclude.list.1)
      cat("\n", length(exclude.list), " subjects are removed due to missing values in data", "\n");
      
    }
    
    dim(aspirin[!(aspirin$id %in% exclude.list),])[1]
    (n0 - length(exclude.list)) * ntimes # now the inconsistency disappears
    
    # Do we need all.y? 
    
    # id_lb = aspirin[aspirin$live_birth == 1 & aspirin$last == 1,1]
    # id_fl = aspirin[aspirin$fetal_loss == 1 & aspirin$last == 1,1]
    # id_wd = aspirin[aspirin$withdraw == 1 & aspirin$last == 1,1]
    # id_ewop = aspirin[aspirin$efuwp == 1 & aspirin$last == 1,1]
    # 
    # # only subjects whose outcome is either livebirth or fetalloss
    # # aspirin <- aspirin[aspirin$id %in% union(id_lb, id_fl),]
    # 
    # aspirin["all.y"] <- NA 
    # aspirin$all.y <- aspirin$fetal_loss + aspirin$efuwp + aspirin$live_birth + aspirin$withdraw
    
    
    # output: aspirin, exclude.list
    list.output <- list(aspirin, exclude.list)
    names(list.output) <- c("data.long", "exclude.list")
    return(list.output)
  }


# -----------------------------------------------------------------------------------
# Fiting G-formula with the appropriate data structure
# input: data info & covariates vectors; outcome; treatment sequence
# output = pointwise estimates via g-formula at each time point 
# -----------------------------------------------------------------------------------


gformula.detinvt <- 
  function(aspirin, dat.names, cov.bs.names, cov.td.names, 
           no.last.tdv = FALSE, lag = 0, outcome.v,
           exclude.list, 
           boot = FALSE, B = 1, 
           a.bar,
           reg.method = "RF")
  {
    
    dat = aspirin[,dat.names]; 
    x.cov.bs = aspirin[aspirin$study_month == 1, cov.bs.names]; 
    cat("\n", dim(x.cov.bs)[2] - 1, " baseline covariates detected", "\n");
    x.cov.td = aspirin[,cov.td.names]; # dim(x.cov.td)
    cat("\n", dim(x.cov.td)[2] - 2, " time-varying covariates detected", "\n");
    
    # effective ids
    uniq_ids <- unique(aspirin$id)
    eff.ids <- uniq_ids[!(uniq_ids %in% exclude.list)]
    
    # remove such irregularities
    dat = dat[dat$id %in% eff.ids,]; # dim(dat)
    x.cov.bs = x.cov.bs[x.cov.bs$id %in% eff.ids,]; # dim(x.cov.bs)
    x.cov.td = x.cov.td[x.cov.td$id %in% eff.ids,]; # dim(x.cov.td)
    
    pnt.est <- estimation.proc(dat, x.cov.bs, x.cov.td, 
                               outcome = outcome.v, 
                               td.cov.names = NULL, no.last.tdv = no.last.tdv,
                               eff.ids = eff.ids, boot = boot, B = B,
                               aBar = a.bar,
                               method = reg.method)
    
    cat("\n", "g-formula estimation for cumulative risk function of", outcome.v, ": ", "\n");
    print(pnt.est[1,])
    plot(pnt.est[1,], type = "l", lwd = 2.5, col = "orange", xlab = "Time", ylab = paste("Probability of", outcome.v))
    return(pnt.est)
    
  }



estimation.proc <- 
  function(dat, x.cov.bs, x.cov.td, td.cov.names = NULL, no.last.tdv = FALSE, lag = 0, outcome,
           eff.ids, boot = FALSE, B = 1, aBar, method = "RF")
  {
    
    all.T <- table(dat[dat$last==1, "study_month"])
    all.T.mat <- t(as.matrix(all.T))
    unique.times <- as.numeric(colnames(all.T.mat))
    ntimes <- max(unique.times) # max {T}
    times.vec <- rev(unique.times) # {T}
    
    result.E.Ybar.a.B = matrix(nrow = B, ncol = ntimes)
    
    for (iter in 1:B) {
      # to get only single time pointwise estimate, set boot = FALSE & B = 1
      if (iter > 1 || boot == TRUE) { # bootstrapped set
        cat("\n", "bootstrapping step: ", iter, "\n");
        boot.list <- sort(resample.list(eff.ids))
        dat.B <- resample.datMat(boot.list, dat)
        x.cov.bs.B <- resample.datMat(boot.list, x.cov.bs)
        x.cov.td.B <- resample.datMat(boot.list, x.cov.td)
      } else { # original set (default)
        dat.B <- dat
        x.cov.bs.B <- x.cov.bs
        x.cov.td.B <- x.cov.td
      }
      # Fit sequential regression for g-formula:
      result.E.Ybar <- g.formula.outcome(times.vec, dat.B, x.cov.bs.B, x.cov.td.B, td.cov.names = td.cov.names, 
                                         no.last.tdv = no.last.tdv, lag = lag,
                                         maxT= ntimes, outcome = outcome,
                                         aBar = aBar,
                                         method = method)
      result.E.Ybar.a.B[iter,] <- result.E.Ybar
    }
    
    return(result.E.Ybar.a.B)
    
  }




#####################################################################################  
# 2) CI estimates of the mean function with bootstrapping
#    when \bar{a}_T is assigned (deterministic static intervention)
##################################################################################### 
# 
# 
# # main loop for g-formula estimation and bootstrapping
# boot = FALSE
# n.B = 10 # to get only single time pointwise estimate, set boot = FALSE & B = 1
# 
# aspirin.boot.ex <- gformula.detinv(aspirin, dat.names, cov.bs.names, cov.td.names, 
#                                    outcome.v= outcome.v,
#                                    exclude.list=exclude.list, 
#                                    a.bar=a.bar,
#                                    B = n.B)
# 
# all.T <- table(aspirin[aspirin$last==1, "study_month"])
# all.T.mat <- t(as.matrix(all.T))
# unique.times <- as.numeric(colnames(all.T.mat))
# times.vec <- rev(unique.times) # {T}
# 
# result.E.Ybar.a.B = aspirin.boot.ex
# 
# bootstrap.plot(times.vec, result.E.Ybar.a.B, ylim.up = 1, outcome.v, alpha = 0.05, method = "percentile", lgd.position = "topleft")

bootstrap.ci <- 
  function(aspirin, dat.names, cov.bs.names, cov.td.names, 
           outcome.v, exclude.list, a.bar, reg.method = "RF", 
           B, 
           alpha = 0.05, method = "pivotal", lgd.position = "topleft", ylim.up = 1)
  {
    
    boot.data <- gformula.detinvt(aspirin, dat.names, cov.bs.names, cov.td.names, 
                                  outcome.v= outcome.v,
                                  exclude.list=exclude.list, 
                                  a.bar=a.bar,
                                  reg.method = reg.method,
                                  B = B)
    
    all.T <- table(aspirin[aspirin$last==1, 2])
    all.T.mat <- t(as.matrix(all.T))
    unique.times <- as.numeric(colnames(all.T.mat))
    times.vec <- rev(unique.times) # {T}
    
    ci <- bootstrap.plot(times.vec, boot.data, ylim.up = 1, outcome.v, alpha, method, lgd.position)
    
    return(ci)
  }


#####################################################################################  
# 3) Point estimate of mean function of Y with complete compliance and complete
#    non-complianece
##################################################################################### 


gformula.detinvt.c.nc <- 
  function(aspirin, dat.names, cov.bs.names, cov.td.names, 
           no.last.tdv = FALSE, lag = 0, outcome.v,
           exclude.list, 
           boot = FALSE, B = 1, 
           reg.method = "RF")
  {
    
    dat = aspirin[,dat.names]; 
    x.cov.bs = aspirin[aspirin$study_month == 1, cov.bs.names]; 
    cat("\n", dim(x.cov.bs)[2] - 1, " baseline covariates detected", "\n");
    x.cov.td = aspirin[,cov.td.names]; # dim(x.cov.td)
    cat("\n", dim(x.cov.td)[2] - 2, " time-varying covariates detected", "\n");
    
    # effective ids
    uniq_ids <- unique(aspirin$id)
    eff.ids <- uniq_ids[!(uniq_ids %in% exclude.list)]
    
    # remove such irregularities
    dat = dat[dat$id %in% eff.ids,]; # dim(dat)
    x.cov.bs = x.cov.bs[x.cov.bs$id %in% eff.ids,]; # dim(x.cov.bs)
    x.cov.td = x.cov.td[x.cov.td$id %in% eff.ids,]; # dim(x.cov.td)
    
    pnt.est <- estimation.proc.c.nc(dat, x.cov.bs, x.cov.td, 
                                    outcome = outcome.v, 
                                    td.cov.names = NULL, no.last.tdv = no.last.tdv,
                                    eff.ids = eff.ids, boot = boot, B = B,
                                    method = reg.method)
    result.E.Ybar.a1.B <- pnt.est[[1]]
    result.E.Ybar.a0.B <- pnt.est[[2]]
    cat("\n", "g-formula estimation for cumulative risk function of", outcome.v)
    cat("\n",  "for both perfect compliance and non-compliance: ", "\n")
    cat("\n",  "- Perfect compliance ", "\n")
    print(result.E.Ybar.a1.B[1,])
    cat("\n",  "- Perfect non-compliance ", "\n")
    print(result.E.Ybar.a0.B[1,])
    
    y.max <- max(max(result.E.Ybar.a1.B[1,]), max(result.E.Ybar.a0.B[1,]))
    
    # par(xpd = T, mar = par()$mar + c(0,0,0,7))
    plot(result.E.Ybar.a1.B[1,], type = "l", lwd = 2.5, col = "Blue", ylim = c(0,y.max), xlab = "Time", ylab = paste("Probability of", outcome.v))
    par(new=T)
    plot(result.E.Ybar.a0.B[1,], type = "l", lwd = 2.5, col = "Red", ylim = c(0,y.max),  xlab = "Time", ylab = paste("Probability of", outcome.v))
    # legend(0.03, 0.015,
    #        c("Compliance","Non-compliance"), 
    #        lty=c(1,1), 
    #        lwd=c(2.5,2.5),col=c("blue","red")) 
    # par(mar=c(5, 4, 4, 2) + 0.1)
    legend("bottomright", c("Compliance","Non-compliance"),
           lty=c(1,1),
           lwd=c(2.5,2.5),col=c("blue","red"))
    
    return(pnt.est)
    
  }

#
estimation.proc.c.nc <- 
  function(dat, x.cov.bs, x.cov.td, td.cov.names = NULL, no.last.tdv = FALSE, lag = 0, outcome,
           eff.ids, boot = FALSE, B = 1, method = "RF")
  {
    
    all.T <- table(dat[dat$last==1, "study_month"])
    all.T.mat <- t(as.matrix(all.T))
    unique.times <- as.numeric(colnames(all.T.mat))
    ntimes <- max(unique.times) # max {T}
    times.vec <- rev(unique.times) # {T}
    
    result.E.Ybar.a0.B = result.E.Ybar.a1.B = matrix(nrow = B, ncol = ntimes)
    
    for (iter in 1:B) {
      
      if (iter > 1 || boot == TRUE) { # bootstrapped set
        cat("\n", "bootstrapping step: ", iter, "\n");
        boot.list <- resample.list(eff.ids)
        dat.B <- resample.datMat(boot.list, dat)
        x.cov.bs.B <- resample.datMat(boot.list, x.cov.bs)
        x.cov.td.B <- resample.datMat(boot.list, x.cov.td)
      } else { # original set
        dat.B <- dat
        x.cov.bs.B <- x.cov.bs
        x.cov.td.B <- x.cov.td
      }
      # Fit sequential regression for g-formula:
      result.E.Ybar <- g.formula.outcome.c.nc(times.vec, dat.B, x.cov.bs.B, x.cov.td.B, td.cov.names = td.cov.names, 
                                              no.last.tdv = no.last.tdv, lag = lag,
                                              maxT = ntimes, outcome = outcome,
                                              method = method)
      result.E.Ybar.a1.B[iter,] <- result.E.Ybar[["E.Ybar.a1"]]
      result.E.Ybar.a0.B[iter,] <- result.E.Ybar[["E.Ybar.a0"]]
    }
    
    boot.result <- list(result.E.Ybar.a1.B, result.E.Ybar.a0.B)
    return(boot.result)
    
  }


#####################################################################################  
# 4) CI estimates of mean function of Y with complete compliance and complete
#    non-complianece with bootstrapping
##################################################################################### 

#
bootstrap.ci.c.nc <- 
  function(aspirin, dat.names, cov.bs.names, cov.td.names, 
           outcome.v, exclude.list, reg.method = "RF", 
           B, 
           alpha = 0.05, method = "pivotal", lgd.position = "topleft", ylim.up = 1)
  {
    
    boot.data <- gformula.detinvt.c.nc(aspirin, dat.names, cov.bs.names, cov.td.names, 
                                       outcome.v= outcome.v,
                                       exclude.list=exclude.list, 
                                       reg.method = reg.method,
                                       B = B)
    
    all.T <- table(aspirin[aspirin$last==1, 2])
    all.T.mat <- t(as.matrix(all.T))
    unique.times <- as.numeric(colnames(all.T.mat))
    times.vec <- rev(unique.times) # {T}
    
    # compute bootstrap CI
    ci <- bootstrap.plot.c.nc(times.vec, boot.data[[1]], boot.data[[2]], ylim.up, outcome.v, alpha, lgd.position)
    
    return(ci)
    
  }
