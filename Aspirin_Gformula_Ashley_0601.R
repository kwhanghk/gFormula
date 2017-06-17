#####################################################################################  
#                   R Code for "G-formula with EAGeR data"                          #                
#                   Kwhangho Kim (kwhanghk@stat.cmu.edu)                            #
#                                                                                   #
#                                                      Final Edit : Mar/28/2017     #
#####################################################################################  

library(ranger)

# setwd("C:/Users/J Kim/Dropbox/R_Working/2017/ADA")
setwd("~/Documents/R_Working/ADA")
source("function_aspirin.R")
load("aspirin.rda") # Ashley's version
head(aspirin)
summary(aspirin)
dim(aspirin)
aspirin[aspirin$id ==1148, ]


# -----------------------------------------------------------------------------------
# setup: manipulate data structure so that it becomes balanced across id
# output = aspirin(elongated)
# -----------------------------------------------------------------------------------

uniq_ids <- unique(aspirin$id)
length(uniq_ids)  #[1] 1228

n0 <- length(unique(aspirin$id)) 
all.T <- table(aspirin[aspirin$last==1,"study_month"])
barplot(cumsum(all.T))
all.T.mat <- t(as.matrix(all.T))
unique.times <- as.numeric(colnames(all.T.mat))
ntimes <- max(unique.times)
times.vec <- rev(unique.times) # {T}

# elongate each data up to study_month = 15 so that the appended data become
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

dim(aspirin_ext)
1228 * ntimes # inconsistency maybe due to some of missing data (4 subjects)

aspirin <- aspirin_ext[order(aspirin_ext$id),]

id_lb = aspirin[aspirin$live_birth == 1 & aspirin$last == 1,1]
id_fl = aspirin[aspirin$fetal_loss == 1 & aspirin$last == 1,1]
id_wd = aspirin[aspirin$withdraw == 1 & aspirin$last == 1,1]
id_ewop = aspirin[aspirin$efuwp == 1 & aspirin$last == 1,1]

# only subjects whose outcome is either livebirth or fetalloss
# aspirin <- aspirin[aspirin$id %in% union(id_lb, id_fl),]

aspirin["all.y"] <- NA 
aspirin$all.y <- aspirin$fetal_loss + aspirin$efuwp + aspirin$live_birth + aspirin$withdraw


# -----------------------------------------------------------------------------------
# Fiting G-formula with the appropriate data structure
# output = pointwise estimates via g-formula at each time point (or bootstrapped CI)
# -----------------------------------------------------------------------------------

# Specify data structure: 
# subject info & outcome of interest
dat.names = c("id","study_month", "compliance", "conceived", "bleeding",
              "live_birth", "fetal_loss", "all.y", "last", 
              "withdraw") # to take censorship process into account, specify "withdraw" variable

# covariates (baseline (bs) & time-dependent (td)):
cov.bs.names = c("id", "site", "rand_date", "eligibility", "age" ,"income","education", 
                 "white","marital", "employed")
cov.td.names = c("id", "study_month", "compliance", "conceived", "bleeding")

dat = aspirin[,dat.names]; dim(dat)
x.cov.bs = aspirin[aspirin$study_month == 1, cov.bs.names]; dim(x.cov.bs)
x.cov.td = aspirin[,cov.td.names]; dim(x.cov.td)

# identify potential outliers:
## identify {id} with missing values before its final value
exclude.list.1 = NULL
for (m in times.vec) {
  id.m.vec <- aspirin[dat$last == 1 & dat$study_month == m, "id"]
  for (id.i in id.m.vec) {
    if (dim(aspirin[aspirin$id == id.i,])[1] != ntimes) {
      exclude.list.1 <- c(exclude.list.1, id.i)
    }
  }
}
exclude.list.1 <- unique(exclude.list.1)

## identify {id} whose outcome do not fall into one of the four categories until the end of the study
exclude.list.2 = NULL
for (m in times.vec) {
  id.m.vec <- aspirin[dat$last == 1 & dat$all.y != 1, "id"]
  for (id.i in id.m.vec) {
    if (dim(aspirin[aspirin$id == id.i,])[1] != m) {
      exclude.list.2 <- c(exclude.list.2, id.i)
    }
  }
}
exclude.list.2 <- unique(exclude.list.2)

# exclude.list <- c(exclude.list.1, exclude.list.2)
exclude.list <- c(exclude.list.1)
length(exclude.list)

dim(aspirin[!(aspirin$id %in% exclude.list),])
(1215 + length(exclude.list.2)) * ntimes # now the inconsistency disappears

# effective ids
eff.ids <- uniq_ids[!(uniq_ids %in% exclude.list)]

# remove such irregularities
dat = dat[!(dat$id %in% exclude.list),]; dim(dat)
x.cov.bs = x.cov.bs[!(x.cov.bs$id %in% exclude.list),]; dim(x.cov.bs)
x.cov.td = x.cov.td[!(x.cov.td$id %in% exclude.list),]; dim(x.cov.td)

# main loop for g-formula estimation and bootstrapping
boot = FALSE
B = 500 # to get only single time pointwise estimate, set boot = FALSE & B = 1
result.E.Ybar.a0.B <- result.E.Ybar.a1.B <- matrix(nrow = B, ncol = ntimes)

for (iter in 1:B) {
  cat("\n", "bootstrapping step: ", iter, "\n");
  
  if (iter > 1 || boot == TRUE) { # bootstrapped set
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
  result.E.Ybar <- g.formula.outcome(times.vec, dat.B, x.cov.bs.B, x.cov.td.B, NULL,
                                     FALSE, 0, ntimes, "fetal_loss")
  result.E.Ybar.a1.B[iter,] <- result.E.Ybar[["E.Ybar.a1"]]
  result.E.Ybar.a0.B[iter,] <- result.E.Ybar[["E.Ybar.a0"]]
}

# Results:

boot.result <- list(result.E.Ybar.a1.B, result.E.Ybar.a0.B)
save(boot.result, file = paste("bootstrap_",B,"_fl.Rda"))

print(result.E.Ybar.a1.B[1,] - result.E.Ybar.a0.B[1,])
plot(result.E.Ybar.a1.B[1,], type = "l", lwd = 2.5, col = "Blue", ylim = c(0,0.3), xlab = "T", ylab="Probability of fetal loss")
par(new=T)
plot(result.E.Ybar.a0.B[1,], type = "l", lwd = 2.5, col = "Red", ylim = c(0,0.3),  xlab = "T", ylab="Probability of fetal loss")
legend("topleft", c("Compliance","Non-compliance"), 
       lty=c(1,1), 
       lwd=c(2.5,2.5),col=c("blue","red")) 

# compute bootstrap CI
bootstrap.plot(times.vec, result.E.Ybar.a1.B, result.E.Ybar.a0.B, 0.5, "fetal loss", 0.05)
