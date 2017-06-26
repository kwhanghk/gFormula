#####################################################################################  
#                   R Code for "G-formula with longitudinal data"                   #                
#                   Kwhangho Kim (kwhanghk@stat.cmu.edu)                            #
#                                                                                   #
#                                                      Final Edit : June/23/2017    #
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
setwd("~/Documents/R_Working/ADA")
source("function_gFormula_class1.R")
source("function_gFormula_class2.R")

# load data
load("aspirin.rda") # Ashley's version
head(aspirin)
summary(aspirin)
dim(aspirin)
aspirin[aspirin$id ==1148, ]



#####################################################################################  
# 1) Point estimate of mean function of Y when we assign a sequence of treatment 
#    \bar{a}_T (deterministic static intervention)
##################################################################################### 


# -----------------------------------------------------------------------------------
# setup: manipulate data structure so that it becomes balanced across id
# intput: aspirin (raw data)
# output = aspirin (elongated), exclude list
# -----------------------------------------------------------------------------------

datastr <- dataStr(aspirin)

aspirin <- datastr$data.long # balanced form
exclude.list <- datastr$exclude.list # any irregularities to be excluded

# -----------------------------------------------------------------------------------
# Fiting G-formula with the appropriate data structure
# input: data info & covariates vectors; outcome; treatment sequence
#        balanced-form data; 
# output = pointwise estimates via g-formula at each time point 
# -----------------------------------------------------------------------------------

# input: 
# id, time, Y, W, other basic info
dat.names = c("id","study_month", "last", # key
              "live_birth", "fetal_loss", "efuwp", "conceived",   # possible outcomes
              "withdraw") # to take censorship process into account, specify this "withdraw" variable

# non time-dependent covariates (baseline) :
cov.bs.names = c("id", # key
                 "site", "rand_date", "eligibility", "age" ,"income","education", 
                 "white","marital", "employed")
# time-dependent covariates :
cov.td.names = c("id", "study_month", # key
                 "compliance", "bleeding", "conceived")
# specify outcome variable of interest:
outcome.v = "live_birth"

# sequence of treatment (deterministic intervention)
a.bar = rep(1,15)


# g-formula estimation
aspirin.a1 <- gformula.detinvt(aspirin, dat.names, cov.bs.names, cov.td.names, 
                               outcome.v = outcome.v,
                               exclude.list =exclude.list, 
                               a.bar=a.bar, reg.method = "RF")



a.bar = rep(0,15)

aspirin.a0 <- gformula.detinvt(aspirin, dat.names, cov.bs.names, cov.td.names, 
                               outcome.v = outcome.v,
                               exclude.list =exclude.list, 
                               a.bar=a.bar, reg.method = "RF")


plot(aspirin.ex[1,], type = "l", lwd = 2.5, col = "orange", xlab = "T", ylab = paste("Probability of", outcome.v))


#
a.bar = c(rep(0,8), rep(1,7))

# a.bar = c(rep(1,8), rep(0,7))
outcome.v = "fetal_loss"
aspirin.a <- gformula.detinvt(aspirin, dat.names, cov.bs.names, cov.td.names, 
                               outcome.v = outcome.v,
                               exclude.list =exclude.list, 
                               a.bar=a.bar, reg.method = "RF")

#####################################################################################  
# 2) CI estimates of the mean function with bootstrapping
#    when \bar{a}_T is assigned (deterministic static intervention)
##################################################################################### 

a.bar = c(rep(0,8), rep(1,7))
n.B = 100 

bootstrap.ci(aspirin, dat.names, cov.bs.names, cov.td.names,
             outcome.v, exclude.list, a.bar, reg.method = "RF", 
             B = n.B, 
             alpha = 0.05, method = "pivotal", lgd.position = "topleft", ylim.up = 1)



#####################################################################################  
# 3) Point estimate of mean function of Y with complete compliance and complete
#    non-complianece
##################################################################################### 


aspirin.lb <- gformula.detinvt.c.nc(aspirin, dat.names, cov.bs.names, cov.td.names, 
                                    outcome.v = "live_birth",
                                    exclude.list =exclude.list, 
                                    reg.method = "RF")


aspirin.fl <- gformula.detinvt.c.nc(aspirin, dat.names, cov.bs.names, cov.td.names, 
                                    outcome.v = "fetal_loss",
                                    exclude.list =exclude.list, 
                                    reg.method = "RF")


aspirin.cr <- gformula.detinvt.c.nc(aspirin, dat.names, cov.bs.names, cov.td.names, 
                                    outcome.v = "conceived",
                                    no.last.tdv = TRUE,
                                    exclude.list =exclude.list, 
                                    reg.method = "RF")


#####################################################################################  
# 4) CI estimates of mean function of Y with complete compliance and complete
#    non-compliance with bootstrapping
##################################################################################### 

outcome.v = "live_birth"
n.B = 100 

bootstrap.ci.c.nc(aspirin, dat.names, cov.bs.names, cov.td.names, 
                  outcome.v, exclude.list, reg.method = "RF", 
                  B = n.B, 
                  alpha = 0.05, method = "pivotal", lgd.position = "topleft", ylim.up = 1)