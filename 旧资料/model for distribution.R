###################################################
# model for distribution   wangbinzjcc 2014-3-17
###################################################
AIC = - 2 ln(L) + 2 k   
BIC = - 2 ln(L) + ln(n)*k     
###########################################################################
###########################################################################
# rm(list=ls())
require(LaplacesDemon)     
#######################################################################
###### ~~~~~~~~~~~ Model for possion distribution ~~~~~~~~~~~~~~~~#####
#######################################################################
Model <- function(parm, Data){
  # parameters
  lambda <- exp(parm[1])
  parm[2] <- lambda
  # log-likelihood
  LL <- sum(dpois(x=Data$y, lambda=lambda,log=T))
  # log-Posterior
  LP <- LL
  #    
  Modelout <- list(LP=LP, Dev=-2*LL,  Monitor=c(LP, lambda),  
                   yhat= rpois(n=Data$N, lambda = lambda) ,
                   parm=parm)
  return(Modelout)  
}
#####################################
# N = 375
# J=1
# PO <- rpois(n=N, lambda = 400)
# plot(PO)
# hist(PO)
# ####################################################
# Initial.Values <- parm <- c(lambda.log=log(1), lambda.00=log(1))
# MyData <- Data <- list(J=J, N=N, 
#                        mon.names=c("lambda"),
#                        parm.names=names(parm), y=PO)
####################################

setwd("F:\\DataW\\DeadStandingTrees")
##
dea.20 <- read.csv("DeadData-20m-abun-aoto.csv")
head(dea.20)
Abu.20 <- dea.20$Abu
hist(Abu.20)
mean(Abu.20)
#
N <- length(Abu.20)
J <- 2
################
Initial.Values <- parm <- c(lambda.log=log(1), lambda.00= 1 )
MyData <- Data <- list(J=J, N=N, 
                       mon.names=c("LogPosterior.wb", "lambda"),
                       parm.names=names(parm), y=Abu.20)
#
Iterations <- 50000 
Status <- 2000
Thinning <- 50 
##########################  Adaptive Metropolis  ##########################
possion.Fit2 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                      Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
                      Algorithm="AM", Specs=list(Adaptive=500, Periodicity=10))

######################  Robust Adaptive Metropolis  #######################
possion.Fit18 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                       Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
                       Algorithm="RAM", Specs=list(alpha.star=0.234, Dist="N", gamma=0.66,
                                                   Periodicity=10))

###########################################################################
#
(FF2 <- possion.Fit2)
#
possion.Fit18
###
FF2$Summary2

#################
#   write.csv(FF2$Summary2, 'possion.distribution.csv')
############################################################################
#   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     #
#   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     #


#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     #
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     #
############################################################################
#  rm(list=ls())
#
require(LaplacesDemon)     
########################################################################
###### Model for Binomial distribution ~~~~~~~~~~~~~~~~#####
Model <- function(parm, Data){
  # parameters
  Size <- parm[1] <- abs(round(parm[1]))
  Prob <- invlogit(parm[2])
  # priors distribution
  Size.priors <- dunif(Size, 0, 100, log=T)
  Prob.priors <- dunif(Prob, 0, 1, log=T)
  # log-likelihood
  LL <- sum(dbinom(x=Data$y, size=Size, prob=Prob, log=T))
  # log-Posterior
  LP <- LL + Prob.priors + Size.priors
  #    
  Modelout <- list(LP=LP, Dev=-2*LL,  Monitor=c(LP, Prob),  
                   yhat=rbinom(length(Data$y), size=Size, prob=Prob),
                   parm=parm)
  return(Modelout)  
}
#####################################################################
# N = 375
# J = 2
# rr <- rbinom(N, size=50, prob=0.3)
# hist(rr) 
# ####################################################
# Initial.Values <- parm <- c(size=10, prob.logit=logit(0.3))
# MyData <- Data <- list(J=J, N=N, 
#                        mon.names=c("LogPosterior", "Prob"),
#                        parm.names=names(parm), y=rr)
######################################################################
setwd("F:\\DataW\\DeadStandingTrees")
##
dea.20 <- read.csv("DeadData-20m-abun-aoto.csv")
head(dea.20)
Abu.20 <- dea.20$Abu
hist(Abu.20)
mean(Abu.20)
#
N <- length(Abu.20)
J <- 2
################
Initial.Values <- parm <- c(size=10, prob.logit=logit(0.3))
MyData <- Data <- list(J=J, N=N, 
                       mon.names=c("LogPosterior.wb", "Prob"),
                       parm.names=names(parm), y=Abu.20)
#
Iterations <- 50000 
Status <- 2000
Thinning <- 50

##########################  Adaptive Metropolis  ##########################
binomial.Fit2 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                      Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
                      Algorithm="AM", Specs=list(Adaptive=500, Periodicity=10))

######################  Robust Adaptive Metropolis  #######################
binomial.Fit18 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                       Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
                       Algorithm="RAM", Specs=list(alpha.star=0.234, Dist="N", gamma=0.66,
                                                   Periodicity=10))

###########################################################################
#
(BF <- binomial.Fit2)
# Fit6
# Fit7
# Fit13
binomial.Fit18
########
###
BF$Summary2

#################
#   write.csv(BF$Summary2, 'binomial.distribution.csv')
############################
Fit <- Fit2
caterpillar.plot(Fit,Parms= c('BETA','RHO'))
#
plot(Fit, BurnIn=50, MyData, PDF=F)
Importance(Fit, Model, MyData, Discrep="Chi-Square")
#

############################################################################# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   #
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   #


#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     #
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     #
############################################################################
#  rm(list=ls())
#
require(LaplacesDemon)     
########################################################################
###### Model for one-zero Binomial distribution ~~~~~~~~~~~~~~~~#####
Model <- function(parm, Data){
  # parameters
  Size <- parm[1] <- 1
  Prob <- invlogit(parm[2])
  # priors distribution
   Prob.priors <- dunif(Prob, 0, 1, log=T)
  # log-likelihood
  LL <- sum(dbinom(x=Data$y, size=Size, prob=Prob, log=T))
  # log-Posterior
  LP <- LL + Prob.priors  
  #    
  Modelout <- list(LP=LP, Dev=-2*LL,  Monitor=c(LP, Prob),  
                   yhat=rbinom(length(Data$y), size=Size, prob=Prob),
                   parm=parm)
  return(Modelout)  
}
#####################################################################
# N = 375
# J = 2
# rr <- rbinom(N, size=50, prob=0.3)
# hist(rr) 
# ####################################################
# Initial.Values <- parm <- c(size=10, prob.logit=logit(0.3))
# MyData <- Data <- list(J=J, N=N, 
#                        mon.names=c("LogPosterior", "Prob"),
#                        parm.names=names(parm), y=rr)
######################################################################
setwd("F:\\DataW\\DeadStandingTrees")
##
dea.20 <- read.csv("DeadData-20m-abun-aoto.csv")
head(dea.20)
Abu.20 <- as.numeric(dea.20$Abu>0)
hist(Abu.20)
mean(Abu.20)
#
N <- length(Abu.20)
J <- 2
################
Initial.Values <- parm <- c(size=10, prob.logit=logit(0.3))
MyData <- Data <- list(J=J, N=N, 
                       mon.names=c("LogPosterior.wb", "Prob"),
                       parm.names=names(parm), y=Abu.20)
#
Iterations <- 50000 
Status <- 2000
Thinning <- 50

##########################  Adaptive Metropolis  ##########################
one.zero.binomial.Fit2 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                               Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
                               Algorithm="AM", Specs=list(Adaptive=500, Periodicity=10))

######################  Robust Adaptive Metropolis  #######################
binomial.Fit18 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                                Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
                                Algorithm="RAM", Specs=list(alpha.star=0.234, Dist="N", gamma=0.66,
                                                            Periodicity=10))

###########################################################################
#
# (OZ <-one.zero.binomial.Fit2 <- binomial.Fit2)
# Fit6
# Fit7
# Fit13
binomial.Fit18
########
###
OZ$Summary2

#################
#   write.csv(OZ$Summary2, 'binomial.distribution.csv')
############################
Fit <- Fit2
caterpillar.plot(Fit,Parms= c('BETA','RHO'))
#
plot(Fit, BurnIn=50, MyData, PDF=F)
Importance(Fit, Model, MyData, Discrep="Chi-Square")
#

############################################################################# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   #
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   #

#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   #
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   #
########################################################################
############################################################################
# rm(list=ls())
#
require(LaplacesDemon)     
########################################################################
###### Model for Negative-Binomial distribution ~~~~~~~~~~~~~~~~#####
Model <- function(parm, Data){
  # parameters
  Size <- parm[1] <- abs(parm[1])
  Mu <-  parm[2] <- abs(parm[2])
  # priors distribution
  Size.priors <- dunif(Size, 0, 50, log=T)
  Mu.priors <- dunif(Mu, 0, 100, log=T)
  # log-likelihood
  LL <- sum(dnbinom(Data$y, size=Size, mu=Mu, log=T))
  # log-Posterior
  LP <- LL + Size.priors + Mu.priors 
  #    
  Modelout <- list(LP=LP, Dev=-2*LL,  Monitor=c(LP),  
                   yhat=rnbinom(length(Data$y), size=Size, mu=Mu),
                   parm=parm)
  return(Modelout)  
}
#####################################################################
# N = 375
# J = 2
# NB <- rnbinom(N, size=2, mu=6)
# plot(NB)
# hist(NB)
# ####################################################
# Initial.Values <- parm <- c(size=2, mu=6)
# MyData <- Data <- list(J=J, N=N, 
#                        mon.names=c("LogPosterior"),
#                        parm.names=names(parm), y=NB)
######################################################################
setwd("F:\\DataW\\DeadStandingTrees")
##
dea.20 <- read.csv("DeadData-20m-abun-aoto.csv")
head(dea.20)
Abu.20 <- dea.20$Abu
hist(Abu.20)
mean(Abu.20)
#
N <- length(Abu.20)
J <- 2
################
Initial.Values <- parm <- c(size=2, mu=6)
MyData <- Data <- list(J=J, N=N, 
                       mon.names=c("LogPosterior.wb"),
                       parm.names=names(parm), y=Abu.20)
######################################################################
Iterations <- 50000 
Status <- 2000
Thinning <- 50

##########################  Adaptive Metropolis  ##########################
NegativeBinomial.Fit2 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                               Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
                               Algorithm="AM", Specs=list(Adaptive=500, Periodicity=10))

######################  Robust Adaptive Metropolis  #######################
NegativeBinomial.Fit18 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                                Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
                                Algorithm="RAM", Specs=list(alpha.star=0.234, Dist="N", gamma=0.66,
                                                            Periodicity=10))

###########################################################################
#
(NB <- NegativeBinomial.Fit2 )
# Fit6
# Fit7
# Fit13
NegativeBinomial.Fit18 
###
NB$Summary2
#################
#   write.csv(NB$Summary2, 'NegativeBinomial.distribution20m.csv')
############################
Fit <- Fit2
caterpillar.plot(Fit,Parms= c('BETA','RHO'))
#
plot(Fit, BurnIn=50, MyData, PDF=F)
Importance(Fit, Model, MyData, Discrep="Chi-Square")

############################################################################# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   #
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   #

#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   #
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   #
########################################################################
############################################################################
# rm(list=ls())
#
require(LaplacesDemon)     
########################################################################
###### ~~~~~~~~~~ Model for nomial distribution ~~~~~~~~~~~~~~~~#####
########################################################################
Model <- function(parm, Data){
  # parameters
  Mu <-  parm[1] <- abs(parm[1])
  Sd <- parm[2] <- abs(parm[2])
  # priors distribution
  Sd.priors <- dunif(Sd, 0, 100, log=T)
  Mu.priors <- dunif(Mu, 0, 100, log=T)
  # log-likelihood
  LL <- sum(dnorm(Data$y, mean=Mu, sd=Sd, log=T))
  # log-Posterior
  LP <- LL + Sd.priors + Mu.priors 
  #    
  Modelout <- list(LP=LP, Dev=-2*LL,  Monitor=c(LP),  
                   yhat=rnorm(length(Data$y), mean=Mu, sd=Sd),
                   parm=parm)
  return(Modelout)  
}
#####################################################################
# N = 375
# J = 2
# Nor <- rnorm(N, mean=6.2, sd=10)
# plot(Nor)
# hist(Nor) 
 
######################################################################
setwd("F:\\DataW\\DeadStandingTrees")
##
dea.20 <- read.csv("DeadData-20m-abun-aoto.csv")
head(dea.20)
Abu.20 <- dea.20$Abu
hist(Abu.20)
mean(Abu.20)
#
N <- length(Abu.20)
J <- 2
################
Initial.Values <- parm <- c(mean=0, sd=1)
MyData <- Data <- list(J=J, N=N, 
                       mon.names=c("LogPosterior"),
                       parm.names=names(parm), y=Abu.20)
######################################################################
Iterations <- 50000 
Status <- 2000
Thinning <- 50

##########################  Adaptive Metropolis  ##########################
norm.Fit2 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                                       Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
                                       Algorithm="AM", Specs=list(Adaptive=500, Periodicity=10))

######################  Robust Adaptive Metropolis  #######################
norm.Fit18 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                                        Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
                                        Algorithm="RAM", Specs=list(alpha.star=0.234, Dist="N", gamma=0.66,
                                                                    Periodicity=10))

###########################################################################
#
(NF <- norm.Fit2)
#
norm.Fit18
############################
###
NF$Summary2
#################
#   write.csv(NF$Summary2, 'nomial.distribution20m.csv')

#################################################################################
# save.image("F:/DataW/DeadStandingTrees/models for distribution2014-3-17.RData")
#################################################################################
# load("F:/DataW/DeadStandingTrees/models for distribution2014-3-17.RData")
#################################################################################

norm.Fit2

possion.Fit2

binomial.Fit2

NegativeBinomial.Fit2 

###############################################


