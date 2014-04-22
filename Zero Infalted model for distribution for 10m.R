###################################################
# model for distribution   wangbinzjcc 2014-3-17
###################################################
# AIC = - 2 ln(L) + 2 k    
# BIC = - 2 ln(L) + ln(n)*k      
###########################################################################
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     #
#      Model for Binomial distribution                               #
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     #
############################################################################
#  rm(list=ls())
#
require(LaplacesDemon)     
########################################################################
# #####################################################################
# # testing data ~~~~~~~~
# # ####################################
# N = 375
# prob11=0.3
# #
# size=10
# prob22=0.3
# #
# ww <- rbinom(N, size=1, prob=prob11); ww
# hist(ww)
# #
# prob22.ww <- ww*prob22
# yy <- rbinom(N,size=size,prob=prob22.ww); yy
# hist(yy)  
# ###################################
# Initial.Values <- parm <- c(Prob11.logiT=logit(0.2911), 
#                             Size=10, 
#                             Prob22.logiT=logit(0.2911),
#                             PropZeroF.logiT=logit(0.02))
# MyData <- Data <- list( parm.names=names(parm), y=yy,
#     mon.names=c('LL.ww.zero', 'LL.cc.abu', "Prob11", 'Prob22', 'PrZeF')
#                       )
###### Model for Binomial distribution ~~~~~~~~~~~~~~~~#####
Model <- function(parm, Data){
  # parameters
  prob11.logit <- parm[1]; prob11 <- invlogit(prob11.logit)
  size <- parm[2] <- abs(round(parm[2]))
  prob22.logit <- parm[3]; prob22 <- invlogit(prob22.logit)
  prz.logit <- parm[4];    prZeF <- invlogit(prz.logit)
  # make false zeros 
  ww <- rep(NA, length(Data$y))
  ww[Data$y > 0] <- 1
  ww[Data$y==0] <- rbinom(sum(Data$y==0), size=1, prZeF)
  ww.prob22 <- ww * prob22
  # priors distribution
  prob11.logit.priors <- dnorm(prob11.logit, 0, 2, log=T)
  size.priors <- dunif(size, 0, 50, log=T)
  prob22.logit.priors <- dnorm(prob22.logit, 0, 2, log=T)
  prz.logit.prior <- dnorm(prz.logit, 0, 2, log=T)
  #
  # log-likelihood
  LL.ww <- sum(dbinom(x=ww, size=1, prob=prob11, log=T))
  LL.cc <- sum(dbinom(x=Data$y, size=size, prob=ww.prob22, log=T))
  # log-Posterior
  LP <- LL.ww + LL.cc + prob11.logit.priors + 
         size.priors + prob22.logit.priors + prz.logit.prior
  #    
  Modelout <- list(LP=LP, Dev=-2*(LL.ww+LL.cc), 
                   Monitor=c(-2*LL.ww, -2*LL.cc, prob11, prob22, prZeF),  
                   yhat=rbinom(length(Data$y), size=size, prob=ww.prob22),
                   parm=parm)
  return(Modelout)  
}
#####################################################################

######################################################################
setwd("F:\\DataW\\DeadStandingTrees")
##
dea.10 <- read.csv("DeadData-10m-abun.csv")
head(dea.10)
Abu.10 <- dea.10$Abu
hist(Abu.10)
yy=Abu.10
#
################
Initial.Values <- parm <- c(Prob11.logiT=logit(0.2911), 
                           Size=10, Prob22.logiT=logit(0.2911),
                           PropZeroF.logiT=logit(0.02))
#
MyData <- Data <- list( parm.names=names(parm), y=yy,
                     mon.names=c('LL.ww.zero', 'LL.cc.abu', 
                     "Prob11", 'Prob22', 'PrZeF') )
#
Iterations <- 50000 
Status <- 2000
Thinning <- 50

##########################  Adaptive Metropolis  ##########################
ZeroInflatedBinomial.Fit2 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                      Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
                      Algorithm="AM", Specs=list(Adaptive=500, Periodicity=10))

# ######################  Robust Adaptive Metropolis  #######################
#  binomial.Fit18 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#                        Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
#                        Algorithm="RAM", Specs=list(alpha.star=0.234, Dist="N", gamma=0.66,
#                                                    Periodicity=10))
# 
# ###########################################################################
#
(BF <- ZeroInflatedBinomial.Fit2)
# 
########
###
BF$Summary2

#################
#   write.csv(BF$Summary2, 'ZeroInflatedBinomial.distribution 10 m.csv')
############################ 

############################################################################# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   #
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   #
#     Model for ZeroInflatedNegative-Binomial distribution ~~~~ #
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   #
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   #
########################################################################
############################################################################
# rm(list=ls())
#
require(LaplacesDemon)     
#
# # #####################################################################
# # # testing data ~~~~~~~~  
# # ####################################
# N = 375
# prob11=0.3
# #
# size=2
# Mu=3
# #
# ww <- rbinom(N, size=1, prob=prob11); ww
# hist(ww)
# #
# Mu.ww <- ww*Mu
# yy <- rnbinom(N,size=size,mu=Mu.ww); yy
# hist(yy)  
# # ####################################################
# 
# Initial.Values <- parm <- c(Prob11.logiT=logit(0.2911), 
#                             Size=2, 
#                             Mu.loG=log(7),
#                             PropZeroF.logiT=logit(0.02))
# MyData <- Data <- list( parm.names=names(parm), y=yy,
#             mon.names=c('LL.ww.zero', 'LL.cc.abu', "Prob11",
#                         'mu', 'PrZeF'))
########################################################################
###### Model for Negative-Binomial distribution ~~~~~~~~~~~~~~~~#####
Model <- function(parm, Data){
  # parameters
  prob11.logit <- parm[1]; prob11 <- invlogit(prob11.logit)
  size <- parm[2] <- abs(parm[2])
  mu.log <- parm[3]; mu <- exp(mu.log)
  prz.logit <- parm[4];    prZeF <- invlogit(prz.logit)
  # make false zeros 
  ww <- rep(NA, length(Data$y))
  ww[Data$y > 0] <- 1
  ww[Data$y==0] <- rbinom(sum(Data$y==0), size=1, prZeF)
  ww.mu <- ww * mu
  # priors distribution
  prob11.logit.prior <- dnorm(prob11.logit, 0, 2, log=T)
  size.prior <- dunif(size, 0, 20, log=T)
  mu.prior <- dgamma(mu, 5, log=T)
  prz.logit.prior <- dnorm(prz.logit, 0, 2, log=T)
  #
  # log-likelihood
  LL.ww <- sum(dbinom(x=ww, size=1, prob=prob11, log=T))
  LL.cc <- sum(dnbinom(x=Data$y, size=size, mu=ww.mu, log=T))
  # log-Posterior
  LP <- LL.ww + LL.cc + prob11.logit.prior + 
        size.prior + mu.prior + prz.logit.prior
  #   
  Modelout <- list(LP=LP, Dev=-2*(LL.ww+LL.cc), 
         Monitor=c(-2*LL.ww, -2*LL.cc, prob11, mu, prZeF),  
           yhat=rnbinom(length(Data$y), size=size, mu=mu),
                   parm=parm)
  return(Modelout)  
}
#####################################################################

######################################################################
setwd("F:\\DataW\\DeadStandingTrees")
##
dea.10 <- read.csv("DeadData-10m-abun.csv")
head(dea.10)
Abu.10 <- dea.10$Abu
hist(Abu.10)
yy=Abu.10
#

################

Initial.Values <- parm <- c(Prob11.logiT=logit(0.2911), 
                            Size=2, 
                            Mu.loG=log(7),
                            PropZeroF.logiT=logit(0.02))
MyData <- Data <- list( parm.names=names(parm), y=yy,
            mon.names=c('LL.ww.zero', 'LL.cc.abu', "Prob11",
                        'mu', 'PrZeF'))
######################################################################
Iterations <- 50000
Status <- 2000
Thinning <- 50

##########################  Adaptive Metropolis  ##########################
zeroInflatedNegativeBinomial.Fit2 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                               Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
                               Algorithm="AM", Specs=list(Adaptive=500, Periodicity=10))

# ######################  Robust Adaptive Metropolis  #######################
# NegativeBinomial.Fit18 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#                                 Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
#                                 Algorithm="RAM", Specs=list(alpha.star=0.234, Dist="N", gamma=0.66,
#                                                             Periodicity=10))
# 
###########################################################################
# 
(NB <- zeroInflatedNegativeBinomial.Fit2)
###
NB$Summary2
#################
# write.csv(NB$Summary2, 'ZeroInflatedNegativeBinomial.distribution 10m.csv')
############################################################################# 

############################################################################# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   #
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   #
#   Model for zero Infalted for possion distribution   #
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   #
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   #
########################################################################
#  Model for zero Inflated  for possion distribution ~~~~###
######################################################################
#  rm(list=ls()) 
#####################################################################
# testing data ~~~~~~~~
# # ####################################
# N = 375
# prob=0.8
# logit(prob)
# lambda=5
# log(lambda)
# #
# ww <- rbinom(N, size=1, prob=prob); ww
# hist(ww)
# lamb.ww <- ww*lambda
# yy <- rpois(N,lambda=lamb.ww); yy
# rpois(10,lambda=1)
# hist(yy)  
# ###################################
 
# #################################
require(LaplacesDemon)  
#
setwd("F:\\DataW\\DeadStandingTrees")
##
dea.10 <- read.csv("DeadData-10m-abun.csv")
head(dea.10)
hist(dea.10$Abu)
yy <- dea.10$Abu
##################################################
Initial.Values <- parm <- c(lambda.loG=log(1.29), prob.logiT=logit(0.2911),
                                  propZeroF.logiT=logit(0.02))
MyData <- Data <- list(mon.names=c('LL.ww.zero', 'LL.cc.abu', "Lambda", 'Prob', 'PrZeF'),
                       parm.names=names(parm), y=yy)
####################################################
# # 
# plot(rnorm(1000,0,1))
# hist(rnorm(1000,0.5,0.25))
# hist(rgamma(1000,20))
# 
# invlogit(-5:5)
# plot(invlogit(runif(100, -5, 5)))
# plot(invlogit(rnorm(100, 0, 2)))
#  
#   
# 
# plot(rgamma(100, 25))
#######################################################################
###### ~~~ Model for zero Infalted possion distribution ~~~~~~~~~~~~~~~~#####
#######################################################################
 Model <- function(parm, Data){
  # parameters
   lambda <- exp(parm[1])
   prob.logit <- parm[2];   prob <- invlogit(prob.logit)
   prz.logit <- parm[3];    prZeF <- invlogit(prz.logit)
  #
   ww <- rep(NA, length(Data$y))
   ww[Data$y > 0] <- 1
   ww[Data$y==0] <- rbinom(sum(Data$y==0), size=1, prZeF)
  #
   ww.lambda <- ww * lambda
  # priors distribution
   lambda.prior <- dgamma(lambda, 5, log=T)
   prob.logit.prior <- dnorm(prob.logit, 0, 2, log=T)
   prz.logit.prior <- dnorm(prz.logit, 0, 2, log=T)
  #
   LL.ww <- sum(dbinom(x=ww, size=1, prob=prob, log=T))
   LL.cc <- sum(dpois(x=Data$y, lambda=ww.lambda, log=T))
  # log-Posterior
   LP <- LL.ww + LL.cc + lambda.prior + prob.logit.prior + prz.logit.prior
  #    
   Modelout <- list(LP=LP, Dev=-2*(LL.ww+LL.cc), 
                   Monitor=c(-2*LL.ww, -2*LL.cc, lambda, prob, prZeF),  
                   yhat=rpois(n=length(Data$y), lambda=ww.lambda) ,
                   parm=parm)
  return(Modelout)  
   }
#########################################

#########################################
Iterations <- 50000
Status <- 2000
Thinning <- 50 
##########################  Adaptive Metropolis  ##########################
possion.Fit2 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                 Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
                 Algorithm="AM", Specs=list(Adaptive=500, Periodicity=10))
#
#

(pp <- possion.Fit2)

ZeroPois <- pp$Summary2

#  write.csv(ZeroPois, "ZeroInflatedPoission.distribution.10m.csv")

####################################################

# ######################  Robust Adaptive Metropolis  #######################
# ZeroInflated.possion.Fit18 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#                                 Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
#                                 Algorithm="RAM", Specs=list(alpha.star=0.234, Dist="N", gamma=0.66,
#                                                             Periodicity=10))
# 
# ###########################################################################
# 

