###################################################
# zero inflated model for distribution   wangbinzjcc 2014-4-21
###################################################

###############################################################
#  rm(list=ls())
###############################################################
setwd("F:\\DataW\\DeadStandingTrees")
dir()
##
dea.5 <- read.csv("DeadData-5m-abun.csv")
head(dea.5)
#
XX <- dea.5
for(i in 1:dim(XX)[1]){
  x0 <- XX[i, 'x']
  y0 <- XX[i, 'y']
  logi.0 <- abs(XX$x-x0)<=1 & abs(XX$y-y0)<=1
  logi.0[i] <- FALSE
  XX[i,"mean.neigh"] <- round(mean(XX[logi.0,'Abu']), 2)
}
dea.5 <- XX
Abu.5 <- dea.5$Abu
hist(Abu.5)
yy=Abu.5
###########################################################
topo.5 <- read.csv("topograCal-5m.csv" )
head(topo.5)
t5.p3 <- topo.5[, c('meanelev', 'slope', 'convex')]
t5.p3.cs <- apply(t5.p3, 2, CenterScale, Binary="centerscale")
t5.p3.cs <- as.matrix(t5.p3.cs)
t5.p3.poly3 <- poly(t5.p3.cs, degree=3,raw =T )
t5.p3.poly3 <- as.data.frame(t5.p3.poly3[,c('1.0.0','2.0.0','3.0.0','0.1.0','0.2.0','0.3.0','0.0.1','0.0.2','0.0.3')])
names(t5.p3.poly3) <- paste(rep(names(t5.p3),each=3), rep(1:3,times=3), sep='')
t5.asp <- data.frame(cos.asp = cos(topo.5$aspect*pi/180), 
                     sin.asp = sin(topo.5$aspect*pi/180))
X.t5.cs <- cbind(t5.p3.poly3, t5.asp)
###############################################################
#
###############################################################
require(LaplacesDemon)     
###############################################################
# testing data ~~~~~~~~
# ####################################
set.seed(6)
N=10000
beta.zero <- c(-2,runif(4, -1, 1))
beta.one  <- c(3, runif(11, -1, 1))
X1 <- matrix(runif(N*12),N,12) ; X1[,1] <- 1
prob.zero.logit <- X1[,1:5] %*% beta.zero + rnorm(N,0,0.1) 
prob.zero <- invlogit(prob.zero.logit)
mean(prob.zero)
rr.zero <- rbinom(N, size=1, prob=1-prob.zero)
hist(rr.zero)
#
lamb.log.one <- X1[,1:12] %*% beta.one + rnorm(N,0,0.1)
lamb.one <- exp(lamb.log.one)
mean(lamb.one)
hist(lamb.one)
#
rr.lamb.one <- rr.zero * lamb.one 
yy <- rpois(N,lambda=rr.lamb.one); yy
hist(yy)  
#
#############################################################
Initial.Values <- parm <- c(beta.zero=rep(0,5), beta.one=rep(0,12))
MyData <- Data <- list( parm.names=names(parm), y=yy,X=X1,
                        mon.names=c("Deviance")  )
###### Model for Binomial distribution ~~~~~~~~~~~~~~~~#####
Model <- function(parm, Data){
  # parameter
  beta.zero <- parm[1:5]
  beta.one <- parm[5+ 1:12]
  #
  prob.zero.logit <- Data$X[,1:5] %*% beta.zero + rnorm(N,0,0.1) 
  prob.zero <- invlogit(prob.zero.logit)
  rr.zero <- rbinom(N, size=1, prob=1-prob.zero)
  #
  lamb.log.one <- X1[,1:12] %*% beta.one + rnorm(N,0,0.1)
  lamb.one <- exp(lamb.log.one)
  ### #   priors distribution
  beta.prior <- sum(dnorm(beta.zero,beta.one), 0, 100,log=T)
  prob.zero.logit.prior <- sum(dnorm(prob.zero.logit, 0, 2, log=T))
  lamb.one.prior <- sum(dgamma(lamb.one, 5, log=T))
  #
  LL <- sum(dpois(x=Data$y, lambda=rr.zero * lamb.one , log=T))
  # log-Posterior
  LP <- LL + beta.prior + prob.zero.logit.prior + lamb.one.prior 
  #    
  Modelout <- list(LP=LP, Dev=-2*LL, 
                   Monitor=c(-2*LL),  
                   yhat=rpois(length(Data$y),lambda=rr.zero * lamb.one ) ,
                   parm=parm)
  return(Modelout)  
}
#####################################################################
#
Iterations <- 5000
Status <- 1000
Thinning <- 50

##########################  Adaptive Metropolis  ##########################
 ZIP.Fit2 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                      Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
                      Algorithm="AM", Specs=list(Adaptive=500, Periodicity=10))

######################  Robust Adaptive Metropolis  #######################
 ZIP.Fit18 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                       Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
                       Algorithm="RAM", Specs=list(alpha.star=0.234, Dist="N", gamma=0.66,
                                                   Periodicity=10))

###########################################################################
#
(BF <- ZeroInflatedBinomial.Fit2)
# 
########
###
BF$Summary2

#################
#   write.csv(BF$Summary2, 'ZeroInflatedBinomial.distribution 5 m.csv')
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
dea.5 <- read.csv("DeadData-5m-abun.csv")
head(dea.5)
yy <- Abu.5 <- dea.5$Abu
hist(Abu.5)
mean(Abu.5)
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
#   write.csv(NB$Summary2, 'ZeroInflatedNegativeBinomial.distribution 5m.csv')
############################################################################# 

############################################################################
#  rm(list=ls())
#
require(LaplacesDemon)     
########################################################################
###### Model for logist distribution  of 0ne-Zero data ~~~~~~~~~~~~~~~~#####
Model <- function(parm, Data){
  # parameters
  Size <- parm[1] <- 1 #abs(round(parm[1]))
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
# #####################################################################
# N = 375
# J = 2
# rr <- rbinom(N, size=1, prob=0.3)
# hist(rr) 
# ####################################################
# Initial.Values <- parm <- c(size=10, prob.logit=logit(0.3))
# MyData <- Data <- list(J=J, N=N, 
#                        mon.names=c("LogPosterior", "Prob"),
#                        parm.names=names(parm), y=rr)
# ######################################################################
setwd("F:\\DataW\\DeadStandingTrees")
##
dea.5 <- read.csv("DeadData-5m-abun.csv")
head(dea.5)
hist(dea.5$Abu)
Abu.5.zero <- as.numeric(dea.5$Abu>0)
hist(Abu.5.zero)
mean(Abu.5.zero)
#
N <- length(Abu.5.zero)
J <- 2
################
Initial.Values <- parm <- c(size=1, prob.logit=logit(0.3))
MyData <- Data <- list(J=J, N=N, 
                       mon.names=c("LogPosterior.wb", "Prob"),
                       parm.names=names(parm), y=Abu.5.zero)
###############################################################
Iterations <- 50000 
Status <- 2000
Thinning <- 50

##########################  Adaptive Metropolis  ##########################
binomial.Fit2 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                               Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
                               Algorithm="AM", Specs=list(Adaptive=500, Periodicity=10))

# ######################  Robust Adaptive Metropolis  #######################
# binomial.Fit18 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#                                 Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
#                                 Algorithm="RAM", Specs=list(alpha.star=0.234, Dist="N", gamma=0.66,
#                                                             Periodicity=10))
# 
# ###########################################################################
#
(OZ05 <- binomial.Fit2)
# Fit6
# Fit7
# Fit13
binomial.Fit18
########
###
OZ05$Summary2

#################
#   write.csv(OZ05$Summary2, 'one-zero-binomial.distribution 5 m.csv')
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
#   Model for zero Infalted  of logist and possion distribution   #
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   #
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   #
########################################################################
#  Model for zero Infalted  of logist and possion distribution ~~~~###
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
dea.5 <- read.csv("DeadData-5m-abun.csv")
head(dea.5)
hist(dea.5$Abu)
yy <- dea.5$Abu
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
#      plot(rgamma(100, 25))
#######################################################################
###### ~~~ Model for zero Infalted  of logist and possion distribution ~~~~~~~~~~~~~~~~#####
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
Iterations <- 100000
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

#  write.csv(ZeroPois, "AIC.ZeroPois.distribution.5m.csv")

####################################################

# ######################  Robust Adaptive Metropolis  #######################
# ZeroInflated.possion.Fit18 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#                                 Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
#                                 Algorithm="RAM", Specs=list(alpha.star=0.234, Dist="N", gamma=0.66,
#                                                             Periodicity=10))
# 
# ###########################################################################
# 

