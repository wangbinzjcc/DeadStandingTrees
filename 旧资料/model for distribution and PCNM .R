###################################################
# model for distribution   wangbinzjcc 2014-3-17
###################################################
rm(list=ls())
#
setwd("F:\\DataW\\DeadStandingTrees")
##
dea.20 <- read.csv("DeadData-20m-abun.csv")
head(dea.20)
dea20 <- dea.20$Abu
t20 <- table(dea20)
t20
mean(t20);var(t20);sd(t20)
plot(t20)
###
image(x=1:25, y=1:15, z=matrix(dea.20$Abu,nrow=25,ncol=15,byrow=T))
###########################################################################
dir()
PCNM.20 <- read.csv("Pcnm.summ 20 m vectors .csv")
PCNM.20 <- PCNM.20[,-1] 
PCNM.posit <- read.csv("Pcnm.summ 20 m Moran_I .csv")
pcnm1 <- PCNM.20[, PCNM.posit$Moran_I.Positive ]
dim(pcnm1)
image(x=1:25, y=1:15, z=matrix(unlist(pcnm1[,30]),nrow=25,ncol=15,byrow=T))

###########################################################################
# rm(list=ls())
require(LaplacesDemon)     
#######################################################################
###### ~~~~~~~~~~~ Model for possion distribution ~~~~~~~~~~~~~~~~#####
#######################################################################
Model <- function(parm, Data){
  # parameters
  alpha= 1
  beta <- parm[1:J]
  #
  beta.priors <- sum(dnorm(beta,0,10,log=T))
  # log-likelihood
  mu <- exp(alpha + Data$X %*% beta)
  LL <- sum(dpois(x=Data$y, lambda=mu, log=T))
  # log-Posterior
  LP <- LL + beta.priors
  #    
  Modelout <- list(LP=LP, Dev=-2*LL,  Monitor=c(LP),  
                   yhat= rpois(n=Data$N, lambda = mu),
                   parm=parm)
  return(Modelout)  
}
#####################################
###  test data ~~~~~~~~~~~~~~
N <- 375
J <- 187
x.test <- matrix(1,N,J)
for(j in 1:J){x.test[,j] <- rnorm(N, runif(1,-1,1), runif(1,0.1,3))}
beta.test <- runif(J, -2, 2)
alpha=  -11
mu.test <- exp((alpha + x.test %*% beta.test)/100) 
summary(mu.test)

y.test <- rpois(n=N, lambda=mu.test)
hist(y.test)
#
#
Initial.Values <- parm <- c(beta.w=rep(1,J))
MyData <- Data <- list(J=J, N=N, X=x.test,
                       mon.names=c("Log-POST"),
                       parm.names=names(parm), y=y.test)
#############################################################
# 
# setwd("F:\\DataW\\DeadStandingTrees")
# ##
# dea.20 <- read.csv("DeadData-20m-abun-aoto.csv")
# head(dea.20)
# Abu.20 <- dea.20$Abu
# hist(Abu.20)
# mean(Abu.20)
# #
# N <- length(Abu.20)
# J <- 2
# ################
# Initial.Values <- parm <- c(lambda.log=log(1), lambda.00= 1 )
# MyData <- Data <- list(J=J, N=N, 
#                        mon.names=c("LogPosterior.wb", "lambda"),
#                        parm.names=names(parm), y=Abu.20)
##############################################################
Iterations <- 50000 
Status <- 2000
Thinning <- 50 
##########################  Adaptive Metropolis  ##########################
possion.PCNM.Fit2 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                      Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
                      Algorithm="AM", Specs=list(Adaptive=500, Periodicity=10))

######################  Robust Adaptive Metropolis  #######################
possion.PCNM.Fit18 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                       Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
                       Algorithm="RAM", Specs=list(alpha.star=0.234, Dist="N", gamma=0.66,
                                                   Periodicity=10))

###########################################################################
#
possion.PCNM.Fit2
#
possion.PCNM.Fit18

Fit <- possion.PCNM.Fit2
caterpillar.plot(Fit,Parms= c('beta'))
#
plot(Fit, BurnIn=50, MyData, PDF=F)
############################################################################
#   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     #
#   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     #


#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     #
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     #
########################################################################
############################################################################
# rm(list=ls())
#
require(LaplacesDemon)     
########################################################################
###### Model for Negative-Binomial distribution ~~~~~~~~~~~~~~~~#####
Model <- function(parm, Data){
  # parameters
  beta <- parm[1:Data$J] 
  Size <- parm[1] <- abs(parm[1])
   # priors distribution
  beta.priors <- sum(dnormv(beta, 0, 500, log=T))
  Size.priors <- dunif(Size, 0, 50, log=T)
  #
  Mu <- exp(Data$X %*% beta)
  # log-likelihood
  LL <- sum(dnbinom(Data$y, size=Size, mu=Mu, log=T))
  # log-Posterior
  LP <- LL + beta.priors + Size.priors
  #    
  Modelout <- list(LP=LP, Dev=-2*LL,  Monitor=c(LP),  
                   yhat=rnbinom(length(Data$y), size=Size, mu=Mu),
                   parm=parm)
  return(Modelout)  
}
#####################################################################
# # test data ~~~~~~~~~~~~~~~~~~~~~~~~~~
# N <- 375
# J <- 187
# x.test <- matrix(1,N,J)
# for(j in 1:J){x.test[,j] <- rnorm(N, runif(1,-1,1), runif(1,0.1,3))}
# beta.test <- runif(J, -2, 2)
# mu.test <- exp((x.test %*% beta.test)/100) 
# summary(mu.test)
# 
# y.test <- rpois(n=N, lambda=mu.test)
# hist(y.test)
# #
# #
# Initial.Values <- parm <- c(beta.w=rep(1,J), Size=1)
# MyData <- Data <- list(J=J, N=N, X=x.test,
#                        mon.names=c("Log-POST"),
#                        parm.names=names(parm), y=y.test)
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
#######################################################
PCNM.20 <- read.csv("Pcnm.summ 20 m vectors .csv")
PCNM.20 <- PCNM.20[,-1] 
PCNM.posit <- read.csv("Pcnm.summ 20 m Moran_I .csv")
pcnm1 <- PCNM.20[, PCNM.posit$Moran_I.Positive ]
J <- dim(pcnm1)[2]
mode(pcnm1)
pcnm1 <- apply(pcnm1, 2 ,as.numeric)
################
Initial.Values <- parm <- c(beta.w=rep(1,J), Size=1)
MyData <- Data <- list(J=J, N=N, X=pcnm1,
                       mon.names=c("Log-Post.wb"),
                       parm.names=names(parm), y=Abu.20)
######################################################################
Iterations <- 50000 * 4
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
(FF <- NegativeBinomial.Fit2)
names(FF)
Initial.Values <- FF$Summary1[,1][1:188]
#
#write.csv(FF$Summary1,'dataForModelofPCNMsummary12014-3-19.csv')
# Fit6
# Fit7
# Fit13
(NN <- NegativeBinomial.Fit18)
Initial.Values <- NN$Summary2[,1][1:188]
# ############################
# Fit <-FF
# Fit <- NN
# caterpillar.plot(Fit,Parms= c('beta'))
# #
# plot(Fit, BurnIn=50, MyData, PDF=F)
# Importance(Fit, Model, MyData, Discrep="Chi-Square")

############################################################################# 
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   #
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   #

#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   #
#  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~   #
########################################################################

#################################################################################
# save.image("F:/DataW/DeadStandingTrees/models for distribution and PCNM 2014-3-19 033033.RData")
#################################################################################
