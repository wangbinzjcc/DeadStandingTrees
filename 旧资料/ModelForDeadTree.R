############################################################
#   model for deadtree    wangbinzjcc 2014-1-24 9:03:43   
############################################################

setwd("F:\\DataW\\DeadStandingTrees")
##
dea.10 <- read.csv("DeadData-10m-abun.csv")
topo.10 <- read.csv("LG.topographic 10m poly4 2013-9-4.csv")
LGtab <- read.csv('LGdata_table.csv')
pi0 <- LGtab$X00枯立木
pi1 <- rowSums(LGtab[-c(1:3)])
head(pi0)
head(LGtab)
lt10 <- log(table(pi0))
plot(lt10)
plot(table(pi0))
#
dir() 
dea.5 <- read.csv("DeadData-5m-abun.csv")
head(dea.5)     
t00 <- table(dea.5$Abu)
plot(t00)
lt05 <- log(t00)
plot(lt05)

#
dea.20 <- read.csv("DeadData-20m-abun.csv")
head(dea.20)
t20 <- table(dea.20$Abu)
lt20 <- log(t20)
mean(lt20);var(lt20)
plot(log(t20))
#######################################################
# de0 <-topo.10[,-c(1,2,3)]
# de0$ymat <- cbind(pi0, pi1 - pi0)
# glm(ymat~., family=binomial(link="logit"), data=de0)
####################################################### 
#  111111
########################################################
require(LaplacesDemon)     
#
Model <- function(parm, Data){
  # parameters
  beta0 <- parm[1] 
  r <- Data$r0
  n <- Data$n0
  # log-likelihood
  logit.pi <- beta0 
  pi <- invlogit(logit.pi)
  LL <- sum(dbinom(x=r, size=n, prob=pi, log=T))
  # priors
  beta0.priors <- dnormv(beta0, 0.0, 1.0E3, log=T)
  #
  LP <- LL + beta0.priors
  #  
  Modelout <- list(LP=LP, Dev=-2*LL,  Monitor=c(LP),  
                   yhat= rbinom(n=length(pi),size=n, prob=pi), parm=parm)
  return(Modelout)  
}
#
#
# ########################################
pi0 <- LGtab$X00枯立木
pi1 <- rowSums(LGtab[-c(1:3)])
#
parm <- c(beta0 =-3.655905)
Data <- list(N= length(pi0),mon.names=c("LogPosterior"), parm.names=names(parm),
             r0=pi0, n0=pi1)
#   Model(parm, Data)
#
#
out <- LaplacesDemon(Model, Data=Data, Initial.Values=parm,
                     Iterations=10000, Status=2000, Thinning=50)
#
out
#
plot(out, BurnIn=50, Data, PDF=F)
#
#######
beta0.mean <- -3.655905
d <- 1
N <- length(pi0)
exp0 = exp(beta0.mean)
pi.mean = exp0/(1+exp0)
BIC = -2*log(prod(dbinom(x=pi0, size=pi1, prob=pi.mean)))+log(N)*d
print(c("BIC ", BIC))

###############################################################################
#  2222
####################################################### 
rm(list=ls())
#
setwd("F:\\DataW\\DeadStandingTrees")
##
dea.10 <- read.csv("DeadData-10m-abun.csv")
topo.10 <- read.csv("LG.topographic 10m poly4 2013-9-4.csv")
LGtab <- read.csv('LGdata_table.csv')
pi0 <- LGtab$X00枯立木
pi1 <- rowSums(LGtab[-c(1:3)])
de0 <-topo.10[, c("meanelev","meanelev2", "convex", "convex2", 
                  "slope","slope2", "asp.cos", "asp.sin")]
de0 <- apply(de0, 2, as.numeric) 
#
require(LaplacesDemon)     
#
Model <- function(parm, Data){
  # parameters
  beta0 <- parm[1] 
  beta1 <- parm[2]
  r <- Data$pi0
  n <- Data$pi1
  X1 <- Data$x1^0.5; X1 <- X1- mean(X1); X1 <- X1/max(abs(X1))
  # log-likelihood
  logit.pi <- beta0 + beta1*X1
  pi <- invlogit(logit.pi)
  LL <- sum(dbinom(x=r, size=n, prob=pi, log=T))
  # priors
  beta0.priors <- dnormv(beta0, 0.0, 1.0E3, log=T)
  beta1.priors <- dnormv(beta1, 0.0, 1.0E3, log=T)
  #
  LP <- LL + beta0.priors + beta1.priors
  #  
  Modelout <- list(LP=LP, Dev=-2*LL,  Monitor=c(LP),  
                   yhat= rbinom(n=length(pi),size=n, prob=pi), parm=parm)
  return(Modelout)  
}
#
#
# ########################################
#
#
#
parm <- c(beta0 =-3.655905, beta1=-0.473)
Data <- list(N= length(pi0),mon.names=c("LogPosterior"), parm.names=names(parm),
             pi0=pi0, pi1=pi1, x1=de0[,"meanelev"])
#
Model(parm, Data)
#
out <- LaplacesDemon(Model, Data=Data, Initial.Values=parm,
                     Iterations=30000, Status=2000, Thinning=50)
#
out
#
plot(out, BurnIn=50, Data, PDF=F)
#
#######
mean00 <- as.data.frame(out$Summary2)$Mean
beta0.mean <- mean00[1]
beta1.mean <- mean00[2]
d <- length(mean00)-2
N <- length(pi0)
exp0 = exp(beta0.mean)
pi.mean = exp0/(1+exp0)
BIC = -2*log(prod(dbinom(x=pi0, size=pi1, prob=pi.mean)))+log(N)*d
print(c("BIC ", BIC))
######################################################
###############################################################################
#  333
####################################################### 
rm(list=ls())
#
setwd("F:\\DataW\\DeadStandingTrees")
##
dea.10 <- read.csv("DeadData-10m-abun.csv")
topo.10 <- read.csv("LG.topographic 10m poly4 2013-9-4.csv")
LGtab <- read.csv('LGdata_table.csv')
pi0 <- LGtab$X00枯立木
pi1 <- rowSums(LGtab[-c(1:3)])
de0 <-topo.10[, c("meanelev","meanelev2", "convex", "convex2", 
                  "slope","slope2", "asp.cos", "asp.sin")]
de0 <- apply(de0, 2, as.numeric) 
#
require(LaplacesDemon)     
#
Model <- function(parm, Data){
  # parameters
  beta0 <- parm[1] 
  beta1 <- parm[2]
  sigma <- exp(parm[3])
  r <- Data$pi0
#  n <- Data$pi1
  X1 <- Data$x1; X1 <- X1- mean(X1); X1 <- X1/max(abs(X1))
  # log-likelihood
  mu <- logit.pi <- beta0 + beta1*X1
#  pi <- invlogit(logit.pi)
#  LL <- sum(dbinom(x=r, size=n, prob=pi, log=T))
  LL <- sum(dnorm(log(r),mu,sigma))
  # priors
  beta0.priors <- dnormv(beta0, 0.0, 1.0E3, log=T)
  beta1.priors <- dnormv(beta1, 0.0, 1.0E3, log=T)
  sigma.priors <- dgamma(sigma, 25, log=T)
  #
  LP <- LL + beta0.priors + beta1.priors + sigma.priors
  #  
  Modelout <- list(LP=LP, Dev=-2*LL,  Monitor=c(LP),  
               #    yhat= rbinom(n=length(pi),size=n, prob=pi), 
                   yhat= rnorm(length(mu),mu,sigma),
                   parm=parm)
  return(Modelout)  
}
#
# #########################################
#
parm <- c(beta0 =0, beta1=0,sigma=log(1))
Data <- list(N= length(pi0),mon.names=c("LogPosterior"), parm.names=names(parm),
             pi0=pi0, pi1=pi1, x1=de0[,"convex"])
#
Model(parm, Data)
#
out <- LaplacesDemon(Model, Data=Data, Initial.Values=parm,
                     Iterations=100000, Status=2000, Thinning=50)
#
out
#
#
plot(out, BurnIn=100, Data, PDF=F)
#
#######
mean00 <- as.data.frame(out$Summary2)$Mean
beta0.mean <- mean00[1]
beta1.mean <- mean00[2]
d <- length(mean00)-2
N <- length(pi0)
exp0 = exp(beta0.mean)
pi.mean = exp0/(1+exp0)
BIC = -2*log(prod(dbinom(x=pi0, size=pi1, prob=pi.mean)))+log(N)*d
print(c("BIC ", BIC))
######################################################