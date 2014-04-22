##################################################################
# models for Negative Binomial , Possion ,Linear, Auto 
#                     wangbinzjcc 2014-3-12 8:31:40
##################################################################
rm(list=ls())
#
require(LaplacesDemon)     
########################################################################
###### Model for Negative Binomial ~~~~~~~~~~~~~~~~#####
 Model <- function(parm, Data){
   # parameters
    Beta <- parm[1:Data$J]
    Rho <- parm[Data$J+1]
    parm[Data$J+2] <- Size <- abs(parm[Data$J+2])
   # priors distribution
    beta.priors <- sum(dnormv(Beta, 0, 1000, log=T))
    rho.priors <- dnormv(Rho, 0, 1000, log=T)
    size.priors <- dgamma(Size, 5, log=T)
 # size.priors <- dhalfcauchy(Size, 25, log=T)
  #
    Mu <- exp(tcrossprod(Data$X, t(Beta))+Rho*Data$aoto)
  # log-likelihood
    LL <- sum(dnbinom(x=Data$y, size=Size, mu=Mu, log=T))
  # log-Posterior
    LP <- LL +beta.priors + rho.priors + size.priors
  #    
    Modelout <- list(LP=LP, Dev=-2*LL,  Monitor=c(LP),  
                   yhat=rnbinom(length(Data$y), size=Size, mu=Mu),
                   parm=parm)
  return(Modelout)  
  }
#####################################
#
N <- 375 
J <- 5
X.test <- replicate(J, sample(1:N/100))
X.test <- x.wb <-  cbind(alpha=1, X.test)
beta.test <- runif(1+J, 0,2)
mu.test <- exp(tcrossprod(X.test, t(beta.test)))
y.test <- y.wb <-  rnbinom(n=N, mu = mu.test, size = 5)
hist(mu.test)
hist(log(mu.test))
####################################################
############################################################
#
setwd("F:\\DataW\\DeadStandingTrees")
##
dea.20 <- read.csv("DeadData-20m-abun-aoto.csv")
head(dea.20)
Abu.20 <- dea.20$Abu
aoto.20 <- CenterScale(dea.20$mean.neigh, Binary="centerscale") 
topo.C <- read.csv('topograCal-20m-CenterSCale.csv' )
topo.C$X <- 1
topo.C <- apply(topo.C, 2, as.numeric)
#
head(topo.C)
#
N <- length(Abu.20)
J <- ncol(topo.C)
#################
#
Initial.Values <- parm0 <- c(BETA=rep(0,1+J), RHO=0, SIZE=0)
MyData <- Data <- list(J=J, N=N, 
               mon.names=c("LogPosterior"),
               parm.names=names(parm0), y=Abu.20, X=topo.C, aoto=aoto.20)
#
Iterations <- 50000
Status <- 2000
Thinning <- 50
 
##########################  Adaptive Metropolis  ##########################
Fit2 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                      Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
                      Algorithm="AM", Specs=list(Adaptive=500, Periodicity=10))

#################  Componentwise Hit-And-Run Metropolis  ##################
Fit6 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                      Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
                      Algorithm="CHARM", Specs=NULL)

###########  Componentwise Hit-And-Run (Adaptive) Metropolis  #############
Fit7 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                      Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
                      Algorithm="CHARM", Specs=list(alpha.star=0.44))

##################  Hit-And-Run (Adaptive) Metropolis  ####################
Fit13 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                       Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
                       Algorithm="HARM", Specs=list(alpha.star=0.234))

######################  Robust Adaptive Metropolis  #######################
Fit18 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                       Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
                       Algorithm="RAM", Specs=list(alpha.star=0.234, Dist="N", gamma=0.66,
                                                   Periodicity=10))

###########################################################################
beta.test
#
Fit2
Fit6
Fit7
Fit13
Fit18
######################################
Initial.Values <- parm <-  c(1.75170100, -6.70155200, 13.67185000, -6.89492400,
             -1.62681300,  4.11091800, -2.64457200, -0.08273443,
              0.14726310, -0.01203643,  0.10797780, -0.12376320,
             -0.11647820,  0.43313350, 13.25728000, 3.78717000)
 ################

Fit <- Fit6
caterpillar.plot(Fit,Parms= c('BETA','RHO'))
#
plot(Fit, BurnIn=50, MyData, PDF=F)
Importance(Fit, Model, MyData, Discrep="Chi-Square")
#
plot(rgamma(1000, 5))

rr <- replicate(1000, interval((-100:-1)/100, 1E-20, Inf))
plot(rr)
