#################################################################
# zero inflated model for distribution   wangbinzjcc 2014-4-21
#################################################################
#
#  rm(list=ls())
############################# 
setwd("F:\\DataW\\DeadStandingTrees")
dir()
##
#dea.5 <- read.csv("DeadData-5m-abun.csv")
dea.5 <- read.csv('DeadData-5m-abun-mrt2014-4-22.csv')
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
head(X.t5.cs)
###############################################################
#
###############################################################
require(LaplacesDemon)     
# ###############################################################
# # testing data ~~~~~~~~
# # ####################################
# set.seed(1)
# N=100
# X1 <- matrix(runif(N*12),N,12) ; X1[,1] <- 1
# 
# beta.zero <- c(-2,runif(4, -2, 2))
# beta.one  <- c(1, runif(11, -1, 1))
# #
# lamb.log.one <- X1[,1:12] %*% beta.one
# lamb.one <- exp(lamb.log.one)
# mean(lamb.one)
# hist(lamb.one)
# #
# prob.zero.logit <- X1[,1:5] %*% beta.zero
# prob.ze <- invlogit(prob.zero.logit)
# y <- ifelse(test= prob.ze>0.5, yes=0, no=lamb.one)
# yy <- rpois(N,lambda=y)
# z <- ifelse(test= y>0, yes=1, no=0)
# #
# hist(yy)  
# #
# #############################################################
# PGF <- function(Data) return(c(rnorm(17,0,1)))
# MyData <- Data <- list( parm.names=names(c(beta.zero=rep(0,5), beta.one=rep(0,12))),
#                         y=yy,X=X1,
#                         PGF=PGF,mon.names=c("Deviance")  )
#####################################################################
yy=dea.5$Abu
Auto=dea.5$mean.neigh
xx=apply(X.t5.cs,2,as.numeric)
head(xx)
mode(xx)
#
 PGF <- function(Data) return(c(rnorm(1+12,0,1), rnorm(4,0,1)))
#
 MyData <- Data <- list( parm.names=names(c(auto=0, beta.one=rep(0,12), beta.zero=rep(0,4))),
                         y=yy,X=xx,
                         PGF=PGF,mon.names=c("loglikelihood")  )
parm <- PGF(Data)
#
###### Model for Binomial distribution ~~~~~~~~~~~~~~~~#####
Model <- function(parm, Data){
  # parameter
  auto <- parm[1]
  beta.one <- parm[2:13]
  beta.zero <- parm[13+ 1:4]
  #
  lamb.log.one <- cbind(alpha=1,Data$X)  %*% beta.one
  lamb.one <- exp(lamb.log.one)
  #
  nam0 <- c('meanelev1','slope1','convex1')
  prob.zero.logit <- cbind(alpha=1,Data$X[,nam0]) %*% beta.zero 
  prob.zero  <- invlogit(prob.zero.logit)
  #
  y.lamb <- ifelse(test= prob.zero>0.5, yes=0, no=lamb.one)
  #
  ### #   priors distribution
  beta.prior <- sum(dnorm(c(auto,beta.one,beta.zero), 0, 5,log=T))
  lamb.one.prior <- sum(dgamma(lamb.one, 3, log=T))
  prob.zero.logit.prior <- sum(dnorm(prob.zero.logit, 0, 2, log=T))
  #
  LL <- sum(dbern(x=as.numeric(Data$y==0), prob=prob.zero , log=T),
            dpois(x=Data$y, lambda=y.lamb, log=T))
  # log-Posterior
  LP <- LL + beta.prior + lamb.one.prior + prob.zero.logit.prior 
  #    
  Modelout <- list(LP=LP, Dev=-2*LL, 
                   Monitor=c(-2*LL),  
                   yhat=rpois(length(Data$y),lambda=y.lamb),
                   parm=parm)
  return(Modelout)  
}
#####################################################################
#
Initial.Values <- parm <- GIV(Model, MyData, n=100000, PGF=TRUE)
#
Initial.Values <- as.initial.values(Fit6)
#
Iterations <- 200000
Status <- 2000
Thinning <- 10
# ##########################  Adaptive Metropolis  ##########################
#  Fit2 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#                       Covar=Fit2$Covar,, Iterations=Iterations, Status=Status, Thinning=Thinning,
#                       Algorithm="AMWG", Specs=list(Periodicity=128))

#################  Componentwise Hit-And-Run Metropolis  ##################
 Fit6 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                       Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
                       Algorithm="CHARM", Specs=NULL)

dput(summary(Fit6),'Fit6')
#dget('Fit6')
#  
# ##########  Componentwise Hit-And-Run (Adaptive) Metropolis  #############
#  Fit7 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#                        Covar=NULL, Iterations=Iterations, Status=Status, Thinning=Thinning,
#                       Algorithm="CHARM", Specs=list(alpha.star=0.44))
#  
# #####################  Robust Adaptive Metropolis  #######################
#  Fit18 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#                        Covar=Fit18$Covar, Iterations=Iterations, Status=Status, Thinning=Thinning,
#                        Algorithm="RAM", Specs=list(20))
# 
# ##########################################################################
#
Juxt <- Juxtapose(list(Fit2,Fit6,Fit7,Fit18 )); Juxt
plot(Juxt, Style="ISM")
#
Consort(Fit2)
Consort(Fit6)
Consort(Fit7)
Consort(Fit18)
Fit  <- Fit6
Consort(Fit)
plot(BMK.Diagnostic(Fit, batches=10))
PosteriorChecks(Fit)
caterpillar.plot(Fit, Parms="beta")
BurnIn <- Fit$Rec.BurnIn.Thinned
plot(Fit, BurnIn, MyData, PDF=FALSE)
Pred <- predict(Fit, Model, MyData)
#
summary(Pred, Discrep="Chi-Square")
summary(Pred, Discrep="Kurtosis")
plot(Pred, Style="Covariates", Data=MyData)
plot(Pred, Style="Density", Rows=1:9)
plot(Pred, Style="ECDF")
plot(Pred, Style="Fitted")
plot(Pred, Style="Jarque-Bera")
plot(Pred, Style="Predictive Quantiles")
plot(Pred, Style="Residual Density")
plot(Pred, Style="Residuals")
Levene.Test(Pred)
Importance(Fit, Model, MyData, Discrep="Chi-Square")
############################################################################


#