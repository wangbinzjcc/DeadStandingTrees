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
dea.5 <- read.csv('DeadData-5m-abun-mrt2014-4-23.csv')
head(dea.5)
#########################################################
# XX <- dea.5
# for(i in 1:dim(XX)[1]){
#   x0 <- XX[i, 'x']
#   y0 <- XX[i, 'y']
#   logi.0 <- abs(XX$x-x0)<=1 & abs(XX$y-y0)<=1
#   logi.0[i] <- FALSE
#   XX[i,"mean.neigh"] <- round(mean(XX[logi.0,'Abu']), 2)
# }
# dea.5 <- XX
#head(dea.5)
#write.csv(dea.5, 'DeadData-5m-abun-mrt2014-4-23.csv')
###########################################################
# topo.5 <- read.csv("topograCal-5m.csv" )
# head(topo.5)
# t5.p3 <- topo.5[, c('meanelev', 'slope', 'convex')]
# t5.p3.cs <- apply(t5.p3, 2, CenterScale, Binary="centerscale")
# t5.p3.cs <- as.matrix(t5.p3.cs)
# t5.p3.poly3 <- poly(t5.p3.cs, degree=3,raw =T )
# t5.p3.poly3 <- as.data.frame(t5.p3.poly3[,c('1.0.0','2.0.0','3.0.0','0.1.0','0.2.0','0.3.0','0.0.1','0.0.2','0.0.3')])
# names(t5.p3.poly3) <- paste(rep(names(t5.p3),each=3), rep(1:3,times=3), sep='')
# t5.asp <- data.frame(cos.asp = cos(topo.5$aspect*pi/180), 
#                      sin.asp = sin(topo.5$aspect*pi/180))
# X.t5.cs <- cbind(t5.p3.poly3, t5.asp)
# head(X.t5.cs)
# write.csv(X.t5.cs, "topograCal-5m-poly3-2014-04-23.csv")
###############################################################

topo.5 <- read.csv("topograCal-5m-poly3-2014-04-23.csv")

###############################################################
require(LaplacesDemon)     
###############################################################
# testing data ~~~~~~~~
# # ####################################
# set.seed(1)
# N=100
# N.tab <- 8
# n.beta <- 3
# X1 <- matrix(runif(N*n.beta), N, n.beta) ; X1[,1] <- 1
# head(X1)
# Tab.num <- sample(1:8, 100,replace=T)
# Tab.lett <- letters[Tab]
# table(Tab.let)
# #
# beta.one.data <- replicate(8, c(runif(1), runif(2, -1, 1)) )
# beta.zero <- c(-2,runif(2, -2, 2))
# #
# beta.one.matrix <- apply(beta.one.data,1,function(xx){xx[Tab.num]})
# lamb.log.one <- rowSums(X1 * beta.one.matrix)
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
# PGF <- function(Data) return(c(rnorm(33,0,1)))
# MyData <- Data <- list( parm.names=names(c(beta.one=rep(0, 3*8), 
#               beta.zero=rep(0,3),mean.beta=rep(0,3),sd.beta.log=rep(0,3))),
#                         y=yy,X=X1,Tab=Tab.num,
#                         PGF=PGF,mon.names=c("sd.beta1",'sd.beta2','sd.beta3')  )
# parm <- PGF(Data)
####################################################################

dat.5 <- dea.5
head(dat.5)
Xdat0 <- topo.5
Xdat <- cbind(alpha=1,Xdat0)
head(Xdat)
yy=dat.5$Abu
Auto=dat.5$mean.neigh
Tab=dat.5$Tab
xx0=apply(Xdat,2,as.numeric)
head(xx0)
#
n.tab <- length(unique(Tab))
tab.num <- match(Tab, unique(Tab))
env.nam <- c('alpha','meanelev1', 'slope1', 'convex1')
n.env <- length(env.nam)
xx <- xx0[,env.nam]
# ##########################################################
# parm <- parm.00 <- c(beta.auto=rnorm(1),beta.one=rnorm(n.tab*n.env),beta.zero=rnorm(n.env))
# beta.auto <- parm[1]
# beta.one <- parm[1 + 1:n.tab*n.env]
# beta.one.matrix<- matrix(beta.one,nrow=n.tab,ncol=n.env)
# beta.one.data <- apply(beta.one.matrix, 2, function(xx){xx[tab.num]})
# lamb.log.one <- rowSums(beta.one.data*xx) + beta.auto * Auto
# lamb.one <- exp(lamb.log.one)
# #
# beta.zero <- parm[(1+n.tab*n.env)+ 1:n.env]
# prob.zero.logit <- xx %*% beta.zero
# prob.ze <- invlogit(prob.zero.logit)
# y.lam <- ifelse(test= prob.ze>0.5, yes=0, no=lamb.one)
# yy <- rpois(N,lambda=y.lam)
# z <- ifelse(test= y>0, yes=1, no=0)
# #
##################################################################
#

#
 parm.00 <- c(beta.auto=rnorm(1),beta.one=rnorm(n.tab*n.env),beta.zero=rnorm(n.env),
             mean.beta.one=rnorm(n.env),sdlog.beta.one=rnorm(n.env))
 PGF <- function(Data) return(rnorm(length(parm.00)))
 MyData <- Data <- list( parm.names=names(parm.00),
                        y=yy,X=xx,
                        Auto=Auto,tab.num=tab.num,
                        PGF=PGF,mon.names=paste("sd.beta",1:n.env)  )
 parm <- PGF(Data)
#
###### Model for Binomial distribution ~~~~~~~~~~~~~~~~#####
Model <- function(parm, Data){
  Ntab <- length(unique(Data$tab.num))
  Nenv <- ncol(Data$X)
  #
  beta.auto <- parm[1]
  beta.one <- parm[1 + 1:(Ntab*Nenv)]
  beta.zero <- parm[(1+Ntab*Nenv)+ 1:Nenv]
  mean.beta.one <- parm[(length(parm.00)-2*Nenv)+1:Nenv]
  sd.beta.one <- exp(parm[(length(parm.00)-Nenv)+ 1:Nenv])
  ##
  beta.one.matrix<- matrix(beta.one,nrow=Ntab,ncol=Nenv)
  beta.one.data <- apply(beta.one.matrix, 2, function(xx){xx[Data$tab.num]})
  lamb.log.one <- rowSums(beta.one.data* Data$X) + beta.auto * Data$Auto
  lamb.one <- exp(lamb.log.one)
  #
  beta.zero <- parm[(1+Ntab*Nenv)+ 1:Nenv]
  prob.zero.logit <- Data$X %*% beta.zero
  prob.zero <- invlogit(prob.zero.logit)
  y.lamb <- ifelse(test= prob.zero>0.5, yes=0, no=lamb.one)
  ### #   priors distribution
  beta.prior <- sum(dnorm(beta.one,
                    mean=rep(mean.beta.one, each=Ntab),
                    sd=rep(sd.beta.one, each=Ntab),log=T))
  mean.beta.prior <- sum(dnorm(mean.beta.one, 0, 50,log=T))
  sd.beta.prior <- sum(dgamma(sd.beta.one,5,log=T))
  lamb.one.prior <- sum(dgamma(lamb.one, 3, log=T))
  prob.zero.logit.prior <- sum(dnorm(prob.zero.logit, 0, 2, log=T))
  #
  LL <- sum(dbern(x=as.numeric(Data$y==0), prob=prob.zero , log=T),
            dpois(x=Data$y, lambda=y.lamb, log=T))
  # log-Posterior
  LP <- LL + beta.prior + mean.beta.prior +
        sd.beta.prior + lamb.one.prior + prob.zero.logit.prior 
  #    
  Modelout <- list(LP=LP, Dev=-2*LL, 
                   Monitor=c(sd.beta.one),  
                   yhat=rpois(length(Data$y),lambda=y.lamb),
                   parm=parm)
  return(Modelout)  
}
#####################################################################
#
Initial.Values <- parm <- GIV(Model, MyData, n=100000, PGF=TRUE)
#
Iterations <-  200000
Status <- 2000
Thinning <- 10
                                
##################################################
for(i in 1:10){
  print(rep(i,10))
Initial.Values <- as.initial.values(Fit6)
#
Fit6  <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                      Covar=NULL, Iterations, Status, Thinning,
                      Algorithm="CHARM", Specs=NULL)
dput(summary(Fit6),"FIT6.Summary.2014.0423.csv")
}
####################################################

Juxt <- Juxtapose(list(Fit2,Fit6,Fit7)); Juxt
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