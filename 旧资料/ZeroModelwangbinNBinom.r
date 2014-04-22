############################################################
#  Zero model for deadtree    wangbinzjcc 2014-3-10 20:54:12
############################################################
rm(list=ls())
ls()
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
lt20 <- log(t20)
mean(lt20);var(lt20);sd(lt20)
plot(log(t20))
###
image(x=1:25, y=1:15, z=matrix(dea.20$Abu,nrow=25,ncol=15,byrow=T))
###############################################################

####################################################### 
#  model Relationship between dea.tree abun and elev.c by 20m*20m
#    dea20 = alpha + beta*elev.c + ee
########################################################

################
require(LaplacesDemon)     
#
Model <- function(parm, Data){
  # parameters
  mu.abs <- abs(parm[1])
  size.abs <- abs(parm[2])
  # priors distribution
  mu0.priors <- dhalfcauchy(mu.abs, 25, log=T)
  size0.priors <- dhalfcauchy(size.abs, 25, log=T)
  # log-likelihood
  LL <- sum(dnbinom(x=Data$y, size=size.abs, mu=mu.abs, log=T))
  # log-Posterior
  LP <- LL +mu0.priors + size0.priors
  #    
  Modelout <- list(LP=LP, Dev=-2*LL,  Monitor=c(LP,mu.abs, size.abs),  
           yhat=rnbinom(length(Data$y), size=size.abs, mu=mu.abs), parm=parm)
  return(Modelout)  
}
#
# ########################################
#
###  test data ~~~~~~~~~~~~~~
N <- 375 # length(dea20)
J <- 2
x.test <- rnbinom(n=N, mu = 4, size = 10)
hist(x.test)
y.wb <- x.test
#
###################################
N <- 375
J <- 2
y.wb <- dea20
#
parm.wb <- c(mu0=0, size0=0)
PGF <- function(MyData){rhalfcauchy(MyData$J, scale=25)}
MyData <- list(J=J, N=length(y.wb), PGF=PGF,
               mon.names=c("LogPosterior",'mu','size'),
               parm.names=names(parm.wb), y=y.wb)
#
Initial.Values <- GIV(Model, MyData, n=10000, PGF=TRUE)
Iterations <- 10000  
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
#
# 100; beta.test
# 
Fit2
Fit6
Fit7
Fit13
Fit18
#
Fit <- Fit7
caterpillar.plot(Fit,Parms=c('mu','size') )
#
points(beta,(5:1),col=2,pch=17)
plot(Fit, BurnIn=50, MyData, PDF=F)
Importance(Fit, Model, MyData, Discrep="Chi-Square")

#####################################################
 
## Alternative parametrization
x1 <- rnbinom(500, mu = 4, size = 1)

x3 <- rnbinom(500, mu = 4, size = 100)

########################################################################
#
head(dea.20)
dea20 <- dea.20$Abu
rr <- rnbinom(length(dea20), mu= 6.02, size=4.76)
#
hist(dea20, ylim=c(0,100),xlim=c(0,20))
hist(rr, ylim=c(0,100),xlim=c(0,20))
#
#



