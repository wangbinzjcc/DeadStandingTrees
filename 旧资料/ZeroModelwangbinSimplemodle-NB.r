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
###########################################
topo <- read.csv("topograCal-20m.csv")
head(topo)
hist(topo$meanelev)
hist(topo$slope)
hist(log1p(topo$convex))
hist(topo$aspect)
image(x=1:25, y=1:15, z=matrix(topo$meanelev,nrow=25,ncol=15,byrow=T))
elev.c <- (topo$meanelev-mean(topo$meanelev))/(2*sd(topo$meanelev))
elev.c <- CenterScale(topo$meanelev, Binary="centerscale") 
elev2.c <- CenterScale(topo$meanelev^2, Binary="centerscale") 
elev3.c <- CenterScale(topo$meanelev^3, Binary="centerscale") 
#
slope.c <- CenterScale(topo$slope, Binary="centerscale") 
slope2.c <- CenterScale(topo$slope^2, Binary="centerscale") 
slope3.c <- CenterScale(topo$slope^3, Binary="centerscale") 
# #
X.c <-  data.frame(Intercept=1,elev.c, elev2.c, elev3.c, slope.c, slope2.c, slope3.c)
X.c <- apply(X.c,2,as.numeric)
####################################################### 
#  model Relationship between dea.tree abun and elev.c by 20m*20m
#    dea20 = alpha + beta*elev.c + ee
########################################################
 
################
require(LaplacesDemon)     
#
Model <- function(parm, Data){
  # parameters
  beta <- parm[1:7] 
  parm[8] <- kappa <- interval(parm[8], 1E-20, Inf)
  # priors distribution
  beta.priors <- sum(dnormv(beta, 0, 1000, log=T))
  kappa.priors <- dhalfcauchy(kappa, 25, log=T)
  # log-likelihood
  mu <- tcrossprod(Data$X,t(beta))
  LL <- sum(dnbinom(Data$y, size=kappa, mu=mu, log=TRUE))
  # log-Posterior
  LP <- LL + beta.priors + kappa.priors
  #    
  Modelout <- list(LP=LP, Dev=-2*LL,  Monitor=c(LP),  
              yhat=rnbinom(length(mu), size=kappa, mu=mu), parm=parm)
  return(Modelout)  
}
#
# ########################################
#

## true data ~~~~~~~~~~~
#
x.wb <- X.c
y.wb <- dea20 
#
N <- length(y.wb)
J <- ncol(X.c)
#
hist(dea20)
###############################################################
#
###################################
#
Initial.Values <- parm <- c(beta=rep(0,J), kappa=0)
MyData　<- Data <- list(J=J, N=N, 
               mon.names=c("LogPosterior"),
               parm.names=names(parm), X=x.wb, y=y.wb)
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
#
# 100; beta.test
#
Fit2
Fit6
Fit7
Fit13
Fit18
#
Fit <- Fit2
caterpillar.plot(Fit, Parms="beta")
#
points(beta,(5:1),col=2,pch=17)
plot(Fit, BurnIn=50, MyData, PDF=F)
Importance(Fit, Model, MyData, Discrep="Chi-Square")
#
#
a <- lm(y.wb~x.wb)
 summary(a) 
a$coefficients
#
#######
Data0 <- list(x=x.wb, y=y.wb)

Fit <- Fit18
p0 <- as.data.frame(Fit$Summary2)$Mean[1:6]


Model(p0, Data0)$LP
Model(p0, Data0)$Dev
###################################
p0 <- c(a$coefficients, log(0.530))

Model(p0, Data0)$LP
Model(p0, Data0)$Dev
###################################
p0 <- rep(0,6)

Model(p0, Data0)$LP
Model(p0, Data0)$Dev
###############################################################################

