############################################################
#  Zero model for deadtree    wangbinzjcc 2014-3-7 8:33:35   
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
elev.c <- (topo$meanelev-mean(topo$meanelev))/sd(topo$meanelev)
hist(elev.c)
####################################################### 
#  model Relationship between dea.tree abun and elev.c by 20m*20m
#    dea.20 = alpha + beta*elev.c + ee
########################################################

################
require(LaplacesDemon)     
#
Model <- function(parm, Data){
  # parameters
  alpha <- parm[1]
  beta <- parm[2] 
  sigma <- exp(parm[3])
  # priors distribution
  alpha.priors <- dnormv(alpha, 0, 100, log=T)
  beta.priors <- dnormv(beta, 0, 100, log=T)
  sigma.priors <- dgamma(sigma, 25, log=T)
  # log-likelihood
  mu <- alpha + beta*Data$x + sigma
  LL <- sum(dnorm(Data$y, mu, sigma, log=T))
  # log-Posterior
  LP <- LL + alpha.priors + beta.priors + sigma.priors
  #    
  Modelout <- list(LP=LP, Dev=-2*LL,  Monitor=c(LP),  
                   yhat= rnorm(length(mu), mu, sigma), parm=parm)
  return(Modelout)  
}
#
# ########################################
#
###  test data ~~~~~~~~~~~~~~
x.test = runif(100, -100, 100)
ee = rnorm(length(x.test),0,0.1)
y.test = 180 - 320*x.test +ee
### ture data ~~~~~~~~~~~~~~~
x.wb <- elev.c 
y.wb <- dea20

#
parm0 <- c(alpha=0, beta=0, sigma=log(1))
PGF <- function(MyData){c(rnormv(MyData$J,0,500), log(rhalfcauchy(1,15)))}
MyData <- list(J=2,
             N=length(x.wb), PGF=PGF,
             mon.names=c("LogPosterior"),
             parm.names=names(parm0), x=x.wb, y=y.wb)
#
Initial.Values <- GIV(Model, MyData,PGF=TRUE)
Iterations <- 100000 
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
#
Fit2
Fit6
Fit7
Fit13
Fit18
Fit <- Fit2
p0 <- as.data.frame(Fit$Summary2)$Mean
#
out <- Fit
#
plot(out, BurnIn=50, Data, PDF=F)
#
a <- lm(y.wb~x.wb)
summary(a)
#
#######
Fit <- Fit6
p0 <- as.data.frame(Fit$Summary2)$Mean

p0 <- c(6.0107,0.4468,log(3.63))

p0 <- c(0,0,0)

parm0 <-  p0
Data0 <- list(x=x.wb, y=y.wb)

Model(parm0, Data0)$LP
Model(parm0, Data0)$Dev
###############################################################################

