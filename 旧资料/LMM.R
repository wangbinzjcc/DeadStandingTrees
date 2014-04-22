########################################################################
#  model only for Negative Binomial, Auto, and Hierarchical regression 
########################################################################
#  chapter 12
####################
rm(list=ls())
n.groups <- 3
n.sample <-20
n <- n.groups * n.sample
pop <- gl(n=n.groups, k=n.sample)
#
original.length <- runif(n, 45, 70)
mu <- mean(original.length)
sd <- sd(original.length)
cat('Mean and sd used to mormalise.original length:', mu, sd, "\n\n")
length <- (original.length-mu)/sd
hist(length, col='grey')
#
#  Xmat <- model.matrix(~pop*length-1-length)
#  head(Xmat)
#
intercept.mean <- 230
intercept.sd <- 20
slope.mean <- 60
slope.sd <- 30
#
rho=0.8
intercept.slope.covariance  <- rho*intercept.sd*slope.sd
#
mu.vector <- c(intercept.mean, slope.mean)
var.cova.matrix <- matrix(c(intercept.sd^2, intercept.slope.covariance, 
                            intercept.slope.covariance, slope.sd^2),2,2)
# # # ~~~~~~~~~~~~~~~~~~~~~~~~~~
require(LearnBayes)
effects <- rmnorm(n=n.groups, mean=mu.vector, varcov=var.cova.matrix)
#
f <- dmnorm(effects,  mu.vector, var.cova.matrix)
#
#effects ~~~~~~~~~~~~~~~~~~~~~~~
apply(effects, 2, mean)
vv <- c(var(effects));vv
rho <- {vv[2]*vv[3] / (vv[1]*vv[4])} ^0.5 ;rho
rho^2
cor(effects[,1], effects[,2])
#
intercept.effects <- effects[,1]
slope.effects <- effects[,2]
#
Alpha <- rep(intercept.effects, times=table(pop))
Beta <- rep(slope.effects, times=table(pop))
# 
lin.pred <- length*Beta +Alpha
# lin.pred <- Xmat %*% all.effects
eps <- rnorm(n=n, mean=0, sd=30)
mass <- lin.pred +eps
#
hist(mass, col='grey')
#

library(lme4)
lme.fit1 <- lmer(mass~length + (length|pop))
lme.fit1
#
summary(lme.fit1)
#
#
library('lattice')
xyplot(mass~length|pop)
#
#
####################################################################
require(LaplacesDemon)
Initial.Values <- parm <- c(alpha=rep(0,n.groups), beta=rep(0,n.groups),
                    mu.alpha=0, mu.beta=0, sigma.alpha=0.1, sigma.beta=0.1, sigma=0,
                    rho=0.1)
length(parm)

Initial.Values  <- c(229, 246, 254, 25, 113, 107, 241, 89, 55, 128, 29, 0.3)

MyData <- Data <- list(J=n.groups, X=length, Pop=pop, mon.names="LP.wb", 
                       parm.names=names(parm), y=mass, DMnorm=LearnBayes::dmnorm)
#
Model <- function(parm, Data)
{
  ### Parameters
  alpha <- parm[1:Data$J]
  beta <- parm[(Data$J+1):(2*Data$J)] 
  mu.alpha <- parm[(2*Data$J)+1]
  mu.beta <- parm[(2*Data$J)+2]
  parm[(2*Data$J)+3] <- sigma.alpha <- abs(parm[(2*Data$J)+3])
  parm[(2*Data$J)+4] <- sigma.beta <- abs(parm[(2*Data$J)+4])
  parm[(2*Data$J)+5] <- sigma <- abs(parm[(2*Data$J)+5])
  parm[2*Data$J+6] <- rho <- interval(parm[2*Data$J+6],-1,1)
#
  covari <- rho* sigma.alpha * sigma.beta
  varcov <- matrix(c(sigma.alpha^2, covari, covari, sigma.beta^2),2,2)
  ### Log(Hyperprior Densities)
  alpha.beta.varcov.prior <- sum(Data$DMnorm(cbind(alpha, beta),  
                            mean=c(mu.alpha, mu.beta), varcov=varcov,log=T))
  sigma.alpha.prior <- dgamma(sigma.alpha, 25, log=T)
  sigma.beta.prior <- dgamma(sigma.beta, 25, log=T)
  rho.prior <- dunif(rho, -1, 1, log=T)
#
  sigma.prior <- dgamma(sigma, 25, log=T)
  ### Log-Likelihood
  Alpha <- rep(alpha, times=table(Data$Pop))
  Beta <- rep(beta, times=table(Data$Pop))
# 
  mu <- Alpha + Data$X*Beta
  LL <- sum(dnorm(Data$y, mu, sigma, log=T))
  ### Log-Posterior
  LP <- LL + alpha.beta.varcov.prior + rho.prior
       + sigma.alpha.prior + sigma.beta.prior + sigma.prior 
  #
Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP),
                   yhat=rnorm(length(mu), mu, sigma), parm=parm)
  return(Modelout)
}
###
Iterations <- 100000 
Status <- 500
Thinning <- 50
###
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

############################################################################################
#
Fit2
Fit6
Fit7
Fit13
Fit18
#
############################################################################################


 



###########################################################################


Fit <- Fit2
caterpillar.plot(Fit,c('alpha','beta','rho') )
#
plot(Fit, BurnIn=50, MyData, PDF=F)
#
#
#
#
#
#
#
#
#
###########################################################################
#                                                                         #
###########################################################################

#

#

###########################################################################
head(sleepstudy)
#
un0 <- unique(sleepstudy$Subject)
with(sleepstudy, plot(Days, Reaction, type='n'))
for(ii in 1:length(un0)){
  with(subset(sleepstudy, Subject==un0[ii]), points(Days, Reaction,lwd=2, cex=1.5,col=1+ii))
}
######
(lm0 <- lm(Reaction~Days, sleepstudy))
abline(lm0, lwd=2)

##############################################################
pop <- gl(n=56,k=10)
pop
#
require(lme4)
## linear mixed models - reference values from older code
(fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))
summary(fm1)# (with its own print method)
(fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy))
(fm3 <- lmer(Reaction ~ Days + (1|Subject) , sleepstudy))
(fm4 <- lmer(Reaction ~ Days + (0+Days|Subject) , sleepstudy))
anova(fm1, fm2)
sm2 <- summary(fm2)
print(fm2, digits=7, ranef.comp="Var") # the print.merMod()         method
print(sm2, digits=3, corr=FALSE)       # the print.summary.merMod() method
#
(vv <- vcov.merMod(fm2, corr=TRUE))
as(vv, "corMatrix")# extracts the ("hidden") 'correlation' entry in @factors
# 
##########################################################################
head(sleepstudy,20)
un0 <- length(unique(sleepstudy$Subject))

##
require(LaplacesDemon)    
#
data(demonsnacks)
y <- log(demonsnacks$Calories)
X <- cbind(1, as.matrix(demonsnacks[,c(7,8,10)]))
J <- ncol(X)
for (j in 2:J) {X[,j] <- CenterScale(X[,j])}
mon.names <- c("LP","delta","sigma","tau")
parm.names <- as.parm.names(list(beta=rep(0,J), gamma=0, log.delta=0,
                                 log.sigma=0, log.tau=0))
MyData <- list(J=J, X=X, mon.names=mon.names, parm.names=parm.names, y=y)
#
Model <- function(parm, Data)
{
  ### Hyperparameters
  gamma <- parm[Data$J+1]
  delta <- exp(parm[Data$J+2])
  tau <- exp(parm[Data$J+4])
  ### Parameters
  beta <- parm[1:Data$J]
  sigma <- exp(parm[Data$J+3])
  ### Log(Hyperprior Densities)
  gamma.prior <- dnorm(gamma, 0, 1000, log=TRUE)
  delta.prior <- dhalfcauchy(delta, 25, log=TRUE)
  tau.prior <- dhalfcauchy(tau, 25, log=TRUE)
  ### Log(Prior Densities)
  beta.prior <- sum(dnormv(beta, gamma, delta, log=TRUE))
  sigma.prior <- dhalfcauchy(sigma, tau, log=TRUE)
  ### Log-Likelihood
  mu <- tcrossprod(Data$X, t(beta))
  LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
  ### Log-Posterior
  LP <- LL + beta.prior + gamma.prior + delta.prior + sigma.prior +
    tau.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,delta,sigma,tau),
                   yhat=rnorm(length(mu), mu, sigma), parm=parm)
  return(Modelout)
}
#
Initial.Values <- c(rep(0,J), 0, rep(1,3))
###
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
############################################################################################


Fit2
Fit6
Fit7
Fit13
Fit18
###########################################################################





###########################################################################













