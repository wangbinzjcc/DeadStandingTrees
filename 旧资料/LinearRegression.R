########################
#
rm(list=ls())
require(LaplacesDemon)   
N <- 10000
J <- 5
X <- matrix(1,N,J)
for (j in 2:J) {X[,j] <- rnorm(N,runif(1,-3,3),runif(1,0.1,1))}
beta <- runif(J,-3,3)
e <- rnorm(N,0,0.1)
y <- tcrossprod(X, t(beta)) + e
mon.names <- c("LP", "sigma")
parm.names <- as.parm.names(list(beta=rep(0,J), log.sigma=0))
PGF <- function(Data){c(rnormv(Data$J,0,10), log(rhalfcauchy(1,5)))}
MyData <- list(J=J, X=X, PGF=PGF,mon.names=mon.names, parm.names=parm.names, y=y)

#
#################################################################
#
Model <- function(parm, Data)
{
  ### Parameters
  beta <- parm[1:Data$J]
  sigma <- exp(parm[Data$J+1])
  ### Log(Prior Densities)
  beta.prior <- sum(dnormv(beta, 0, 500, log=TRUE))
  sigma.prior <- dgamma(sigma, 10, log=TRUE)
  ### Log-Likelihood
  mu <- tcrossprod(Data$X, t(beta))
  LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
  ### Log-Posterior
  LP <- LL + beta.prior + sigma.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP, sigma),
                   yhat=rnorm(length(mu), mu, sigma), parm=parm)
  return(Modelout)
}

############################  Initial Values  #############################
Initial.Values <- GIV(Model, MyData, PGF=TRUE)
Initial.Values <- beta
###########################################################################
# Examples of MCMC Algorithms                                             #
###########################################################################

########################  Hit-And-Run Metropolis  #########################
Fit0 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
                     Algorithm="HARM", Specs=NULL)
#Fit0
#beta 
 
#plot(Fit0, BurnIn=50, MyData, PDF=F)
#  Consort(Fit)
#  plot(BMK.Diagnostic(Fit))
#  PosteriorChecks(Fit)
#  caterpillar.plot(Fit, Parms="beta")
#  BurnIn <- Fit$Rec.BurnIn.Thinned
#  plot(Fit, BurnIn, MyData, PDF=FALSE)
#  Pred <- predict(Fit, Model, MyData)
#  summary(Pred, Discrep="Chi-Square")
#  plot(Pred, Style="Covariates", Data=MyData)
#  plot(Pred, Style="Density", Rows=1:9)
#  plot(Pred, Style="ECDF")
#  plot(Pred, Style="Fitted")
#  plot(Pred, Style="Jarque-Bera")
#  plot(Pred, Style="Predictive Quantiles")
#  plot(Pred, Style="Residual Density")
#  plot(Pred, Style="Residuals")
#  Levene.Test(Pred)
#  Importance(Fit, Model, MyData, Discrep="Chi-Square")

##################  Adaptive Hamiltonian Monte Carlo  #####################
Fit1 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=2000, Status=100, Thinning=1,
      Algorithm="AHMC", Specs=list(epsilon=rep(0.02, length(Initial.Values)),
      L=2, Periodicity=10))

##########################  Adaptive Metropolis  ##########################
Fit2 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=2000, Status=100, Thinning=1,
     Algorithm="AM", Specs=list(Adaptive=500, Periodicity=10))

###################  Adaptive Metropolis-within-Gibbs  ####################
Fit3 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=2000, Status=100, Thinning=1,
     Algorithm="AMWG", Specs=list(Periodicity=50))

######################  Adaptive-Mixture Metropolis  ######################
Fit4 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=2000, Status=100, Thinning=1,
     Algorithm="AMM", Specs=list(Adaptive=500, B=NULL, Periodicity=10,
     w=0.05))

###################  Affine-Invariant Ensemble Sampler  ###################
Fit5 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=2000, Status=100, Thinning=1,
     Algorithm="AIES", Specs=list(Nc=2*length(Initial.Values), Z=NULL,
     beta=2, CPUs=1, Packages=NULL, Dyn.lib=NULL))

#################  Componentwise Hit-And-Run Metropolis  ##################
Fit6 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=2000, Status=100, Thinning=1,
     Algorithm="CHARM", Specs=NULL)

###########  Componentwise Hit-And-Run (Adaptive) Metropolis  #############
Fit7 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=2000, Status=100, Thinning=1,
     Algorithm="CHARM", Specs=list(alpha.star=0.44))

#################  Delayed Rejection Adaptive Metropolis  #################
Fit8 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=2000, Status=100, Thinning=1,
     Algorithm="DRAM", Specs=list(Adaptive=500, Periodicity=10))

#####################  Delayed Rejection Metropolis  ######################
Fit9 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=2000, Status=100, Thinning=1,
     Algorithm="DRM", Specs=NULL)

##################  Differential Evolution Markov Chain  ##################
Fit10 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=2000, Status=100, Thinning=1,
     Algorithm="DEMC", Specs=list(Nc=3, Z=NULL, gamma=NULL, w=0.1))

#######################  Hamiltonian Monte Carlo  #########################
Fit11 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=2000, Status=100, Thinning=1,
     Algorithm="HMC", Specs=list(epsilon=rep(0.02, length(Initial.Values)),
     L=2))

#############  Hamiltonian Monte Carlo with Dual-Averaging  ###############
Fit12 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=2000, Status=100, Thinning=1,
     Algorithm="HMCDA", Specs=list(A=500, delta=0.65, epsilon=NULL,
     Lmax=1000, lambda=0.1))

##################  Hit-And-Run (Adaptive) Metropolis  ####################
Fit13 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                     Covar=NULL, Iterations=2000, Status=100, Thinning=1,
                     Algorithm="HARM", Specs=list(alpha.star=0.234))

########################  Independence Metropolis  ########################
### Note: the mu and Covar arguments are populated from a previous Laplace
### Approximation.
Fit14 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=Fit$Covar, Iterations=2000, Status=100, Thinning=1,
     Algorithm="IM",
     Specs=list(mu=Fit$Summary1[1:length(Initial.Values),1]))

#########################  Interchain Adaptation  #########################
Initial.Values <- rbind(Initial.Values, GIV(Model, MyData, PGF=TRUE))
Fit15 <- LaplacesDemon.hpc(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=2000, Status=100, Thinning=1,
     Algorithm="INCA", Specs=list(Adaptive=500, Periodicity=10),
     Chains=2, CPUs=2, Packages=NULL, Dyn.libs=NULL)

#######################  Metropolis-within-Gibbs  #########################
Fit16 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=2000, Status=100, Thinning=1,
     Algorithm="MWG", Specs=NULL)

##########################  No-U-Turn Sampler  ############################
Fit17 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=2000, Status=10, Thinning=1,
     Algorithm="NUTS", Specs=list(A=50, delta=0.6, epsilon=NULL))

######################  Robust Adaptive Metropolis  #######################
Fit18 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=2000, Status=100, Thinning=1,
     Algorithm="RAM", Specs=list(alpha.star=0.234, Dist="N", gamma=0.66,
     Periodicity=10))

###########################  Reversible-Jump  #############################
bin.n <- J-1
bin.p <- 0.2
parm.p <- c(1, rep(1/(J-1),(J-1)), 1)
selectable <- c(0, rep(1,J-1), 0)
Fit19 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=2000, Status=100, Thinning=1,
     Algorithm="RJ", Specs=list(bin.n=bin.n, bin.p=bin.p,
          parm.p=parm.p, selectable=selectable,
          selected=c(0,rep(1,J-1),0)))

########################  Random-Walk Metropolis  #########################
Fit20 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=2000, Status=100, Thinning=1,
     Algorithm="RWM", Specs=NULL)

##############  Sequential Adaptive Metropolis-within-Gibbs  ##############
#NOTE: The SAMWG algorithm is only for state-space models (SSMs)
Fit21 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=2000, Status=100, Thinning=1,
     Algorithm="SAMWG", Specs=list(Dyn=Dyn, Periodicity=50))

##################  Sequential Metropolis-within-Gibbs  ###################
#NOTE: The SMWG algorithm is only for state-space models (SSMs)
Fit22 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=2000, Status=100, Thinning=1,
     Algorithm="SMWG", Specs=list(Dyn=Dyn))

# #############################  Slice Sampler  #############################
# m <- Inf; w <- 1
# Fit23 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
#      Covar=NULL, Iterations=2000, Status=100, Thinning=1,
#      Algorithm="Slice", Specs=list(m=m, w=w))
# 
###################  Tempered Hamiltonian Monte Carlo  ####################
Fit24 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=2000, Status=100, Thinning=1,
     Algorithm="THMC", Specs=list(epsilon=rep(0.05,length(Initial.Values)),
     L=2, Temperature=2))

###############################  t-walk  #################################
Fit25 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=2000, Status=100, Thinning=1,
     Algorithm="twalk", Specs=list(SIV=NULL, n1=4, at=6, aw=1.5))
 
#End
############################################################################
Fit0
Fit1
Fit2
Fit3
Fit4
Fit5
Fit6
Fit7
Fit8
Fit9
Fit10
Fit11
Fit12#!
Fit13
Fit14#!
Fit15#!
Fit16
Fit17#!
Fit18
Fit19
Fit20
Fit21#!
Fit22#!
Fit23#!
Fit24
 
####################
#
op <- par(no.readonly = TRUE) 
Fit = Fit0
beta 
par(op)
caterpillar.plot(Fit, Parms="beta")
#
points(beta,(5:1),col=2,pch=17)
plot(Fit, BurnIn=50, MyData, PDF=F)
#
#
#
#
Fit = Fit1
beta 
par(op)
caterpillar.plot(Fit, Parms="beta")
#
points(beta,(5:1),col=2,pch=17)
plot(Fit, BurnIn=50, MyData, PDF=F)
#
#
##
#
Fit = Fit2
beta 
par(op)
caterpillar.plot(Fit, Parms="beta")
#
points(beta,(5:1),col=2,pch=17)
plot(Fit, BurnIn=50, MyData, PDF=F)
#
#
##
#
Fit = Fit3
beta 
par(op)
caterpillar.plot(Fit, Parms="beta")
#
points(beta,(5:1),col=2,pch=17)
plot(Fit, BurnIn=50, MyData, PDF=F)
#
#
##
#
Fit = Fit4
beta 
par(op)
caterpillar.plot(Fit, Parms="beta")
#
points(beta,(5:1),col=2,pch=17)
plot(Fit, BurnIn=50, MyData, PDF=F)
#
#
##
#
Fit = Fit5
beta 
par(op)
caterpillar.plot(Fit, Parms="beta")
#
points(beta,(5:1),col=2,pch=17)
plot(Fit, BurnIn=50, MyData, PDF=F)
#
#
##
#
Fit = Fit6
beta 
par(op)
caterpillar.plot(Fit, Parms="beta")
#
points(beta,(5:1),col=2,pch=17)
plot(Fit, BurnIn=50, MyData, PDF=F)
#
#
##
#
Fit = Fit7
beta 
par(op)
caterpillar.plot(Fit, Parms="beta")
#
points(beta,(5:1),col=2,pch=17)
plot(Fit, BurnIn=50, MyData, PDF=F)
#
#
##
#
Fit = Fit8
beta 
par(op)
caterpillar.plot(Fit, Parms="beta")
#
points(beta,(5:1),col=2,pch=17)
plot(Fit, BurnIn=50, MyData, PDF=F)
#
#
##
#
Fit = Fit9
beta 
par(op)
caterpillar.plot(Fit, Parms="beta")
#
points(beta,(5:1),col=2,pch=17)
plot(Fit, BurnIn=50, MyData, PDF=F)
#
#
##
# # 
#
Fit = Fit13
beta 
par(op)
caterpillar.plot(Fit, Parms="beta")
#
points(beta,(5:1),col=2,pch=17)
plot(Fit, BurnIn=50, MyData, PDF=F)
#
#
 
#
Fit = Fit16
beta 
par(op)
caterpillar.plot(Fit, Parms="beta")
#
points(beta,(5:1),col=2,pch=17)
plot(Fit, BurnIn=50, MyData, PDF=F)
#
#
 
#
Fit = Fit18 
beta 
par(op)
caterpillar.plot(Fit, Parms="beta")
#
points(beta,(5:1),col=2,pch=17)
plot(Fit, BurnIn=50, MyData, PDF=F)
#
#
##
#
Fit = Fit19
beta 
par(op)
caterpillar.plot(Fit, Parms="beta")
#
points(beta,(5:1),col=2,pch=17)
plot(Fit, BurnIn=50, MyData, PDF=F)
#
#
##
#
Fit = Fit20
beta 
par(op)
caterpillar.plot(Fit, Parms="beta")
#
points(beta,(5:1),col=2,pch=17)
plot(Fit, BurnIn=50, MyData, PDF=F)
#
#
 