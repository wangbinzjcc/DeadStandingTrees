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
###########################################################################
# Examples of MCMC Algorithms                                             #
###########################################################################
Iterations=200000
Status=1000
Thinning=50
########################  Hit-And-Run Metropolis  #########################

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

################################################################################

#plot(Fit, BurnIn=50, MyData, PDF=F)
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


#End
############################################################################

Fit2  # AM
Fit6  # CHARM  
Fit7  # CHARM
Fit13   # HARM
Fit18   # RAM

####################
#
par(op)
op <- par(no.readonly = TRUE) 

par(op)
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
##
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
#######################################################################

