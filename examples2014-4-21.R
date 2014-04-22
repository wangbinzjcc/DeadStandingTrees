#############################################################
# The accompanying Examples vignette is a compendium of examples.
####################  Load the LaplacesDemon Library  #####################
##
Bi <- rbinom(n=400,size=10,prob=0.3)
Pi <- rpois(n=400, lambda = 3)
NB1 <- rnbinom(n=400,size=1,mu=3)
NB2 <- rnbinom(n=400,size=3,mu=3)
#
plot(Bi)
hist(Bi)
plot(Pi)
hist(Pi)
plot(NB1)
hist(NB1)
plot(NB2)
hist(NB2)
#
image(matrix(Bi,20,20))
contour(matrix(Bi,20,20),add=T)

image(matrix(Pi,20,20))
contour(matrix(Pi,20,20),add=T)

image(matrix(NB1,20,20))
contour(matrix(NB1,20,20),add=T)

image(matrix(NB2,20,20))
contour(matrix(NB2,20,20),add=T)
######################################################

library(LaplacesDemon)

##############################  Demon Data  ###############################
data(demonsnacks)
y <- log(demonsnacks$Calories)
X <- cbind(1, as.matrix(log(demonsnacks[,c(1,4,10)]+1)))
J <- ncol(X)
for (j in 2:J) {X[,j] <- CenterScale(X[,j])}
mon.names <- c("LP","sigma")
parm.names <- as.parm.names(list(beta=rep(0,J), log.sigma=0))
PGF <- function(Data) return(c(rnormv(Data$J,0,10),
                               log(rhalfcauchy(1,5))))
MyData <- list(J=J, PGF=PGF, X=X, mon.names=mon.names,
               parm.names=parm.names, y=y)

##########################  Model Specification  ##########################
Model <- function(parm, Data)
{
  ### Parameters
  beta <- parm[1:Data$J]
  sigma <- exp(parm[Data$J+1])
  ### Log of Prior Densities
  beta.prior <- sum(dnormv(beta, 0, 1000, log=TRUE))
  sigma.prior <- dhalfcauchy(sigma, 25, log=TRUE)
  ### Log-Likelihood
  mu <- tcrossprod(Data$X, t(beta))
  LL <- sum(dnorm(Data$y, mu, sigma, log=TRUE))
  ### Log-Posterior
  LP <- LL + beta.prior + sigma.prior
  Modelout <- list(LP=LP, Dev=-2*LL, Monitor=c(LP,sigma),
                   yhat=rnorm(length(mu), mu, sigma), parm=parm)
  return(Modelout)
}

set.seed(666)

############################  Initial Values  #############################
Initial.Values <- GIV(Model, MyData, PGF=TRUE)

###########################################################################
# Examples of MCMC Algorithms                                             #
###########################################################################

########################  Hit-And-Run Metropolis  #########################
Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
                     Covar=NULL, Iterations=200000, Status=500, Thinning=36,
                     Algorithm="HARM", Specs=NULL)
Consort(Fit)
plot(BMK.Diagnostic(Fit, batches=10))
PosteriorChecks(Fit)
caterpillar.plot(Fit, Parms="beta")
BurnIn <- Fit$Rec.BurnIn.Thinned
plot(Fit, BurnIn, MyData, PDF=FALSE)
Pred <- predict(Fit, Model, MyData)
??summary 
??predict
??plot
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
# Juxtapose

########################  Adaptive Griddy-Gibbs  ##########################
Fit1 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
     Algorithm="AGG", Specs=list(Grid=GaussHermiteQuadRule(3)$nodes,
     dparm=NULL, smax=Inf, CPUs=1, Packages=NULL, Dyn.libs=NULL))

##################  Adaptive Hamiltonian Monte Carlo  #####################
Fit2 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
     Algorithm="AHMC", Specs=list(epsilon=rep(0.02, length(Initial.Values)),
     L=2, Periodicity=10))

##########################  Adaptive Metropolis  ##########################
Fit3 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
     Algorithm="AM", Specs=list(Adaptive=500, Periodicity=10))

###################  Adaptive Metropolis-within-Gibbs  ####################
Fit4 <- LaplacesDemon(Model, Data=MyData, Initial.Values,
     Covar=NULL, Iterations=1000, Status=100, Thinning=1,
     Algorithm="AMWG", Specs=list(Periodicity=50))
#########################################


 Juxt <- Juxtapose(list(Fit1=Fit1, Fit2=Fit2, Fit3=Fit3, Fit4=Fit4)); Juxt
 plot(Juxt, Style="ISM")

