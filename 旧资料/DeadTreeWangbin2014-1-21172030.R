#########################################################
#    DeadStandingTrees      wangbinzjcc 2014-1-21 17:17:49
#######################################################
#
setwd("F:\\DataW\\DeadStandingTrees")
##
dea.10 <- read.csv("DeadData-10m-abun.csv")
head(dea.10)
dim(dea.10)
#
topo.10 <- read.csv("LG.topographic 10m poly4 2013-9-4.csv")
# 
LGtab <- read.csv('LGdata_table.csv')
head(LGtab)
pi0 <- LGtab$X00枯立木
pi1 <- rowSums(LGtab[-c(1:3)])
#######################################################
head(topo.10)
de0 <-topo.10[,-c(1,2,3)]
de0$ymat <- cbind(pi0, pi1 - pi0)
glm(ymat~., family=binomial(link="logit"), data=de0)
####################################################### 

require(LaplacesDemon)     
#
Model <- function(parm, Data){
# parameters
  beta0 <- parm[1] 
  beta1 <- parm[2]  
  beta2 <- parm[3]
 # beta3 <- parm[4]
 # beta4 <- parm[5]
 # beta5 <- parm[6]
 # beta6 <- parm[7]
 # beta7 <- parm[8]
 # beta8 <- parm[9]
#
  x1 <- Data$meanelev
  x2 <- Data$meanelev2
#  x3 <- Data$convex
#  x4 <- Data$convex2
#  x5 <- Data$slope
#  x6 <- Data$slope2
#  x7 <- Data$asp.cos
#  x8 <- Data$asp.sin
#  
  r <- Data$r0
  n <- Data$n0
  X1 <- x1 - mean(x1)
  X2 <- x2 - mean(x2)
  # X3 <- x3 - mean(x3)
  # X4 <- x4 - mean(x4)
  # X5 <- x5 - mean(x5)
  # X6 <- x6 - mean(x6)
  # X7 <- x7 - mean(x7)
  # X8 <- x8 - mean(x8)
# log-likelihood
  logit.pi <- beta0 + beta1*X1 +  beta2*X2
     # + beta3*X3 + beta4*X4 + beta5*X5  +  beta6*X6 + beta7*X7 + beta8*X8
  pi <- invlogit(logit.pi)
  LL <- sum(dbinom(x=r, size=n, prob=pi, log=T))
  # priors
  beta0.priors <- dnormv(beta0, 0.0, 1.0E3, log=T)
  beta1.priors <- dnormv(beta1, 0.0, 1.0E3, log=T)
  beta2.priors <- dnormv(beta2, 0.0, 1.0E3, log=T)
#  beta3.priors <- dnormv(beta3, 0.0, 1.0E3, log=T)
#  beta4.priors <- dnormv(beta4, 0.0, 1.0E3, log=T)
#  beta5.priors <- dnormv(beta5, 0.0, 1.0E3, log=T)
#  beta6.priors <- dnormv(beta6, 0.0, 1.0E3, log=T)
#  beta7.priors <- dnormv(beta7, 0.0, 1.0E3, log=T)
#  beta8.priors <- dnormv(beta8, 0.0, 1.0E3, log=T)
  #
  LP <- LL + beta0.priors + beta1.priors +  beta2.priors 
  # + beta3.priors + beta4.priors + beta5.priors + beta6.priors + beta7.priors + beta8.priors   
  #  
  Modelout <- list(LP=LP, Dev=-2*LL,  Monitor=c(LP),  
                   yhat= rbinom(n=length(pi),size=n, prob=pi), parm=parm)
  return(Modelout)  
}
#
#
# ########################################
pi0 <- LGtab$X00枯立木
pi1 <- rowSums(LGtab[-c(1:3)])
de0 <-topo.10[, c("meanelev","meanelev2", "convex", "convex2", 
                  "slope","slope2", "asp.cos", "asp.sin")]
de0 <- apply(de0, 2, as.numeric)
#
parm <- c(beta0 =-3.65, beta1 =0, beta2 =0
          # , beta3 =0, beta4 =0, beta5 =0, beta6 =0, beta7 =0, beta8 =0
          )
Data <- list(N= length(pi0),mon.names=c("LogPosterior"), parm.names=names(parm),
             r0=pi0, n0=pi1, meanelev=de0[1], meanelev2=de0[2]
           #  , convex=de0[3], convex2=de0[4], 
           #   slope=de0[5], slope2=de0[6], asp.cos=de0[7], asp.sin=de0[8]
             )
#   Model(parm, Data)
# Run LaplacesDemon
out1 <- LaplaceApproximation(Model, Data=Data, parm=parm, Iterations=2000, Method="HAR")
#
parm <- as.data.frame(out1$Summary2)$Mean
#
#
out <- LaplacesDemon(Model, Data=Data, Initial.Values=parm,
                     Iterations=500000, Status=2000, Thinning=50)
#
out
#
plot(out, BurnIn=50, Data, PDF=F)
parm <- as.data.frame(out$Summary1)$Mean[1:5]
#
#
out2 <- LaplacesDemon.hpc(Model, Data=Data, Initial.Values=parm, Chains=5, CPUs=2,
                          Iterations=500000, Status=2000, Thinning=50)
#
#
plot(out2, BurnIn=50, Data, PDF=F)

as.data.frame(out2[[3]]$Summary2)$Mean[1:5]
#
# 
#
plot(out, BurnIn=50, Data, PDF=F)
#



######################################################################################
