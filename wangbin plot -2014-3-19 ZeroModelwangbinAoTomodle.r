############################################################
#  Zero model for deadtree    wangbinzjcc 2014-3-10 20:54:12
############################################################
rm(list=ls())
ls()
####################
dir()
dea.5 <- read.csv("DeadData-5m-abun.csv")
head(dea.5)
dea5 <- dea.5$Abu
t5 <- table(dea5)
#
windowsFonts(Times = windowsFont("Times New Roman"))
tiff('dead trees Frequency 5m.tiff',
     family="Times", pointsize=8,
     width = 80, height = 60, units = "mm",
     res=600,compression = "lzw")
#
op00 <- par(mar=c(5,5,5,1), mex=0.43)
plot(t5,xlab="Abundance",ylab="Frequency", main = "Abundance by 5m * 5m")
#
dev.off()

#####################

windowsFonts(Times = windowsFont("Times New Roman"))
tiff('dead trees abundance 5m.tiff',
     family="Times", pointsize=8,
     width = 80, height = 60, units = "mm",
     res=600,compression = "lzw")
#
op00 <- par(mar=c(3,3,5,1), mex=0.43)
image(x=1:100, y=1:60, z=matrix(dea10,nrow=100,ncol=60,byrow=T),
      main = "Dead trees' abundance by 5m * 5m")
dev.off()

#########################################################################
dir()
dea.10 <- read.csv("DeadData-10m-abun.csv")
head(dea.10)
dea10 <- dea.10$Abu
t10 <- table(dea10)
#
windowsFonts(Times = windowsFont("Times New Roman"))
tiff('dead trees Frequency 10m.tiff',
     family="Times", pointsize=8,
     width = 80, height = 60, units = "mm",
     res=600,compression = "lzw")
#
op00 <- par(mar=c(5,5,5,1), mex=0.43)
plot(t10,xlab="Abundance",ylab="Frequency", main = "Abundance by 10m * 10m")
#
dev.off()

#####################

windowsFonts(Times = windowsFont("Times New Roman"))
tiff('dead trees abundance 10m.tiff',
     family="Times", pointsize=8,
     width = 80, height = 60, units = "mm",
     res=600,compression = "lzw")
#
op00 <- par(mar=c(3,3,5,1), mex=0.43)
image(x=1:50, y=1:30, z=matrix(dea10,nrow=50,ncol=30,byrow=T),
      main = "Dead trees' abundance by 10m * 10m")
dev.off()

#######################################################################################

#
setwd("F:\\DataW\\DeadStandingTrees")
##
dea.20 <- read.csv("DeadData-20m-abun.csv")
head(dea.20)
dea20 <- dea.20$Abu
t20 <- table(dea20)

########################################################

windowsFonts(Times = windowsFont("Times New Roman"))
tiff('dead trees Frequency 20m.tiff',
     family="Times", pointsize=8,
     width = 80, height = 60, units = "mm",
     res=600,compression = "lzw")
#
op00 <- par(mar=c(5,5,5,1), mex=0.43)
plot(t20,xlab="Abundance",ylab="Frequency", main = "Abundance by 20m * 20m")
#
dev.off()

########################################################


mean(t20);var(t20);sd(t20)
plot(t20)
lt20 <- log(t20)
mean(lt20);var(lt20);sd(lt20)
plot(log(t20))
###
windowsFonts(Times = windowsFont("Times New Roman"))
tiff('dead trees abundance 20m.tiff',
     family="Times", pointsize=8,
     width = 80, height = 60, units = "mm",
     res=600,compression = "lzw")
#
op00 <- par(mar=c(3,3,5,1), mex=0.43)
image(x=1:25, y=1:15, z=matrix(dea.20$Abu,nrow=25,ncol=15,byrow=T),
      main = "Dead trees' abundance by 20m * 20m")
dev.off()
###############################################################
#####  aoto regression  ~~~~~~~~
head(dea.20)
i=5
XX <- dea.20
for(i in 1:dim(dea.20)[1]){
  x0 <- XX[i, 'x']
  y0 <- XX[i, 'y']
  logi.0 <- abs(XX$x - x0)<=1 & abs(XX$y-y0)<=1
  logi.0[i] <- FALSE
  XX[i,"mean.neigh"] <- round(mean(XX[logi.0,'Abu']), 2)
}
dea.20 <- XX
write.csv(dea.20, "DeadData-20m-abun-aoto.csv")
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
#
convex.c <- CenterScale(topo$convex, Binary="centerscale") 
convex2.c <- CenterScale(topo$convex^2, Binary="centerscale") 
convex3.c <- CenterScale(topo$convex^3, Binary="centerscale") 
#
aspect.c <- CenterScale(cos(topo$aspect*pi/180), Binary="centerscale") 
aspect2.c <- CenterScale(cos(topo$aspect*pi/180)^2, Binary="centerscale") 
aspect3.c <- CenterScale(cos(topo$aspect*pi/180)^3, Binary="centerscale") 
 #
X.c <-  data.frame(elev.c, elev2.c, elev3.c, slope.c, slope2.c, slope3.c,
            convex.c, convex2.c, convex3.c, aspect.c, aspect2.c, aspect3.c)
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
  alpha <- parm[1]
  beta <- parm[c(2,3,4,5,6,7,8,9,10,11,12,13)] 
  rho <- parm[14]
  sigma <- exp(parm[15])
  # priors distribution
  alpha.priors <- dnormv(alpha, 0, 1000, log=T)
  beta.priors <- sum(dnormv(beta, 0, 1000, log=T))
  rho.priors <- dnormv(rho, 0, 1000, log=T)
  sigma.priors <- dgamma(sigma, 25, log=T)
  # log-likelihood
  mu <- alpha + tcrossprod(Data$X,t(beta)) + rho*Data$NI+ sigma
  LL <- sum(dnorm(Data$y, mu, sigma, log=T))
  # log-Posterior
  LP <- LL + alpha.priors + beta.priors + rho.priors + sigma.priors
  #    
  Modelout <- list(LP=LP, Dev=-2*LL,  Monitor=c(LP),  
                   yhat= rnorm(length(mu), mu, sigma), parm=parm)
  return(Modelout)  
}
#
# ########################################
#
###  test data ~~~~~~~~~~~~~~
N <- 375 # length(dea20)
J <- 12
x.test <- matrix(1,N,J)
for(j in 1:J){x.test[,j] <- rnorm(N, runif(1,-10,10), runif(1,0.1,3))}
NI.test <- rnorm(N,10,3)
alpha.test <- 100
beta.test <- runif(J, -100, 100)
rho.test <- 50
ee <- rnorm(N, 0, 0.5) 
y.test <- alpha.test + tcrossprod(x.test,t(beta.test)) +rho.test*NI.test + ee
hist(y.test)
#
 
#
x.wb <- x.test
y.wb <- y.test
NI.wb <- NI.test
###############################################################
### ture data ~~~~~~~~~~~~~~~
#
#
x.wb <-  X.c 
y.wb <- dea.20$Abu
NI.wb <- dea.20$mean.neigh
#
 N <- length(y.wb)
 J <- 12  # beta numbers
#
###################################
#
parm.wb <- c(alpha=0, beta=rep(0,J),rho=0, sigma=log(1))
PGF <- function(MyData){c(rnormv(1+MyData$J+1,0,1000), log(rhalfcauchy(1, 25)))}
MyData <- list(J=J, N=length(y.wb), PGF=PGF,
               mon.names=c("LogPosterior"),
               parm.names=names(parm.wb), X=x.wb, NI=NI.wb, y=y.wb)
#
Initial.Values <- GIV(Model, MyData, n=10000, PGF=TRUE)
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
# 100; beta.test
# 
Fit2
Fit6
Fit7
Fit13
Fit18
#
Fit <- Fit7
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
#   autoregression
###################################################
Data <- list(X=x.wb, y=y.wb, NI=NI.wb)
Fit <- Fit7
beta.data <- as.data.frame(Fit$Summary2)$Mean[1:15]
alpha <- beta.data[1]
beta <- beta.data[2:13]
rho <- beta.data[14]
sigma <- beta.data[15]
mu <- alpha + tcrossprod(Data$X,t(beta)) + rho * Data$NI + sigma
yhat <- rnorm(length(mu), mu, sigma)
y.residual <- Data$y - yhat -3

image(x=1:25, y=1:15, z=matrix(Data$y,nrow=25,ncol=15,byrow=T))
image(x=1:25, y=1:15, z=matrix(y.residual,nrow=25,ncol=15,byrow=T))

hist(yhat)
hist(Data$y)
hist(y.residual)
plot(yhat, y.residual)
abline(lm(yhat~y.residual),col=3)
abline(h=0, col=2)
plot(NI.wb,y.residual)
abline(h=0, col=2)
################################
ls()
dea20.xy <- dea.20
head(dea20.xy)
dea.20
library(gstat)
coordinates(dea20.xy) = ~x+y
v0 = variogram(Abu ~1, dea20.xy)
plot(v0)
v1= variogram(y.residual~1,dea20.xy)
plot(v1, col=2)
plot(variogram(y.residual~1, dea20.xy, cloud=TRUE))
##
require(spdep)
bci.xy <- expand.grid(x=1:10,y=1:20)
bci.xy <- dea.20[c('x','y')]
bci.nb <- dnearneigh(as.matrix(bci.xy),0,1)
plot(bci.nb, bci.xy)
 #
bcisp.I <- sp.correlogram(bci.nb,dea20,order=25,method="I",zero.policy=T)
bciYres.I <- sp.correlogram(bci.nb,y.residual,order=25,method="I",zero.policy=T)
plot(bcisp.I, ylim=c(-0.1,0.4))
plot(bciYres.I, ylim=c(-0.1,0.4))
#
bcisp.C <- sp.correlogram(bci.nb,dea20,order=25,method="C",zero.policy=T)
bciYres.C <- sp.correlogram(bci.nb,y.residual,order=25,method="C",zero.policy=T)
plot(bcisp.C, ylim=c(0.65,1.1))
plot(bciYres.C, ylim=c(0.65,1.1))
#
#####################################################






AIC
BIC



