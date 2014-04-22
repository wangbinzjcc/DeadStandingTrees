############################################################## 
# DeadStandingTree  Bayesian model   wangbinzjcc2014-1-9
##############################################################
setwd("F:\\DataW\\DeadStandingTrees")
dat <- read.csv('deadtree_data.csv')
dat0 <- dat[is.na(dat$bra),]
head(dat0)
##########################
gx <- gsub('[N|n|m][G|a|c]([0-9]{2})[0-9]{5}','\\1',dat0$no.)
gy <- gsub('[N|n|m][G|a|c][0-9]{2}([0-9]{2})[0-9]{3}','\\1',dat0$no.)
#
gxy <- gsub('[N|n|m][G|a|c]([0-9]{4})[0-9]{3}','\\1',dat0$no.)
#
speAbu <- table(gxy)
#
exp0 <- expand.grid(x0=c(paste(0,1:9,sep=''), 10:50), y0=c(paste(0,1:9,sep=''), 10:30)) 
xy00 <- paste(exp0$x0, exp0$y0, sep='')
no.0 <- match(x=names(speAbu), table=xy00)
Are00 <- rep(0,length(xy00));names(Are00) <- xy00
Are00[no.0] <- speAbu; SpeAbu.dat <-  Are00 
Dead.dat <- as.data.frame(cbind(exp0, SpeAbu.dat))
rownames(Dead.dat) <- xy00; colnames(Dead.dat) <- c('x', 'y', 'Abu')
head(Dead.dat)
write.csv(Dead.dat, 'DeadData-10m-abun.csv')
##########################################################
############################################################

setwd("F:\\DataW\\DeadStandingTrees")
dat <- read.csv('deadtree_data.csv')
dat0 <- dat[is.na(dat$bra),]
head(dat0)

gx0 <- gsub('[N|n|m][G|a|c]([0-9]{2})[0-9]{5}','\\1',dat0$no.)
gy0 <- gsub('[N|n|m][G|a|c][0-9]{2}([0-9]{2})[0-9]{3}','\\1',dat0$no.)
#
gx <- ceiling(as.numeric(gx0)/2) 
gx[nchar(gx)==1] <- paste('0', gx[nchar(gx)==1], sep='')
gy <- ceiling(as.numeric(gy0)/2) 
gy[nchar(gy)==1] <- paste('0', gy[nchar(gy)==1], sep='')
gxy <- paste(gx,gy,sep='')
#
speAbu <- table(gxy)
length(speAbu)
#
exp0 <- expand.grid(x0=c(paste(0,1:9,sep=''), 10:25), y0=c(paste(0,1:9,sep=''), 10:15)) 
xy00 <- paste(exp0$x0, exp0$y0, sep='')
no.0 <- match(x=names(speAbu), table=xy00)
Are00 <- rep(0,length(xy00));names(Are00) <- xy00
Are00[no.0] <- speAbu; SpeAbu.dat <-  Are00 
Dead.dat <- as.data.frame(cbind(exp0, SpeAbu.dat))
rownames(Dead.dat) <- xy00; colnames(Dead.dat) <- c('x', 'y', 'Abu')
head(Dead.dat)
write.csv(Dead.dat, 'DeadData-20m-abun.csv')
########################################################

setwd("F:\\DataW\\DeadStandingTrees")
dat <- read.csv('deadtree_data.csv')
dat0 <- dat[is.na(dat$bra),]
head(dat0)
#
dat0$y[1974] <- dat0$y[1974]+0.01
gx <- ceiling(as.numeric(dat0$x)/5)-1 
gx[nchar(gx)==1] <- paste('0', gx[nchar(gx)==1], sep='')
gy <- ceiling(as.numeric(dat0$y)/5)-1
gy[nchar(gy)==1] <- paste('0', gy[nchar(gy)==1], sep='')
#
speAbu <- table(paste(gx,gy)) 
#
exp0 <- expand.grid(x0=c(paste(0,0:9,sep=''), 10:99), y0=c(paste(0,0:9,sep=''), 10:59)) 
head(exp0)
xy00 <- paste(exp0$x0,exp0$y0)
 
no.0 <- match(x=names(speAbu), table=xy00)
Are00 <- rep(0,length(xy00));names(Are00) <- xy00
Are00[no.0] <- speAbu; SpeAbu.dat <-  Are00 
Dead.dat <- as.data.frame(cbind(exp0, SpeAbu.dat))
rownames(Dead.dat) <- xy00; colnames(Dead.dat) <- c('x', 'y', 'Abu')
head(Dead.dat)
write.csv(Dead.dat, 'DeadData-5m-abun.csv')
#######################################################

dir()
##
dea.05 <- read.csv("DeadData-5m-abun.csv")
dea.10 <- read.csv("DeadData-10m-abun.csv")
dea.20 <- read.csv("DeadData-20m-abun.csv")
#
topo.10 <- read.csv("LG.topographic 10m poly4 2013-9-4.csv")
topo.20 <- read.csv("topograCal-20m.csv")
topo.05 <- read.csv("topograCal-5m.csv")
# 
#######################################################
dead0 <- c(dea.05$Abu, dea.10$Abu, dea.20$Abu)
dead1 <- cbind(dea.abu=dead0, 
        qua=c(rep('a',dim(dea.05)[1]), rep('b',dim(dea.10)[1]), rep('c',dim(dea.20)[1]))
               ) 
######################################################################

top.05.0 <- (topo.05[,c('slope', 'convex')])^3
top.10.0 <- topo.10[,c('slope3', 'convex3')]; names(top.10.0) <- c('slope', 'convex')
top.20.0 <- (topo.20[,c('slope', 'convex')])^3
topo11 <- rbind(top.05.0, top.10.0, top.20.0)
topo12 <- cbind(topo11, 
         qua1=c(rep('a',dim(top.05.0)[1]), rep('b', dim(top.10.0)[1]), rep('c',dim(top.20.0)[1]))
                )      
dim(topo12)
dead.env.dat <- cbind(dead1, topo12)
write.csv(dead.env.dat, 'dead.env.dat.wb.2014.csv')
#######################################################################
##############################################################
#######################################################################

require(LaplacesDemon) 

#
枯立木多度： Y ~ Posi(λ)
地形因子：  λ =  exp(α + X β);    i = 1,2;   j = 1,2,3
参数分布：
α.j   ~  Normal( μ.a, σ.a)        # 3个尺度
β.ij  ~  Normal( μ.b, σ.b)       # 3个尺度×2个参数
误差分布：       
μ.a ~ Normal(0,1000)
μ.b ~ Normal(0,1000)
#
σ.a ~ Gamma(25)
σ.a ~ Gamma(25)
#############
setwd("F:\\DataW\\DeadStandingTrees")
#
colnames(Xmat);y0=De;y0
write.csv(De,"yDead.csv")
write.csv(Xmat,"Xmat.csv")
Xmat0 <- read.csv("Xmat.csv")
yDe <- read.csv("yDead.csv")
head(Xmat0)
head(yDe)
########################################################################
########################################################################
########################################################################
#
setwd('F:\\GitHub\\DeadStandingTrees')

Xmat0 <- read.csv("Xmat.csv")

yDe <- read.csv("yDead.csv")
head(Xmat0)
head(yDe)
####################
require(LaplacesDemon)          #  R package
#

parm=c(alpha=1:3,beta=4:9, miu.a=10,miu.b=11,sigm.a_log=log(12),sigm.b_log=log(13))


parm=c(alpha=c(-0.05,-0.07,-0.1),beta=c(-0.06,-0.1,-0.06,-0.09,-0.1,-0.08), miu.a=0,miu.b=3.00,sigm.a_log=3,sigm.b_log=8)
parm.names <- as.parm.names(parm)
Data=list(x=Xmat0[,-1], parm.names=parm.names, mon.names=c("sigma.a","sigma.b"),
          y=yDe$x)

Model=function(parm, Data){
  # Priors
  alpha = parm[1:3]    # 3
  beta  = parm[4:9]    # 6
  #
  miu.a = parm[10]
  miu.b = parm[11]
  sigm.a = exp(parm[12])    # 1
  sigm.b = exp(parm[13])   # 1
  
  #
  alpha.prior <- sum(dnormv(alpha, miu.a, sigm.a, log=T))
  beta.prior <- sum(dnormv(beta, miu.b, sigm.b, log=T))
  miu.a.prior <- dnormv(miu.a, 0, 1000, log=T) 
  miu.b.prior <- dnormv(miu.b, 0, 1000, log=T) 
  sigm.a.prior <- dgamma(sigm.a, 25, log=T)
  sigm.b.prior <- dgamma(sigm.b, 25, log=T)  
  #
  mu= exp(rowSums(Data$x * c(alpha,beta)))
  LL <- sum(dpois(Data$y,mu,log=T))
  #
  LP <- LL + alpha.prior + beta.prior + miu.a.prior+ miu.b.prior + sigm.a.prior + sigm.b.prior
  
  # Model out
  yhat <- rpois(length(mu) ,mu)
  Modelout <- list(LP=LP, Dev= -2*LL, Monitor=c(sigm.a, sigm.b),
                   yhat=yhat, parm=parm)
  return(Modelout)              
}
###
out <- LaplacesDemon(Model, Initial.Values=parm, Data=Data,
                     Iterations=1000000, Status=200, Thinning=30)


#
# save.image("F:/GitHub/DeadStandingTrees/1000000timesBayes.RData")
# load("F:/GitHub/DeadStandingTrees/1000000timesBayes.RData")
out
# dput(out,"10000000timesLaplacesDenon")
# out0 <- dget("10000000timesLaplacesDenon")
summary(out)
rownames(summary(out))
#
o1 <- out[['Posterior1']]
o1 <- as.data.frame(o1)
head(o1)
aaa <- o1$alpha1
plot(aaa, type='o')
hist(aaa)
summary(aaa)
#
apply(o1, 2, quantile, probs=c(0.5, 0.05, 0.95))
out[['Summary2']]
###
quantile(aaa,probs=c(0.05,0.95))

out[[i]];i=i+1
summary(out[[16]] )
#
#
for(i in rownames(summary(out))){
  aa <- out[[i]]
  try(write.csv(aa, paste(i, '.csv'))) 
}
#

out[['Summary2']]


####################################################################



