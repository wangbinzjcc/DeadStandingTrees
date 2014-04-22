###############################################################
#pcnm build , xianqianxuanze ,fangchafenjie; wangbinzjcc
###############
setwd('F:\\DataW\\DeadStandingTrees')
require(vegan)
#require(ade4)
dir()
#```
pcn0 <- read.csv("Pcnm.summ 10 m vectors .csv")
pcn01 <- cbind(expand.grid(y=1:30, x=1:50),pcn0[, -1])
p_Mor <- read.csv("Pcnm.summ 10 m Moran_I .csv")  
pcn1 <- pcn01[, p_Mor$Moran_I.Positive]
head(pcn1)[,1:10] ## ````
pcn11 <-pcn1[,8]
pp <- pcn1[pcn11>mean(pcn11),]
with(pp, plot(x, y,xlim=c(0,50),ylim=c(0,30))) 

#```
env00 <- read.csv("LG.topographic 10m poly4 2013-9-4.csv")
env11 <- env00[env00$meanelev>250,]
with(env11,plot(x,y))

#
LGtab <- read.csv('LGdata_table.csv')
head(LGtab)
la0 <- LGtab$X00枯立木 / rowSums(LGtab[-c(1:3)])
summary(la0)
dea11 <- LGtab[la0>0.04,]
with(dea11, plot(ox, oy))
#```
DeadDat <- read.csv("Deadtree-dataframe-dbhAre-abun.csv")
DeadDat <- DeadDat[order(DeadDat$x, DeadDat$y), ]
head(DeadDat)
 
summary(DeadDat$dbh.Are)
LL <- DeadDat
plot(LL$x, LL$y, cex=(LL$dbh.Are)^0.2)
#
#################################################################
#################################################################
#
### xiangqianxuanze diandiantu :  ~~~~ 
#  install.packages("PCNM", repos="http://R-Forge.R-project.org")
library(PCNM)    
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
# ################################### 
##################################################################
####################################################################
LGtab <- read.csv('LGdata_table.csv')
head(LGtab)
tab0 <- LGtab[, -c(1:5)] 
#########################
head(env01)
env00 <- read.csv("LG.topographic 10m poly4 2013-9-4.csv")
env001 <- env00[, -c(1:3)]
envtab01 <- cbind(env001, tab0 ,Abu=rowSums(tab0))
head(envtab01)
#
DeadDat <- read.csv("Deadtree-dataframe-dbhAre-abun.csv")
DeadDat <- DeadDat[order(DeadDat$x, DeadDat$y), ]
head(DeadDat)
LG0a <- log1p(DeadDat[, -c(1:3)])$dbh.Are
LG0a <- LGtab$X00枯立木/ rowSums(LGtab[-c(1:3)])
head(LG0a)
#
s.anova <- c() ; PCNM.sign <- list()
#mite.h.det <- resid(lm(as.matrix(mite.h) ~ ., data=mite.xy))
env.rda <- vegan::rda(LG0a, envtab01)
env.R2a <- RsquareAdj(env.rda)$adj.r.squared 
(LGa.env.fwd <- forward.sel(LG0a, envtab01, 
                            adjR2thresh=env.R2a))
#  write.csv(LG.PCNM.fwd, "PCNM.fwd_DeadDataFrame.csv") 
#   LG.PCNM.fwd <- read.csv("PCNM.fwd_DeadDataFrame.csv")
###############################################################
PCNM.pos <- pcn1[, -c(1:2)]
LG0a <- log1p(DeadDat[, -c(1:3)])$dbh.Are
head(LG0a)
#
s.anova <- c() ; PCNM.sign <- list()
#mite.h.det <- resid(lm(as.matrix(mite.h) ~ ., data=mite.xy))
PCNM.rda <- vegan::rda(LG0a, as.matrix(PCNM.pos))
LGa.R2a <- RsquareAdj(PCNM.rda)$adj.r.squared 
(LGa.PCNM.fwd <- forward.sel(LG0a, as.matrix(PCNM.pos), 
                            adjR2thresh=LGa.R2a))
#
 
hist(LGa.PCNM.fwd$order, breaks=20)
#  write.csv(LG.PCNM.fwd, "PCNM.fwd_DeadDataFrame.csv") 
#  LG.PCNM.fwd <- read.csv("PCNM.fwd_DeadDataFrame.csv")
###################
env_a_fwd <- envtab01[, LGa.env.fwd$variables]
head(env_a_fwd)
pcn_a_fwd <- PCNM.pos[, LGa.PCNM.fwd$variables]
head(pcn_a_fwd)
#
(LGa.varpart <- varpart(LG0a, env_a_fwd, pcn_a_fwd) )
plot(LGa.varpart, digits=3)
# ########################################################################
##########################################################################
###############################################################

head(env01)
env001 <- env00[, -c(1:3)]
env01 <- cbind(env001, tab0 ,Abu=rowSums(tab0))
#
head(DeadDat)
LG0b <- log1p(DeadDat[, -c(1:3)])$Abu
head(LG0b)
#
s.anova <- c() ; PCNM.sign <- list()
#mite.h.det <- resid(lm(as.matrix(mite.h) ~ ., data=mite.xy))
env.rda <- vegan::rda(LG0b, env01)
env.R2a <- RsquareAdj(env.rda)$adj.r.squared 
(LG.env.fwd <- forward.sel(LG0b, env01, 
                           adjR2thresh=env.R2a))
LGb.env.fwd <- LG.env.fwd
#  write.csv(LG.PCNM.fwd, "PCNM.fwd_DeadDataFrame.csv") 
#   LG.PCNM.fwd <- read.csv("PCNM.fwd_DeadDataFrame.csv")
######### 
PCNM.pos <- pcn1[, -c(1:2)]
LG0b <- log1p(DeadDat[, -c(1:3)])$Abu
head(LG0b)
#
s.anova <- c() ; PCNM.sign <- list()
#mite.h.det <- resid(lm(as.matrix(mite.h) ~ ., data=mite.xy))
PCNM.rda <- vegan::rda(LG0b, as.matrix(PCNM.pos))
res.anova <- anova.cca(PCNM.rda)
LG.R2a <- RsquareAdj(PCNM.rda)$adj.r.squared 
(LG.PCNM.fwd <- forward.sel(LG0b, as.matrix(PCNM.pos), 
                            adjR2thresh=LG.R2a))
LGb.PCNM.fwd <- LG.PCNM.fwd
#
 hist(LGb.PCNM.fwd$order, breaks=20)
#  write.csv(LG.PCNM.fwd, "PCNM.fwd_DeadDataFrame.csv") 
#   LG.PCNM.fwd <- read.csv("PCNM.fwd_DeadDataFrame.csv")
env_fwd <- env01[, LGb.env.fwd$variables]
head(env_fwd)
pcn_fwd <- PCNM.pos[, LGb.PCNM.fwd$variables]
head(pcn_fwd)
#
(LGb.varpart <- varpart(LG0b, env_fwd, pcn_fwd) )
plot(LGb.varpart, digits=3)
#########################################################################

13.68+23.31-27

5.17+11.75-13.89

#################################################################
# ###############################################################
### fangchafenjiequanquantu : ~~~~~
head(LG)
dim(LG)
dim(env00)
head(env00)

head(env01)
head(env_fwd)
head(LG0)

head(tab_fwd)
LG.varpart <- varpart(LG0$dbh.Are, env_fwd, pcn_fwd) 

LG.varpart <- varpart(LG0, env_fwd, pcn_fwd, tab_fwd) 

varpart(LG, env01, pcn_fwd) 

windows(title="Mite - environment - PCNM variation partitioning",6,3)
par(mfrow=c(1,2),mex= 0.39,mar=c(1,1, 1, 1) , pty="m")
showvarparts(4)
plot(LG.varpart, digits=4)


###############################################
#   绘制维恩图    wangbin 2012-5-27 17:55:41    #
###############################################
x11(4.1,3.1)
par(mex=0.1)

plot(1:10,type="n", axes="" )  
box() 
symbols(c(4.5,6.2), c(6,6), circles = c(4.5,6.8),  inches = 1.0
        ,fg =c(2,4), add=T)
a=0.3
b=10.2
c=43.2

text(3.4, 6.5, "[a] ")
text(3.4, 5.8, paste(a,'%',sep=''))

text(5, 6.5, "[b] ")
text(5, 5.8, paste(b,'%',sep=''))

text(7, 6.5, "[c] ")
text(7, 5.8, paste(c,'%',sep=''))

text(7.8, 1.5, paste("不可解释部分[d] =",100-a-b-c,"%",sep=''))

text(2,8,paste('地形因子\n解释部分=\n',a+b,'%',sep=''))

text(9.4, 8, paste("PCNM\n解释部分\n=",b+c,'%',sep=''),col=4)
 


text(7.8, 1.5, paste("不可解释部分[d] ",sep=''))

text(2,7,paste('地形因子\n 解释部分=',sep=''))

text(9.4, 7, paste("PCNM\n =解释部分",sep=''),col=4)



#################################

########### 

# #######################################################
#  CA PaiXu BankuaiTuͼ            2012??5??27??20:36:42           # 
# #################################################################
# ------------
head(LG); head(env11) ; head(pcn_fwd[,1:10])

dat.ca <-  summary(vegan::cca(LG))$sites 
LG_ca.sum <-summary(vegan::cca(LG))
LG_ca.sum[6]
summary(summary(vegan::cca(LG)))
head(dat.ca)
# ~~~~~~~~~~
XY.10 <-  env00[,c("x","y")]
x11()
par(mfrow=c(3,2), mex=0.05 ,mar=c(0, 0, 20 , 0) , pty="m" )

for( j in 1:6 ){
  s.value(XY.10, dat.ca[,j] , method="greylevel" , csize=0.3, clegend=1, grid=FALSE, include.origin=FALSE, addaxes = FALSE)
  title(main=paste("CA",j),cex.main=1 )
}
dim(dat.ca)

Rda_all <- vegan::cca(LG, cbind(pcn_fwd, env11)) 
summary(summary(Rda_all))
names(summary(Rda_all))
x11()
plot(Rda_all)

for(i in 1:15){
write.csv(summary(Rda_all)[[i]],paste(names(summary(Rda_all))[i],".csv")) 
}

## ~~~~~

names00 <- colnames(LG)[1:50]
spe.rda <- vegan::cca(LG)
x11()
  par( mex= 0.39,mar=c(5.6, 5.6, 5, 0.5) , pty="m" )
  
  plot(spe.rda, type='n',scaling=1); grid()
  abline(h=0,v=0, lty=2)
  spe.sc <- scores(spe.rda, choices=1:2, scaling=1, display="sp")
  points( spe.sc[,1], spe.sc[, 2], pch=3, col=4, cex=0.3 ) 
  
  text(spe.sc[names00,1]-0.09, spe.sc[names00, 2]-0.04,names00, cex=0.8, col=2)  
  title(main="物种组成 CA排序图", cex=0.8) 


plot(LG.rda)

# #######################################-------------------------------
#CCA a+b+c KeJieShi Tu : ~~~~~
nam00 <- colnames(dat.ca) 
png( "LGCABanKuaiTu22.png" ,  width = 600, height =600)
x11()
XY.10 <- cbind(as.numeric(substr(rownames(LG),1,2)), as.numeric(substr(rownames(LG),3,5)))
par(mfrow=c(6,3), mex=0.05 ,mar=c(0, 0, 20 , 0) , pty="m" )
for( i in nam00[1:6]){
  dat.abc <- rda(dat.ca[,i], cbind(pcn_fwd, env11)) 
  dat.ab <- rda(dat.ca[,i], env11)
  dat.c <- rda(dat.ca[,i], pcn_fwd, env11)
  names(summary(dat.c))
  RsquareAdj(dat.abc)
  dat.000 <- list(dat.abc, dat.ab, dat.c)
  dat.001 <- sapply(dat.000,function(x)RsquareAdj(x)$adj.r)
  if(any(is.na(dat.001))) {dat.001[which(is.na(dat.001))] <- 
                             as.numeric(RsquareAdj(dat.000[[which(is.na(dat.001))]])[1])*0.95}
  d001 <- round(dat.001,4)*100
  for( j in 1:3 ){
    su.d <- summary(dat.000[[j]])$sites[,1]
    s.value(XY.10, su.d, method="greylevel" , csize=0.3, clegend=1, grid=FALSE, include.origin=FALSE, addaxes = FALSE)
    title(main=paste(i, "explained by", c("[a+b+c]", "[a+b]", "[c]")[j],"= ",d001[j],"%"),cex.main=1 )
  }  }

dev.off()

# #################################################################
#   R2 PingFangTuͼ                 wangbin 2012??5??27??21:13:47           #
# #################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rd.c00 <- rda(LG, pcn_fwd, env11)
rd.ab00 <- rda(LG, env11)
summary(summary(rd.c00))

rd.c01 <- summary(rd.c00)$constraints 

rd.ab01 <- summary(rd.ab00)$constraints

r.squ <- c()
for(i in 1:50){
  r.squ[i] <- RsquareAdj(rda(rd.ab01, PCNM.pos[, 10*i+(-9:0)]))$adj.r.squared 
}

r.squ0 <- c()
for(i in 1:50){
  r.squ0[i] <- RsquareAdj(rda(rd.c01, PCNM.pos[, 10*i+(-9:0)]))$adj.r.squared 
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
x11(6,3.5)
par( mex=0.5, mar=c(5,6,2,1))
plot(r.squ,type="o",pch=2, col=4, cex=0.8
     , xlab="PCNM block" , ylab= expression(Adjusted*"  "*R^2)) 
points(r.squ0,type="o",pch=16, col=2, cex=0.8)
abline(h=0, type="l", lty=2)

leg.txt <-   expression("Adjusted  " * R^2 * "  for  [a+b]" , 
                        "Adjusted  " * R^2 * "  for  [c]")  

legend(20, 0.25,leg.txt , pch = c(2,16), col = c(4, 2),
       bty="n", cex = 0.8,lty=1, lwd=1 )

# ######################################################




