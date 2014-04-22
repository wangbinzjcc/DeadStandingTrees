###############################################################
#pcnm build , xianqianxuanze ,fangchafenjie; wangbinzjcc
###############
require(vegan)
require(ade4)
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
env00 <- read.csv("地形因子 10 m  2012-10-10144056 .csv")
env00$asp.sin <- sin(env00$aspect*pi/180)
env00$asp.cos <- cos(env00$aspect*pi/180)
head(env00)  ## ````
env11 <- env00[env00$meanelev>250,]
with(env11,plot(x,y))

#```
LGdat <- read.csv("lgdat2013.csv")
LGdat$x[LGdat$x>=500] <- 499.999
LGdat$y[LGdat$y>=300] <- 299.999
LG_tab <- table(paste(LGdat$x %/% 10 + 1,  LGdat$y %/% 10 + 1)
                , LGdat$sp)
a <- rownames(LG_tab)
bb <- order(ox <- as.numeric(substr(a, 1, 2)), oy <-as.numeric(substr(a, 3, 5)))
LG_tab2 <- cbind(ox, oy, LG_tab)[bb,]
head(LG_tab2)[,1:10] ## ````
################################################
#write.csv(LG_tab2, "LGdata_table.csv")
################################################
LL <- LG_tab2[LG_tab2[,'蚬木'] > 0,]
with(as.data.frame(LL), plot(ox, oy))
#
### xiangqianxuanze diandiantu :  ~~~~ 
# 
library(PCNM)         
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
PCNM.pos <- pcn1[, -c(1:2)]
dim(PCNM.pos)
head(PCNM.pos)[1:10]  # !
#  PCNM FangkuanTu~~~~~
 length(dat.ord) 

head(pcn1)[,1:10]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  windows(title="PCNM variables (grid)",30,65)
par(mfrow=c(4,2),oma=c(0, 0, 0, 0), mex=2,mar=c(1,1, 1, 1) , pty="m" )
somePCNM2 <- c(1,5, 10, 20, 50, 100, 400,750 )
for(i in somePCNM2 ){
s.value(pcn1[,c('x','y')], PCNM.pos[, i], method="greylevel", 
        csub=0.5, csize=0.3, clegend=2, grid=FALSE, include.origin=FALSE, addaxes = FALSE)
title(main=paste("PCNM", i),cex.main=1.2 )
                    }
# ##~~~
data <-  LG_tab2[,-c(1:2)]
IV.nam <- rev(sort(table(LGdat$sp)))
Nam00 <- names(IV.nam)[IV.nam>15]
LG0 <- data[, match(Nam00, colnames(data))] 
LG <- LG0[, colSums(LG0)>0]
#
dim(LG)
head(LG)    # !
#
env11 <- env00[,-c(1:4)]
head(env11)  # !
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
windows(title="env11 variables (grid)",30,65)
par(mfrow=c(3,2),oma=c(0, 0, 0, 0), mex=1,mar=c(1,1, 1, 1) , pty="m" )
for(i in 1:dim(env11)[2]){
  s.value(env00[,c('x','y')], env11[, i], method="greylevel", 
          csub=0.5, csize=0.3, clegend=2, grid=FALSE, include.origin=FALSE, addaxes = FALSE)
  title(main=names(env11)[i],cex.main=1.2 )
}
# ~~~

windows(title="env11 hist",30,65)
par(mfrow=c(3,2),oma=c(0, 0, 0, 0), mex=1,mar=c(5,5, 1, 1) , pty="m" )
for(i in 1:dim(env11)[2]){
hist(env11[,i],col='light green',right=F,xlab = names(env11)[i], main=NULL)
}
# ################################### 

s.anova <- c() ; PCNM.sign <- list()
#mite.h.det <- resid(lm(as.matrix(mite.h) ~ ., data=mite.xy))
PCNM.rda <- vegan::rda(LG, as.matrix(PCNM.pos))
res.anova <- anova.cca(PCNM.rda)
LG.R2a <- RsquareAdj(PCNM.rda)$adj.r.squared 
(LG.PCNM.fwd <- forward.sel(LG, as.matrix(PCNM.pos), 
                              adjR2thresh=LG.R2a))
#  write.csv(LG.PCNM.fwd, "PCNM.fwd_LGLG.csv") 
#   LG.PCNM.fwd <- read.csv("PCNM.fwd_LGLG.csv")
# ########################################################################

head(LG.PCNM.fwd)
dat.ord <-  sort(LG.PCNM.fwd$order)

pcn_fwd <- PCNM.pos[,dat.ord]
names(pcn_fwd) <- gsub('vectors', 'pcnm', names(pcn_fwd))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d_c <- cut(dat.ord,c(0,200,400,600,800))

ord1 <- dat.ord[d_c==unique(d_c)[1]]
ord2 <- dat.ord[d_c==unique(d_c)[2]]
ord3 <- dat.ord[d_c==unique(d_c)[3]]
ord4 <- dat.ord[d_c==unique(d_c)[4]]

at00 <- seq(1,201, by=10) 
#Xiangqianxuanze diandiantu:  ~~~~~~~~~~

x11(6,3)
par(mfrow=c(4, 1),oma=c(0, 0, 0, 0), mex= 0.25, mar=c(5, 1, 1, 1)
    ,xaxt="s", pty="m" )
# 1:200 ~~~~~~~~ 
plot(1:200, rep(0,200), type="n",axes=F, xlab="" ,ylab="") ;box()
axis( 1,at00 , paste( at00-1 ),cex.axis=0.9 )
points(ord1+1, rep(0,length(ord1)), pch =20, cex=0.55,col=4)
# 201:400 ~~~~~~~ 
plot(1:200, rep(0,200), type="n",axes=F, xlab="" ,ylab="" ) ;box()
axis(1, at00 , paste(at00-1+200 ),cex.axis=0.9)
points(ord2-200+1, rep(0,length(ord2)), pch =20, cex=0.55,col=4)

# 401:600 ~~~~~~~~
plot(1:200, rep(0,200), type="n",axes=F, xlab="" ,ylab="" ) ;box()
axis(1, at00 , paste(at00-1+400 ),cex.axis=0.9)
points(ord3-400+1, rep(0,length(ord3)), pch =20, cex=0.55,col=4)

# 601:800 ~~~~~~~~
plot(1:200, rep(0,200), type="n",axes=F, xlab="" ,ylab="" ) ;box()
axis(1, at00 , paste(at00-1+600 ),cex.axis=0.9)
points(ord4-600+1, rep(0,length(ord4)), pch =20, cex=0.55,col=4)

# ###############################################################
### fangchafenjiequanquantu : ~~~~~

LG.varpart <- varpart(LG, env11, pcn_fwd) 

windows(title="Mite - environment - PCNM variation partitioning",6,3)
par(mfrow=c(1,2),mex= 0.39,mar=c(1,1, 1, 1) , pty="m")
x11()
showvarparts(2)
plot(LG.varpart, digits=2)


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




