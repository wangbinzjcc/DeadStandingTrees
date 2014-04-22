###################################################################
#   nongang data for 20m * 20m abundance   wangbinzjcc 2014-1-25
###################################################################
#
setwd('F:\\DataW\\sapling-LiDX')
#
dir()
#
dat <- read.csv('lgdat2013.csv')
dat0 <- dat[is.na(dat$bra), ]
#
gx0 <- gsub('[N|n|m][G|a|c]([0-9]{2})[0-9]{5}','\\1',dat0$no.)
gy0 <- gsub('[N|n|m][G|a|c][0-9]{2}([0-9]{2})[0-9]{3}','\\1',dat0$no.)
#
gx <- ceiling(as.numeric(gx0)/2) 
gx[nchar(gx)==1] <- paste('0', gx[nchar(gx)==1], sep='')
gy <- ceiling(as.numeric(gy0)/2) 
gy[nchar(gy)==1] <- paste('0', gy[nchar(gy)==1], sep='')
gxy <- paste(gx,gy,sep='')
#
speAbu <- table(gxy, dat0$sp)
#
xx <- gsub('([0-9]{2})[0-9]{2}', '\\1', rownames(speAbu))
yy <- gsub('[0-9]{2}([0-9]{2})', '\\1', rownames(speAbu))
#
spe20Abu <- cbind(x=as.numeric(xx), y=as.numeric(yy), speAbu)
#
# write.csv(spe20Abu, 'LgData-20m-abun.csv')
####################################################
dir()
env <- read.csv("topograCal-20m.csv")
#
env$cos.asp <- cos(env$aspect*pi/180)
#
head(env); dim(env)
attach(env)
Tab.0 <- rep(NA,length(X))
elev0 <- c(meanelev-min(meanelev))
logi.1 <- elev0<76.3
logi.2 <- slope<31.52 & elev0<44.6
logi.3 <- elev0>=12.85
logi.4 <- cos.asp>0.008238
logi.5 <- elev0<135.3
logi.6 <- cos.asp>0.356
logi.7 <- slope<41.17
Tab.0[logi.1 & logi.2 & logi.3]   <- "A"
Tab.0[logi.1 & logi.2 & !logi.3]  <- "B"
Tab.0[logi.1 & !logi.2 & logi.4]  <- "D"
Tab.0[logi.1 & !logi.2 & !logi.4] <- "C"
Tab.0[!logi.1 & logi.5 & logi.6]  <- "F"
Tab.0[!logi.1 & logi.5 & !logi.6 & logi.7]  <- "G"
Tab.0[!logi.1 & logi.5 & !logi.6 & !logi.7] <- "E"
Tab.0[!logi.1 & !logi.5 ]  <- "H"
detach(env)


ha0 <- Tab.0
head(env)
xy.mrt <- cbind(x=env$x, y=env$y, tab=Tab.0)
xy.mrt <- as.data.frame(xy.mrt)
head(xy.mrt)
write.csv(xy.mrt,"NonggangXY.Mrt.csv")

##
############################# 
setwd("F:\\DataW\\DeadStandingTrees")
dea.5 <- read.csv("DeadData-5m-abun.csv")
head(dea.5)
head(xy.mrt)

aa <- unique(dea.5$y)
dea.5$y20 <- (dea.5$y %/% 4 + 1) *20
dea.5$x20 <- (dea.5$x %/% 4 + 1) *20
head(dea.5)
m00 <- match(paste(dea.5$x20,dea.5$y20,sep=''), paste(xy.mrt$x,xy.mrt$y,sep=''))
dea.5$Tab <- xy.mrt$tab[m00]
write.csv(dea.5,'DeadData-5m-abun-mrt2014-4-22.csv')
head(dea.5)
aa<-aov(Abu~Tab, data=dea.5[dea.5$Abu!=0,])
summary(aa)
aa<-aov(Abu~Tab, data=dea.5)
summary(aa)
aa
plot(dea.5$Abu~dea.5$Tab)
with(dea.5[dea.5$Abu!=0,], tapply(Abu, Tab, summary))
with(dea.5, tapply(Abu, Tab, summary))
################################################



unique(xy.mrt$y)
unique(dea.5$x+1)

dea.5[1:20,]






ha1 <- paste('a',c(ha0))
hab <- match(ha1, unique(ha1))
rn.ha <- 1:length(unique(c(hab)))
dim(hab) <- c(15, 25) 
image(1:25, 1:15,  t(hab), col=n.ha+1)
xy <- expand.grid(y=1:15, x=1:25)
text(xy$x, xy$y, ha0)



##########
require(mvpart)
#
# pick cv size and do PCA
#fit <- mvpart(data.matrix(spe20Abu[,-(1:2)])~meanelev+aspect+slope+ convex+cos.asp+
#                sin.asp, env,xv='p')
fit <- mvpart(data.matrix(spe20Abu[,-(1:2)])~meanelev + slope + convex + cos.asp +
                sin.asp, env,xv="1se")
str(fit)
ha0 <- fit$where
ha1 <- paste('a',c(ha0))
hab <- match(ha1, unique(ha1))
n.ha <- 1:length(unique(c(hab)))
dim(hab) <- c(15, 25) 
image(1:25, 1:15, t(hab), col=n.ha+1)
xy <- expand.grid(x=1:25, y=1:15)
text(xy$x, xy$y, c(t(hab)))
# write.csv(cbind(xy,c(t(hab))),'LG.hat.20m.csv')
########################################## 
dir()
dat00 <- read.csv("种子雨坐标.csv")
head(dat00)
gx0 <- dat00$xx
gy0 <- dat00$yy
gx <- ceiling(as.numeric(gx0)/20) 
gy <- ceiling(as.numeric(gy0)/20)  
#####
dat01 <- read.csv('LG.hat.20m.csv')
head(dat01)
m00 <- match(paste(gx,gy),paste(dat01$x,dat01$y))
dat02 <- dat01[match(paste(gx,gy),paste(dat01$x,dat01$y)), ]
head(dat02)
# write.csv(dat02,'hat.sapling.20m.csv')
###############
datayf <- read.csv('datayf-sp3.csv')
datayf <- apply(datayf,2,as.numeric)
head(datayf)
datyf.new <- matrix(0, 15*25,length(colnames(datayf))); colnames(datyf.new)<-colnames(datayf)
datyf.new[m00, ] <- datayf
head(datyf.new)
dim(datyf.new)
data003 <- cbind(dat01, total=rowSums(datyf.new), datyf.new)
head(data003)
#
with(subset(data003,苹婆>0),plot(x,y))
with(subset(data003,蚬木>0),plot(x,y))
#
# write.csv(data003,"LG.hat.spe.20m-wb-2014-1-25.csv")
########################################################
data(spider)
mvpart(data.matrix(spider[,1:12])~herbs+reft+moss+sand+twigs+
         water,spider)       # defaults
mvpart(data.matrix(spider[,1:12])~herbs+reft+moss+sand+twigs+
         water,spider,xv="p")  # pick the tree size
# pick cv size and do PCA
fit <- mvpart(data.matrix(spider[,1:12])~herbs+reft+moss+sand+
                twigs+water,spider,xv="1se",pca=TRUE)
rpart.pca(fit,interact=TRUE,wgt.ave=TRUE) # interactive PCA plot of saved multivariate tree











