###############################################################
#  Dead standing trees ppp points patterns 王斌 2013-12-3 10:27:30            #
###############################################################
setwd("F:\\DataW\\DeadStandingTrees")
dir()
ls()
dat <- read.csv('lgdat2013.csv')
########
if(any(dat$dbh < 1)){dat$dbh[dat$dbh<1]=1}
w00　<- duplicated(paste(dat$x, dat$y))
dat$x[w00] <- dat$x[w00] + runif(sum(w00),0,0.1)
head(dat)
###
deadtre <- subset(dat, sp=='00枯立木' & is.na(bra))
head(deadtre)
dim(deadtre)
#write.csv(deadtre,'deadtree_data.csv')

#
dat01 <- subset(deadtre, sp=='00枯立木' & dbh>=1)
head(dat01)

#
dat02 <- subset(deadtre, sp=='00枯立木' & dbh>=1 & dbh<5)

#
dat020 <- subset(deadtre, sp=='00枯立木' & dbh>=5 & dbh<10)

#
dat03 <- subset(deadtre, sp=='00枯立木' & dbh>=10 & dbh<20)

#
dat04 <- subset(deadtre, sp=='00枯立木' & dbh>=20)

###
tiff('DeadStandscattersplot.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")

plot(dat01$x, dat01$y, type='p', xlab='', ylab='')

dev.off()
###

tiff('DeadStandscattersplot-0-inf.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")

plot(dat01$x, dat01$y, type='n', xlab='', ylab='')
points(dat04$x, dat04$y, pch=7, col=2)
points(dat03$x, dat03$y, pch=6, col=3)
points(dat020$x, dat020$y, pch=1, col=4)
points(dat02$x, dat02$y, pch=3, col=1)

dev.off()
###################################################
tiff('DeadStandscattersplot＞20.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
plot(dat04$x, dat04$y, pch=7, col=2, main="dbh＞20cm")
dev.off()
###
tiff('DeadStandscattersplot-10-20.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
plot(dat03$x, dat03$y, pch=6, col=3, main='dbh-10-20cm')
dev.off()
###
tiff('DeadStandscattersplot-5-10.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
plot(dat020$x, dat020$y, pch=1, col=4, main='dbh-5-10cm')
dev.off()
###
tiff('DeadStandscattersplot-1-5.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
plot(dat02$x, dat02$y, pch=3, col=1, main='dbh-1-5cm')
dev.off()
############
###
tiff('DeadStandsHistPlot.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
hist(dat01$dbh)
dev.off()

#
hist(dat01$dbh)
hist(dat02$dbh)
hist(dat03$dbh)
hist(dat04$dbh)
###############################################

###############################################
require(spatstat)
##
xm1=dat01
xm2<-ppp(xm1$x,xm1$y,marks=as.factor(xm1$sp),c(0,500),c(0,300))
xm3.3<-envelope(xm2,pcf,nsim=199,nrank=5)
##
tiff('xm.dat01.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
par( mar=c(5.0,5.0,0.5,0.5 ),pty="m",mex=0.5) 
plot(xm3.3,legend=FALSE,xlab="",ylab="", 
     xlim=c(0,70), ylim=c(0,2), main="",bty="l")
legend("top", "枯立木 dbh>1 cm", inset =0.01,bty="n")
dev.off()
##
xm1=dat02
xm2<-ppp(xm1$x,xm1$y,marks=as.factor(xm1$sp),c(0,500),c(0,300))
xm3.3<-envelope(xm2,pcf,nsim=199,nrank=5)
##
tiff('xm.dat02.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
par( mar=c(5.0,5.0,0.5,0.5 ),pty="m",mex=0.5) 
plot(xm3.3,legend=FALSE,xlab="",ylab="", 
     xlim=c(0,70), ylim=c(0,2), main="",bty="l")
legend("top", "枯立木 1cm < dbh < 5cm", inset =0.01,bty="n")
dev.off()
##
xm1=dat020
xm2<-ppp(xm1$x,xm1$y,marks=as.factor(xm1$sp),c(0,500),c(0,300))
xm3.3<-envelope(xm2,pcf,nsim=199,nrank=5)
##
tiff('xm.dat03.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
par( mar=c(5.0,5.0,0.5,0.5 ),pty="m",mex=0.5) 
plot(xm3.3,legend=FALSE,xlab="",ylab="", 
     xlim=c(0,70), ylim=c(0,2), main="",bty="l")
legend("top", "枯立木 5cm < dbh < 10cm", inset =0.01,bty="n")
dev.off()
##
xm1=dat03
xm2<-ppp(xm1$x,xm1$y,marks=as.factor(xm1$sp),c(0,500),c(0,300))
xm3.3<-envelope(xm2,pcf,nsim=199,nrank=5)
##
tiff('xm.dat04.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
par( mar=c(5.0,5.0,0.5,0.5 ),pty="m",mex=0.5) 
plot(xm3.3,legend=FALSE,xlab="",ylab="", 
     xlim=c(0,70), ylim=c(0,2), main="",bty="l")
legend("top", "枯立木 10cm < dbh < 20cm", inset =0.01,bty="n")
dev.off()
##
xm1=dat04
xm2<-ppp(xm1$x,xm1$y,marks=as.factor(xm1$sp),c(0,500),c(0,300))
xm3.3<-envelope(xm2,pcf,nsim=199,nrank=5)
##
tiff('xm.dat05.tiff',
     width = 3000, height = 2800,res=600,compression = "lzw")
par( mar=c(5.0,5.0,0.5,0.5 ),pty="m",mex=0.5) 
plot(xm3.3,legend=FALSE,xlab="",ylab="", 
     xlim=c(0,70), ylim=c(0,2), main="",bty="l")
legend("top", "枯立木 dbh > 20cm", inset =0.01,bty="n")
dev.off()
##



#####################################################

