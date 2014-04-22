#####################################################
# from deadtree data to data-DBH.area-Abundance   
#                            wangbinzjcc  2013-12-30
######################################################
setwd("F:\\DataW\\DeadStandingTrees")
dat <- read.csv('deadtree_data.csv')
dat0 <- dat[is.na(dat$bra),]
head(dat0)
dim(dat0)
mean(dat0$dbh)
max(dat0$dbh)
###################
gx <- gsub('[N|n|m][G|a|c]([0-9]{2})[0-9]{5}','\\1',dat0$no.)
gy <- gsub('[N|n|m][G|a|c][0-9]{2}([0-9]{2})[0-9]{3}','\\1',dat0$no.)
#
gxy <- gsub('[N|n|m][G|a|c]([0-9]{4})[0-9]{3}','\\1',dat0$no.)
#
dbhArea <- tapply(dat0$dbh, gxy, function(x){sum(pi*(x/2)^2)})
speAbu <- table(gxy)
#
exp0 <- expand.grid(x0=c(paste(0,1:9,sep=''), 10:50), y0=c(paste(0,1:9,sep=''), 10:30)) 
xy00 <- paste(exp0$x0, exp0$y0, sep='')
no.0 <- match(x=names(dbhArea), table=xy00)
Are00 <- rep(0,length(xy00));names(Are00) <- xy00
Are00[no.0] <- dbhArea; DbhAre.dat <- Are00
length(DbhAre.dat)
Are00[no.0] <- speAbu; SpeAbu.dat <-  Are00 
Dead.dat <- as.data.frame(cbind(exp0, DbhAre.dat,SpeAbu.dat))
rownames(Dead.dat) <- xy00; colnames(Dead.dat) <- c('x', 'y', 'dbh.Are', 'Abu')
head(Dead.dat)
write.csv(Dead.dat, 'DeadData-dbhAre-abun.csv')
hist(log(DbhAre.dat))
hist(log(SpeAbu.dat))
#
############################
 