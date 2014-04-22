#############################################################
#   统计绘制 二维图 斑块图 CA轴   wangbin2012年5月18日10:56:44    #
############################################################# 
# ******************************************
#        程序 Plot.Cca  Function          ~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  Plot.Cca.WB<-function(spe.rda, names00 = names(ngdata01)){
     par( mex= 0.39,mar=c(5.6, 5.6, 5, 0.5) , pty="m" )

     plot(spe.rda, type='n',scaling=1); grid()
     abline(h=0,v=0, lty=2)
     spe.sc <- scores(spe.rda, choices=1:2, scaling=1, display="sp")
     points( spe.sc[,1], spe.sc[, 2], pch=3, col=4, cex=0.3 ) 
 
     text(spe.sc[names00,1]-0.09, spe.sc[names00, 2]-0.04,names00, cex=1, col=1)  
     text(spe.sc[4,1]-0.09, spe.sc[4, 2]-0.04,names00[4], cex=1, col=2)  
     title(main="物种组成 CA排序图", cex=0.8)                 }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~end
 
# **************************************************
   require(ade4)

   require(vegan)
# ~~~~~~~~~~~~~~~~~~
   data <- read.csv("小重值数据2012-5-27 142333.csv")
   IV.nam <- read.csv("IV.4全部 2012-4-28224922.csv")
   head(data )
  dim(data )   
    head(IV.nam)
# ~~~~~~~~~~~~~~~~~~
  x <- as.numeric(substr(data[,"X.1"],1,2))
  y <- as.numeric(substr(data[,"X.1"],3,5))
  logi.a <- x>14 & x<=34 & y>10 & y<=30

   Nam00 <- IV.nam$X[IV.nam$IV>0.1]
# ~~~~~~~~~~~~~~~~~~~

   ngdata01 <- data[logi.a, match(Nam00, colnames(data))]+0.00001 
  head(ngdata01)
   ngda.ran <-decostand(ngdata01,"range")
   ngda.hel <-decostand(ngda.ran,"hellinger",MARGIN=2)
   dim(ngda.hel)
   head(ngda.hel)
 
# ~~~~~~~~~~~~~~~~~~~
   dec0 <- vegan::cca(ngda.hel, scale=T)

#     
  # write.csv(summary(dec0)$species,"cca.summ.species.csv")
  # write.csv(summary(dec0)$sites,"cca.summ.sites.csv")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#        二维图 CA轴
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
png("二维图 CA排序.png",  width = 600, height = 600 ) 
 
  Plot.Cca.WB(spe.rda=dec0, names00 = names(ngdata01)[1:50]) 

dev.off()

# ************************************************************* 
#        斑块图 CA轴
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 summary(summary(dec0)) 

  datas.all <- summary(dec0)$site 
  datas.sp <- summary(dec0)$species 

head(datas.sp)
head(datas.all)
names(datas.all)
 ma00 <- which.max(as.data.frame(datas.all)[,4])
 mi00 <- which.min(as.data.frame(datas.all)[,4])
datas.all[,4][c(ma00,mi00)] <- datas.all[,4][c(ma00,mi00)]/3

   XY.10 <- data.frame(x[logi.a], y[logi.a])

    names.0 <- paste("物种组成 CA轴",1:6, sep="")

# ----------------------------
   png("斑块图-物种 CA轴 .png", width = 600, height = 600 ) 
    
  par(mfrow=c(2,2), mex= 0.05,mar=c(0, 0, 23,0) , pty="m" ) 
   for(i in 1:4){
      s.value(XY.10, datas.all[,i] , method="greylevel", csize=0.7, clegend=1.8 , grid=FALSE, include.origin=FALSE, addaxes = FALSE)
      title(main=names.0[i])  
                }

     dev.off()           

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




#######################################################
#     空间自相关图  CA轴
#   ################################################### 

#************************************************** 
# FUNCTION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
  
 Plot.Maron.WB <- function(subs.cor.i,names.i){
  
    aa.c <- subs.cor.i
    plot(aa.c, main=names.i)
    aa <- as.data.frame(print(aa.c))                           

    t.a <- aa[,5]
    t.y <- aa[,1]
    t.ad <- aa[,3]  
    text.0 <-  character()
    length(text.0) <- length(t.a)
    text.0[t.a<100] <- ""
    text.0[t.a<0.05] <- "*"
    text.0[t.a<0.001] <- "**"
    text.0[t.a<0.001] <- "***"
    
    T.Y <- t.y+(1+abs(t.ad/max(t.ad)))*(1+(1-max(t.y))*0.00014)^10000*max(t.y)*0.07-0.03
    text(1:length(t.y), T.Y, text.0 , col=2 , cex=0.9)
                                  }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~END
library(PCNM)	 
require(vegan)

#*********************************************************** 

nb1 <- dnearneigh(as.matrix(XY.10), 0, 1)  ### 尺度为 0~1 个10m
mite.env <- round(datas.all,3)
 
subs.correlog <- list()
for( i in 5:6){
subs.dens <- mite.env[,i]
subs.correlog[[i]] <- sp.correlogram(nb1,subs.dens, order=25, method="I", 
  zero.policy=TRUE)     
              }
 
#***********************************************************
#   绘图 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

png("自相关图-物种CA轴.png", width = 600, height = 600 )
par(mfrow=c(2,2), mex= 0.39,mar=c(5.6, 5.2, 5, 0.5) , pty="m" )

 for( i in 1:4){
           subs.cor.i <- subs.correlog[[i]]
           names.i <- names.0[[i]]
           Plot.Maron.WB(subs.cor.i, names.i)
               }

dev.off()

# ***********************************************************



