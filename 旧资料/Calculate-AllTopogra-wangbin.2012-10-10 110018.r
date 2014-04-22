###########################################################################################
#  calculate AllTopogra,eg. slope aspect  convex elev                                         #
#                                             wangbinzjcc 2013-1-5 17:25:30
#                                             wangbinzjcc 2014-1-9 2:37:36
###########################################################################################
           
# one slope of a grid-plot ; 
# elev3 is dataframe of coordinates and elevation at 4 corner.
 SlopeCalc0 <- function(i, elev3 = elev3, gridsize = gridsize){
       z  <- sort(elev3[-i, 3])
       slope0 <- acos(gridsize / 
                    sqrt(gridsize^2 + (z[3] - z[1])^2  + (z[2] - z[1])^2))
       result.slope <- 180 * slope0 / pi         
       return(result.slope)                                   }

# aspect of a grid-plot
 AspectCalc0 <- function(elev3 = elev3){
       fx <- sum(elev3[elev3$x == max(elev3$x), 3]  - 
                   elev3[elev3$x == min(elev3$x), 3] )
       fy <- sum(elev3[elev3$y == max(elev3$y), 3] - 
                   elev3[elev3$y == min(elev3$y), 3] )
       result.aspect <- 180 - atan(fx / fy) * 180 / pi + 90 * fx / abs(fx)
       return(result.aspect)           }

# calulate elev， slop，aspect of a grid-plot  ;
# elev1 is dataframe of coordinates and elevation of a grid-plot
 TopograResult0 <- function(no.i = no.i, elev1=elev1, gridsize = gridsize){     
       elev2 <- elev1[elev1$x != min(elev1$x) & elev1$y != min(elev1$y), ]            # Remove the two edges
       corn0.i <- elev2[no.i, 1:2]                                            # Extract the maximum coordinates of quadrats 
       xy0 <- expand.grid( rbind(corn0.i - gridsize, corn0.i) )               # Coordinates of the four corners of the quadrat  
       elev3 <- elev1[match(paste(xy0$x, xy0$y), paste(elev1$x, elev1$y)), 1:3]       # Extract quadrat corners of altitude
       m.ele1 <- mean(elev3[, 3])                                              # Calculate the average altitude    
       slop1 <- mean(sapply(1:4, SlopeCalc0, elev3 = elev3, gridsize = gridsize))      # Calculate the slope  
       aspec1 <- AspectCalc0(elev3 = elev3)                                    # Calculate the aspect 
       result.topogra <- data.frame(meanelev = m.ele1, slope = slop1, aspect = aspec1)  # 3 terrain factor results 
       return(result.topogra)                                           }
                                                                 
# Calculated quadrats irregularities (1); 
# result1 three terrain factors result in the quadrats
 ConvexCalc0 <- function(n1.i = n1.i, result1 = result1, gridsize = gridsize){
       result1 <- result1[, c("x", "y", "meanelev")]             # Organize data
       qua0.i <- as.numeric(result1[n1.i, 1:2])                  # Extract quadrats coordinates
       size1 <- seq(-gridsize, gridsize, by = gridsize)
       xy1 <- expand.grid(data.frame(x = qua0.i[1] + size1, y = qua0.i[2] + size1))     # 9 adjacent quadrats coordinates 
       xy2 <- xy1[-match(paste(qua0.i[1], qua0.i[2]), paste(xy1$x, xy1$y)), ]           #  8 adjacent quadrats coordinates
       neighb.elev <- result1[match(paste(xy2$x, xy2$y), paste(result1$x, result1$y)), 3]    # 8 quadrats adjacent average altitude
       result.convex <- as.numeric(result1[n1.i, 3]) - mean(unlist(neighb.elev))        # Bump degrees ( 1 )        
       return(result.convex)                                                 }

##  ##  ##

# Calculation the four topographical factors within the entire sample ; 
# elevfile plots the original coordinates and altitude data
 AllTopograCalcResult <- function(elevfile = elevfile, gridsize = 20){ 
       elev0 <- elevfile[order(elevfile$x, elevfile$y), c("x", "y", "elev")]     # Organize data

       elev1 <- elev0[(elev0$x-min(elev0$x)) %% gridsize == 0 & 
                        (elev0$y-min(elev0$y)) %% gridsize == 0, ]               # Converted into the corresponding scale data

       elev2 <- elev1[elev1$x != min(elev1$x) & elev1$y != min(elev1$y), ]       # Remove the two edges of each sample ( maximum internal coordinate square number)


       result0 <- t(sapply(1 : dim(elev2)[1], 
                           TopograResult0, elev1 = elev1, gridsize = gridsize))     # Calculate three topographic factors within the entire sample

       result1 <- cbind(no. = 1:nrow(elev2), elev2[, 1:2], result0)
 
       convex.dat0 <- sapply(1:dim(elev2)[1], 
                             ConvexCalc0, result1 = result1, gridsize = gridsize)    # Bump degrees ( 1 )
       
       logic.1 <- result1$x == min(result1$x) | result1$x == max(result1$x) | 
                      result1$y == min(result1$y) | result1$y == max(result1$y)        # The edge of the plot
       siz.cut <- cut(gridsize, c(0, 15, 30, 45, 90, "inf"), c(5, 10, 20, 25, 50))     # Different scales corresponding to different values
       res.xy  <- result1[logic.1, c("x", "y")] - as.numeric(as.character(siz.cut))    # Coordinates of the center of the edge quadrat
      
       logic.2 <- match(paste(res.xy$x, res.xy$y), paste(elev0$x, elev0$y))              # Edge quadrats No.       
       convex.dat0[logic.1] <- elev0[logic.2, 3] - unlist(result1[logic.1, "meanelev"])  # Plot edge irregularities (2)
 

       result.al0 <- cbind(result1, convex = convex.dat0)                    # All terrain factor results
       result.all <- matrix(unlist(result.al0), nrow = nrow(result.al0))
       colnames(result.all) <- names(result.al0)     

       return(round(result.all, 3))                                 }                        
                                                            
##############################################################################

###############################################################################
#    setwd("F:\\DataW\\DeadStandingTrees")
#    elevfile <- read.csv("elev5m2013.csv")   
#    head(elevfile)  # included "x","y","elev"
#    length( unique(elevfile$y) ) ;length( unique(elevfile$x) ) ;dim(elevfile)  

#   All.result <- AllTopograCalcResult(elevfile=elevfile ,gridsize=20)   
#    write.csv(All.result, 'topograCal-20m.csv')
###
#    All.result <- AllTopograCalcResult(elevfile=elevfile ,gridsize=5) 
#     write.csv(All.result, 'topograCal-5m.csv')
################################################################################

   
