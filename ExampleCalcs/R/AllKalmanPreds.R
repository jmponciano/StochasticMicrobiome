load("allhmpdata3.0.RData")

# Creating a list with time series data that mirrors the list 'integrated.data2'
# but only contains time series of obs.-error corrected time series of abundances
# These corrected time series will be the Kalman-filter predictions computed
# using the Dennis and Ponciano 2014 paper. All the code for this paper is in 
# the folder ROUSS

#----------------------------------------------------------------------
#        PARAMETER ESTIMATION, PARAMETRIC BOOTSTRAP AND PREDICTIONS
#----------------------------------------------------------------------
# Before doing the calculations, the user has to specify ONLY the following 4 options:

# 1. Do you want to compute the ML estimates or the REML estimates?

method <- "REML" # alternatively, set method <- "ML"

# 2. Do you want to plot the predictions?
pred.plot <- "FALSE" # Set it to "FALSE" if you do not want to plot the predictions

# 3. Do you want to plot the parametric bootstrap distribution of the estimates?
pboot.plot <- "FALSE" # Set it to "FALSE" if you do not want to plot the bootstrap distribution of the estimates

# 4. How many bootstrap replicates?
NBoot <- 2; # 10 just to try, for formal results use 1000 or more 



#----------------------------------------------------------------------
#  Load all the OUSS model functions and needed packages
#----------------------------------------------------------------------
library("MASS")
source("ROUSS/ROUSSE-2.0.R")

all.Short.KalmanXs <- list()
all.Short.KalmanXs.CI1 <- list()
all.Short.KalmanXs.CI2 <- list()
all.Long.KalmanXs <- list()
all.Long.KalmanXs.CI1 <- list()
all.Long.KalmanXs.CI2 <- list()
all.pboot.cis <- list()

#nshort <- length(integrated.data2)

# Check:
#for(i in 1:nshort){
# 	one.womansdat <- integrated.data2[[i]]
# 	Tvec <- one.womansdat$raw.days
# 	
# 	print(sum(table(Tvec)>1))
#}



for(i in 1:nshort){
  
  one.womansdat <- integrated.data2[[i]]
  
  one.womdims <- dim(one.womansdat)
  last.spcol  <- one.womdims[2]-2
  first.spcol <- 9
  cols2extract <- first.spcol:last.spcol
  nspp <- length(cols2extract)

  # getting the time series time data
  Time.t <- correct.time(one.womansdat$raw.days)
  #correct.time(Time.t)
  print(i)
  print(rbind(one.womansdat$raw.days,Time.t))
   
  tt <- as.numeric(Time.t)-as.numeric(Time.t[1])
  long.t <- tt[1]:max(tt)
  
  ithlens     <- length(as.numeric(Time.t))
  ithlong.lens <- length(long.t)
  Short.KalmanXs     <- matrix(0,nrow=ithlens,ncol=nspp);colnames(Short.KalmanXs) <- names(one.womansdat[,cols2extract]);
  Short.KalmanXs.CI1 <- matrix(0,nrow=ithlens,ncol=nspp);colnames(Short.KalmanXs.CI1) <- names(one.womansdat[,cols2extract]);
  Short.KalmanXs.CI2 <- matrix(0,nrow=ithlens,ncol=nspp);colnames(Short.KalmanXs.CI2) <- names(one.womansdat[,cols2extract]);
  Long.KalmanXs     <- matrix(0,nrow=ithlong.lens,ncol=nspp);colnames(Long.KalmanXs) <- names(one.womansdat[,cols2extract]);
  Long.KalmanXs.CI1 <- matrix(0,nrow=ithlong.lens,ncol=nspp);colnames(Long.KalmanXs.CI1) <- names(one.womansdat[,cols2extract]);
  Long.KalmanXs.CI2 <- matrix(0,nrow=ithlong.lens,ncol=nspp);colnames(Long.KalmanXs.CI2) <- names(one.womansdat[,cols2extract]);
  ith.pboot.cis<- list() 
  
  for(j in 1:nspp){
           
    # Extracting each time series, one by one
    Observed.t <- one.womansdat[,cols2extract[j]] 

    #-------- Log-transform the observations to carry all the calculations in this program--------------------#
    log.obs    <- log(Observed.t);
    #---------------------------------------------------------------------------------------------------------#

    #-------------------------------------------------------------------------------------------------#
    #  5. RUN THE FOLLOWING LINE OF CODE TO COMPUTE THE ESTIMATES, PREDICTIONS,
    #     AND RUN A PARAMETRIC BOOTSTRAP. THE USER DOES NOT NEED TO MODIFY THIS LINE OF CODE.
    #     THE OUTPUT OF THE FUNCTION 'ROUSS.CALCS' IS A LIST AND THE USER CAN RETRIEVE EACH OF THE   
    #     LIST ELEMENTS PRINTED AND SAVED IN THE OBJECT "all.results". 
    #     THE 95% PARAMETRIC BOOTSTRAP FOR BOTH, THE PARAMETERS AND THE PREDICTIONS ARE COMPUTED 
    #     BY THE FUNCTION "ROUSS.CALCS" 	
    #     
    #     THE OUTCOME OF THE FUNCTION: A LIST WITH THE OBJECTS:
    #     $parms.est = mles of the OU state space model
    #     $lnLhat = the maximized log-likelihood
    #     $pbootmat = the parametric bootstrap mles matrix
    #     $pboot.cis = the parametric bootstrap confidence intervals of the OU model parameters
    #     $pboot.preds1 = the matrix of kalman estimates, columns are 1, time; 2, 2.5% pctle, 3, REMLES; 4, 97.5 pclte.
    #     $pboot.preds2 = the matrix of kalman estimates with the added missing time step predictions, 
    #       columns are 1, time; 2, 2.5% pctle, 3, REMLES; 4, 97.5 pclte.    
    #-------------------------------------------------------------------------------------------------#
    
    jthts.results <- ROUSS.CALCS(Yobs=log.obs,Tvec=as.numeric(Time.t), pmethod=method, nboot=NBoot, 
                                 plot.pred=pred.plot, plot.bootdists = pboot.plot);    
      
    ith.pboot.cis[[j]] <- jthts.results$pboot.cis
    Short.KalmanXs[,j] <- jthts.results$pboot.preds1[,3]
    Short.KalmanXs.CI1[,j] <- jthts.results$pboot.preds1[,2]
    Short.KalmanXs.CI2[,j] <- jthts.results$pboot.preds1[,4]

    Long.KalmanXs[,j] <- jthts.results$pboot.preds2[,3]
    Long.KalmanXs.CI1[,j] <- jthts.results$pboot.preds2[,2]
    Long.KalmanXs.CI2[,j] <- jthts.results$pboot.preds2[,4]

  }

  all.Short.KalmanXs[[i]] <- cbind(tt,Short.KalmanXs);
  all.Short.KalmanXs.CI1[[i]] <- cbind(tt, Short.KalmanXs.CI1);
  all.Short.KalmanXs.CI2[[i]] <- cbind(tt, Short.KalmanXs.CI2);
  all.Long.KalmanXs[[i]] <- cbind(long.t, Long.KalmanXs);
  all.Long.KalmanXs.CI1[[i]] <- cbind(long.t, Long.KalmanXs.CI1);
  all.Long.KalmanXs.CI2[[i]] <- cbind(long.t,Long.KalmanXs.CI2);
  all.pboot.cis[[i]] <- ith.pboot.cis
    
}

names(all.Short.KalmanXs) <- names(integrated.data2)
names(all.Short.KalmanXs.CI1) <- names(integrated.data2)
names(all.Short.KalmanXs.CI2) <- names(integrated.data2)
names(all.Long.KalmanXs) <- names(integrated.data2)
names(all.Long.KalmanXs.CI1) <- names(integrated.data2)
names(all.Long.KalmanXs.CI2) <- names(integrated.data2)
names(all.pboot.cis) <- names(integrated.data2)

save.image("allhmpdata5.0.RData")

 for(i in 1:nshort){
	
	 fname <- paste0("longkalman-",names(integrated.data2)[i] ,".txt")
	 write.matrix(x=all.Long.KalmanXs[[i]], file=fname)
	
 }



#save.image("allhmpdata4.0.RData")




