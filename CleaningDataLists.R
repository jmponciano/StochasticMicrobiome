###  Functions to aid in the data loading process

### This function takes a list of lists
### where each "sub-list" contains three matrices.
### Each one of these three matrices has 1 rep of the counts
### for every bacteria.  Their first column are the days (1 to 70)
### Returns two lists: the first one is another list of lists with 
### the same matrices, but without the rows with NA's 
### (that is, when no data was gathered) and the second one
### is another list of lists where each sublist contains three
### vectors, the vectors of removed rows so that I can go
### back to the metadata and pull the rows interesting for the
### analysis.

nas.cleaning <- function(womentsdat){

	n <- length(womentsdat)
	out.ts <- list()
	out.narows <- list()  
  
  for(i in 1:n){
    
    ts.mats <- womentsdat[[i]]
  	orig.tvec <- ts.mats[[1]][,1]
  	# These three vectors will be used for the "emptyness vec", which happens if all NA's
  	na.emptyvec1 <- rep(0,length(orig.tvec))
  	na.emptyvec2 <- rep(0,length(orig.tvec))
  	na.emptyvec3 <- rep(0,length(orig.tvec)) 	

	out.narows[[i]] <- list()

    # make three time series matrices (rows=days, cols=spp counts)
    # that DO NOT have the rows with na's
    has.nas1 <- sum(is.na(ts.mats[[1]][,2]))
    if(has.nas1>0){
      na.ind1 <- which(is.na(ts.mats[[1]][,2])==TRUE,arr.ind=TRUE);
      tsmat1 <- ts.mats[[1]][-na.ind1,];
      out.narows[[i]][[1]] <- na.ind1;
      na.emptyvec1[na.ind1] <- orig.tvec[na.ind1];
      rm(na.ind1);
    		}else{tsmat1 <- ts.mats[[1]]
    			out.narows[[i]][[1]] <- 0;
            }
    
    has.nas2 <- sum(is.na(ts.mats[[2]][,2]))
    if(has.nas2>0){
      na.ind2 <- which(is.na(ts.mats[[2]][,2])==TRUE,arr.ind=TRUE);
      tsmat2 <- ts.mats[[2]][-na.ind2,];
      out.narows[[i]][[2]] <- na.ind2;
      na.emptyvec2[na.ind2] <- orig.tvec[na.ind2]
      rm(na.ind2);	
    	}else{tsmat2 <- ts.mats[[2]];
    		out.narows[[i]][[2]] <- 0;
    	}
    
    has.nas3 <- sum(is.na(ts.mats[[3]][,2]))
    if(has.nas3>0){
      na.ind3 <- which(is.na(ts.mats[[3]][,2])==TRUE,arr.ind=TRUE)
      tsmat3 <- ts.mats[[3]][-na.ind3,];
      out.narows[[i]][[3]] <- na.ind3;
      na.emptyvec3[na.ind3] <- orig.tvec[na.ind3]
      rm(na.ind3)
   	 	}else{tsmat3 <- ts.mats[[3]]
    		out.narows[[i]][[3]] <- 0
    	}

    emptyness.test <- (sum((orig.tvec!=na.emptyvec1)&(orig.tvec!=na.emptyvec2)&(orig.tvec!=na.emptyvec3), na.rm=TRUE)==0)
	
	
    if(emptyness.test==TRUE){out.ts[[i]] <- list(ltsmat1=as.matrix(0),ltsmat2=as.matrix(0),ltsmat3=as.matrix(0)); 
    		
    		}else{
    
	    t1 <- tsmat1[,1];t2 <- tsmat2[,1];t3 <- tsmat3[,1];
    	mint1 <- min(t1,na.rm=TRUE);t1 <- t1-mint1;
    	lent1 <- length(t1);
   	 
    	if(sum(is.na(t1))>0){
    		nainds <- which(is.na(t1)==TRUE,arr.ind=TRUE)
    		for(j in 1:length(nainds)){
    			bad.ind <- nainds[j]
    			t1[bad.ind] <- t1[(bad.ind-1)]+1
    		}
    		rm(nainds)
    	}
    
    	tsmat1[,1] <- t1

    	mint2 <- min(t2,na.rm=TRUE);t2 <- t2-mint2;
    	lent2 <- length(t2);
    
    	if(sum(is.na(t2))>0){
    		nainds <- which(is.na(t2)==TRUE,arr.ind=TRUE)
    		for(j in 1:length(nainds)){
    			bad.ind <- nainds[j]
    			t2[bad.ind] <- t2[(bad.ind-1)]+1
    		}
   		    	rm(nainds)
    	}

     	tsmat2[,1] <- t2


    	mint3 <- min(t3,na.rm=TRUE);t3 <- t3-mint3;
    	lent3 <- length(t3);
    
    	if(sum(is.na(t3))>0){
    		nainds <- which(is.na(t3)==TRUE,arr.ind=TRUE)
    		for(j in 1:length(nainds)){
    			bad.ind <- nainds[j]
    			t3[bad.ind] <- t3[(bad.ind-1)]+1
    		}
    	    	rm(nainds)
    	}
    
    	tsmat3[,1] <- t3
  
  
    
    	## Next step:
    	## If there is one species that has all zeros throughout, eliminate it
    	## from the matrix of time series data
    	ncols <-ncol(tsmat1)
    	# Are there columns with all 0's in each matrix?
    	zero.onespp1 <- sum(apply(tsmat1[,-c(1,ncols)],2,sum)==0)
    	zero.onespp2 <- sum(apply(tsmat2[,-c(1,ncols)],2,sum)==0)        
    	zero.onespp3 <- sum(apply(tsmat3[,-c(1,ncols)],2,sum)==0)    
   	 
   	 	max.zeros <- max(c(zero.onespp1,zero.onespp2,zero.onespp3))
		where.max <- which(c(zero.onespp1,zero.onespp2,zero.onespp3)==max.zeros, arr.ind=TRUE)   	 
   	 	ndif <- length(where.max)
		
		tsmat11 <- tsmat1
		tsmat22 <- tsmat2
		tsmat33 <- tsmat3		
   	 	
   	 	if(ndif==1){
   	 		if(where.max==1){all.zero.ind <- which(apply(tsmat1,2,sum)==0,arr.ind=TRUE)
   	 			}else if(where.max==2){
   	 				all.zero.ind <- which(apply(tsmat2,2,sum)==0,arr.ind=TRUE)
   	 			}else if(where.max==3){
   	 				all.zero.ind <- which(apply(tsmat3,2,sum)==0,arr.ind=TRUE)
   	 		}
		tsmat11 <- tsmat1[,-all.zero.ind]
		tsmat22 <- tsmat2[,-all.zero.ind]
		tsmat33 <- tsmat3[,-all.zero.ind]
		
		}else{   	 
   	 
    		if(zero.onespp1>0){
      			zero.ind1 <- which(apply(tsmat1,2,sum)==0,arr.ind=TRUE)
      			tsmat11   <- tsmat1[,-zero.ind1]
      			print(paste0("1st matrix for woman id # ", i))
      			print(focal.spp[zero.ind1])
    		}
    
    		if(zero.onespp2>0){
      			zero.ind2 <- which(apply(tsmat2,2,sum)==0,arr.ind=TRUE)
      			tsmat22   <- tsmat2[,-zero.ind2]
      			print(paste0("2nd matrix for woman id # ", i))
      			print(focal.spp[zero.ind2])
    		}
    
    		if(zero.onespp3>0){
      			zero.ind3 <- which(apply(tsmat3,2,sum)==0,arr.ind=TRUE)
      			tsmat33   <- tsmat3[,-zero.ind3]
      			print(paste0("3rd matrix for woman id # ", i))
      			print(focal.spp[zero.ind3])
    		}
		
		}
		    
    	test <- sum(length(tsmat11[1,]) == c(length(tsmat22[1,]),length(tsmat33[1,])))
    	if(test!=2){print("Watch it! A species is missing in one replicate but not in the others");
                print(paste0("This happens for woman id # ", i))}
    
    	tsmat11[tsmat11==0] <- .Machine$double.eps;tsmat11[1,1] <- 0;
    	tsmat22[tsmat22==0] <- .Machine$double.eps;tsmat22[1,1] <- 0;    
    	tsmat33[tsmat33==0] <- .Machine$double.eps;tsmat33[1,1] <- 0;
    
    	ltsmat1 <- cbind(tsmat11[,1],log(tsmat11[,-1]))
    	ltsmat2 <- cbind(tsmat22[,1],log(tsmat22[,-1]))
    	ltsmat3 <- cbind(tsmat33[,1],log(tsmat33[,-1]))

    	out.ts[[i]] <- list(ltsmat1=ltsmat1, ltsmat2 = ltsmat2, ltsmat3 = ltsmat3)
	}
  }
  
  names(out.ts) <- names(womentsdat)
  names(out.narows) <- names(womentsdat)
  

  return(list(out.ts=out.ts,out.narows=out.narows))
  
}