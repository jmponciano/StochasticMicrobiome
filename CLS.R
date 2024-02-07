##### Program for CLS estimation of the MARSS(1) model
# statevar is the time series data species as columns and time as rows
# covariates are the covariates for the model
library(MASS)
library("rjags")


mars.cls	<-	function(comm.mat,covariate=NULL){
	
	#if(is.null(covariate)==FALSE){no.covar	<-	ncol(covariate)
	#NAs			<-	matrix(0,nrow=nrow(covariate),ncol=no.covar)
	#for(i in 1:no.covar){NAs.col	<-	which(is.na(covariate[,i]))
	#					 NAs[NAs.col,i]<-1}
	#which.nas	<-	rowSums(NAs)
	#comm.mat	<-	comm.mat[-which(which.nas>0),]	
	#covariate	<-	covariate[-which(which.nas>0),]}		 
	
	statevar<-comm.mat[2:nrow(comm.mat),]
	lagstate<-comm.mat[1:(nrow(comm.mat)-1),]
	Q	<-	length(statevar[,1])	# length of the time series
	P	<-	length(statevar[1,])	# Number of variates in the model
	
	if(is.null(covariate)==FALSE){R	<-	length(covariate[1,])}else {R <- 0}	# Number of covariates in the model 
	

	#Initialize variates
	E	<-	matrix(0,nrow=Q,ncol=P)	# residuals for least squares
	A	<-	rep(0,P)				# intercepts fo the variates
	B	<-	matrix(0,nrow=P,ncol=P)	# parameters for the variates
	if(is.null(covariate)==FALSE){C	<-	matrix(0,nrow=P,ncol=R)}else{C	<-	NULL}	# parameters for the covariates
	
	Yhat<-	matrix(0,ncol=P,nrow=Q)	# estimates
	varY<-	matrix(0,ncol=P,nrow=Q)	# variates for total variance
	varDY<-	matrix(0,ncol=P,nrow=Q)	# variates for per capita variance

	# loop over each variate

		for (dv in 1:P){
			Y	<-	statevar[,dv]
			# predictors
			if(is.null(covariate)==FALSE){X<-as.matrix(cbind(rep(1,Q),lagstate,covariate[-nrow(covariate),]))
				}else{X	<-	as.matrix(cbind(rep(1,Q),lagstate))}
		
			# CLS estimates
			beta	<-	ginv((t(X) %*% X)) %*% (t(X) %*% Y)
			# calculate residuals for that dependent variate
			Yhat [,dv]	<-	X %*% beta
			E[,dv]		<-	Y - Yhat[,dv]
		
			#parameter estimates
			A[dv]	<-	beta[1,1]
			B[dv,]	<-	beta[2:(P+1),1]
			if (is.null(covariate)==FALSE){
				C[dv,]	<-	t(beta[(P+2):nrow(beta),])}
			else{C	<- NULL}
		
			# accumulate Y variates for calculation of R^2 for observed vs. predicted
			varY[,dv]	<-	Y - mean(Y)
		
			# accumulate Y variates for calculation of R^2 for observed change in
			# x vs. predicted
			varDY[,dv]	<-	Y-lagstate[,dv]
		}
	
	# Assess the CLS fits

	# calculate CLS log-likelihood

	sigma	<-	t(E) %*% E/Q
	lnlike 	<-	-Q*(P/2)*log(2*pi)-Q/2*log(det(sigma)+0.001)-Q*P/2

	# calculate CLS explained variance

	varMatrix	<-	t(varY) %*% varY/Q
	R2	<-	1-(diag(sigma)/diag(varMatrix))

	# calculate CLS explained variance in change in x
	varMatrix_D	<-	t(varDY) %*% varDY/Q
	R2_D = 1-diag(sigma)/diag(varMatrix_D)

	# correlation matrix for process error
	d	<-	diag(1/sqrt(diag(sigma)))
	corrmatrix	<-	d*sigma*d

	# parameter count = # non-zero parameters in A,B,C and some of sigma
	par	<- sum(A!=0)+sum(sum(B!=0))+P*(P+1)/2 +sum(sum(C!=0))

	# AIC and BIC
	lsAIC	<-	-2*lnlike + 2*par
	lsBIC	<-	-2*lnlike + par*log(Q)

	spp.names <- colnames(comm.mat);
	colnames(B) <- spp.names;
	row.names(B) <- spp.names;

	results	<-	list(A=A,B=B,sigma=sigma,C=C,E=E,Yhat=Yhat,R2=R2,R2_D=R2_D,AIC=lsAIC,BIC=lsBIC)
	return(results)
}

#### Function stability() calculates the three measures of stability of a community as 
# proposed by Ives et al 2003. 
# B is a community matrix as an output from the mars.cls() function $B
# sigma is the sigma matrix from mars.cls() function

stability	<- function(Bmat,Sig){

	Bmat	<-	as.matrix(Bmat)
	nspecies <- nrow(Bmat)
	stab1 <- det(Bmat)^(2)
	lams.vec <- Re(eigen(Bmat)$values)
	mean.return <- max(lams.vec) # return rate for the mean
	var.return <- max(Re(eigen(kronecker(Bmat,Bmat))$values)) # return rate for the variance
	Vinf.true  <- matrix(solve(diag(nspecies*nspecies) - kronecker(Bmat,Bmat))%*%as.vector(Sig) , nrow=nspecies,ncol=nspecies)
	stab3 <- -sum(diag(Bmat))/sum(diag(Vinf.true))
	spp.contribs <- (1/lams.vec^2)/sum(1/lams.vec^2) 
	names(spp.contribs) <- colnames(Bmat)

	
	return(list(Var.prop=stab1,Return.times=c(mean.return=mean.return,var.return=var.return),Reactivity=stab3, Spp.contribs = spp.contribs))
}

plot.marts <- function(tsmat, loc.legend){
	
	tsmat <- exp(tsmat);
	max.abund <- max(tsmat)
	min.abund <- min(tsmat)
	spp.names <- colnames(tsmat)
	len <- nrow(tsmat)
	mycols <- topo.colors(length(spp.names))
	plot(x=1:len,y=tsmat[,1],xlab="Time (days)", ylab="Estimated abundance (qpcr reads)", type="b", lwd=2, col=mycols[1], ylim=c(min.abund,max.abund), main="Total abundances", cex.lab=1.5)
	for(i in 2:length(spp.names)){
		
		points(x=1:len, y=tsmat[,i], type="b", lwd=2, col=mycols[i])
		
	}
	legend(loc.legend, legend=spp.names, lty=1, lwd=2, col=mycols)
	
	
}


#par(mfrow=c(3,4))
#for(i in 1:12){
#	
#	pred.ts <- exp(kalman.mat[i,])
#	plot(days.kalman,pred.ts, type="b", col="blue")
#	points(days1,tsmat1[i,], type="b", col="red")
#}


plot.pmarts <- function(fname, loc.legend){
	
	tsmat <- as.matrix(read.delim(fname));
	tsmat <- exp(tsmat);
	pmat <- t(apply(tsmat,1,FUN=function(x){x/sum(x)}));
	spp.names <- colnames(tsmat)
	len <- nrow(tsmat)
	mycols <- topo.colors(length(spp.names))
	plot(x=1:len,y=pmat[,1],xlab="Time (days)", ylab="Estimated relative abundance", type="b", lwd=2, col=mycols[1], ylim=c(0,1), main="Relative abundances")
	for(i in 2:length(spp.names)){
		
		points(x=1:len, y=pmat[,i], type="b", lwd=2, col=mycols[i])
		
	}
	#legend(loc.legend, legend=spp.names, lty=1, lwd=2, col=mycols)
	
	
}

#This function creates a color scale for use with e.g. the image()
#function. Input parameters should be consistent with those
#used in the corresponding image plot. The "horiz" argument
#defines whether the scale is horizonal(=TRUE) or vertical(=FALSE).
#Depending on the orientation, x- or y-limits may be defined that
#are different from the z-limits and will reduce the range of
#colors displayed.
 
image.scale <- function(z, zlim, col = heat.colors(12),
breaks, horiz=TRUE, ylim=NULL, xlim=NULL, ...){
 if(!missing(breaks)){
  if(length(breaks) != (length(col)+1)){stop("must have one more break than colour")}
 }
 if(missing(breaks) & !missing(zlim)){
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1)) 
 }
 if(missing(breaks) & missing(zlim)){
  zlim <- range(z, na.rm=TRUE)
  zlim[2] <- zlim[2]+c(zlim[2]-zlim[1])*(1E-3)#adds a bit to the range in both directions
  zlim[1] <- zlim[1]-c(zlim[2]-zlim[1])*(1E-3)
  breaks <- seq(zlim[1], zlim[2], length.out=(length(col)+1))
 }
 poly <- vector(mode="list", length(col))
 for(i in seq(poly)){
  poly[[i]] <- c(breaks[i], breaks[i+1], breaks[i+1], breaks[i])
 }
 xaxt <- ifelse(horiz, "s", "n")
 yaxt <- ifelse(horiz, "n", "s")
 if(horiz){YLIM<-c(0,1); XLIM<-range(breaks)}
 if(!horiz){YLIM<-range(breaks); XLIM<-c(0,1)}
 if(missing(xlim)) xlim=XLIM
 if(missing(ylim)) ylim=YLIM
 plot(1,1,t="n",ylim=ylim, xlim=xlim, xaxt=xaxt, yaxt=yaxt,xlab="",ylab="", yaxs="i", ...)  
 for(i in seq(poly)){
  if(horiz){
   polygon(poly[[i]], c(0,0,1,1), col=col[i], border=NA)
  }
  if(!horiz){
   polygon(c(0,0,1,1), poly[[i]], col=col[i], border=NA)
  }
 }
}

# Extracts time series data from a big time series matrix with all spp
# Arguments are filename and an array of species names to extract
# This time series data contains all the women's records
ssp.select <- function(fname,spp.extract){
	
	# Assumes original time series matrix has time in rows and spp names in cols
	alldata <- as.matrix(read.table(file=fname, header=TRUE));
	all.names <- colnames(alldata);
	nspp <- length(spp.extract)
	col.indexs <- rep(0,nspp);
	for(i in 1:nspp){col.indexs[i]<-which(all.names==spp.extract[i], arr.ind=TRUE)}
	out <- alldata[,c(1,2,col.indexs)];
	return(out)
}

# From the time series returned by ssp.select, extract the ts data for one woman
one.woman <- function(allts, woman.ID){
	
	out <- allts[allts[,1]==woman.ID,]
	return(out)
}

# Takes as argument three filenames, one per replicated sample and an array
# of the desired species to analize, and creates
	# 1. A list which is a set of three-dimensional matrices, one for every woman, 
	# Dim1 = species 
	# Dim2 = day
	# Dim3 = rep
	# 2. A vector of days for which there is data, output as the colnames of the ts matrices
all.women <- function(reps.fnames, focal.spp){
	
	rep1tsdata <- ssp.select(fname=reps.fnames[1], spp.extract=focal.spp)
	rep2tsdata <- ssp.select(fname=reps.fnames[2], spp.extract=focal.spp)
	rep3tsdata <- ssp.select(fname=reps.fnames[3], spp.extract=focal.spp)
	wom.ids <- names(table(rep1tsdata[,1]));
	nwomen  <- length(wom.ids);
	nspp    <- length(focal.spp);
	women.list <- list()
	for(i in 1:nwomen){
		
		rep1 <- one.woman(allts=rep1tsdata, woman.ID=wom.ids[i]);
		rep2 <- one.woman(allts=rep2tsdata, woman.ID=wom.ids[i]);
		rep3 <- one.woman(allts=rep3tsdata, woman.ID=wom.ids[i]);
		ndays <- nrow(rep1);
		submat1 <- matrix(as.numeric(rep1[,-c(1,2)]),nrow=nspp,ncol=ndays, byrow=TRUE);
		submat2 <- matrix(as.numeric(rep2[,-c(1,2)]),nrow=nspp,ncol=ndays, byrow=TRUE);		
		submat3 <- matrix(as.numeric(rep3[,-c(1,2)]),nrow=nspp,ncol=ndays, byrow=TRUE);
		submat <- array(0,dim=c(nspp,ndays,3));
		dimnames(submat) <- list(focal.spp, rep1[,2],c("rep1", "rep2", "rep3"))
		submat[,,1] <- submat1;
		submat[,,2] <- submat2;		
		submat[,,3] <- submat3;		
		women.list[[i]] <- submat; 
	}
	names(women.list) <- wom.ids;
	return(women.list)
}


GSSkalman.3reps <- function(ts.mats, K, my.thin, nchains, niter, nadapt){

		has.nas1 <-sum(is.na(ts.mats[1,,1]))
		if(has.nas1>0){		 
			na.ind1 <- which(is.na(ts.mats[1,,1])==TRUE, arr.ind=TRUE);
			tsmat1 <- ts.mats[,-na.ind1,1]; 
		}else{tsmat1 <- ts.mats[,,1]};

		has.nas2 <-sum(is.na(ts.mats[1,,2]))
		if(has.nas2>0){		 
			na.ind2 <- which(is.na(ts.mats[1,,2])==TRUE, arr.ind=TRUE);
			tsmat2 <- ts.mats[,-na.ind2,2]; 
		}else{tsmat2 <- ts.mats[,,2]};
		
		has.nas3 <-sum(is.na(ts.mats[1,,3]))
		if(has.nas3>0){		 
			na.ind3 <- which(is.na(ts.mats[1,,3])==TRUE, arr.ind=TRUE);
			tsmat3 <- ts.mats[,-na.ind3,3];
		}else{tsmat3 <- ts.mats[,,3]};
	
	tsmat1[tsmat1==0] <-.Machine$double.eps;
	tsmat2[tsmat2==0] <-.Machine$double.eps;
	tsmat3[tsmat3==0] <-.Machine$double.eps;
	ltsmat1 <- log(tsmat1);
	ltsmat2 <- log(tsmat2);
	ltsmat3 <- log(tsmat3);

	days1 <- as.numeric(colnames(tsmat1))
	days1 <- days1-min(days1) + 1;
	lendays1 <- length(days1)
	days2 <- as.numeric(colnames(tsmat2))
	days2 <- days2-min(days2) + 1;
	lendays2 <- length(days2)
	days3 <- as.numeric(colnames(tsmat3))
	days3 <- days3-min(days3) + 1;
	lendays3 <- length(days3)
	all.days <- c(days1,days2,days3)
	last.day <- max(all.days)
	days.kalman <- 1:last.day
	ndayskal <- length(days.kalman)
	nspp <- nrow(tsmat1)
	kalman.mat <- matrix(0,nrow=nspp, ncol=ndayskal)

	for(j in 1:nspp){

		Y1 <- matrix(0,nrow=K,ncol=lendays1);
		Y2 <- matrix(0,nrow=K,ncol=lendays2);
		Y3 <- matrix(0,nrow=K,ncol=lendays3);
		for(k in 1:K){
			Y1[k,] <- ltsmat1[j,];
			Y2[k,] <- ltsmat2[j,];
			Y3[k,] <- ltsmat3[j,];						
		}				
		
		gss.data <- list(Y1=Y1,Y2=Y2,Y3=Y3, K=K, len=ndayskal, days1=days1, days2=days2, 
						  days3=days3, lendays1=lendays1,lendays2=lendays2, lendays3=lendays3);
		mardc.model <- jags.model('GSS3repsjags.txt', data=gss.data, n.chains=nchains, n.adapt=nadapt);
		mar2.samples <- coda.samples(model=mardc.model, variable.names=c('a', 'cc','sig','tau'),n.iter=niter, thin=my.thin)
		mar.chain <- as.matrix(mar2.samples)
		mles <- apply(mar.chain,2,mean)
		
		xgy.data <- list(Y1=Y1[1,], Y2=Y2[1,], Y3=Y3[1,], a=mles[1],cc=mles[2], sig=mles[3], tau=mles[4], 
			             len=ndayskal,len1=lendays1,len2=lendays2,len3=lendays3,days1=days1, days2=days2,days3=days3)
		xgy.model <- jags.model('hofXgYGSS3reps.txt', data=xgy.data, n.chains=nchains, n.adapt=nadapt);		
		xgy.samples <- coda.samples(model=xgy.model, variable.names='X1', n.iter=niter,thin=my.thin);	
		xgy.chain <- as.matrix(xgy.samples)
		kalman.preds <- apply(xgy.chain,2,mean)	
		kalman.mat[j,] <- kalman.preds		
	}
	return(kalman.mat)
}
	

# 5. Kalman estimation for the GSS model for a set of time series.
#	comm.mat is the community matrix with species as rows and time as columns
#	K defines the numebr of times for cloning the data to improve the MLEs
#	Full.results stores the full mcmcchain estimations
#	ncahins number of chains to run
#	nadapt number of itertions to sample from
#	niter number of iterations tu run
#	If plot.it=TRUE, then specify number of columns and number of rows in which you want to plot the graphs by
#	specifying cols = Number of Columns, rows = Number of Rows

GSSkalman.1rep	<-	function(comm.mat,K,full.results=FALSE,nchains,nadapt,niter,plot.it=FALSE,cols=1,rows=1){

	len			<-	ncol(comm.mat)
	nspp		<-	nrow(comm.mat)
	MLEs		<-	matrix(NA,nrow=nspp,ncol=3,dimnames=list(rownames(comm.mat),c("a","cc","sig")))	
	Xt.pred		<-	matrix(NA,nrow=nspp,ncol=len,dimnames=list(rownames(comm.mat),1:len))
	mcmc.full	<-	list()
	nchains		<-	nchains
	nadapt		<-	nadapt
	niter		<-	niter
	for (i in 1:nspp){
		tskreps1	<-	matrix(rep(comm.mat[i,],K),nrow=K,ncol=len,byrow=T)
		# Setting up the jags model:
		## DEFINE lendays and days AS ELEMENTS OF THE DATA MATRIX (SWAP len FOR lendays)
		Gmodel <- jags.model('Gompertz-State_Space_Mod.txt', data= list('Y1'= tskreps1, 'K'= K,'len'=len), n.chains = nchains, n.adapt = nadapt)

		#Doing the MCMC
		G.samples <- coda.samples(Gmodel, c('a', 'cc', 'sig1'),n.iter=niter)
		G.mcmcchain <- as.matrix(G.samples)
		MLES <- apply(G.mcmcchain,2,mean)
		VCOV <- K*var(G.mcmcchain)

		#####----- After getting the MLes, sample from h(X|Y), tilting at the mles
		#####----- The model file for the conditional posterior of X|Y is "hofXgYGompertz1rep.txt"

		Y1		<-	as.vector(comm.mat[i,])
		## DEFINE lendays and days AS ELEMENTS OF THE DATA MATRIX (SWAP len FOR lendays)
		hofx.data <- list(len=len,cc=MLES[2], a=MLES[1], sig=MLES[3],Y1=Y1)

		# Setting up the jags model:
		hofxgy.model <- jags.model('hofXgYGompertz.txt', data= hofx.data, n.chains = nchains, n.adapt = nadapt)

		#Doing the MCMC
		hofxgy.samples <- coda.samples(hofxgy.model, c('X1'),niter)
		hofxgy.mcmcchain <- as.matrix(hofxgy.samples)

		# Then, the MLEs of where the process truly was (the Kalman estimates) is given by the mean
		# of the posterior samples of the X's given the Y's:
		Xt.predicted <- apply(hofxgy.mcmcchain, 2, mean)

		MLEs[i,]		<-	MLES
		Xt.pred[i,]		<-	Xt.predicted
		mcmc.full[[i]]	<-	G.mcmcchain
	}
	# Nomas para ver:
	if(plot.it==TRUE){
		
		par(mfrow=c(rows,cols),mar=c(3,3,1,1),mgp=c(2,0.5,0))
		
		for (i in 1:nspp){
		
			plot(0:(len-1), comm.mat[i,], pch=16, xlab="Time", ylab="Counts",main=rownames(comm.mat)[i])
			points(0:(len-1), exp(Xt.pred[i,]), type="l", col="red")
		}
	}
	if(full.results==FALSE){
		results	<-	list(MLEs=MLEs,Xt.pred=Xt.pred)
		}else{results<- list(MLEs=MLEs,Xt.pred=Xt.pred,mcmc.full=mcmc.full)}
	return(results)
}

mars.cls.boot	<-	function(nboot,comm.mat,cls.fit,covariate=NULL){
	
	lagstate	<-	comm.mat[1:(nrow(comm.mat)-1),]
	statevar	<-	comm.mat[2:nrow(comm.mat),]
	
	Q	<-	nrow(statevar)
	P	<-	ncol(statevar)
	bootA	<-	matrix(NA,nrow=P,ncol=nboot)
	bootB	<-	array(NA,dim=c(P,P,nboot))
	if(is.null(covariate)==TRUE){bootC	<-	NULL
		}else{bootC	<-	array(NA,dim=c(P,R,nboot))}
	bootSigma<-	array(NA,dim=c(P,P,nboot))
		
	
	if(is.null(covariate)==TRUE){R	<-	NULL
		}else{R	<-	ncol(covariate)}
	
	for (i in 1:nboot){
		# Randomize quartets of residuals
		rand.num	<-	sample(1:Q,Q,replace=TRUE) # get "new" residuals
		Estar		<-	cls.fit$E[rand.num,]
	
		s1		<-	length(comm.mat[,1])
		Xstar	<-	matrix(0,ncol=P,nrow=s1)
		Xstar[1,]<-	comm.mat[1,]	
	
			for(j in 2:s1){
				if(is.null(covariate)==TRUE){Xstar[j,]	<-	cls.fit$A + cls.fit$B%*%Xstar[j-1,] + Estar[j-1,]
					}else if(R==1){Xstar[j,]	<-	cls.fit$A + cls.fit$B%*%Xstar[j-1,] + cls.fit$C*covariate[j-1] + Estar[j-1,]
					}else if(R>1){Xstar[j,]	<-	cls.fit$A + cls.fit$B%*%Xstar[j-1,] + cls.fit$C%*%covariate[j-1] + Estar[j-1,]}
			}
		
	
		fit	<-	mars.cls(comm.mat=Xstar,covariate=covariate)
		bootA[,i]	<-	fit$A
		bootB[,,i]	<-	fit$B
		bootC[,,i]	<-	fit$C
		bootSigma[,,i]<-	fit$sigma
	}
	A.ci	<-	matrix(NA,ncol=2,nrow=P)
	B.ciL	<-	matrix(NA,ncol=P,nrow=P)
	B.ciU	<-	matrix(NA,ncol=P,nrow=P)
	if(is.null(covariate)==FALSE){C.ciL	<- matrix(NA,ncol=R,nrow=P)
									C.ciU	<- matrix(NA,ncol=R,nrow=P)
	}else{C.ciL=NULL
			C.ciU=NULL}
	Sigma.ciU	<-	matrix(NA,ncol=P,nrow=P)
	Sigma.ciL	<-	matrix(NA,ncol=P,nrow=P)
	
		for (l in 1:P){
			A.ci[l,]	<-	quantile(bootA[l,],probs=c(0.025,0.975))
			for(k in 1:P){
				B.ciL[l,k]	<-	quantile(bootB[l,k,],probs=0.025)
				B.ciU[l,k]	<-	quantile(bootB[l,k,],probs=0.975)
				Sigma.ciL[l,k]<-	quantile(bootSigma[l,k,],probs=0.025)
				Sigma.ciU[l,k]<-	quantile(bootSigma[l,k,],probs=0.975)
			}
			if(is.null(covariate)==FALSE){
			for(s in 1:R){
				C.ciL[l,s]	<-	quantile(bootC[l,s,],probs=0.025)
				C.ciU[l,s]	<-	quantile(bootC[l,s,],probs=0.975)
			}}
		}
		results	<-	list(A=A.ci,B.L=B.ciL,B.U=B.ciU,C.U=C.ciU,C.L=C.ciL,Sigma.U=Sigma.ciU,Sigma.L=Sigma.ciL)
		return(results)
}




