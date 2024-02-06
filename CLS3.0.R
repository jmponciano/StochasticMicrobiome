##### Program for CLS estimation of the MARSS(1) model
# statevar is the time series data species as columns and time as rows
# covariates are the covariates for the model
library(MASS)


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

	results	<-	list(A=A,B=B,sigma=sigma,C=C,E=E,Yhat=Yhat,R2=R2,R2_D=R2_D,AIC=lsAIC,BIC=lsBIC, lnlike=lnlike)
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


B.colors <- function(B){
	
	my.cols <- 1:4#grey.colors(n=4)
	wh.col1 <- which(B <= -0.5)
	wh.col2 <- which((-0.5<=B)&(B<=0))
	wh.col4 <- which((0<=B)&(B<=0.5))
	wh.col4 <- which((0.5<=B)&(B<=1))
	wh.col5 <- which(B>=1)
	
	B[wh.col1] <- my.cols[1]
	B[wh.col2] <- my.cols[2]
	B[wh.col3] <- my.cols[3]		
	B[wh.col4] <- my.cols[4]	
	B[wh.col5] <- my.cols[5]	
		
	return(B)

}

mycols.ftn <- function(n){
	d <- 360/n
	h <- cumsum(c(15,rep(d,(n-1))))
	return(hcl(h=h, c=100,l=65))
}

unique.cols <- function(mycols,local.sppnames, all.sppnames){
	
	# mycols length = all.sppnames' length
	
	cols.inds <- match(local.sppnames,all.sppnames)
	return(mycols[cols.inds])
	
}



