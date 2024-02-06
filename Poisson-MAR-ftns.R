# Functions for the Poisson-MAR model:
library(MASS)
library("rjags")

# 1. rand.MVN:  Multivariate Normal random number generator
#    n = number of random samples of a MVN vector
#    mu = mean vector of the MVN distribution to sample from
#    cov.mat = Variance-covariance matrix of the MVN distribution to sample from
randmvn <- function(n,mu.vec, cov.mat){
	
	p <- length(mu.vec);
	Tau <- chol(cov.mat);
	Zmat <- matrix(rnorm(n=p*n,mean=0,sd=1),nrow=p,ncol=n); #generate normal deviates outside loop
	out <- matrix(0,nrow=p,ncol=n);
	for(i in 1:n){
		
		Z <- Zmat[,i];
		out[,i] <- t(Tau)%*%Z + mu.vec
		
		}
	
	return(out)
	
}

# Returns the population abundance (NOT log)
# for a 2 spp MAR model
Twospp.mar <- function(A,B,Sigma,len,Xo,rnd.init=FALSE){
	
	# if rnd.int==FALSE, provide ann initial Xo value
	Xmat <- matrix(0,ncol=len,nrow=2)
	X.start <- Xo
	if(rnd.init==TRUE){
		Vec.V <- ginv(diag(4) - kronecker(B,B))%*%as.vector(Sigma) # Eq. 17th Ives et al 2003
		V <- matrix(Vec.V, nrow=2,ncol=2,byrow=FALSE)
		#if(V[2,2]<0){V[2,2] <-Sigma[2,2]/(1-B[2,2]^2)}
		muvec <- ginv(diag(2)-B)%*%A
		zero.mean.tst <- muvec<=log(.Machine$double.xmin)
		if(sum(zero.mean.tst)>0){muvec[zero.mean.tst==1]<- log(.Machine$double.xmin)}
		X.start <- randmvn(n=1,mu.vec= muvec , cov.mat=V)

	}
	Xmat[,1] <- X.start
	
	for( i in 2:len){
		
		im1 <- i-1;
		Xim1  <- matrix(Xmat[,im1], nrow=2,ncol=1)
		mui <- A + B%*%Xim1
		rand.trans <- randmvn(n=1,mu.vec= mui, cov.mat=Sigma)
		Xmat[,i] <- rand.trans
	}
	
	return(t(exp(Xmat)))	
}



# 2. A 2 spp. simulation , compositional sim. ftn, Gompertz map & Cobweb, MAR(1) likelihood function  ---------#######
Twospp.mar1 <- function(A,B,Sigma,len, plot.it=FALSE){
	
	Xmat <- matrix(0,ncol=len,nrow=2)
	Nmat <- Xmat	
	Vec.V <- ginv(diag(4) - kronecker(B,B))%*%as.vector(Sigma) # Eq. 17th Ives et al 2003
	V <- matrix(Vec.V, nrow=2,ncol=2,byrow=FALSE)
	Xo <- randmvn(n=1,mu.vec= ginv(diag(2)-B)%*%A, cov.mat=V)
	Xmat[,1] <- Xo
	Nmat[,1] <- rpois(n=2, lambda=exp(Xo))
	
	for( i in 2:len){
		
		im1 <- i-1;
		Xim1  <- matrix(Xmat[,im1], nrow=2,ncol=1)
		mui <- A + B%*%Xim1
		rand.trans <- randmvn(n=1,mu.vec= mui, cov.mat=Sigma)
		Xmat[,i] <- rand.trans
		Nmat[,i] <- rpois(n=2,lambda=exp(rand.trans))
	}
	
	if(plot.it==TRUE){
		plot(1:len, Nmat[1,], pch=16, type="b", col="red", ylim=c(0,max(Nmat)), main ="True Abundances", xlab="Time", ylab="2-species population sizes")
		points(1:len, Nmat[2,],pch=16, type="b", col="blue")
	}
	return(list(Xmat=Xmat, Nmat=Nmat))	
}


# 3. Three species Simulation

Threespp.mar1 <- function(A,B,Sigma,len, plot.it=FALSE){
	
	Xmat <- matrix(0,ncol=len,nrow=3)
	Nmat <- Xmat	
	Vec.V <- ginv(diag(9) - kronecker(B,B))%*%as.vector(Sigma) # Eq. 17th Ives et al 2003
	V <- matrix(Vec.V, nrow=3,ncol=3,byrow=FALSE)
	Xo <- randmvn(n=1,mu.vec= ginv(diag(3)-B)%*%A, cov.mat=V)
	Xmat[,1] <- Xo
	Nmat[,1] <- rpois(n=3, lambda=exp(Xo))
	
	for( i in 2:len){
		
		im1 <- i-1;
		Xim1  <- matrix(Xmat[,im1], nrow=3,ncol=1)
		mui <- A + B%*%Xim1
		rand.trans <- randmvn(n=1,mu.vec= mui, cov.mat=Sigma)
		Xmat[,i] <- rand.trans
		Nmat[,i] <- rpois(n=3,lambda=exp(rand.trans))
	}
	
	if(plot.it==TRUE){
		plot(1:len, Nmat[1,], pch=16, type="b", col="red", ylim=c(0,max(Nmat)), main ="True Abundances", xlab="Time", ylab="3-species population sizes")
		points(1:len, Nmat[2,],pch=16, type="b", col="blue")
		points(1:len, Nmat[3,],pch=16, type="b", col="black")
	}
	return(list(Xmat=Xmat, Nmat=Nmat))	
}

Threespp.marN <- function(A,B,Sigma,tausq,len,nreps, plot.it=FALSE){
	
	tau  <- sqrt(tausq);
	Xmat <- matrix(0,ncol=len,nrow=3)
	Nmat <- array(0,dim=c(3,len,nreps))	
	Vec.V <- ginv(diag(9) - kronecker(B,B))%*%as.vector(Sigma) # Eq. 17th Ives et al 2003
	V <- matrix(Vec.V, nrow=3,ncol=3,byrow=FALSE)
	Xo <- randmvn(n=1,mu.vec= ginv(diag(3)-B)%*%A, cov.mat=V)
	Xmat[,1] <- Xo

	for(n in 1:nreps){

		Nmat[,1,n] <- rnorm(n=3, mean=exp(Xo), sd=tau)		
	}

	
	for( i in 2:len){
		
		im1 <- i-1;
		Xim1  <- matrix(Xmat[,im1], nrow=3,ncol=1)
		mui <- A + B%*%Xim1
		rand.trans <- randmvn(n=1,mu.vec= mui, cov.mat=Sigma)
		Xmat[,i] <- rand.trans
		for(n in 1:nreps){	
			Nmat[,i,n] <- rnorm(n=3,mean=exp(rand.trans), sd=tau)
		}
	}
	
	if(plot.it==TRUE){
	
		par(mfrow=c(1,nreps))
		for(n in 1:nreps){	
		
			plot(1:len, Nmat[1,,n], pch=16, type="b", col="red", ylim=c(0,max(Nmat)), main ="True Abundances", xlab="Time", ylab="3-species population sizes")
			points(1:len, Nmat[2,,n],pch=16, type="b", col="blue")
			points(1:len, Nmat[3,,n],pch=16, type="b", col="black")
		}
	
	}
	return(list(Xmat=Xmat, Nmat=Nmat))	
}


Threespp.marZLN <- function(A,B,Sigma,tausq,ndetect, len,nreps, plot.it=FALSE){
	
	tau  <- sqrt(tausq);
	Xmat <- matrix(0,ncol=len,nrow=3)
	Nmat <- array(0,dim=c(3,len,nreps))	
	Vec.V <- ginv(diag(9) - kronecker(B,B))%*%as.vector(Sigma) # Eq. 17th Ives et al 2003
	V <- matrix(Vec.V, nrow=3,ncol=3,byrow=FALSE)
	Xo <- randmvn(n=1,mu.vec= ginv(diag(3)-B)%*%A, cov.mat=V)
	Xmat[,1] <- Xo
	if(length(ndetect)==1){ndetect <- rep(ndetect,3)}
	
	for(n in 1:nreps){

		Nmat[,1,n] <- rnorm(n=3, mean=exp(Xo), sd=tau)		
	}

	
	for( i in 2:len){
		
		im1 <- i-1;
		Xim1  <- matrix(Xmat[,im1], nrow=3,ncol=1)
		mui <- A + B%*%Xim1
		rand.trans <- randmvn(n=1,mu.vec= mui, cov.mat=Sigma)
		Xmat[,i] <- rand.trans
		for(n in 1:nreps){	
			simn <- exp(rnorm(n=3,mean=rand.trans, sd=tau))
			is.detected <- simn>ndetect
			Nmat[,i,n] <- is.detected*simn
		}
	}
	
	if(plot.it==TRUE){
	
		par(mfrow=c(1,nreps))
		for(n in 1:nreps){	
		
			plot(1:len, Nmat[1,,n], pch=16, type="b", col="red", ylim=c(0,max(Nmat)), main ="True Abundances", xlab="Time", ylab="3-species population sizes")
			points(1:len, Nmat[2,,n],pch=16, type="b", col="blue")
			points(1:len, Nmat[3,,n],pch=16, type="b", col="black")
		}
	
	}
	return(list(Xmat=Xmat, Nmat=Nmat))	
}

Pspp.marZLN <- function(A,B,Sigma,tausq,ndetect, len,nreps, plot.it=FALSE){
	
	p    <- length(A);
	tau  <- sqrt(tausq);
	Xmat <- matrix(0,ncol=len,nrow=p)
	Nmat <- array(0,dim=c(p,len,nreps))	
	Vec.V <- ginv(diag(p*p) - kronecker(B,B))%*%as.vector(Sigma) # Eq. 17th Ives et al 2003
	V <- matrix(Vec.V, nrow=p,ncol=p,byrow=FALSE)
	Xo <- randmvn(n=1,mu.vec= ginv(diag(p)-B)%*%A, cov.mat=V)
	Xmat[,1] <- Xo
	if(length(ndetect)==1){ndetect <- rep(ndetect,p)}
	
	for(n in 1:nreps){

		Nmat[,1,n] <- rnorm(n=p, mean=exp(Xo), sd=tau)		
	}

	
	for( i in 2:len){
		
		im1 <- i-1;
		Xim1  <- matrix(Xmat[,im1], nrow=p,ncol=1)
		mui <- A + B%*%Xim1
		rand.trans <- randmvn(n=1,mu.vec= mui, cov.mat=Sigma)
		Xmat[,i] <- rand.trans
		for(n in 1:nreps){	
			simn <- exp(rnorm(n=p,mean=rand.trans, sd=tau))
			is.detected <- simn>ndetect
			Nmat[,i,n] <- is.detected*simn
		}
	}
	
	if(plot.it==TRUE){
	
		par(mfrow=c(1,nreps))
		for(n in 1:nreps){	
		
			plot(1:len, Nmat[1,,n], pch=16, type="b", col="red", ylim=c(0,max(Nmat)), main =paste0("Replicate ", n), xlab="Time", ylab=paste0(p," species population sizes"))
			points(1:len, Nmat[2,,n],pch=16, type="b", col="blue")
			points(1:len, Nmat[3,,n],pch=16, type="b", col="black")
		}
	
	}
	return(list(Xmat=Xmat, Nmat=Nmat))	
}



# 4. A rjags Poisson-MARS estimation function:

Pois.mars <- function(ts.mat, K, my.thin, nchains, niter, nadapt, MARS.model = "Full"){
	
	Nsim1 <- ts.mat;
	len <- length(Nsim1[1,])
	nspp <- length(Nsim1[,1])
	mu0 <- rep(0,nspp)
	for(i in 1:nspp){mu0[i] <- mean(log(Nsim1[i,Nsim1[i,]>0]))}
	S0 <- diag(length(mu0))
	S1 <- diag(length(mu0))

	# Cloned data set!
	Nsim1K <- array(0,dim=c(nspp,len,K))
	for(i in 1:K){
		Nsim1K[, ,i] <- Nsim1;	
	}

	if(MARS.model=="Full"){
		
		mardc.data <- list(N = Nsim1K, K = K, mu0=mu0, S0=S0, S1=S1, len=len, nspp=nspp)
		
		# Setting up the jags model
		mardc.model <- jags.model('MARSDC.txt', data=mardc.data, n.chains=nchains, n.adapt=nadapt);
		# Sampling from the MCMC 
		mar1.samples <- coda.samples(model=mardc.model, variable.names=c('A', 'B','Sigma'),n.iter=niter, thin=my.thin)
		mar.chain <- as.matrix(mar1.samples)
		mles <- apply(mar.chain,2,mean)

		
	}else if(MARS.model=="Diag. Sigma"){
		
		
		mardc.data <- list(N = Nsim1K, K = K, mu0=mu0, S0=S0, len=len, nspp=nspp)
		
		# Setting up the jags model
		mardc2.model <- jags.model('MARSDC-simple.txt', data=mardc.data, n.chains=nchains, n.adapt=nadapt);

		# Sampling from the MCMC 
		mar2.samples <- coda.samples(model=mardc2.model, variable.names=c('A', 'B','sig'),n.iter=niter, thin=my.thin)
		mar.chain <- as.matrix(mar2.samples)
		mles <- apply(mar.chain,2,mean)
		
	}

	out <- list(outchain = mar.chain, mles=mles)
	
	return(out)
	
}


Pois.mars.nreps <- function(ts.mat, K, my.thin, nchains, niter, nadapt, nreps,days,totdays, mars.model="Full"){

	#  totdays is the real length of the time series, some time steps are unobserved

	Nsim1 <- ts.mat;
	lendays <- length(Nsim1[1,,1]) # = number of days for which we have data
	nspp <- length(Nsim1[,1,1])
	mu0 <- rep(0,nspp)
	for(i in 1:nspp){mu0[i] <- mean(log(Nsim1[i,Nsim1[i,,1]>0,1]))}
	S0 <- diag(length(mu0))
	S1 <- diag(length(mu0))

	# Checking for 0's
	for(i in 1:nspp){
		for(j in 1:lendays){
			for(ell in 1:nreps){
				
				if(Nsim1[i,j,ell]==0){Nsim1[i,j,ell] <- .Machine$double.eps}
				
			}
		}
		
	}
	# Cloned data set!
	Nsim1K <- array(0,dim=c(nspp,lendays,nreps,K))
	for(i in 1:K){
		Nsim1K[, , ,i] <- Nsim1;	
	}
	YK <- log(Nsim1K)
	
	if(mars.model=="Full"){
		# Data for jags
		mardcn.data <- list(Y = YK, K = K, mu0=mu0, S0=S0, S1=S1, len=totdays, nspp=nspp, nreps=nreps, lendays=lendays, days=days)
	
		# Setting up the jags model
		mardcn.model <- jags.model('MARSDCnreps.txt', data=mardcn.data, n.chains=nchains, n.adapt=nadapt);
	
		# Sampling from the MCMC 
		mar.samples <- coda.samples(model=mardcn.model, variable.names=c('A', 'B','Sigma','tau'),n.iter=niter, thin=my.thin)
	}else if(mars.model=="Diag. Sigma"){
		# Data for jags
		mardcn.data <- list(Y = YK, K = K, mu0=mu0, S0=S0, len=totdays, nspp=nspp, nreps=nreps, lendays=lendays, days=days)
	
		# Setting up the jags model
		mardcn.model <- jags.model('MARSDC-simple-nreps.txt', data=mardcn.data, n.chains=nchains, n.adapt=nadapt);
	
		# Sampling from the MCMC 
		mar.samples <- coda.samples(model=mardcn.model, variable.names=c('A', 'B','sig','tau'),n.iter=niter, thin=my.thin)
		
	}
	
	mar.chain <- as.matrix(mar.samples)
	mles <- apply(mar.chain,2,mean)
	fishinv <- K*var(mar.chain)
	out <- list(outchain = mar.chain, mles=mles, fishinv)
	
	return(out)
	
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
		Gmodel <- jags.model('Gompertz-State_Space_Mod.txt', data= list('Y1'= tskreps1, 'K'= K,'len'=len), n.chains = nchains, n.adapt = nadapt)

		#Doing the MCMC
		G.samples <- coda.samples(Gmodel, c('a', 'cc', 'sig1'),n.iter=niter)
		G.mcmcchain <- as.matrix(G.samples)
		MLES <- apply(G.mcmcchain,2,mean)
		VCOV <- K*var(G.mcmcchain)

		#####----- After getting the MLes, sample from h(X|Y), tilting at the mles
		#####----- The model file for the conditional posterior of X|Y is "hofXgYGompertz1rep.txt"

		Y1		<-	as.vector(comm.mat[i,])
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
		results	<-	list(MLEs,Xt.pred)
		}else{results<- list(MLEs,Xt.pred,mcmc.full)}
	return(results)
}

# 6. Probability density function for the Multivariate Normal Distribution 

dmvnorm	<-	function(x,mu,sigma){
	p	<-	length(x)
	log.lik	<- (-p/2)*log(2*pi)-(1/2)*log(prod(diag(sigma)))-(1/2)*(t(x-mu)%*%ginv(sigma)%*%(x-mu))
	
	return(log.lik)	
	
}



