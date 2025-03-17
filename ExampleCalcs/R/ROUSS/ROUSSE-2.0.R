######  R Appendix accompanying the manuscript:  
######  Brian Dennis and Jose Miguel Ponciano. 2013. Density dependent state space model for population abundance data with unequal time intervals.
######  To be submitted to Ecology

######  PART I)  OUSS model R functions
######  PART II) EGSS model R functions
######  PART III) Density dependence vs. Density independence Parametric Bootstrap Likelihood Ratio Test (PBLRT) functions 

######  All EGSS Functions originally published in 
######  Humbert, J.-Y., Mills, S., Horne, J. S. and Dennis, B. 2009. A better way to estimate population trend. Oikos 118: 1940–1946.


correct.time <- function(Tvec){
	
	are.there.reps <- sum(table(Tvec)>1)
	len.Tvec <- length(Tvec)
	
	if(are.there.reps>0){
		
		where.reps <- which(table(Tvec)>1, arr.ind=TRUE)
		rep.nums <- as.numeric(names(table(Tvec))[where.reps]) # vec. of all numbers that are repeated
		nreps    <- table(Tvec)[where.reps] # vec. of how many times each of the repeated nums. are repeated
		len.reps <- length(nreps)
		where.inTvec <- match(x=rep.nums, table=Tvec) # only position of first
	
		for(h in 1:len.reps){
			
			hth.num   <- rep.nums[h] #the offending repeating number
			hth.nreps <- nreps[h] # how many times that number is repeated
			hth.ind   <- where.inTvec[h] # position in the original vector of offending number
			next.left <- Tvec[(hth.ind-1)] # number immediately to the left of offending num.
			
			### Idea of what's next: substitute sequence of repeating numbers by 
			### numbers differing only by one.
			absdiff.left <- abs(hth.num-next.left) # difference with num. to the left
			next.right.ind <- (hth.ind+hth.nreps) # index of the number to the right, after reps.
			if(next.right.ind>len.Tvec){
				next.right <- next.left+1
				absdiff.right <- absdiff.left+1
				}else{	
					next.right<- Tvec[next.right.ind]
					absdiff.right<- abs(hth.num-next.right)				
				}
			

			if(absdiff.right>absdiff.left){
				sum.quant <- 0:(hth.nreps-1)
			}else if(absdiff.right<absdiff.left){
				sum.quant <- -(hth.nreps-1):0
			}else{sum.quant <- 0:(hth.nreps-1)}
			new.ts    <- rep(hth.num,hth.nreps) + sum.quant

			Tvec[hth.ind:(hth.ind+hth.nreps-1)] <- new.ts			
			
		}		
		
		
	}
	return(Tvec)
	
}


## Left or Right, to which side is the difference of more than 1


#----------------- PART I):  The OUSS functions --------------------------------# 
# 1. Negative log-likelihood score, for ML estimation or Model Selection (eq. 11)

#  ML objective function "negloglike.OU.ml" is negative of log-likelihood;
#  The function uses a multivariate normal log-
#  likelihood (Eq. XX).  The three function arguments are:
#  fguess, vector of parameters (transformed to the real line), 
#  yt, vector of time series of log observed abundances (cannot have 0's ),
#  tt, vector of observation times.
negloglike.OU.ml <- function(fguess,yt,tt){
	
	mu <- fguess[1];
	guess  <- exp(fguess[2:4]); # Constrains parameters > 0
	theta  <- guess[1]; 
	betasq <- guess[2];
	tausq  <- guess[3];
	q      <- length(yt) - 1;
	qp1    <- q+1;
	Var.inf<- betasq/(2*theta); # Better than Dr. Seuss' thing 1 and thing 2!
	ss <- tt[2:qp1] - tt[1:q];
	t.cols     <- matrix(rep(tt,each=qp1),nrow=qp1,ncol=qp1, byrow=FALSE);
	t.rows     <- t(t.cols);
	abs.diffs  <- abs(t.rows-t.cols);
	V <- Var.inf*exp(-theta*abs.diffs);
	diag(V)  <- diag(V) + rep(tausq,qp1);
	mu.vec <- rep(mu,qp1);
	
	neglogl <- (qp1/2)*log(2*pi)+(1/2)*log(det(V))+ (1/2)*(yt-mu.vec)%*%ginv(V)%*%(yt-mu.vec);
	
	return(neglogl)
	
	}

#  2. Negative log-likelihood for REML estimation or Model Selection (eq. 14)
#  REML objective function "negloglike.OU.reml" is negative of log-likelihood
#  for first differences of the log-scale observations.  The REML objective
#  function uses equations ??.  The three
#  function arguments are: 
#  phi, vector of parameters (transformed to the real line), 
#  yt, vector of time series observations (log scale),
#  tt, vector of observation times.
#  The function performs the differencing.
negloglike.OU.reml=function(phi,yt,tt)
{
   theta  <- exp(phi[1]);          #  Constrains th > 0.
   betasq <- exp(phi[2]);      #  Constrains betasq > 0. 
   tausq  <- exp(phi[3]);       #  Constrains tausq > 0.
   Var.inf  <- betasq/(2*theta);       #  Recurring quantity.
   q   <- length(yt)-1;
   qp1 <- q+1;
   ss     <- tt[2:qp1]-tt[1:q];         #  Time intervals.

	t.cols     <- matrix(rep(tt,each=qp1),nrow=qp1,ncol=qp1, byrow=FALSE);
	t.rows     <- t(t.cols);
	abs.diffs  <- abs(t.rows-t.cols);
	Sigma.mat <- Var.inf*exp(-theta*abs.diffs);
   Itausq <- matrix(0,qp1,qp1);
   diag(Itausq) <- rep(tausq,qp1);
   V  <- Sigma.mat+Itausq;
   Dmat <- cbind(-diag(1,q),matrix(0,q,1)) + cbind(matrix(0,q,1),diag(1,q));    #  Differencing matrix.
   Phi.mat <- Dmat%*%V%*%t(Dmat);           #  REML var-cov mattrix.
   wt <- yt[2:qp1]-yt[1:q];
   ofn <- (q/2)*log(2*pi)+0.5*log(det(Phi.mat)) + 0.5*wt%*%ginv(Phi.mat)%*%wt;
   
   if(is.infinite(ofn)==TRUE){return(50000)}else{
   
    return(ofn)}
}

# 3. rand.MVN:  Multivariate Normal random number generator
#    n = number of random samples of a MVN vector
#    mu = mean vector of the MVN distribution to sample from
#    cov.mat = Variance-covariance matrix of the MVN distribution to sample from
randmvn <- function(n,mu.vec, cov.mat){
	
	p <- length(mu.vec);
	Tau <- chol(cov.mat, pivot=TRUE);
	Zmat <- matrix(rnorm(n=p*n,mean=0,sd=1),nrow=p,ncol=n); #generate normal deviates outside loop
	out <- matrix(0,nrow=p,ncol=n);
	for(i in 1:n){
		
		Z <- Zmat[,i];
		out[,i] <- t(Tau)%*%Z + mu.vec
		
		}
	
	return(out)
	
}


# 4. Simulation function: need to provide parameter values and a vector of observation times
#    The multivariate Normal model (see eq. 11) is used to simulate the data
#    nsims = number of bootstrap replicates to simulate
#    parms = vector of parameter values = c(mu_hat, theta_hat, betasq_hat,tausq_hat), where 'hat' denotes
#            either the ML or the REML estimates.
#    Tvec = vector of ORIGINAL observation times (t_0, t_1, t_2, ..., t_q)
ROUSS.sim <- function(nsims,parms,Tvec){

	tt <- Tvec-Tvec[1]
	mu <- parms[1];
	theta  <- parms[2]; 
	betasq <- parms[3];
	tausq  <- parms[4];
	q      <- length(tt) - 1;
	qp1    <- q+1;
	Var.inf<- betasq/(2*theta); 
	ss <- tt[2:qp1] - tt[1:q];
	t.cols     <- matrix(rep(tt,each=qp1),nrow=qp1,ncol=qp1, byrow=FALSE);
	t.rows     <- t(t.cols);
	abs.diffs  <- abs(t.rows-t.cols);
	V <- Var.inf*exp(-theta*abs.diffs);
	diag(V)  <- diag(V) + rep(tausq,qp1);
	m.vec <- rep(mu,qp1);
	out <- randmvn(n=nsims,mu.vec=m.vec, cov.mat = V)
	return(out)
}

# 5. Computing rough initial guesses for ML estimation
#     Yobs = log(Observed time series of abundances)
#     Tvec = vector of sampling times

guess.calc <- function(Yobs,Tvec){
	
	T.t <-Tvec-Tvec[1]; #  For calculations, time starts at zero.
	q <- length(Yobs)-1;      #  Number of time series transitions, q.
	qp1 <- q+1;              #  q+1 gets used a lot, too.
	S.t <- T.t[2:qp1]-T.t[1:q];  #  Time intervals.
	Ybar <- mean(Yobs);
	Yvar <- sum((Yobs-Ybar)*(Yobs-Ybar))/q;
	mu1 <- Ybar;

	# Kludge an initial value for theta based on mean of Y(t+s) given Y(t).
	th1<- -mean(log(abs((Yobs[2:qp1]-mu1)/(Yobs[1:q]-mu1)))/S.t);            
	bsq1<- 2*th1*Yvar/(1+2*th1);         # Moment estimate using stationary
	tsq1<- bsq1;                         #   variance, with betasq=tausq.

	#three 0's 
	three0s <- sum(c(th1,bsq1,tsq1))
	if(three0s==0|is.na(three0s)){th1 <- 0.5;bsq1 <- 0.09; tsq1 <- 0.23;}

	
	out1 <- c(th1,bsq1,tsq1);
	if(sum(out1<1e-7)>=1){out1 <- c(0.5,0.09,0.23)}
	out <- c(mu1,out1);
	return(abs(out))
	
}

# 6.  Computing the ML estimates of the OUSS model, maximized log-likelihood and the AIC score 
#     Yobs = log(Observed time series of abundances)
#     Tvec = vector of sampling times
#     parms.guess = vector of initial guesses for the parameters c(mu_hat, theta_hat, betasq_hat,tausq_hat)
# 
ROUSS.ML <- function(Yobs,Tvec, parms.guess){
	
	tt  <- Tvec-Tvec[1];
	q   <- length(tt) -1; 
	qp1 <- q+1;
	guess.optim <- c(parms.guess[1], log(parms.guess[2:4])) 
	optim.out <- optim(par=guess.optim,fn=negloglike.OU.ml,method="Nelder-Mead",yt=Yobs,tt=tt)
	mles <- c(optim.out$par[1], exp(optim.out$par[2:4]))
	lnL.hat <- - optim.out$value[1]
	AIC     <- -2*lnL.hat + 2*4 #where 4 = length(mles)... 
	
	out <- list(mles=mles, lnL.hat = lnL.hat, AIC=AIC)
	return(out)
}

# 7.  Computing the REML estimates of the OUSS model, maximized log-likelihood and the AIC score 
#     Yobs = log(Observed time series of abundances)
#     Tvec = vector of sampling times
#     parms.guess = vector of initial guesses for the parameters c(mu_hat, theta_hat, betasq_hat,tausq_hat)
# 
ROUSS.REML <- function(Yobs, Tvec, parms.guess){

	tt  <- Tvec-Tvec[1];
	q   <- length(tt) -1; 
	qp1 <- q+1;
	ss  <- tt[2:qp1]-tt[1:q]; 
	guess.optim <- log(parms.guess[2:4]);
	optim.out <- optim(par = guess.optim, fn=negloglike.OU.reml,method="Nelder-Mead", yt=Yobs,tt=tt);
	remls <- exp(optim.out$par);
	lnL.hat <- -optim.out$value[1];
	theta.reml <- remls[1];
	betasq.reml <- remls[2];
	tausq.reml <- remls[3];
	Var.inf <- betasq.reml/(2*theta.reml)
	vx=matrix(1,qp1,qp1);
	for (ti in 1:q){
   		vx[(ti+1):qp1,ti]=exp(-theta.reml*cumsum(ss[ti:q]));
   		vx[ti,(ti+1):qp1]=vx[(ti+1):qp1,ti];
	}
	Sigma.mat=vx*Var.inf;
	Itausq=matrix(0,qp1,qp1);
	diag(Itausq)=rep(tausq.reml,qp1);
	V.reml=Sigma.mat+Itausq;
	j=matrix(1,qp1,1);
	Vinv=ginv(V.reml);
	mu.reml=(t(j)%*%Vinv%*%Yobs)/(t(j)%*%Vinv%*%j);  #  REML estimate of mu.

	out <- list(remls = c(mu.reml, theta.reml, betasq.reml,tausq.reml), lnLhat = lnL.hat)
	return(out)	
}


# 8. Parameteric Bootstrap function

# B = number of bootstrap replicates
# parms = c(mu.mle,theta.mle,betasq.mle,tausq.mle); ML estimates of the parms. with the original time series
# tt = vector of observation times (t_0, t_1, t_2, ..., t_q) from the original time series
# if REML="TRUE", then the parametric bootsrap is computed for REML estimation.  In that case, 'parms' has to
# contain the REML estimates for the original time series

ROUSS.pboot <- function(B=2,parms,Yobs, Tvec, REML="FALSE", plot.it="FALSE"){
	
	tt <- Tvec-Tvec[1]
	long.t <- tt[1]:max(tt)
	nparms    <- length(parms);
	preds.boot1<- matrix(0,nrow=B,ncol=length(tt))
	preds.boot2<- matrix(0,nrow=B,ncol=length(long.t))
	
	if(REML=="TRUE"){

		boot.remles <- matrix(0,nrow=B,ncol=nparms+1); 
		all.sims  <- ROUSS.sim(nsims=B,parms=parms,Tvec=Tvec);
		all.preds <- ROUSS.predict(parms=parms, Yobs=Yobs,Tvec=Tvec,plot.it="FALSE")
		reml.preds<- all.preds[[1]][,2]
		reml.longpreds <- all.preds[[2]][,2]		

		for(b in 1:B ){
		
			bth.timeseries <- all.sims[,b];
			remles.out <- ROUSS.REML(Yobs=bth.timeseries,Tvec=Tvec, parms.guess=parms);
			boot.remles[b,] <- c(remles.out$remls, remles.out$lnLhat);
			#all.bootpreds <- ROUSS.predict(parms=remles.out$remls, Yobs=bth.timeseries,Tvec=Tvec,plot.it="FALSE");			
			# Bootstrap prediction conditioned on observed time series: 
			all.bootpreds <- ROUSS.predict(parms=remles.out$remls, Yobs=Yobs,Tvec=Tvec,plot.it="FALSE");
			preds.boot1[b,] <- all.bootpreds[[1]][,2]
			preds.boot2[b,] <- all.bootpreds[[2]][,2]
			
		} 
	
		CIs.mat <- apply(boot.remles,2,FUN=function(x){quantile(x,probs=c(0.025,0.975))});
		CIs.mat <- rbind(CIs.mat[1,1:4],parms,CIs.mat[2,1:4]);
		rownames(CIs.mat) <- c("2.5%","REMLE","97.5%");
		colnames(CIs.mat) <- c("mu", "theta","betasq","tausq");
		
		#preds.CIs1 <- apply(preds.boot1,2,FUN=function(x){quantile(x,probs=c(0.025,0.975))});
		#mean.boots <- apply(preds.boot1,2,FUN=function(x){quantile(x,probs=0.50)})
		#preds.CIs1 <- t(rbind(Tvec,reml.preds-(mean.boots-preds.CIs1[1,]), reml.preds, reml.preds+(preds.CIs1[2,]-mean.boots)));
		
		#preds.CIs1 <- apply(preds.boot1,2,FUN=function(x){quantile(x,probs=c(0.025,0.5,0.975))});
		#preds.CIs1 <- t(rbind(Tvec, preds.CIs1))
		preds.se <- apply(preds.boot1,2,FUN=function(x){sqrt(var(x))}) 
		preds.CIs1 <- cbind(Tvec, reml.preds-1.95*preds.se, reml.preds, reml.preds+1.95*preds.se)
		preds.CIs1[preds.CIs1[,2]<0, 2] <- 0.01 
		colnames(preds.CIs1) <- c("Year","2.5%","REMLE","97.5%");

		#preds.CIs2 <- apply(preds.boot2,2,FUN=function(x){quantile(x,probs=c(0.025,0.975))});
		#mean.boots2 <- apply(preds.boot2,2,FUN=function(x){quantile(x,probs=0.50)})
		#preds.CIs2 <- t(rbind(long.t+Tvec[1],reml.longpreds-(mean.boots2-preds.CIs2[1,]), reml.longpreds, reml.longpreds+(preds.CIs2[2,]-mean.boots2)));

		#preds.CIs2 <- apply(preds.boot2,2,FUN=function(x){quantile(x,probs=c(0.025,0.5,0.975))});
		#preds.CIs2 <- t(rbind(long.t+Tvec[1],preds.CIs2));
		preds.se2 <- apply(preds.boot2,2,FUN=function(x){sqrt(var(x))}) 
		preds.CIs2 <- cbind(long.t+Tvec[1], reml.longpreds-1.95*preds.se2, reml.longpreds, reml.longpreds+1.95*preds.se2)
		preds.CIs2[preds.CIs2[,2]<0,2] <- 0.01
		colnames(preds.CIs2) <- c("Year","2.5%","REMLE","97.5%");
		#pred.CIs2 <- cbind(preds.CIs2[,1], reml.longpreds-(preds.CIs2[,3]-preds.CIs2[,2]),reml.longpreds,reml.longpreds+(preds.CIs2[,4]-preds.CIs2[,3]))
		
		boot.list <- list(boot.remles = boot.remles, CIs.mat = CIs.mat, preds.CIs1 = preds.CIs1, 
					 preds.CIs2=preds.CIs2);
		
		if(plot.it=="TRUE"){
			
			X11();
			par(mfrow=c(2,2));
			hist(boot.remles[,1],main=expression(hat(mu)),xlab="");abline(v=parms[1],lwd=2,col="red");
			hist(boot.remles[,2],main=expression(hat(theta)),xlab="");abline(v=parms[2],lwd=2,col="red");
			hist(boot.remles[,3],main=expression(hat(beta^2)),xlab="");abline(v=parms[3],lwd=2,col="red");
			hist(boot.remles[,4],main=expression(hat(tau^2)),xlab="");abline(v=parms[4],lwd=2,col="red");				
		}
		return(boot.list)
		
		}else{

		boot.mles <- matrix(0,nrow=B,ncol=nparms+2);
		all.sims  <- ROUSS.sim(nsims=B,parms=parms,Tvec=Tvec);
		all.preds <- ROUSS.predict(parms=parms, Yobs=Yobs,Tvec=Tvec,plot.it="FALSE")
		ml.preds<- all.preds[[1]][,2]
		ml.longpreds <- all.preds[[2]][,2]		
		
		for(b in 1:B ){
		
			bth.timeseries <- all.sims[,b];
			mles.out <- ROUSS.ML(Yobs= bth.timeseries,Tvec=Tvec, parms.guess=parms);
			boot.mles[b,] <- c(mles.out$mles, mles.out$lnL.hat,mles.out$AIC);
			all.bootpreds <- ROUSS.predict(parms=mles.out$mles, Yobs=Yobs,Tvec=Tvec,plot.it="FALSE");
			preds.boot1[b,] <- all.bootpreds[[1]][,2]
			preds.boot2[b,] <- all.bootpreds[[2]][,2]

			} 
		
		CIs.mat <- apply(boot.mles,2,FUN=function(x){quantile(x,probs=c(0.025,0.975))});
		CIs.mat <- rbind(CIs.mat[1,1:4],parms,CIs.mat[2,1:4]);
		rownames(CIs.mat) <- c("2.5%","MLE","97.5%");
		colnames(CIs.mat) <- c("mu", "theta","betasq","tausq");
		
		#preds.CIs1 <- apply(preds.boot1,2,FUN=function(x){quantile(x,probs=c(0.025,0.975))});
		#preds.CIs1 <- t(rbind(Tvec,preds.CIs1[1,], ml.preds, preds.CIs1[2,]));
		preds.CIs1 <- apply(preds.boot1,2,FUN=function(x){quantile(x,probs=c(0.025,0.5,0.975))});
		preds.CIs1 <- t(rbind(Tvec,preds.CIs1));
		colnames(preds.CIs1) <- c("Year","2.5%","MLE","97.5%");
		preds.CIs1 <- cbind(preds.CIs1[,1], ml.preds-(preds.CIs1[,3]-preds.CIs1[,2]),ml.preds,ml.preds+(preds.CIs1[,4]-preds.CIs1[,3]))
		preds.CIs1[pred.CIs1[,2]<0, 2] <- 0.01 

		#preds.CIs2 <- apply(preds.boot2,2,FUN=function(x){quantile(x,probs=c(0.025,0.975))});
		#preds.CIs2 <- t(rbind(Tvec,preds.CIs2[1,], ml.preds, preds.CIs2[2,]));
		preds.CIs2 <- apply(preds.boot2,2,FUN=function(x){quantile(x,probs=c(0.025,0.5,0.975))});
		preds.CIs2 <- t(rbind(long.t+Tvec[1],preds.CIs2));
		preds.CIs2 <- cbind(preds.CIs2[,1], ml.longpreds-(preds.CIs2[,3]-preds.CIs2[,2]),ml.longpreds,ml.longpreds+(preds.CIs2[,4]-preds.CIs2[,3]))
		preds.CIs2[pred.CIs2[,2]<0,2] <- 0.01
		
		colnames(preds.CIs2) <- c("Year","2.5%","MLE","97.5%");
		
		boot.list <- list(boot.mles = boot.mles, CIs.mat = CIs.mat, preds.CIs1 = preds.CIs1
						, preds.CIs2 = preds.CIs2)

		if(plot.it=="TRUE"){
			
			X11();
			par(mfrow=c(2,2));
			hist(boot.mles[,1],main=expression(hat(mu)),xlab="");abline(v=parms[1],lwd=2,col="red");
			hist(boot.mles[,2],main=expression(hat(theta)),xlab="");abline(v=parms[2],lwd=2,col="red");
			hist(boot.mles[,3],main=expression(hat(beta^2)),xlab="");abline(v=parms[3],lwd=2,col="red");
			hist(boot.mles[,4],main=expression(hat(tau^2)),xlab="");abline(v=parms[4],lwd=2,col="red");				
		}
		
		return(boot.list)
		
	}# End if/else
		
}

#ROUSS.pboot(B=20,parms=linx.remles,Yobs=log(Observed.t), Tvec=Time.t, REML="TRUE", plot.it="FALSE")


# 9. Function to compute the predicted trajectory of the unobserved process
#    The arguments are:
#    parms = ML or REML estimates of c(mu,theta,betasq,tausq), which ever you prefer
#    Yobs  = Log observations
#    Tvec  = vector of original observation times (t_0,t_1,...,t_q)
ROUSS.predict <- function(parms, Yobs,Tvec,plot.it="TRUE"){
	
	qp1    <- length(Yobs);
	q      <- qp1-1;
	tt     <- Tvec - Tvec[1];
	ss     <- tt[2:qp1] - tt[1:q];
	nmiss  <- ss-1;
	long.nmiss <- c(0,nmiss);
	Nmiss <- sum(nmiss)
	#cbind(long.nmiss,Tvec,Yobs)

	t.cols     <- matrix(rep(tt,each=qp1),nrow=qp1,ncol=qp1, byrow=FALSE);
	t.rows     <- t(t.cols);
	abs.diffs  <- abs(t.rows-t.cols);

	long.t <- tt[1]:max(tt)
	where.miss <- which(is.na(match(x=long.t,table=tt)),arr.ind=TRUE)
	lt.cols <- matrix(rep(long.t),nrow=(qp1+Nmiss),ncol=(qp1+Nmiss), byrow=FALSE);
	lt.rows <- t(lt.cols);
	labs.diffs <- abs(lt.rows-lt.cols);
	

	mu     <- parms[1];
	theta  <- parms[2];
	betasq <- parms[3];
	tausq  <- parms[4];
	Var.inf<- betasq/(2*theta); 
	Sigma.mat  <- Var.inf*exp(-theta*abs.diffs);
	Itausq     <- matrix(0,qp1,qp1);
	diag(Itausq)<- rep(tausq,qp1);
	V <- Sigma.mat+Itausq;

	long.V <- Var.inf*exp(-theta*labs.diffs) + diag(rep(tausq,(qp1+Nmiss)))

	Predict.t <- rep(0,qp1);
	Muvec  <- rep(mu,q);
	miss.predict <- list()
	Muvec.miss <- rep(mu,qp1);
	start.miss <- 1
	stop.miss <- 0
	for (tj in 1:qp1){
		Y.omitj <- Yobs[-tj];    #  Omit observation at time tj.
		V.omitj <- V[-tj,-tj];  #  Omit row tj and col tj from var-cov matrix.
		V12     <- V[tj,-tj];       #  Submatrix:  row tj without col tj.
		Predict.t[tj] <- mu+V12%*%ginv(V.omitj)%*%(Y.omitj-Muvec);  #  Graybill's 1976 Thm.
		
		if(long.nmiss[tj]==0){miss.predict[[tj]]<-Predict.t[tj]}else if(long.nmiss[tj]>0){

			start.miss <- stop.miss+1
			ntjmiss <- long.nmiss[tj]
			mu.miss <- rep(mu,ntjmiss);
			ind.tjmiss <- where.miss[start.miss:(start.miss+(ntjmiss-1))]
			stop.miss  <- stop.miss+ntjmiss
	
			longV12 <- long.V[ind.tjmiss,-where.miss]		
			
			miss.predict[[tj]] <- c(mu.miss + longV12%*%ginv(V)%*%(Yobs-Muvec.miss),Predict.t[tj])
		}
	}

	#cbind(long.nmiss,Tvec,Yobs,Predict.t)
	#cbind(long.t, as.vector(unlist(miss.predict)))


	Predict.t <- exp(Predict.t);
	LPredict.t <- exp(as.vector(unlist(miss.predict)))

	
	isinf <- sum(is.infinite(Predict.t))
	if(isinf>0){
	  
	  where.infs <- which(is.infinite(Predict.t)==TRUE, arr.ind=TRUE)
	  Predict.t[where.infs] <- .Machine$double.xmax
	  
	}

	isinf2 <- sum(is.infinite(LPredict.t))
	if(isinf2>0){
	  
	  where.infs <- which(is.infinite(LPredict.t)==TRUE, arr.ind=TRUE)
	  LPredict.t[where.infs] <- .Machine$double.xmax
	  
	}

	

	if(plot.it=="TRUE"){	
		#  Plot the data & model-fitted values
		X11()
		plot(Tvec,exp(Yobs),xlab="Time",ylab="Population abundance",type="b",cex=1.5, 
			main="Predicted (--) and observed (-o-) abundances");#  Population data are circles.
		par(lty="dashed"); #  Predicted abundances are dashed line.
		points(Tvec,Predict.t, type="l", lwd=1);
	}
	return(list(cbind(Tvec,Predict.t), cbind(long.t,LPredict.t) ))

}

# 10. Function to run the estimation, compute the predictions and run a parametric bootstrap
ROUSS.CALCS <- function(Yobs,Tvec,pmethod="ML",nboot, plot.pred="TRUE",plot.bootdists = "TRUE",log.out=FALSE){

	# 10.1 Compute a rough guess of the parameter estimates to initialize the search:
	guesscalc <- guess.calc(Yobs = Yobs,Tvec=Tvec)
	

	# 10.2 Compute either the ML or the REML estimates, according to what you specified in point 1. above
	if(pmethod=="ML"){
		best.guess <- ROUSS.ML(Yobs= Yobs,Tvec=Tvec, parms.guess=guesscalc);
		AIC <- best.guess[[3]];
		reml.option <- "FALSE"}else if (pmethod=="REML"){
			best.guess <- ROUSS.REML(Yobs=Yobs,Tvec=Tvec, parms.guess=guesscalc);
			reml.option <- "TRUE"}else{
				print("Error: ML and REML are the only options allowed for 'method' ")}

	# 10.3 Parameter estimates and maximized log-likelihood (we will print these at the end)
	parms.est <- best.guess[[1]];
	lnLhat    <- best.guess[[2]];

	# 10.4 Computing the predictions
	#rouss.preds <- ROUSS.predict(parms=parms.est, Yobs=Yobs,Tvec=Tvec,plot.it=plot.pred);
	#print("These are the predictions of what the true, unobserved population abundances were");
	#print(rouss.preds);

	# 10.5 Parametric bootstrap: computing both, parameters and predictions 95 % CI's
	pboot.calcs <- ROUSS.pboot(B=nboot,parms=parms.est,Yobs=Yobs, Tvec=Tvec, REML=reml.option, plot.it=plot.bootdists);
	#print("These are the predictions of what the true, unobserved population abundances were, along with 95% CI's");
	#print(pboot.calcs$preds.CIs)
	
	if(plot.pred=="TRUE"){
		points(Tvec,pboot.calcs$preds.CIs1[,2],type="l",col="blue");
		points(Tvec,pboot.calcs$preds.CIs1[,4],type="l",col="blue");
		}
	
	#print("This is the matrix of Parametric Bootstrap 95% CI's along with the estimates");
	#print(pboot.calcs$CIs.mat)
	#print("Maximized log-likelihood score")
	#print(lnLhat)
	# Next line is from ROUSS.pboot, copied here just to check object names that are needed here
	#boot.list <- list(boot.mles = boot.mles, CIs.mat = CIs.mat, preds.CIs1 = preds.CIs1, preds.CIs2 = preds.CIs2)
	
	
	if(pmethod=="ML"){
		#print("AIC score");print(AIC);
		out <- list(parms.est = parms.est, lnLhat = lnLhat, AIC = AIC, pbootmat = pboot.calcs[[1]], 
		            pboot.cis = pboot.calcs[[2]], pboot.preds1 = pboot.calcs$preds.CIs1, 
		            pboot.preds2 = pboot.calcs$preds.CIs2)}else{
		              
			out <- list(parms.est = parms.est, lnLhat = lnLhat, pbootmat = pboot.calcs[[1]], 
			            pboot.cis = pboot.calcs[[2]], pboot.preds1 = pboot.calcs$preds.CIs1, 
			            pboot.preds2 = pboot.calcs$preds.CIs2)}
	
	if(log.out==TRUE){
	  
	  out$pbootmat <- log(out$pbootmat)
	  out$pboot.cis <- log(out$pboot.cis)
	  out$pboot.preds1 <- log(out$pboot.preds1)
	  out$pboot.preds2 <- log(out$pboot.preds2)
	  
	}
	
	
	return(out);
	

}


#----------------- PART II):  The EGSS functions --------------------------------# 

# 11.  Negative Log-Lilkelihood function for the EGSS model
negloglike.EGSS.ml <- function(theta,yt,tt){
	
	mu <- theta[1];
	sigmasq <- exp(theta[2]);
	tausq   <- exp(theta[3]);
	xo      <- theta[4];
	q       <- length(yt) - 1;
	qp1     <- q+1;
	yt      <- matrix(yt,nrow=qp1,ncol=1); # makes data a matrix object
	vx      <- matrix(0,qp1,qp1);
	for(i in 1:q){ # In the next line, I preferred to create the matrix with the correct dimensions
		           # instead of relying on R's automatic recycling-to-match-dimensions property
		vx[((i+1):qp1),((i+1):qp1)] <- matrix(1,(qp1-i),(qp1-i))*tt[(i+1)];
	}
	Sigma.mat    <- sigmasq*vx;
	Itausq       <- matrix(rep(0,(qp1*qp1)), nrow=qp1, ncol=qp1);
	diag(Itausq) <- rep(tausq,qp1);
	V            <- Sigma.mat + Itausq;
	mu.vec <- matrix((xo+mu*tt), nrow=qp1,ncol=1);
	
	return((qp1/2)*log(2*pi) + 0.5*log(det(V)) + 0.5*(t(yt-mu.vec)%*%ginv(V)%*%(yt-mu.vec)))	
}


#  12. Function to propose an initial guess of the parameters to feed the parameter estimation function under the EGSS model 

init.egss <- function(Yobs,Tvec){

	q   <- length(Yobs)-1;       #  Number of time series transitions, q.
	qp1 <- q+1;                  #  q+1 gets used a lot, too.
	T.t <-Tvec-Tvec[1];          #  For calculations, time starts at zero.
	Ybar <- mean(Yobs);
	Tbar <- mean(T.t);
	mu.egoe <- sum((T.t-Tbar)*(Yobs-Ybar))/sum((T.t-Tbar)*(T.t-Tbar));
	x0.egoe <- Ybar-mu.egoe*Tbar;
	ssq.egoe <- 0;
	Yhat.egoe <- x0.egoe+mu.egoe*T.t;
	tsq.egoe <- sum((Yobs-Yhat.egoe)*(Yobs-Yhat.egoe))/(q-1);
	
	S.t <- T.t[2:qp1]-T.t[1:q];  #  Time intervals.
	Ttr <- sqrt(S.t);
	Ytr <- (Yobs[2:qp1]- Yobs[1:q])/Ttr;
	mu.egpn <- sum(Ttr*Ytr)/sum(Ttr*Ttr);
	Ytrhat <- mu.egpn*Ttr;
	ssq.egpn <- sum((Ytr-Ytrhat)*(Ytr-Ytrhat))/(q-1);
	tsq.egpn <- 0;
	x0.egpn  <- Yobs[1];

	mu0 <- (mu.egoe+mu.egpn)/2;
	ssq0 <- ssq.egpn/2;
	tsq0 <- tsq.egoe/2;
	xo.out <- x0.egoe;
	
	return(c(mu0, ssq0, tsq0, xo.out))
	
}


#  13. Function to optimize the negative log-likelihood and return parameter estimates, likelihood score and AIC, BIC value
EGSS.ML <- function(Yobs,Tvec, parms.guess){
	
	tt  <- Tvec-Tvec[1];
	q   <- length(tt) -1; 
	qp1 <- q+1;
	guess.optim <- c(parms.guess[1], log(parms.guess[2:3]), parms.guess[4]) 
	optim.out <- optim(par=guess.optim,fn=negloglike.EGSS.ml,method="Nelder-Mead",yt=Yobs,tt=tt)
	mles <- c(optim.out$par[1], exp(optim.out$par[2:3]), optim.out$par[4])
	lnL.hat <- - optim.out$value[1]
	AIC     <- -2*lnL.hat + 2*4 #where 4 = length(mles)... 
	
	out <- list(mles=mles, lnL.hat = lnL.hat, AIC=AIC)
	return(out)
}

#guess.egss   <- init.egss(Yobs=yt,Tvec=Time.t)
#trial.egssml <- EGSS.ML(Yobs=yt, Tvec=Time.t, parms.guess=guess.egss)


#  14. Simulation from the EGSS model for parametric bootstrap:  simulate data using the observed time intervals
#      and the ML estimates.  The simulation is greatly eased by the fact that the log observations are
#      multivariate normally distributed

EGSS.sim <- function(nsims,parms,Tvec){
	
	tt  <- Tvec-Tvec[1];
	mu  <- parms[1];
	sigmasq <- parms[2];
	tausq   <- parms[3];
	xo      <- parms[4];
	q       <- length(tt) - 1;
	qp1     <- q+1;
	vx      <- matrix(0,qp1,qp1);
	for(i in 1:q){ # In the next line, I preferred to create the matrix with the correct dimensions
		           # instead of relying on R's automatic recycling-to-match-dimensions property
		vx[((i+1):qp1),((i+1):qp1)] <- matrix(1,(qp1-i),(qp1-i))*tt[(i+1)];
	} 
	Sigma.mat    <- sigmasq*vx;
	Itausq       <- matrix(rep(0,(qp1*qp1)), nrow=qp1, ncol=qp1);
	diag(Itausq) <- rep(tausq,qp1);
	V            <- Sigma.mat + Itausq;
	mu.vec <- matrix((xo+mu*tt), nrow=qp1,ncol=1);
	out <- randmvn(n=nsims,mu.vec=mu.vec, cov.mat=V);
	
	return(out)
}

# PART III) Parametric bootstrap Likelihood Ratio Test
# 15. Function to do a PBLRT between the EGSS model (Null) and the OUSS model (alternative)

PBLRT.ddpcont <- function(B=10,Yobs,Tvec, alpha=0.05, plot.PBLR=TRUE){
	
	#  Estimation under the null Ho: the EGSS model
	guess.egss   <- init.egss(Yobs=Yobs,Tvec=Tvec);
	egssml       <- EGSS.ML(Yobs=Yobs, Tvec=Tvec, parms.guess=guess.egss);
	lnL.Ho       <- egssml$lnL.hat
	
	#  Estimation under the alternative H1:  the OUSS model;
	guess.rouss  <- guess.calc(Yobs=Yobs, Tvec=Tvec);
	roussml      <- ROUSS.ML(Yobs=Yobs, Tvec=Tvec, parms.guess=guess.rouss);
	lnL.H1       <- roussml$lnL.hat
	
	#  Computing the observed Lam.obs = -2*ln(L.Ho/L.H1) 
	Lam.obs <- -2*(lnL.Ho-lnL.H1)
	
	#  Simulating under Ho using the ML estimates 
	mlsforboot <- egssml$mles;
	sims.mat   <- EGSS.sim(nsims=B,parms=mlsforboot,Tvec=Tvec)
	
	#  Looping over the simulations and maximizing the likelihood for both models each time	
	#  and computing the bootstrapped Lam.boot = -2*ln(L.Ho/L.H1)
	Lam.vec <- rep(0,B);
	
	for(b in 1:B){
		
		# Estimating parameters and maximizing under the null
		Yobs.b <- sims.mat[,b];
		egssml.b       <- EGSS.ML(Yobs=Yobs.b, Tvec=Tvec, parms.guess=mlsforboot);
		lnL.Ho.b       <- egssml.b$lnL.hat

		# Estimating parameters and maximizing under the alternative		
		guess.rouss.b  <- guess.calc(Yobs=Yobs.b, Tvec=Tvec);
		roussml.b      <- ROUSS.ML(Yobs=Yobs.b, Tvec=Tvec, parms.guess=guess.rouss.b);
		lnL.H1.b       <- roussml.b$lnL.hat
		
		# Computing the boostrapped Likelihood Ratio
		Lam.boot <- -2*(lnL.Ho.b-lnL.H1.b)
		Lam.vec[b] <- Lam.boot
		
	}
    
	#  Computing the proportion of the simulations that have a more extreme Lamb.boot value
	#  than the observed Lamb.obs
	pval <- sum(Lam.vec>Lam.obs)/B;
	Decision.rule <- "Fail to Reject Density Independence"
	if(pval<alpha){	Decision.rule <- "Reject Density Independence"}
	
	#  returning out the vector of Lam.boot values and a plot if plot.PBLR=TRUE
	out <- list(egssml = egssml, roussml = roussml, Lam.obs = Lam.obs, Lam.vec = Lam.vec, pvalue = pval, Decision.rule=Decision.rule)

	if(plot.PBLR==TRUE){hist(Lam.vec, main="Pboot distribution of -2*ln(L.Ho/L.H1)"); abline(v=Lam.obs, col="red");}
	
	return(out)
	
}


# PART IV) ML estimation for the nonstationary case

negloglike.OU.ml2 <- function(fguess,yt,tt,statio=1){
	
	# if statio=1, then use the stationary likelihood. else if statio = 0, use the nonstationary case
	
	if(statio==1){
		mu <- fguess[1];
		guess  <- exp(fguess[2:4]); # Constrains parameters > 0
		theta  <- guess[1]; 
		betasq <- guess[2];
		tausq  <- guess[3];
		q      <- length(yt) - 1;
		qp1    <- q+1;
		Var.inf<- betasq/(2*theta); # Better than Dr. Seuss' thing 1 and thing 2!
		ss <- tt[2:qp1] - tt[1:q];
		t.cols     <- matrix(rep(tt,each=qp1),nrow=qp1,ncol=qp1, byrow=FALSE);
		t.rows     <- t(t.cols);
		abs.diffs  <- abs(t.rows-t.cols);
		V          <- Var.inf*exp(-theta*abs.diffs);
		diag(V)    <- diag(V) + rep(tausq,qp1);
		mu.vec     <- rep(mu,qp1);
	
		neglogl <- (qp1/2)*log(2*pi)+(1/2)*log(det(V))+ (1/2)*(yt-mu.vec)%*%ginv(V)%*%(yt-mu.vec);
		
	}else if(statio==0){
		
		mu <- fguess[1];
		guess  <- exp(fguess[2:5]); # Constrains parameters > 0
		theta  <- guess[1]; 
		betasq <- guess[2];
		tausq  <- guess[3];
		xo     <- guess[4];	
		q      <- length(yt) - 1;
		qp1    <- q+1;
		Var.inf<- betasq/(2*theta); # Better than Dr. Seuss' thing 1 and thing 2!
		ss <- tt[2:qp1] - tt[1:q];
		t.cols     <- matrix(rep(tt,each=qp1),nrow=qp1,ncol=qp1, byrow=FALSE);
		t.rows     <- t(t.cols);
		abs.diffs  <- abs(t.rows-t.cols);
		min.titj   <- pmin(t.rows,t.cols); #elementwise minima between two matrices
		V          <- Var.inf*(1-exp(-2*theta*min.titj))*exp(-theta*abs.diffs);
		diag(V)    <- diag(V) + rep(tausq,qp1);
		mu.vec     <- rep(mu,qp1) - (mu-xo)*exp(-theta*tt);
	
		neglogl    <- (qp1/2)*log(2*pi)+(1/2)*log(det(V))+ (1/2)*(yt-mu.vec)%*%ginv(V)%*%(yt-mu.vec);
		
	}
	
	return(neglogl)
	
}


ROUSS.ML2 <- function(Yobs,Tvec, parms.guess,statio=0){
	
	tt  <- Tvec-Tvec[1];
	q   <- length(tt) -1; 
	qp1 <- q+1;
	if(statio==1){
		guess.optim <- c(parms.guess[1], log(parms.guess[2:4])) 
		optim.out <- optim(par=guess.optim,fn=negloglike.OU.ml,method="Nelder-Mead",yt=Yobs,tt=tt)
		mles <- c(optim.out$par[1], exp(optim.out$par[2:4]))
		lnL.hat <- - optim.out$value[1]
		AIC     <- -2*lnL.hat + 2*4 #where 4 = length(mles)... 
	}else if(statio==0){
		
		guess.optim <- c(parms.guess[1], log(parms.guess[2:4]), Yobs[1]);
		optim.out <- optim(par=guess.optim, fn=negloglike.OU.ml2, method="Nelder-Mead", yt=Yobs, tt=tt,statio=statio);
		mles      <- c(optim.out$par[1], exp(optim.out$par[2:5]));
		lnL.hat   <- - optim.out$value[1]
		AIC       <- -2*lnL.hat + 2*5 # where 5= length(mles) 
	}
	out <- list(mles=mles, lnL.hat = lnL.hat, AIC=AIC)
	return(out)
}

joint.negloglike.OU <- function(fguess,Yobs.mat,tt,statio=0){
	
	nts <- ncol(Yobs.mat);
	nllikes <- rep(0,nts);
		for(i in 1:nts){
			Yobs <- Yobs.mat[,i]
			nllikes[i] <- negloglike.OU.ml2(fguess=fguess,yt=Yobs,tt=tt,statio=statio)
		}
	return(sum(nllikes))
	
}


JOINT.ROUSS.ML <- function(Yobs.mat, Tvec, parms.guess,statio=0){

	tt  <- Tvec-Tvec[1];
	q   <- length(tt) -1; 
	qp1 <- q+1;

	if(statio==1){
		guess.optim <- c(parms.guess[1], log(parms.guess[2:4])) 
		optim.out <- optim(par=guess.optim,fn=joint.negloglike.OU, method="Nelder-Mead",Yobs.mat=Yobs.mat,tt=tt,statio=1)
		mles <- c(optim.out$par[1],exp(optim.out$par[2:4]))
		lnL.hat <- - optim.out$value[1]
		AIC     <- -2*lnL.hat + 2*4 #where 4 = length(mles)... 
	}else if(statio==0){
		
		guess.optim <- c(parms.guess[1], log(parms.guess[2:4]), Yobs.mat[1,1]);
		optim.out <- optim(par=guess.optim, fn=joint.negloglike.OU, method="Nelder-Mead", Yobs.mat=Yobs.mat, tt=tt,statio=0);
		mles      <- c(optim.out$par[1],exp(optim.out$par[2:5]));
		lnL.hat   <- - optim.out$value[1]
		AIC       <- -2*lnL.hat + 2*5 # where 5= length(mles) 
	}
	out <- list(mles=mles, lnL.hat = lnL.hat, AIC=AIC)
	return(out)
	
}


####-------- 2 spp. simulation functions, compositional sim. ftn, Gompertz map & Cobweb, MAR(1) likelihood function  ---------#######
Twospp.mar1 <- function(A,B,Sigma,Xo,len){
	
	Xmat <- matrix(0,ncol=len,nrow=2)
	Xmat[,1] <- Xo
	
	
	for( i in 2:len){
		
		im1 <- i-1;
		Xim1  <- matrix(Xmat[,im1], nrow=2,ncol=1)
		mui <- A + B%*%Xim1
		rand.trans <- randmvn(n=1,mu.vec= mui, cov.mat=Sigma)
		Xmat[,i] <- rand.trans
	}
	return(exp(Xmat))	
}

##transformations used to constrain parameterss
atanh	<- function(r) { (1/2)*log((1+r)/(1-r)) }
tanh		<- function(x) { (exp(2*x)-1)/(exp(2*x)+1) }

mar2sppllike <- function(guess,Xmat,model="Full"){
	
	qp1 <- dim(Xmat)[2];
	q   <- qp1-1;
	A   <- matrix(exp(guess[1:2]), nrow=2,ncol=1);

	if(model=="Full"){
		Bparms <- guess[3:6]
		sigpars<- guess[7:9]
		c11  <- tanh(Bparms[1]) # constrains c11 <1
		c21  <- Bparms[2]
		c12  <- Bparms[3]
		c22  <- tanh(Bparms[4]) # constrains c22 <2
		B    <- matrix(c(c11, c21, c12, c22), nrow=2, ncol=2, byrow=TRUE)
	}else if(model=="No Competition"){

		Bparms <- guess[3:4]
		sigpars<- guess[5:7]
		c11  <- tanh(Bparms[1]) # constrains c11 <1
		c21  <- 0
		c12  <- 0
		c22  <- tanh(Bparms[2]) # constrains c22 <2
		B    <- matrix(c(c11, c21, c12, c22), nrow=2, ncol=2, byrow=TRUE)
	}
	sigsq1 <- exp(sigpars[1])
	sig12  <- sigpars[2]
	sigsq2 <- exp(sigpars[3])				
	detSig <- sigsq1*sigsq2 - sig12*sig12;
	if(((c11-1)^2<0.001) | ((c22-1)^2<0.001)){return(.Machine$double.xmax)}
	if(detSig<=0){return(.Machine$double.xmax)}
	Siginv <- (1/detSig)*matrix(c(sigsq2,-sig12,-sig12,sigsq1), nrow=2, ncol=2, byrow=TRUE)
	Part1  <- -q*log(2*pi) - (q/2)*log(detSig)
	Part2v <- rep(0,q)

	for(i in 2:qp1){

		Xt  <- Xmat[,i];
		Xtm1<- Xmat[,(i-1)];
		mut <- A + B%*%Xtm1;
		Part2v[i] <- -0.5*t(Xt-mut)%*%Siginv%*%(Xt-mut)		
	}	

	neglogl <- -Part1 -sum(Part2v)

	return(neglogl)

}

kalman.predict2spp <- function(Fracs.obs, Totvec.obs, guess.oup){ # Specification of tvec is missing!!!
	
	lents  <- length(Totvec.obs)
	tvec   <- 0:(lents-1)
	Totmat <- matrix(rep(Totvec.obs,2), nrow=2, ncol=lents, byrow=TRUE)
	Ymat.obs <- log(Totmat*Fracs.obs)

	ts1.obs <- Ymat.obs[1,]
	ts2.obs <- Ymat.obs[2,]
	
	ou.mles1 <-ROUSS.ML2(Yobs=ts1.obs, Tvec=tvec, parms.guess= guess.oup, statio=0)
	ou.hats1  <- ou.mles1$mles
	ts1.preds <- log(ROUSS.predict(parms=ou.hats1, Yobs=ts1.obs, Tvec=tvec, plot.it="FALSE")[,2])
	ts1.preds[1] <- ts1.obs[1]
	
	ou.mles2 <-ROUSS.ML2(Yobs=ts2.obs, Tvec=tvec, parms.guess= guess.oup, statio=0)
	ou.hats2  <- ou.mles2$mles
	ts2.preds <- log(ROUSS.predict(parms=ou.hats2, Yobs=ts2.obs, Tvec=tvec, plot.it="FALSE")[,2])
	ts2.preds[1] <- ts2.obs[1]
	
	out <- matrix(rbind(ts1.preds,ts2.preds), nrow=2, ncol=lents)
	return(out)
	
}

MAR1.mles <- function(guess.mar, Xmat, model="Full"){
	
	# use the function mar2sppllike
	#mar2sppllike(guess=guess.mar, Xmat= trial)
	out <- optim(par=guess.mar,fn=mar2sppllike, method="Nelder-Mead", Xmat=Xmat, model=model)
	guess <- out$par
	A   <- matrix(exp(guess[1:2]), nrow=2,ncol=1);
	
	if(model=="Full"){
		Bparms <- guess[3:6]
		sigpars<- guess[7:9]
		c11  <- tanh(Bparms[1]) # constrains c11 <1
		c21  <- Bparms[2]
		c12  <- Bparms[3]
		c22  <- tanh(Bparms[4]) # constrains c22 <2
		B    <- matrix(c(c11, c21, c12, c22), nrow=2, ncol=2, byrow=TRUE)
	}else if(model=="No Competition"){
		Bparms <- guess[3:4]
		sigpars<- guess[5:7]
		c11  <- tanh(Bparms[1]) # constrains c11 <1
		c21  <- 0
		c12  <- 0
		c22  <- tanh(Bparms[2]) # constrains c22 <2
		B    <- matrix(c(c11, c21, c12, c22), nrow=2, ncol=2, byrow=TRUE)
	}
	
	sigsq1 <- exp(sigpars[1])
	sig12  <- sigpars[2]
	sigsq2 <- exp(sigpars[3])	
	Sigma <- matrix(c(sigsq1, sig12, sig12,sigsq2), nrow=2, ncol=2, byrow=TRUE)
	
	BIC   <- 2*out$val + length(guess)*log(2*dim(Xmat)[2])
	
	outlist <- list(A, B, Sigma, negLL.hat = out$val, BIC=BIC)
	
	return(outlist)
	
}

composit.sim <- function(Truth, tausq){
	
	len <- dim(Truth)[2]
	Itausq <- tausq*diag(2)
	zeromean <- matrix(0, nrow=2,ncol=1)
	noise <- randmvn(n=len, mu.vec=zeromean, cov.mat=Itausq)
	Y.error <- exp(log(Truth) + noise)
	fracs   <- apply(Y.error, 2, FUN= function(x){x/sum(x)})
	
	return(list(Truth=Truth, Yt = Y.error, fracs=fracs))

}




########--------- Gompertz---------########
Gompertz<-function(a,b,len,no){
	
	pop.vec     <-rep(0,len);
	pop.vec[1]  <- no;
	for(i in 1:(len-1)){
		pop.vec[(i+1)] <- pop.vec[i]*exp(a+(b*log(pop.vec[i])));
	} 
	
	return(pop.vec);
}


Gompertz.map<- function(a,b,nt.vec){
	
	len<- length(nt.vec);
	ntp1.vec <- rep(0,len);
	
	for(i in 1:len){
		
		ntp1.vec[i] <- nt.vec[i]*exp(a+(b*log(nt.vec[i])));
	}
	
	return(ntp1.vec)
	
}

####- Gompertz Cobweb:

Gompertz.cobweb<-function(a,b,start.val,len){

	# len has to be an odd number	
	x.web <- rep(0,len);
	y.web <- x.web;
	x.web[1]<- start.val;
	y.web[1]<- 0;
	
	for(i in seq(2,(len-1),by=2)){
		
		next.pred <- x.web[(i-1)]*exp(a+(b*log(x.web[(i-1)])));
		x.web[i]  <- x.web[(i-1)];
		y.web[i]  <- next.pred;
		x.web[(i+1)] <- next.pred;
		y.web[(i+1)] <- next.pred;
		}		
	
	return(cbind(x.web,y.web));
}







