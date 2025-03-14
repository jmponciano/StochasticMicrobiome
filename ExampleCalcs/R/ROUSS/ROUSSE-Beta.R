## OUSS model R functions

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
   ss     <- tt[2:qp1]-tt[1:q];         #  Time intervals.
   q   <- length(yt)-1;
   qp1 <- q+1;

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
   return(ofn);
}

# 3. rand.MVN:  Multivariate Normal random number generator
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

	out <- c(mu1,th1,bsq1,tsq1);
	return(out)
	
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
#         if 
# tt = vector of observation times (t_0, t_1, t_2, ..., t_q) from the original time series
# if REML="TRUE", then the parametric bootsrap is computed for REML estimation.  In that case, 'parms' has to
# contain the REML estimates for the original time series

ROUSS.pboot <- function(B=2,parms, Tvec, REML="FALSE", plot.it="FALSE"){
	
	tt <- Tvec-Tvec[1]
	nparms    <- length(parms);
	
	if(REML=="TRUE"){

		boot.remles <- matrix(0,nrow=B,ncol=nparms+1); 
		all.sims  <- ROUSS.sim(nsims=B,parms=parms,Tvec=Tvec);

		for(b in 1:B ){
		
			bth.timeseries <- all.sims[,b];
			remles.out <- ROUSS.REML(Yobs=bth.timeseries,Tvec=Tvec, parms.guess=parms);
			boot.remles[b,] <- c(remles.out$remls, remles.out$lnLhat);
		} 
	
		CIs.mat <- apply(boot.remles,2,FUN=function(x){quantile(x,probs=c(0.025,0.975))});
		CIs.mat <- rbind(CIs.mat[1,1:4],parms,CIs.mat[2,1:4]);
		rownames(CIs.mat) <- c("2.5%","REMLE","97.5%");
		colnames(CIs.mat) <- c("mu", "theta","betasq","tausq")
		
		boot.list <- list(boot.remles = boot.remles, CIs.mat = CIs.mat)
		
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
		for(b in 1:B ){
		
			bth.timeseries <- all.sims[,b];
			mles.out <- ROUSS.ML(Yobs= bth.timeseries,Tvec=Tvec, parms.guess=parms);
			boot.mles[b,] <- c(mles.out$mles, mles.out$lnL.hat,mles.out$AIC);
			} 
		
		CIs.mat <- apply(boot.mles,2,FUN=function(x){quantile(x,probs=c(0.025,0.975))});
		CIs.mat <- rbind(CIs.mat[1,1:4],parms,CIs.mat[2,1:4]);
		rownames(CIs.mat) <- c("2.5%","MLE","97.5%");
		colnames(CIs.mat) <- c("mu", "theta","betasq","tausq")
		
		boot.list <- list(boot.mles = boot.mles, CIs.mat = CIs.mat)

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
	mu     <- parms[1];
	theta  <- parms[2];
	betasq <- parms[3];
	tausq  <- parms[4];
	
	Var.inf<- betasq/(2*theta); 
	t.cols     <- matrix(rep(tt,each=qp1),nrow=qp1,ncol=qp1, byrow=FALSE);
	t.rows     <- t(t.cols);
	abs.diffs  <- abs(t.rows-t.cols);
	Sigma.mat  <- Var.inf*exp(-theta*abs.diffs);
	Itausq     <- matrix(0,qp1,qp1);
	diag(Itausq)<- rep(tausq,qp1);
	V <- Sigma.mat+Itausq;

	Predict.t <- rep(0,qp1);
	Muvec=rep(mu,q);
	for (tj in 1:qp1){
		Y.omitj <- Yobs[-tj];    #  Omit observation at time tj.
		V.omitj <- V[-tj,-tj];  #  Omit row tj and col tj from var-cov matrix.
		V12     <- V[tj,-tj];       #  Submatrix:  row tj without col tj.
		Predict.t[tj] <- mu+V12%*%ginv(V.omitj)%*%(Y.omitj-Muvec);  #  Usual expression for conditional MV normal mean.
	}

	Predict.t <- exp(Predict.t);

	if(plot.it=="TRUE"){	
		#  Plot the data & model-fitted values
		plot(Time.t,exp(Yobs),xlab="Time",ylab="Population abundance",type="o",pch=1,cex=1.5, main="Predicted (--) and observed (-o-) abundances");#  Population data are circles.
		par(lty="dashed"); #  Predicted abundances are dashed line.
		points(Tvec,Predict.t, type="l", lwd=1);	}
	return(cbind(Tvec,Predict.t))

}


ROUSS.CALCS <- function(Yobs,Tvec,pmethod="ML",nboot, plot.pred="TRUE",plot.bootdists = "TRUE"){

	# 5.1 Compute a rough guess of the parameter estimates to initialize the search:
	guesscalc <- guess.calc(Yobs = log.obs,Tvec=Tvec)

	# 5.2 Compute either the ML or the REML estimates, according to what you specified in point 1. above
	if(pmethod=="ML"){
		best.guess <- ROUSS.ML(Yobs= Yobs,Tvec=Tvec, parms.guess=guesscalc);
		AIC <- best.guess[[3]];
		reml.option <- "FALSE"}else if (pmethod=="REML"){
			best.guess <- ROUSS.REML(Yobs=Yobs,Tvec=Tvec, parms.guess=guesscalc);
			reml.option <- "TRUE"}else{
				print("Error: ML and REML are the only options allowed for 'method' ")}

	# 5.3 Parameter estimates and maximized log-likelihood (we will print these at the end)
	parms.est <- best.guess[[1]];
	lnLhat    <- best.guess[[2]];

	# 5.4 Computing the predictions
	rouss.preds <- ROUSS.predict(parms=parms.est, Yobs=Yobs,Tvec=Tvec,plot.it=plot.pred);
	print("These are the predictions of what the true, unobserved population abundances were");
	print(rouss.preds);

	# 5.5 Parametric bootstrap
	pboot.calcs <- ROUSS.pboot(B=nboot,parms=parms.est, Tvec=Tvec, REML=reml.option, plot.it=plot.bootdists);	print("This is the matrix of Parametric Bootstrap 95% CI's along with the estimates");
	print(pboot.calcs$CIs.mat)
	print("Maximized log-likelihood score")
	print(lnLhat)
	

	
	if(pmethod=="ML"){
		print("AIC score");print(AIC);
		out <- list(parms.est = parms.est, lnLhat = lnLhat, AIC = AIC, rouss.preds=rouss.preds, pbootmat = pboot.calcs[[1]], pboot.cis = pboot.calcs[[2]])}else{
			out <- list(parms.est = parms.est, lnLhat = lnLhat, rouss.preds=rouss.preds, pbootmat = pboot.calcs[[1]], pboot.cis = pboot.calcs[[2]])}
	
	return(out);
	

}



######################## MAIN ##############################

# Trial
library(MASS)
Observed.t=c(346,675,802,1478,1173,756,861,972,854,1161,1318,901,901,1173,
   608,811,903,584,1179,1020,1129,966);

Time.t <- c(1956,1957,1958,1959,1960,1961,1962,1963,1964,1965,1970,1971,1972,
   1973,1974,1975,1976,1977,1978,1979,1980,1981);
tt  <- Time.t-Time.t[1];
q   <- length(tt) -1; 
qp1 <- q+1;
theta <- 1.5
mu <- 1.4
betasq <- 0.09
tausq <- 0.23
fguess <- c(mu,log(c(theta,betasq,tausq)))
log.obs <- log(Observed.t)

p1tm <- proc.time()
negloglike.OU.ml(fguess=fguess,yt=log.obs,tt=tt)
my.nll.time <- proc.time()-p1tm
print(my.nll.time)

ROUSS.sim(nsims=2,parms=c(mu,theta,betasq,tausq),Tvec=Time.t)
guesscalc <- guess.calc(Yobs = log.obs,Tvec=Time.t)
ROUSS.ML(Yobs= log.obs,Tvec=Time.t, parms.guess=guesscalc)
$mles
[1] 6.7947258819 1.6925283890 0.3339970731 0.0001871666

$lnL.hat
[1] -5.394046

$AIC
[1] 18.78809

ROUSS.REML(Yobs=log.obs,Tvec=Time.t, parms.guess=guesscalc)
$remls
[1] 6.792985294 1.259852634 0.271751547 0.000748084

$lnLhat
[1] -6.94946

linx.remles <- c(6.792985294, 1.259852634, 0.271751547, 0.000748084);
linx.mles <- c(6.7947258819, 1.6925283890, 0.3339970731, 0.0001871666);

ROUSS.predict(parms=linx.remles, Yobs=log.obs,Tvec=Time.t,plot.it="TRUE")

trial.pboot0 <- ROUSS.pboot(B=200,parms=linx.remles, Tvec=Time.t, REML="TRUE", plot.it="TRUE")
trial.pboot0

trial.pboot <- ROUSS.pboot(B=1000,parms=linx.mles, Tvec=Time.t,REML="FALSE")
par(mfrow=c(2,2))
hist(trial.pboot[,1],main="mu",xlab="");abline(v=linx.mles[1],lwd=2,col="red");
hist(trial.pboot[,2],main="theta",xlab="");abline(v=linx.mles[2],lwd=2,col="red");
hist(trial.pboot[,3],main="betasq",xlab="");abline(v=linx.mles[3],lwd=2,col="red");
hist(trial.pboot[,4],main="tausq",xlab="");abline(v=linx.mles[4],lwd=2,col="red");

trialrml.pboot <- ROUSS.pboot(B=1000,parms=linx.remles, Tvec=Time.t,REML="FALSE")
par(mfrow=c(2,2))
hist(trialrml.pboot[,1],main="mu",xlab="");abline(v=linx.remles[1],lwd=2,col="red");
hist(trialrml.pboot[,2],main="theta",xlab="");abline(v=linx.remles[2],lwd=2,col="red");
hist(trialrml.pboot[,3],main="betasq",xlab="");abline(v=linx.remles[3],lwd=2,col="red");
hist(trialrml.pboot[,4],main="tausq",xlab="");abline(v=linx.remles[4],lwd=2,col="red");


# 5.1 Compute a rough guess of the parameter estimates to initialize the search:
guesscalc <- guess.calc(Yobs = log.obs,Tvec=Time.t)

# 5.2 Compute either the ML or the REML estimates, according to what you specified in point 1. above
if(method=="ML"){
	best.guess <- ROUSS.ML(Yobs= log.obs,Tvec=Time.t, parms.guess=guesscalc);
	AIC <- best.guess[[3]];
	reml.option <- "FALSE"}else if (method=="REML"){
		best.guess <- ROUSS.REML(Yobs=log.obs,Tvec=Time.t, parms.guess=guesscalc);
		reml.option <- "TRUE"}else{
			print("Error: ML and REML are the only options allowed for 'method' ")}

# 5.3 Parameter estimates and maximized log-likelihood (we will print these at the end)
parms.est <- best.guess[[1]]
lnLhat    <- best.guess[[2]]

# 5.4 Computing the predictions
rouss.preds <- ROUSS.predict(parms=parms.est, Yobs=log.obs,Tvec=Time.t,plot.it=pred.plot)
print("These are the predictions of what the true, unobserved population abundances were")
print(rouss.preds)

# 5.5 Parametric bootstrap
pboot.calcs <- ROUSS.pboot(B=NBoot,parms=parms.est, Tvec=Time.t, REML=reml.option, plot.it=pboot.plot)
print("This is the matrix of Parametric Bootstrap 95% CI's along with the estimates")
print(pboot.calcs$CIs.mat)
print("Maximized log-likelihood score")
print(lnLhat)
if(method=="ML"){print("AIC score");print(AIC)}


