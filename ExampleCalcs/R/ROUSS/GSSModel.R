#  This is an R file that contains functions to estimate the GSS model parameters and simulate
#  data.  Here are the contents:

#1. Likelihood functions for
# a. ML, univariate time series
# b. REML, univariate time series
# c. Replicated samples MLEs 
#
#2. Functions to simulate data for
# a. replicated samples time series, starting from the stationary distribution
# b. univariate time series, starting from the stationary distribution

# 3. Examples for each of the methods
# All programs written by Jose Miguel Ponciano

library(MASS);


# 1.
################ Likelihood functions for ML, REML, RSML

# a.  First the Multivariate matrix notation to get initial parameter estimates;
#     guess = = [log(a) 1-log(c) log(sigsq) log(tausq)]'
#     yt = a vector of data
#     Assumes first observation comes from the stationary distribution
kalman.iter2<- function(guess,yt){
    a            <- exp(guess[1]);   #constrains a >0;
    cc           <- 1-exp(guess[2]); #constrains c < 1 (from negloglike.m); 
    #if(cc==1){cc <- 0.9999};   
    sigmasq      <- exp(guess[3]);   #constrains sigmasq > 0 (from negloglike.m); 
    tausq        <- exp(guess[4]);   #constrains tausq > 0 (from negloglike.m)
    q            <- length(yt)-1;
    qp1          <- q + 1;
    yt           <- matrix(log(yt),nrow= qp1,ncol=1);
    sigmat       <- toeplitz(c(1,cumprod(rep(cc,q))))*(sigmasq/(1-(cc^2)));
    Itausq       <- matrix(rep(0,(qp1*qp1)),nrow = q+1,ncol = q+1);
    diag(Itausq) <- rep(tausq,q+1);
    V            <- sigmat + Itausq; 
    mu           <- matrix((a/(1-cc))*rep(1,qp1),nrow=qp1,ncol=1);
    ofn <- ((qp1)/2)*log(2*pi) + (0.5*log(det(V))) + (0.5*(t(yt-mu)%*%ginv(V)%*%(yt-mu)));
    return(ofn);

}

# b. The REML estimates:
# Theta is an initial geuss of parameter estimates:
# theta = [1-log(c) log(sigsq) log(tausq)]'
negloglike <- function(theta,w){
  cc     <-  1-exp(theta[1]);
  sigsq  <- exp(theta[2]);
  tausq  <- exp(theta[3]);
  q      <- length(w);
  sigmat <- toeplitz(c(1,cumprod(rep(cc,q))))*(sigsq/(1-(cc^2)));
  Dmat   <- toeplitz(c(-1,1,rep(0,q)))
  for(i in 1:q){Dmat[i+1,i]<-0};
  qp1    <- q+1;
  Dmat   <- Dmat[1:q,1:qp1];
  Identi <- matrix(rep(0,(qp1*qp1)),nrow = q+1,ncol = q+1);
  diag(Identi) <- rep(tausq,q+1);
  phimat <- Dmat%*%(sigmat + Identi)%*%t(Dmat);

  out <- (q/2)*log(2*pi) + (0.5*log(det(phimat))) + (0.5*(w%*%ginv(phimat)%*%w));
  return(out);
}


# c. This function calculates the Kalman Iterated Form  for REPLICATED SAMPLING
# of the negative-loglikelihood ,eqs. 14-23:
# data is a matrix, with replicates in rows and time in columns, 
# guess is a vector of parameters guess, statio refers to whether
# or not the stationary assumption is applicable (1=yes, 0=no)
kalman.nrep <- function(guess,yt,statio){
# if statio = 1, then guess =[a, c, sigmasq, tausq], IF statio=0, guess=[x0,a,c, sigmasq, tausq]
    yt      <-  log(yt);
    len     <- ncol(yt);
    pt      <- nrow(yt);
if (statio==1){
    a          <- exp(guess[1]);   #constrains a >0;
    cc         <- 1-exp(guess[2]); #constrains c < 1 (from negloglike.m); 
    sigmasq    <- exp(guess[3]);   #constrains sigmasq > 0 (from negloglike.m); 
    tausq      <- exp(guess[4]);   #constrains tausq > 0 (from negloglike.m)
    ut         <- matrix(rep(0,len),nrow = len,ncol=1);
    psitsq     <- ut;
    mt         <- matrix(rep(0,pt*len),nrow=pt,ncol=len);
    Vt         <- matrix(rep(0,pt*pt*len),nrow=pt,ncol=pt*len);
    ut[1,1]    <- a/(1-cc);
    psitsq[1,1]<- (sigmasq/(1-cc^2));
    Itausq     <- diag(rep(tausq,pt));
    mt[,1]     <- rep(ut[1,1],pt);
    Vt0        <- matrix(rep(psitsq[1,1],pt*pt),nrow=pt,ncol=pt) + Itausq;
    Vt[1:pt,1:pt]<- Vt0; 
        for (i in 2:len){
            jt.1        <- matrix(rep(1,pt),nrow=1,ncol=pt);#row vect. (1xpt);
            ut[i,1]     <- a + cc*(ut[i-1] +  psitsq[i-1]*jt.1%*%ginv(Vt[,(pt*(i-2)+1):((i-1)*pt)])%*%(yt[,i-1]-mt[,i-1]));
            psitsq[i,1] <- (cc^2)*psitsq[i-1]*(1- psitsq[i-1]*jt.1%*%ginv(Vt[,(pt*(i-2)+1):((i-1)*pt)])%*%t(jt.1)) + sigmasq;
            mt[,i]      <- t(jt.1)*ut[i,1];
            Vt[,(pt*(i-1) +1):(i*pt)]  <- matrix(rep(psitsq[i,1],pt*pt),nrow=pt,ncol=pt) + Itausq;}
  }
          
    else {
    Xo         <- guess[1];
    a          <- exp(guess[2]);   #constrains a >0;
    cc         <- 1-exp(guess[3]); #constrains c < 1 (from negloglike.m); 
    sigmasq    <- exp(guess[4]);   #constrains sigmasq > 0 (from negloglike.m); 
    tausq      <- exp(guess[5]);   #constrains tausq > 0 (from negloglike.m)
    ut         <- matrix(rep(0,len),nrow = len,ncol=1);
    psitsq     <- ut;
    mt         <- matrix(rep(0,pt*len),nrow=pt,ncol=len);
    Vt         <- matrix(rep(0,pt*pt*len),nrow=pt,ncol=pt*len);
    ut[1,1]    <- Xo;
    #psitsq[1,1]<- (sigmasq/(1-c^2));
    Itausq     <- diag(rep(tausq,pt));
    mt[,1]     <- rep(ut[1,1],pt);
    Vt0        <- Itausq;
    Vt[1:pt,1:pt]<- Vt0; 
        for (i in 2:len){
            jt.1        <- matrix(rep(1,pt),nrow=1,ncol=pt);#row vect. (1xpt);
            ut[i,1]     <- a + cc*(ut[i-1,1] +  psitsq[i-1,1]*jt.1%*%ginv(Vt[,(pt*(i-2)+1):((i-1)*pt)])%*%(yt[,i-1]-mt[,i-1]));
            psitsq[i,1] <- (cc^2)*psitsq[i-1,1]*(1- psitsq[i-1,1]*jt.1%*%ginv(Vt[,(pt*(i-2)+1):((i-1)*pt)])%*%t(jt.1)) + sigmasq;
            mt[,i]      <- t(jt.1)*ut[i,1];
            Vt[,(pt*(i-1) +1):(i*pt)]  <- matrix(rep(psitsq[i,1],pt*pt),nrow=pt,ncol=pt) + Itausq;}
  }
   
    detvect  <- rep(0,len);
    quadform <- detvect;
    
    for (i in 1:len){
        detvect[i] <- det(Vt[,(pt*(i-1)+1):(i*pt)]);
        quadform[i]<- t(yt[,i]-mt[,i])%*%ginv(Vt[,(pt*(i-1)+1):(i*pt)])%*%(yt[,i]-mt[,i]);}

    p     <- pt*len;
    logL <- (p/2)*log(2*pi) + 0.5*sum(log(detvect)) + 0.5*sum(quadform);
    return(logL)
  }


# 2.
########## Functions used to simulate data:############
# a. simulating a replicated samples time series with an even number (pt) of replicates
# and of length 'len'.  "guess" is the initial guess for the parameter estimates [a,c,sigmasq,tausq]'

spsimmv3 <- function(guess,pt,len){
a         <- guess[1]; 
cc        <- guess[2]; 
sigmasq   <- guess[3]; 
tausq     <- guess[4];    
npred     <- matrix(rep(0,len*pt),nrow=pt,ncol=len);
Et        <- rnorm(len,0,sqrt(sigmasq));
Ftprime   <- matrix(rnorm(len*pt,0,1),nrow=pt,ncol=len);
Itausq    <- diag(rep(tausq,pt));
Te        <- chol(Itausq);
F1        <- t(Te)%*%Ftprime[,1];
statiom   <- a/(1-cc);
statiovar <- (sigmasq/(1-cc^2)); #From the X process
statiosd  <- sqrt(statiovar);
X0        <- statiom + rnorm(1,mean=0, sd=statiosd);
npred[,1] <- rep(X0,pt) + F1;
jt        <- rep(1,pt);
for (i in 2:len){
  npred[,i] <- a*jt + cc*npred[,(i-1)] + Et[i]*jt + t(Te)%*%Ftprime[,i] - cc*(t(Te)%*%Ftprime[,(i-1)]);
}

Yt        <- exp(npred);
return(Yt);

}


#b. Now to simulate just a univariate time series, no replicas:
spsim3 <- function(guess,len){
a         <- guess[1]; 
cc        <- guess[2]; 
sigmasq   <- guess[3]; 
tausq     <- guess[4];    
npred     <- rep(0,len);
Et        <- rnorm(len,0,sqrt(sigmasq));
Ft        <- rnorm(len,0,sqrt(tausq));
statiom   <- a/(1-cc);
statiovar <- (sigmasq/(1-cc^2)) + tausq;
statiosd  <- sqrt(statiovar);
npred[1]  <- rnorm(1,mean = statiom,sd = statiosd);

for (i in 2:len){
  npred[i] <- (a + cc*npred[(i-1)] + Et[i] + Ft[i] - cc*Ft[(i-1)]); 
}

Yt        <- exp(npred);
return(Yt);

}

# spsim3(c(0.39,0.79,0.10,0.23), 30)

# Function to compute the Kalman predictions
# The arguments are: 1) parms.hat = ml or reml estimates of the GSS model parameters
#                    2) yt: original log-time series of observations
# Returns: 1) The one-step predicted OBSERVED values: mt, the vector of means of the Yts
#          2) The estimated mean vector of the true (Xs) population log abundances
#          3) The estimated variance vector of the true (Xs) population log abundances
#             Estimates in 2 and 3 correspond to Kalman-updated estimates
pred.kalman <- function(parms.hat,yt){
	
	a     <- parms.hat[1];
	cc    <- parms.hat[2];
	sigsq <- parms.hat[3];
	tausq <- parms.hat[4];
	
	qp1    <- length(yt); # This is q + 1, an inheritance of the Dennis and Taper 1994 notation.
	mvec   <- rep(0,qp1);
	muvec  <- rep(0,qp1);
	vsqvec <- rep(0,qp1);
	psisqv <- rep(0,qp1);
	alphav <- rep(0,qp1);
	betasqv<- rep(0,qp1);
	
	mvec[1]   <- a/(1-cc);
	vsqvec[1] <- sigsq/(1-cc*cc) + tausq;
	
	muvec[1]  <- mvec[1];
	psisqv[1] <- sigsq/(1-cc*cc);  
	
	alphav[1] <- mvec[1] + psisqv[1]*(yt[1]-mvec[1])/vsqvec[1];
	betasqv[1]<- (psisqv[1]/vsqvec[1])*tausq;
	
	for(i in 2:qp1){
		
		im1       <- i-1;
		mvec[i]   <- a+cc*(mvec[im1] + ((vsqvec[im1]-tausq)/vsqvec[im1])*(yt[im1]-mvec[im1]));
		vsqvec[i] <- cc*cc*((vsqvec[im1]-tausq)/vsqvec[im1])*tausq + sigsq + tausq;
		muvec[i]  <- mvec[i];
		psisqv[i] <- vsqvec[i] - tausq;
		alphav[i] <- muvec[i] + psisqv[i]*((yt[i] - mvec[i])/vsqvec[i]);
		betasqv[i]<- psisqv[1]*(1-(psisqv[i]/vsqvec[i]));
	}
	
	one.step.obs <- mvec;
	true.xs.mean <- alphav;
	true.xs.var   <- betasqv;
	
	out <- list(one.step.obs = mvec, true.xs.mean = alphav, true.xs.var = betasqv);
	
	return(out)
	
}


# 3. Examples for each of the methods:

# The ML estimates:
Setophaga <- c(18,10,9,14,17,14,5,10,9,5,11,11,4,5,4,8,2,3,9,2,4,7,4,1,2,4,11,11,9,6);
SERUmles <- optim(par=c(log(0.39),log(1-0.79),log(0.09),log(0.2)),kalman.iter2,NULL,method="Nelder-Mead",yt = Setophaga);
SERUmles <- c(exp(SERUmles$par[1]),1-exp(SERUmles$par[2]),exp(SERUmles$par[3]),exp(SERUmles$par[4]));
SERUmles

kalmanpreds <- pred.kalman(parms.hat=SERUmles, yt= log(Setophaga))
pred.means  <- exp(kalmanpreds$true.xs.mean)

plot(1:length(Setophaga), Setophaga, type="b", pch=1);
points(1:length(Setophaga), pred.means, type="b", pch=16);

# The REML estimates:

#Guetting the first differences:
startmean <- log(Setophaga);
q  <- length(startmean) -1;
y1 <- startmean[1:q];
y2 <- startmean[2:(q+1)];
w  <- y2-y1;
qp1 <-q+1;
#negloglike(c((1-log(0.79)),(log(0.09)),(log(0.2315))), w)
# Optimizing the negative log-likelihood
sr <- optim(par=c(log(1-0.793),log(0.097),log(0.231)),negloglike,NULL,method="Nelder-Mead",control = list(maxit=1000, trace=F),w=w);
rml<- c((1-exp(sr$par[1])),exp(sr$par[2]),exp(sr$par[3]));

#Calculating the estimated var-cov matrix and then calculating the estimate of 'a';
cc           <- rml[1];sigsq <- rml[2]; tausq <-rml[3]; 
sigmat       <- toeplitz(c(1,cumprod(rep(cc,q))))*(sigsq/(1-(cc^2)));
Identi       <- matrix(rep(0,(qp1*qp1)),nrow = q+1,ncol = q+1);
diag(Identi) <- rep(tausq,qp1);
psimat       <- sigmat + Identi;
psimatinv    <- ginv(psimat);
jvect        <- rep(1,qp1);
a            <- (jvect%*%psimatinv%*%startmean)/(jvect%*%psimatinv%*%jvect);
a            <- a*(1-cc); 

SERUremls <- c(a,cc,sigsq,tausq);
SERUremls
# The RSML estimates:

kalmanpreds <- pred.kalman(parms.hat=SERUremls, yt= log(Setophaga))
pred.means2  <- exp(kalmanpreds$true.xs.mean)

####  This corresponds to Figure 1 in Dennis et al 2006.
plot(1:length(Setophaga), Setophaga, type="b", pch=1);
points(1:length(Setophaga), pred.means2, type="b", pch=16);


# Simulating data, 2 replicates, 30 lenght :
trial <- spsimmv3(c(0.3929,0.7934,0.09725,0.2315), 2,30)

kalman.nrep(c(log(0.3929),log(1-0.7934),log(0.09725),log(0.2315)),trial,statio=1)
#Guetting the RSMLEs for the simulated data 'trial':
SERUrsmles <- optim(par=c(log(0.39),log(1-0.79),log(0.09),log(0.2)),kalman.nrep,NULL,method="Nelder-Mead",yt = trial,statio=1);
SERUrsmles <- c(exp(SERUrsmles$par[1]),1-exp(SERUrsmles$par[2]),exp(SERUrsmles$par[3]),exp(SERUrsmles$par[4]));
SERUrsmles
