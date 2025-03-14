## OUSS model R functions

# 1. Negative log-likelihood score, for ML estimation or Model Selection (eq. 11)
negloglike.OU.ml=function(phi,yt,tt)  
{
   mu=phi[1];
   theta=exp(phi[2]);         #  Constrains theta > 0.
   betasq=exp(phi[3]);      #  Constrains betasq > 0. 
   tausq=exp(phi[4]);       #  Constrains tausq > 0.
   thing=betasq/(2*theta);   #  Recurring quantity: variance of statio. distrib.
   q=length(yt)-1; 			# This line was below and should have been defined here
   ss=tt[2:qp1]-tt[1:q];      #  Time intervals.
   qp1=q+1;
   vx=matrix(1,qp1,qp1);  # Preallocate matrix for autocorrelations.
   #  Following loop calculates the model autocorrelations and puts 
   #    them in vx.
   for (ti in 1:q)
   {
      #print(cumsum(ss[ti:q]))
      vx[(ti+1):qp1,ti]= exp(-theta*cumsum(ss[ti:q]));
      vx[ti,(ti+1):qp1]=vx[(ti+1):qp1,ti];
   }
   
   Sigma.mat=vx*thing;   #  Variance-covariance matrix for X(t).
   Itausq=matrix(0,qp1,qp1);
   diag(Itausq)=rep(tausq,qp1);
   V=Sigma.mat+Itausq;   #  Variance-covariance matrix for Y(t).
   mm=rep(mu,qp1);
   ofn=(qp1/2)*log(2*pi)+0.5*log(det(V))+
      0.5*(yt-mm)%*%ginv(V)%*%(yt-mm);
   return(ofn);
}

mynegloglike.OU.ml <- function(fguess,yt,tt){
	
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
negloglike.OU.ml(phi=fguess,yt=log.obs,tt=tt)
brians.nll <- proc.time()-p1tm

p2tm <- proc.time()
mynegloglike.OU.ml(fguess=fguess,yt=log.obs,tt=tt)
my.nll <- proc.time()-p2tm

print(brians.nll)
print(my.nll)


# 2. Negative log-likelihood for REML estimation or Model Selection (eq. 14)
# 3. Optimization of either Negative-loglikelihood function
# 4. For REML estimation, take the output of the optimization, read it and provide estimate of mu along with CI's
# 5. Simulation function: need to provide parameter values and length of simulations
# 6. Parametric Bootstrap Confidence Intervals for all the parameters

