#### This is the 'Master File'.  It contains the main functions regarding simulation and estimation.  
#### This file needs to be sourced in order to run all the other R files in this folder
#### All programs written by Jose Miguel Ponciano at University of Florida, Spring 2014

library(MASS)


# Function to draw random vectors from a Multivariate normal distribution with mean vector 'mu.vec' and variance 'cov.mat'
my.rmvn <- function(n,mu.vec, cov.mat){
	
	p <- length(mu.vec);
	Tau <- chol(cov.mat);
	Zmat <- matrix(rnorm(n=p*n,mean=0,sd=1),nrow=p,ncol=n);
	out <- matrix(0,nrow=p,ncol=n);
	for(i in 1:n){
		
		Z <- Zmat[,i];
		out[,i] <- t(Tau)%*%Z + mu.vec
		
		}
	
	return(out)
	
}


# Calculator function to compute the mean and the variance-covariance matrix of the process over time

mu.vcov.calc <- function(a1,c1,sigmasq1,a2,c2,sigmasq2,tau,len){
	
	taupqp1 <- len;
	q <- taupqp1-tau-1;
	qp1 <- q + 1;

	#set up matrix calcs.
	c1.taucum   <- cumprod(rep(c1,tau));
	c2.qcum     <- cumprod(rep(c2,q));
	
	tau.ones      <- rep(1,tau);
	c1.taudenom   <- rep((1-c1),tau);
	g1.part1      <- (tau.ones-c1.taucum)/c1.taudenom;
	tauth.term    <- g1.part1[tau]; 
	g1.part2      <- c2.qcum*tauth.term; 
	g1            <- matrix(c(0,g1.part1,g1.part2), nrow=taupqp1,ncol=1,byrow=TRUE);
	
	q.ones        <- rep(1,q);
	c2.qdenom     <- rep((1-c2),q);
	g2.part2      <- (q.ones-c2.qcum)/c2.qdenom;
    g2.part1      <- rep(0, (tau+1));
    g2            <- matrix(c(g2.part1,g2.part2), nrow=taupqp1, ncol=1, byrow=TRUE);
    
    Gamma11       <- toeplitz(c(1,c1.taucum));
    zind          <- which(upper.tri(Gamma11)==TRUE, arr.ind=TRUE);
    Gamma11[zind] <- Gamma11[zind]*0; # dim=(tau+1) x (tau+1)
    
    Gamma12       <- matrix(0,nrow=(tau+1), ncol=q); # dim = (tau+1) x q
    
    Gamma22       <- toeplitz(c(1,c2.qcum[1:(q-1)])); 
	zind2         <- which(upper.tri(Gamma22)==TRUE, arr.ind=TRUE);
    Gamma22[zind2]<- Gamma22[zind2]*0; # dim = q x q
    
    rev.c1s       <- c(rev(c1.taucum),1); # length =  tau+1
    c1.taumat     <- matrix( rep(rev.c1s,each=q), nrow=q, ncol=(tau+1), byrow=FALSE);
    c2.qmat       <- matrix(rep(c2.qcum, (tau+1)), nrow=q, ncol=(tau+1), byrow=FALSE);
    
    Gamma21       <-  c2.qmat*c1.taumat; # dim q x (tau + 1) 

	Gamma         <- rbind(cbind(Gamma11, Gamma12),cbind(Gamma21,Gamma22));

	# MVN mean vector
	mvec          <- matrix( c(a1/(1-c1), rep(0,tau+q)), nrow=taupqp1, ncol=1, byrow=TRUE);
	mu.tau        <- a1*g1 + a2*g2 + Gamma%*%mvec

	# MVN Covariance matrix Sigma.tau = Gamma x Phi x Gamma'
	Phi           <- matrix(0, nrow=taupqp1, ncol=taupqp1)
	diag(Phi)     <- c(sigmasq1/(1-c1^2), rep(sigmasq1,tau), rep(sigmasq2,q))
	Sigma.tau     <- Gamma%*%Phi%*%t(Gamma)
	
	return(list(mu.vec=mu.tau, Sigma = Sigma.tau))
	
}

## Multivariate Normal likelihood calculation with an added change point in the dynamics:
## Model flag refers to one of the 8 possible sub-models contained in the general model that 
## says that the dynamics is different before and after the breakpoint: 
## only one parameter may change after the breakpoint 'M', or any two, or any three, etc...
## sampling error, however, is assumed to remain the same.

model.profmle <- function(par.vec, data.vec, M, model.flag){
	
	tau     <- M;
	taupqp1 <- length(data.vec); # tau + q + 1
	q       <- taupqp1-tau-1;	
	qp1     <- q + 1;	

	
	if(model.flag=="m1"){
		
		a1      <- exp(par.vec[1]); # constrains a>0
		c1      <- tanh(par.vec[2]);# constrains c<1
		sigmasq1  <- exp(par.vec[3]); # constrains sigmasq>0
		vsq.val <- exp(par.vec[4]); # constrains nusq >0
		vsq.vec <- rep(vsq.val, taupqp1);
		a2 <- a1; c2 <- c1; sigmasq2 <- sigmasq1;
		
	}else if(model.flag=="m2"){
		
		a1       <- exp(par.vec[1]);
		c1       <- tanh(par.vec[2]);
		sigmasq1   <- exp(par.vec[3]);
		sigmasq2   <- exp(par.vec[4]);
		a2 <- a1; c2 <- c1;
		
	}else if(model.flag=="m3"){
		
		a1       <- exp(par.vec[1]);
		c1       <- tanh(par.vec[2]);
		c2       <- tanh(par.vec[3]);
		sigmasq1 <- exp(par.vec[4]); 
		a2 <- a1; sigmasq2 <- sigmasq1;
		
	}else if(model.flag=="m4"){
		
		a1       <- exp(par.vec[1]);
		a2       <- exp(par.vec[2]); 
		c1       <- tanh(par.vec[3]); 
		sigmasq1 <- exp(par.vec[4]); 
		c2 <- c1; sigmasq2 <- sigmasq1;
		
	}else if(model.flag=="m5"){
		
		a1       <- exp(par.vec[1]); 
		a2       <- exp(par.vec[2]); 
		c1       <- tanh(par.vec[3]);
		c2       <- tanh(par.vec[4]);		
		sigmasq1 <-  exp(par.vec[5]); 
		sigmasq2 <- sigmasq1; 
		
	}else if(model.flag=="m6"){
		
		a1       <- exp(par.vec[1]) ;
		c1       <- tanh(par.vec[2]);
		c2       <- tanh(par.vec[3]);
		sigmasq1 <- exp(par.vec[4]);
		sigmasq2 <- exp(par.vec[5]);
		a2       <- a1;

	}else if(model.flag=="m7"){

		a1       <- exp(par.vec[1]); 
		a2       <- exp(par.vec[2]); 
		c1       <- tanh(par.vec[3]); 
		sigmasq1 <- exp(par.vec[4]);
		sigmasq2 <- exp(par.vec[5]);
		c2       <- c1;		
		
	}else if(model.flag=="m8"){
		
		a1       <- exp(par.vec[1]); 
		a2       <- exp(par.vec[2]); 
		c1       <- tanh(par.vec[3]); 
		c2       <- tanh(par.vec[4]);
		sigmasq1 <- exp(par.vec[5]);
		sigmasq2 <- exp(par.vec[6]);
		
	}else{print("Model number not specified")}

	#if(sigmasq1 < sigma.thresh) {sigma.thresh.flag <<-TRUE; sigmasq1 <- sigma.thresh}
	#if(sigmasq2 < sigma.thresh) {sigma.thresh.flag <<-TRUE; sigmasq2 <- sigma.thresh}
	if(any(!is.finite(c(a1, a2, c1, c2, sigmasq1, sigmasq2)))) { return(500000) }else{
	
	#set up matrix calcs.
	c1.taucum   <- cumprod(rep(c1,tau));
	c2.qcum     <- cumprod(rep(c2,q));
	
	tau.ones      <- rep(1,tau);
	c1.taudenom   <- rep((1-c1),tau);
	g1.part1      <- (tau.ones-c1.taucum)/c1.taudenom;
	tauth.term    <- g1.part1[tau]; 
	g1.part2      <- c2.qcum*tauth.term; 
	g1            <- matrix(c(0,g1.part1,g1.part2), nrow=taupqp1,ncol=1,byrow=TRUE);
	
	q.ones        <- rep(1,q);
	c2.qdenom     <- rep((1-c2),q);
	g2.part2      <- (q.ones-c2.qcum)/c2.qdenom;
    g2.part1      <- rep(0, (tau+1));
    g2            <- matrix(c(g2.part1,g2.part2), nrow=taupqp1, ncol=1, byrow=TRUE);
    
    Gamma11       <- toeplitz(c(1,c1.taucum));
    zind          <- which(upper.tri(Gamma11)==TRUE, arr.ind=TRUE);
    Gamma11[zind] <- Gamma11[zind]*0; # dim=(tau+1) x (tau+1)
    
    Gamma12       <- matrix(0,nrow=(tau+1), ncol=q); # dim = (tau+1) x q
    
    Gamma22       <- toeplitz(c(1,c2.qcum[1:(q-1)])); 
	zind2         <- which(upper.tri(Gamma22)==TRUE, arr.ind=TRUE);
    Gamma22[zind2]<- Gamma22[zind2]*0; # dim = q x q
    
    rev.c1s       <- c(rev(c1.taucum),1); # length =  tau+1
    c1.taumat     <- matrix( rep(rev.c1s,each=q), nrow=q, ncol=(tau+1), byrow=FALSE);
    c2.qmat       <- matrix(rep(c2.qcum, (tau+1)), nrow=q, ncol=(tau+1), byrow=FALSE);
    
    Gamma21       <-  c2.qmat*c1.taumat; # dim q x (tau + 1) 

	Gamma         <- rbind(cbind(Gamma11, Gamma12),cbind(Gamma21,Gamma22));

	# MVN mean vector
	mvec          <- matrix( c(a1/(1-c1), rep(0,tau+q)), nrow=taupqp1, ncol=1, byrow=TRUE);
	mu.tau        <- a1*g1 + a2*g2 + Gamma%*%mvec

	# MVN Covariance matrix Sigma.tau = Gamma x Phi x Gamma'
	Phi           <- matrix(0, nrow=taupqp1, ncol=taupqp1)
	diag(Phi)     <- c(sigmasq1/(1-c1^2), rep(sigmasq1,tau), rep(sigmasq2,q))
	Sigma.tau     <- Gamma%*%Phi%*%t(Gamma)
	
	missing.vec   <- which(is.na(data.vec), arr.ind=TRUE)
	nmiss         <- length(missing.vec)
	yt            <- matrix(data.vec, nrow=taupqp1, ncol=1)
	new.len       <- taupqp1 - nmiss;
	 
	if(nmiss>0){ 
		
		mu.tau    <- mu.tau[-missing.vec];
		Sigma.tau <- Sigma.tau[-missing.vec,];
		Sigma.tau <- Sigma.tau[,-missing.vec];
		yt        <- yt[-missing.vec]; 
		vsq.vec   <- vsq.vec[-missing.vec];
	}

	# Observation error var-covar
	vsq.vec  <- rep(exp(par.vec[7]), nrow(Sigma.tau));
	V.tau         <- Sigma.tau + diag(vsq.vec);
    if(any(diag(V.tau) < 0 )) {return(NA)}    
    Vinv			<- try(solve(V.tau,tol=1e-25),silent=T);


	Negloglike    <- ((new.len)/2)*log(2*pi) + (0.5*log(det(V.tau))) + (0.5*(t(yt-mu.tau)%*%Vinv%*%(yt-mu.tau)))	
	
	return(Negloglike)
		
	}
}


###  Function to simulate under the demographic + two forms of environmental noise!!
negbinde.sim <-function(no,len,scalep,shapep,ddp=list("Ricker",b=-0.0004)){
	#  ATT.!: b< 0
	alpha      <- scalep;
	k          <- shapep;
	pop.vec    <- rep(0,len);
	pop.vec[1] <- no;
	lam.vec <- rgamma(n=len,shape=k,rate=alpha);
	
	model.form <- ddp[[1]]
	
	for(i in 2:len){
		
		lam <- lam.vec[i]
		nt <- pop.vec[(i-1)];
		if(model.form=="Ricker"){
			b  <- ddp[[2]]
			pt <- exp(b*nt);}
		else if(model.form=="Below"){
			K  <- ddp[[2]]
			beta <- ddp[[3]]
			pt <- 1/(1+ (lam-1)*(nt/K)^(beta))}
		else{print("Crap! You chose a non-specified model!");break;}
			
		s  <- rpois(n=1,lambda=nt*pt*lam);
		pop.vec[i] <- s # s, for survivors...
		
		}
	return(pop.vec)
		
}


#############---------- Auxiliary functions ---------------######


##transformations used to constrain parameterss
atanh	<- function(r) { (1/2)*log((1+r)/(1-r)) }
tanh		<- function(x) { (exp(2*x)-1)/(exp(2*x)+1) }


# Function to create curly braces
# x, y position where to put the braces
# range is the widht
# position: 1 vertical, 2 horizontal
# direction: 1 left/down, 2 right/up
CurlyBraces <- function(x, y, mrange, direction = 1 ,my.lwd=2) {

    a=c(1,2,3,48,50)    # set flexion point for spline
    b=c(0,.2,.28,.7,.8) # set depth for spline flexion point

    curve = spline(a, b, n = 50, method = "natural")$y / 2 

    curve = c(curve,rev(curve))

    a_sequence = rep(x,100)
    b_sequence = seq(y-mrange/2,y+mrange/2,length=100)  

    # direction
    if(direction==1)
    a_sequence = a_sequence+curve
    if(direction==2)
    a_sequence = a_sequence-curve

    # pos
    #if(pos==1)
    #points(a_sequence,b_sequence,type="l",lwd=my.lwd) # vertical
    #if(pos==2)
    #points(b_sequence,a_sequence,type="l",lwd=my.lwd) # horizontal
	
	return(list(a_sequence=a_sequence,b_sequence=b_sequence))
	
    }

#plot(0,0,ylim=c(-10,10),xlim=c(-30,30))
#CurlyBraces(2, 0, 10, pos = 1, direction = 1 )
#CurlyBraces(2, 0, 10, pos = 2, direction = 2 )
#CurlyBraces(2, 0, 5,  pos = 1, direction = 2 )
#CurlyBraces(1, 0, 10, pos = 2, direction = 1 )
#CurlyBraces(1, 0, 5,  pos = 2, direction = 2 )


