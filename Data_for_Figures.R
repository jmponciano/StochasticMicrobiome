### Program to plot the different scenarios in strength of inter and intraspecifc interactions in time series of community abundances
### In the system there is two competing species and one mutualistic species. 
# Species 1 competes with species 2: Negative coefficients
# Species 3 is mutualistic with species 1: Positive coefficients
# Speices 2 and species 3 compete: Negative coefficients
# in all examples, environmental variance is kept constant
# A is also kept constant
# In all cases intraspecific competition is negative reflecting density dependence


source("Poisson-MAR-ftns.R")
library(MASS)
source('CLS.R')

nspecies <- 3

sigsq1 <- 0.05; sig12  <- 0.005; sig13 <- sig12;
sig21  <- sig12; sigsq2 <- 0.05; sig23 <- sig13;
sig31 <- sig13; sig32 <- sig13; sigsq3 <- 0.05; 

sig	<-	matrix(c(sigsq1,sig12,sig13,sig21,sigsq2,sig23,sig31,sig32,sigsq3),nrow=3,ncol=3,byrow=T)

#####################---------Example 1: weak intra and interspecific interactions everywhere ------################

c11  <- 0.75 ;  c12 <- -0.06; c13 <- 0.04;
c21 <-  -0.1;  c22  <- 0.75; c23 <- -0.05;
c31 <- 0.07;   c32 <- -0.02;   c33 <- 0.75;

true.B	<-	matrix(c(c11,c12,c13,c21,c22,c23,c31,c32,c33),ncol=3,nrow=3,byrow=T) 

a1     <-  1.9#1.78;#(c11-1)*log(3)
a2     <-  1.3#1.5;#(c22-1)*log(2)
a3     <-  1.1#1.2;#(c33-1)*log(2)
Avec   <- matrix(c(a1,a2,a3), nrow=3,ncol=1)

Bmat1   <- matrix(c(c11, c12,c13,c21,c22,c23,c31,c32,c33), nrow=3, ncol=3, byrow=TRUE); #print(eigen(Bmat)$val)
Sig <- matrix(c(sigsq1, sig12,sig13,sig21,sigsq2,sig23,sig31,sig32,sigsq3), nrow=3,ncol=3, byrow=TRUE)
tausq <- 0.2315
true.parms <- c(Avec, as.vector(Bmat1), as.vector(Sig), 1/sqrt(tausq)) # length=22
mu.vec <- exp(solve(diag(nspecies)-Bmat1)%*%Avec)

#############------------------- Simulating data -------------############################################
# Setting up parameters for simulation
lensim <- 70;
nreps <- 1;
ts.mat.t <- Pspp.marZLN(A=Avec,B=Bmat1,Sigma=Sig,tausq=tausq,ndetect=0.25, len=lensim,nreps=nreps, plot.it=F)$Xmat
#colnames(ts.mat.rep)	<-	1:70

#ts.mat	<- GSSkalman.3reps(ts.mats=ts.mat.rep,K=10,my.thin=5,nchains=2,niter=10000,nadapt=1000)

ts.mat	<-	exp(ts.mat.t)
prop.mat	<-	matrix(NA,ncol=lensim,nrow=3)
dev	<-	colSums(ts.mat)

		for(j in 1:3){
	prop.mat[j,]	<-	ts.mat[j,]/dev
	}
	
	
#colnames(prop.mat)	<-	1:70
#prop.mat	<-	GSSkalman.3reps(ts.mats=prop.mat.rep,K=10,my.thin=5,nchains=2,niter=10000,nadapt=1000)

#####################---------Example 2: strong intraspecific interactions all over, weak inter-specific elsewhere------################

c11  <- 0.01 ;  c12 <- -0.06; c13 <- 0.04;
c21 <-  -0.01; c22  <- 0.01; c23 <- -0.05;
c31 <- 0.07;   c32 <- -0.02;   c33 <- 0.01; 

true.B2	<-	matrix(c(c11,c12,c13,c21,c22,c23,c31,c32,c33),ncol=3,nrow=3,byrow=T)

a1     <-  1.9;#(c11-1)*log(3)
a2     <-  1.3;#(c22-1)*log(2)
a3     <-  1.1;#(c33-1)*log(2)
Avec   <- matrix(c(a1,a2,a3), nrow=3,ncol=1)
Bmat2   <- matrix(c(c11, c12,c13,c21,c22,c23,c31,c32,c33), nrow=3, ncol=3, byrow=TRUE); #print(eigen(Bmat)$val)
Sig <- matrix(c(sigsq1, sig12,sig13,sig21,sigsq2,sig23,sig31,sig32,sigsq3), nrow=3,ncol=3, byrow=TRUE)
tausq <- 0.2315
true.parms <- c(Avec, as.vector(Bmat2), as.vector(Sig), 1/sqrt(tausq)) # length=22
mu.vec <- exp(solve(diag(nspecies)-Bmat2)%*%Avec)
#############------------------- Simulating data -------------############################################
# Setting up parameters for simulation
lensim <- 70;
nreps <- 1;
ts.mat2.t <- Pspp.marZLN(A=Avec,B=Bmat2,Sigma=Sig,tausq=tausq,ndetect=0.25, len=lensim,nreps=nreps, plot.it=F)$Xmat
#colnames(ts.mat2.rep)	<-	1:70
#ts.mat2	<-	GSSkalman.3reps(ts.mats=ts.mat2.rep,K=10,my.thin=5,nchains=2,niter=10000,nadapt=1000)

ts.mat2	<-	exp(ts.mat2.t)
prop.mat2	<-	matrix(NA,ncol=lensim,nrow=3)
dev2	<-	colSums(ts.mat2)

	for(j in 1:3){
		prop.mat2[j,]	<-	ts.mat2[j,]/dev2
		}


#colnames(prop.mat.rep2)	<-	1:70
#prop.mat2	<-	GSSkalman.3reps(ts.mats=prop.mat.rep2,K=10,my.thin=5,nchains=2,niter=10000,nadapt=1000)


#####################---------Example 3: Strong intraspecific in C, low interspecific competitions elswhere------################

c11  <- 0.75 ;  c12 <- -0.06; c13 <- 0.04;
c21 <-  -0.01;  c22  <- 0.75; c23 <- -0.05;
c31 <- 0.07;   c32 <- -0.02;   c33 <- -0.5; 

true.B3	<-	matrix(c(c11,c12,c13,c21,c22,c23,c31,c32,c33),ncol=3,nrow=3,byrow=T)

a1     <-  1.9#1.78;#(c11-1)*log(3)
a2     <-  1.3#1.5;#(c22-1)*log(2)
a3     <-  1.1#1.2;#(c33-1)*log(2)
Avec   <- matrix(c(a1,a2,a3), nrow=3,ncol=1)
Bmat3   <- matrix(c(c11, c12,c13,c21,c22,c23,c31,c32,c33), nrow=3, ncol=3, byrow=TRUE); #print(eigen(Bmat)$val)
Sig <- matrix(c(sigsq1, sig12,sig13,sig21,sigsq2,sig23,sig31,sig32,sigsq3), nrow=3,ncol=3, byrow=TRUE)
tausq <- 0.2315
true.parms <- c(Avec, as.vector(Bmat3), as.vector(Sig), 1/sqrt(tausq)) # length=22
mu.vec <- exp(solve(diag(nspecies)-Bmat3)%*%Avec)
#############------------------- Simulating data -------------############################################
# Setting up parameters for simulation
lensim <- 70;
nreps <- 1;
ts.mat3.t <- Pspp.marZLN(A=Avec,B=Bmat3,Sigma=Sig,tausq=tausq,ndetect=0.25, len=lensim,nreps=nreps, plot.it=FALSE)$Xmat
#colnames(ts.mat3.rep)	<-	1:70
#ts.mat3	<-	GSSkalman.3reps(ts.mats=ts.mat3.rep,K=10,my.thin=5,nchains=2,niter=10000,nadapt=1000)

ts.mat3	<-	exp(ts.mat3.t)
prop.mat3	<-	matrix(NA,ncol=lensim,nrow=3)
dev3	<-	colSums(ts.mat3)
	for(j in 1:3){
		prop.mat3[j,]	<-	ts.mat3[j,]/dev3
	}

	
#colnames(prop.mat.rep3)	<-	1:70
#prop.mat3	<-	GSSkalman.3reps(ts.mats=prop.mat.rep3,K=10,my.thin=5,nchains=2,niter=10000,nadapt=1000)


#####################---------Example 4: High negative interspecific from B to A and from B to C, low intraspecific interactions ------################

#c11  <- 0.55 ;  c12 <- -0.06; c13 <- 0.04;
#c21 <-  -0.6;  c22  <- 0.55; c23 <- -0.75;
#c31 <- 0.07;   c32 <- -0.02;   c33 <- 0.55; 

c11  <- 0.55 ;  c12 <- -0.6; c13 <- 0.07;
c21 <-  -0.06;  c22  <- 0.55; c23 <- -0.02;
c31 <- 0.04;   c32 <- -0.75;   c33 <- 0.55; 



true.B4	<-	matrix(c(c11,c12,c13,c21,c22,c23,c31,c32,c33),ncol=3,nrow=3,byrow=T)

a1     <-  1.9#1.78;#(c11-1)*log(3)
a2     <-  1.3#1.5;#(c22-1)*log(2)
a3     <-  1.1#1.2;#(c33-1)*log(2)
Avec   <- matrix(c(a1,a2,a3), nrow=3,ncol=1)
Bmat4   <- matrix(c(c11, c12,c13,c21,c22,c23,c31,c32,c33), nrow=3, ncol=3, byrow=TRUE); #print(eigen(Bmat)$val)
Sig <- matrix(c(sigsq1, sig12,sig13,sig21,sigsq2,sig23,sig31,sig32,sigsq3), nrow=3,ncol=3, byrow=TRUE)
tausq <- 0.2315
true.parms <- c(Avec, as.vector(Bmat4), as.vector(Sig), 1/sqrt(tausq)) # length=22
mu.vec <- exp(solve(diag(nspecies)-Bmat4)%*%Avec)
#############------------------- Simulating data -------------############################################
# Setting up parameters for simulation
lensim <- 70;
nreps <- 1;
ts.mat4.t <- Pspp.marZLN(A=Avec,B=Bmat4,Sigma=Sig,tausq=tausq,ndetect=1, len=lensim,nreps=nreps, plot.it=FALSE)$Xmat
#colnames(ts.mat4.rep)	<-	1:70
#ts.mat4.rep[which(ts.mat4.rep<0)]<-0
#ts.mat4		<-	GSSkalman.3reps(ts.mats=ts.mat4.rep,K=10,my.thin=5,nchains=2,niter=10000,nadapt=1000)


ts.mat4	<-	exp(ts.mat4.t)
prop.mat4	<-	matrix(NA,ncol=lensim,nrow=3)
dev4	<-	colSums(ts.mat4)
	for(j in 1:3){
		prop.mat4[j,]	<-	ts.mat4[j,]/dev4
	}


#colnames(prop.mat.rep4)	<-	1:70
#prop.mat4	<-	GSSkalman.3reps(ts.mats=prop.mat.rep4,K=10,my.thin=5,nchains=2,niter=10000,nadapt=1000)

####### Stability measurements #########

ex1.stab	<-	stability(Bmat1,sig)
ex2.stab	<-	stability(Bmat2,sig)
ex3.stab	<-	stability(Bmat3,sig)
ex4.stab	<-	stability(Bmat4,sig)


save.image("FourStabilityScenarios.RData")


####### Estimate the parameters #######
#### Example 1

N.mles	<-	mars.cls(t(ts.mat.t))
N.boot	<-	mars.cls.boot(nboot=1000,comm.mat = t(ts.mat.t),cls.fit = N.mles)
upper.B	<-	c(N.boot$B.U[1,],N.boot$B.U[2,],N.boot$B.U[3,])
lower.B	<-	c(N.boot$B.L[1,],N.boot$B.L[2,],N.boot$B.L[3,])

p.mles	<-	mars.cls(t(prop.mat))
p.boot	<-	mars.cls.boot(nboot=1000,comm.mat = t(prop.mat),cls.fit = p.mles)

upper.Bp	<-	c(p.boot$B.U[1,],p.boot$B.U[2,],p.boot$B.U[3,])
lower.Bp	<-	c(p.boot$B.L[1,],p.boot$B.L[2,],p.boot$B.L[3,])

#### Example 2


N.mles2	<-	mars.cls(t(ts.mat2.t))
N.boot2	<-	mars.cls.boot(nboot=1000,comm.mat = t(ts.mat2.t),cls.fit = N.mles2)
upper.B2	<-	c(N.boot2$B.U[1,],N.boot2$B.U[2,],N.boot2$B.U[3,])
lower.B2	<-	c(N.boot2$B.L[1,],N.boot2$B.L[2,],N.boot2$B.L[3,])

p.mles2	<-	mars.cls(t(prop.mat2))
p.boot2	<-	mars.cls.boot(nboot=1000,comm.mat = t(prop.mat2),cls.fit = p.mles2)

upper.Bp2	<-	c(p.boot2$B.U[1,],p.boot2$B.U[2,],p.boot2$B.U[3,])
lower.Bp2	<-	c(p.boot2$B.L[1,],p.boot2$B.L[2,],p.boot2$B.L[3,])

#### Example 3

N.mles3	<-	mars.cls(t(ts.mat3.t))
N.boot3	<-	mars.cls.boot(nboot=1000,comm.mat = t(ts.mat3.t),cls.fit = N.mles3)
upper.B3	<-	c(N.boot3$B.U[1,],N.boot3$B.U[2,],N.boot3$B.U[3,])
lower.B2	<-	c(N.boot3$B.L[1,],N.boot3$B.L[2,],N.boot3$B.L[3,])

p.mles3	<-	mars.cls(t(prop.mat3))
p.boot3	<-	mars.cls.boot(nboot=1000,comm.mat = t(prop.mat3),cls.fit = p.mles3)

upper.Bp3	<-	c(p.boot3$B.U[1,],p.boot3$B.U[2,],p.boot3$B.U[3,])
lower.Bp3	<-	c(p.boot3$B.L[1,],p.boot3$B.L[2,],p.boot3$B.L[3,])

#### Example 4

N.mles4	<-	mars.cls(t(ts.mat4.t))
N.boot4	<-	mars.cls.boot(nboot=1000,comm.mat = t(ts.mat4.t),cls.fit = N.mles4)
upper.B4	<-	c(N.boot4$B.U[1,],N.boot4$B.U[2,],N.boot4$B.U[3,])
lower.B2	<-	c(N.boot4$B.L[1,],N.boot4$B.L[2,],N.boot4$B.L[3,])

p.mles4	<-	mars.cls(t(prop.mat4))
p.boot4	<-	mars.cls.boot(nboot=1000,comm.mat = t(prop.mat4),cls.fit = p.mles4)

upper.Bp4	<-	c(p.boot4$B.U[1,],p.boot4$B.U[2,],p.boot4$B.U[3,])
lower.Bp4	<-	c(p.boot4$B.L[1,],p.boot4$B.L[2,],p.boot4$B.L[3,])


######## -------------Alternative data for figure 2------------------ ##########
# Setting up parameters for simulation
lensim <- 70;
nreps <- 1;
nsims <- 1000

####### ------------- Example 1: --------------- ##################

Bmattots <- matrix(0,nrow=nsims,ncol=9)
Bmatprops <- matrix(0,nrow=nsims,ncol=9)

for(i in 1:1000){
	
	ts.mat.t <- Pspp.marZLN(A=Avec,B=Bmat1,Sigma=Sig,tausq=tausq,ndetect=0.25, len=lensim,nreps=nreps, plot.it=F)$Xmat
	ts.mat	<-	ts.mat.t# exp(ts.mat.t)
	prop.mat	<-	matrix(NA,ncol=lensim,nrow=3)
	dev	<-	colSums(exp(ts.mat))

	for(j in 1:3){
		prop.mat[j,] <-	exp(ts.mat[j,])/dev
	}

	N.mles	<-	mars.cls(t(ts.mat))
	p.mles	<-	mars.cls(t(prop.mat))

	Bmattots[i,] <- as.vector((N.mles$B))
	Bmatprops[i,] <- as.vector((p.mles$B))

}

truthmat1 <- matrix(rep(as.vector(Bmat1),nsims),nrow=nsims,ncol=9,byrow=TRUE)

####### ---------- Example 2 --------------- ####################

Bmattots2 <- matrix(0,nrow=nsims,ncol=9)
Bmatprops2 <- matrix(0,nrow=nsims,ncol=9)

for(i in 1:1000){
	
	ts.mat2.t <- Pspp.marZLN(A=Avec,B=Bmat2,Sigma=Sig,tausq=tausq,ndetect=0.25, len=lensim,nreps=nreps, plot.it=F)$Xmat
	ts.mat2	<-	ts.mat2.t
	prop.mat2	<-	matrix(NA,ncol=lensim,nrow=3)
	dev	<-	colSums(exp(ts.mat2))

	for(j in 1:3){
		prop.mat2[j,] <-	exp(ts.mat2[j,])/dev
	}


	N.mles2	<-	mars.cls(t(ts.mat2))
	p.mles2	<-	mars.cls(t(prop.mat2))

	Bmattots2[i,] <- as.vector((N.mles2$B))
	Bmatprops2[i,] <- as.vector((p.mles2$B))

}

truthmat2 <- matrix(rep(as.vector(Bmat2),nsims),nrow=nsims,ncol=9,byrow=TRUE)

########## ---------- Example 3 ------------- ####################

Bmattots3 <- matrix(0,nrow=nsims,ncol=9)
Bmatprops3 <- matrix(0,nrow=nsims,ncol=9)

for(i in 1:1000){
	
	ts.mat3.t <- Pspp.marZLN(A=Avec,B=Bmat3,Sigma=Sig,tausq=tausq,ndetect=0.25, len=lensim,nreps=nreps, plot.it=F)$Xmat
	ts.mat3	<-	ts.mat3.t
	prop.mat3	<-	matrix(NA,ncol=lensim,nrow=3)
	dev	<-	colSums(exp(ts.mat3))

	for(j in 1:3){
		prop.mat3[j,] <-	exp(ts.mat3[j,])/dev
	}


	N.mles3	<-	mars.cls(t(ts.mat3))
	p.mles3	<-	mars.cls(t(prop.mat3))

	Bmattots3[i,] <- as.vector((N.mles3$B))
	Bmatprops3[i,] <- as.vector((p.mles3$B))

}

truthmat3 <- matrix(rep(as.vector(Bmat3),nsims),nrow=nsims,ncol=9,byrow=TRUE)

######### ------- Example 4 --------------- #################

Bmattots4 <- matrix(0,nrow=nsims,ncol=9)
Bmatprops4 <- matrix(0,nrow=nsims,ncol=9)

for(i in 1:1000){
	
	ts.mat4.t <- Pspp.marZLN(A=Avec,B=Bmat4,Sigma=Sig,tausq=tausq,ndetect=0.25, len=lensim,nreps=nreps, plot.it=F)$Xmat
	ts.mat4	<-	ts.mat4.t
	prop.mat4	<-	matrix(NA,ncol=lensim,nrow=3)
	dev	<-	colSums(exp(ts.mat4))

	for(j in 1:3){
		prop.mat4[j,] <-	exp(ts.mat4[j,])/dev
	}


	N.mles4	<-	mars.cls(t(ts.mat4))
	p.mles4	<-	mars.cls(t(prop.mat4))

	Bmattots4[i,] <- as.vector((N.mles4$B))
	Bmatprops4[i,] <- as.vector((p.mles4$B))

}

truthmat4 <- matrix(rep(as.vector(Bmat4),nsims),nrow=nsims,ncol=9,byrow=TRUE)

save.image("examples_for_plot.RData")