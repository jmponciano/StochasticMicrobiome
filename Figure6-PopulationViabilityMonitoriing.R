# Figure 6:  The "Hurricane plot" Population viability monitoring and estimating 
# the temporally varying cances of Lactobacillus persistence
load("Stability4TwoSppData1.0.RData")

### First we need these functions:

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

library(Rfast)


thres <- 0.45#0.25
thres.times <- 0:10#c(3,5,10,20,50,75,100)
len.sim <- max(thres.times)+1
nthreshs <- length(thres.times)
labs.thresh <- paste0("PT",thres.times)
# PFall.mat <- matrix(0,nrow=nshort,ncol=nthreshs)
# colnames(PFall.mat) <- labs.thresh
# row.names(PFall.mat) <- short.labs
nchains <- 100

ith.ind <- 3#46#10#25#
ith.commat <- Two.spp.commats[[ith.ind]]
ith.len <- nrow(ith.commat)
ith.A <- mles.list[[ith.ind]]$A
print(ith.A)
ith.B <- mles.list[[ith.ind]]$B
print(ith.B)
# Now increase the strength of DDP for Lactobacillus
#ith.B[1,1] <- 0.95
#ith.B[1,2] <- -ith.B[1,2]
#ith.B[2,1] <- -ith.B[2,1] 
ith.S <- mles.list[[ith.ind]]$sigma
#ith.muvec <- ginv(diag(2)-ith.B)%*%ith.A

Lfrac.Obstraj <- exp(ith.commat)[,1]/apply(exp(ith.commat),1,sum)



start.times <- seq(30,69,by=10)
tiff("ExtinctProbsSequence-W-3-Figure6.tiff", width=8,height=8, units="in", res=600, compression="lzw",type="cairo", family="times")
par(mfrow=c(2,2),mar=c(3,3,3,1),  oma=c(2,2,2,1), mgp=c(2,0.5,0))

for(m in start.times){
  
  xinit <- ith.commat[m,]
  plot(1:m,Lfrac.Obstraj[1:m], type="b", lwd=2, xlim=c(0,80),ylim=c(0,1), ylab="Lactobacillus fraction", xlab="Time (days)", cex.lab=1.25)
  last.points <- rep(0,nchains)
  for(j in 1:nchains){
    sim.mat <- Twospp.mar(A=ith.A, B=ith.B, Sigma=ith.S,len=len.sim,Xo=xinit,rnd.init=FALSE)
    jth.frac <- sim.mat[,1]/apply(sim.mat,1,sum)
    points((m+thres.times), jth.frac, type="l", col="grey")
    last.points[j] <- jth.frac[len.sim]
  }
  abline(h=0.5, lty=2, lwd=1)		
  
  #mean.proj <- mean(last.points)
  #sd.proj <- sqrt(var(last.points))
  props4proj <- seq(from=0.001,to=0.999, by=0.001)
  mles.4beta <- beta.mle(last.points)	
  alpha.hat <- mles.4beta$param[1]
  beta.hat <- mles.4beta$param[2]
  proj.dens <- dbeta(x=props4proj,shape1=alpha.hat, shape2=beta.hat)	
  #proj.dens <- dnorm(x=props4proj,mean=mean.proj,sd=sd.proj)	
  where.half <- which(props4proj==0.5,arr.ind=TRUE)
  
  props4prob <- props4proj[1:where.half]
  proj.dens2 <- proj.dens[1:where.half]
  nptsproj <- length(props4prob)
  bottom.xs <- seq(from=m+10, to=m+10+proj.dens2[1]*5,by=0.01)
  top.xs <- seq(from=m+10+proj.dens2[nptsproj]*5, to=m+10,by=-0.01)
  poly.xs <- c(rep(m+10, nptsproj),bottom.xs,m+10+proj.dens2*5,top.xs)
  
  nbottom <- length(bottom.xs)
  ntop <- length(top.xs)
  
  bottom.ys <- rep(0,nbottom)
  top.ys <- rep(0.5,ntop)
  
  poly.ys <- c(rev(props4prob),bottom.ys,props4prob,top.ys)
  
  polygon(x=poly.xs,y=poly.ys, col="lightpink",border=NA)
  points(m+10+proj.dens*5,props4proj, type="l", col="blue", lwd=2)	
  
  p.below <- pbeta(q=0.5,shape1=alpha.hat, shape2=beta.hat)
  title(main=paste0("P(Falling under 0.5) = ",signif(p.below,digits=2)), cex=1.5)
}
dev.off()

#matplot(cbind(1:ith.len,1:ith.len), ith.commat,type="b", pch=1, col=mycols.ftn(2), lwd=2, xlab="Days (t)", ylab="log abundances")

###### Only 1 plot, large:
m <- 27

tiff("ExtinctProbsAt-t27-W-3.tiff", width=8,height=8, units="in", res=600, compression="lzw",type="cairo", family="times")
par(mar=c(3,3,3,1),  oma=c(2,2,2,1), mgp=c(2,0.5,0))

xinit <- ith.commat[m,]
plot(1:m,Lfrac.Obstraj[1:m], type="b", lwd=2, xlim=c(0,60),ylim=c(0,1), ylab="Lactobacillus fraction", xlab="Time (days)", cex.lab=1.25)
last.points <- rep(0,nchains)
for(j in 1:nchains){
  sim.mat <- Twospp.mar(A=ith.A, B=ith.B, Sigma=ith.S,len=len.sim,Xo=xinit,rnd.init=FALSE)
  jth.frac <- sim.mat[,1]/apply(sim.mat,1,sum)
  points((m+thres.times), jth.frac, type="l", col="grey")
  last.points[j] <- jth.frac[len.sim]
}
abline(h=0.5, lty=2, lwd=1)		

#mean.proj <- mean(last.points)
#sd.proj <- sqrt(var(last.points))
props4proj <- seq(from=0.001,to=0.999, by=0.001)
mles.4beta <- beta.mle(last.points)	
alpha.hat <- mles.4beta$param[1]
beta.hat <- mles.4beta$param[2]
proj.dens <- dbeta(x=props4proj,shape1=alpha.hat, shape2=beta.hat)	
#proj.dens <- dnorm(x=props4proj,mean=mean.proj,sd=sd.proj)	
where.half <- which(props4proj==0.5,arr.ind=TRUE)

props4prob <- props4proj[1:where.half]
proj.dens2 <- proj.dens[1:where.half]
nptsproj <- length(props4prob)
bottom.xs <- seq(from=m+10, to=m+10+proj.dens2[1]*5,by=0.01)
top.xs <- seq(from=m+10+proj.dens2[nptsproj]*5, to=m+10,by=-0.01)
poly.xs <- c(rep(m+10, nptsproj),bottom.xs,m+10+proj.dens2*5,top.xs)

nbottom <- length(bottom.xs)
ntop <- length(top.xs)

bottom.ys <- rep(0,nbottom)
top.ys <- rep(0.5,ntop)

poly.ys <- c(rev(props4prob),bottom.ys,props4prob,top.ys)

polygon(x=poly.xs,y=poly.ys, col="lightpink",border=NA)
points(m+10+proj.dens*5,props4proj, type="l", col="blue", lwd=2)	

p.below <- pbeta(q=0.5,shape1=alpha.hat, shape2=beta.hat)
title(main=paste0("P(Falling under 0.5) = ",signif(p.below,digits=2)), cex=1.5)
dev.off()


############# Next, applying this concept to all women in the data set:

# Now, computing the probability of falling under clinical threshold

thres <- 0.45#0.25
thres.times <- 2:30#c(3,5,10,20,50,75,100)
len.sim <- max(thres.times)
nthreshs <- length(thres.times)
labs.thresh <- paste0("PT",thres.times)
PFall.mat <- matrix(0,nrow=nshort,ncol=nthreshs)
colnames(PFall.mat) <- labs.thresh
row.names(PFall.mat) <- short.labs
nchains <- 200
NtimesNo <- 1

# Idea: track percent of times pop. dropped below threshold, not only if it did
#       threshold is defined as 45% of the total abundance

for(i in 1:nshort){
  ith.ind <- i#10#25#
  ith.commat <- Two.spp.commats[[ith.ind]]
  ith.len <- nrow(ith.commat)
  ith.A <- mles.list[[ith.ind]]$A
  ith.B <- mles.list[[ith.ind]]$B
  ith.S <- mles.list[[ith.ind]]$sigma
  ith.muvec <- ginv(diag(2)-ith.B)%*%ith.A
  ith.half <- 0.5*sum(ith.muvec)
  #print(exp(ith.muvec))
  ith.xo <- rep(ith.half,2) #ith.muvec #log(NtimesNo) +  ith.commat[ith.len,]
  #ith.initot <- sum(exp(ith.muvec))
  #Nthres <- thres*ith.initot
  for(k in 1:nthreshs){
    code.line1 <- paste0("cross.thres",k," <-rep(0,nchains)")
    eval(parse(text=code.line1))
  }
  
  for(j in 1:nchains){
    sim.mat <- Twospp.mar(A=ith.A, B=ith.B, Sigma=ith.S,len=len.sim,Xo=ith.xo,rnd.init=FALSE)
    #print(sim.mat)
    for(h in 1:nthreshs){
      #			test <- sum((sim.mat[1:thres.times[h],1]/apply(sim.mat[1:thres.times[h],],1,sum))<thres)>0			
      test <- sum((sim.mat[thres.times[h],1]/sum(sim.mat[thres.times[h],]))<thres)>0
      code.line2 <- paste0("cross.thres",h,"[j] <- test")
      eval(parse(text=code.line2))
    }
  }
  for(ell in 1:nthreshs){
    code.line3 <- paste0("PFall.mat[i,",ell,"] <- sum(cross.thres",ell,")/nchains")
    eval(parse(text=code.line3))
  }  
  
}	

Tthres.mat <- matrix(rep(c(1,thres.times),nshort), nrow=nthreshs+1, ncol=nshort,byrow=FALSE)
tPFall.mat <- t(cbind(rep(0,nshort),PFall.mat))

thisdir <- getwd()


############## Figure 7 from the main text:  look under the folder "/Wi-2sppPlots/Stable/PFall0.45-W-3.tiff" and
############## /Wi-2sppPlots/Stable/PFall0.45-W-60.tiff"                                                        

############## To plot results for woman 3 just run the code inside the for loop setting i <- 3

for(i in 1:nshort){
  stab.level <- belongingsF[i]
  plotfname <- paste0(thisdir,"/Wi-2sppPlots/",stab.level,"/PFall",thres,"-", short.labs[i],".tiff")	
  tiff(plotfname, width=6,height=10, units="in", res=600, compression="lzw",type="cairo", family="times")
  par(mar=c(3,3,3,1),  oma=c(2,2,2,1), mgp=c(2,0.5,0))
  matplot(x=Tthres.mat,y=1-tPFall.mat, type="b", pch=1, col=mycols.ftn(nshort), lwd=2, xlab="Days (t)", ylab=paste0("P(Lactobacillus is above threshold)"), 
          cex.lab=1.25, main=paste0(short.labs[i],": Pr(persisting above ",thres*100,"% of total on day 't')"),bty="n")
  points(Tthres.mat[,i], 1-tPFall.mat[,i], type="b", pch=16, col="black", lwd=2)
  dev.off()
}

###  Grab one example, say W-60, and project changes in population parameters to see how the P(Persistence changes)

New.PFallvec <- rep(0,nthreshs)
# Set i= 46 to work with woman 60
ith.ind <- 46#10#25#
ith.commat <- Two.spp.commats[[ith.ind]]
ith.len <- nrow(ith.commat)
ith.A <- mles.list[[ith.ind]]$A
print(ith.A)
ith.B <- mles.list[[ith.ind]]$B
print(ith.B)
# Now increase the strength of DDP for Lactobacillus
#ith.B[1,1] <- 0.95
ith.B[1,2] <- -ith.B[1,2]
ith.B[2,1] <- -ith.B[2,1] 
ith.S <- mles.list[[ith.ind]]$sigma
ith.muvec <- ginv(diag(2)-ith.B)%*%ith.A
ith.half <- 0.5*sum(ith.muvec)
ith.xo <- rep(ith.half,2)
for(k in 1:nthreshs){
  code.line1 <- paste0("cross.thres",k," <-rep(0,nchains)")
  eval(parse(text=code.line1))
}
for(j in 1:nchains){
  sim.mat <- Twospp.mar(A=ith.A, B=ith.B, Sigma=ith.S,len=len.sim,Xo=ith.xo,rnd.init=FALSE)
  #print(sim.mat)
  for(h in 1:nthreshs){
    test <- sum((sim.mat[thres.times[h],1]/sum(sim.mat[thres.times[h],]))<thres)>0
    code.line2 <- paste0("cross.thres",h,"[j] <- test")
    eval(parse(text=code.line2))
  }
}
for(ell in 1:nthreshs){
  code.line3 <- paste0("New.PFallvec[",ell,"] <- sum(cross.thres",ell,")/nchains")
  eval(parse(text=code.line3))
}  

New.PPersist <- c(1,1-New.PFallvec)

rbind(New.PPersist,c(1,1-PFall.mat[46,]))


stab.level <- belongingsF[46]
plotfname <- paste0(thisdir,"/Wi-2sppPlots/",stab.level,"/PFallHyp",thres,"-", short.labs[46],"-swapBsigns.tiff")	
tiff(plotfname, width=6,height=10, units="in", res=600, compression="lzw",type="cairo", family="times")
par(mar=c(3,3,3,1),  oma=c(2,2,2,1), mgp=c(2,0.5,0))
plot(x=Tthres.mat[,1],y=c(1,1-PFall.mat[46,]), type="b", pch=1, col="black", lwd=2, xlab="Days (t)", ylab=paste0("P(Lactobacillus is above threshold)"), 
     cex.lab=1.25, main=paste0(short.labs[46],": Pr(persisting above ",thres*100,"% of total on day 't')"),bty="n")
points(Tthres.mat[,1], New.PPersist, type="b", pch=16, col="lightpink", lwd=2)
legend(x=15,y=0.8, legend=c("Estimated dynamics","Restored dynamics"), pch=c(16,16), col=c("black","lightpink"), cex=1.25, bty="n")
dev.off()

save.image("Stability4TwoSppData4.0.RData") 



