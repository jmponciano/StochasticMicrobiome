########## Functions used saved in BPGSSToolkit.R ##########
source("BPGSSToolkit.R")

# What is the half-life of the density dependent effect?


len1 <- 230;
tau1 <- 130;
a1 <- 3.5;
c1 <- 0.75;
sigmasq1 <- sigmasq2 <- 0.11#0.09725;
a2 <- 0.30;
c2 <- (1/2)^(1/10);
mu1 <- a1/(1-c1);
mu2 <- a2/(1-c2);
t.half <- log(1/2)/log(abs(c2))
mid.thalf <- tau1+t.half/2
mean.mean <- (mu1+mu2)/2

e.scaling <- 1/5

lowleft.m <- mu.vcov.calc(a1 = a1, c1 = c1, sigmasq1 = sigmasq1, a2 = a2, c2 = c2, sigmasq2=sigmasq2, tau=tau1,len=len1);
mean.t <- lowleft.m$mu.vec;
var.t <- lowleft.m$Sigma;
ts.reps <- my.rmvn(n=5, mu.vec=mean.t, cov.mat = var.t)
ts.reps.e <- my.rmvn(n=5, mu.vec=mean.t, cov.mat = var.t*e.scaling)
ts.reps.unstable <- my.rmvn(n=5, mu.vec=mean.t, cov.mat = var.t*3)


# Plot of 5 time series replicates stopping at the breakpoint
tiff("Stable-Unstable.tiff",height=6,width=8,units="in",res=600,compression="lzw",type="cairo")
par(mfrow=c(1,2))

# Stable dynamics
par(bty="n", xaxt="n", yaxt="n", oma=c(1,1,1,1), mar=c(4,4,4,2),mgp=c(2,-2,-2))
matplot(1:tau1, ts.reps[1:tau1,], type="l", col=rep("darkgrey",5), lty=rep(1,5), ylim=c(mu1-2.5, mu1+2.5), xlim=c(-5,len1), xlab="Time", ylab="Log population abundance", cex.lab=0.9, main="Stable", cex.main=0.9)
arrows(x0=c(0,0),y0=c(mu1-2.5,mu1-2.5),x1=c(0,len1),y1=c(1.05*max(ts.reps),mu1-2.5),lwd=1, code=2,length=0.1)
matpoints(1:tau1, ts.reps.e[1:tau1,], type="l", col=rep("blue",5), lty=rep(1,5))

# Dotted line showing mean of first statio pdf
points(1:tau1, rep(mu1,tau1), type="l", lty=2, lwd=1)
text(x=-10,y=mu1, labels=expression(mu[1]),cex=1.05, srt=90)

# pdf of the statio distrib before the break point
sd1 <- sqrt(sigmasq1/(1-c1^2));
supp1 <- seq(mu1-3*sd1, mu1+3*sd1, by=0.01);
dens1 <- 50*dnorm(supp1,mean=mu1,sd=sd1);
supp2 <- seq(mu1-4*(sd1*e.scaling), mu1+4*(sd1*e.scaling), by=0.01);
dens2 <- 25*dnorm(supp2, mean=mu1, sd=e.scaling*sd1)
segments(x0=tau1,y0=mu1-3*sd1,x1=tau1,y1=mu1+3*sd1,col="black",lwd=1);
points(tau1+dens1,supp1,type="l",lty=1,lwd=1,col="darkgrey")
points(tau1+dens2,supp2,type="l",lty=1,lwd=1,col="blue")

# Unstable dynamics
#par(oma=c(1,1,1,2))
matplot(1:tau1, ts.reps.unstable[1:tau1,], type="l", col=rep("darkgrey",5), lty=rep(1,5), ylim=c(mu1-2.3, mu1+2.3), xlim=c(-5,len1), xlab="Time", ylab="Log population abundance", cex.lab=0.9,bty="n", xaxt="n", yaxt="n", main="Unstable", cex.main=0.9)
arrows(x0=c(0,0),y0=c(mu1-2.4,mu1-2.4),x1=c(0,len1),y1=c(1.05*max(ts.reps),mu1-2.4),lwd=1, code=2,length=0.1)
matpoints(1:tau1, ts.reps.e[1:tau1,], type="l", col=rep("blue",5), lty=rep(1,5))

# Dotted line showing mean of first statio pdf
points(1:tau1, rep(mu1,tau1), type="l", lty=2, lwd=1)
text(x=-10,y=mu1, labels=expression(mu[1]),cex=1.05, srt=90)

# pdf of the statio distrib before the break point
sdu <- sqrt(sigmasq1/(1-c1^2))*2.4;
suppu <- seq(mu1-2.4*sdu, mu1+2.4*sdu, by=0.01);
densu <- 50*dnorm(suppu,mean=mu1,sd=sdu);
segments(x0=tau1,y0=mu1-2.4*sdu,x1=tau1,y1=mu1+2.4*sdu,col="black",lwd=1);
points(tau1+densu,suppu,type="l",lty=1,lwd=1,col="darkgrey")
points(tau1+dens2,supp2,type="l",lty=1,lwd=1,col="blue")
dev.off()