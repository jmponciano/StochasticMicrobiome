###  Figure 8, persistence probabilities for woman-60

load("Stability4TwoSppData4.0.RData")
library("MASS")


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
plotfname <- paste0(thisdir,"/Wi-2sppPlots/",stab.level,"/PFallHyp",thres,"-", short.labs[46],"-Examples4paper2023.tiff")	
tiff(plotfname, width=8,height=10, units="in", res=600, compression="lzw",type="cairo", family="times")
par(mfrow=c(1,2), mar=c(4,4,4,1),  oma=c(2,2,2,1), mgp=c(2,0.5,0))

plot(x=Tthres.mat[,1],y=c(1,1-PFall.mat[46,]), type="b", pch=1, col="black", lwd=2, ylab="", 
main="",bty="n", xlab="")
points(Tthres.mat[,1], New.PPersist, type="b", pch=16, col="lightpink", lwd=2)
legend(x=10,y=0.8, legend=c("Estimated dynamics","Restored dynamics"), pch=c(16,16), col=c("black","lightpink"), cex=1.05, bty="n")
text(x=1.5,y=1.026,labels="(a)", cex=1.35)
#dev.off()

#####  Now increasing the growth rate of lactobacillus

New.PFallvec <- rep(0,nthreshs)
# Set i= 46 to work with woman 60
ith.ind <- 46#10#25#
ith.commat <- Two.spp.commats[[ith.ind]]
ith.len <- nrow(ith.commat)
ith.A <- mles.list[[ith.ind]]$A
print(ith.A)
ith.B <- mles.list[[ith.ind]]$B
print(ith.B)
ith.A[1] <- ith.A[1]*1.25

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
plot(x=Tthres.mat[,1],y=c(1,1-PFall.mat[46,]), type="b", pch=1, col="black", lwd=2, ylab="", 
main="",bty="n", xlab="")
points(Tthres.mat[,1], New.PPersist, type="b", pch=16, col="lightpink", lwd=2)
legend(x=10,y=0.8, legend=c("Estimated dynamics","Restored dynamics"), pch=c(16,16), col=c("black","lightpink"), cex=1.05, bty="n")
text(x=1.5,y=1.026,labels="(b)",cex=1.35)

mtext(text=paste0(short.labs[46],": Pr(persisting above ",thres*100,"% of total on day 't')"), side=3, outer=TRUE, adj=0.5, cex=1.5)

mtext(text="P(Lactobacillus is above threshold)", side=2,outer=TRUE, adj=0.5, cex=1.15)

mtext(text="Days (t) since final observation", side=1,outer=TRUE, adj=0.5, cex=1.15)
dev.off()

