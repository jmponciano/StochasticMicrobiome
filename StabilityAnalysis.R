# R files work flow:

# 1. Start with the R file HMP_metadata_updated_taxonomy_StR_CSTs_020619.R
# Run it without running all the lines saving RData files or text files
# The final product from that file is the list "integrated.data2"

# 2. continue with the file "AllKalmanPreds.R" that feeds on the integrated.data2 list

# Stability analysis of the Kalman-estimated time series of abundances
source("CLS2.0.R")

#load("integrateddata2.RData")
nshort <- length(integrated.data2)
mles.list <- list()
stab.mles <- list()

for(i in 1:nshort){
  
  fname <- paste0("longkalman-",names(integrated.data2)[i] ,".txt")
  kalman.mat <- read.table(file=fname, header=TRUE)
  comm.mat <- log(kalman.mat[,-1])
  comm.mat.mles <- mars.cls(comm.mat)
  mles.list[[i]] <- comm.mat.mles
  stab.mles[[i]] <- stability(comm.mat.mles$B, comm.mat.mles$sigma)	
  
}

names(mles.list) <- names(all.Long.KalmanXs)
names(stab.mles) <- names(all.Long.KalmanXs)

Stab.stats <-	matrix(0,ncol=6,nrow=nshort)
n.spp	<-	rep(NA,nshort)
colnames(Stab.stats)	<-	c("Var.prop","Mean.Ret","Var.Ret","React","Mean DDP", "Var DDP")
rownames(Stab.stats)	<-	names(mles.list)

for(i in 1:nshort){
  B.tot   <- mles.list[[i]]$B;
  ave.diag <- mean(diag(B.tot));
  var.diag <- var(diag(B.tot));
  stab.wom	<-	stab.mles[[i]];
  Stab.stats[i,1]	<-	stab.wom[[1]];
  Stab.stats[i,2:3]<-	stab.wom[[2]];
  Stab.stats[i,4]	<-	stab.wom[[3]];
  Stab.stats[i,5] <- ave.diag;
  Stab.stats[i,6] <- var.diag;
  n.spp[i]		<-	nrow(B.tot)
}

write.table(Stab.stats,file="Stabstats-May2019.txt", sep="	", row.names=TRUE)


######################  PCA #########################################
# Data for cluster according to stability measures:
stab.clust.data <- Stab.stats[,1:4]
Xnp <- stab.clust.data;
dimX <- dim(Xnp);
n    <- dimX[1];
p    <- dimX[2];
S    <- var(Xnp);
Dm12       <- matrix(0,nrow=p,ncol=p);
diag(Dm12) <- 1/sqrt(diag(S));
R          <- Dm12%*%S%*%Dm12;
x.bar      <- apply(Xnp,2,mean);
Znp        <- matrix(0,nrow=n,ncol=p);
rownames(Znp) <- rownames(stab.clust.data);
colnames(Znp) <- colnames(stab.clust.data)
for(i in 1:n){
  Znp[i,] <- Dm12%*%(Xnp[i,]-x.bar);	
}

eigs.R <- eigen(R);
princomp.scoresR <- matrix(0,nrow=n,ncol=p); # Componentes principales ('scores')
corr.pcompi.varkR <- matrix(0,nrow=p,ncol=p);
prop.varsR        <- rep(0,p);
for(i in 1:p){
  princomp.scoresR[,i]  <- Znp%*%eigs.R$vectors[,i];
  corr.pcompi.varkR[,i] <- sqrt(eigs.R$values[i])*eigs.R$vectors[,i];
  prop.varsR[i]         <- eigs.R$values[i]/p;
}
# Row and col names of the resulting matrices:
names.corrs <- list(colnames(Xnp),paste0("Princomp ",1:p));
names.scores <- list(row.names(Xnp),paste0("Princomp ",1:p));
row.names(R) <- colnames(Xnp);
colnames(R) <- colnames(Xnp);
dimnames(corr.pcompi.varkR) <- names.corrs;
dimnames(princomp.scoresR) <- names.scores;

pos.pos	<-	which(princomp.scoresR[,1]>0&princomp.scoresR[,2]>0)
pos.neg	<-	which(princomp.scoresR[,1]>0&princomp.scoresR[,2]<0)
rest.pca<-	which(princomp.scoresR[,1]<0)

scaled.scores1 <- princomp.scoresR[,1]/sqrt(sum(princomp.scoresR[,1]^2))
scaled.scores2 <- princomp.scoresR[,2]/sqrt(sum(princomp.scoresR[,2]^2))

cumsum(prop.varsR) # 85% at two dims, that's good
# [1] 0.4432911 0.8534923 0.9438577 1.0000000
#> corr.pcompi.varkR
#         Princomp 1 Princomp 2 Princomp 3 Princomp 4
#Var.prop -0.5282259  0.7553914 -0.3218372 -0.2162915
#Mean.Ret -0.6985602 -0.5999748 -0.3190253  0.2242025
#Var.Ret  -0.7561569 -0.5384687  0.2610437 -0.2648291
#React     0.6590768 -0.6482825 -0.2965833 -0.2395536


pca.scores <- prop.varsR[1]*princomp.scoresR[,1] + prop.varsR[2]*princomp.scoresR[,2] + prop.varsR[3]*princomp.scoresR[,3] + prop.varsR[4]*princomp.scoresR[,4]

##  this plot shows that the Mean and Variance return times 
##  are the most correlated with the overall PCA score (can be seen from corr.pcompi.varkR)
par(mfrow=c(2,2))
plot(pca.scores, Znp[,1], pch=19)
plot(pca.scores, Znp[,2], pch=19)
plot(pca.scores, Znp[,3], pch=19)
plot(pca.scores, Znp[,4], pch=19)

### This plot shows each stability measure  ploted as a function of 
### the PC axis scores with which it has the highest correlation (from corr.pcompi.varkR)
par(mfrow=c(2,2))
plot(princomp.scoresR[,2],Znp[,1], pch=19,main="Std. Variance proportion")
plot(princomp.scoresR[,1],Znp[,2], pch=19,main="Std. Mean Return Time")
plot(princomp.scoresR[,1],Znp[,3], pch=19,main="Std. Var. Return Time")
plot(princomp.scoresR[,1],Znp[,4], pch=19,main="Std. Reactivity")



################### Doing the PCA plot
short.labs <- paste0("W-",substring(row.names(Xnp), first=6,last=nchar(row.names(Xnp))))

tiff("PCA-88women-May2019-StrengthDDPinPCA.tiff", width=10,height=6, units="in", res=600, compression="lzw",
     type="cairo", family="times")
ones <- matrix(1,nrow=3,ncol=3)
twos <- matrix(2,nrow=3,ncol=3)
threes <- matrix(3,nrow=3,ncol=3)
fours <- matrix(4,nrow=3,ncol=9)	
matlay <- rbind(cbind(ones,twos,threes),fours)
layout(mat=matlay)
par(mar=c(3,3,3,1),  oma=c(2,2,2,1), mgp=c(2,0.5,0))
## Zoomed-out figure
plot(scaled.scores1,scaled.scores2, pch=19, main="PCA, full scale",xlab="P.C. I",ylab="P. C. II", ylim=c(-0.65, max(scaled.scores2)),bty="n")
abline(h=0, lty=3);abline(v=0, lty=3)
scale <- sqrt(2/p) #sqrt(d/p), where d is the num. of dimensions plotted
palette("default")
scaled.corrs1 <- 0.1*corr.pcompi.varkR[,1]/sqrt(sum(corr.pcompi.varkR[,1]^2))
scaled.corrs2 <- 0.2*corr.pcompi.varkR[,2]/sqrt(sum(corr.pcompi.varkR[,2]^2))

for (i in 1:4){
  arrows(x0=0,y0=0, x1= scaled.corrs1[i],y1=scaled.corrs2[i], col=i,length=0.10, lwd=2)
}
legend(-0.65,0.6, legend=c("Var. prop.", "Mean Ret.", "Var. Ret.", "Reactivity"), col=1:4, lty=1,lwd=2,cex=0.95, bty="n")
some.labels <- 1:nshort 
some.labels.names <- short.labs[some.labels]
text(scaled.scores1[some.labels],scaled.scores2[some.labels],some.labels.names,pos=1, cex=0.85) 

## First Zoomed-in figure

plot(scaled.scores1,scaled.scores2, pch=19, main="PCA, small scale",xlab="P.C. I",ylab="P. C. II", ylim=c(-0.15, 0.15),xlim=c(-0.1,0.10),bty="n")
abline(h=0, lty=3);abline(v=0, lty=3)
scale <- sqrt(2/p) #sqrt(d/p), where d is the num. of dimensions plotted
palette("default")
scaled.corrs1 <- 0.1*corr.pcompi.varkR[,1]/sqrt(sum(corr.pcompi.varkR[,1]^2))
scaled.corrs2 <- 0.2*corr.pcompi.varkR[,2]/sqrt(sum(corr.pcompi.varkR[,2]^2))

for (i in 1:4){
  arrows(x0=0,y0=0, x1= scaled.corrs1[i],y1=scaled.corrs2[i], col=i,length=0.10, lwd=2)
}
legend("topright", legend=c("Var. prop.", "Mean Ret.", "Var. Ret.", "Reactivity"), col=1:4, lty=1,lwd=2,cex=0.95, bty="n")
some.labels <- 1:nshort 
some.labels.names <- short.labs[some.labels]
text(scaled.scores1[some.labels],scaled.scores2[some.labels],some.labels.names,pos=1, cex=0.85) 

## Second Zoomed-in figure
plot(scaled.scores1,scaled.scores2, pch=19, main="PCA, smaller scale",xlab="P.C. I",ylab="P. C. II", ylim=c(-0.05, 0.05),xlim=c(-0.02,0.055),bty="n")
abline(h=0, lty=3);abline(v=0, lty=3)
scale <- sqrt(2/p) #sqrt(d/p), where d is the num. of dimensions plotted
palette("default")
scaled.corrs1 <- 0.04*corr.pcompi.varkR[,1]/sqrt(sum(corr.pcompi.varkR[,1]^2))
scaled.corrs2 <- 0.04*corr.pcompi.varkR[,2]/sqrt(sum(corr.pcompi.varkR[,2]^2))

for (i in 1:4){
  arrows(x0=0,y0=0, x1= scaled.corrs1[i],y1=scaled.corrs2[i], col=i,length=0.10, lwd=2)
}
legend(0.03,0.05, legend=c("Var. prop.", "Mean Ret.", "Var. Ret.", "Reactivity"), col=1:4, lty=1,lwd=2,cex=0.95, bty="n")
some.labels <- 1:nshort 
some.labels.names <- short.labs[some.labels]
text(scaled.scores1[some.labels],scaled.scores2[some.labels],some.labels.names,pos=1, cex=0.85) 

# And the boxplot of some women
labs2look <- c("W-55","W-121", "W-13","W-27","W-26","W-58","W-17","W-88","W-7", "W-112","W-29", "W-22", "W-122", "W-116", "W-31", "W-11", "W-28", "W-96","W-134", "W-66", "W-82", "W-76", "W-92", "W-101", "W-71","W-60")
nlabs2look <- length(labs2look)

plot(0,0, type="n", xlim=c(1,nlabs2look), xaxt="n", ylab="Strength of density-dependence", xlab="Sample subjects to examine with PCA", bty="l", cex.lab=1.25,cex.axis=0.8)
polygon(x=c(0,(nlabs2look+1),(nlabs2look+1),0), y=c(-0.5,-0.5,0.5,0.5), col="pink", border=FALSE)

for(i in 1:nlabs2look){
  j <- which(short.labs==labs2look[i],arr.ind=TRUE)
  jth.B <- as.vector(diag(mles.list[[j]]$B))
  boxplot(jth.B, at=i, add=TRUE, axes=FALSE, outline=FALSE, boxwex=0.9,col="grey")
}
axis(1,1:nlabs2look, labels=labs2look, tick=TRUE, cex.lab=1.25,cex.axis=0.8)

dev.off()




####  PCA plots with the PCA scores and strength of density dependence of Lactobacillus (L) and non-Lactobacillus (NL)

tiff("PCA-88women-May2019-BothStrengthsDDPinPCA.tiff", width=10,height=6, units="in", res=600, compression="lzw",
     type="cairo", family="times")
ones <- matrix(1,nrow=3,ncol=3)
twos <- matrix(2,nrow=3,ncol=3)
threes <- matrix(3,nrow=3,ncol=3)
fours <- matrix(4,nrow=3,ncol=9)	
matlay <- rbind(cbind(ones,twos,threes),fours)
layout(mat=matlay)
par(mar=c(3,3,3,1),  oma=c(2,2,2,1), mgp=c(2,0.5,0))
## Zoomed-out figure
plot(scaled.scores1,scaled.scores2, pch=19, main="PCA, full scale",xlab="P.C. I",ylab="P. C. II", ylim=c(-0.65, max(scaled.scores2)),bty="n")
abline(h=0, lty=3);abline(v=0, lty=3)
scale <- sqrt(2/p) #sqrt(d/p), where d is the num. of dimensions plotted
palette("default")
scaled.corrs1 <- 0.1*corr.pcompi.varkR[,1]/sqrt(sum(corr.pcompi.varkR[,1]^2))
scaled.corrs2 <- 0.2*corr.pcompi.varkR[,2]/sqrt(sum(corr.pcompi.varkR[,2]^2))

for (i in 1:4){
  arrows(x0=0,y0=0, x1= scaled.corrs1[i],y1=scaled.corrs2[i], col=i,length=0.10, lwd=2)
}
legend(-0.65,0.6, legend=c("Var. prop.", "Mean Ret.", "Var. Ret.", "Reactivity"), col=1:4, lty=1,lwd=2,cex=0.95, bty="n")
some.labels <- 1:nshort 
some.labels.names <- short.labs[some.labels]
text(scaled.scores1[some.labels],scaled.scores2[some.labels],some.labels.names,pos=1, cex=0.85) 

## First Zoomed-in figure

plot(scaled.scores1,scaled.scores2, pch=19, main="PCA, small scale",xlab="P.C. I",ylab="P. C. II", ylim=c(-0.15, 0.15),xlim=c(-0.1,0.10),bty="n")
abline(h=0, lty=3);abline(v=0, lty=3)
scale <- sqrt(2/p) #sqrt(d/p), where d is the num. of dimensions plotted
palette("default")
scaled.corrs1 <- 0.1*corr.pcompi.varkR[,1]/sqrt(sum(corr.pcompi.varkR[,1]^2))
scaled.corrs2 <- 0.2*corr.pcompi.varkR[,2]/sqrt(sum(corr.pcompi.varkR[,2]^2))

for (i in 1:4){
  arrows(x0=0,y0=0, x1= scaled.corrs1[i],y1=scaled.corrs2[i], col=i,length=0.10, lwd=2)
}
legend("topright", legend=c("Var. prop.", "Mean Ret.", "Var. Ret.", "Reactivity"), col=1:4, lty=1,lwd=2,cex=0.95, bty="n")
some.labels <- 1:nshort 
some.labels.names <- short.labs[some.labels]
text(scaled.scores1[some.labels],scaled.scores2[some.labels],some.labels.names,pos=1, cex=0.85) 

## Second Zoomed-in figure
plot(scaled.scores1,scaled.scores2, pch=19, main="PCA, smaller scale",xlab="P.C. I",ylab="P. C. II", ylim=c(-0.05, 0.05),xlim=c(-0.02,0.055),bty="n")
abline(h=0, lty=3);abline(v=0, lty=3)
scale <- sqrt(2/p) #sqrt(d/p), where d is the num. of dimensions plotted
palette("default")
scaled.corrs1 <- 0.04*corr.pcompi.varkR[,1]/sqrt(sum(corr.pcompi.varkR[,1]^2))
scaled.corrs2 <- 0.04*corr.pcompi.varkR[,2]/sqrt(sum(corr.pcompi.varkR[,2]^2))

for (i in 1:4){
  arrows(x0=0,y0=0, x1= scaled.corrs1[i],y1=scaled.corrs2[i], col=i,length=0.10, lwd=2)
}
legend(0.03,0.05, legend=c("Var. prop.", "Mean Ret.", "Var. Ret.", "Reactivity"), col=1:4, lty=1,lwd=2,cex=0.95, bty="n")
some.labels <- 1:nshort 
some.labels.names <- short.labs[some.labels]
text(scaled.scores1[some.labels],scaled.scores2[some.labels],some.labels.names,pos=1, cex=0.85) 

# And the boxplot of some women
labs2look <- c("W-55","W-121", "W-13","W-27","W-26","W-58","W-17","W-88","W-7", "W-112","W-29", "W-22", "W-122", "W-116", "W-31", "W-11", "W-28", "W-96","W-134", "W-66", "W-82", "W-76", "W-92", "W-101", "W-71","W-60")
nlabs2look <- length(labs2look)
some.pcascores <- rep(0, nlabs2look)
for(i in 1:nlabs2look){
  j <- which(short.labs==labs2look[i],arr.ind=TRUE)
  some.pcascores[i] <- pca.scores[j]
}

ordered.inds <- order(some.pcascores,decreasing=TRUE)
labs2look <- labs2look[ordered.inds]
some.pcascores <- some.pcascores[ordered.inds]


#par(mar=c(3,3,3,1),  oma=c(2,2,2,1), mgp=c(2,0.5,0))
plot(0,0, type="n", xlim=c(1,nlabs2look+1), xaxt="n", ylab="Strength of density-dependence", xlab="Sample subjects to examine with PCA", bty="l", cex.lab=1.25,cex.axis=0.8,ylim=c(-1,1.5))
polygon(x=c(0,(nlabs2look+1),(nlabs2look+1),0), y=c(-0.5,-0.5,0.5,0.5), col="pink", border=FALSE)


for(i in 1:nlabs2look){
  j <- which(short.labs==labs2look[i],arr.ind=TRUE)
  Btest <- mles.list[[j]]$B
  jth.B <- as.vector(diag(Btest))
  Lacto.cols <- which(grepl("Lactobacillus", colnames(Btest))==TRUE, arr.ind=TRUE)
  NonLacto.cols <- (1:ncol(Btest))[-Lacto.cols]
  boxplot(jth.B[Lacto.cols], at=i-0.125, add=TRUE, axes=FALSE, outline=FALSE, boxwex=0.4,col="blue")
  boxplot(jth.B[NonLacto.cols], at=i+0.125, add=TRUE, axes=FALSE, outline=FALSE, boxwex=0.4,col="grey")
}
axis(1,1:nlabs2look, labels=labs2look, tick=TRUE, cex.lab=1.25,cex.axis=0.8)
abline(h=1,lty=2,lwd=1.5)
legend("bottomleft",legend=c("Lactobacillus", "Non-Lactobacillus"),col=c("blue","grey"),lty=1,lwd=2,cex=0.95, bty="n")
text(x=1:nlabs2look,y=rep(1.4,nlabs2look), round(some.pcascores,digits=2), cex=0.8)

dev.off()


####  PCA plots with the PCA scores and interaction strength of of Lactobacillus to non Lactobacillus (L->NL) and non-Lactobacillus to Lactobacillus (NL->L)

tiff("PCA-88women-May2019-CrossInteractStrInPCA.tiff", width=10,height=6, units="in", res=600, compression="lzw",
     type="cairo", family="times")
ones <- matrix(1,nrow=3,ncol=3)
twos <- matrix(2,nrow=3,ncol=3)
threes <- matrix(3,nrow=3,ncol=3)
fours <- matrix(4,nrow=3,ncol=9)	
matlay <- rbind(cbind(ones,twos,threes),fours)
layout(mat=matlay)
par(mar=c(3,3,3,1),  oma=c(2,2,2,1), mgp=c(2,0.5,0))
## Zoomed-out figure
plot(scaled.scores1,scaled.scores2, pch=19, main="PCA, full scale",xlab="P.C. I",ylab="P. C. II", ylim=c(-0.65, max(scaled.scores2)),bty="n")
abline(h=0, lty=3);abline(v=0, lty=3)
scale <- sqrt(2/p) #sqrt(d/p), where d is the num. of dimensions plotted
palette("default")
scaled.corrs1 <- 0.1*corr.pcompi.varkR[,1]/sqrt(sum(corr.pcompi.varkR[,1]^2))
scaled.corrs2 <- 0.2*corr.pcompi.varkR[,2]/sqrt(sum(corr.pcompi.varkR[,2]^2))

for (i in 1:4){
  arrows(x0=0,y0=0, x1= scaled.corrs1[i],y1=scaled.corrs2[i], col=i,length=0.10, lwd=2)
}
legend(-0.65,0.6, legend=c("Var. prop.", "Mean Ret.", "Var. Ret.", "Reactivity"), col=1:4, lty=1,lwd=2,cex=0.95, bty="n")
some.labels <- 1:nshort 
some.labels.names <- short.labs[some.labels]
text(scaled.scores1[some.labels],scaled.scores2[some.labels],some.labels.names,pos=1, cex=0.85) 

## First Zoomed-in figure

plot(scaled.scores1,scaled.scores2, pch=19, main="PCA, small scale",xlab="P.C. I",ylab="P. C. II", ylim=c(-0.15, 0.15),xlim=c(-0.1,0.10),bty="n")
abline(h=0, lty=3);abline(v=0, lty=3)
scale <- sqrt(2/p) #sqrt(d/p), where d is the num. of dimensions plotted
palette("default")
scaled.corrs1 <- 0.1*corr.pcompi.varkR[,1]/sqrt(sum(corr.pcompi.varkR[,1]^2))
scaled.corrs2 <- 0.2*corr.pcompi.varkR[,2]/sqrt(sum(corr.pcompi.varkR[,2]^2))

for (i in 1:4){
  arrows(x0=0,y0=0, x1= scaled.corrs1[i],y1=scaled.corrs2[i], col=i,length=0.10, lwd=2)
}
legend("topright", legend=c("Var. prop.", "Mean Ret.", "Var. Ret.", "Reactivity"), col=1:4, lty=1,lwd=2,cex=0.95, bty="n")
some.labels <- 1:nshort 
some.labels.names <- short.labs[some.labels]
text(scaled.scores1[some.labels],scaled.scores2[some.labels],some.labels.names,pos=1, cex=0.85) 

## Second Zoomed-in figure
plot(scaled.scores1,scaled.scores2, pch=19, main="PCA, smaller scale",xlab="P.C. I",ylab="P. C. II", ylim=c(-0.05, 0.05),xlim=c(-0.02,0.055),bty="n")
abline(h=0, lty=3);abline(v=0, lty=3)
scale <- sqrt(2/p) #sqrt(d/p), where d is the num. of dimensions plotted
palette("default")
scaled.corrs1 <- 0.04*corr.pcompi.varkR[,1]/sqrt(sum(corr.pcompi.varkR[,1]^2))
scaled.corrs2 <- 0.04*corr.pcompi.varkR[,2]/sqrt(sum(corr.pcompi.varkR[,2]^2))

for (i in 1:4){
  arrows(x0=0,y0=0, x1= scaled.corrs1[i],y1=scaled.corrs2[i], col=i,length=0.10, lwd=2)
}
legend(0.03,0.05, legend=c("Var. prop.", "Mean Ret.", "Var. Ret.", "Reactivity"), col=1:4, lty=1,lwd=2,cex=0.95, bty="n")
some.labels <- 1:nshort 
some.labels.names <- short.labs[some.labels]
text(scaled.scores1[some.labels],scaled.scores2[some.labels],some.labels.names,pos=1, cex=0.85) 

# And the boxplot of some women
labs2look <- c("W-55","W-121", "W-13","W-27","W-26","W-58","W-17","W-88","W-7", "W-112","W-29", "W-22", "W-122", "W-116", "W-31", "W-11", "W-28", "W-96","W-134", "W-66", "W-82", "W-76", "W-92", "W-101", "W-71","W-60")
nlabs2look <- length(labs2look)
some.pcascores <- rep(0, nlabs2look)
for(i in 1:nlabs2look){
  j <- which(short.labs==labs2look[i],arr.ind=TRUE)
  some.pcascores[i] <- pca.scores[j]
}

ordered.inds <- order(some.pcascores,decreasing=TRUE)
labs2look <- labs2look[ordered.inds]
some.pcascores <- some.pcascores[ordered.inds]


#par(mar=c(3,3,3,1),  oma=c(2,2,2,1), mgp=c(2,0.5,0))
plot(0,0, type="n", xlim=c(1,nlabs2look+1), xaxt="n", ylab="Strength of inter-specific effect", xlab="Sample subjects to examine with PCA", bty="l", cex.lab=1.25,cex.axis=0.8,ylim=c(-1.5,2))
polygon(x=c(0,(nlabs2look+1),(nlabs2look+1),0), y=c(-1.0,-1.0,0,0), col="pink", border=FALSE)


for(i in 1:nlabs2look){
  j <- which(short.labs==labs2look[i],arr.ind=TRUE)
  Btest <- mles.list[[j]]$B
  Lacto.cols <- which(grepl("Lactobacillus", colnames(Btest))==TRUE, arr.ind=TRUE)
  NonLacto.cols <- (1:ncol(Btest))[-Lacto.cols]
  L2NLvec <- as.vector(Btest[NonLacto.cols,Lacto.cols])
  NL2Lvec <- as.vector(Btest[Lacto.cols,NonLacto.cols])
  boxplot(L2NLvec, at=i-0.125, add=TRUE, axes=FALSE, outline=FALSE, boxwex=0.4,col="blue")
  boxplot(NL2Lvec, at=i+0.125, add=TRUE, axes=FALSE, outline=FALSE, boxwex=0.4,col="grey")
}
axis(1,1:nlabs2look, labels=labs2look, tick=TRUE, cex.lab=1.25,cex.axis=0.8)
abline(h=1,lty=2,lwd=1.5)
legend("bottomleft",legend=c("Lactobacillus on Non-Lactobacillus", "Non-Lactobacillus on Lactobacillus"),col=c("blue","grey"),lty=1,lwd=2,cex=0.95, bty="n")
text(x=1:nlabs2look,y=rep(1.4,nlabs2look), round(some.pcascores,digits=2), cex=0.8)

dev.off()


### Plotting only the intra and the inter specific competition effects in two rows
# And the boxplot of some women
labs2look <- c("W-55","W-121", "W-13","W-27","W-26","W-58","W-17","W-88","W-7", "W-112","W-29", "W-22", "W-122", "W-116", "W-31", "W-11", "W-28", "W-96","W-134", "W-66", "W-82", "W-76", "W-92", "W-101", "W-71","W-60")
some.pcascores <- rep(0, nlabs2look)
for(i in 1:nlabs2look){
  j <- which(short.labs==labs2look[i],arr.ind=TRUE)
  some.pcascores[i] <- pca.scores[j]
}


tiff("Boxplots-88women-May2019-AllInteractStrInPCA.tiff", width=10,height=6, units="in", res=600, compression="lzw", type="cairo", family="times")
ones <- matrix(1,nrow=3,ncol=9)
twos <- matrix(2,nrow=3,ncol=9)	
matlay <- rbind(ones,twos)
layout(mat=matlay)
par(mar=c(3,3,3,1),  oma=c(2,2,2,1), mgp=c(2,0.5,0))


ordered.inds <- order(some.pcascores,decreasing=TRUE)
labs2look <- labs2look[ordered.inds]
some.pcascores <- some.pcascores[ordered.inds]
nlabs2look <- length(labs2look)

### First row
plot(0,0, type="n", xlim=c(1,nlabs2look+1), xaxt="n", ylab="Strength of density-dependence", xlab="Sample subjects", bty="l", cex.lab=1.25,cex.axis=0.8,ylim=c(-1,1.5))
polygon(x=c(0,(nlabs2look+1),(nlabs2look+1),0), y=c(-0.5,-0.5,0.5,0.5), col="pink", border=FALSE)


for(i in 1:nlabs2look){
  j <- which(short.labs==labs2look[i],arr.ind=TRUE)
  Btest <- mles.list[[j]]$B
  jth.B <- as.vector(diag(Btest))
  Lacto.cols <- which(grepl("Lactobacillus", colnames(Btest))==TRUE, arr.ind=TRUE)
  NonLacto.cols <- (1:ncol(Btest))[-Lacto.cols]
  boxplot(jth.B[Lacto.cols], at=i-0.125, add=TRUE, axes=FALSE, outline=FALSE, boxwex=0.4,col="blue")
  boxplot(jth.B[NonLacto.cols], at=i+0.125, add=TRUE, axes=FALSE, outline=FALSE, boxwex=0.4,col="grey")
  
}
axis(1,1:nlabs2look, labels=labs2look, tick=TRUE, cex.lab=1.25,cex.axis=0.8)
abline(h=1,lty=2,lwd=1)
legend("bottomleft",legend=c("Lactobacillus", "Non-Lactobacillus"),col=c("blue","grey"),lty=1,lwd=2,cex=0.95, bty="n")
text(x=1:nlabs2look,y=rep(1.4,nlabs2look), round(some.pcascores,digits=2), cex=0.8)


### Second row
plot(0,0, type="n", xlim=c(1,nlabs2look+1), xaxt="n", ylab="Strength of inter-specific effect", xlab="Sample subjects", bty="l", cex.lab=1.25,cex.axis=0.8,ylim=c(-1.5,2))
polygon(x=c(0,(nlabs2look+1),(nlabs2look+1),0), y=c(-1.0,-1.0,0,0), col="pink", border=FALSE)

for(i in 1:nlabs2look){
  j <- which(short.labs==labs2look[i],arr.ind=TRUE)
  Btest <- mles.list[[j]]$B
  Lacto.cols <- which(grepl("Lactobacillus", colnames(Btest))==TRUE, arr.ind=TRUE)
  NonLacto.cols <- (1:ncol(Btest))[-Lacto.cols]
  L2NLvec <- as.vector(Btest[NonLacto.cols,Lacto.cols])
  NL2Lvec <- as.vector(Btest[Lacto.cols,NonLacto.cols])
  boxplot(L2NLvec, at=i-0.125, add=TRUE, axes=FALSE, outline=FALSE, boxwex=0.4,col="blue")
  boxplot(NL2Lvec, at=i+0.125, add=TRUE, axes=FALSE, outline=FALSE, boxwex=0.4,col="grey")
  
}
axis(1,1:nlabs2look, labels=labs2look, tick=TRUE, cex.lab=1.25,cex.axis=0.8)
abline(h=0,lty=2,lwd=1)
legend("bottomleft",legend=c("Lactobacillus on Non-Lactobacillus", "Non-Lactobacillus on Lactobacillus"),col=c("blue","grey"),lty=1,lwd=2,cex=0.95, bty="n")
text(x=1:nlabs2look,y=rep(1.4,nlabs2look), round(some.pcascores,digits=2), cex=0.8)

dev.off()




##### The individual women plots with ALL bells and whistles:

ones <- matrix(1,nrow=2,ncol=4)
twos <- matrix(2,nrow=4,ncol=4)
ones.twos <- rbind(ones,twos)
threes <- matrix(3,nrow=3,ncol=2)
fours <- matrix(rep(c(4,5),each=3),nrow=3,ncol=2,byrow=FALSE)
threes.fours <- rbind(threes,fours)
matlay <- cbind(ones.twos,threes.fours)


dens4plot <- density(pca.scores)
xscores <- dens4plot$x
lenx <- length(xscores)
range.x <- quantile(xscores,probs=c(0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1))


for(i in 1:nshort){
  
  Wom.num <- as.numeric(substring(row.names(Xnp)[i], first=6,last=nchar(row.names(Xnp)[i])))
  fname <- paste0("longkalman-woman",Wom.num,".txt")
  kalman.mat <- read.table(file=fname, header=TRUE)
  comm.mat <- log(kalman.mat[,-1])
  tmat <- matrix(rep(kalman.mat[,1],ncol(comm.mat)), ncol=ncol(comm.mat), nrow=length(kalman.mat[,1]),byrow=FALSE)	
  mycols <- mycols.ftn(n=ncol(comm.mat))
  V.prop <- round(Stab.stats[i,1],digits=3)
  MRT    <- round(Stab.stats[i,2], digits=3)
  main.lab <- paste0(short.labs[i], ", PCA score = ",round(pca.scores[i], digits=3))
  
  plotfname <- paste0(short.labs[i],"Allplots.tiff")
  tiff(plotfname, width=8,height=6, units="in", res=600, compression="lzw",type="cairo", family="times")
  
  layout(mat=matlay)
  par(mar=c(3,3,3,1),  oma=c(2,2,2,2), mgp=c(2,0.5,0))
  
  # First plot
  plot(0,0,xlim=range(xscores),ylim=c(0,25), type="n",xaxt="n", yaxt="n", bty="n",xlab="",ylab="")
  points(xscores, dens4plot$y, type="l", lwd=3, col="blue")
  polygon(x=c(0,xscores,0,rev(xscores),0), y=c(0,rep(0,lenx),0,rev(dens4plot$y),0), col="blue")
  abline(v=pca.scores[i],lwd=2,col="red", lty=1)
  axis(1,at=range.x,labels=round(range.x,digits=2), tick=TRUE, cex.axis=0.8)
  legend("topleft", legend=c(paste0("PCA score for ",short.labs[i]),"Density of all PCA scores"), col=c("red","blue"), lty=1, lwd=2,bty="n",cex=1.1)
  
  # Second plot
  matplot(x=tmat,y=comm.mat,type="b", pch=16, col=mycols, xlab="Time (days)", ylab="Log-abundances", cex.lab=1.25, main=main.lab)
  spp.names <- colnames(comm.mat)
  
  # Third plot
  plot(0,0,type="n",axes=FALSE,xlab="", ylab="")
  legend("topright", legend=spp.names, col=mycols, lty=1,lwd=2, bty="n", cex=1.1)
  
  # Fourth and fifth plots
  Btest <- mles.list[[i]]$B
  
  # Just the fourth plot now
  plot(0,0, type="n", xlim=c(0,2), xaxt="n", ylab="Strength of density-dependence", xlab="", bty="l", 
       cex.lab=1.25,cex.axis=0.8,ylim=c(-1,1.5))
  polygon(x=c(0,2,2,0), y=c(-0.5,-0.5,0.5,0.5), col="pink", border=FALSE)
  
  jth.B <- as.vector(diag(Btest))
  Lacto.cols <- which(grepl("Lactobacillus", colnames(Btest))==TRUE, arr.ind=TRUE)
  NonLacto.cols <- (1:ncol(Btest))[-Lacto.cols]
  boxplot(jth.B[Lacto.cols], at=1-0.25, add=TRUE, axes=FALSE, outline=FALSE, boxwex=0.4,col="blue")
  boxplot(jth.B[NonLacto.cols], at=1+0.25, add=TRUE, axes=FALSE, outline=FALSE, boxwex=0.4,col="grey")
  axis(1,1, labels=short.labs[i], tick=TRUE, cex.lab=1.25,cex.axis=0.8)
  abline(h=1,lty=2,lwd=1)
  legend("topleft",legend=c("L", "Non-L"),col=c("blue","grey"),lty=1,lwd=2,cex=0.8, bty="n")
  
  
  # Just the Fifth plot
  plot(0,0, type="n", xlim=c(0,2), xaxt="n", ylab="Strength of inter-specific effects", xlab="", bty="l", 	
       cex.lab=1.25,cex.axis=0.8,ylim=c(-1.5,2))	
  polygon(x=c(0,2,2,0), y=c(-1.0,-1.0,0,0), col="pink", border=FALSE)
  
  Lacto.cols <- which(grepl("Lactobacillus", colnames(Btest))==TRUE, arr.ind=TRUE)
  NonLacto.cols <- (1:ncol(Btest))[-Lacto.cols]
  L2NLvec <- as.vector(Btest[NonLacto.cols,Lacto.cols])
  NL2Lvec <- as.vector(Btest[Lacto.cols,NonLacto.cols])
  boxplot(L2NLvec, at=1-0.25, add=TRUE, axes=FALSE, outline=FALSE, boxwex=0.4,col="blue")
  boxplot(NL2Lvec, at=1+0.25, add=TRUE, axes=FALSE, outline=FALSE, boxwex=0.4,col="grey")
  axis(1,1, labels=short.labs[i], tick=TRUE, cex.lab=1.25,cex.axis=0.8)
  abline(h=0,lty=2,lwd=1)
  legend("topleft",legend=c("L->Non-L", "Non-L->L"),col=c("blue","grey"),lty=1,lwd=2,cex=0.8, bty="n")
  
  dev.off()
  
  
  
}

#### Looking for explanations:
order.inds.scores <- order(pca.scores,decreasing=TRUE)
ordered.stabstats <- cbind(Stab.stats[order.inds.scores,],pca.scores[order.inds.scores])
colnames(ordered.stabstats) <- c(colnames(Stab.stats), "PCA Score")
write.matrix(x=ordered.stabstats, file="OrderedStabStats.txt")



save.image("MobiluncusRedux2.0.RData")


###### Thrashing around to prep for next round of search:
for(i in 1:nshort){
  B.tot   <- mles.list[[i]]$B;
  print(i)
  print("Gardnerella")
  is.Garden <- sum(grepl(pattern="Gardnerella", x=colnames(B.tot)))
  print(is.Garden)
  if(is.Garden==0){print(colnames(B.tot))}
}	

short.labs[34]


for(i in 1:nshort){
  B.tot   <- mles.list[[i]]$B;
  num.lactob <- sum(grepl(pattern="Lactobacillus", x=colnames(B.tot)))
  num.nonlact <- dim(B.tot)[1]-num.lactob
  print(i)
  print(c(num.lactob,num.nonlact))
}	



