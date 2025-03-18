# Lactobacillus vs. the rest of the world

source("CLS3.0.R") #### Warning!!!! The stability() function in CLS3.0.R outputs 
                   #### results in a slightly different format than the 
                   #### stability() function in the ExampleCalcs folder. The  
                   #### numbers in the output are identical but the way the 
                   #### list of statistics is returned is slightly different.
                   #### Hence, the results from this stability() function are
                   #### read differently. Follow the code in the Examplecalcs
                   #### folder to replicate the paper results.

load("integrateddata2.RData")
nshort <- length(integrated.data2)
mles.list <- list()
stab.mles <- list()
Two.spp.commats <- list()


for(i in 1:nshort){
	
	fname <- paste0("longkalman-",names(integrated.data2)[i] ,".txt")
	kalman.mat <- read.table(file=fname, header=TRUE)
	comm.mat1 <- kalman.mat[,-1]
	Lactocols <- which(grepl("Lactobacillus", colnames(comm.mat1))==TRUE,arr.ind=TRUE)
	NonLactocols <- (1:ncol(comm.mat1))[-Lactocols]
	if(length(Lactocols)==1){Lacto.vec <- log(comm.mat1[,Lactocols])}else{
		Lacto.vec <- log(apply(comm.mat1[,Lactocols],1,sum))}
	if(length(NonLactocols)==1){NonLac.vec <- log(comm.mat1[,NonLactocols])}else{
		NonLac.vec <- log(apply(comm.mat1[,NonLactocols],1,sum))}
	comm.mat <- cbind(Lacto.vec,NonLac.vec)
	colnames(comm.mat) <- c("Lactobacillus", "Non-Lactobacillus")	
	Two.spp.commats[[i]] <- comm.mat
	comm.mat.mles <- mars.cls(comm.mat)
	mles.list[[i]] <- comm.mat.mles
	stab.mles[[i]] <- stability(comm.mat.mles$B, comm.mat.mles$sigma)	
	
}

names(mles.list) <- names(integrated.data2)
names(stab.mles) <- names(integrated.data2)
names(Two.spp.commats) <- names(integrated.data2)


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

#write.table(Stab.stats,file="Stabstats2spp-May2019.txt", sep="	", row.names=TRUE)

######################  PCA #########################################
# Data for cluster according to stability measures:
stab.clust.data <- Stab.stats[,1:4]

write.table(x=stab.clust.data, file="TableS1.txt")
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

write.table(corr.pcompi.varkR, file="TableS2.txt")

pos.pos	<-	which(princomp.scoresR[,1]>0&princomp.scoresR[,2]>0)
pos.neg	<-	which(princomp.scoresR[,1]>0&princomp.scoresR[,2]<0)
rest.pca<-	which(princomp.scoresR[,1]<0)

scaled.scores1 <- princomp.scoresR[,1]/sqrt(sum(princomp.scoresR[,1]^2))
scaled.scores2 <- princomp.scoresR[,2]/sqrt(sum(princomp.scoresR[,2]^2))

cumsum(prop.varsR)
#[1] 0.6034070 0.8482549 0.9974194 1.0000000

#> corr.pcompi.varkR
#         Princomp 1  Princomp 2  Princomp 3   Princomp 4
#Var.prop  0.7291374 -0.09614292  0.67757406  0.002931069
#Mean.Ret  0.9556163 -0.06774068 -0.27780469  0.070945933
#Var.Ret   0.9657255 -0.05357288 -0.24336227 -0.072656147
#React    -0.1901533 -0.98116728 -0.03392658 -0.001218271

pca.scores <- prop.varsR[1]*princomp.scoresR[,1] + prop.varsR[2]*princomp.scoresR[,2] + prop.varsR[3]*princomp.scoresR[,3] + prop.varsR[4]*princomp.scoresR[,4];

#### Looking for explanations:
order.inds.scores <- order(pca.scores,decreasing=TRUE)
ordered.stabstats <- cbind(Stab.stats[order.inds.scores,],pca.scores[order.inds.scores])
colnames(ordered.stabstats) <- c(colnames(Stab.stats), "PCA Score")
write.matrix(x=ordered.stabstats, file="OrderedStabStats2spp.txt")



##  this plot shows that the Variance proportion and the Mean and Variance return times 
##  are the most correlated with the overall PCA score (can be seen from corr.pcompi.varkR)
tiff("StabStatsVsPCAscore.tiff", width=8,height=8, units="in", res=600, compression="lzw",type="cairo", family="times")
par(mfrow=c(2,2),mar=c(3,3,3,1),  oma=c(2,2,2,1), mgp=c(2,0.5,0))
plot(pca.scores, Znp[,1], pch=19, xlab="PCA scores", ylab="Std. Variance proportion", bty="n",cex.lab=1.25)
plot(pca.scores, Znp[,2], pch=19, xlab="PCA scores", ylab="Std. Mean Return Time", bty="n",cex.lab=1.25)
plot(pca.scores, Znp[,3], pch=19, xlab="PCA scores", ylab="Std. Var. Return Time", bty="n",cex.lab=1.25)
plot(pca.scores, Stab.stats[,"Mean DDP"], pch=19, xlab="PCA scores", ylab="Mean Density dependence", bty="n",cex.lab=1.25)
dev.off()


##### Modeling stability metrics as a function of the pca scores-  code not needed for the PCA plot
df4fit <- data.frame(vprop=Znp[-30,1],mrt=Znp[-30,2],vrt=Znp[-30,3],ddp = Stab.stats[-30,"Mean DDP"],pcascore = pca.scores[-30],pcascore2 = pca.scores[-30]^2)
pcasc4fit <- seq(min(pca.scores[-30]),max(pca.scores[-30]),by=0.005)

df4fit <- data.frame(vprop=Znp[,1],mrt=Znp[,2],vrt=Znp[,3],ddp = Stab.stats[,"Mean DDP"],pcascore = pca.scores,pcascore2 = pca.scores^2)
pcasc4fit <- seq(min(pca.scores),max(pca.scores),by=0.005)

df4pred <-  data.frame(pcascore=pcasc4fit,pcascore2=pcasc4fit^2)
lm1 <- lm(vprop~pcascore+pcascore2, data=df4fit)
lm2 <- lm(mrt~pcascore+pcascore2, data=df4fit)
lm3 <- lm(vrt~pcascore+pcascore2, data=df4fit)
lm4 <- lm(ddp~pcascore+pcascore2, data=df4fit)
pred1fit <- predict.lm(lm1,newdata=df4pred,se.fit=TRUE)
pred2fit <- predict.lm(lm2,newdata=df4pred,se.fit=TRUE)
pred3fit <- predict.lm(lm3,newdata=df4pred,se.fit=TRUE)
pred4fit <- predict.lm(lm4,newdata=df4pred,se.fit=TRUE)

pred1fit$se.fit

tiff("StabStatsVsPCAscore-fancy.tiff", width=8,height=8, units="in", res=600, compression="lzw",type="cairo", family="times")
par(mfrow=c(2,2),mar=c(3,3,3,1),  oma=c(2,2,2,1), mgp=c(2,0.5,0))
plot(pca.scores[-30], Znp[-30,1], pch=19, xlab="PCA scores", ylab="Std. Variance proportion", bty="n",cex.lab=1.25,type="n")
yvec1 <- c(pred1fit$fit +1.96*pred1fit$se.fit,rev(pred1fit$fit -1.96*pred1fit$se.fit))
xvec <- c(pcasc4fit,rev(pcasc4fit))
polygon(x=xvec,y=yvec1,col="grey", border="grey")
points(pca.scores[-30], Znp[-30,1], pch=19)
points(pcasc4fit,pred1fit$fit,type="l",col="black",lwd=2)

plot(pca.scores[-30], Znp[-30,2], pch=19, xlab="PCA scores", ylab="Std. Mean Return Time", bty="n",cex.lab=1.25,type="n")
yvec2 <- c(pred2fit$fit +1.96*pred2fit$se.fit,rev(pred2fit$fit -1.96*pred2fit$se.fit))
polygon(x=xvec,y=yvec2,col="grey", border="grey")
points(pca.scores[-30], Znp[-30,2], pch=19)
points(pcasc4fit,pred2fit$fit,type="l",col="black",lwd=2)

plot(pca.scores[-30], Znp[-30,3], pch=19, xlab="PCA scores", ylab="Std. Var. Return Time", bty="n",cex.lab=1.25,type="n")
yvec3 <- c(pred3fit$fit +1.96*pred3fit$se.fit,rev(pred3fit$fit -1.96*pred3fit$se.fit))
polygon(x=xvec,y=yvec3,col="grey", border="grey")
points(pca.scores[-30], Znp[-30,3], pch=19)
points(pcasc4fit,pred3fit$fit,type="l",col="black",lwd=2)

plot(pca.scores[-30], Stab.stats[-30,"Mean DDP"], pch=19, xlab="PCA scores", ylab="Mean Density dependence", bty="n",cex.lab=1.25,type="n")
yvec4 <- c(pred4fit$fit +1.96*pred4fit$se.fit,rev(pred4fit$fit -1.96*pred4fit$se.fit))
polygon(x=xvec,y=yvec4,col="grey", border="grey")
points(pca.scores[-30], Stab.stats[-30,"Mean DDP"], pch=19)
points(pcasc4fit,pred4fit$fit,type="l",col="black",lwd=2)

dev.off()
################### End of modeling stability metrics as a function of the pca scores



#plot(pca.scores, Znp[,4], pch=19,ylim=c(-1,0.3), xlab="PCA scores", ylab="Std. Reactivity")
#plot(pca.scores, Stab.stats[,"Var DDP"],pch=19, xlab="PCA scores", ylab="Var Density dependence")

### This plot shows each stability measure  plotted as a function of 
### the PC axis scores with which it has the highest correlation (from corr.pcompi.varkR)
par(mfrow=c(2,2))
plot(princomp.scoresR[,1],Znp[,1], pch=19,main="Std. Variance proportion")
plot(princomp.scoresR[,1],Znp[,2], pch=19,main="Std. Mean Return Time")
plot(princomp.scoresR[,1],Znp[,3], pch=19,main="Std. Var. Return Time")
plot(princomp.scoresR[,2],Znp[,4], pch=19,main="Std. Reactivity")


#################### This cluster analysis is needed to obtain the stability classification
# Scores for cluster
scores.cluster <- cbind(scaled.scores1,scaled.scores2,pca.scores, mean.ddp=Stab.stats[,"Mean DDP"])
kmeans.clust <- kmeans(x=scores.cluster,centers=4)

cluster.means <- kmeans.clust$centers

write.table(cluster.means, file="TableS3.txt")

# Looking at the correlation of the prin comps with each variable:
#Princomp 1  Princomp 2  Princomp 3   Princomp 4
#Var.prop -0.7291374 -0.09614292  0.67757406  0.002931069
#Mean.Ret -0.9556163 -0.06774068 -0.27780469  0.070945933
#Var.Ret  -0.9657255 -0.05357288 -0.24336227 -0.072656146
#React     0.1901533 -0.98116728 -0.03392658 -0.001218271

# Tells you that: the smaller the value of Var.prop, Mean.Ret, Var.Ret, the
# higher the PC I score, which also means that the more stable a dynamics is.
# PC I explains 60% of the variation
# As well, Reactivity scores are highly negatively correlated with PC II
# so the smaller the reactivity, the more stable the dynamics.  PC II explains 24% 
# of the variations

### Looking at the means of the cluster groups tells you that
### the clusters 1 and 3, have on average positive scaled.scores1 
### and the two lowest values of the density dependent coefficient (hence stronger
### self regulation via intra-specific competition).  The highest the mean score
### in PC I and the lowest the value of the ddp, the more stable the dynamics

# scaled.scores1 scaled.scores2 pca.scores  mean.ddp
# 1     0.22339761     0.02708272  2.0710675 0.5205932
# 2    -0.01622555    -0.01845970 -0.2262513 0.7864128
# 3     0.05608588     0.05714975  0.5858070 0.6901813
# 4    -0.08951403    -0.01986391 -0.7799067 0.9059187

### Toether, the above observations suggest using:
# 1 = very stable
# 2 = unstable
# 3 = stable
# 4 = highly unstable



#names(kmeans.clust)
belongings <- kmeans.clust$cluster
kcols <- mycols.ftn(max(belongings)+4)
cols4plot <- rep(kcols[1],nshort)
for(i in 1:nshort){
	
	cols4plot[i] <- kcols[belongings[i]]
}

belongingsF <- as.factor(belongings)
levels(belongingsF)

#pcascores.belongings.df <- data.frame(scores.cluster, belongings)
#View(pcascores.belongings.df)


#belongings.NumF <- as.factor(belongings)

levels(belongingsF) <- c("Highly stable","Unstable","Stable","Highly unstable") # check that this is so

belongings.df <- data.frame(short.labs,belongings, belongingsF,scores.cluster)

belongings.df[belongings.df$belongingsF=="Highly stable",]
belongings.df[belongings.df$belongingsF=="Stable",]
belongings.df[belongings.df$belongingsF=="Unstable",]
belongings.df[belongings.df$belongingsF=="Highly unstable",]
belongings.df


######## Figure 5 in the November 2, 2023 Manuscript ###############
tiff("PCA-88women-2spp-Nov2023-noWlabs.tiff", width=8,height=8, units="in", res=600, compression="lzw",
	type="cairo", family="times")

par(mar=c(3,3,3,1),  oma=c(2,2,2,1), mgp=c(2,0.5,0))
#short.labs <- paste0("W-",substring(row.names(Xnp), first=6,last=nchar(row.names(Xnp))))
#. ATTENTION! below, you might need to multiply by -1 the x-axis coordinates (it's just
#  a matter of how the function "eigen" solves for the eigenvalues, this sign is
#  arbitrary in the PCA)
plot(scaled.scores1,scaled.scores2, pch=19, main="",xlab="P.C. I",ylab="P. C. II",bty="n", 
     xlim=c(-0.2,0.6), ylim=c(-0.075,0.10), col=cols4plot,cex=1.25)
abline(h=0, lty=3);abline(v=0, lty=3)
scale <- sqrt(2/p) #sqrt(d/p), where d is the num. of dimensions plotted
palette("default")
scaled.corrs1 <- 0.1*corr.pcompi.varkR[,1]/sqrt(sum(corr.pcompi.varkR[,1]^2))
scaled.corrs2 <- 0.05*corr.pcompi.varkR[,2]/sqrt(sum(corr.pcompi.varkR[,2]^2))

for (i in 1:4){
	arrows(x0=0,y0=0, x1= scaled.corrs1[i],y1=scaled.corrs2[i], col=kcols[(i+4)],length=0.10, lwd=2)
}
legend("bottomright", legend=c("Var. prop.", "Mean Ret.", "Var. Ret.", "Reactivity"), col=kcols[5:8],lty=1,lwd=2,cex=1.25, bty="n")
legend(x=0.43,y=-0.01, legend=c("Highly stable","Stable","Unstable","Highly unstable"), 
       col=c(kcols[1],kcols[3], kcols[2],kcols[4]), pch=rep(19,4), cex=1.25, bty="n")

# Uncomment to plot the women labels
#some.labels <- 1:nshort 
#some.labels.names <- short.labs[some.labels]
#text(scaled.scores1[some.labels],scaled.scores2[some.labels],some.labels.names,pos=1, cex=0.8) 

dev.off()


save.image("Stability4TwoSppData1.0.RData") # Saved on Nov 3rd 2023


################  Figure 4 in the Nov 2 manuscript ##############


########### ROUSS SHOW ##############
# 1. Do you want to compute the ML estimates or the REML estimates?
method <- "REML" # alternatively, set method <- "ML"
# 2. Do you want to plot the predictions?
pred.plot <- "FALSE" # Set it to "FALSE" if you do not want to plot the predictions
# 3. Do you want to plot the parametric bootstrap distribution of the estimates?
pboot.plot <- "FALSE" # Set it to "FALSE" if you do not want to plot the bootstrap distribution of the estimates


library("MASS")
source("ROUSS/ROUSSE-2.0.R")
# 4. How many bootstrap replicates?
NBoot <- 1000; # 10 just to try, for formal results use 1000 or more 

samp.tss <- integrated.data2[[25]] # W-30
Time.t <- as.numeric(samp.tss$raw.days)
tt <- Time.t-Time.t[1]
long.t <- tt[1]:max(tt)
ithlens     <- length(Time.t)
ithlong.lens <- length(long.t)
all.lacto <- apply(samp.tss[,9:12],1,sum)
log.obs <- log(all.lacto)

ts.results <- ROUSS.CALCS(Yobs=log.obs,Tvec=Time.t, pmethod=method, nboot=NBoot, plot.pred=pred.plot, plot.bootdists = pboot.plot);   

ts.pboot.cis <- ts.results$pboot.cis
longpred.mle <- ts.results$pboot.preds2[,3]
longpred.lci <- ts.results$pboot.preds2[,2]
longpred.uci <- ts.results$pboot.preds2[,4]

shortpred.mle <- ts.results$pboot.preds1[,3]
shortpred.lci <- ts.results$pboot.preds1[,2]
shortpred.uci <- ts.results$pboot.preds1[,4]

tiff("ExampleKalman.tiff", width=8,height=6, units="in", res=600, compression="lzw",type="cairo", family="times")
par(mar=c(3,3,3,1),  oma=c(2,2,2,1), mgp=c(2,0.5,0))
plot(tt,all.lacto, pch=1,type="n", bty="n", xlab="Time (days)", ylab="Abundance of Lactobacilli")
polygon(x=c(long.t,rev(long.t)), y=c(longpred.lci,rev(longpred.uci)), col="grey", border=NA)
points(tt,all.lacto, pch=1, type="b")
points(long.t,longpred.mle, pch=16, type="b")
legend("topleft", legend=c("Observed", "Estimated without sampling error"), pch=c(1,16),bty="n")
dev.off()

####################. End of Figure 4 in the code. ################


dens4plot <- density(pca.scores)
xscores <- dens4plot$x
lenx <- length(xscores)
range.x <- quantile(xscores,probs=c(0,0.125,0.25,0.375,0.5,0.625,0.75,0.875,1))


matlay <- rbind(matrix(1,nrow=2,ncol=3),matrix(2,nrow=3,ncol=3))
mycols <- mycols.ftn(n=2)
thisdir <- getwd()

for(i in 1:nshort){
	
	fname <- paste0("longkalman-",names(integrated.data2)[i] ,".txt")
	kalman.mat <- read.table(file=fname, header=TRUE)
	comm.mat <- Two.spp.commats[[i]]
	tmat <- matrix(rep(kalman.mat[,1],ncol(comm.mat)), ncol=ncol(comm.mat), nrow=length(kalman.mat[,1]),byrow=FALSE)	
	
	B <- mles.list[[i]]$B
	B11 <- signif(B[1,1],3);
	B12 <- signif(B[1,2],3);
	B21 <- signif(B[2,1],3);
	B22 <- signif(B[2,2],3);
	pcasc <- signif(pca.scores[i], digits=3)
	bquote.bit <- bquote(paste(.(short.labs[i]), ', PCA score = ',.(pcasc),', B'[11]*'=',.(B11),', B'[12]*'=',.(B12),', B'[21]*'=',.(B21),', B'[22]*'=',.(B22)))

	## plotfname <- paste0(short.labs[i],"-TwoSppAllplots.tiff")
	## coment out next two lines if you want the plots in log scale 

	#levels(belongingsF) <- c("Unstable","Highly unstable", "Stable","Highly stable") # check that this is so
	stab.level <- belongingsF[i]
	plotfname <- paste0(thisdir,"/Wi-2sppPlots/",stab.level,"/N-",short.labs[i],"-TwoSppAllplots.tiff")	
	comm.mat <- exp(comm.mat)
	tiff(plotfname, width=8,height=6, units="in", res=600, compression="lzw",type="cairo", family="times")
	
	layout(mat=matlay)
	par(mar=c(3,3,3,1),  oma=c(2,2,2,2), mgp=c(2,0.5,0))

	# First plot
	plot(0,0,xlim=range(xscores),ylim=c(0,1), type="n",xaxt="n", yaxt="n", bty="n",xlab="",ylab="")
	points(xscores, dens4plot$y, type="l", lwd=3, col="blue")
	polygon(x=c(0,xscores,0,rev(xscores),0), y=c(0,rep(0,lenx),0,rev(dens4plot$y),0), col="blue")
	abline(v=pca.scores[i],lwd=2,col="red", lty=1)
	axis(1,at=range.x,labels=round(range.x,digits=2), tick=TRUE, cex.axis=0.8)
	legend("topleft", legend=c(paste0("PCA score for ",short.labs[i]),"Density of all PCA scores"), col=c("red","blue"), 
	lty=1, lwd=2, bty="n",cex=1.15)

	# Second plot
	#matplot(x=tmat,y=comm.mat,type="b", pch=16, col=mycols, xlab="Time (days)", ylab="Log-abundances", cex.lab=1.25, main=main.lab,bty="n")
	# comment out to plot in log scale
	matplot(x=tmat,y=comm.mat,type="b", pch=16, col=mycols, xlab="Time (days)", ylab="Abundances", cex.lab=1.25, main=bquote.bit,bty="n")
	spp.names <- colnames(comm.mat)
	legend("topright", legend=spp.names, col=mycols, lty=1,lwd=2, bty="n", cex=1.15)
	dev.off()
}

#### Stopping times


# W-34 (stable)
 mles.list[[87]]$B
 
# W-19 (Highly Unstable)
 mles.list[[17]]$B

# W-11 (Highly unstable)
 mles.list[[10]]$B

# W-30 (highly stable) 
 mles.list[[25]]$B 

# W-34 (stable)
# >  mles.list[[87]]$B
#                  Lactobacillus Non-Lactobacillus
#Lactobacillus         0.2032576       -0.03518137
#Non-Lactobacillus    -0.2253270        0.88273072
#>  

#> # W-11 (Highly unstable)
#>  mles.list[[10]]$B
#                  Lactobacillus Non-Lactobacillus
#Lactobacillus      9.205391e-01        0.09986235
#Non-Lactobacillus -6.788782e-06        0.99654138

# W-19 (Highly Unstable)
#>  mles.list[[17]]$B
#                  Lactobacillus Non-Lactobacillus
#Lactobacillus         0.9988011      -0.001645383
#Non-Lactobacillus     0.1531567       0.754865813

# W-30 (highly stable) 
#>  mles.list[[25]]$B 
#                  Lactobacillus Non-Lactobacillus
#Lactobacillus        0.29999003        -0.2870767
#Non-Lactobacillus   -0.03207672         0.1570908


source("Poisson-MAR-ftns.R")

len.sim <- 100
#pct.drop <- 0.60

ith.ind <- 19#25#10#
# retrieving the time vector
fname <- paste0("longkalman-",names(integrated.data2)[ith.ind] ,".txt")
kalman.mat <- read.table(file=fname, header=TRUE)
tmat <- matrix(rep(kalman.mat[,1],2), ncol=2, nrow=length(kalman.mat[,1]),byrow=FALSE)	
tmat.sim <- matrix(rep(1:len.sim,2), ncol=2, nrow=len.sim, byrow=FALSE)

ith.commat <- Two.spp.commats[[ith.ind]]
ith.len <- nrow(ith.commat)
ith.xo <- ith.commat[ith.len,]
ith.A <- mles.list[[ith.ind]]$A
ith.B <- mles.list[[ith.ind]]$B
ith.S <- mles.list[[ith.ind]]$sigma
ith.muvec <- exp(ginv(diag(2)-ith.B)%*%ith.A)
sim.mat <- Twospp.mar(A=ith.A, B=ith.B, Sigma=ith.S,len=len.sim,Xo=ith.xo)


# Comparing stopping times of an unstable community vs. a stable community
# Doesn't work well...from here till line 531

ith.cuttoff <- min(sim.mat[,1])
ith.cuttoffpct <- (ith.cuttoff/sim.mat[1,1])*100
main.lab <- paste0(short.labs[ith.ind],", PCA score = ",signif(pca.scores[ith.ind],digits=2), ", min(Lactobacillus) = ", 
                   signif(ith.cuttoffpct,digits=2), "% of Initial")
plotfname <- paste0("PFUCT-",short.labs[ith.ind],".tiff")
tiff(plotfname, width=8,height=8, units="in", res=600, compression="lzw",type="cairo", family="times")
par(mar=c(3,3,3,1),  oma=c(2,2,2,2), mgp=c(2,0.5,0))
matplot(x=tmat.sim,y=sim.mat,type="b", pch=16, col=mycols, xlab="Time (days)", ylab="Abundances", cex.lab=1.25) 
abline(h=ith.cuttoff, col="red", lwd=2, lty=1)
abline(h=sim.mat[1,1], lty=2,col="black")
spp.names <- colnames(comm.mat)
legend("topright", legend=spp.names, col=mycols, lty=1,lwd=2, bty="n", cex=1.15)
dev.off()


len.sim <- 100
nchains <- 6

statio.stats1 <- list()
statio.stats2 <- rep(0,nshort)

for(i in 1:nshort){
	ith.ind <- i#10#25#
	ith.commat <- Two.spp.commats[[ith.ind]]
	ith.len <- nrow(ith.commat)
	ith.xo <- ith.commat[ith.len,]
	ith.A <- mles.list[[ith.ind]]$A
	ith.B <- mles.list[[ith.ind]]$B
	ith.S <- mles.list[[ith.ind]]$sigma
	ith.muvec <- exp(ginv(diag(2)-ith.B)%*%ith.A)
	
	Vec.V <- ginv(diag(4) - kronecker(ith.B,ith.B))%*%as.vector(ith.S) # Eq. 17th Ives et al 2003
	V <- matrix(Vec.V, nrow=2,ncol=2,byrow=FALSE)
	nzeros <- sum(V<=0)
	if(nzeros>0){rnd.start=FALSE}else{rnd.start=TRUE}
	nzeros2 <- sum(ith.muvec==0)
	if(nzeros2>0){rnd.start=FALSE}
	chains.list <- list()
	for(j in 1:nchains){
		sim.mat <- Twospp.mar(A=ith.A, B=ith.B, Sigma=ith.S,len=len.sim,Xo=ith.xo,rnd.init=rnd.start)
  		chains.list[[j]] <- mcmc(sim.mat)
	}  

	simsasmcmc <- mcmc.list(chains.list)
	diag.stats <- gelman.diag(simsasmcmc)
	statio.stats1[[i]] <- diag.stats[[1]]
	statio.stats2[i] <- diag.stats[[2]]
}	

tiff("GelmansRvsPCAscore.tiff", width=8,height=8, units="in", res=600, compression="lzw",type="cairo", family="times")
par(mar=c(3,3,3,1),  oma=c(2,2,2,1), mgp=c(2,0.5,0))
plot(pca.scores, log(statio.stats2),pch=19, col=cols4plot,ylab="Gelman's R statistic", xlab="PCA score", bty="n",yaxt="n",cex.lab=1.25)
abline(h=log(1.25),col="red", lwd=2)
legend("topleft", legend=c("Highly stable","Stable", "Unstable","Highly unstable"), col=c(kcols[4], kcols[3],kcols[1],kcols[2]), pch=rep(19,4), cex=1.25, bty="n")
axis(side=2, at=c(0, 0.5,1,1.5,2), labels=signif(exp(c(0,0.5,1,1.5,2)),digits=2) )
dev.off()


# Saving the belonging stats in the data frame:

belongings.df$nshort <- 1:nshort
belongings.df$GelmansR <- statio.stats2


#### instability does not mean pop is going to crash, it may grow exponentially like here!!!
ith.ind <- 17#25#10#
fname <- paste0("longkalman-",names(integrated.data2)[ith.ind] ,".txt")
kalman.mat <- read.table(file=fname, header=TRUE)
tmat <- matrix(rep(kalman.mat[,1],2), ncol=2, nrow=length(kalman.mat[,1]),byrow=FALSE)	
ith.commat <- exp(Two.spp.commats[[ith.ind]])

matplot(x=tmat,y=ith.commat,type="b", pch=16, col=mycols, xlab="Time (days)", ylab="Abundances", cex.lab=1.25, main=short.labs[ith.ind],bty="n")
spp.names <- colnames(ith.commat)
legend("topright", legend=spp.names, col=mycols, lty=1,lwd=2, bty="n", cex=1.15)

tiff(paste0(thisdir,"/Wi-2sppPlots/Highly unstable/N-W-19-LactobExpGrowth.tiff"), width=8,height=8, units="in", res=600, compression="lzw",type="cairo", family="times")
par(mar=c(3,3,3,1),  oma=c(2,2,2,1), mgp=c(2,0.5,0))
plot(x=tmat[,1], y=ith.commat[,1], type="b", pch=16, col=mycols[1], xlab="Time (days)", ylab="Abundances", cex.lab=1.25, main=short.labs[ith.ind],bty="n")
dev.off()

save.image("Stability4TwoSppData2.0.RData")	


# Now, computing the probability of falling under clinical threshold

thres <- 0.25
thres.times <- 2:20#c(3,5,10,20,50,75,100)
len.sim <- max(thres.times)
nthreshs <- length(thres.times)
labs.thresh <- paste0("PT",thres.times)
Pfuct.mat <- matrix(0,nrow=nshort,ncol=nthreshs)
colnames(Pfuct.mat) <- labs.thresh
row.names(Pfuct.mat) <- short.labs
nchains <- 100
NtimesNo <- 1

# Idea: track percent of times pop. dropped below threshold, not only if it did

for(i in 1:nshort){
	ith.ind <- i#10#25#
	ith.commat <- Two.spp.commats[[ith.ind]]
	ith.len <- nrow(ith.commat)
	ith.A <- mles.list[[ith.ind]]$A
	ith.B <- mles.list[[ith.ind]]$B
	ith.S <- mles.list[[ith.ind]]$sigma
	ith.muvec <- ginv(diag(2)-ith.B)%*%ith.A
	ith.xo <- log(NtimesNo) +  ith.commat[ith.len,]
	#ith.initot <- sum(exp(ith.xo))
	#Nthres <- thres*ith.initot
	for(k in 1:nthreshs){
		code.line1 <- paste0("cross.thres",k," <-rep(0,nchains)")
		eval(parse(text=code.line1))
	}

	for(j in 1:nchains){
		sim.mat <- Twospp.mar(A=ith.A, B=ith.B, Sigma=ith.S,len=len.sim,Xo=ith.xo,rnd.init=FALSE)
		for(h in 1:nthreshs){
#			test <- sum((sim.mat[1:thres.times[h],1]/apply(sim.mat[1:thres.times[h],],1,sum))<thres)>0			
			test <- sum((sim.mat[thres.times[h],1]/sum(sim.mat[thres.times[h],]))<thres)>0
			code.line2 <- paste0("cross.thres",h,"[j] <- test")
			eval(parse(text=code.line2))
		}
	}
	for(ell in 1:nthreshs){
		code.line3 <- paste0("Pfuct.mat[i,",ell,"] <- sum(cross.thres",ell,")/nchains")
		eval(parse(text=code.line3))
	}  
				
}	

Tthres.mat <- matrix(rep(thres.times,nshort), nrow=nthreshs, ncol=nshort,byrow=FALSE)
tPfuct.mat <- t(Pfuct.mat)

for(i in 1:nshort){
	stab.level <- belongingsF[i]
	plotfname <- paste0(thisdir,"/Wi-2sppPlots/",stab.level,"/PFUCT",thres,"-", short.labs[i],".tiff")	
	tiff(plotfname, width=6,height=10, units="in", res=600, compression="lzw",type="cairo", family="times")
	par(mar=c(3,3,3,1),  oma=c(2,2,2,1), mgp=c(2,0.5,0))
	matplot(x=Tthres.mat,y=1-tPfuct.mat, type="b", pch=1, col=mycols.ftn(nshort), lwd=2, xlab="Days (t)", ylab=paste0("P(Lactobacillus is above threshold)"), 
	cex.lab=1.25, main=paste0(short.labs[i],": Pr(persisting above ",thres*100,"% of total on day 't')"),bty="n")
	points(Tthres.mat[,i], 1-tPfuct.mat[,i], type="b", pch=16, col="black", lwd=2)
	dev.off()
}

