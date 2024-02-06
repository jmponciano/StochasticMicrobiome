# Lactobacillus vs. "Gardnerella_vaginalis" vs. the rest 

source("CLS3.0.R")

load("integrateddata2.RData")
nshort <- length(integrated.data2)
mles.list <- list()
stab.mles <- list()
Three.spp.commats <- list()

# i=26, woman31 does not have gardnerella, use the only other two species: Prevotella_bivia and gDialister
# i=34, woman44 does not have gardnerella, use Prevotella_bivia (possibly +g_Dialister), g_Mobiluncus
# i=
for(i in 1:nshort){
	
	print(paste("community number ",i, "case name: ",names(integrated.data2)[i] ))
	
	fname <- paste0("longkalman-",names(integrated.data2)[i] ,".txt")
	kalman.mat <- read.table(file=fname, header=TRUE)
	comm.mat1 <- kalman.mat[,-1]
	print(colnames(comm.mat1))
	print(grepl("Gardnerella_vaginalis", colnames(comm.mat1)))
	
	Lactocols <- which(grepl("Lactobacillus", colnames(comm.mat1))==TRUE,arr.ind=TRUE)
	if(i==26){
		Prevocol <- which(grepl("Prevotella_bivia", colnames(comm.mat1))==TRUE,arr.ind=TRUE)
		Dialistcol <-  which(grepl("g_Dialister", colnames(comm.mat1))==TRUE,arr.ind=TRUE)
		Lacto.vec <- log(apply(comm.mat1[,Lactocols],1,sum))
		Prevovec <- log(comm.mat1[,Prevocol])
		Dialistvec <- log(comm.mat1[,Dialistcol])
		comm.mat <- cbind(Lacto.vec, Prevovec,Dialistvec)
		colnames(comm.mat) <- c("Lactobacillus", "Prevotella_bivia", "g_Dialister")
		}else if(i==34){
			Prevocol <- which(grepl("Prevotella_bivia", colnames(comm.mat1))==TRUE,arr.ind=TRUE)
			Prevovec <- log(comm.mat1[,Prevocol])			
			Mobilcol <- which(grepl("g_Mobiluncus", colnames(comm.mat1))==TRUE,arr.ind=TRUE)
			Mobilvec <- log(comm.mat1[,Mobilcol])			
			Lacto.vec <- log(apply(comm.mat1[,Lactocols],1,sum))
			comm.mat <- cbind(Lacto.vec, Prevovec,Mobilvec)
			colnames(comm.mat) <- c("Lactobacillus", "Prevotella_bivia", "g_Mobiluncus")
		}else{
		Gardnercol<- which(grepl("Gardnerella_vaginalis", colnames(comm.mat1))==TRUE, arr.ind=TRUE)
		Gardnervec <- log(comm.mat1[,Gardnercol])
		Othercols <- (1:ncol(comm.mat1))[-c(Lactocols,Gardnercol)]
		if(length(Lactocols)==1){Lacto.vec <- log(comm.mat1[,Lactocols])}else{
			Lacto.vec <- log(apply(comm.mat1[,Lactocols],1,sum))}
		if(length(Othercols)==1){
			Other.vec <- log(comm.mat1[,Othercols])
			other.name <- colnames(comm.mat1)[Othercols]
			comm.mat <- cbind(Lacto.vec,Gardnervec,Other.vec)
			colnames(comm.mat) <- c("Lactobacillus","Gardnerella_vaginalis", other.name)
			}else{
				Other.vec <- log(apply(comm.mat1[,Othercols],1,sum))
				comm.mat <- cbind(Lacto.vec,Gardnervec, Other.vec)
				colnames(comm.mat) <- c("Lactobacillus", "Gardnerella_vaginalis", "All other spp.")	
			}
		}	
			
	Three.spp.commats[[i]] <- comm.mat
	comm.mat.mles <- mars.cls(comm.mat)
	mles.list[[i]] <- comm.mat.mles
	stab.mles[[i]] <- stability(comm.mat.mles$B, comm.mat.mles$sigma)	
	
}

names(mles.list) <- names(integrated.data2)
names(stab.mles) <- names(integrated.data2)
names(Three.spp.commats) <- names(integrated.data2)

# Print results just to check
for(i in 1:nshort){
	
	print(mles.list[[i]]$B)
	
}


# Now, store the effects of every spp on every other
# Except for communities 26 and 34 that don't have Gardnerella
all.interactions <- matrix(0,nrow=nshort,ncol=9)
row.names(all.interactions) <- names(mles.list)
colnames(all.interactions) <- c("L2L", "L2G", "L2R", "G2L", "G2G", "G2R", "R2L", "R2G","R2R")
for(i in 1:nshort){
	all.interactions[i,] <- as.vector(mles.list[[i]]$B)
}


main.labs <- c("Lactobacillus on Lactobacillus", "Lactobacillus on Gardnerella", "Lactobacillus on the rest", 
"Gardnerella on Lactobacillus", "Gardnerella on Gardnerella", "Gardnerella on the rest", "All others on Lactobacillus", "All others on Gardnerella","All others on themselves")

par(mfrow=c(2,3),mar=c(3,4,3,1),  oma=c(3,4,3,1), mgp=c(2,0.5,0))
hist(all.interactions[,1], main = "Lactobacillus on Lactobacillus", xlab="Strength of Density Dependence", cex.lab=1.5, ylab="Frequency", cex.main=1.5)
hist(all.interactions[,5], main = "Gardnerella on Gardnerella", xlab="Strength of Density Dependence", cex.lab=1.5, ylab="Frequency", cex.main=1.5)
hist(all.interactions[,9], main = "All others on themselves", xlab="Strength of Density Dependence", cex.lab=1.5, ylab="Frequency", cex.main=1.5)
plot(all.interactions[,2], all.interactions[,4], xlab="Lactobacillus on Gardnerella", ylab="Gardnerella on Lactobacillus", pch=19,xlim=c(-1.5,1.5), ylim=c(-1.5,1.5), cex.lab=1.5)
plot(all.interactions[,3], all.interactions[,7], xlab="Lactobacillus on the rest", ylab="All others on Lactobacillus", pch=19,cex.lab=1.5)
plot(all.interactions[,6], all.interactions[,8], xlab="Gardnerella on the rest", ylab="All others on Gardnerella", pch=19,cex.lab=1.5)

save.image(file="ThreeSppWomensInteractions-November2023.RData")


