# Meeting with Larry: March 5 2019:
#
# 1. In the manuscript NIH manuscript 151018.docx:
#
# Figure 1: is now figure 2, figure 1 should show the temporal variation of abundances AND why 
# a categorical classsification like NUGENT falls short.

# Focus of the figure is to say:  look , under the same amount of environmental temporal
# variability, the variability in abundances over time changes dramatically as a function of the
# ecological community structure.

# Now think about which set of interactions to illustrate (passenger-driver idea, 
# does including 4 species changes qualitatively the patters)

#### 
source("CleaningDataLists.R")
all.hmpdata <- read.csv("HMP_metadata_updated_taxonomy_StR_CSTs_020619.csv", header=TRUE)

names(all.hmpdata)
all.names <- names(all.hmpdata)
women.num <- unique(all.hmpdata$PID)
id.week.day <- all.hmpdata$UID
patient.id <- all.hmpdata$PID
allweeks <- all.hmpdata$WEEK
alldays <- all.hmpdata$DAY
nugent.score <- all.hmpdata$NUGENT_SCORE
nugent.class <- all.hmpdata$NUGENT_CLASS
ph <- all.hmpdata$PH
menses.norm <- all.hmpdata$MENSTRUATION_NORMALIZED
menses.norm.int <- all.hmpdata$MENSTRUATION_NORM_INT 
total.reads <- all.hmpdata$X16s_total_reads
last.bact.sp <- "g_Pediococcus"
first.bact.sp <- "Lactobacillus_iners"
all.bact.reads <- cbind(all.hmpdata[,which(all.names==first.bact.sp,arr.ind=TRUE):which(all.names==last.bact.sp,arr.ind=TRUE)],total.reads)
nbact.spp <- ncol(all.bact.reads)
colnames(all.bact.reads)[nbact.spp] <- "Total_sum"
total.reads.mat <- matrix(rep(total.reads,nbact.spp),ncol=nbact.spp, nrow=length(total.reads),byrow=FALSE)
all.bact.pcts <- all.bact.reads/total.reads.mat
qPCR.Ct.rep1 <- all.hmpdata$qPCR_Ct_rep1
qPCR.Ct.rep2 <- all.hmpdata$qPCR_Ct_rep2
qPCR.Ct.rep3 <- all.hmpdata$qPCR_Ct_rep3
qPCR.copies.rep1 <- all.hmpdata$qPCR_copies_swab_rep1
qPCR.copies.rep2 <- all.hmpdata$qPCR_copies_swab_rep2
qPCR.copies.rep3 <- all.hmpdata$qPCR_copies_swab_rep3
qPCR.copies.ave  <- all.hmpdata$qPCR_copies_swab

# just a little visualization of the Ct number vs the number of copies
#plot(qPCR.Ct.rep1[qPCR.Ct.rep1>10], log(qPCR.copies.rep1[qPCR.Ct.rep1>10]), pch=16)


total.PCR.Ct1.mat <- matrix(rep(qPCR.Ct.rep1,nbact.spp),ncol=nbact.spp, nrow=length(qPCR.Ct.rep1),byrow=FALSE)
total.PCR.Ct2.mat <- matrix(rep(qPCR.Ct.rep2,nbact.spp),ncol=nbact.spp, nrow=length(qPCR.Ct.rep2),byrow=FALSE)
total.PCR.Ct3.mat <- matrix(rep(qPCR.Ct.rep3,nbact.spp),ncol=nbact.spp, nrow=length(qPCR.Ct.rep3),byrow=FALSE)

total.PCR.copies1.mat <- matrix(rep(qPCR.copies.rep1,nbact.spp),ncol=nbact.spp, nrow=length(qPCR.copies.rep1),byrow=FALSE)
total.PCR.copies2.mat <- matrix(rep(qPCR.copies.rep2,nbact.spp),ncol=nbact.spp, nrow=length(qPCR.copies.rep2),byrow=FALSE)
total.PCR.copies3.mat <- matrix(rep(qPCR.copies.rep3,nbact.spp),ncol=nbact.spp, nrow=length(qPCR.copies.rep3),byrow=FALSE)

PCR.Ct1.mat <- all.bact.pcts*total.PCR.Ct1.mat
PCR.Ct2.mat <- all.bact.pcts*total.PCR.Ct2.mat
PCR.Ct3.mat <- all.bact.pcts*total.PCR.Ct3.mat

PCR.copies1.mat <- all.bact.pcts*total.PCR.copies1.mat/4000 # 4000 adjustment per Michael France
PCR.copies2.mat <- all.bact.pcts*total.PCR.copies2.mat/4000
PCR.copies3.mat <- all.bact.pcts*total.PCR.copies3.mat/4000




hmpdata <- data.frame(id.week.day,patient.id,allweeks,alldays,nugent.score,nugent.class,ph,menses.norm, menses.norm.int, total.reads, all.bact.reads, qPCR.Ct.rep1, qPCR.Ct.rep2, qPCR.Ct.rep3, qPCR.copies.rep1, qPCR.copies.rep2, qPCR.copies.rep3) 

PCR.Ct.data1 <- data.frame(id.week.day,patient.id,allweeks,alldays,nugent.score,nugent.class,ph,menses.norm, menses.norm.int,PCR.Ct1.mat)
PCR.Ct.data2 <- data.frame(id.week.day,patient.id,allweeks,alldays,nugent.score,nugent.class,ph,menses.norm, menses.norm.int,PCR.Ct2.mat)
PCR.Ct.data3 <- data.frame(id.week.day,patient.id,allweeks,alldays,nugent.score,nugent.class,ph,menses.norm, menses.norm.int,PCR.Ct3.mat)

PCR.copies.data1 <- data.frame(id.week.day,patient.id,allweeks,alldays,nugent.score,nugent.class,ph,menses.norm, menses.norm.int,PCR.copies1.mat)
PCR.copies.data2 <- data.frame(id.week.day,patient.id,allweeks,alldays,nugent.score,nugent.class,ph,menses.norm, menses.norm.int,PCR.copies2.mat)
PCR.copies.data3 <- data.frame(id.week.day,patient.id,allweeks,alldays,nugent.score,nugent.class,ph,menses.norm, menses.norm.int,PCR.copies3.mat)



###  Creating a list of lists with all the data separated per woman:

nwomen <- length(unique(patient.id)) # or length(women.num)

all.woments <- list()

for(i in 1:nwomen){
	
	ith.patientid <- women.num[i]
	ith.rows <- which(hmpdata$patient.id ==ith.patientid,arr.ind=TRUE)
	ith.woman.metadata <- cbind(hmpdata[ith.rows,1:9],all.hmpdata[ith.rows,c(3,150)])
	ith.woman.PCR.data1 <- PCR.copies1.mat[ith.rows,]
	ith.woman.PCR.data2 <- PCR.copies2.mat[ith.rows,]
	ith.woman.PCR.data3 <- PCR.copies3.mat[ith.rows,]	
	
	all.woments[[i]] <- list(ith.woman.metadata,ith.woman.PCR.data1, ith.woman.PCR.data2, ith.woman.PCR.data3)
}

names(all.woments) <- paste0("woman",women.num)

#### Choosing the set of species to look at:  we ran these 14 spp
focal.spp <- c("Lactobacillus_iners", "Lactobacillus_crispatus", "Lactobacillus_jensenii", "Lactobacillus_gasseri", "Gardnerella_vaginalis", "g_Megasphaera", "g_Gemella", "g_Leptotrichia", "g_Parvimonas", "Atopobium_vaginae", "Prevotella_buccalis", "Prevotella_bivia","g_Dialister","g_Mobiluncus", "Total_sum")  

#####  Checking if all of these spp are in the new data set
for(i in 1:(length(focal.spp)-1)){ # the -1 is there because we don't want to check if the total is 
								   # the original data set, it doesn't make sense...
	
	print(paste(which(all.names==focal.spp[i], arr.ind=TRUE), focal.spp[i]))	
	
}

#Printing all the "wanted genus" spp:
all.bact.spp <- colnames(all.bact.reads)
nbact <- length(all.bact.spp)
wanted.bact <-  "Mobiluncus"#"Total_sum"#"Dialister"#"Parvimonas"#"Parvimonas" #"Eggerthella" # "Leptotrichia"#"Mycoplasma" #"Gemella"#"Megasphaera" #"Porphyromonas" # "Dialister" #"Leptotrichia" #"Lactobacillus"
is.wanted.in <- grepl(wanted.bact,all.bact.spp, fixed=TRUE)
all.bact.spp[is.wanted.in>0]
 
#### Creating a new list with one list per woman.  Each one of the woman's list contains
#### three matrices (one of the three reps) of the selected time series

all.woments.short <- list()

cols.focalspp <- match(x=focal.spp,table=names(all.woments[[1]][[2]]))

for(i in 1:nwomen){
	################	START HERE, INSERT PH AFTER TIME AND DRAG IT ALL OVER !!!!!!!!!!!!!!!
		
	datarep1 <- cbind(all.woments[[i]][[1]][,10], all.woments[[i]][[2]][,cols.focalspp],all.woments[[i]][[1]][,7])
	datarep2 <- cbind(all.woments[[i]][[1]][,10],all.woments[[i]][[3]][,cols.focalspp],all.woments[[i]][[1]][,7])
	datarep3 <- cbind(all.woments[[i]][[1]][,10],all.woments[[i]][[4]][,cols.focalspp],all.woments[[i]][[1]][,7])		

	all.woments.short[[i]] <- list(datarep1,datarep2,datarep3)
}

names(all.woments.short) <- paste0("woman",women.num)

##  OK: cleaning up the raw data (0's and NA's):
clean.data <- nas.cleaning(womentsdat=all.woments.short)

#save.image(file="allhmpdata.RData")
#load(file="allhmpdata.RData")



### Saving an indicator variable testing whether there is enough data for estimation
### for each woman, and for each replicate per woman

dims.womenmats <- list()
estim.test <- rep(0,nwomen)
nestim.women <- 0


estim.woments <- list()
#short.miss.inds <- list()
short.dim.womats <- list()
estim.womencovars <- list()

for(i in 1:nwomen){
	
	np1 <- dim(clean.data$out.ts[[i]][[1]]);
	np2 <- dim(clean.data$out.ts[[i]][[2]]);
	np3 <- dim(clean.data$out.ts[[i]][[3]]);

	dims.womenmats[[i]] <- rbind(np1,np2,np3);

	n1 <- np1[1];n2 <- np2[1];n3 <- np3[1];
	p <- np1[2]-1; # minus 1 because there is one column with the totals for every spp

	test.val <- ((n1-1)*p + (n2-1)*p +(n3-1)*p)>p*(2*p+1)
	estim.test[i] <- test.val 
	if(test.val==TRUE){
		nestim.women <- nestim.women+1;
		estim.woments[[nestim.women]] <- clean.data$out.ts[[i]]
		cov.list <- list()
		for(k in 1:3){

			missk <- clean.data$out.narows[[i]][[k]];
			sum.miss <- sum(missk)
			if(sum.miss>0){
				cov.list[[k]] <- all.woments[[i]][[1]][-missk,]			
			}else if(sum.miss==0){
				cov.list[[k]] <- all.woments[[i]][[1]]
			}
			
		}

		estim.womencovars[[nestim.women]] <- cov.list
		#short.miss.inds[[nestim.women]] <- all.miss.inds[[i]]
		short.dim.womats[[nestim.women]] <- rbind(np1,np2,np3);
		print(c(nestim.women,i))
	}

}

names(dims.womenmats) <- names(all.woments.short)
names(estim.woments) <- names(all.woments.short)[which(estim.test==1,arr.ind=TRUE)]
#names(short.miss.inds) <- names(estim.woments)
names(short.dim.womats) <- names(estim.woments)
names(estim.womencovars) <- names(estim.woments) 


print(nestim.women)

save.image(file="allhmpdata.RData")

save(estim.woments, file="estim.woments.RData")
save(short.miss.inds, file="short.miss.inds.RData")
save(short.dim.womats, file="short.dim.womats.RData")
save(estim.womencovars, file="estim.womencovars.RData")


### Run from here onwards
load(file="estim.woments.RData")
#load(file="short.miss.inds.RData")
load(file="short.dim.womats.RData")
load(file="estim.womencovars.RData")


nshort <- length(names(estim.womencovars))
short.labs <- rep(0,nshort)
mediangt4.5 <- rep(0,nshort)



###############################  pH's figures ##########################
tiff("VaginalpH-pinkNgrey.tiff", width=6,height=8, units="in", res=600, compression="lzw",
	type="cairo", family="times")

par(mfrow=c(2,1), oma=c(3,1.5,1.5,1.5),mar=c(3.5,3.5,1.5,1), mgp=c(2,0.75,0))
plot(0,0, type="n", xlim=c(1,44), ylim=c(3.8,7.2), xaxt="n", ylab="Vaginal pH", xlab="Subjects 1 to 44   ", bty="l", cex.lab=1.5)
polygon(x=c(0,46,46,0), y=c(4.5,4.5,7.2,7.2), col="pink", border=FALSE)

for(i in 1:44){
	
	ith.phvec <- estim.womencovars[[i]][[1]][,7]
	boxplot(ith.phvec, at=i, add=TRUE, axes=FALSE, outline=FALSE, boxwex=0.9,col="grey")
	short.labs[i] <- estim.womencovars[[i]][[1]][1,2]
	mediangt4.5[i] <- median(ith.phvec, na.rm=TRUE) > 4.5
	
}
axis(1,1:44,labels=FALSE, tick=FALSE, cex.lab=1.5)

plot(0,0, type="n", xlim=c(1,45), ylim=c(3.8,7.2), xaxt="n", ylab="Vaginal pH", xlab="Subjects 45 to 88   ", bty="l", cex.lab=1.5)
polygon(x=c(0,46,46,0), y=c(4.5,4.5,7.2,7.2), col="pink", border=FALSE)
for(i in 45:nshort){
	
	j <- i-44
	ith.phvec <- estim.womencovars[[i]][[1]][,7]
	boxplot(ith.phvec, at=j, add=TRUE, axes=FALSE, outline=FALSE, boxwex=0.9,col="grey")
	short.labs[i] <- estim.womencovars[[i]][[1]][1,2]
	mediangt4.5[i] <- median(ith.phvec, na.rm=TRUE) > 4.5
		
}
axis(1,1:44,labels=FALSE, tick=FALSE, cex.lab=1.5)
#mtext("Boxplots of the vaginal pH for 130 subjects during 70 days", side=1, outer=TRUE, cex=1.5)

dev.off()


############### Long Horizontal single pH figure. ################
tiff("VaginalpH-pinkNgrey-long.tiff", width=8,height=6, units="in", res=600, compression="lzw",
     type="cairo", family="times")
par(oma=c(1.5,1.5,1.5,1.5),mar=c(4.5,4.5,1.5,1))
plot(0,0, type="n", xlim=c(1,nshort), ylim=c(3.8,7.2), xaxt="n", ylab="Vaginal pH", 
     xlab="Women with complete time series data ", bty="l", cex.lab=1.5)
polygon(x=c(0,nshort+1,nshort+1,0), y=c(4.5,4.5,7.2,7.2), col="pink", border=FALSE)

for(i in 1:nshort){
  
  ith.phvec <- estim.womencovars[[i]][[1]][,7]
  boxplot(ith.phvec, at=i, add=TRUE, axes=FALSE, outline=FALSE, boxwex=0.9,col="grey")
  short.labs[i] <- estim.womencovars[[i]][[1]][1,2]
  mediangt4.5[i] <- median(ith.phvec, na.rm=TRUE) > 4.5
  
}
axis(1,1:nshort,labels=1:nshort, tick=TRUE, cex.lab=1.5,cex.axis=0.85,tcl=-0.35)
dev.off()




####  Printing and saving "cbind" Lactobacillus_iners with pH for all women

estim.phNLactob1 <- list()
for(i in 1:nshort){
	last.col <- ncol(estim.woments[[i]][[1]])
	uni.mat <- cbind(estim.womencovars[[i]][[1]][,c(10,5:9,11)],estim.woments[[i]][[1]])			
	colnames(uni.mat)[c(1:8,(last.col+7))] <- c("raw.days","nugent.score","nugent.class","pH","menses.norm","menses.int", "CST","math.days", "log_pH") 
	estim.phNLactob1[[i]]<- uni.mat
	print(uni.mat)
}
names(estim.phNLactob1) <- names(estim.woments)

estim.phNLactob2 <- list()
for(i in 1:nshort){

	last.col <- ncol(estim.woments[[i]][[2]])
	uni.mat <- cbind(estim.womencovars[[i]][[2]][,c(10,5:9,11)],estim.woments[[i]][[2]])			
	colnames(uni.mat)[c(1:8,(last.col+7))] <- c("raw.days","nugent.score","nugent.class","pH","menses.norm","menses.int", "CST","math.days", "log_pH") 
	estim.phNLactob2[[i]]<- uni.mat
	#print(uni.mat)
}
names(estim.phNLactob2) <- names(estim.woments)


estim.phNLactob3 <- list()
for(i in 1:nshort){
	
	last.col <- ncol(estim.woments[[i]][[3]])
	uni.mat <- cbind(estim.womencovars[[i]][[3]][,c(10,5:9,11)],estim.woments[[i]][[3]])			
	colnames(uni.mat)[c(1:8,(last.col+7))] <- c("raw.days","nugent.score","nugent.class","pH","menses.norm","menses.int", "CST","math.days", "log_pH") 
	estim.phNLactob3[[i]]<- uni.mat
	#print(uni.mat)
}
names(estim.phNLactob3) <- names(estim.woments)

#### PICK UP HERE: REGULARIZE DIMENSIONS (NUMBER OF ROWS) ACROSS THE THREE REPS
#### AND GET AVERAGED COUNTS OVER THE THREE REPS.

integrated.data <- list()
for(i in 1:nshort){
	
	nspp <- short.dim.womats[[i]][1,2] -1

	bact.lcounts1 <- estim.phNLactob1[[i]];bact.lcounts2 <- estim.phNLactob2[[i]];bact.lcounts3 <- estim.phNLactob3[[i]];
	counts1 <- exp(bact.lcounts1[,9:(8+nspp)]);counts2 <- exp(bact.lcounts2[,9:(8+nspp)]);counts3 <- exp(bact.lcounts3[,9:(8+nspp)]);	
	    

	list.counts <- list()
	list.counts[[1]] <- counts1;list.counts[[2]] <- counts2;list.counts[[3]] <- counts3;

	covars.list <- list()
	covars.list[[1]] <- bact.lcounts1[,1:8];covars.list[[2]] <- bact.lcounts2[,1:8];covars.list[[3]] <- bact.lcounts3[,1:8]	
	
	imperfect.days1 <- bact.lcounts1[,8]+1;
	imperfect.days2 <- bact.lcounts2[,8]+1;
	imperfect.days3 <- bact.lcounts3[,8]+1;	
	
	imperfect.list <- list();
	imperfect.list[[1]] <- imperfect.days1; 
	imperfect.list[[2]] <- imperfect.days2; 
	imperfect.list[[3]] <- imperfect.days3;
	max.days <- max(c(imperfect.days1,imperfect.days2,imperfect.days3))
	min.days <- min(c(imperfect.days1,imperfect.days2,imperfect.days3))
	perfect.days <- min.days:max.days
	counter <- 0

	super.avecounts <- matrix(NA, nrow=length(perfect.days),ncol=nspp)
	colnames(super.avecounts) <- colnames(counts1)
	super.covsmat <- matrix(NA,nrow=length(perfect.days),ncol=8)
	colnames(super.covsmat) <- colnames(bact.lcounts1)[1:8]
	for(j in perfect.days){
		
		counter <- counter+1;
		match1 <- match(x=j,table=imperfect.days1)
		match2 <- match(x=j,table=imperfect.days2)		
		match3 <- match(x=j,table=imperfect.days3)		
		where.matches <- c(match1,match2,match3)
		not.nas <- 1-c(is.na(match1), is.na(match2), is.na(match3))
		num.present <- sum(not.nas)
		#print(c(where.matches, num.present))
		pos.present <- which(not.nas==1,arr.ind=TRUE)
		
		if(num.present>0){
		
			new.where.matches <- where.matches[pos.present]
			rows2average <- matrix(0,nrow=num.present,ncol=nspp)
			for(k in 1:num.present){
				
				#kth.imp.days <- imperfect.list[[pos.present[k]]] 
				kth.row.num <- new.where.matches[k]
				rows2average[k,] <- matrix(unlist(list.counts[[pos.present[k]]][kth.row.num,]),nrow=1,ncol=nspp)
				
			}
			ave.entry  <- apply(rows2average,2,mean)
			super.avecounts[counter,] <- ave.entry
			
			#covs.kth.days <- imperfect.list[[pos.present[1]]]
			covs.kth.rownum <- new.where.matches[1]
			covs.row <- unlist(covars.list[[pos.present[1]]][covs.kth.rownum,])
			super.covsmat[counter,] <- covs.row		
		
		}
		
		
	} # end 1:max.days loop
	super.covsmat[,1] <- perfect.days 
	
	integrated.ithdf <- data.frame(super.covsmat,super.avecounts)
	integrated.ithdf$nugent.class <- as.factor(integrated.ithdf$nugent.class)
	levels(integrated.ithdf$nugent.class) <- levels(as.factor(bact.lcounts1[,3]))

	integrated.ithdf$CST <- as.factor(integrated.ithdf$CST)
	levels(integrated.ithdf$CST) <- levels(as.factor(bact.lcounts1[,7]))

	integrated.data[[i]] <- integrated.ithdf

}
names(integrated.data) <- names(estim.woments)


save(integrated.data, file="integrated-data.RData")

## Now removing all the rows with NA's with every integrated data file:

integrated.data2 <- list()
for(i in 1:nshort){
	
	
	integ.data <- integrated.data[[i]]
	
	rows2rm <- which(is.na(integ.data$nugent.score)==TRUE,arr.ind=TRUE)
	
	if(length(rows2rm)==0){integrated.data2[[i]] <- integ.data}else{
		
		integrated.data2[[i]] <- integ.data[-rows2rm,]
		}
	
}
names(integrated.data2) <- names(estim.woments)
length(integrated.data2)

save(integrated.data2, file="integrated-data2.RData")
save.image(file="allhmpdata3.0.RData")
# Now upload "allhmpdata3.0.RData and play with it!


