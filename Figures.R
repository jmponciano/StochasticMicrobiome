# The following code plots the figure 1 and figure 2 of the NIH data, using simulations
# The data follows the code "Data_for_Figures.R" and loads the data from the resulting code

library(MASS)
library(shape)
source('CLS.R')
source("Poisson-MAR-ftns.R")
source("error.barfunc.R")
load("FourStabilityScenarios.RData")

tiff("FourStabilityScenarios-2023.tiff",height=25,width=35,units="cm",res=300,compression="lzw",type="cairo")
mat	<-	matrix(c(1:8),nrow=2,ncol=4,byrow=T)
layout(mat)
par(mar=c(0,0,0,0),oma=c(0.5,1,0.5,0.5))
plot(0,0,pch=".",col="white",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1),axes=F)
###### Hypothetical examples 
##### Scenario 1: All interactions weak

text(0.5,0.775,"1",cex=4, col="black")
text(0.2,0.225,"2",cex=4, col="blue")
text(0.8,0.225,"3",cex=4, col="red")
arrows(0.2,0.3,0.45,0.7,lwd=2,length=0.15)
arrows(0.4,0.7,0.15,0.30,lwd=2,length=0.15)
arrows(0.3,0.25,0.69,0.25,lwd=2,length=0.15)
arrows(0.69,0.2,0.3,0.2,lwd=2,length=0.15)
arrows(0.8,0.3,0.55,0.7,lwd=2,length=0.15)
arrows(0.6,0.7,0.85,0.3,lwd=2,length=0.15)
plotcircle(mid=c(0.5,0.85),r=0.1,from=315*pi/180,to=225*pi/180,arrow=F, lcol="black")
Arrowhead(x0=0.43,y0=0.8,angle=317,arr.lwd=3,lcol="black")
plotcircle(mid=c(0.12,0.15),r=0.1,from=90*pi/180,to=360*pi/180,arrow=F,lcol="blue")
Arrowhead(x0=0.22,y0=0.15,angle=90,arr.lwd=3,lcol="blue")
plotcircle(mid=c(0.88,0.15),r=0.1,from=180*pi/180,to=90*pi/180,arrow=F,lcol="red")
Arrowhead(x0=0.78,y0=0.15,angle=90,arr.lwd=3,lcol="red")
coord.min	<-	c(0.5,0.85,0.12,0.15,0.88,0.15,0.4,0.5,0.2,0.5,0.5,0.28,0.5,0.15)
mat.min		<-	matrix(coord.min,nrow=7,ncol=2,byrow=T)
	for (i in 1:7){
		if(i==2){my.col<-"blue"}else if(i==3){my.col<-"red"}else{my.col<-"black"}
	  text(mat.min[i,1],mat.min[i,2],"-",cex=4, col=my.col)
	}
coord.plus	<-	c(0.6,0.5,0.8,0.5)	
mat.plus	<-	matrix(coord.plus,nrow=2,ncol=2,byrow=T)
	for (i in 1:2){
		text(mat.plus[i,1],mat.plus[i,2],"+",cex=4)
	}
#box()
mtext("(a)",adj=0,cex=2,line=-5)

###### Scenario 2: Only C has strong intraspecific competition
plot(0,0,pch=".",col="white",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1),axes=F)
text(0.5,0.775,"1",cex=4, col="black")
text(0.2,0.225,"2",cex=4, col="blue")
text(0.8,0.225,"3",cex=4, col="red")
arrows(0.2,0.3,0.45,0.7,lwd=2,length=0.15)
arrows(0.4,0.7,0.15,0.30,lwd=2,length=0.15)
arrows(0.3,0.25,0.69,0.25,lwd=2,length=0.15)
arrows(0.69,0.2,0.3,0.2,lwd=2,length=0.15)
arrows(0.8,0.3,0.55,0.7,lwd=2,length=0.15)
arrows(0.6,0.7,0.85,0.3,lwd=2,length=0.15)
plotcircle(mid=c(0.5,0.85),r=0.1,from=315*pi/180,to=225*pi/180,arrow=F)
Arrowhead(x0=0.43,y0=0.8,angle=317,arr.lwd=3)
plotcircle(mid=c(0.12,0.15),r=0.1,from=90*pi/180,to=360*pi/180,arrow=F,lcol="blue")
Arrowhead(x0=0.22,y0=0.15,angle=90,arr.lwd=3,lcol="blue")
plotcircle(mid=c(0.88,0.15),r=0.1,from=180*pi/180,to=90*pi/180,arrow=F,lwd=6,lcol="red")
Arrowhead(x0=0.78,y0=0.15,angle=90,arr.lwd=3,lcol="red")
coord.min	<-	c(0.5,0.85,0.12,0.15,0.88,0.15,0.4,0.5,0.2,0.5,0.5,0.28,0.5,0.15)
mat.min		<-	matrix(coord.min,nrow=7,ncol=2,byrow=T)
	for (i in 1:7){
	  if(i==2){my.col<-"blue"}else if(i==3){my.col<-"red"}else{my.col<-"black"}
	  text(mat.min[i,1],mat.min[i,2],"-",cex=4, col=my.col)
	}
coord.plus	<-	c(0.6,0.5,0.8,0.5)	
mat.plus	<-	matrix(coord.plus,nrow=2,ncol=2,byrow=T)
	for (i in 1:2){
		text(mat.plus[i,1],mat.plus[i,2],"+",cex=4)
	}
#box()
mtext("(b)",adj=0,cex=2,line=-5)


####### Scenario 3: Weak intraspecific competittion but B has strong negative effect on A and C
plot(0,0,pch=".",col="white",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1),axes=F)
text(0.5,0.775,"1",cex=4, col="black")
text(0.2,0.225,"2",cex=4, col="blue")
text(0.8,0.225,"3",cex=4, col="red")
arrows(0.2,0.3,0.45,0.7,lwd=6,length=0.15)
arrows(0.4,0.7,0.15,0.30,lwd=2,length=0.15)
arrows(0.3,0.25,0.69,0.25,lwd=6,length=0.15)
arrows(0.69,0.2,0.3,0.2,lwd=2,length=0.15)
arrows(0.8,0.3,0.55,0.7,lwd=2,length=0.15)
arrows(0.6,0.7,0.85,0.3,lwd=2,length=0.15)
plotcircle(mid=c(0.5,0.85),r=0.1,from=315*pi/180,to=225*pi/180,arrow=F)
Arrowhead(x0=0.43,y0=0.8,angle=317,arr.lwd=3)
plotcircle(mid=c(0.12,0.15),r=0.1,from=90*pi/180,to=360*pi/180,arrow=F,lcol="blue")
Arrowhead(x0=0.22,y0=0.15,angle=90,arr.lwd=3,lcol="blue")
plotcircle(mid=c(0.88,0.15),r=0.1,from=180*pi/180,to=90*pi/180,arrow=F,lcol="red")
Arrowhead(x0=0.78,y0=0.15,angle=90,arr.lwd=3,lcol="red")
coord.min	<-	c(0.5,0.85,0.12,0.15,0.88,0.15,0.4,0.5,0.2,0.5,0.5,0.28,0.5,0.15)
mat.min		<-	matrix(coord.min,nrow=7,ncol=2,byrow=T)
	for (i in 1:7){
	  if(i==2){my.col<-"blue"}else if(i==3){my.col<-"red"}else{my.col<-"black"}
	  text(mat.min[i,1],mat.min[i,2],"-",cex=4, col=my.col)
	}
coord.plus	<-	c(0.6,0.5,0.8,0.5)	
mat.plus	<-	matrix(coord.plus,nrow=2,ncol=2,byrow=T)
	for (i in 1:2){
		text(mat.plus[i,1],mat.plus[i,2],"+",cex=4)
	}
#box()
mtext("(c)",adj=0,cex=2,line=-5)

###### Scenario 4: Strong Intraspecific competition 
plot(0,0,pch=".",col="white",xlab="",ylab="",xlim=c(0,1),ylim=c(0,1),axes=F)
text(0.5,0.775,"1",cex=4, col="black")
text(0.2,0.225,"2",cex=4, col="blue")
text(0.8,0.225,"3",cex=4, col="red")
arrows(0.2,0.3,0.45,0.7,lwd=2,length=0.15)
arrows(0.4,0.7,0.15,0.30,lwd=2,length=0.15)
arrows(0.3,0.25,0.69,0.25,lwd=2,length=0.15)
arrows(0.69,0.2,0.3,0.2,lwd=2,length=0.15)
arrows(0.8,0.3,0.55,0.7,lwd=2,length=0.15)
arrows(0.6,0.7,0.85,0.3,lwd=2,length=0.15)
plotcircle(mid=c(0.5,0.85),r=0.1,from=315*pi/180,to=225*pi/180,arrow=F,lwd=6)
Arrowhead(x0=0.43,y0=0.8,angle=317,arr.lwd=3)
plotcircle(mid=c(0.12,0.15),r=0.1,from=90*pi/180,to=360*pi/180,arrow=F,lwd=6,lcol="blue")
Arrowhead(x0=0.22,y0=0.15,angle=90,arr.lwd=3,lcol="blue")
plotcircle(mid=c(0.88,0.15),r=0.1,from=180*pi/180,to=90*pi/180,arrow=F,lwd=6,lcol="red")
Arrowhead(x0=0.78,y0=0.15,angle=90,arr.lwd=3,lcol="red")
coord.min	<-	c(0.5,0.85,0.12,0.15,0.88,0.15,0.4,0.5,0.2,0.5,0.5,0.28,0.5,0.15)
mat.min		<-	matrix(coord.min,nrow=7,ncol=2,byrow=T)
for (i in 1:7){
  if(i==2){my.col<-"blue"}else if(i==3){my.col<-"red"}else{my.col<-"black"}    
  text(mat.min[i,1],mat.min[i,2],"-",cex=4, col=my.col)
}
coord.plus	<-	c(0.6,0.5,0.8,0.5)	
mat.plus	<-	matrix(coord.plus,nrow=2,ncol=2,byrow=T)
for (i in 1:2){
  text(mat.plus[i,1],mat.plus[i,2],"+",cex=4)
}
#box()
mtext("(d)",adj=0,cex=2,line=-5)

par(mar=c(3,3,1,1),mgp=c(1.5,0.5,0))
plot(1:70,ts.mat[1,],xlab="Time",ylab="3 species population sizes",cex.axis=1.2,cex.lab=1.2,pch=19,type="b",ylim=c(0,max(ts.mat)))
points(1:70,ts.mat[2,],pch=19,col="blue",type="b")
points(1:70,ts.mat[3,],pch=19,col="red",type="b")
mtext(paste0("VP=",round(ex1.stab[[1]],2),"\n","MR=",round(ex1.stab[[2]][1],2),"\n","VR=",round(ex1.stab[[2]][2],2),"\n","R=",round(ex1.stab[[3]],2)),at=55,line=-6,adj=0)

plot(1:70,ts.mat3[1,],xlab="Time",ylab="3 species population sizes",cex.axis=1.2,cex.lab=1.2,pch=19,type="b",ylim=c(0,max(ts.mat3)))
points(1:70,ts.mat3[2,],pch=19,col="blue",type="b")
points(1:70,ts.mat3[3,],pch=19,col="red",type="b")
mtext(paste0("VP=",round(ex3.stab[[1]],2),"\n","MR=",round(ex3.stab[[2]][1],2),"\n","VR=",round(ex3.stab[[2]][2],2),"\n","R=",round(ex3.stab[[3]],2)),at=55,line=-6,adj=0)
plot(1:70,ts.mat4[1,],xlab="Time",ylab="3 species population sizes",cex.axis=1.2,cex.lab=1.2,pch=19,type="b",ylim=c(0,max(ts.mat4)))
points(1:70,ts.mat4[2,],pch=19,col="blue",type="b")
points(1:70,ts.mat4[3,],pch=19,col="red",type="b")
mtext(paste0("VP=",round(ex4.stab[[1]],2),"\n","MR=",round(ex4.stab[[2]][1],2),"\n","VR=",round(ex4.stab[[2]][2],2),"\n","R=",round(ex4.stab[[3]],2)),at=55,line=-6,adj=0)

plot(1:70,ts.mat2[1,],xlab="Time",ylab="3 species population sizes",cex.axis=1.2,cex.lab=1.2,pch=19,type="b",ylim=c(0,max(ts.mat2)))
points(1:70,ts.mat2[2,],pch=19,col="blue",type="b")
points(1:70,ts.mat2[3,],pch=19,col="red",type="b")
mtext(paste0("VP=",round(ex2.stab[[1]],2),"\n","MR=",round(ex2.stab[[2]][1],2),"\n","VR=",round(ex2.stab[[2]][2],2),"\n","R=",round(ex2.stab[[3]],2)),at=55,line=-2,adj=0,side=1)


dev.off()




