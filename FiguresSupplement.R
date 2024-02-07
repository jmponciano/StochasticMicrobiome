rm(list=ls())

library(MASS)
library(shape)
source('CLS.R')
source("Poisson-MAR-ftns.R")
source("error.barfunc.R")
load("examples_for_plot.RData")


#################  THERE ARE TWO VERSIONS OF THIS FIGURE, I USED THE SECOND ONE #########
#################  THE CODE FOR THE SECOND VERSION STARTS ON LINE 145 BELOW #############
################# ---------------- Figures for the supplement ###########################
################# ---------------- Example 1 ----------------  ##########################
days	<-	70
cols	<-	c("black","blue","red")

tiff("Scenario1.tiff",height=25,width=25,units="cm",res=300,compression="lzw",type="cairo")
par(mfrow=c(2,2),mar=c(4,4,1,1),mgp=c(2.5,1,0))

plot(1:days,ts.mat[1,],pch=19,type="b",col=cols[1],xlab="Time",ylab="3 Spp abundances",ylim=c(0,max(ts.mat)),cex.axis=1.5,cex.lab=1.5)
for (i in 2:3){
  points(1:days,ts.mat[i,],pch=19,type="b",col=cols[i])
}
#legend("topleft",legend=c("Species A","Species B","Species C"),pch=19,col=cols)

plot(1:days,prop.mat[1,],pch=19,type="b",col=cols[1],xlab="Time",ylab="3 Spp proportions",ylim=c(0,max(prop.mat)),cex.axis=1.5,cex.lab=1.5)
for (i in 2:3){
  points(1:days,prop.mat[i,],pch=19,type="b",col=cols[i])
}

plot(1:9,true.B,pch=19,xlab="Interaction Coefficient",ylab="B",xaxt="n",ylim=c(-1.5,1.5),cex.lab=1.5,cex.axis=1.5,cex=2)
error.bar(1:9,N.mles$B,upper=N.boot$B.U-N.mles$B,lower=N.mles$B-N.boot$B.L)
points(1:9,N.mles$B,pch=19,col="grey",cex=2)
axis(side=1,at=1:9,labels=expression("B"[AA],"B"[BA],"B"[CA],"B"[AB],"B"[BB],"B"[CB],"B"[AC],"B"[BC],"B"[CC]),cex.axis=1.5)
abline(h=0,lty=2,lwd=2,col="gray")
plot(1:9,true.B,pch=19,xlab="Interaction Coefficient",xaxt="n",ylim=c(-1.5,1.5),cex.lab=1.5,cex.axis=1.5,ylab="B",cex=2)
error.bar(1:9,p.mles$B,upper=p.boot$B.U-p.mles$B,lower=p.mles$B-p.boot$B.L)
points(1:9,p.mles$B,pch=19,col="grey",cex=2)
axis(side=1,at=1:9,labels=expression("B"[AA],"B"[BA],"B"[CA],"B"[AB],"B"[BB],"B"[CB],"B"[AC],"B"[BC],"B"[CC]),cex.axis=1.5)
abline(h=0,lty=2,lwd=2,col="gray")
legend("bottomright",legend=c("True B","Estimated B"),pch=19,col=c("black","grey"))
dev.off()

################# ---------------- Example 2 ----------------  ###########################

days	<-	70
cols	<-	c("black","blue","red")

tiff("Scenario-old2.tiff",height=25,width=25,units="cm",res=300,compression="lzw",type="cairo")
par(mfrow=c(2,2),mar=c(4,4,1,1),mgp=c(2.5,1,0))

plot(1:days,ts.mat2[1,],pch=19,type="b",col=cols[1],xlab="Time",ylab="3 Spp abundances",ylim=c(0,max(ts.mat2)),cex.axis=1.5,cex.lab=1.5)
for (i in 2:3){
  points(1:days,ts.mat2[i,],pch=19,type="b",col=cols[i])
}
#legend("topleft",legend=c("Species A","Species B","Species C"),pch=19,col=cols)

plot(1:days,prop.mat2[1,],pch=19,type="b",col=cols[1],xlab="Time",ylab="3 Spp proportions",ylim=c(0,1),cex.axis=1.5,cex.lab=1.5)
for (i in 2:3){
  points(1:days,prop.mat2[i,],pch=19,type="b",col=cols[i])
}

plot(1:9,true.B2,pch=19,xlab="Interaction Coefficient",ylab="B",xaxt="n",cex=2,ylim=c(-1.5,1.5),cex.lab=1.5,cex.axis=1.5)
points(1:9,N.mles2$B,pch=19,col="grey",cex=2)
axis(side=1,at=1:9,labels=expression("B"[AA],"B"[AB],"B"[AC],"B"[BA],"B"[BB],"B"[BC],"B"[CA],"B"[CB],"B"[CC]),cex.axis=1.5)
abline(h=0,lty=2,lwd=2,col="gray")
error.bar(1:9,N.mles2$B,upper=N.boot2$B.U-N.mles2$B,lower=N.mles2$B-N.boot2$B.L)
plot(1:9,true.B2,pch=19,xlab="Interaction Coefficient",xaxt="n",cex=2,ylim=c(-1.5,1.5),cex.lab=1.5,cex.axis=1.5,ylab="B")
points(1:9,p.mles2$B,pch=19,col="grey",cex=2)
axis(side=1,at=1:9,labels=expression("B"[AA],"B"[AB],"B"[AC],"B"[BA],"B"[BB],"B"[BC],"B"[CA],"B"[CB],"B"[CC]),cex.axis=1.5)
abline(h=0,lty=2,lwd=2,col="gray")
error.bar(1:9,p.mles2$B,upper=p.boot2$B.U-p.mles2$B,lower=p.mles2$B-p.boot2$B.L)
legend("bottomright",legend=c("True B","Estimated B"),pch=19,col=c("black","grey"))
dev.off()

################# ---------------- Example 3 ----------------  ###########################
days	<-	70
cols	<-	c("black","blue","red")

tiff("Senario-old3.tiff",height=25,width=25,units="cm",res=300,compression="lzw",type="cairo")
par(mfrow=c(2,2),mar=c(4,4,1,1),mgp=c(2.5,1,0))

plot(1:days,ts.mat3[1,],pch=19,type="b",col=cols[1],xlab="Time",ylab="3 Spp abundances",ylim=c(0,max(ts.mat3)),cex.axis=1.5,cex.lab=1.5)
for (i in 2:3){
  points(1:days,ts.mat3[i,],pch=19,type="b",col=cols[i])
}
#legend("topleft",legend=c("Species A","Species B","Species C"),pch=19,col=cols)

plot(1:days,prop.mat3[1,],pch=19,type="b",col=cols[1],xlab="Time",ylab="3 Spp proportions",ylim=c(0,1),cex.axis=1.5,cex.lab=1.5)
for (i in 2:3){
  points(1:days,prop.mat3[i,],pch=19,type="b",col=cols[i])
}

plot(1:9,true.B3,pch=19,xlab="Interaction Coefficient",ylab="B",xaxt="n",cex=2,ylim=c(-1.5,1.5),cex.lab=1.5,cex.axis=1.5)
points(1:9,N.mles3$B,pch=19,col="grey",cex=2)
axis(side=1,at=1:9,labels=expression("B"[AA],"B"[AB],"B"[AC],"B"[BA],"B"[BB],"B"[BC],"B"[CA],"B"[CB],"B"[CC]),cex.axis=1.5)
abline(h=0,lty=2,lwd=2,col="gray")
error.bar(1:9,N.mles3$B,upper=N.boot3$B.U-N.mles3$B,lower=N.mles3$B-N.boot3$B.L)
plot(1:9,true.B3,pch=19,xlab="Interaction Coefficient",xaxt="n",cex=2,ylim=c(-1.5,1.5),cex.lab=1.5,cex.axis=1.5,ylab="B")
points(1:9,p.mles3$B,pch=19,col="grey",cex=2)
axis(side=1,at=1:9,labels=expression("B"[AA],"B"[AB],"B"[AC],"B"[BA],"B"[BB],"B"[BC],"B"[CA],"B"[CB],"B"[CC]),cex.axis=1.5)
abline(h=0,lty=2,lwd=2,col="gray")
error.bar(1:9,p.mles3$B,upper=p.boot3$B.U-p.mles3$B,lower=p.mles3$B-p.boot3$B.L)
legend("bottomright",legend=c("True B","Estimated B"),pch=19,col=c("black","grey"))
dev.off()

################# ---------------- Example 4 ----------------  ###########################
days	<-	70
cols	<-	c("black","blue","red")

tiff("Scenario-old4.tiff",height=25,width=25,units="cm",res=300,compression="lzw",type="cairo")
par(mfrow=c(2,2),mar=c(4,4,1,1),mgp=c(2.5,1,0))

plot(1:days,ts.mat4[1,],pch=19,type="b",col=cols[1],xlab="Time",ylab="3 Spp abundances",ylim=c(0,max(ts.mat4)),cex.axis=1.5,cex.lab=1.5)
for (i in 2:3){
  points(1:days,ts.mat4[i,],pch=19,type="b",col=cols[i])
}
#legend("topleft",legend=c("Species A","Species B","Species C"),pch=19,col=cols)

plot(1:days,prop.mat4[1,],pch=19,type="b",col=cols[1],xlab="Time",ylab="3 Spp proportions",ylim=c(0,1),cex.axis=1.5,cex.lab=1.5)
for (i in 2:3){
  points(1:days,prop.mat4[i,],pch=19,type="b",col=cols[i])
}

plot(1:9,true.B4,pch=19,xlab="Interaction Coefficient",ylab="B",xaxt="n",cex=2,ylim=c(-1.5,1.5),cex.lab=1.5,cex.axis=1.5)
points(1:9,N.mles4$B,pch=19,col="grey",cex=2)
axis(side=1,at=1:9,labels=expression("B"[AA],"B"[AB],"B"[AC],"B"[BA],"B"[BB],"B"[BC],"B"[CA],"B"[CB],"B"[CC]),cex.axis=1.5)
abline(h=0,lty=2,lwd=2,col="gray")
error.bar(1:9,N.mles4$B,upper=N.boot4$B.U-N.mles4$B,lower=N.mles4$B-N.boot4$B.L)
plot(1:9,true.B4,pch=19,xlab="Interaction Coefficient",xaxt="n",cex=2,ylim=c(-1.5,1.5),cex.lab=1.5,cex.axis=1.5,ylab="B")
points(1:9,p.mles4$B,pch=19,col="grey",cex=2)
axis(side=1,at=1:9,labels=expression("B"[AA],"B"[AB],"B"[AC],"B"[BA],"B"[BB],"B"[BC],"B"[CA],"B"[CB],"B"[CC]),cex.axis=1.5)
abline(h=0,lty=2,lwd=2,col="gray")
error.bar(1:9,p.mles4$B,upper=p.boot4$B.U-p.mles4$B,lower=p.mles4$B-p.boot4$B.L)
legend("bottomright",legend=c("True B","Estimated B"),pch=19,col=c("black","grey"))
dev.off()







######## -------- Supplementary Figure Alternative with the arrow diagrams inside!!!! ------###########
######## -------- Supplementary Figure Alternative with the arrow diagrams inside!!!! ------###########
######## -------- Supplementary Figure Alternative with the arrow diagrams inside!!!! ------###########
######## -------- Supplementary Figure Alternative with the arrow diagrams inside!!!! ------###########
######## -------- Supplementary Figure Alternative with the arrow diagrams inside!!!! ------###########
######## -------- Supplementary Figure Alternative with the arrow diagrams inside!!!! ------###########


################# ---------------- Example 1 ----------------  ###########################

days	<-	70
cols	<-	c("black","blue","red")


tiff("Scenario-old1-new1.tiff",height=25,width=25,units="cm",res=300,compression="lzw",type="cairo")
par(mfrow=c(2,2),mar=c(4,4,1,1),mgp=c(2.5,1,0))

plot(1:days,exp(ts.mat[1,]),pch=19,type="b",col=cols[1],xlab="Time",ylab="3 Spp abundances",ylim=c(0,max(exp(ts.mat))),cex.axis=1.5,cex.lab=1.5)
for (i in 2:3){
  points(1:days,exp(ts.mat[i,]),pch=19,type="b",col=cols[i])
}
#legend("topleft",legend=c("Species A","Species B","Species C"),pch=19,col=cols)

plot(1:days,prop.mat[1,],pch=19,type="b",col=cols[1],xlab="Time",ylab="3 Spp proportions",ylim=c(0,1),cex.axis=1.5,cex.lab=1.5)
for (i in 2:3){
  points(1:days,prop.mat[i,],pch=19,type="b",col=cols[i])
}

boxplot(Bmattots/truthmat1, ylim=c(-10,10),xlab="Interaction Coefficient", ylab=expression("Estimated B"/"True B"),xaxt="n",cex.axis=1.5,cex.lab=1.5)
abline(h=1,col="grey",lty=2,lwd=2)
axis(side=1,at=1:9,labels=expression("B"["11"],"B"["21"],"B"["31"],"B"["12"],"B"["22"],"B"["32"],"B"["13"],"B"["23"],"B"["33"]),cex.axis=1.5)
boxplot(Bmatprops/truthmat1, ylim=c(-10,10),xlab="Interaction Coefficient", ylab=expression("Estimated B"/"True B"),xaxt="n",cex.axis=1.5,cex.lab=1.5)
abline(h=1,col="grey",lty=2,lwd=2)
axis(side=1,at=1:9,labels=expression("B"["11"],"B"["21"],"B"["31"],"B"["12"],"B"["22"],"B"["32"],"B"["13"],"B"["23"],"B"["33"]),cex.axis=1.5)

par(new=T,fig=c(0.675,0.875,0.69,0.9),mar=c(0,0,0,0),oma=c(0,1,0,0))
plot(1,1,pch=".",col="white",xaxt="n",yaxt="n",xlim=c(0,1),ylim=c(0,1),bty="n")
text(0.5,0.775,"1",cex=1.5, col="black")
text(0.18,0.225,"2",cex=1.5, col="blue")
text(0.82,0.225,"3",cex=1.5, col="red")
arrows(0.2,0.3,0.45,0.7,lwd=1,length=0.13)
arrows(0.4,0.7,0.15,0.30,lwd=1,length=0.13)
arrows(0.3,0.25,0.69,0.25,lwd=1,length=0.13)
arrows(0.69,0.2,0.3,0.2,lwd=1,length=0.13)
arrows(0.8,0.3,0.55,0.7,lwd=1,length=0.13)
arrows(0.6,0.7,0.85,0.3,lwd=1,length=0.13)
plotcircle(mid=c(0.5,0.85),r=0.1,from=315*pi/180,to=225*pi/180,arrow=F,lwd=1)
Arrowhead(x0=0.4275,y0=0.7875,angle=310,arr.lwd=0.5)
plotcircle(mid=c(0.12,0.15),r=0.1,from=90*pi/180,to=360*pi/180,arrow=F,lwd=1,lcol="blue")
Arrowhead(x0=0.22,y0=0.15,angle=90,arr.lwd=0.5,lcol="blue")
plotcircle(mid=c(0.88,0.15),r=0.1,from=180*pi/180,to=90*pi/180,arrow=F,lwd=1,lcol="red")
Arrowhead(x0=0.78,y0=0.15,angle=90,arr.lwd=0.5,lcol="red")
coord.min	<-	c(0.5,0.85,0.12,0.15,0.88,0.15,0.4,0.5,0.2,0.5,0.5,0.28,0.5,0.15)
mat.min		<-	matrix(coord.min,nrow=7,ncol=2,byrow=T)
for (i in 1:7){
  if(i==2){my.col<-"blue"}else if(i==3){my.col<-"red"}else{my.col<-"black"}
  text(mat.min[i,1],mat.min[i,2],"-",cex=4, col=my.col)
}
coord.plus	<-	c(0.6,0.5,0.8,0.5)	
mat.plus	<-	matrix(coord.plus,nrow=2,ncol=2,byrow=T)
for (i in 1:2){
  text(mat.plus[i,1],mat.plus[i,2],"+",cex=2)
}

dev.off()



################# ---------------- Example 2 ----------------  ###########################

tiff("Scenario-old2-new4.tiff",height=25,width=25,units="cm",res=300,compression="lzw",type="cairo")
par(mfrow=c(2,2),mar=c(4,4,1,1),mgp=c(2.5,1,0))

plot(1:days,exp(ts.mat2[1,]),pch=19,type="b",col=cols[1],xlab="Time",ylab="3 Spp abundances",ylim=c(0,max(exp(ts.mat2))),cex.axis=1.5,cex.lab=1.5)
for (i in 2:3){
  points(1:days,exp(ts.mat2[i,]),pch=19,type="b",col=cols[i])
}
#legend("topleft",legend=c("Species A","Species B","Species C"),pch=19,col=cols)

plot(1:days,prop.mat2[1,],pch=19,type="b",col=cols[1],xlab="Time",ylab="3 Spp proportions",ylim=c(0,1),cex.axis=1.5,cex.lab=1.5)
for (i in 2:3){
  points(1:days,prop.mat2[i,],pch=19,type="b",col=cols[i])
}

boxplot(Bmattots2/truthmat2, ylim=c(-10,10),xlab="Interaction Coefficient", ylab=expression("Estimated B"/"True B"),xaxt="n",cex.axis=1.5,cex.lab=1.5)
abline(h=1,col="grey",lty=2,lwd=2)
axis(side=1,at=1:9,labels=expression("B"["11"],"B"["21"],"B"["31"],"B"["12"],"B"["22"],"B"["32"],"B"["13"],"B"["23"],"B"["33"]),cex.axis=1.5)
boxplot(Bmatprops2/truthmat2, ylim=c(-10,10),xlab="Interaction Coefficient", ylab=expression("Estimated B"/"True B"),xaxt="n",cex.axis=1.5,cex.lab=1.5)
abline(h=1,col="grey",lty=2,lwd=2)
axis(side=1,at=1:9,labels=expression("B"["11"],"B"["21"],"B"["31"],"B"["12"],"B"["22"],"B"["32"],"B"["13"],"B"["23"],"B"["33"]),cex.axis=1.5)

par(new=T,fig=c(0.575,0.775,0.785,0.995),mar=c(0,0,0,0),oma=c(0,1,0,0))
plot(1,1,pch=".",col="white",xaxt="n",yaxt="n",xlim=c(0,1),ylim=c(0,1),bty="n")
text(0.5,0.775,"1",cex=1.5,col="black")
text(0.18,0.225,"2",cex=1.5,col="blue")
text(0.82,0.225,"3",cex=1.5,col="red")
arrows(0.2,0.3,0.45,0.7,lwd=1,length=0.13)
arrows(0.4,0.7,0.15,0.30,lwd=1,length=0.13)
arrows(0.3,0.25,0.69,0.25,lwd=1,length=0.13)
arrows(0.69,0.2,0.3,0.2,lwd=1,length=0.13)
arrows(0.8,0.3,0.55,0.7,lwd=1,length=0.13)
arrows(0.6,0.7,0.85,0.3,lwd=1,length=0.13)
plotcircle(mid=c(0.5,0.85),r=0.1,from=315*pi/180,to=225*pi/180,arrow=F,lwd=3)
Arrowhead(x0=0.4275,y0=0.7875,angle=310,arr.lwd=0.5)
plotcircle(mid=c(0.12,0.15),r=0.1,from=90*pi/180,to=360*pi/180,arrow=F,lwd=3,lcol="blue")
Arrowhead(x0=0.22,y0=0.15,angle=90,arr.lwd=0.5,lcol="blue")
plotcircle(mid=c(0.88,0.15),r=0.1,from=180*pi/180,to=90*pi/180,arrow=F,lwd=3,lcol="red")
Arrowhead(x0=0.78,y0=0.15,angle=90,arr.lwd=0.5,lcol="red")
coord.min	<-	c(0.5,0.85,0.12,0.15,0.88,0.15,0.4,0.5,0.2,0.5,0.5,0.28,0.5,0.15)
mat.min		<-	matrix(coord.min,nrow=7,ncol=2,byrow=T)
for (i in 1:7){
  if(i==2){my.col<-"blue"}else if(i==3){my.col<-"red"}else{my.col<-"black"}
  text(mat.min[i,1],mat.min[i,2],"-",cex=4, col=my.col)
}
coord.plus	<-	c(0.6,0.5,0.8,0.5)	
mat.plus	<-	matrix(coord.plus,nrow=2,ncol=2,byrow=T)
for (i in 1:2){
  text(mat.plus[i,1],mat.plus[i,2],"+",cex=2)
}

dev.off()

################# ---------------- Example 3 ----------------  ###########################
tiff("Scenario-old3-new2.tiff",height=25,width=25,units="cm",res=300,compression="lzw",type="cairo")
par(mfrow=c(2,2),mar=c(4,4,1,1),mgp=c(2.5,1,0))

plot(1:days,exp(ts.mat3[1,]),pch=19,type="b",col=cols[1],xlab="Time",ylab="3 Spp abundances",ylim=c(0,max(exp(ts.mat3))),cex.axis=1.5,cex.lab=1.5)
for (i in 2:3){
  points(1:days,exp(ts.mat3[i,]),pch=19,type="b",col=cols[i])
}
#legend("topleft",legend=c("Species A","Species B","Species C"),pch=19,col=cols)

plot(1:days,prop.mat3[1,],pch=19,type="b",col=cols[1],xlab="Time",ylab="3 proportions",ylim=c(0,1),cex.axis=1.5,cex.lab=1.5)
for (i in 2:3){
  points(1:days,prop.mat3[i,],pch=19,type="b",col=cols[i])
}

boxplot(Bmattots3/truthmat3, ylim=c(-10,10),xlab="Interaction Coefficient", ylab=expression("Estimated B"/"True B"),xaxt="n",cex.axis=1.5,cex.lab=1.5)
abline(h=1,col="grey",lty=2,lwd=2)
axis(side=1,at=1:9,labels=expression("B"["11"],"B"["21"],"B"["31"],"B"["12"],"B"["22"],"B"["32"],"B"["13"],"B"["23"],"B"["33"]),cex.axis=1.5)
boxplot(Bmatprops3/truthmat3, ylim=c(-10,10),xlab="Interaction Coefficient", ylab=expression("Estimated B"/"True B"),xaxt="n",cex.axis=1.5,cex.lab=1.5)
abline(h=1,col="grey",lty=2,lwd=2)
axis(side=1,at=1:9,labels=expression("B"["11"],"B"["21"],"B"["31"],"B"["12"],"B"["22"],"B"["32"],"B"["13"],"B"["23"],"B"["33"]),cex.axis=1.5)

par(new=T,fig=c(0.575,0.775,0.69,0.9),mar=c(0,0,0,0),oma=c(0,1,0,0))
plot(1,1,pch=".",col="white",xaxt="n",yaxt="n",xlim=c(0,1),ylim=c(0,1),bty="n")
text(0.5,0.775,"1",cex=1.5, col="black")
text(0.18,0.225,"2",cex=1.5, col="blue")
text(0.82,0.225,"3",cex=1.5, col="red")
arrows(0.2,0.3,0.45,0.7,lwd=1,length=0.13)
arrows(0.4,0.7,0.15,0.30,lwd=1,length=0.13)
arrows(0.3,0.25,0.69,0.25,lwd=1,length=0.13)
arrows(0.69,0.2,0.3,0.2,lwd=1,length=0.13)
arrows(0.8,0.3,0.55,0.7,lwd=1,length=0.13)
arrows(0.6,0.7,0.85,0.3,lwd=1,length=0.13)
plotcircle(mid=c(0.5,0.85),r=0.1,from=315*pi/180,to=225*pi/180,arrow=F,lwd=1)
Arrowhead(x0=0.4275,y0=0.7875,angle=310,arr.lwd=0.5)
plotcircle(mid=c(0.12,0.15),r=0.1,from=90*pi/180,to=360*pi/180,arrow=F,lwd=1, lcol="blue")
Arrowhead(x0=0.22,y0=0.15,angle=90,arr.lwd=0.5,lcol="blue")
plotcircle(mid=c(0.88,0.15),r=0.1,from=180*pi/180,to=90*pi/180,arrow=F,lwd=3,lcol="red")
Arrowhead(x0=0.78,y0=0.15,angle=90,arr.lwd=0.5,lcol="red")
coord.min	<-	c(0.5,0.85,0.12,0.15,0.88,0.15,0.4,0.5,0.2,0.5,0.5,0.28,0.5,0.15)
mat.min		<-	matrix(coord.min,nrow=7,ncol=2,byrow=T)
for (i in 1:7){
  if(i==2){my.col<-"blue"}else if(i==3){my.col<-"red"}else{my.col<-"black"}
  text(mat.min[i,1],mat.min[i,2],"-",cex=4, col=my.col)
}
coord.plus	<-	c(0.6,0.5,0.8,0.5)	
mat.plus	<-	matrix(coord.plus,nrow=2,ncol=2,byrow=T)
for (i in 1:2){
  text(mat.plus[i,1],mat.plus[i,2],"+",cex=2)
}


dev.off()


################# ---------------- Example 4 ----------------  ###########################

tiff("Scenario-old4-new3.tiff",height=25,width=25,units="cm",res=300,compression="lzw",type="cairo")
par(mfrow=c(2,2),mar=c(4,4,1,1),mgp=c(2.5,1,0))

plot(1:days,exp(ts.mat4[1,]),pch=19,type="b",col=cols[1],xlab="Time",ylab="3 Spp abundances",ylim=c(0,max(exp(ts.mat4))),cex.axis=1.5,cex.lab=1.5)
for (i in 2:3){
  points(1:days,exp(ts.mat4[i,]),pch=19,type="b",col=cols[i])
}
#legend("topleft",legend=c("Species A","Species B","Species C"),pch=19,col=cols)

plot(1:days,prop.mat4[1,],pch=19,type="b",col=cols[1],xlab="Time",ylab="3 Spp proportions",ylim=c(0,1),cex.axis=1.5,cex.lab=1.5)
for (i in 2:3){
  points(1:days,prop.mat4[i,],pch=19,type="b",col=cols[i])
}

boxplot(Bmattots4/truthmat4, ylim=c(-10,10),xlab="Interaction Coefficient", ylab=expression("Estimated B"/"True B"),xaxt="n",cex.axis=1.5,cex.lab=1.5)
abline(h=1,col="grey",lty=2,lwd=2)
axis(side=1,at=1:9,labels=expression("B"["11"],"B"["21"],"B"["31"],"B"["12"],"B"["22"],"B"["32"],"B"["13"],"B"["23"],"B"["33"]),cex.axis=1.5)
boxplot(Bmatprops4/truthmat4, ylim=c(-10,10),xlab="Interaction Coefficient", ylab=expression("Estimated B"/"True B"),xaxt="n",cex.axis=1.5,cex.lab=1.5)
abline(h=1,col="grey",lty=2,lwd=2)
axis(side=1,at=1:9,labels=expression("B"["11"],"B"["21"],"B"["31"],"B"["12"],"B"["22"],"B"["32"],"B"["13"],"B"["23"],"B"["33"]),cex.axis=1.5)

par(new=T,fig=c(0.575,0.775,0.69,0.9),mar=c(0,0,0,0),oma=c(0,1,0,0))
plot(1,1,pch=".",col="white",xaxt="n",yaxt="n",xlim=c(0,1),ylim=c(0,1),bty="n")
text(0.5,0.775,"1",cex=1.5, col="black")
text(0.18,0.225,"2",cex=1.5, col="blue")
text(0.82,0.225,"3",cex=1.5, col="red")
arrows(0.2,0.3,0.45,0.7,lwd=3,length=0.13)
arrows(0.4,0.7,0.15,0.30,lwd=1,length=0.13)
arrows(0.3,0.25,0.69,0.25,lwd=3,length=0.13)
arrows(0.69,0.2,0.3,0.2,lwd=1,length=0.13)
arrows(0.8,0.3,0.55,0.7,lwd=1,length=0.13)
arrows(0.6,0.7,0.85,0.3,lwd=1,length=0.13)
plotcircle(mid=c(0.5,0.85),r=0.1,from=315*pi/180,to=225*pi/180,arrow=F,lwd=1)
Arrowhead(x0=0.4275,y0=0.7875,angle=310,arr.lwd=0.5)
plotcircle(mid=c(0.12,0.15),r=0.1,from=90*pi/180,to=360*pi/180,arrow=F,lwd=1,lcol="blue")
Arrowhead(x0=0.22,y0=0.15,angle=90,arr.lwd=0.5,lcol="blue")
plotcircle(mid=c(0.88,0.15),r=0.1,from=180*pi/180,to=90*pi/180,arrow=F,lwd=1,lcol="red")
Arrowhead(x0=0.78,y0=0.15,angle=90,arr.lwd=0.5, lcol="red")
coord.min	<-	c(0.5,0.85,0.12,0.15,0.88,0.15,0.4,0.5,0.2,0.5,0.5,0.28,0.5,0.15)
mat.min		<-	matrix(coord.min,nrow=7,ncol=2,byrow=T)
for (i in 1:7){
  if(i==2){my.col<-"blue"}else if(i==3){my.col<-"red"}else{my.col<-"black"}
  text(mat.min[i,1],mat.min[i,2],"-",cex=4, col=my.col)
}
coord.plus	<-	c(0.6,0.5,0.8,0.5)	
mat.plus	<-	matrix(coord.plus,nrow=2,ncol=2,byrow=T)
for (i in 1:2){
  text(mat.plus[i,1],mat.plus[i,2],"+",cex=2)
}

dev.off()




