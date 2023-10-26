#lib
library(stringr)
library(tidyr)
library(dplyr)
library(plyr)
library(ggplot2)
library(readxl)
library(reshape2)
library(tables)
library(doBy)
library(readxl)



gt.dat <- read_excel(path="IMPACC_plot_06072022.xlsx", sheet="eq5d5l")

gt.dat$visitnum <- rep(NA,nrow(gt.dat))
gt.dat$visitnum[gt.dat$VisitTypeID==30] <- 1
gt.dat$visitnum[gt.dat$VisitTypeID==40] <- 2
gt.dat$visitnum[gt.dat$VisitTypeID==50] <- 3
gt.dat$visitnum[gt.dat$VisitTypeID==60] <- 4

gt.dat$month <- rep(NA,nrow(gt.dat))
gt.dat$month[gt.dat$VisitTypeID==30] <- 3
gt.dat$month[gt.dat$VisitTypeID==40] <- 6
gt.dat$month[gt.dat$VisitTypeID==50] <- 9
gt.dat$month[gt.dat$VisitTypeID==60] <- 12


gt.dat$score <- gt.dat$eq5d5l
gt.dat$cluster <- gt.dat$class_eq5d5l


gt.dat$scorejit <- gt.dat$score + rnorm(nrow(gt.dat),mean=0,sd=0.1)
gt.dat$monthjit <- gt.dat$month + rnorm(nrow(gt.dat),mean=0,sd=0.1)
vtimes <- c(3,6,9,12)
#NA,i,j i=number of clusters, j=number of visits=4
medianmat <- matrix(NA,4,4)
for (i in 1:3) for (j in 1:4) medianmat[i,j] <- median(gt.dat$score[gt.dat$visitnum==j & gt.dat$cluster==i],na.rm=T)

upmat <- matrix(NA,4,4)
for (i in 1:3) for (j in 1:4) {
  x <- gt.dat$score[gt.dat$visitnum==j & gt.dat$cluster==i]
  upmat[i,j] <- quantile(x,probs=0.75,na.rm=T) 
}
dnmat <- matrix(NA,3,4)
for (i in 1:3) for (j in 1:4) {
  x <- gt.dat$score[gt.dat$visitnum==j & gt.dat$cluster==i]
  dnmat[i,j] <- quantile(x,probs=0.25,na.rm=T)

}


grpcol3 <- c("blue","orange","red","black")



##### Figure 1. EQ-5D-5L trajectories by group

jpeg("traj_eq5d5l.tiff",quality=100,height=480,width=960)
plot(gt.dat$month,gt.dat$score,type="n",ylab="Score",xlab="Months",main=paste())

for (i in 1:3) for (j in 1:3) {
  polygon(x=c(vtimes[j],vtimes[j+1],vtimes[j+1],vtimes[j],vtimes[j]),
          y=c(upmat[i,j],upmat[i,j+1],dnmat[i,j+1],dnmat[i,j],upmat[i,j]),
          col=grpcol3[i],density=30,border=NA)
}
for (i in 1:3) lines(vtimes,medianmat[i,],col=grpcol3[i],lwd=5)
#add legend
legend(x="bottomleft",pch=19,cex=1.5,col=grpcol3,legend=c("High","Medium","Low"))
dev.off()



##########
gt.dat <- read_excel(path="IMPACC_plot_06072022.xlsx", sheet="health")

gt.dat$visitnum <- rep(NA,nrow(gt.dat))
gt.dat$visitnum[gt.dat$VisitTypeID==30] <- 1
gt.dat$visitnum[gt.dat$VisitTypeID==40] <- 2
gt.dat$visitnum[gt.dat$VisitTypeID==50] <- 3
gt.dat$visitnum[gt.dat$VisitTypeID==60] <- 4

gt.dat$month <- rep(NA,nrow(gt.dat))
gt.dat$month[gt.dat$VisitTypeID==30] <- 3
gt.dat$month[gt.dat$VisitTypeID==40] <- 6
gt.dat$month[gt.dat$VisitTypeID==50] <- 9
gt.dat$month[gt.dat$VisitTypeID==60] <- 12


gt.dat$score <- gt.dat$health_score
gt.dat$cluster <- gt.dat$class_health


gt.dat$scorejit <- gt.dat$score + rnorm(nrow(gt.dat),mean=0,sd=0.1)
gt.dat$monthjit <- gt.dat$month + rnorm(nrow(gt.dat),mean=0,sd=0.1)
vtimes <- c(3,6,9,12)
#NA,i,j i=number of clusters, j=number of visits=4
medianmat <- matrix(NA,3,4)
for (i in 1:3) for (j in 1:4) medianmat[i,j] <- median(gt.dat$score[gt.dat$visitnum==j & gt.dat$cluster==i],na.rm=T)

upmat <- matrix(NA,3,4)
for (i in 1:3) for (j in 1:4) {
  x <- gt.dat$score[gt.dat$visitnum==j & gt.dat$cluster==i]
  upmat[i,j] <- quantile(x,probs=0.75,na.rm=T) 
}
dnmat <- matrix(NA,3,4)
for (i in 1:3) for (j in 1:4) {
  x <- gt.dat$score[gt.dat$visitnum==j & gt.dat$cluster==i]
  dnmat[i,j] <- quantile(x,probs=0.25,na.rm=T) 
}




##### Figure 2. Health recovery score trajectories by group


jpeg("traj_health.tiff",quality=100,height=480,width=960)
plot(gt.dat$month,gt.dat$score,type="n",ylab="Score",xlab="Months",main=paste())

for (i in 1:3) for (j in 1:3) {
  polygon(x=c(vtimes[j],vtimes[j+1],vtimes[j+1],vtimes[j],vtimes[j]),
          y=c(upmat[i,j],upmat[i,j+1],dnmat[i,j+1],dnmat[i,j],upmat[i,j]),
          col=grpcol3[i],density=30,border=NA)
}
for (i in 1:3) lines(vtimes,medianmat[i,],col=grpcol3[i],lwd=5)
#add legend
legend(x="bottomleft",pch=19,cex=1.5,col=grpcol3,legend=c("High","Medium","Low"))
dev.off()



###########
gt.dat <- read_excel(path="/IMPACC_plot_06072022.xlsx", sheet="psycho")

gt.dat$visitnum <- rep(NA,nrow(gt.dat))
gt.dat$visitnum[gt.dat$VisitTypeID==30] <- 1
gt.dat$visitnum[gt.dat$VisitTypeID==40] <- 2
gt.dat$visitnum[gt.dat$VisitTypeID==50] <- 3
gt.dat$visitnum[gt.dat$VisitTypeID==60] <- 4

gt.dat$month <- rep(NA,nrow(gt.dat))
gt.dat$month[gt.dat$VisitTypeID==30] <- 3
gt.dat$month[gt.dat$VisitTypeID==40] <- 6
gt.dat$month[gt.dat$VisitTypeID==50] <- 9
gt.dat$month[gt.dat$VisitTypeID==60] <- 12


gt.dat$score <- gt.dat$PROMIS_Psycho
gt.dat$cluster <- gt.dat$class_cognitive


gt.dat$scorejit <- gt.dat$score + rnorm(nrow(gt.dat),mean=0,sd=0.1)
gt.dat$monthjit <- gt.dat$month + rnorm(nrow(gt.dat),mean=0,sd=0.1)
vtimes <- c(3,6,9,12)
#NA,i,j i=number of clusters, j=number of visits=4
medianmat <- matrix(NA,3,4)
for (i in 1:3) for (j in 1:4) medianmat[i,j] <- median(gt.dat$score[gt.dat$visitnum==j & gt.dat$cluster==i],na.rm=T)

upmat <- matrix(NA,3,4)
for (i in 1:3) for (j in 1:4) {
  x <- gt.dat$score[gt.dat$visitnum==j & gt.dat$cluster==i]
  upmat[i,j] <- quantile(x,probs=0.75,na.rm=T)
}
dnmat <- matrix(NA,3,4)
for (i in 1:3) for (j in 1:4) {
  x <- gt.dat$score[gt.dat$visitnum==j & gt.dat$cluster==i]
  dnmat[i,j] <- quantile(x,probs=0.25,na.rm=T)
}


##### Figure 4. PROMIS cognitive score trajectories by group

jpeg("traj_cognitive.tiff",quality=100,height=480,width=960)
plot(gt.dat$month,gt.dat$score,type="n",ylab="Score",xlab="Months",main=paste())

for (i in 1:3) for (j in 1:3) {
  polygon(x=c(vtimes[j],vtimes[j+1],vtimes[j+1],vtimes[j],vtimes[j]),
          y=c(upmat[i,j],upmat[i,j+1],dnmat[i,j+1],dnmat[i,j],upmat[i,j]),
          col=grpcol3[i],density=30,border=NA)
}
for (i in 1:3) lines(vtimes,medianmat[i,],col=grpcol3[i],lwd=5)
#add legend
legend(x="bottomleft",pch=19,cex=1.5,col=grpcol3,legend=c("High","Medium","Low"))
dev.off()


#######
gt.dat <- read_excel(path="IMPACC_plot_06072022.xlsx", sheet="dyspnea2")

gt.dat$visitnum <- rep(NA,nrow(gt.dat))
gt.dat$visitnum[gt.dat$VisitTypeID==30] <- 1
gt.dat$visitnum[gt.dat$VisitTypeID==40] <- 2
gt.dat$visitnum[gt.dat$VisitTypeID==50] <- 3
gt.dat$visitnum[gt.dat$VisitTypeID==60] <- 4

gt.dat$month <- rep(NA,nrow(gt.dat))
gt.dat$month[gt.dat$VisitTypeID==30] <- 3
gt.dat$month[gt.dat$VisitTypeID==40] <- 6
gt.dat$month[gt.dat$VisitTypeID==50] <- 9
gt.dat$month[gt.dat$VisitTypeID==60] <- 12


gt.dat$score <- gt.dat$PROMIS_Dyspnea
gt.dat$cluster <- gt.dat$class_dyspnea

gt.dat$scorejit <- gt.dat$score + rnorm(nrow(gt.dat),mean=0,sd=0.1)
gt.dat$monthjit <- gt.dat$month + rnorm(nrow(gt.dat),mean=0,sd=0.1)
vtimes <- c(3,6,9,12)
#NA,i,j i=number of clusters, j=number of visits=4
medianmat <- matrix(NA,3,4)
for (i in 1:3) for (j in 1:4) medianmat[i,j] <- median(gt.dat$score[gt.dat$visitnum==j & gt.dat$cluster==i],na.rm=T)

upmat <- matrix(NA,3,4)
for (i in 1:3) for (j in 1:4) {
  x <- gt.dat$score[gt.dat$visitnum==j & gt.dat$cluster==i]
  upmat[i,j] <- quantile(x,probs=0.75,na.rm=T)
}
dnmat <- matrix(NA,3,4)
for (i in 1:3) for (j in 1:4) {
  x <- gt.dat$score[gt.dat$visitnum==j & gt.dat$cluster==i]
  dnmat[i,j] <- quantile(x,probs=0.25,na.rm=T)
}


##### Figure 5. PROMIS dyspnea score trajectories by group


jpeg("traj_dyspnea.tiff",quality=100,height=480,width=960)
plot(gt.dat$month,gt.dat$score,type="n",ylab="Score",xlab="Months",main=paste())

for (i in 1:3) for (j in 1:3) {
  polygon(x=c(vtimes[j],vtimes[j+1],vtimes[j+1],vtimes[j],vtimes[j]),
          y=c(upmat[i,j],upmat[i,j+1],dnmat[i,j+1],dnmat[i,j],upmat[i,j]),
          col=grpcol3[i],density=30,border=NA)
}
for (i in 1:3) lines(vtimes,medianmat[i,],col=grpcol3[i],lwd=5)
#add legend
legend(x="bottomleft",pch=19,cex=1.5,col=grpcol3,legend=c("Low","Medium","High"))
dev.off()




##no traj
#PROMIS Physical Function Score

#PROMIS Global Mental Health Score
#PROMIS Psychosocial Illness Impact Positive Score

#impact
gt.dat <- read_excel(path="PROMIS_scores_05272022.xlsx", sheet="Impact")

gt.dat$visitnum <- rep(NA,nrow(gt.dat))
gt.dat$visitnum[gt.dat$VisitTypeID==30] <- 1
gt.dat$visitnum[gt.dat$VisitTypeID==40] <- 2
gt.dat$visitnum[gt.dat$VisitTypeID==50] <- 3
gt.dat$visitnum[gt.dat$VisitTypeID==60] <- 4

gt.dat$month <- rep(NA,nrow(gt.dat))
gt.dat$month[gt.dat$VisitTypeID==30] <- 3
gt.dat$month[gt.dat$VisitTypeID==40] <- 6
gt.dat$month[gt.dat$VisitTypeID==50] <- 9
gt.dat$month[gt.dat$VisitTypeID==60] <- 12


gt.dat$score <- gt.dat$PROMIS_Impact
gt.dat$cluster <- 1

gt.dat$scorejit <- gt.dat$score + rnorm(nrow(gt.dat),mean=0,sd=0.1)
gt.dat$monthjit <- gt.dat$month + rnorm(nrow(gt.dat),mean=0,sd=0.1)
vtimes <- c(3,6,9,12)
#NA,i,j i=number of clusters, j=number of visits=4
medianmat <- matrix(NA,1,4)
for (i in 1:1) for (j in 1:4) medianmat[i,j] <- median(gt.dat$score[gt.dat$visitnum==j & gt.dat$cluster==i],na.rm=T)

upmat <- matrix(NA,1,4)
for (i in 1:1) for (j in 1:4) {
  x <- gt.dat$score[gt.dat$visitnum==j & gt.dat$cluster==i]
  upmat[i,j] <- quantile(x,probs=0.75,na.rm=T)
}
dnmat <- matrix(NA,1,4)
for (i in 1:1) for (j in 1:4) {
  x <- gt.dat$score[gt.dat$visitnum==j & gt.dat$cluster==i]
  dnmat[i,j] <- quantile(x,probs=0.25,na.rm=T)
}


#jpeg("/Users/sonicpumpkin/Documents/BCH/impacc/final2023/traj_plots/notraj_impact.tiff",quality=100,height=480,width=960)
plot(gt.dat$month,gt.dat$score,type="n",ylab="Score",xlab="Months",main=paste(),ylim=c(0,100))

for (i in 1:1) for (j in 1:3) {
  polygon(x=c(vtimes[j],vtimes[j+1],vtimes[j+1],vtimes[j],vtimes[j]),
          y=c(upmat[i,j],upmat[i,j+1],dnmat[i,j+1],dnmat[i,j],upmat[i,j]),
          col=grpcol3[i],density=30,border=NA)
}
for (i in 1:1) lines(vtimes,medianmat[i,],col=grpcol3[i],lwd=5)
#add legend
#dev.off()




#mental

gt.dat <- read_excel(path="PROMIS_scores_05272022.xlsx", sheet="Mental")

gt.dat$visitnum <- rep(NA,nrow(gt.dat))
gt.dat$visitnum[gt.dat$VisitTypeID==30] <- 1
gt.dat$visitnum[gt.dat$VisitTypeID==40] <- 2
gt.dat$visitnum[gt.dat$VisitTypeID==50] <- 3
gt.dat$visitnum[gt.dat$VisitTypeID==60] <- 4

gt.dat$month <- rep(NA,nrow(gt.dat))
gt.dat$month[gt.dat$VisitTypeID==30] <- 3
gt.dat$month[gt.dat$VisitTypeID==40] <- 6
gt.dat$month[gt.dat$VisitTypeID==50] <- 9
gt.dat$month[gt.dat$VisitTypeID==60] <- 12


gt.dat$score <- gt.dat$PROMIS_Mental
gt.dat$cluster <- 1

gt.dat$scorejit <- gt.dat$score + rnorm(nrow(gt.dat),mean=0,sd=0.1)
gt.dat$monthjit <- gt.dat$month + rnorm(nrow(gt.dat),mean=0,sd=0.1)
vtimes <- c(3,6,9,12)
#NA,i,j i=number of clusters, j=number of visits=4
medianmat <- matrix(NA,1,4)
for (i in 1:1) for (j in 1:4) medianmat[i,j] <- median(gt.dat$score[gt.dat$visitnum==j & gt.dat$cluster==i],na.rm=T)

upmat <- matrix(NA,1,4)
for (i in 1:1) for (j in 1:4) {
  x <- gt.dat$score[gt.dat$visitnum==j & gt.dat$cluster==i]
  upmat[i,j] <- quantile(x,probs=0.75,na.rm=T)
}
dnmat <- matrix(NA,1,4)
for (i in 1:1) for (j in 1:4) {
  x <- gt.dat$score[gt.dat$visitnum==j & gt.dat$cluster==i]
  dnmat[i,j] <- quantile(x,probs=0.25,na.rm=T)
}


jpeg("notraj_mental.tiff",quality=100,height=480,width=960)
plot(gt.dat$month,gt.dat$score,type="n",ylab="Score",xlab="Months",main=paste(),ylim=c(0,100))

for (i in 1:1) for (j in 1:3) {
  polygon(x=c(vtimes[j],vtimes[j+1],vtimes[j+1],vtimes[j],vtimes[j]),
          y=c(upmat[i,j],upmat[i,j+1],dnmat[i,j+1],dnmat[i,j],upmat[i,j]),
          col=grpcol3[i],density=30,border=NA)
}
for (i in 1:1) lines(vtimes,medianmat[i,],col=grpcol3[i],lwd=5)
#add legend
dev.off()



#PROMIS physical


gt.dat <- read_excel(path="PROMIS_scores_05272022.xlsx", sheet="Physical")

gt.dat$visitnum <- rep(NA,nrow(gt.dat))
gt.dat$visitnum[gt.dat$VisitTypeID==30] <- 1
gt.dat$visitnum[gt.dat$VisitTypeID==40] <- 2
gt.dat$visitnum[gt.dat$VisitTypeID==50] <- 3
gt.dat$visitnum[gt.dat$VisitTypeID==60] <- 4

gt.dat$month <- rep(NA,nrow(gt.dat))
gt.dat$month[gt.dat$VisitTypeID==30] <- 3
gt.dat$month[gt.dat$VisitTypeID==40] <- 6
gt.dat$month[gt.dat$VisitTypeID==50] <- 9
gt.dat$month[gt.dat$VisitTypeID==60] <- 12


gt.dat$score <- gt.dat$PROMIS_Physical
gt.dat$cluster <- 1

gt.dat$scorejit <- gt.dat$score + rnorm(nrow(gt.dat),mean=0,sd=0.1)
gt.dat$monthjit <- gt.dat$month + rnorm(nrow(gt.dat),mean=0,sd=0.1)
vtimes <- c(3,6,9,12)
#NA,i,j i=number of clusters, j=number of visits=4
medianmat <- matrix(NA,1,4)
for (i in 1:1) for (j in 1:4) medianmat[i,j] <- median(gt.dat$score[gt.dat$visitnum==j & gt.dat$cluster==i],na.rm=T)

upmat <- matrix(NA,1,4)
for (i in 1:1) for (j in 1:4) {
  x <- gt.dat$score[gt.dat$visitnum==j & gt.dat$cluster==i]
  upmat[i,j] <- quantile(x,probs=0.75,na.rm=T)
}
dnmat <- matrix(NA,1,4)
for (i in 1:1) for (j in 1:4) {
  x <- gt.dat$score[gt.dat$visitnum==j & gt.dat$cluster==i]
  dnmat[i,j] <- quantile(x,probs=0.25,na.rm=T)
}


jpeg("notraj_physical.tiff",quality=100,height=480,width=960)
plot(gt.dat$month,gt.dat$score,type="n",ylab="Score",xlab="Months",main=paste(),ylim=c(0,100))

for (i in 1:1) for (j in 1:3) {
  polygon(x=c(vtimes[j],vtimes[j+1],vtimes[j+1],vtimes[j],vtimes[j]),
          y=c(upmat[i,j],upmat[i,j+1],dnmat[i,j+1],dnmat[i,j],upmat[i,j]),
          col=grpcol3[i],density=30,border=NA)
}
for (i in 1:1) lines(vtimes,medianmat[i,],col=grpcol3[i],lwd=5)
#add legend
dev.off()








