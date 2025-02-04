##Figure 1
svg(paste0("output/cd4_fig1.svg"),width=10.4,height=14.9) 
layout(matrix(c(1,1,1,1,2,2,3,3,4,4,5,5,1,1,1,1,2,2,3,3,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8,7,7,7,7,7,7,8,8,8,8,8,8,9,9,10,10,10,10,11,11,12,12,12,12,9,9,10,10,10,10,11,11,12,12,12,12,13,13,13,13,13,13,14,14,14,15,15,15,13,13,13,13,13,13,16,16,16,17,17,17),nrow=8,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))

#PCA and clusters
i=is[2]
tmp=tcell.final.acute[[i]][apply(tcell.final.acute[[i]],1,function(x){sum(is.na(x))/length(x)*100})==0,]
tmpp=prcomp((tmp[,js]>0)+0,scale=F)
tmpc=rowSums(tmp[,js[1:12]]>0)

ii=cut(tmpc,breaks=seq(min(tmpc),max(tmpc),len=100),include.lowest=T)
tcolors=colorRampPalette(spctrl[rev(c(2:5,8:10))])(99)[ii]
plot(-tmpp$x[,1],tmpp$x[,2],xlab=paste0("PC1 (",round(summary(tmpp)[[6]][2,1]*100,1),"%)"),ylab=paste0("PC2 (",round(summary(tmpp)[[6]][2,2]*100,1),"%)"),xaxt='n',yaxt='n',pch=21,bg=tcolors,cex=1)
axis(1,at=seq(-3,2,1),labels=seq(-3,2,1))
axis(2,at=seq(-3,2,1),labels=seq(-3,2,1),las=2)
abline(h=0,v=0,lty=3)
legend("topleft",title="No. +",legend=seq(0,12,2),pch=21,pt.bg=tcolors[match(seq(0,12,2),tmpc)],bty='n')

par(bty='n',mai=c(0.62,0.22,0.22,0.12))

b=barplot(-tmpp$rotation[,1],horiz=T,xlab="PC1 loading",yaxt='n',xaxt='n',col=c(rep("cornflowerblue",6),rep("chartreuse3",6)))
axis(1,at=seq(0,0.4,0.1),labels=seq(0,0.4,0.1))
text(rep(0,length(b)),b,labels=rownames(tmpp$rotation),pos=4,cex=0.75)

b=barplot(tmpp$rotation[,2],horiz=T,xlab="PC2 loading",yaxt='n',xaxt='n',col=c(rep("cornflowerblue",6),rep("chartreuse3",6)),xlim=c(-0.4,0.4))
axis(1,at=seq(-0.4,0.4,0.2),labels=seq(-0.4,0.4,0.2))
#text(rep(0,length(b[1:6])),b[1:6],labels=rownames(tmpp$rotation)[1:6],pos=2)
#text(rep(0,length(b[7:12])),b[7:12],labels=rownames(tmpp$rotation)[7:12],pos=4)

#totalk2.4=kmeans(tmpp$x[,1:2],centers=2)
par(bty='n',mai=c(0.62,0.52,0.22,0.22))
z=2
plot(-tmpp$x[totalk2.4$cluster==z,1],tmpp$x[totalk2.4$cluster==z,2],xlab="PC1",ylab="PC2",xaxt='n',yaxt='n',pch=21,bg="orange",cex=0.8,xlim=c(-2,2),ylim=c(-2,2),main=paste0("Cluster-",z,": Low"))
axis(1,at=seq(-3,2,1),labels=seq(-3,2,1),pos=-2)
axis(2,at=seq(-3,2,1),labels=seq(-3,2,1),las=2,pos=-2)
z=1
plot(-tmpp$x[totalk2.4$cluster==z,1],tmpp$x[totalk2.4$cluster==z,2],xlab="PC1",ylab="PC2",xaxt='n',yaxt='n',pch=21,bg="magenta",cex=0.8,xlim=c(-2,2),ylim=c(-2,2),main=paste0("Cluster-",z,": High"))
axis(1,at=seq(-3,2,1),labels=seq(-3,2,1),pos=-2)
axis(2,at=seq(-3,2,1),labels=seq(-3,2,1),las=2,pos=-2)

#tt=totalk2.4$cluster[match(tcell.info.acute$event_id,rownames(tmpp$x))]
#tcell.info.acute$responder_c2.4=NA
#tcell.info.acute$responder_c2.4[which(tt==1)]="high"
#tcell.info.acute$responder_c2.4[which(tt==2)]="low"
#tcell.info.acute$responder_c2.4=factor(tcell.info.acute$responder_c2.4)
x
plot(density(log(tcell.info.acute$event_date[which(tcell.info.acute$responder_c2.4=="low")])),xlim=c(0,4.5),main="",xaxt='n',col="orange",xlab="Days post-admission",las=2)
polygon(density(log(tcell.info.acute$event_date[which(tcell.info.acute$responder_c2.4=="low")])),col=adjustcolor("orange",alpha.f=0.5))
lines(density(log(tcell.info.acute$event_date[which(tcell.info.acute$responder_c2.4=="high")])),col="magenta")
polygon(density(log(tcell.info.acute$event_date[which(tcell.info.acute$responder_c2.4=="high")])),col=adjustcolor("magenta",alpha.f=0.5))
axis(1,at=log(c(1,2,5,10,20,50)),labels=c(1,2,5,10,20,50),pos=0)
m=summary(glm(tcell.info.acute$responder_c2.4~log(tcell.info.acute$event_date),family="binomial"))$coefficients
legend("right",legend=paste0("P = ",round(m[2,4],3)),bty='n')
legend("topright",legend=c("Low","High"),pch=22,pt.bg=c("orange","magenta"),pt.cex=1.5,bty='n')

#RBD/vrld
j=js[1]
for(k in virvars[c(3,1)])
{
  temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tcell.final.acute[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin)) #get data
  temp$var=tcell.info.acute[match(temp$event_id,tcell.info.acute$event_id),k]	#get variable
  temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
  temp$livecd3=livecd3.final[[i]][match(temp$event_id,names(livecd3.final[[i]]))] #get stim live count
  temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
  temp$gate=as.numeric(temp$gate)
  temp$days=as.numeric(temp$days)
  temp$days_symptom=as.numeric(temp$days_symptom)
  
  temp$gate2=log2(temp$gate+0.001)  #make log-transformed variable
  temp$gate3=log2(temp$gate)  #as above but don't add minimum value
  temp$gate3[is.infinite(temp$gate3)]=NA
  
  temp$responder=temp$gate>0 #make responder
  
  temp$days_binned=tcell.info.acute$event_date_binned[match(temp$event_id,tcell.info.acute$event_id)]
  temp$days_symptom_binned=tcell.info.acute$event_date_symptom_binned[match(temp$event_id,tcell.info.acute$event_id)]
  
  temp$logdays=log(temp$days) #needed for lme
  temp$logdays_symptom=log(temp$days_symptom) #needed for lme
  temp=temp[!is.na(temp$var),]	#remove missing
  temp$responder_c2=factor(tcell.info.acute$responder_c2.4[match(temp$event_id,tcell.info.acute$event_id)],c("low","high"))
  
  lgth=length(levels(temp$responder_c2))
  tmpcols=c("orange","magenta")
  temp$col="orange"
  temp$col[temp$responder_c2=="high"]="magenta"
  
  m=gam(var~responder_c2+s(logdays,bs="cr",by=responder_c2),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$responder_c2),],random =list(id=~1),family="gaussian")
  a1=anova(m)
  aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
  aa2=sum(a1$s.table[-c(1), 2])
  
  pints=as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))
  newpvs=paste0(c("Slope","High"),"-P=",sprintf("%.1e",c(pchisq(aa1,aa2, lower.tail = F),pints)))
  plot(temp$logdays,temp$var,xlab="Days post-admission",pch=21,bg=adjustcolor("orange",alpha.f = 0.5),ylab=virnams[k],main="",type='n',xaxt='n',las=2)
  axis(1,at=log(c(1,2,3,5,10,20,30,50,100)),labels=c(1,2,3,5,10,20,30,50,100))
  points(temp$logdays,temp$var,bg=adjustcolor(temp$col,alpha.f=0.5),pch=21,cex=1)
  
  for(z in 1:lgth)
  {
    tryCatch({
      sp.cis=spline.cis(cbind(temp$logdays[which(!is.na(temp$logdays)&!is.na(temp$var)&temp$responder_c2==levels(temp$responder_c2)[z])],temp$var[which(!is.na(temp$logdays)&!is.na(temp$var)&temp$responder_c2==levels(temp$responder_c2)[z])]),B=100,alpha=0.05)
      lines(x=sp.cis$x,y=sp.cis$main.curve,col=tmpcols[z],lwd=2)
      lines(x=sp.cis$x,y=sp.cis$upper.ci,col=tmpcols[z],lty=3,lwd=1)
      lines(x=sp.cis$x,y=sp.cis$lower.ci,col=tmpcols[z],lty=3,lwd=1)
    },error=function(e){})
  }
  if(k==virvars[3])
  {
    legend("bottomright",legend=c(newpvs,"Low"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",c("magenta","orange")),pt.cex=1.5,title="Responder")		
  } else {
    legend("right",legend=c(newpvs,"Low"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",c("magenta","orange")),pt.cex=1.5,title="Responder")
  }
}		  

#TG/d28 mortality
for(k in virvars[c(12,16)])
{
  temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tcell.final.acute[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin)) #get data
  temp$var=tcell.info.acute[match(temp$event_id,tcell.info.acute$event_id),k]	#get variable
  temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
  temp$livecd3=livecd3.final[[i]][match(temp$event_id,names(livecd3.final[[i]]))] #get stim live count
  temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
  temp$gate=as.numeric(temp$gate)
  temp$days=as.numeric(temp$days)
  temp$days_symptom=as.numeric(temp$days_symptom)
  
  temp$gate2=log2(temp$gate+0.001)  #make log-transformed variable
  temp$gate3=log2(temp$gate)  #as above but don't add minimum value
  temp$gate3[is.infinite(temp$gate3)]=NA
  
  temp$responder=temp$gate>0 #make responder
  
  temp$days_binned=tcell.info.acute$event_date_binned[match(temp$event_id,tcell.info.acute$event_id)]
  temp$days_symptom_binned=tcell.info.acute$event_date_symptom_binned[match(temp$event_id,tcell.info.acute$event_id)]
  
  temp$logdays=log(temp$days) #needed for lme
  temp$logdays_symptom=log(temp$days_symptom) #needed for lme
  temp=temp[!is.na(temp$var),]	#remove missing
  temp$responder_c2=factor(tcell.info.acute$responder_c2.4[match(temp$event_id,tcell.info.acute$event_id)],c("low","high"))
  
  temp$var=factor(temp$var)
  lgth=length(levels(temp$var))
  temp$col=vircols[[k]][match(temp$var,levels(temp$var))]
  
  tmp=table(temp$responder_c2,temp$var)
  tttt=colSums(tmp)
  tmp2=tmp[2,]/colSums(tmp)*100
  m=glmer(responder_c2~var+(1|id),temp,family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),nAGQ=20)
  newpvs=sprintf("%.1e",summary(m)$coefficients[2,4])
  ymax=100
  tmp=barplot(tmp2,main="",ylab="High (%)",xaxt='n',col=unlist(lapply(vircols[[k]],function(x){rep(x,1)})),las=2,xlab=virnams[k],ylim=c(0,ymax))
  axis(1,tmp,c("FALSE","TRUE"),tick=F,las=1)
  legend("topright",legend=paste0("P=",newpvs),bty='n')
  text(tmp,rep(5,length(tmp)),labels=unlist(tttt))
  
  
  m=gam(responder_c2~var+s(logdays,bs="cr",by=var),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$responder_c2),],random =list(id=~1),family="binomial")
  a1=anova(m)
  aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
  aa2=sum(a1$s.table[-c(1), 2])
  
  pints=as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))
  newpvs=paste0(c("Slope","High"),"-P=",sprintf("%.1e",c(pchisq(aa1,aa2, lower.tail = F),pints)))
  tt=vector("list",lgth)
  ttt=vector("list",lgth)
  tttt=vector("list",lgth)
  ymax=100
  tmp=table("gate"=temp$responder_c2,temp$days_binned)
  tnams=colnames(tmp)
  for(z in 1:lgth)
  {
    tt[[z]]=table("gate"=temp$responder_c2[which(temp$var==levels(temp$var)[z])],temp$days_binned[which(temp$var==levels(temp$var)[z])])		
    tt[[z]]=tt[[z]][match(c("low","high"),rownames(tt[[z]])),]
    rownames(tt[[z]])=c("low","high")
    tt[[z]][is.na(tt[[z]])]=0
    ttt[[z]]=tt[[z]][2,]/colSums(tt[[z]])*100
    tttt[[z]]=colSums(tt[[z]])
  }
  tmp=barplot(unlist(ttt),main="",ylab="High (%)",xaxt='n',col=unlist(lapply(vircols[[k]],function(x){rep(x,3)})),las=2,xlab="Days post-admission",ylim=c(0,ymax))
  axis(1,tmp,rep(tnams,lgth),tick=F,las=1)
  legend("topright",legend=c(newpvs,"TRUE","FALSE"),bty='n',pch=c(NA,rep(22,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=virnams[k])
  text(tmp,rep(5,length(tmp)),labels=unlist(tttt))
}

dev.off()

##Figure 1
svg(paste0("output/cd8_fig1.svg"),width=10.4,height=14.9) 
layout(matrix(c(1,1,1,1,2,2,3,3,4,4,5,5,1,1,1,1,2,2,3,3,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8,7,7,7,7,7,7,8,8,8,8,8,8,9,9,10,10,10,10,11,11,12,12,12,12,9,9,10,10,10,10,11,11,12,12,12,12,13,13,13,13,13,13,14,14,14,15,15,15,13,13,13,13,13,13,16,16,16,17,17,17),nrow=8,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))

#PCA and clusters
i=is[3]
tmp=tcell.final.acute[[i]][apply(tcell.final.acute[[i]],1,function(x){sum(is.na(x))/length(x)*100})==0,]
tmpp=prcomp((tmp[,js]>0)+0,scale=F)
tmpc=rowSums(tmp[,js[1:12]]>0)

ii=cut(tmpc,breaks=seq(min(tmpc),max(tmpc),len=100),include.lowest=T)
tcolors=colorRampPalette(spctrl[rev(c(2:5,8:10))])(99)[ii]
plot(-tmpp$x[,1],tmpp$x[,2],xlab=paste0("PC1 (",round(summary(tmpp)[[6]][2,1]*100,1),"%)"),ylab=paste0("PC2 (",round(summary(tmpp)[[6]][2,2]*100,1),"%)"),xaxt='n',yaxt='n',pch=21,bg=tcolors,cex=1)
axis(1,at=seq(-3,2,1),labels=seq(-3,2,1))
axis(2,at=seq(-3,2,1),labels=seq(-3,2,1),las=2)
abline(h=0,v=0,lty=3)
legend("topleft",title="No. +",legend=seq(0,12,2),pch=21,pt.bg=tcolors[match(seq(0,12,2),tmpc)],bty='n')

par(bty='n',mai=c(0.62,0.22,0.22,0.12))

b=barplot(-tmpp$rotation[,1],horiz=T,xlab="PC1 loading",yaxt='n',xaxt='n',col=c(rep("cornflowerblue",6),rep("chartreuse3",6)))
axis(1,at=seq(0,0.4,0.1),labels=seq(0,0.4,0.1))
text(rep(0,length(b)),b,labels=rownames(tmpp$rotation),pos=4,cex=0.75)

b=barplot(tmpp$rotation[,2],horiz=T,xlab="PC2 loading",yaxt='n',xaxt='n',col=c(rep("cornflowerblue",6),rep("chartreuse3",6)),xlim=c(-0.4,0.4))
axis(1,at=seq(-0.4,0.4,0.2),labels=seq(-0.4,0.4,0.2))
#text(rep(0,length(b[1:6])),b[1:6],labels=rownames(tmpp$rotation)[1:6],pos=2)
#text(rep(0,length(b[7:12])),b[7:12],labels=rownames(tmpp$rotation)[7:12],pos=4)

#totalk2.8=kmeans(tmpp$x[,1:2],centers=2)
par(bty='n',mai=c(0.62,0.52,0.22,0.22))
z=1
plot(-tmpp$x[totalk2.8$cluster==z,1],tmpp$x[totalk2.8$cluster==z,2],xlab="PC1",ylab="PC2",xaxt='n',yaxt='n',pch=21,bg="orange",cex=0.8,xlim=c(-2,2),ylim=c(-2,2),main=paste0("Cluster-",z,": Low"))
axis(1,at=seq(-3,2,1),labels=seq(-3,2,1),pos=-2)
axis(2,at=seq(-3,2,1),labels=seq(-3,2,1),las=2,pos=-2)
z=2
plot(-tmpp$x[totalk2.8$cluster==z,1],tmpp$x[totalk2.8$cluster==z,2],xlab="PC1",ylab="PC2",xaxt='n',yaxt='n',pch=21,bg="magenta",cex=0.8,xlim=c(-2,2),ylim=c(-2,2),main=paste0("Cluster-",z,": High"))
axis(1,at=seq(-3,2,1),labels=seq(-3,2,1),pos=-2)
axis(2,at=seq(-3,2,1),labels=seq(-3,2,1),las=2,pos=-2)

#tt=totalk2.8$cluster[match(tcell.info.acute$event_id,rownames(tmpp$x))]
#tcell.info.acute$responder_c2.8=NA
#tcell.info.acute$responder_c2.8[which(tt==2)]="high"
#tcell.info.acute$responder_c2.8[which(tt==1)]="low"
#tcell.info.acute$responder_c2.8=factor(tcell.info.acute$responder_c2.8)

plot(density(log(tcell.info.acute$event_date[which(tcell.info.acute$responder_c2.8=="low")])),xlim=c(0,4.5),main="",xaxt='n',col="orange",xlab="Days post-admission",las=2)
polygon(density(log(tcell.info.acute$event_date[which(tcell.info.acute$responder_c2.8=="low")])),col=adjustcolor("orange",alpha.f=0.5))
lines(density(log(tcell.info.acute$event_date[which(tcell.info.acute$responder_c2.8=="high")])),col="magenta")
polygon(density(log(tcell.info.acute$event_date[which(tcell.info.acute$responder_c2.8=="high")])),col=adjustcolor("magenta",alpha.f=0.5))
axis(1,at=log(c(1,2,5,10,20,50)),labels=c(1,2,5,10,20,50),pos=0)
m=summary(glm(tcell.info.acute$responder_c2.8~log(tcell.info.acute$event_date),family="binomial"))$coefficients
legend("right",legend=paste0("P = ",round(m[2,4],3)),bty='n')
legend("topright",legend=c("Low","High"),pch=22,pt.bg=c("orange","magenta"),pt.cex=1.5,bty='n')

#RBD/vrld
j=js[1]
for(k in virvars[c(3,1)])
{
  temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tcell.final.acute[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin)) #get data
  temp$var=tcell.info.acute[match(temp$event_id,tcell.info.acute$event_id),k]	#get variable
  temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
  temp$livecd3=livecd3.final[[i]][match(temp$event_id,names(livecd3.final[[i]]))] #get stim live count
  temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
  temp$gate=as.numeric(temp$gate)
  temp$days=as.numeric(temp$days)
  temp$days_symptom=as.numeric(temp$days_symptom)
  
  temp$gate2=log2(temp$gate+0.001)  #make log-transformed variable
  temp$gate3=log2(temp$gate)  #as above but don't add minimum value
  temp$gate3[is.infinite(temp$gate3)]=NA
  
  temp$responder=temp$gate>0 #make responder
  
  temp$days_binned=tcell.info.acute$event_date_binned[match(temp$event_id,tcell.info.acute$event_id)]
  temp$days_symptom_binned=tcell.info.acute$event_date_symptom_binned[match(temp$event_id,tcell.info.acute$event_id)]
  
  temp$logdays=log(temp$days) #needed for lme
  temp$logdays_symptom=log(temp$days_symptom) #needed for lme
  temp=temp[!is.na(temp$var),]	#remove missing
  temp$responder_c2=factor(tcell.info.acute$responder_c2.8[match(temp$event_id,tcell.info.acute$event_id)],c("low","high"))
  
  lgth=length(levels(temp$responder_c2))
  tmpcols=c("orange","magenta")
  temp$col="orange"
  temp$col[temp$responder_c2=="high"]="magenta"
  
  m=gam(var~responder_c2+s(logdays,bs="cr",by=responder_c2),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$responder_c2),],random =list(id=~1),family="gaussian")
  a1=anova(m)
  aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
  aa2=sum(a1$s.table[-c(1), 2])
  
  pints=as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))
  newpvs=paste0(c("Slope","High"),"-P=",sprintf("%.1e",c(pchisq(aa1,aa2, lower.tail = F),pints)))
  plot(temp$logdays,temp$var,xlab="Days post-admission",pch=21,bg=adjustcolor("orange",alpha.f = 0.5),ylab=virnams[k],main="",type='n',xaxt='n',las=2)
  axis(1,at=log(c(1,2,3,5,10,20,30,50,100)),labels=c(1,2,3,5,10,20,30,50,100))
  points(temp$logdays,temp$var,bg=adjustcolor(temp$col,alpha.f=0.5),pch=21,cex=1)
  
  for(z in 1:lgth)
  {
    tryCatch({
      sp.cis=spline.cis(cbind(temp$logdays[which(!is.na(temp$logdays)&!is.na(temp$var)&temp$responder_c2==levels(temp$responder_c2)[z])],temp$var[which(!is.na(temp$logdays)&!is.na(temp$var)&temp$responder_c2==levels(temp$responder_c2)[z])]),B=100,alpha=0.05)
      lines(x=sp.cis$x,y=sp.cis$main.curve,col=tmpcols[z],lwd=2)
      lines(x=sp.cis$x,y=sp.cis$upper.ci,col=tmpcols[z],lty=3,lwd=1)
      lines(x=sp.cis$x,y=sp.cis$lower.ci,col=tmpcols[z],lty=3,lwd=1)
    },error=function(e){})
  }
  if(k==virvars[3])
  {
    legend("bottomright",legend=c(newpvs,"Low"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",c("magenta","orange")),pt.cex=1.5,title="Responder")		
  } else {
    legend("right",legend=c(newpvs,"Low"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",c("magenta","orange")),pt.cex=1.5,title="Responder")
  }
}		  

#TG/d28 mortality
for(k in virvars[c(12,16)])
{
  temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tcell.final.acute[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin)) #get data
  temp$var=tcell.info.acute[match(temp$event_id,tcell.info.acute$event_id),k]	#get variable
  temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
  temp$livecd3=livecd3.final[[i]][match(temp$event_id,names(livecd3.final[[i]]))] #get stim live count
  temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
  temp$gate=as.numeric(temp$gate)
  temp$days=as.numeric(temp$days)
  temp$days_symptom=as.numeric(temp$days_symptom)
  
  temp$gate2=log2(temp$gate+0.001)  #make log-transformed variable
  temp$gate3=log2(temp$gate)  #as above but don't add minimum value
  temp$gate3[is.infinite(temp$gate3)]=NA
  
  temp$responder=temp$gate>0 #make responder
  
  temp$days_binned=tcell.info.acute$event_date_binned[match(temp$event_id,tcell.info.acute$event_id)]
  temp$days_symptom_binned=tcell.info.acute$event_date_symptom_binned[match(temp$event_id,tcell.info.acute$event_id)]
  
  temp$logdays=log(temp$days) #needed for lme
  temp$logdays_symptom=log(temp$days_symptom) #needed for lme
  temp=temp[!is.na(temp$var),]	#remove missing
  temp$responder_c2=factor(tcell.info.acute$responder_c2.8[match(temp$event_id,tcell.info.acute$event_id)],c("low","high"))
  
  temp$var=factor(temp$var)
  lgth=length(levels(temp$var))
  temp$col=vircols[[k]][match(temp$var,levels(temp$var))]
  
  tmp=table(temp$responder_c2,temp$var)
  tttt=colSums(tmp)
  tmp2=tmp[2,]/colSums(tmp)*100
  m=glmer(responder_c2~var+(1|id),temp,family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),nAGQ=20)
  newpvs=sprintf("%.1e",summary(m)$coefficients[2,4])
  ymax=100
  tmp=barplot(tmp2,main="",ylab="High (%)",xaxt='n',col=unlist(lapply(vircols[[k]],function(x){rep(x,1)})),las=2,xlab=virnams[k],ylim=c(0,ymax))
  axis(1,tmp,c("FALSE","TRUE"),tick=F,las=1)
  legend("topright",legend=paste0("P=",newpvs),bty='n')
  text(tmp,rep(5,length(tmp)),labels=unlist(tttt))
  
  
  m=gam(responder_c2~var+s(logdays,bs="cr",by=var),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$responder_c2),],random =list(id=~1),family="binomial")
  a1=anova(m)
  aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
  aa2=sum(a1$s.table[-c(1), 2])
  
  pints=as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))
  newpvs=paste0(c("Slope","High"),"-P=",sprintf("%.1e",c(pchisq(aa1,aa2, lower.tail = F),pints)))
  tt=vector("list",lgth)
  ttt=vector("list",lgth)
  tttt=vector("list",lgth)
  ymax=100
  tmp=table("gate"=temp$responder_c2,temp$days_binned)
  tnams=colnames(tmp)
  for(z in 1:lgth)
  {
    tt[[z]]=table("gate"=temp$responder_c2[which(temp$var==levels(temp$var)[z])],temp$days_binned[which(temp$var==levels(temp$var)[z])])		
    tt[[z]]=tt[[z]][match(c("low","high"),rownames(tt[[z]])),]
    rownames(tt[[z]])=c("low","high")
    tt[[z]][is.na(tt[[z]])]=0
    ttt[[z]]=tt[[z]][2,]/colSums(tt[[z]])*100
    tttt[[z]]=colSums(tt[[z]])
  }
  tmp=barplot(unlist(ttt),main="",ylab="High (%)",xaxt='n',col=unlist(lapply(vircols[[k]],function(x){rep(x,3)})),las=2,xlab="Days post-admission",ylim=c(0,ymax))
  axis(1,tmp,rep(tnams,lgth),tick=F,las=1)
  legend("topright",legend=c(newpvs,"TRUE","FALSE"),bty='n',pch=c(NA,rep(22,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=virnams[k])
  text(tmp,rep(5,length(tmp)),labels=unlist(tttt))
}

dev.off()

##Figure 6
svg(paste0("output/cd4cd8_new_fig6.svg"),width=10.4,height=14.9) 
layout(matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,1,1,1,1,2,2,2,2,4,4,4,4,5,5,5,5,8,8,8,8,11,11,11,11,5,5,5,5,8,8,8,8,11,11,11,11,6,6,6,6,9,9,9,9,12,12,12,12,6,6,6,6,9,9,9,9,12,12,12,12,7,7,7,7,10,10,10,10,13,13,13,13,7,7,7,7,10,10,10,10,13,13,13,13),nrow=8,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))

plot.new()
plot.new()
plot.new()
plot.new()

k=virvars[1]
l=js[5]

for(j in js[c(8,10,12)])
{
  i=is[3]
  temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tcell.final.acute[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin)) #get data
  temp$var=tcell.info.acute[match(temp$event_id,tcell.info.acute$event_id),k]	#get variable
  temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
  temp$livecd3=livecd3.final[[i]][match(temp$event_id,names(livecd3.final[[i]]))] #get stim live count
  temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
  temp$gate=as.numeric(temp$gate)
  temp$days=as.numeric(temp$days)
  temp$days_symptom=as.numeric(temp$days_symptom)
  
  temp$gate2=log2(temp$gate+0.001)  #make log-transformed variable
  temp$gate3=log2(temp$gate)  #as above but don't add minimum value
  temp$gate3[is.infinite(temp$gate3)]=NA
  
  temp$responder=as.factor(temp$gate>0) #make responder
  
  temp$days_binned=tcell.info.acute$event_date_binned[match(temp$event_id,tcell.info.acute$event_id)]
  temp$days_symptom_binned=tcell.info.acute$event_date_symptom_binned[match(temp$event_id,tcell.info.acute$event_id)]
  
  temp$logdays=log(temp$days) #needed for lme
  temp$logdays_symptom=log(temp$days_symptom) #needed for lme
  temp=temp[!is.na(temp$var),]	#remove missing
  
  lgth=length(levels(temp$responder))
  
  temp$aim=(tcell.final.acute[[2]][match(temp$event_id,tcell.info.acute$event_id),l]>0)+0
  temp$newvar=factor(paste0(temp$responder,"_",temp$aim))
  temp=temp[grep("NA",temp$newvar,invert=T),]
  lgth=length(levels(temp$newvar))
  temp$col=newnewcols[match(temp$newvar,levels(temp$newvar))]
  
  
  m=gam(var~newvar+s(logdays,bs="cr",by=newvar),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$newvar),],random =list(id=~1),family="gaussian")
  a1=anova(m)
  aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
  aa2=sum(a1$s.table[-c(1), 2])
  
  pints=as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))
  newpvs=paste0(c("Slope",paste0(c("R-","NR+","R+"),j)),"-P=",sprintf("%.1e",c(pchisq(aa1,aa2, lower.tail = F),pints)))
  plot(temp$logdays,temp$var,xlab="Days post-admission",pch=21,bg=adjustcolor("orange",alpha.f = 0.5),ylab=virnams[k],main=l,type='n',xaxt='n',las=2)
  axis(1,at=log(c(1,2,3,5,10,20,30,50,100)),labels=c(1,2,3,5,10,20,30,50,100))
  points(temp$logdays,temp$var,bg=adjustcolor(temp$col,alpha.f=0.5),pch=21,cex=1)
  
  for(z in 1:lgth)
  {
    tryCatch({
      sp.cis=spline.cis(cbind(temp$logdays[which(!is.na(temp$logdays)&!is.na(temp$var)&temp$newvar==levels(temp$newvar)[z])],temp$var[which(!is.na(temp$logdays)&!is.na(temp$var)&temp$newvar==levels(temp$newvar)[z])]),B=100,alpha=0.05)
      lines(x=sp.cis$x,y=sp.cis$main.curve,col=newnewcols[z],lwd=2)
      #lines(x=sp.cis$x,y=sp.cis$upper.ci,col=newnewcols[z],lty=3,lwd=1)
      #lines(x=sp.cis$x,y=sp.cis$lower.ci,col=newnewcols[z],lty=3,lwd=1)
    },error=function(e){})
  }
  if(k==virvars[3])
  {
    #		legend("bottomright",legend=c(newpvs,"Non-responder"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=j)		
    legend(2.2,1,legend=c(newpvs,paste0("NR-",j)),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",newnewcols[c(2:4,1)]),pt.cex=1.5,cex=0.75)		
  } else {
    #		legend("topright",legend=c(newpvs,"Non-responder"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=j)
    legend(2.25,4.75,legend=c(newpvs,paste0("NR-",j)),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",newnewcols[c(2:4,1)]),pt.cex=1.5,cex=0.75)
  }	
  
}

k=virvars[12]
l=js[1]

for(j in js[c(8,10,12)])
{
  i=is[3]
  temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tcell.final.acute[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin)) #get data
  temp$var=tcell.info.acute[match(temp$event_id,tcell.info.acute$event_id),k]	#get variable
  temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
  temp$livecd3=livecd3.final[[i]][match(temp$event_id,names(livecd3.final[[i]]))] #get stim live count
  temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
  temp$gate=as.numeric(temp$gate)
  temp$days=as.numeric(temp$days)
  temp$days_symptom=as.numeric(temp$days_symptom)
  
  temp$gate3=log2(temp$gate)  #as above but don't add minimum value
  temp$gate3[is.infinite(temp$gate3)]=NA
  
  temp$logdays=log(temp$days) #needed for lme
  temp$logdays_symptom=log(temp$days_symptom)
  temp=temp[!is.na(temp$var)&!is.na(temp$gate),]	#remove missing
  
  temp$var=factor(temp$var)
  lgth=length(levels(temp$var))
  temp$col=vircols[[k]][match(temp$var,levels(temp$var))]
  
  temp$aim=(tcell.final.acute[[2]][match(temp$event_id,tcell.info.acute$event_id),l]>0)+0
  temp$newvar=factor(paste0(temp$var,"_",temp$aim))
  temp$newvar=relevel(temp$newvar,"FALSE_1")
  temp=temp[grep("NA",temp$newvar,invert=T),]
  
  lgth=length(levels(temp$newvar))
  temp$col=newnewcols[match(temp$newvar,levels(temp$newvar))]
  m=gam(gate3~newvar+s(logdays,bs="cr",by=newvar),data=temp[!is.na(temp$logdays)&!is.na(temp$newvar)&!is.na(temp$gate3),],random =list(id=~1))
  a1=anova(m)
  aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
  aa2=sum(a1$s.table[-c(1), 2])
  
  pints=sprintf("%.1e",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
  newpvs=paste0(c("Slope",c("TG123+NR","TG45+NR","TG45+R")),"-P=",c(sprintf("%.1e",pchisq(aa1,aa2, lower.tail = F)),pints))		
  
  plot(temp$logdays,temp$gate3,xlab="Days post-admission",yaxt='n',xaxt='n',pch=21,bg=adjustcolor("azure4",alpha.f = 0.5),ylab=paste0(j," (%)"),type='n',main=l)
  points(temp$logdays,temp$gate3,bg=temp$col,pch=21,cex=1)
  axis(2,at=log2(c(0.001,0.01,0.10,1.00,10.00,100.00)),labels=c(0.001,0.01,0.1,1,10,100),las=2)
  if(j%in%js[c(10,12)])
  {
    legend("topright",legend=c(newpvs,"TG123+R"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",newnewcols[c(2:4,1)]),pt.cex=1.5,cex=0.75)		
  } else {
    legend("bottomright",legend=c(newpvs,"TG123+R"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",newnewcols[c(2:4,1)]),pt.cex=1.5,cex=0.75)
  }
  for(z in 1:lgth)
  {
    tryCatch({
      sp.cis=spline.cis(cbind(temp$logdays[which(!is.na(temp$logdays)&!is.na(temp$gate3)&temp$newvar==levels(temp$newvar)[z])],temp$gate3[which(!is.na(temp$logdays)&!is.na(temp$gate3)&temp$newvar==levels(temp$newvar)[z])]),B=100,alpha=0.05)
      lines(x=sp.cis$x,y=sp.cis$main.curve,col=newnewcols[z],lwd=2)
      #lines(x=sp.cis$x,y=sp.cis$upper.ci,col=newnewcols[z],lty=3,lwd=2)
      #lines(x=sp.cis$x,y=sp.cis$lower.ci,col=newnewcols[z],lty=3,lwd=2)
    },error=function(e){})
  }
  axis(1,at=log(c(1,2,3,5,10,20,30)),labels=c(1,2,3,5,10,20,30))  
}


k=virvars[12]
j=js[1]

for(l in js[c(8,10,12)])
{
  i=is[2]
  temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tcell.final.acute[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin)) #get data
  temp$var=tcell.info.acute[match(temp$event_id,tcell.info.acute$event_id),k]#get variable
  temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
  temp$livecd3=livecd3.final[[i]][match(temp$event_id,names(livecd3.final[[i]]))] #get stim live count
  temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
  temp$gate=as.numeric(temp$gate)
  temp$days=as.numeric(temp$days)
  temp$days_symptom=as.numeric(temp$days_symptom)
  temp$gate2=log2(temp$gate+0.001)  #make log-transformed variable
  temp$gate3=log2(temp$gate)  #as above but don't add minimum value
  temp$gate3[is.infinite(temp$gate3)]=NA
  temp$responder=as.factor(temp$gate>0) #make responder
  temp$days_binned=tcell.info.acute$event_date_binned[match(temp$event_id,tcell.info.acute$event_id)]
  temp$days_symptom_binned=tcell.info.acute$event_date_symptom_binned[match(temp$event_id,tcell.info.acute$event_id)]
  temp$logdays=log(temp$days) #needed for lme
  temp$logdays_symptom=log(temp$days_symptom) #needed for lme
  temp=temp[!is.na(temp$var),]#remove missing
  temp$var=factor(temp$var)
  lgth=length(levels(temp$var))
  temp$col=vircols[[k]][match(temp$var,levels(temp$var))]
  temp$ics=tcell.final.acute[[3]][match(temp$event_id,tcell.info.acute$event_id),l]
  
  temp=temp[which(temp$ics>0),]
  temp$ics=temp$ics>median(temp$ics)
  temp$newvar=factor(paste0(temp$responder,"_",temp$ics))
  tmp=table(temp$var,temp$newvar)
  tttt=colSums(tmp)
  tmp2=tmp[2,]/colSums(tmp)*100
  m=glm(var~relevel(temp$newvar,"TRUE_FALSE"),temp,family="binomial")
  oldnewpvs=sprintf("%.1e",summary(m)$coefficients[-1,4])
  ymax=100
  tmp=barplot(tmp2,main="",ylab="TG45 (%)",xaxt='n',col=newnewcols[c(1,3,2,4)],las=2,xlab=paste0(j," & ",l),ylim=c(0,ymax))
  axis(1,tmp,c("NR+ICSlow","NR+ICShigh","R+ICSlow","R+ICShigh"),tick=F,las=1,cex.axis=0.8)
  text(tmp,rep(5,length(tmp)),labels=unlist(tttt),col="white")
  text(tmp[-3],rep(70,length(tmp[-3])),labels=paste0("P=",oldnewpvs))
}

dev.off()

##Figure 2
svg(paste0("output/cd4cd8_new_fig2.svg"),width=10.4,height=14.9) 
layout(matrix(c(1,1,1,1,5,5,5,5,9,9,9,9,1,1,1,1,5,5,5,5,9,9,9,9,2,2,2,2,6,6,6,6,10,10,10,10,2,2,2,2,6,6,6,6,10,10,10,10,3,3,3,3,7,7,7,7,11,11,11,11,3,3,3,3,7,7,7,7,11,11,11,11,4,4,4,4,8,8,8,8,12,12,12,12,4,4,4,4,8,8,8,8,12,12,12,12),nrow=8,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))

#i=is[2]

#AIM responder vs rbd, vrld, TG, d28 mortality
for(k in virvars[c(3,1)])
{
  for(j in js[c(5,2,4,6)])
  {
    
    if(j==js[6])
    {
      i=is[3]
    } else {
      i=is[2]
    }
    temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tcell.final.acute[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin)) #get data
    temp$var=tcell.info.acute[match(temp$event_id,tcell.info.acute$event_id),k]	#get variable
    temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
    temp$livecd3=livecd3.final[[i]][match(temp$event_id,names(livecd3.final[[i]]))] #get stim live count
    temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
    temp$gate=as.numeric(temp$gate)
    temp$days=as.numeric(temp$days)
    temp$days_symptom=as.numeric(temp$days_symptom)
    
    temp$gate2=log2(temp$gate+0.001)  #make log-transformed variable
    temp$gate3=log2(temp$gate)  #as above but don't add minimum value
    temp$gate3[is.infinite(temp$gate3)]=NA
    
    temp$responder=as.factor(temp$gate>0) #make responder
    
    temp$days_binned=tcell.info.acute$event_date_binned[match(temp$event_id,tcell.info.acute$event_id)]
    temp$days_symptom_binned=tcell.info.acute$event_date_symptom_binned[match(temp$event_id,tcell.info.acute$event_id)]
    
    temp$logdays=log(temp$days) #needed for lme
    temp$logdays_symptom=log(temp$days_symptom) #needed for lme
    temp=temp[!is.na(temp$var),]	#remove missing
    
    lgth=length(levels(temp$responder))
    tmpcols=vircols[[k]]
    temp$col=vircols[[k]][1]
    temp$col[temp$responder==T]=vircols[[k]][2]
    
    m=gam(var~responder+s(logdays,bs="cr",by=responder),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$responder),],random =list(id=~1),family="gaussian")
    a1=anova(m)
    aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
    aa2=sum(a1$s.table[-c(1), 2])
    
    pints=as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))
    newpvs=paste0(c("Slope","Responder"),"-P=",sprintf("%.1e",c(pchisq(aa1,aa2, lower.tail = F),pints)))
    plot(temp$logdays,temp$var,xlab="Days post-admission",pch=21,bg=adjustcolor("orange",alpha.f = 0.5),ylab=virnams[k],main="",type='n',xaxt='n',las=2)
    axis(1,at=log(c(1,2,3,5,10,20,30,50,100)),labels=c(1,2,3,5,10,20,30,50,100))
    points(temp$logdays,temp$var,bg=adjustcolor(temp$col,alpha.f=0.5),pch=21,cex=1)
    
    for(z in 1:lgth)
    {
      tryCatch({
        sp.cis=spline.cis(cbind(temp$logdays[which(!is.na(temp$logdays)&!is.na(temp$var)&temp$responder==levels(temp$responder)[z])],temp$var[which(!is.na(temp$logdays)&!is.na(temp$var)&temp$responder==levels(temp$responder)[z])]),B=100,alpha=0.05)
        lines(x=sp.cis$x,y=sp.cis$main.curve,col=tmpcols[z],lwd=2)
        lines(x=sp.cis$x,y=sp.cis$upper.ci,col=tmpcols[z],lty=3,lwd=1)
        lines(x=sp.cis$x,y=sp.cis$lower.ci,col=tmpcols[z],lty=3,lwd=1)
      },error=function(e){})
    }
    if(k==virvars[3])
    {
      #		legend("bottomright",legend=c(newpvs,"Non-responder"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=j)		
      legend(1.9,1,legend=c(newpvs,"Non-responder"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=j)		
    } else {
      #		legend("topright",legend=c(newpvs,"Non-responder"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=j)
      legend(1.8,4.75,legend=c(newpvs,"Non-responder"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=j)
    }
  }
}

dev.off()	

##Figure 3
svg(paste0("output/cd4cd8_new_fig3.svg"),width=10.4,height=14.9) 
layout(matrix(c(1,1,1,1,5,5,5,5,9,9,9,9,1,1,1,1,5,5,5,5,9,9,9,9,2,2,2,2,6,6,6,6,10,10,10,10,2,2,2,2,6,6,6,6,10,10,10,10,3,3,3,3,7,7,7,7,11,11,11,11,3,3,3,3,7,7,7,7,11,11,11,11,4,4,4,4,8,8,8,8,12,12,12,12,4,4,4,4,8,8,8,8,13,13,13,13),nrow=8,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))

i=is[1]

for(k in virvars[c(12,16)])
{
  for(j in js[c(2,1,3,6)])
  {
    
    if(j==js[6])
    {
      i=is[3]
    } else {
      i=is[2]
    }
    
    temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tcell.final.acute[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin)) #get data
    temp$var=tcell.info.acute[match(temp$event_id,tcell.info.acute$event_id),k]	#get variable
    temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
    temp$livecd3=livecd3.final[[i]][match(temp$event_id,names(livecd3.final[[i]]))] #get stim live count
    temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
    temp$gate=as.numeric(temp$gate)
    temp$days=as.numeric(temp$days)
    temp$days_symptom=as.numeric(temp$days_symptom)
    
    temp$gate2=log2(temp$gate+0.001)  #make log-transformed variable
    temp$gate3=log2(temp$gate)  #as above but don't add minimum value
    temp$gate3[is.infinite(temp$gate3)]=NA
    
    temp$responder=as.factor(temp$gate>0) #make responder
    
    temp$days_binned=tcell.info.acute$event_date_binned[match(temp$event_id,tcell.info.acute$event_id)]
    temp$days_symptom_binned=tcell.info.acute$event_date_symptom_binned[match(temp$event_id,tcell.info.acute$event_id)]
    
    temp$logdays=log(temp$days) #needed for lme
    temp$logdays_symptom=log(temp$days_symptom) #needed for lme
    temp=temp[!is.na(temp$var),]	#remove missing
    
    temp$var=factor(temp$var)
    lgth=length(levels(temp$var))
    temp$col=vircols[[k]][match(temp$var,levels(temp$var))]
    
    tmp=table(temp$responder,temp$var)
    tttt=colSums(tmp)
    tmp2=tmp[2,]/colSums(tmp)*100
    m=glmer(responder~var+(1|id),temp,family="binomial")
    oldnewpvs=sprintf("%.1e",summary(m)$coefficients[2,4])
    ymax=100
    #tmp=barplot(tmp2,main="",ylab="Responder (%)",xaxt='n',col=unlist(lapply(vircols[[k]],function(x){rep(x,1)})),las=2,xlab=virnams[k],ylim=c(0,ymax))
    #axis(1,tmp,c("FALSE","TRUE"),tick=F,las=1)
    #legend("topright",legend=paste0("P=",newpvs),bty='n')
    #text(tmp,rep(5,length(tmp)),labels=unlist(tttt))
    
    m=gam(responder~var+s(logdays,bs="cr",by=var),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$responder),],random =list(id=~1),family="binomial")
    a1=anova(m)
    aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
    aa2=sum(a1$s.table[-c(1), 2])
    
    pints=as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))
    newpvs=paste0(c("Slope"),"-P=",sprintf("%.1e",pchisq(aa1,aa2, lower.tail = F)))
    tt=vector("list",lgth)
    ttt=vector("list",lgth)
    tttt=vector("list",lgth)
    ymax=100
    tmp=table("gate"=temp$responder,temp$days_binned)
    tnams=colnames(tmp)
    for(z in 1:lgth)
    {
      tt[[z]]=table("gate"=temp$responder[which(temp$var==levels(temp$var)[z])],temp$days_binned[which(temp$var==levels(temp$var)[z])])		
      tt[[z]]=tt[[z]][match(c("FALSE","TRUE"),rownames(tt[[z]])),]
      rownames(tt[[z]])=c("FALSE","TRUE")
      tt[[z]][is.na(tt[[z]])]=0
      ttt[[z]]=tt[[z]][2,]/colSums(tt[[z]])*100
      tttt[[z]]=colSums(tt[[z]])
    }
    tmp=barplot(unlist(ttt),main="",ylab=paste0(j," Responder (%)"),xaxt='n',col=unlist(lapply(vircols[[k]],function(x){rep(x,3)})),las=2,xlab="Days post-admission",ylim=c(0,ymax))
    axis(1,tmp,rep(tnams,lgth),tick=F,las=1)
    legend("topright",legend=c(newpvs,paste0("TRUE-P=",oldnewpvs),"FALSE"),bty='n',pch=c(NA,rep(22,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=virnams[k])
    text(tmp,rep(5,length(tmp)),labels=unlist(tttt))
  }
}


dev.off()

##Figure 5
svg(paste0("output/cd4cd8_new_fig5.svg"),width=10.4,height=14.9) 
#layout(matrix(c(1,1,1,1,2,2,3,3,4,4,5,5,1,1,1,1,2,2,3,3,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8,7,7,7,7,7,7,8,8,8,8,8,8,9,9,10,10,10,10,11,11,12,12,12,12,9,9,10,10,10,10,11,11,12,12,12,12,13,13,13,13,13,13,14,14,14,15,15,15,13,13,13,13,13,13,16,16,16,17,17,17),nrow=8,byrow=T))
layout(matrix(c(1,1,1,1,5,5,5,5,9,9,9,9,1,1,1,1,5,5,5,5,9,9,9,9,2,2,2,2,6,6,6,6,10,10,10,10,2,2,2,2,6,6,6,6,10,10,10,10,3,3,3,3,7,7,7,7,11,11,11,11,3,3,3,3,7,7,7,7,11,11,11,11,4,4,4,4,8,8,8,8,12,12,12,12,4,4,4,4,8,8,8,8,13,13,13,13),nrow=8,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))

i=is[2]

#responder vs vrld
k=virvars[1]
for(j in js[c(9,10,8,12)])
{
  
  if(j==js[9])
  {
    i=is[2]
  } else {
    i=is[3]
  }
  
  temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tcell.final.acute[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin)) #get data
  temp$var=tcell.info.acute[match(temp$event_id,tcell.info.acute$event_id),k]	#get variable
  temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
  temp$livecd3=livecd3.final[[i]][match(temp$event_id,names(livecd3.final[[i]]))] #get stim live count
  temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
  temp$gate=as.numeric(temp$gate)
  temp$days=as.numeric(temp$days)
  temp$days_symptom=as.numeric(temp$days_symptom)
  
  temp$gate2=log2(temp$gate+0.001)  #make log-transformed variable
  temp$gate3=log2(temp$gate)  #as above but don't add minimum value
  temp$gate3[is.infinite(temp$gate3)]=NA
  
  temp$responder=as.factor(temp$gate>0) #make responder
  
  temp$days_binned=tcell.info.acute$event_date_binned[match(temp$event_id,tcell.info.acute$event_id)]
  temp$days_symptom_binned=tcell.info.acute$event_date_symptom_binned[match(temp$event_id,tcell.info.acute$event_id)]
  
  temp$logdays=log(temp$days) #needed for lme
  temp$logdays_symptom=log(temp$days_symptom) #needed for lme
  temp=temp[!is.na(temp$var),]	#remove missing
  
  lgth=length(levels(temp$responder))
  tmpcols=vircols[[k]]
  temp$col=vircols[[k]][1]
  temp$col[temp$responder==T]=vircols[[k]][2]
  
  m=gam(var~responder+s(logdays,bs="cr",by=responder),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$responder),],random =list(id=~1),family="gaussian")
  a1=anova(m)
  aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
  aa2=sum(a1$s.table[-c(1), 2])
  
  pints=as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))
  newpvs=paste0(c("Slope","Responder"),"-P=",sprintf("%.1e",c(pchisq(aa1,aa2, lower.tail = F),pints)))
  plot(temp$logdays,temp$var,xlab="Days post-admission",pch=21,bg=adjustcolor("orange",alpha.f = 0.5),ylab=virnams[k],main="",type='n',xaxt='n',las=2)
  axis(1,at=log(c(1,2,3,5,10,20,30,50,100)),labels=c(1,2,3,5,10,20,30,50,100))
  points(temp$logdays,temp$var,bg=adjustcolor(temp$col,alpha.f=0.5),pch=21,cex=1)
  
  for(z in 1:lgth)
  {
    tryCatch({
      sp.cis=spline.cis(cbind(temp$logdays[which(!is.na(temp$logdays)&!is.na(temp$var)&temp$responder==levels(temp$responder)[z])],temp$var[which(!is.na(temp$logdays)&!is.na(temp$var)&temp$responder==levels(temp$responder)[z])]),B=100,alpha=0.05)
      lines(x=sp.cis$x,y=sp.cis$main.curve,col=tmpcols[z],lwd=2)
      lines(x=sp.cis$x,y=sp.cis$upper.ci,col=tmpcols[z],lty=3,lwd=1)
      lines(x=sp.cis$x,y=sp.cis$lower.ci,col=tmpcols[z],lty=3,lwd=1)
    },error=function(e){})
  }
  if(k==virvars[3])
  {
    #		legend("bottomright",legend=c(newpvs,"Non-responder"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=j)		
    legend(1.9,1,legend=c(newpvs,"Non-responder"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=j)		
  } else {
    #		legend("topright",legend=c(newpvs,"Non-responder"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=j)
    legend(1.9,4.75,legend=c(newpvs,"Non-responder"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=j)
  }
}	

#vs TG
k=virvars[12]
for(j in js[c(9,10,8,12)])
{
  
  if(j==js[9])
  {
    i=is[2]
  } else {
    i=is[3]
  }
  
  temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tcell.final.acute[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin)) #get data
  temp$var=tcell.info.acute[match(temp$event_id,tcell.info.acute$event_id),k]	#get variable
  temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
  temp$livecd3=livecd3.final[[i]][match(temp$event_id,names(livecd3.final[[i]]))] #get stim live count
  temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
  temp$gate=as.numeric(temp$gate)
  temp$days=as.numeric(temp$days)
  temp$days_symptom=as.numeric(temp$days_symptom)
  
  temp$gate3=log2(temp$gate)  #as above but don't add minimum value
  temp$gate3[is.infinite(temp$gate3)]=NA
  
  temp$logdays=log(temp$days) #needed for lme
  temp$logdays_symptom=log(temp$days_symptom)
  temp=temp[!is.na(temp$var)&!is.na(temp$gate),]	#remove missing
  
  temp$var=factor(temp$var)
  lgth=length(levels(temp$var))
  temp$col=vircols[[k]][match(temp$var,levels(temp$var))]
  
  m=gam(gate3~var+s(logdays,bs="cr",by=var),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$gate3),],random =list(id=~1))
  a1=anova(m)
  aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
  aa2=sum(a1$s.table[-c(1), 2])
  
  pints=sprintf("%.1e",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
  newpvs=paste0(c("Slope",levels(temp$var)[-1]),"-P=",c(sprintf("%.1e",pchisq(aa1,aa2, lower.tail = F)),pints))		
  
  plot(temp$logdays,temp$gate3,xlab="Days post-admission",yaxt='n',xaxt='n',pch=21,bg=adjustcolor("azure4",alpha.f = 0.5),ylab=paste0(j," (%)"),type='n')
  points(temp$logdays,temp$gate3,bg=temp$col,pch=21,cex=1)
  axis(2,at=log2(c(0.001,0.01,0.10,1.00,10.00,100.00)),labels=c(0.001,0.01,0.1,1,10,100),las=2)
  if(j%in%js[c(10,12)])
  {
    legend("topright",legend=c(newpvs,"FALSE"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=virnams[k])		
  } else {
    legend("bottomright",legend=c(newpvs,"FALSE"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=virnams[k])
  }
  for(z in 1:lgth)
  {
    tryCatch({
      sp.cis=spline.cis(cbind(temp$logdays[which(!is.na(temp$logdays)&!is.na(temp$gate3)&temp$var==levels(temp$var)[z])],temp$gate3[which(!is.na(temp$logdays)&!is.na(temp$gate3)&temp$var==levels(temp$var)[z])]),B=100,alpha=0.05)
      lines(x=sp.cis$x,y=sp.cis$main.curve,col=vircols[[k]][z],lwd=2)
      lines(x=sp.cis$x,y=sp.cis$upper.ci,col=vircols[[k]][z],lty=3,lwd=2)
      lines(x=sp.cis$x,y=sp.cis$lower.ci,col=vircols[[k]][z],lty=3,lwd=2)
    },error=function(e){})
  }
  axis(1,at=log(c(1,2,3,5,10,20,30)),labels=c(1,2,3,5,10,20,30))		
}

k=virvars[16]
j=js[9]
i=is[2]

temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tcell.final.acute[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin)) #get data
temp$var=tcell.info.acute[match(temp$event_id,tcell.info.acute$event_id),k]	#get variable
temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
temp$livecd3=livecd3.final[[i]][match(temp$event_id,names(livecd3.final[[i]]))] #get stim live count
temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
temp$gate=as.numeric(temp$gate)
temp$days=as.numeric(temp$days)
temp$days_symptom=as.numeric(temp$days_symptom)

temp$gate3=log2(temp$gate)  #as above but don't add minimum value
temp$gate3[is.infinite(temp$gate3)]=NA

temp$logdays=log(temp$days) #needed for lme
temp$logdays_symptom=log(temp$days_symptom)
temp=temp[!is.na(temp$var)&!is.na(temp$gate),]	#remove missing

temp$var=factor(temp$var=="1")
lgth=length(levels(temp$var))
temp$col=vircols[[k]][match(temp$var,levels(temp$var))]

m=gam(gate3~var+s(logdays,bs="cr",by=var),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$gate3),],random =list(id=~1))
a1=anova(m)
aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
aa2=sum(a1$s.table[-c(1), 2])

pints=sprintf("%.1e",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
newpvs=paste0(c("Slope",levels(temp$var)[-1]),"-P=",c(sprintf("%.1e",pchisq(aa1,aa2, lower.tail = F)),pints))		

plot(temp$logdays,temp$gate3,xlab="Days post-admission",yaxt='n',xaxt='n',pch=21,bg=adjustcolor("azure4",alpha.f = 0.5),ylab=paste0(j," (%)"),type='n')
points(temp$logdays,temp$gate3,bg=temp$col,pch=21,cex=1)
axis(2,at=log2(c(0.001,0.01,0.10,1.00,10.00,100.00)),labels=c(0.001,0.01,0.1,1,10,100),las=2)
if(j%in%js[c(10,12)])
{
  legend("topright",legend=c(newpvs,"FALSE"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=virnams[k])		
} else {
  legend("bottomright",legend=c(newpvs,"FALSE"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=virnams[k])
}
for(z in 1:lgth)
{
  tryCatch({
    sp.cis=spline.cis(cbind(temp$logdays[which(!is.na(temp$logdays)&!is.na(temp$gate3)&temp$var==levels(temp$var)[z])],temp$gate3[which(!is.na(temp$logdays)&!is.na(temp$gate3)&temp$var==levels(temp$var)[z])]),B=100,alpha=0.05)
    lines(x=sp.cis$x,y=sp.cis$main.curve,col=vircols[[k]][z],lwd=2)
    lines(x=sp.cis$x,y=sp.cis$upper.ci,col=vircols[[k]][z],lty=3,lwd=2)
    lines(x=sp.cis$x,y=sp.cis$lower.ci,col=vircols[[k]][z],lty=3,lwd=2)
  },error=function(e){})
}
axis(1,at=log(c(1,2,3,5,10,20,30)),labels=c(1,2,3,5,10,20,30))


dev.off

#ncs.final.acute=ncs[match(rownames(tcell.final.acute[[1]]),rownames(ncs)),]

svg(paste0("output/dmso_fig1.svg"),width=10.4,height=14.9) 
layout(matrix(c(1,1,1,1,2,2,3,3,4,4,5,5,1,1,1,1,2,2,3,3,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8,7,7,7,7,7,7,8,8,8,8,8,8,9,9,10,10,10,10,11,11,12,12,12,12,9,9,10,10,10,10,11,11,12,12,12,12,13,13,13,13,13,13,14,14,14,15,15,15,13,13,13,13,13,13,16,16,16,17,17,17),nrow=8,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))

#PCA and clusters
i=is[1]
#tmp=ncs.final.acute[apply(ncs.final.acute,1,function(x){sum(is.na(x))/length(x)*100})==0,]
#tmpp=prcomp(tmp[,js],scale=T)
#tmpc=rowSums(tmp[,js[1:12]]>0)

#tmpp2=prcomp(tmp[tmpp$x[,1]<8,js],scale=T)

ii=cut(tmpc,breaks=seq(min(tmpc),max(tmpc),len=100),include.lowest=T)
tcolors=colorRampPalette(spctrl[rev(c(2:5,8:10))])(99)[ii]
plot(tmpp2$x[,1],tmpp2$x[,2],xlab=paste0("PC1 (",round(summary(tmpp2)[[6]][2,1]*100,1),"%)"),ylab=paste0("PC2 (",round(summary(tmpp2)[[6]][2,2]*100,1),"%)"),xaxt='n',yaxt='n',pch=21,bg=tcolors,cex=1)
axis(1,at=seq(-5,15,5),labels=seq(-5,15,5))
axis(2,at=seq(-15,10,5),labels=seq(-15,10,5),las=2)
abline(h=0,v=0,lty=3)
legend("bottomleft",title="No. +",legend=seq(0,12,2),pch=21,pt.bg=tcolors[match(seq(0,12,2),tmpc)],bty='n')

par(bty='n',mai=c(0.62,0.22,0.22,0.12))

b=barplot(tmpp2$rotation[,1],horiz=T,xlab="PC1 loading",yaxt='n',xaxt='n',col=c(rep("cornflowerblue",6),rep("chartreuse3",6)))
axis(1,at=seq(0,0.4,0.1),labels=seq(0,0.4,0.1))
text(rep(0,length(b)),b,labels=rownames(tmpp2$rotation),pos=4,cex=0.75)

b=barplot(tmpp2$rotation[,2],horiz=T,xlab="PC2 loading",yaxt='n',xaxt='n',col=c(rep("cornflowerblue",6),rep("chartreuse3",6)),xlim=c(-0.4,0.4))
axis(1,at=seq(-0.4,0.4,0.2),labels=seq(-0.4,0.4,0.2))
#text(rep(0,length(b[1:6])),b[1:6],labels=rownames(tmpp2$rotation)[1:6],pos=2)
#text(rep(0,length(b[7:12])),b[7:12],labels=rownames(tmpp2$rotation)[7:12],pos=4)

#ncsk2=kmeans(tmpp2$x[,1:2],centers=2)
par(bty='n',mai=c(0.62,0.52,0.22,0.22))
z=1
plot(tmpp2$x[ncsk2$cluster==z,1],tmpp2$x[ncsk2$cluster==z,2],xlab="PC1",ylab="PC2",xaxt='n',yaxt='n',pch=21,bg="orange",cex=0.8,xlim=c(-2,10),ylim=c(-15,6),main=paste0("Cluster-",z,": Low"))
axis(1,at=seq(-5,15,5),labels=seq(-5,15,5),pos=-15)
axis(2,at=seq(-15,10,5),labels=seq(-15,10,5),las=2,pos=-2)
z=2
plot(tmpp2$x[ncsk2$cluster==z,1],tmpp2$x[ncsk2$cluster==z,2],xlab="PC1",ylab="PC2",xaxt='n',yaxt='n',pch=21,bg="magenta",cex=0.8,xlim=c(-2,10),ylim=c(-15,6),main=paste0("Cluster-",z,": High"))
axis(1,at=seq(-5,15,5),labels=seq(-5,15,5),pos=-15)
axis(2,at=seq(-15,10,5),labels=seq(-15,10,5),las=2,pos=-2)

#tt=ncsk2$cluster[match(tcell.info.acute$event_id,rownames(tmpp2$x))]
#tcell.info.acute$ncs_responder=NA
#tcell.info.acute$ncs_responder[which(tt==1)]="low"
#tcell.info.acute$ncs_responder[which(tt==2)]="high"
#tcell.info.acute$ncs_responder=factor(tcell.info.acute$ncs_responder)

plot(density(log(tcell.info.acute$event_date[which(tcell.info.acute$ncs_responder=="low")])),xlim=c(0,4.5),main="",xaxt='n',col="orange",xlab="Days post-admission",las=2,ylim=c(0,0.42))
polygon(density(log(tcell.info.acute$event_date[which(tcell.info.acute$ncs_responder=="low")])),col=adjustcolor("orange",alpha.f=0.5))
lines(density(log(tcell.info.acute$event_date[which(tcell.info.acute$ncs_responder=="high")])),col="magenta")
polygon(density(log(tcell.info.acute$event_date[which(tcell.info.acute$ncs_responder=="high")])),col=adjustcolor("magenta",alpha.f=0.5))
axis(1,at=log(c(1,2,5,10,20,50)),labels=c(1,2,5,10,20,50),pos=0)
m=summary(glm(tcell.info.acute$ncs_responder~log(tcell.info.acute$event_date),family="binomial"))$coefficients
legend("right",legend=paste0("P = ",round(m[2,4],3)),bty='n')
legend("topright",legend=c("Low","High"),pch=22,pt.bg=c("orange","magenta"),pt.cex=1.5,bty='n')

#RBD/vrld
j=js[1]
for(k in virvars[c(3,1)])
{
  temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(ncs.final.acute[,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin)) #get data
  temp$var=tcell.info.acute[match(temp$event_id,tcell.info.acute$event_id),k]	#get variable
  temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
  temp$livecd3=livecd3.final[[i]][match(temp$event_id,names(livecd3.final[[i]]))] #get stim live count
  temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
  temp$gate=as.numeric(temp$gate)
  temp$days=as.numeric(temp$days)
  temp$days_symptom=as.numeric(temp$days_symptom)
  
  temp$gate2=log2(temp$gate+0.001)  #make log-transformed variable
  temp$gate3=log2(temp$gate)  #as above but don't add minimum value
  temp$gate3[is.infinite(temp$gate3)]=NA
  
  temp$responder=temp$gate>0 #make responder
  
  temp$days_binned=tcell.info.acute$event_date_binned[match(temp$event_id,tcell.info.acute$event_id)]
  temp$days_symptom_binned=tcell.info.acute$event_date_symptom_binned[match(temp$event_id,tcell.info.acute$event_id)]
  
  temp$logdays=log(temp$days) #needed for lme
  temp$logdays_symptom=log(temp$days_symptom) #needed for lme
  temp=temp[!is.na(temp$var),]	#remove missing
  temp$ncs_responder_c2=factor(tcell.info.acute$ncs_responder[match(temp$event_id,tcell.info.acute$event_id)],c("low","high"))
  
  lgth=length(levels(temp$ncs_responder_c2))
  tmpcols=c("orange","magenta")
  temp$col="orange"
  temp$col[temp$ncs_responder_c2=="high"]="magenta"
  
  m=gam(var~ncs_responder_c2+s(logdays,bs="cr",by=ncs_responder_c2),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$ncs_responder_c2),],random =list(id=~1),family="gaussian")
  a1=anova(m)
  aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
  aa2=sum(a1$s.table[-c(1), 2])
  
  pints=as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))
  newpvs=paste0(c("Slope","High"),"-P=",sprintf("%.1e",c(pchisq(aa1,aa2, lower.tail = F),pints)))
  plot(temp$logdays,temp$var,xlab="Days post-admission",pch=21,bg=adjustcolor("orange",alpha.f = 0.5),ylab=virnams[k],main="",type='n',xaxt='n',las=2)
  axis(1,at=log(c(1,2,3,5,10,20,30,50,100)),labels=c(1,2,3,5,10,20,30,50,100))
  points(temp$logdays,temp$var,bg=adjustcolor(temp$col,alpha.f=0.5),pch=21,cex=1)
  
  for(z in 1:lgth)
  {
    tryCatch({
      sp.cis=spline.cis(cbind(temp$logdays[which(!is.na(temp$logdays)&!is.na(temp$var)&temp$ncs_responder_c2==levels(temp$ncs_responder_c2)[z])],temp$var[which(!is.na(temp$logdays)&!is.na(temp$var)&temp$ncs_responder_c2==levels(temp$ncs_responder_c2)[z])]),B=100,alpha=0.05)
      lines(x=sp.cis$x,y=sp.cis$main.curve,col=tmpcols[z],lwd=2)
      lines(x=sp.cis$x,y=sp.cis$upper.ci,col=tmpcols[z],lty=3,lwd=1)
      lines(x=sp.cis$x,y=sp.cis$lower.ci,col=tmpcols[z],lty=3,lwd=1)
    },error=function(e){})
  }
  if(k==virvars[3])
  {
    legend("bottomright",legend=c(newpvs,"Low"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",c("magenta","orange")),pt.cex=1.5,title="Responder")		
  } else {
    legend("right",legend=c(newpvs,"Low"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",c("magenta","orange")),pt.cex=1.5,title="Responder")
  }
}		  

#TG/d28 mortality
for(k in virvars[c(12,16)])
{
  temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(ncs.final.acute[,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin)) #get data
  temp$var=tcell.info.acute[match(temp$event_id,tcell.info.acute$event_id),k]	#get variable
  temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
  temp$livecd3=livecd3.final[[i]][match(temp$event_id,names(livecd3.final[[i]]))] #get stim live count
  temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
  temp$gate=as.numeric(temp$gate)
  temp$days=as.numeric(temp$days)
  temp$days_symptom=as.numeric(temp$days_symptom)
  
  temp$gate2=log2(temp$gate+0.001)  #make log-transformed variable
  temp$gate3=log2(temp$gate)  #as above but don't add minimum value
  temp$gate3[is.infinite(temp$gate3)]=NA
  
  temp$responder=temp$gate>0 #make responder
  
  temp$days_binned=tcell.info.acute$event_date_binned[match(temp$event_id,tcell.info.acute$event_id)]
  temp$days_symptom_binned=tcell.info.acute$event_date_symptom_binned[match(temp$event_id,tcell.info.acute$event_id)]
  
  temp$logdays=log(temp$days) #needed for lme
  temp$logdays_symptom=log(temp$days_symptom) #needed for lme
  temp=temp[!is.na(temp$var),]	#remove missing
  temp$ncs_responder_c2=factor(tcell.info.acute$ncs_responder[match(temp$event_id,tcell.info.acute$event_id)],c("low","high"))
  
  temp$var=factor(temp$var)
  lgth=length(levels(temp$var))
  temp$col=vircols[[k]][match(temp$var,levels(temp$var))]
  
  tmp=table(temp$ncs_responder_c2,temp$var)
  tttt=colSums(tmp)
  tmp2=tmp[2,]/colSums(tmp)*100
  m=glmer(ncs_responder_c2~var+(1|id),temp,family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)),nAGQ=20)
  newpvs=sprintf("%.1e",summary(m)$coefficients[2,4])
  ymax=100
  tmp=barplot(tmp2,main="",ylab="High (%)",xaxt='n',col=unlist(lapply(vircols[[k]],function(x){rep(x,1)})),las=2,xlab=virnams[k],ylim=c(0,ymax))
  axis(1,tmp,c("FALSE","TRUE"),tick=F,las=1)
  legend("topright",legend=paste0("P=",newpvs),bty='n')
  text(tmp,rep(5,length(tmp)),labels=unlist(tttt))
  
  
  m=gam(ncs_responder_c2~var+s(logdays,bs="cr",by=var),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$ncs_responder_c2),],random =list(id=~1),family="binomial")
  a1=anova(m)
  aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
  aa2=sum(a1$s.table[-c(1), 2])
  
  pints=as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))
  newpvs=paste0(c("Slope","High"),"-P=",sprintf("%.1e",c(pchisq(aa1,aa2, lower.tail = F),pints)))
  tt=vector("list",lgth)
  ttt=vector("list",lgth)
  tttt=vector("list",lgth)
  ymax=100
  tmp=table("gate"=temp$ncs_responder_c2,temp$days_binned)
  tnams=colnames(tmp)
  for(z in 1:lgth)
  {
    tt[[z]]=table("gate"=temp$ncs_responder_c2[which(temp$var==levels(temp$var)[z])],temp$days_binned[which(temp$var==levels(temp$var)[z])])		
    tt[[z]]=tt[[z]][match(c("low","high"),rownames(tt[[z]])),]
    rownames(tt[[z]])=c("low","high")
    tt[[z]][is.na(tt[[z]])]=0
    ttt[[z]]=tt[[z]][2,]/colSums(tt[[z]])*100
    tttt[[z]]=colSums(tt[[z]])
  }
  tmp=barplot(unlist(ttt),main="",ylab="High (%)",xaxt='n',col=unlist(lapply(vircols[[k]],function(x){rep(x,3)})),las=2,xlab="Days post-admission",ylim=c(0,ymax))
  axis(1,tmp,rep(tnams,lgth),tick=F,las=1)
  legend("topright",legend=c(newpvs,"TRUE","FALSE"),bty='n',pch=c(NA,rep(22,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=virnams[k])
  text(tmp,rep(5,length(tmp)),labels=unlist(tttt))
}

dev.off()