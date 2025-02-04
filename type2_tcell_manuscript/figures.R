##Figure 1
svg(paste0("output/fig1.svg"),width=10.4,height=14.9) 
layout(matrix(c(1,1,1,1,2,2,3,3,4,4,5,5,1,1,1,1,2,2,3,3,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8,7,7,7,7,7,7,8,8,8,8,8,8,9,9,10,10,10,10,11,11,12,12,12,12,9,9,10,10,10,10,11,11,12,12,12,12,13,13,13,13,13,13,14,14,14,15,15,15,13,13,13,13,13,13,16,16,16,17,17,17),nrow=8,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))

#PCA and clusters
i=is[1]
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

#totalk2=kmeans(tmpp$x[,1:2],centers=2)
par(bty='n',mai=c(0.62,0.52,0.22,0.22))
z=1
plot(-tmpp$x[totalk2$cluster==z,1],tmpp$x[totalk2$cluster==z,2],xlab="PC1",ylab="PC2",xaxt='n',yaxt='n',pch=21,bg="orange",cex=0.8,xlim=c(-2,2),ylim=c(-1.9,1.7),main=paste0("Low-responders"))
axis(1,at=seq(-3,2,1),labels=seq(-3,2,1),pos=-2)
axis(2,at=seq(-3,2,1),labels=seq(-3,2,1),las=2,pos=-2)
z=2
plot(-tmpp$x[totalk2$cluster==z,1],tmpp$x[totalk2$cluster==z,2],xlab="PC1",ylab="PC2",xaxt='n',yaxt='n',pch=21,bg="magenta",cex=0.8,xlim=c(-2,2),ylim=c(-1.9,1.7),main=paste0("High-responders"))
axis(1,at=seq(-3,2,1),labels=seq(-3,2,1),pos=-2)
axis(2,at=seq(-3,2,1),labels=seq(-3,2,1),las=2,pos=-2)

plot(density(log(tcell.info.acute$event_date[which(tcell.info.acute$responder_c2=="low")])),xlim=c(0,4.5),main="",xaxt='n',col="orange",xlab="Days post-admission",las=2)
polygon(density(log(tcell.info.acute$event_date[which(tcell.info.acute$responder_c2=="low")])),col=adjustcolor("orange",alpha.f=0.5))
lines(density(log(tcell.info.acute$event_date[which(tcell.info.acute$responder_c2=="high")])),col="magenta")
polygon(density(log(tcell.info.acute$event_date[which(tcell.info.acute$responder_c2=="high")])),col=adjustcolor("magenta",alpha.f=0.5))
axis(1,at=log(c(1,2,5,10,20,50)),labels=c(1,2,5,10,20,50),pos=0)
m=summary(glm(tcell.info.acute$responder_c2~log(tcell.info.acute$event_date),family="binomial"))$coefficients
legend("right",legend=paste0("P = ",round(m[2,4],3)),bty='n')
legend("topright",legend=c("Low-responder","High-responder"),pch=22,pt.bg=c("orange","magenta"),pt.cex=1.5,bty='n')

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
  temp$responder_c2=factor(tcell.info.acute$responder_c2[match(temp$event_id,tcell.info.acute$event_id)],c("low","high"))
  
  lgth=length(levels(temp$responder_c2))
  tmpcols=c("orange","magenta")
  temp$col="orange"
  temp$col[temp$responder_c2=="high"]="magenta"
  
  if(as.numeric(unlist(strsplit(resp.acute.pvalues.int[[i]][[2]][1,k],"/")))<0.001)
  {
    pints=sprintf("%.1E",as.numeric(unlist(strsplit(resp.acute.pvalues.int[[i]][[2]][1,k],"/"))))
  } else {
    pints=round(as.numeric(unlist(strsplit(resp.acute.pvalues.int[[i]][[2]][1,k],"/"))),3)
  }
  #newpvs=paste0(c("Slope",levels(temp$var)[-1]),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))		
  if(resp.acute.pvalues[[i]][[2]][1,k]<0.001)
  {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",resp.acute.pvalues[[i]][[2]][1,k]),pints))				
  } else {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(resp.acute.pvalues[[i]][[2]][1,k],3),pints))		
  }
  
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
    legend("bottomright",legend=c(newpvs,"Low-responder","High-responder"),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",c("orange","magenta")),pt.cex=1.5)		
  } else {
    legend("right",legend=c(newpvs,"Low-responder","High-responder"),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",c("orange","magenta")),pt.cex=1.5)
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
  temp$responder_c2=factor(tcell.info.acute$responder_c2[match(temp$event_id,tcell.info.acute$event_id)],c("low","high"))
  
  temp$var=factor(temp$var)
  lgth=length(levels(temp$var))
  temp$col=vircols[[k]][match(temp$var,levels(temp$var))]
  
  tmp=table(temp$responder_c2,temp$var)
  tttt=colSums(tmp)
  tmp2=tmp[2,]/colSums(tmp)*100
  m=glmer(responder_c2~var+(1|id),temp,family="binomial")
  if(summary(m)$coefficients[2,4]<0.001)
  {
    newpvs=sprintf("%.1e",summary(m)$coefficients[2,4])
  } else {
    newpvs=round(summary(m)$coefficients[2,4],3)
  }
  
  ymax=100
  tmp=barplot(tmp2,main="",ylab="High-responder (%)",xaxt='n',col=unlist(lapply(vircols[[k]],function(x){rep(x,1)})),las=2,xlab=virnams[k],ylim=c(0,ymax))
  axis(1,tmp,virlvls[[k]],tick=F,las=1)
  legend("topright",legend=paste0("P=",newpvs),bty='n')
  text(tmp,rep(5,length(tmp)),labels=unlist(tttt))
  
  
  if(resp.acute.pvalues[[i]][[2]][1,k]<0.001)
  {
    newpvs=paste0(c("Slope"),"-P=",c(sprintf("%.1E",resp.acute.pvalues[[i]][[2]][1,k])))				
  } else {
    newpvs=paste0(c("Slope"),"-P=",c(round(resp.acute.pvalues[[i]][[2]][1,k],3)))		
  }
  
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
  tmp=barplot(unlist(ttt),main="",ylab="High-responder (%)",xaxt='n',col=unlist(lapply(vircols[[k]],function(x){rep(x,3)})),las=2,xlab="Days post-admission",ylim=c(0,ymax))
  axis(1,tmp,rep(tnams,lgth),tick=F,las=1)
  legend("topright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,rep(22,lgth)),pt.bg=c("white",vircols[[k]]),pt.cex=1.5)
  text(tmp,rep(5,length(tmp)),labels=unlist(tttt))
}

#omics
plot.new()

#indiv examples
for(k in omics.vars[c(2,11,38,55)])
{
  temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tcell.final.acute[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin,"responder"=tcell.info.acute$responder_c2,"PC1"=tcell.info.acute$PC1)) #get data
  if(omics.type[k]=="cytof")
  {
    tmp=cytof[gsub(".*[-]","",cytof$event_id)=="1",]
    tmpp=cytof.pct[gsub(".*[-]","",cytof$event_id)=="1",]
    temp$bvar=tmpp[match(temp$id,tmp$participant_id),k]
  } else {
    tmp=olink[gsub(".*[-]","",olink$event_id)=="1",]
    temp$bvar=olink[match(temp$id,olink$participant_id),k]
  }
  temp$PC1=as.numeric(temp$PC1)
  temp$bvar2=log10(temp$bvar+min(temp$bvar[which(temp$bvar>0)]))	
  if(omics.type[k]=="cytof")
  {
    temp$var=cytof.pct[match(temp$event_id,cytof$event_id),k]
  } else {
    temp$var=olink[match(temp$event_id,olink$event_id),k]
  }	
  temp$var2=log10(temp$var+min(temp$var[which(temp$var>0)]))
  temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
  temp$livecd3=livecd3.final[[1]][match(temp$event_id,names(livecd3.final[[1]]))] #get stim live count
  temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
  temp$responder=factor(temp$responder=="1")
  temp$logdays=log(as.numeric(temp$days))	
  
  yrng=range(c(temp$var+min(temp$var[which(temp$var>0)]),temp$bvar+min(temp$bvar[which(temp$bvar>0)])),na.rm=T)
  
  if(yrng[1]>0)
  {
    boxplot(var2~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab="",ylab="",type='n',border="white",pch="",col="white",xaxt='n',log='y',las=2)
    par(new=T)
    boxplot(bvar2~responder,temp,xlim=c(0.5,5.5),ylim=log10(yrng),xlab="Responder",ylab=omics.nams2[k],type='n',border="white",pch="",col="white",yaxt='n',xaxt='n')
    stripchart(bvar2~responder,temp[which(temp$responder==F),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("orange",alpha.f=0.5),cex=0.5,lwd=0.2)
    stripchart(bvar2~responder,temp[which(temp$responder==T),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("magenta",alpha.f=0.5),cex=0.5,lwd=0.2)
    par(new=T)
    boxplot(bvar2~responder,temp,xlim=c(0.5,5.5),ylim=log10(yrng),xlab="",ylab="",type='n',medcol="black",pch="",col=adjustcolor(c("orange","magenta"),alpha.f=0.1),xaxt='n',yaxt='n')
    par(new=T)
    boxplot(var2~responder,temp,xlim=c(0.5,5.5),ylim=log10(yrng),xlab="Responder",ylab="",type='n',border="white",pch="",col="white",yaxt='n',xaxt='n',at=4:5)
    stripchart(var2~responder,temp[which(temp$responder==F),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("orange",alpha.f=0.5),cex=0.5,at=4:5,lwd=0.2)
    stripchart(var2~responder,temp[which(temp$responder==T),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("magenta",alpha.f=0.5),cex=0.5,at=4:5,lwd=0.2)
    par(new=T)
    boxplot(var2~responder,temp,xlim=c(0.5,5.5),ylim=log10(yrng),xlab="",ylab="",type='n',medcol="black",pch="",col=adjustcolor(c("orange","magenta"),alpha.f=0.25),xaxt='n',yaxt='n',at=4:5)
    axis(1,at=c(1:2,4:5),labels=rep(c("Low","High"),2))
    
    if(omics.baseline.res.responder.log[k,1]<0.001)
    {
      newpvs1=sprintf("%.1E",omics.baseline.res.responder.log[k,1])				
    } else {
      newpvs1=round(omics.baseline.res.responder.log[k,1],3)		
    }
    if(omics.res.responder.log[k,1]<0.001)
    {
      newpvs2=sprintf("%.1E",omics.res.responder.log[k,1])				
    } else {
      newpvs2=round(omics.res.responder.log[k,1],3)		
    }		
    newpvs=c(newpvs1,newpvs2)
  } else {
    yrng=range(c(temp$var,temp$bvar),na.rm=T)
    boxplot(bvar~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab="Responder",ylab=omics.nams2[k],type='n',border="white",pch="",col="white",xaxt='n',las=2)
    stripchart(bvar~responder,temp[which(temp$responder==F),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("orange",alpha.f=0.5),cex=0.25,lwd=0.2)
    stripchart(bvar~responder,temp[which(temp$responder==T),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("magenta",alpha.f=0.5),cex=0.25,lwd=0.2)
    par(new=T)
    boxplot(bvar~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab="",ylab="",type='n',medcol="black",pch="",col=adjustcolor(c("orange","magenta"),alpha.f=0.1),xaxt='n',yaxt='n')
    par(new=T)
    boxplot(var~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab="Responder",ylab="",type='n',border="white",pch="",col="white",yaxt='n',xaxt='n',at=4:5)
    stripchart(var~responder,temp[which(temp$responder==F),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("orange",alpha.f=0.5),cex=0.25,at=4:5,lwd=0.2)
    stripchart(var~responder,temp[which(temp$responder==T),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("magenta",alpha.f=0.5),cex=0.25,at=4:5,lwd=0.2)
    par(new=T)
    boxplot(var~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab="",ylab="",type='n',medcol="black",pch="",col=adjustcolor(c("orange","magenta"),alpha.f=0.25),xaxt='n',yaxt='n',at=4:5)
    axis(1,at=c(1:2,4:5),labels=rep(c("Low","High"),2))
    newpvs=sprintf("%.1e",c(omics.baseline.res.responder.log[k,1],omics.res.responder.log[k,1]))
  }
  abline(v=3,lty=3)
  axis(3,at=c(1.5,4.5),tick=F,labels=paste0(c("Baseline-P=","P="),newpvs),line=-0.5)
  
}

dev.off()

##all omics
#pheatmap(omics.responder.aves[,apply(omics.responder.aves,2,sd)>0],display_numbers=omics.responder.sigs[,apply(omics.responder.aves,2,sd)>0],scale="column",cluster_rows=F,labels_col=omics.nams2[apply(omics.responder.aves,2,sd)>0],angle=45,labels_row=c("Baseline:Low","Baseline:High","Low","High"),fontsize_number=10,filename="output/fig1_omics.png",width=10.4,height=3.725,fontsize_col=6,fontsize_row=6,legend=F,cellwidth=10,cellheight=10)
#dev.off()
##sig omics
#pheatmap(omics.responder.aves[,apply(omics.responder.sigs,2,function(x){sum(x=="*")})>0],display_numbers=omics.responder.sigs[,apply(omics.responder.sigs,2,function(x){sum(x=="*")})>0],scale="column",cluster_rows=F,labels_col=omics.nams2[apply(omics.responder.sigs,2,function(x){sum(x=="*")})>0],angle=45,labels_row=c("Baseline:Low","Baseline:High","Low","High"),fontsize_number=20,filename="output/fig1_omics_sig.png",width=5.2,height=3.725,fontsize_col=8,fontsize_row=8,legend=F,cellwidth=20,cellheight=20)
#pheatmap(omics.responder.aves[,apply(omics.responder.sigs,2,function(x){sum(x=="*")})>0],display_numbers=omics.responder.sigs[,apply(omics.responder.sigs,2,function(x){sum(x=="*")})>0],scale="column",cluster_rows=F,labels_col=omics.nams2[apply(omics.responder.sigs,2,function(x){sum(x=="*")})>0],angle=45,labels_row=c("Baseline:Low","Baseline:High","Low","High"),fontsize_number=20,filename="output/fig1_omics_sig_legend.png",width=5.2,height=3.725,fontsize_col=8,fontsize_row=8,legend=T,cellwidth=20,cellheight=20)
#dev.off()

##Figure 2
svg(paste0("output/new_fig2.svg"),width=10.4,height=14.9) 
layout(matrix(c(1,1,1,1,5,5,5,5,9,9,9,9,1,1,1,1,5,5,5,5,9,9,9,9,2,2,2,2,6,6,6,6,10,10,10,10,2,2,2,2,6,6,6,6,10,10,10,10,3,3,3,3,7,7,7,7,11,11,11,11,3,3,3,3,7,7,7,7,11,11,11,11,4,4,4,4,8,8,8,8,12,12,12,12,4,4,4,4,8,8,8,8,12,12,12,12),nrow=8,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))

i=is[1]

#AIM responder vs rbd, vrld, TG, d28 mortality
for(k in virvars[c(3,1)])
{
  for(j in js[c(5,2,4,6)])
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
    
    if(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))<0.001)
    {
      pints=sprintf("%.1E",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
    } else {
      pints=round(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))),3)
    }
    if(pchisq(aa1,aa2, lower.tail = F)<0.001)
    {
      newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))				
    } else {
      newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints))		
    }
    
    plot(temp$logdays,temp$var,xlab="Days post-admission",pch=21,bg=adjustcolor("orange",alpha.f = 0.5),ylab=virnams[k],type='n',xaxt='n',las=2,main=j)
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
      legend(2.2,1,legend=c(newpvs,"Non-responder","Responder"),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)		
    } else {
      #		legend("topright",legend=c(newpvs,"Non-responder"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=j)
      legend(2.25,4.75,legend=c(newpvs,"Non-responder","Responder"),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)
    }
  }
}

plot(aims$x,aims$y,pch=20,cex=0.5,col=grps.cols2[match(aims$group,grps)],xaxt='n',yaxt='n',xlab="",ylab="")
lines(c(-30,25),c(10,20))
legend("bottomleft",c("Naive","CM","EM","TEMRA"),pch=20,col=grps.cols2[1:4],bty='n',pt.cex=1.5)
text(c(-28,-28),c(7.5,13),labels=c("CD4","CD8"))
arrows(-31,-32,27,-32,length=0.1)
arrows(-31,-32,-31,33,length=0.1)
mtext("tSNE-1",side=1,line=0.5,cex=0.75)
mtext("tSNE-2",side=2,line=0,cex=0.75)

plot(aims$x,aims$y,pch=20,cex=0.5,col=grps.aims.cols[match(aims$aim_num,grps.aims)],xlab="",ylab="",xaxt='n',yaxt='n')
lines(c(-30,25),c(10,20))
legend("bottomleft",legend=1:5,pch=20,col=grps.aims.cols,bty='n',pt.cex=1.5,title="AIM+ markers")
text(c(-28,-28),c(7.5,13),labels=c("CD4","CD8"))
arrows(-31,-32,27,-32,length=0.1)
arrows(-31,-32,-31,33,length=0.1)

for(k in aim.vars[1:2])
{
  tmp=table(aims$aim_num,aims[[k]])
  tt=barplot(t(tmp)/colSums(tmp)*100,beside=T,las=2,ylab="",xlab="",col=c("azure4","azure2"),xaxt='n')
  axis(1,colMeans(tt),labels=1:5,cex.axis=0.75,pos=2,tick=F)
  legend("topright",horiz=F,legend=aims.vars.labs[[k]],pch=22,pt.bg=c("azure4","azure2"),bty='n',pt.cex=1.5)
  mtext("AIM+ markers",side=1,line=1.5,cex=0.75)
  mtext("%",side=2,line=2.5,cex=0.75)
  if(summary(glm(aims[[k]]~aims$aim_num,family="binomial"))$coefficients[2,4]<0.001)
  {
    newpv=sprintf("%.1E",summary(glm(aims[[k]]~aims$aim_num,family="binomial"))$coefficients[2,4])
  } else {
    newpv=round(summary(glm(aims[[k]]~aims$aim_num,family="binomial"))$coefficients[2,4],3)
  }
  legend("right",legend=paste0("P = ",newpv),bty='n')
}

dev.off()	

##%cyto by memory subsets
#tt=table(aims$group,aims$aim_num)
#barplot(t(tt/rowSums(tt))*100,las=2,col=grps.aims.cols,names.arg=gsub("_","+",rownames(tt)),ylab="%")
#legend("topleft",legend=1:5,pch=22,pt.bg=grps.aims.cols,bty='n',pt.cex=1.5,title="AIM+ markers",horiz=T)

##Figure 3
svg(paste0("output/new_fig3.svg"),width=10.4,height=14.9) 
layout(matrix(c(1,1,1,1,5,5,5,5,9,9,9,9,1,1,1,1,5,5,5,5,9,9,9,9,2,2,2,2,6,6,6,6,10,10,10,10,2,2,2,2,6,6,6,6,10,10,10,10,3,3,3,3,7,7,7,7,11,11,11,11,3,3,3,3,7,7,7,7,11,11,11,11,4,4,4,4,8,8,8,8,12,12,12,12,4,4,4,4,8,8,8,8,13,13,13,13),nrow=8,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))

i=is[1]

for(k in virvars[c(12,16)])
{
  for(j in js[c(2,1,3,6)])
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
    if(summary(m)$coefficients[2,4]<0.001)
    {
      oldnewpvs=sprintf("%.1e",summary(m)$coefficients[2,4])
    } else {
      oldnewpvs=round(summary(m)$coefficients[2,4],3)
    }
    ymax=100
    #tmp=barplot(tmp2,main="",ylab="Responder (%)",xaxt='n',col=unlist(lapply(vircols[[k]],function(x){rep(x,1)})),las=2,xlab=virnams[k],ylim=c(0,ymax))
    #axis(1,tmp,c("FALSE","TRUE"),tick=F,las=1)
    #legend("topright",legend=paste0("P=",newpvs),bty='n')
    #text(tmp,rep(5,length(tmp)),labels=unlist(tttt))
    
    m=gam(responder~var+s(logdays,bs="cr",by=var),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$responder),],random =list(id=~1),family="binomial")
    a1=anova(m)
    aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
    aa2=sum(a1$s.table[-c(1), 2])
    
    if(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))<0.001)
    {
      pints=sprintf("%.1E",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
    } else {
      pints=round(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))),3)
    }
    if(pchisq(aa1,aa2, lower.tail = F)<0.001)
    {
      newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),oldnewpvs))				
    } else {
      newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),oldnewpvs))		
    }
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
    tmp=barplot(unlist(ttt),main=j,ylab=paste0("Responder (%)"),xaxt='n',col=unlist(lapply(vircols[[k]],function(x){rep(x,3)})),las=2,xlab="Days post-admission",ylim=c(0,ymax))
    axis(1,tmp,rep(tnams,lgth),tick=F,las=1)
    legend("topright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(22,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)
    text(tmp,rep(5,length(tmp)),labels=unlist(tttt))
  }
}

for(k in aim.vars[3:4])
{
  tmp=table(aims$aim_num,aims[[k]])
  tt=barplot(t(tmp)/colSums(tmp)*100,beside=T,las=2,ylab="",xlab="",col=c("azure4","azure2"),xaxt='n')
  axis(1,colMeans(tt),labels=1:5,cex.axis=0.75,pos=2,tick=F)
  legend("topright",horiz=F,legend=aims.vars.labs[[k]],pch=22,pt.bg=c("azure4","azure2"),bty='n',pt.cex=1.5)
  mtext("AIM+ markers",side=1,line=1.5,cex=0.75)
  mtext("%",side=2,line=2.5,cex=0.75)
  
  if(summary(glm(aims[[k]]~aims$aim_num,family="binomial"))$coefficients[2,4]<0.001)
  {
    newpv=sprintf("%.1E",summary(glm(aims[[k]]~aims$aim_num,family="binomial"))$coefficients[2,4])
  } else {
    newpv=round(summary(glm(aims[[k]]~aims$aim_num,family="binomial"))$coefficients[2,4],3)
  }
  legend("right",legend=paste0("P = ",newpv),bty='n')
}

cytocol=c("azure2",grps.cytos.cols)[match(aims[["cyto_pos_indiv"]],c("None",grps.cytos))]
cytocol[is.na(aims$group)]="white"
plot(aims$x,aims$y,pch=20,cex=0.5,col=cytocol,xlab="",ylab="",xaxt='n',yaxt='n')
legend("bottomleft",c("None",grps.cytos),pch=20,col=c("azure2",grps.cytos.cols),bty='n',pt.cex=1.5)
lines(c(-30,25),c(10,20))
text(c(-28,-28),c(7.5,13),labels=c("CD4","CD8"))
arrows(-31,-32,27,-32,length=0.1)
arrows(-31,-32,-31,33,length=0.1)
mtext("tSNE-1",side=1,line=0.5,cex=0.75)
mtext("tSNE-2",side=2,line=0,cex=0.75)

par(bty='n',mai=c(0.32,0.42,0.32,0.22))

tmp=table(aims$cyto_pos>0,aims$aim_num,aims$type2)
tmpp=tmp
tmpp[,,1][2,]=tmp[,,1][2,]/colSums(tmp[,,1])*100
tmpp[,,1][1,]=100-tmpp[,,1][2,]
tmpp[,,2][2,]=tmp[,,2][2,]/colSums(tmp[,,2])*100
tmpp[,,2][1,]=100-tmpp[,,2][2,]
tmpp4=tmpp[,,1][,match(grps.aims,colnames(tmpp[,,1]))]
tmpp8=tmpp[,,2][,match(grps.aims,colnames(tmpp[,,2]))]
tmp4=tmp[,,1][,match(grps.aims,colnames(tmp[,,1]))]
tmp8=tmp[,,2][,match(grps.aims,colnames(tmp[,,2]))]

tt=barplot(tmpp4,cex.names=0.75,yaxt='n',ylab="%",xlab="",col=c("azure4","azure2"),xaxt='n')
axis(2,seq(0,100,20),labels=seq(0,100,20),cex.axis=0.75,pos=0,las=2)
legend("topleft",horiz=T,legend=c("Cytokine-","Cytokine+"),pch=22,pt.bg=c("azure4","azure2"),bty='n',pt.cex=1.5)
axis(1,tt,labels=grps.aims,cex.axis=1,pos=10,tick=F)
mtext("CD4+ T cells",side=1,line=1.5,cex=0.75)
mtext("%",side=2,line=1.5,cex=0.75)
#summary(glm(aims$cyto_pos[which(aims$type2=="CD4")]>0~aims$aim_num[which(aims$type2=="CD4")],family="binomial"))$coefficients[2,4]

tt=barplot(tmpp8,cex.names=0.75,yaxt='n',ylab="%",xlab="",col=c("azure4","azure2"),xaxt='n')
axis(2,seq(0,100,20),labels=seq(0,100,20),cex.axis=0.75,pos=0,las=2)
axis(1,tt,labels=grps.aims,cex.axis=1,pos=10,tick=F)
mtext("CD8+ T cells",side=1,line=1.5,cex=0.75)
mtext("%",side=2,line=1.5,cex=0.75)
#summary(glm(aims$cyto_pos[which(aims$type2=="CD8")]>0~aims$aim_num[which(aims$type2=="CD8")],family="binomial"))$coefficients[2,4]

dev.off()

##Supp Figure 2
svg(paste0("output/new_supp_fig2.svg"),width=10.4,height=14.9) 
layout(matrix(c(1,1,1,2,2,2,5,5,5,6,6,6,1,1,1,2,2,2,5,5,5,6,6,6,3,3,3,4,4,4,7,7,7,8,8,8,3,3,3,4,4,4,7,7,7,8,8,8,9,9,9,10,10,10,13,13,13,14,14,14,9,9,9,10,10,10,13,13,13,14,14,14,11,11,11,12,12,12,15,15,15,16,16,16,11,11,11,12,12,12,15,15,15,16,16,16),nrow=8,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))

i=is[1]

#AIM responder vs rbd, vrld, TG, d28 mortality
for(k in virvars[c(3,1)])
{
  for(j in js[c(1,3)])
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
    
    if(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))<0.001)
    {
      pints=sprintf("%.1E",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
    } else {
      pints=round(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))),3)
    }
    if(pchisq(aa1,aa2, lower.tail = F)<0.001)
    {
      newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))				
    } else {
      newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints))		
    }
    
    plot(temp$logdays,temp$var,xlab="Days post-admission",pch=21,bg=adjustcolor("orange",alpha.f = 0.5),ylab=virnams[k],main=j,type='n',xaxt='n',las=2)
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
      legend(2.2,1,legend=c(newpvs,"Non-responder","Responder"),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1,cex=0.6)		
    } else {
      #		legend("topright",legend=c(newpvs,"Non-responder"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=j)
      legend(2.25,4.75,legend=c(newpvs,"Non-responder","Responder"),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1,cex=0.6)
    }
  }
}

for(k in aim.vars[1:2])
{
  for(j in 1:length(levels(aims[[k]])))
  {
    plot(aims$x,aims$y,pch=20,cex=0.5,col=grps.aims.cols[match(aims$aim_num,grps.aims)],xlab="",ylab="",xaxt='n',yaxt='n',type='n',main=aims.vars.labs[[k]][j])
    points(aims$x[which(aims[[k]]==levels(aims[[k]])[j])],aims$y[which(aims[[k]]==levels(aims[[k]])[j])],pch=20,cex=0.5,col=grps.aims.cols[match(aims$aim_num,grps.aims)][which(aims[[k]]==levels(aims[[k]])[j])])
    lines(c(-30,25),c(10,20))
    legend("bottomleft",legend=1:5,pch=20,col=grps.aims.cols,bty='n',pt.cex=1.5,title="AIM+ markers")
    text(c(-28,-28),c(7.5,13),labels=c("CD4","CD8"))
    arrows(-31,-32,27,-32,length=0.1)
    arrows(-31,-32,-31,33,length=0.1)
    
    if(k==aim.vars[1] & j==1)
    {
      mtext("tSNE-1",side=1,line=0.5,cex=0.75)
      mtext("tSNE-2",side=2,line=0,cex=0.75)
    }
    
  }
}

for(k in virvars[c(12,16)])
{
  for(j in js[c(5,4)])
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
    if(summary(m)$coefficients[2,4]<0.001)
    {
      oldnewpvs=sprintf("%.1e",summary(m)$coefficients[2,4])
    } else {
      oldnewpvs=round(summary(m)$coefficients[2,4],3)
    }
    ymax=100
    #tmp=barplot(tmp2,main="",ylab="Responder (%)",xaxt='n',col=unlist(lapply(vircols[[k]],function(x){rep(x,1)})),las=2,xlab=virnams[k],ylim=c(0,ymax))
    #axis(1,tmp,c("FALSE","TRUE"),tick=F,las=1)
    #legend("topright",legend=paste0("P=",newpvs),bty='n')
    #text(tmp,rep(5,length(tmp)),labels=unlist(tttt))
    
    m=gam(responder~var+s(logdays,bs="cr",by=var),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$responder),],random =list(id=~1),family="binomial")
    a1=anova(m)
    aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
    aa2=sum(a1$s.table[-c(1), 2])
    
    if(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))<0.001)
    {
      pints=sprintf("%.1E",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
    } else {
      pints=round(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))),3)
    }
    if(pchisq(aa1,aa2, lower.tail = F)<0.001)
    {
      newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),oldnewpvs))				
    } else {
      newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),oldnewpvs))		
    }
    
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
    tmp=barplot(unlist(ttt),ylab=paste0("Responder (%)"),xaxt='n',col=unlist(lapply(vircols[[k]],function(x){rep(x,3)})),las=2,xlab="Days post-admission",ylim=c(0,ymax),main=j)
    axis(1,tmp,rep(tnams,lgth),tick=F,las=1)
    legend("topright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(22,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1,cex=0.6)
    text(tmp,rep(5,length(tmp)),labels=unlist(tttt))
  }
}

for(k in aim.vars[3:4])
{
  for(j in 1:length(levels(aims[[k]])))
  {
    plot(aims$x,aims$y,pch=20,cex=0.5,col=grps.aims.cols[match(aims$aim_num,grps.aims)],xlab="",ylab="",xaxt='n',yaxt='n',type='n',main=aims.vars.labs[[k]][j])
    points(aims$x[which(aims[[k]]==levels(aims[[k]])[j])],aims$y[which(aims[[k]]==levels(aims[[k]])[j])],pch=20,cex=0.5,col=grps.aims.cols[match(aims$aim_num,grps.aims)][which(aims[[k]]==levels(aims[[k]])[j])])
    lines(c(-30,25),c(10,20))
    legend("bottomleft",legend=1:5,pch=20,col=grps.aims.cols,bty='n',pt.cex=1.5,title="AIM+ markers")
    text(c(-28,-28),c(7.5,13),labels=c("CD4","CD8"))
    arrows(-31,-32,27,-32,length=0.1)
    arrows(-31,-32,-31,33,length=0.1)
  }
}

dev.off()

#Supp Figure 4
svg(paste0("output/new_supp_fig4.svg"),width=10.4,height=14.9) 
layout(matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,15,15,15,16,16,16,17,17,17,18,18,18,19,19,19,20,20,20,21,21,21,22,22,22,23,23,23,24,24,24),nrow=6,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))

i=is[1]
for(j in js[c(1:6)])
{
  for(k in omics.vars[c(55,54,46,29)])
  {
    temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tcell.final.acute[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin)) #get data
    if(omics.type[k]=="cytof")
    {
      tmp=cytof[gsub(".*[-]","",cytof$event_id)=="1",]
      tmpp=cytof.pct[gsub(".*[-]","",cytof$event_id)=="1",]
      temp$bvar=tmpp[match(temp$id,tmp$participant_id),k]
    } else {
      tmp=olink[gsub(".*[-]","",olink$event_id)=="1",]
      temp$bvar=olink[match(temp$id,olink$participant_id),k]
    }
    temp$gate=as.numeric(temp$gate)
    temp$bvar2=log10(temp$bvar+min(temp$bvar[which(temp$bvar>0)]))	
    if(omics.type[k]=="cytof")
    {
      temp$var=cytof.pct[match(temp$event_id,cytof$event_id),k]
    } else {
      temp$var=olink[match(temp$event_id,olink$event_id),k]
    }	
    temp$var2=log10(temp$var+min(temp$var[which(temp$var>0)]))
    temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
    temp$livecd3=livecd3.final[[1]][match(temp$event_id,names(livecd3.final[[1]]))] #get stim live count
    temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
    temp$responder=factor(temp$gate>0)
    temp$logdays=log(as.numeric(temp$days))	
    
    yrng=range(c(temp$var+min(temp$var[which(temp$var>0)]),temp$bvar+min(temp$bvar[which(temp$bvar>0)])),na.rm=T)
    
    if(tcell.bomics.or.pv[[i]][k,j]<0.001)
    {
      newpvs1=sprintf("%.1E",tcell.bomics.or.pv[[i]][k,j])				
    } else {
      newpvs1=round(tcell.bomics.or.pv[[i]][k,j],3)		
    }
    if(tcell.omics.or.pv[[i]][k,j]<0.001)
    {
      newpvs2=sprintf("%.1E",tcell.omics.or.pv[[i]][k,j])				
    } else {
      newpvs2=round(tcell.omics.or.pv[[i]][k,j],3)		
    }		
    newpvs=c(newpvs1,newpvs2)
    
    if(yrng[1]>0)
    {
      boxplot(var2~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab="",ylab="",type='n',border="white",pch="",col="white",xaxt='n',log='y',las=2)
      par(new=T)
      boxplot(bvar2~responder,temp,xlim=c(0.5,5.5),ylim=log10(yrng),xlab=j,ylab=omics.nams2[k],type='n',border="white",pch="",col="white",yaxt='n',xaxt='n')
      stripchart(bvar2~responder,temp[which(temp$responder==F),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("azure4",alpha.f=0.5),cex=0.5,lwd=0.2)
      stripchart(bvar2~responder,temp[which(temp$responder==T),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("coral2",alpha.f=0.5),cex=0.5,lwd=0.2)
      par(new=T)
      boxplot(bvar2~responder,temp,xlim=c(0.5,5.5),ylim=log10(yrng),xlab="",ylab="",type='n',medcol="black",pch="",col=adjustcolor(c("azure4","coral2"),alpha.f=0.1),xaxt='n',yaxt='n')
      par(new=T)
      boxplot(var2~responder,temp,xlim=c(0.5,5.5),ylim=log10(yrng),xlab=j,ylab="",type='n',border="white",pch="",col="white",yaxt='n',xaxt='n',at=4:5)
      stripchart(var2~responder,temp[which(temp$responder==F),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("azure4",alpha.f=0.5),cex=0.5,at=4:5,lwd=0.2)
      stripchart(var2~responder,temp[which(temp$responder==T),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("coral2",alpha.f=0.5),cex=0.5,at=4:5,lwd=0.2)
      par(new=T)
      boxplot(var2~responder,temp,xlim=c(0.5,5.5),ylim=log10(yrng),xlab="",ylab="",type='n',medcol="black",pch="",col=adjustcolor(c("azure4","coral2"),alpha.f=0.25),xaxt='n',yaxt='n',at=4:5)
      axis(1,at=c(1:2,4:5),labels=rep(c("NR","R"),2),cex.axis=0.8)
    } else {
      yrng=range(c(temp$var,temp$bvar),na.rm=T)
      boxplot(bvar~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab=j,ylab=omics.nams2[k],type='n',border="white",pch="",col="white",xaxt='n',las=2)
      stripchart(bvar~responder,temp[which(temp$responder==F),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("azure4",alpha.f=0.5),cex=0.25,lwd=0.2)
      stripchart(bvar~responder,temp[which(temp$responder==T),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("coral2",alpha.f=0.5),cex=0.25,lwd=0.2)
      par(new=T)
      boxplot(bvar~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab="",ylab="",type='n',medcol="black",pch="",col=adjustcolor(c("azure4","coral2"),alpha.f=0.1),xaxt='n',yaxt='n')
      par(new=T)
      boxplot(var~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab=j,ylab="",type='n',border="white",pch="",col="white",yaxt='n',xaxt='n',at=4:5)
      stripchart(var~responder,temp[which(temp$responder==F),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("azure4",alpha.f=0.5),cex=0.25,at=4:5,lwd=0.2)
      stripchart(var~responder,temp[which(temp$responder==T),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("coral2",alpha.f=0.5),cex=0.25,at=4:5,lwd=0.2)
      par(new=T)
      boxplot(var~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab="",ylab="",type='n',medcol="black",pch="",col=adjustcolor(c("azure4","coral2"),alpha.f=0.25),xaxt='n',yaxt='n',at=4:5)
      axis(1,at=c(1:2,4:5),labels=rep(c("NR","R"),2),cex.axis=0.8)
    }
    abline(v=3,lty=3)
    axis(3,at=c(1.5,4.5),tick=F,labels=paste0(c("Baseline-P=","P="),newpvs),line=-0.5)
    
  }
}	

dev.off()

#Supp Figure 10
svg(paste0("output/new_supp_fig10.svg"),width=10.4,height=14.9) 
layout(matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,15,15,15,16,16,16,17,17,17,18,18,18,19,19,19,20,20,20,21,21,21,22,22,22,23,23,23,24,24,24),nrow=6,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))

i=is[1]
for(j in js[c(7:12)])
{
  for(k in omics.vars[c(51,48,15,29)])
  {
    temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tcell.final.acute[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin)) #get data
    if(omics.type[k]=="cytof")
    {
      tmp=cytof[gsub(".*[-]","",cytof$event_id)=="1",]
      tmpp=cytof.pct[gsub(".*[-]","",cytof$event_id)=="1",]
      temp$bvar=tmpp[match(temp$id,tmp$participant_id),k]
    } else {
      tmp=olink[gsub(".*[-]","",olink$event_id)=="1",]
      temp$bvar=olink[match(temp$id,olink$participant_id),k]
    }
    temp$gate=as.numeric(temp$gate)
    temp$bvar2=log10(temp$bvar+min(temp$bvar[which(temp$bvar>0)]))	
    if(omics.type[k]=="cytof")
    {
      temp$var=cytof.pct[match(temp$event_id,cytof$event_id),k]
    } else {
      temp$var=olink[match(temp$event_id,olink$event_id),k]
    }	
    temp$var2=log10(temp$var+min(temp$var[which(temp$var>0)]))
    temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
    temp$livecd3=livecd3.final[[1]][match(temp$event_id,names(livecd3.final[[1]]))] #get stim live count
    temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
    temp$responder=factor(temp$gate>0)
    temp$logdays=log(as.numeric(temp$days))	
    
    yrng=range(c(temp$var+min(temp$var[which(temp$var>0)]),temp$bvar+min(temp$bvar[which(temp$bvar>0)])),na.rm=T)
    if(tcell.bomics.or.pv[[i]][k,j]<0.001)
    {
      newpvs1=sprintf("%.1E",tcell.bomics.or.pv[[i]][k,j])				
    } else {
      newpvs1=round(tcell.bomics.or.pv[[i]][k,j],3)		
    }
    if(tcell.omics.or.pv[[i]][k,j]<0.001)
    {
      newpvs2=sprintf("%.1E",tcell.omics.or.pv[[i]][k,j])				
    } else {
      newpvs2=round(tcell.omics.or.pv[[i]][k,j],3)		
    }		
    newpvs=c(newpvs1,newpvs2)
    if(yrng[1]>0)
    {
      boxplot(var2~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab="",ylab="",type='n',border="white",pch="",col="white",xaxt='n',log='y',las=2)
      par(new=T)
      boxplot(bvar2~responder,temp,xlim=c(0.5,5.5),ylim=log10(yrng),xlab=j,ylab=omics.nams2[k],type='n',border="white",pch="",col="white",yaxt='n',xaxt='n')
      stripchart(bvar2~responder,temp[which(temp$responder==F),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("azure4",alpha.f=0.5),cex=0.5,lwd=0.2)
      stripchart(bvar2~responder,temp[which(temp$responder==T),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("coral2",alpha.f=0.5),cex=0.5,lwd=0.2)
      par(new=T)
      boxplot(bvar2~responder,temp,xlim=c(0.5,5.5),ylim=log10(yrng),xlab="",ylab="",type='n',medcol="black",pch="",col=adjustcolor(c("azure4","coral2"),alpha.f=0.1),xaxt='n',yaxt='n')
      par(new=T)
      boxplot(var2~responder,temp,xlim=c(0.5,5.5),ylim=log10(yrng),xlab=j,ylab="",type='n',border="white",pch="",col="white",yaxt='n',xaxt='n',at=4:5)
      stripchart(var2~responder,temp[which(temp$responder==F),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("azure4",alpha.f=0.5),cex=0.5,at=4:5,lwd=0.2)
      stripchart(var2~responder,temp[which(temp$responder==T),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("coral2",alpha.f=0.5),cex=0.5,at=4:5,lwd=0.2)
      par(new=T)
      boxplot(var2~responder,temp,xlim=c(0.5,5.5),ylim=log10(yrng),xlab="",ylab="",type='n',medcol="black",pch="",col=adjustcolor(c("azure4","coral2"),alpha.f=0.25),xaxt='n',yaxt='n',at=4:5)
      axis(1,at=c(1:2,4:5),labels=rep(c("NR","R"),2),cex.axis=0.8)
    } else {
      yrng=range(c(temp$var,temp$bvar),na.rm=T)
      boxplot(bvar~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab=j,ylab=omics.nams2[k],type='n',border="white",pch="",col="white",xaxt='n',las=2)
      stripchart(bvar~responder,temp[which(temp$responder==F),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("azure4",alpha.f=0.5),cex=0.25,lwd=0.2)
      stripchart(bvar~responder,temp[which(temp$responder==T),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("coral2",alpha.f=0.5),cex=0.25,lwd=0.2)
      par(new=T)
      boxplot(bvar~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab="",ylab="",type='n',medcol="black",pch="",col=adjustcolor(c("azure4","coral2"),alpha.f=0.1),xaxt='n',yaxt='n')
      par(new=T)
      boxplot(var~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab=j,ylab="",type='n',border="white",pch="",col="white",yaxt='n',xaxt='n',at=4:5)
      stripchart(var~responder,temp[which(temp$responder==F),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("azure4",alpha.f=0.5),cex=0.25,at=4:5,lwd=0.2)
      stripchart(var~responder,temp[which(temp$responder==T),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("coral2",alpha.f=0.5),cex=0.25,at=4:5,lwd=0.2)
      par(new=T)
      boxplot(var~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab="",ylab="",type='n',medcol="black",pch="",col=adjustcolor(c("azure4","coral2"),alpha.f=0.25),xaxt='n',yaxt='n',at=4:5)
      axis(1,at=c(1:2,4:5),labels=rep(c("NR","R"),2),cex.axis=0.8)
    }
    abline(v=3,lty=3)
    axis(3,at=c(1.5,4.5),tick=F,labels=paste0(c("Baseline-P=","P="),newpvs),line=-0.5)
    
  }
}	

dev.off()

##Figure 4

svg(paste0("output/new_fig4.svg"),width=10.4,height=14.9) 
#layout(matrix(c(1,1,1,1,2,2,3,3,4,4,5,5,1,1,1,1,2,2,3,3,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8,7,7,7,7,7,7,8,8,8,8,8,8,9,9,10,10,10,10,11,11,12,12,12,12,9,9,10,10,10,10,11,11,12,12,12,12,13,13,13,13,13,13,14,14,14,15,15,15,13,13,13,13,13,13,16,16,16,17,17,17),nrow=8,byrow=T))
layout(matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,1,1,1,1,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8,6,6,6,6,7,7,7,7,8,8,8,8,9,9,9,9,10,10,10,10,11,11,11,11,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,13,13,13,13,14,15,15,16,16,12,12,12,13,13,13,13,14,17,17,18,18),nrow=8,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))

i=is[1]

#TFH vs gates
plot.new()

#tfh vs RBD boxplots
k=virvars[3]
for(j in tfh.js[c(2,7,10,8)])
{
  
  temp=data.frame(cbind("days"=tfh.info.acute$event_date,"days_symptom"=tfh.info.acute$event_date_symptom,"gate"=unlist(tfh.final.acute[[i]][,j]),"id"=tfh.info.acute$participant_id,"event_id"=tfh.info.acute$event_id,"symptom_bin"=tfh.info.acute$symptom_date_bin)) #get data
  temp$var=tfh.info.acute[match(temp$event_id,tfh.info.acute$event_id),k]	#get variable
  temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
  temp$livecd3=livecd3.final[[i]][match(temp$event_id,names(livecd3.final[[i]]))] #get stim live count
  temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
  temp$gate=as.numeric(temp$gate)
  temp$days=as.numeric(temp$days)
  temp$days_symptom=as.numeric(temp$days_symptom)
  
  temp$gate3=log2(temp$gate)  #as above but don't add minimum value
  temp$gate3[is.infinite(temp$gate3)]=NA
  
  temp$logdays=log(temp$days) #needed for lme
  temp$logdays_symptom=log(temp$days_symptom) #needed for lme
  temp=temp[!is.na(temp$var)&!is.na(temp$gate),]
  temp$responder=factor(temp$gate>0)
  
  m=lme(var~responder,data=temp,random=~1|id)
  
  boxplot(temp$var~temp$responder,ylab=virnams[k],xlab=tfh.js3b[j],type='n',pch="",col="white",las=2,xaxt='n',border="white")
  axis(1,at=1:2,labels=c("Non-responder","Responder"))
  stripchart(temp$var[which(temp$responder==F)]~temp$responder[which(temp$responder==F)],add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor(vircols[[k]],alpha.f=0.5)[1],jitter=0.25,cex=0.7)
  stripchart(temp$var[which(temp$responder==T)]~temp$responder[which(temp$responder==T)],add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor(vircols[[k]],alpha.f=0.5)[2],jitter=0.25,cex=0.7)
  par(new=T)
  boxplot(temp$var~temp$responder,xlab="",ylab="",type='n',pch="",col=adjustcolor(vircols[[k]],alpha.f=0.1),yaxt='n',las=2,xaxt='n',medcol="black")
  axis(1,at=c(-1,3))
  if(summary(m)$tTable[2,5]<0.001)
  {
    axis(3,at=1.5,line=-0.5,labels=paste0("P=",sprintf("%.1e",summary(m)$tTable[2,5])),tick=F,cex.axis=0.8)
  } else {
    axis(3,at=1.5,line=-0.5,labels=paste0("P=",round(summary(m)$tTable[2,5],3)),tick=F,cex.axis=0.8)		
  }
  
}

#tfh vs vrld/tg
for(k in virvars[c(2,12)])
{
  for(j in tfh.js[c(2,7,8)])
  {
    temp=data.frame(cbind("days"=tfh.info.acute$event_date,"days_symptom"=tfh.info.acute$event_date_symptom,"gate"=unlist(tfh.final.acute[[i]][,j]),"id"=tfh.info.acute$participant_id,"event_id"=tfh.info.acute$event_id,"symptom_bin"=tfh.info.acute$symptom_date_bin)) #get data
    temp$var=tfh.info.acute[match(temp$event_id,tfh.info.acute$event_id),k]	#get variable
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
    
    if(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))<0.001)
    {
      pints=sprintf("%.1E",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
    } else {
      pints=round(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))),3)
    }
    if(pchisq(aa1,aa2, lower.tail = F)<0.001)
    {
      newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))				
    } else {
      newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints))		
    }
    
    plot(temp$logdays,temp$gate3,xlab="Days post-admission",yaxt='n',xaxt='n',pch=21,bg=adjustcolor("azure4",alpha.f = 0.5),ylab=paste0("(%)"),type='n',main=tfh.js3b[j])
    points(temp$logdays,temp$gate3,bg=temp$col,pch=21,cex=1)
    axis(2,at=log2(c(0.001,0.01,0.10,1.00,10.00,100.00)),labels=c(0.001,0.01,0.1,1,10,100),las=2)
    if(j%in%tfh.js[c(2,10)])
    {
      legend("bottomright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)		
    } else {
      legend("topleft",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)
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
}

#CD8 Tfh vs rbd, vrld and tg
k=virvars[3]
j=tfh.js[11]

temp=data.frame(cbind("days"=tfh.info.acute$event_date,"days_symptom"=tfh.info.acute$event_date_symptom,"gate"=unlist(tfh.final.acute[[i]][,j]),"id"=tfh.info.acute$participant_id,"event_id"=tfh.info.acute$event_id,"symptom_bin"=tfh.info.acute$symptom_date_bin)) #get data
temp$var=tfh.info.acute[match(temp$event_id,tfh.info.acute$event_id),k]	#get variable
temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
temp$livecd3=livecd3.final[[i]][match(temp$event_id,names(livecd3.final[[i]]))] #get stim live count
temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
temp$gate=as.numeric(temp$gate)
temp$days=as.numeric(temp$days)
temp$days_symptom=as.numeric(temp$days_symptom)

temp$gate3=log2(temp$gate)  #as above but don't add minimum value
temp$gate3[is.infinite(temp$gate3)]=NA

temp$logdays=log(temp$days) #needed for lme
temp$logdays_symptom=log(temp$days_symptom) #needed for lme
temp=temp[!is.na(temp$var)&!is.na(temp$gate),]
temp$responder=factor(temp$gate>0)

m=lme(var~responder,data=temp,random=~1|id)

boxplot(temp$var~temp$responder,ylab=virnams[k],xlab=tfh.js3b[j],type='n',pch="",col="white",las=2,xaxt='n',border="white")
axis(1,at=1:2,labels=c("Non-responder","Responder"))
stripchart(temp$var[which(temp$responder==F)]~temp$responder[which(temp$responder==F)],add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor(vircols[[k]],alpha.f=0.5)[1],jitter=0.25,cex=0.7)
stripchart(temp$var[which(temp$responder==T)]~temp$responder[which(temp$responder==T)],add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor(vircols[[k]],alpha.f=0.5)[2],jitter=0.25,cex=0.7)
par(new=T)
boxplot(temp$var~temp$responder,xlab="",ylab="",type='n',pch="",col=adjustcolor(vircols[[k]],alpha.f=0.1),yaxt='n',las=2,xaxt='n',medcol="black")
axis(1,at=c(-1,3))
if(summary(m)$tTable[2,5]<0.001)
{
  axis(3,at=1.5,line=-0.5,labels=paste0("P=",sprintf("%.1e",summary(m)$tTable[2,5])),tick=F,cex.axis=0.8)
} else {
  axis(3,at=1.5,line=-0.5,labels=paste0("P=",round(summary(m)$tTable[2,5],3)),tick=F,cex.axis=0.8)		
}

#k=virvars[12]
for(k in virvars[c(12)])
{
  temp=data.frame(cbind("days"=tfh.info.acute$event_date,"days_symptom"=tfh.info.acute$event_date_symptom,"gate"=unlist(tfh.final.acute[[i]][,j]),"id"=tfh.info.acute$participant_id,"event_id"=tfh.info.acute$event_id,"symptom_bin"=tfh.info.acute$symptom_date_bin)) #get data
  temp$var=tfh.info.acute[match(temp$event_id,tfh.info.acute$event_id),k]	#get variable
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
  
  if(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))<0.001)
  {
    pints=sprintf("%.1E",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
  } else {
    pints=round(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))),3)
  }
  if(pchisq(aa1,aa2, lower.tail = F)<0.001)
  {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))				
  } else {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints))		
  }
  
  plot(temp$logdays,temp$gate3,xlab="Days post-admission",yaxt='n',xaxt='n',pch=21,bg=adjustcolor("azure4",alpha.f = 0.5),ylab=paste0("(%)"),type='n',main=tfh.js3b[j])
  points(temp$logdays,temp$gate3,bg=temp$col,pch=21,cex=1)
  axis(2,at=log2(c(0.001,0.01,0.10,1.00,10.00,100.00)),labels=c(0.001,0.01,0.1,1,10,100),las=2)
  if(j%in%tfh.js[c(2,10)])
  {
    legend("bottomright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)		
  } else {
    legend("topright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)
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

#Tfh vs omics
plot.new()
for(j in tfh.js[c(7,8)])
{
  for(k in omics.vars[c(32,35,38)])
  {
    if(j==tfh.js[8]&k==omics.vars[32]|j==tfh.js[7]&k==omics.vars[35])
    {
    } else {
      
      temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tfh.final.acute.merge[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin)) #get data
      if(omics.type[k]=="cytof")
      {
        tmp=cytof[gsub(".*[-]","",cytof$event_id)=="1",]
        tmpp=cytof.pct[gsub(".*[-]","",cytof$event_id)=="1",]
        temp$bvar=tmpp[match(temp$id,tmp$participant_id),k]
      } else {
        tmp=olink[gsub(".*[-]","",olink$event_id)=="1",]
        temp$bvar=olink[match(temp$id,olink$participant_id),k]
      }
      temp$gate=as.numeric(temp$gate)
      temp$bvar2=log10(temp$bvar+min(temp$bvar[which(temp$bvar>0)]))	
      if(omics.type[k]=="cytof")
      {
        temp$var=cytof.pct[match(temp$event_id,cytof$event_id),k]
      } else {
        temp$var=olink[match(temp$event_id,olink$event_id),k]
      }	
      temp$var2=log10(temp$var+min(temp$var[which(temp$var>0)]))
      temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
      temp$livecd3=livecd3.final[[1]][match(temp$event_id,names(livecd3.final[[1]]))] #get stim live count
      temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
      temp$responder=factor(temp$gate>0)
      temp$logdays=log(as.numeric(temp$days))	
      
      yrng=range(c(temp$var+min(temp$var[which(temp$var>0)]),temp$bvar+min(temp$bvar[which(temp$bvar>0)])),na.rm=T)
      if(tfh.bomics.or.pv[[i]][k,j]<0.001)
      {
        newpvs1=sprintf("%.1E",tfh.bomics.or.pv[[i]][k,j])				
      } else {
        newpvs1=round(tfh.bomics.or.pv[[i]][k,j],3)		
      }
      if(tfh.omics.or.pv[[i]][k,j]<0.001)
      {
        newpvs2=sprintf("%.1E",tfh.omics.or.pv[[i]][k,j])				
      } else {
        newpvs2=round(tfh.omics.or.pv[[i]][k,j],3)		
      }		
      newpvs=c(newpvs1,newpvs2)
      if(yrng[1]>0)
      {
        boxplot(var2~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab="",ylab="",type='n',border="white",pch="",col="white",xaxt='n',log='y',las=2)
        par(new=T)
        boxplot(bvar2~responder,temp,xlim=c(0.5,5.5),ylim=log10(yrng),xlab=tfh.js3b[j],ylab=omics.nams2[k],type='n',border="white",pch="",col="white",yaxt='n',xaxt='n')
        stripchart(bvar2~responder,temp[which(temp$responder==F),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("azure4",alpha.f=0.5),cex=0.5,lwd=0.2)
        stripchart(bvar2~responder,temp[which(temp$responder==T),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("coral2",alpha.f=0.5),cex=0.5,lwd=0.2)
        par(new=T)
        boxplot(bvar2~responder,temp,xlim=c(0.5,5.5),ylim=log10(yrng),xlab="",ylab="",type='n',medcol="black",pch="",col=adjustcolor(c("azure4","coral2"),alpha.f=0.1),xaxt='n',yaxt='n')
        par(new=T)
        boxplot(var2~responder,temp,xlim=c(0.5,5.5),ylim=log10(yrng),xlab=tfh.js3b[j],ylab="",type='n',border="white",pch="",col="white",yaxt='n',xaxt='n',at=4:5)
        stripchart(var2~responder,temp[which(temp$responder==F),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("azure4",alpha.f=0.5),cex=0.5,at=4:5,lwd=0.2)
        stripchart(var2~responder,temp[which(temp$responder==T),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("coral2",alpha.f=0.5),cex=0.5,at=4:5,lwd=0.2)
        par(new=T)
        boxplot(var2~responder,temp,xlim=c(0.5,5.5),ylim=log10(yrng),xlab="",ylab="",type='n',medcol="black",pch="",col=adjustcolor(c("azure4","coral2"),alpha.f=0.25),xaxt='n',yaxt='n',at=4:5)
        axis(1,at=c(1:2,4:5),labels=rep(c("NR","R"),2),cex.axis=0.5)
      } else {
        yrng=range(c(temp$var,temp$bvar),na.rm=T)
        boxplot(bvar~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab=tfh.js3b[j],ylab=omics.nams2[k],type='n',border="white",pch="",col="white",xaxt='n',las=2)
        stripchart(bvar~responder,temp[which(temp$responder==F),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("azure4",alpha.f=0.5),cex=0.25,lwd=0.2)
        stripchart(bvar~responder,temp[which(temp$responder==T),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("coral2",alpha.f=0.5),cex=0.25,lwd=0.2)
        par(new=T)
        boxplot(bvar~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab="",ylab="",type='n',medcol="black",pch="",col=adjustcolor(c("azure4","coral2"),alpha.f=0.1),xaxt='n',yaxt='n')
        par(new=T)
        boxplot(var~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab=tfh.js3b[j],ylab="",type='n',border="white",pch="",col="white",yaxt='n',xaxt='n',at=4:5)
        stripchart(var~responder,temp[which(temp$responder==F),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("azure4",alpha.f=0.5),cex=0.25,at=4:5,lwd=0.2)
        stripchart(var~responder,temp[which(temp$responder==T),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("coral2",alpha.f=0.5),cex=0.25,at=4:5,lwd=0.2)
        par(new=T)
        boxplot(var~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab="",ylab="",type='n',medcol="black",pch="",col=adjustcolor(c("azure4","coral2"),alpha.f=0.25),xaxt='n',yaxt='n',at=4:5)
        axis(1,at=c(1:2,4:5),labels=rep(c("NR","R"),2),cex.axis=0.5)
      }
      abline(v=3,lty=3)
      axis(3,at=c(1.5,4.5),tick=F,labels=paste0(c("Baseline-P=","P="),newpvs),line=-0.5)
    }}
}

dev.off()
##see 2024.09.20_tfh_resoonder_omics.R for details

#bottom-right heatmap
i=is[1]
tmp.or=cbind(tfh.bomics.or[[i]][,c(11,9,2,4,8,3,7,10,5,6)],tfh.omics.or[[i]][,c(11,9,2,4,8,3,7,10,5,6)])
tmp.sigs=cbind(tfh.bomics.sigs[[i]][,c(11,9,2,4,8,3,7,10,5,6)],tfh.omics.sigs[[i]][,c(11,9,2,4,8,3,7,10,5,6)])
lgth=length(which(apply(tmp.sigs,1,function(x){sum(x!="")})>0))
#pheatmap(t(log2(tmp.or[which(apply(tmp.sigs,1,function(x){sum(x!=""&x!="*")})>0),])),display_numbers=t(tmp.sigs[which(apply(tmp.sigs,1,function(x){sum(x!=""&x!="*")})>0),]),cluster_rows=F,cluster_cols=T,labels_row=rep(tfh.js3[-2][c(10,8,1,3,7,2,6,9,4,5)],2),angle=45,labels_col=omics.nams2[which(apply(tmp.sigs,1,function(x){sum(x!=""&x!="*")})>0)],fontsize_number=8,filename=paste0("output/local/tfh/heatmap_combined_omics_median_sigs_flip_",i,"_trim.png"),height=(24.4/round(59/lgth)),width=10.725,fontsize_col=6,fontsize_row=6,legend=T,cellwidth=12,cellheight=12,breaks=seq(-2,2,0.05),color=colorRampPalette(rev(ryb[2:10]))(length(seq(-2,2,0.05))),fontsize=6,gaps_row=10)
pheatmap(t(log2(tmp.or[which(apply(tmp.sigs,1,function(x){sum(x!="")})>0),])),display_numbers=t(tmp.sigs[which(apply(tmp.sigs,1,function(x){sum(x!="")})>0),]),cluster_rows=F,cluster_cols=T,labels_row=rep(tfh.js3[-2][c(10,8,1,3,7,2,6,9,4,5)],2),angle=45,labels_col=omics.nams2[which(apply(tmp.sigs,1,function(x){sum(x!="")})>0)],fontsize_number=8,filename=paste0("output/local/tfh/heatmap_combined_omics_median_sigs_flip_",i,"_trim.png"),height=(24.4/round(59/lgth)),width=10.725,fontsize_col=6,fontsize_row=6,legend=T,cellwidth=12,cellheight=12,breaks=seq(-2,2,0.05),color=colorRampPalette(rev(ryb[2:10]))(length(seq(-2,2,0.05))),fontsize=6,gaps_row=10)


##Supp Figure 5

svg(paste0("output/new_supp_fig5.svg"),width=10.4,height=14.9) 
#layout(matrix(c(1,1,1,1,2,2,3,3,4,4,5,5,1,1,1,1,2,2,3,3,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8,7,7,7,7,7,7,8,8,8,8,8,8,9,9,10,10,10,10,11,11,12,12,12,12,9,9,10,10,10,10,11,11,12,12,12,12,13,13,13,13,13,13,14,14,14,15,15,15,13,13,13,13,13,13,16,16,16,17,17,17),nrow=8,byrow=T))
layout(matrix(c(1,1,1,1,5,5,5,5,9,9,9,9,1,1,1,1,5,5,5,5,9,9,9,9,2,2,2,2,6,6,6,6,10,10,10,10,2,2,2,2,6,6,6,6,10,10,10,10,3,3,3,3,7,7,7,7,11,11,11,11,3,3,3,3,7,7,7,7,11,11,11,11,4,4,4,4,8,8,8,8,12,12,12,12,4,4,4,4,8,8,8,8,12,12,12,12),nrow=8,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))

i=is[1]

k=virvars[16]

for(j in tfh.js[c(2,7,8,11)])
{
  temp=data.frame(cbind("days"=tfh.info.acute$event_date,"days_symptom"=tfh.info.acute$event_date_symptom,"gate"=unlist(tfh.final.acute[[i]][,j]),"id"=tfh.info.acute$participant_id,"event_id"=tfh.info.acute$event_id,"symptom_bin"=tfh.info.acute$symptom_date_bin)) #get data
  temp$var=tfh.info.acute[match(temp$event_id,tfh.info.acute$event_id),k]	#get variable
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
  
  if(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))<0.001)
  {
    pints=sprintf("%.1E",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
  } else {
    pints=round(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))),3)
  }
  if(pchisq(aa1,aa2, lower.tail = F)<0.001)
  {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))				
  } else {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints))		
  }	
  
  plot(temp$logdays,temp$gate3,xlab="Days post-admission",yaxt='n',xaxt='n',pch=21,bg=adjustcolor("azure4",alpha.f = 0.5),ylab=paste0("%"),type='n',main=tfh.js3b[j])
  points(temp$logdays,temp$gate3,bg=temp$col,pch=21,cex=1)
  axis(2,at=log2(c(0.001,0.01,0.10,1.00,10.00,100.00)),labels=c(0.001,0.01,0.1,1,10,100),las=2)
  if(j%in%tfh.js[c(2,10)])
  {
    legend("bottomright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)		
  } else {
    legend("topleft",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)
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

for(k in virvars[c(2,12)])
{
  for(j in tfh.js[c(10,3,6,5)])
  {
    temp=data.frame(cbind("days"=tfh.info.acute$event_date,"days_symptom"=tfh.info.acute$event_date_symptom,"gate"=unlist(tfh.final.acute[[i]][,j]),"id"=tfh.info.acute$participant_id,"event_id"=tfh.info.acute$event_id,"symptom_bin"=tfh.info.acute$symptom_date_bin)) #get data
    temp$var=tfh.info.acute[match(temp$event_id,tfh.info.acute$event_id),k]	#get variable
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
    
    if(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))<0.001)
    {
      pints=sprintf("%.1E",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
    } else {
      pints=round(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))),3)
    }
    if(pchisq(aa1,aa2, lower.tail = F)<0.001)
    {
      newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))				
    } else {
      newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints))		
    }
    
    plot(temp$logdays,temp$gate3,xlab="Days post-admission",yaxt='n',xaxt='n',pch=21,bg=adjustcolor("azure4",alpha.f = 0.5),ylab=paste0("%"),type='n',main=tfh.js3b[j])
    points(temp$logdays,temp$gate3,bg=temp$col,pch=21,cex=1)
    axis(2,at=log2(c(0.001,0.01,0.10,1.00,10.00,100.00)),labels=c(0.001,0.01,0.1,1,10,100),las=2)
    if(j%in%tfh.js[c(2,10)])
    {
      legend("bottomright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)		
    } else {
      legend("topleft",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)
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
}

dev.off()

##Supp Figure 6
svg(paste0("output/new_supp_fig6.svg"),width=10.4,height=14.9) 
layout(matrix(c(1,1,1,1,1,1,1,1,6,2,2,2,1,1,1,1,1,1,1,1,6,3,3,3,1,1,1,1,1,1,1,1,6,4,4,4,1,1,1,1,1,1,1,1,6,5,5,5,rep(6,48)),nrow=8,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))

#heatmap of tfh vs omics
plot.new()

for(j in tfh.js[c(2,10)])
{
  for(k in omics.vars[c(35,38)])
  {
    temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tfh.final.acute.merge[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin)) #get data
    if(omics.type[k]=="cytof")
    {
      tmp=cytof[gsub(".*[-]","",cytof$event_id)=="1",]
      tmpp=cytof.pct[gsub(".*[-]","",cytof$event_id)=="1",]
      temp$bvar=tmpp[match(temp$id,tmp$participant_id),k]
    } else {
      tmp=olink[gsub(".*[-]","",olink$event_id)=="1",]
      temp$bvar=olink[match(temp$id,olink$participant_id),k]
    }
    temp$gate=as.numeric(temp$gate)
    temp$bvar2=log10(temp$bvar+min(temp$bvar[which(temp$bvar>0)]))	
    if(omics.type[k]=="cytof")
    {
      temp$var=cytof.pct[match(temp$event_id,cytof$event_id),k]
    } else {
      temp$var=olink[match(temp$event_id,olink$event_id),k]
    }	
    temp$var2=log10(temp$var+min(temp$var[which(temp$var>0)]))
    temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
    temp$livecd3=livecd3.final[[1]][match(temp$event_id,names(livecd3.final[[1]]))] #get stim live count
    temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
    temp$responder=factor(temp$gate>0)
    temp$logdays=log(as.numeric(temp$days))	
    
    yrng=range(c(temp$var+min(temp$var[which(temp$var>0)]),temp$bvar+min(temp$bvar[which(temp$bvar>0)])),na.rm=T)
    if(tfh.bomics.or.pv[[i]][k,j]<0.001)
    {
      newpvs1=sprintf("%.1E",tfh.bomics.or.pv[[i]][k,j])				
    } else {
      newpvs1=round(tfh.bomics.or.pv[[i]][k,j],3)		
    }
    if(tfh.omics.or.pv[[i]][k,j]<0.001)
    {
      newpvs2=sprintf("%.1E",tfh.omics.or.pv[[i]][k,j])				
    } else {
      newpvs2=round(tfh.omics.or.pv[[i]][k,j],3)		
    }		
    newpvs=c(newpvs1,newpvs2)
    if(yrng[1]>0)
    {
      boxplot(var2~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab="",ylab="",type='n',border="white",pch="",col="white",xaxt='n',log='y',las=2)
      par(new=T)
      boxplot(bvar2~responder,temp,xlim=c(0.5,5.5),ylim=log10(yrng),xlab=tfh.js3b[j],ylab=omics.nams2[k],type='n',border="white",pch="",col="white",yaxt='n',xaxt='n')
      stripchart(bvar2~responder,temp[which(temp$responder==F),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("azure4",alpha.f=0.5),cex=0.5,lwd=0.2)
      stripchart(bvar2~responder,temp[which(temp$responder==T),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("coral2",alpha.f=0.5),cex=0.5,lwd=0.2)
      par(new=T)
      boxplot(bvar2~responder,temp,xlim=c(0.5,5.5),ylim=log10(yrng),xlab="",ylab="",type='n',medcol="black",pch="",col=adjustcolor(c("azure4","coral2"),alpha.f=0.1),xaxt='n',yaxt='n')
      par(new=T)
      boxplot(var2~responder,temp,xlim=c(0.5,5.5),ylim=log10(yrng),xlab=tfh.js3b[j],ylab="",type='n',border="white",pch="",col="white",yaxt='n',xaxt='n',at=4:5)
      stripchart(var2~responder,temp[which(temp$responder==F),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("azure4",alpha.f=0.5),cex=0.5,at=4:5,lwd=0.2)
      stripchart(var2~responder,temp[which(temp$responder==T),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("coral2",alpha.f=0.5),cex=0.5,at=4:5,lwd=0.2)
      par(new=T)
      boxplot(var2~responder,temp,xlim=c(0.5,5.5),ylim=log10(yrng),xlab="",ylab="",type='n',medcol="black",pch="",col=adjustcolor(c("azure4","coral2"),alpha.f=0.25),xaxt='n',yaxt='n',at=4:5)
      axis(1,at=c(1:2,4:5),labels=rep(c("NR","R"),2),cex.axis=0.8)
    } else {
      yrng=range(c(temp$var,temp$bvar),na.rm=T)
      boxplot(bvar~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab=tfh.js3b[j],ylab=omics.nams2[k],type='n',border="white",pch="",col="white",xaxt='n',las=2)
      stripchart(bvar~responder,temp[which(temp$responder==F),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("azure4",alpha.f=0.5),cex=0.25,lwd=0.2)
      stripchart(bvar~responder,temp[which(temp$responder==T),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("coral2",alpha.f=0.5),cex=0.25,lwd=0.2)
      par(new=T)
      boxplot(bvar~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab="",ylab="",type='n',medcol="black",pch="",col=adjustcolor(c("azure4","coral2"),alpha.f=0.1),xaxt='n',yaxt='n')
      par(new=T)
      boxplot(var~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab=tfh.js3b[j],ylab="",type='n',border="white",pch="",col="white",yaxt='n',xaxt='n',at=4:5)
      stripchart(var~responder,temp[which(temp$responder==F),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("azure4",alpha.f=0.5),cex=0.25,at=4:5,lwd=0.2)
      stripchart(var~responder,temp[which(temp$responder==T),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("coral2",alpha.f=0.5),cex=0.25,at=4:5,lwd=0.2)
      par(new=T)
      boxplot(var~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab="",ylab="",type='n',medcol="black",pch="",col=adjustcolor(c("azure4","coral2"),alpha.f=0.25),xaxt='n',yaxt='n',at=4:5)
      axis(1,at=c(1:2,4:5),labels=rep(c("NR","R"),2),cex.axis=0.8)
    }
    abline(v=3,lty=3)
    axis(3,at=c(1.5,4.5),tick=F,labels=paste0(c("Baseline-P=","P="),newpvs),line=-0.5)
  }
}	

dev.off()

##Figure 5
svg(paste0("output/new_fig5.svg"),width=10.4,height=14.9) 
#layout(matrix(c(1,1,1,1,2,2,3,3,4,4,5,5,1,1,1,1,2,2,3,3,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8,7,7,7,7,7,7,8,8,8,8,8,8,9,9,10,10,10,10,11,11,12,12,12,12,9,9,10,10,10,10,11,11,12,12,12,12,13,13,13,13,13,13,14,14,14,15,15,15,13,13,13,13,13,13,16,16,16,17,17,17),nrow=8,byrow=T))
layout(matrix(c(1,1,1,1,5,5,5,5,9,9,9,9,1,1,1,1,5,5,5,5,9,9,9,9,2,2,2,2,6,6,6,6,10,10,10,10,2,2,2,2,6,6,6,6,10,10,10,10,3,3,3,3,7,7,7,7,11,11,11,11,3,3,3,3,7,7,7,7,11,11,11,11,4,4,4,4,8,8,8,8,12,12,12,12,4,4,4,4,8,8,8,8,13,13,13,13),nrow=8,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))

i=is[1]

#responder vs vrld
k=virvars[1]
for(j in js[c(9,10,8,12)])
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
  
  if(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))<0.001)
  {
    pints=sprintf("%.1E",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
  } else {
    pints=round(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))),3)
  }
  if(pchisq(aa1,aa2, lower.tail = F)<0.001)
  {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))				
  } else {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints))		
  }
  
  plot(temp$logdays,temp$var,xlab="Days post-admission",pch=21,bg=adjustcolor("orange",alpha.f = 0.5),ylab=virnams[k],main=j,type='n',xaxt='n',las=2)
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
    legend(2.2,1,legend=c(newpvs,"Non-responder","Responder"),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)		
  } else {
    #		legend("topright",legend=c(newpvs,"Non-responder"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=j)
    legend(2.25,4.75,legend=c(newpvs,"Non-responder","Responder"),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)
  }
}	

#vs TG
k=virvars[12]
for(j in js[c(9,10,8,12)])
{
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
  
  if(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))<0.001)
  {
    pints=sprintf("%.1E",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
  } else {
    pints=round(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))),3)
  }
  if(pchisq(aa1,aa2, lower.tail = F)<0.001)
  {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))				
  } else {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints))		
  }
  
  plot(temp$logdays,temp$gate3,xlab="Days post-admission",yaxt='n',xaxt='n',pch=21,bg=adjustcolor("azure4",alpha.f = 0.5),ylab=paste0("%"),type='n',main=j)
  points(temp$logdays,temp$gate3,bg=temp$col,pch=21,cex=1)
  axis(2,at=log2(c(0.001,0.01,0.10,1.00,10.00,100.00)),labels=c(0.001,0.01,0.1,1,10,100),las=2)
  if(j%in%js[c(10,12)])
  {
    legend("topright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)		
  } else {
    legend("bottomright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)
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

if(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))<0.001)
{
  pints=sprintf("%.1E",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
} else {
  pints=round(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))),3)
}
if(pchisq(aa1,aa2, lower.tail = F)<0.001)
{
  newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))				
} else {
  newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints))		
}

plot(temp$logdays,temp$gate3,xlab="Days post-admission",yaxt='n',xaxt='n',pch=21,bg=adjustcolor("azure4",alpha.f = 0.5),ylab=paste0("%"),type='n',main=j)
points(temp$logdays,temp$gate3,bg=temp$col,pch=21,cex=1)
axis(2,at=log2(c(0.001,0.01,0.10,1.00,10.00,100.00)),labels=c(0.001,0.01,0.1,1,10,100),las=2)
if(j%in%js[c(10,12)])
{
  legend("topright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)		
} else {
  legend("bottomright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)
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

par(bty='n',mai=c(0.32,0.22,0.22,0.22))
##cyto+ memory subsets
grps=c("CD4_N","CD4_CM","CD4_EM","CD4_TEMRA","CD8_N","CD8_CM","CD8_EM","CD8_TEMRA")		
grps.cols2=c("#1F78B4","#33A02C","#FF7F00","#6A3D9A","#1F78B4","#33A02C","#FF7F00","#6A3D9A")	
cytos=readRDS("output/2024.06.21_cytos.obj")

plot(cytos$x,cytos$y,pch=20,cex=0.5,col=grps.cols2[match(cytos$group,grps)],xlab="tSNE-1",ylab="tSNE-2",xaxt='n',yaxt='n')
lines(c(-10,12),c(-25,25),lty=3)
legend("topright",c("Naive","CM","EM","TEMRA"),pch=20,col=grps.cols2[1:4],bty='n',pt.cex=1.5)
text(c(7,16),c(25,25),labels=c("CD4","CD8"))
arrows(-31,-26.5,37,-26.5,length=0.1)
arrows(-31,-26.5,-31,28,length=0.1)
mtext("tSNE-1",side=1,line=0.5,cex=0.75)
mtext("tSNE-2",side=2,line=0,cex=0.75)

#now cytokines
grps.cytos.cols=c("#BEAED4","#8DD3C7","#FFD92F","#666666","#66A61E","#E31A1C","#FF7F00")
grps.cytos=unique(cytos$cytokine2)[c(3,6,1,4,5,7,8)]
plot(cytos$x,cytos$y,pch=20,cex=0.5,col=grps.cytos.cols[match(cytos$cytokine2,grps.cytos)],xlab="tSNE-1",ylab="tSNE-2",xaxt='n',yaxt='n')
lines(c(-10,12),c(-25,25),lty=3)
text(c(7,16),c(25,25),labels=c("CD4","CD8"))
arrows(-31,-26.5,37,-26.5,length=0.1)
arrows(-31,-26.5,-31,28,length=0.1)
legend("topright",grps.cytos[1:4],pch=20,col=grps.cytos.cols[1:4],bty='n',pt.cex=1.5)
legend("bottomright",grps.cytos[5:7],pch=20,col=grps.cytos.cols[5:7],bty='n',pt.cex=1.5)
mtext("tSNE-1",side=1,line=0.5,cex=0.75)
mtext("tSNE-2",side=2,line=0,cex=0.75)

#vs AIM
par(bty='n',mai=c(0.32,0.42,0.32,0.22))

tmp=table(cytos$aim_pos>0,cytos$cytokine2,cytos$type2)
tmpp=tmp
tmpp[,,1][2,]=tmp[,,1][2,]/colSums(tmp[,,1])*100
tmpp[,,1][1,]=100-tmpp[,,1][2,]
tmpp[,,2][2,]=tmp[,,2][2,]/colSums(tmp[,,2])*100
tmpp[,,2][1,]=100-tmpp[,,2][2,]
tmpp4=tmpp[,,1][,match(grps.cytos,colnames(tmpp[,,1]))]
tmpp8=tmpp[,,2][,match(grps.cytos,colnames(tmpp[,,2]))]
tmp4=tmp[,,1][,match(grps.cytos,colnames(tmp[,,1]))]
tmp8=tmp[,,2][,match(grps.cytos,colnames(tmp[,,2]))]

tt=barplot(tmpp4,cex.names=0.75,yaxt='n',ylab="%",xlab="",col=c("azure4","azure2"),xaxt='n')
axis(2,seq(0,100,20),labels=seq(0,100,20),cex.axis=0.75,pos=0,las=2)
legend("topleft",horiz=T,legend=c("AIM-","AIM+"),pch=22,pt.bg=c("azure4","azure2"),bty='n',pt.cex=1.5)
axis(1,tt,labels=grps.cytos,cex.axis=0.5,pos=10,tick=F)
mtext("CD4+ T cells",side=1,line=1.5,cex=0.75)
mtext("%",side=2,line=1.5,cex=0.75)

tt=barplot(tmpp8,cex.names=0.75,yaxt='n',ylab="%",xlab="",col=c("azure4","azure2"),xaxt='n')
axis(2,seq(0,100,20),labels=seq(0,100,20),cex.axis=0.75,pos=0,las=2)
axis(1,tt,labels=grps.cytos,cex.axis=0.5,pos=10,tick=F)
mtext("CD8+ T cells",side=1,line=1.5,cex=0.75)
mtext("%",side=2,line=1.5,cex=0.75)

dev.off()

cytos[["sid"]]=gsub("D0","D1",paste0(cytos[["pid"]],"_",gsub("[-]","",cytos[["visit"]])))
xinfo=readRDS("output/2024.03.25_rajeshids_eventids.obj")
cytos[["event_id"]]=xinfo[[1]]$event_id[match(cytos[["sid"]],gsub("Day","D",xinfo[[1]]$sid))]

cytos[["days"]]=tcell.info$event_date[match(cytos[["event_id"]],tcell.info$event_id)]
cytos[["days_group"]]=NA
cytos[["days_group"]][which(cytos[["days"]]<4)]="1-3"
cytos[["days_group"]][which(cytos[["days"]]>=4)]="4-10"
cytos[["days_group"]][which(cytos[["days"]]>=11)]="11+"
cytos[["days_group"]]=factor(cytos[["days_group"]],c("1-3","4-10","11+"))

cytos[["tg45"]]=tcell.info.acute$tg45[match(cytos$event_id,tcell.info.acute$event_id)]
cytos[["vrld"]]=tcell.info.acute$vrld[match(cytos$event_id,tcell.info.acute$event_id)]

##Supp Figure 8
svg(paste0("output/new_supp_fig8.svg"),width=10.4,height=14.9) 
layout(matrix(c(1,1,1,1,5,5,5,5,9,9,9,9,1,1,1,1,5,5,5,5,9,9,9,9,2,2,2,2,6,6,6,6,9,9,9,9,2,2,2,2,6,6,6,6,9,9,9,9,3,3,3,3,7,7,7,7,9,9,9,9,3,3,3,3,7,7,7,7,9,9,9,9,4,4,4,4,8,8,8,8,9,9,9,9,4,4,4,4,8,8,8,8,9,9,9,9),nrow=8,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))

i=is[1]

k=virvars[1]
for(j in js[c(7,11)])
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
  
  if(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))<0.001)
  {
    pints=sprintf("%.1E",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
  } else {
    pints=round(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))),3)
  }
  if(pchisq(aa1,aa2, lower.tail = F)<0.001)
  {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))				
  } else {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints))		
  }
  
  plot(temp$logdays,temp$var,xlab="Days post-admission",pch=21,bg=adjustcolor("orange",alpha.f = 0.5),ylab=virnams[k],main=j,type='n',xaxt='n',las=2)
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
    legend(2.2,1,legend=c(newpvs,"Non-responder","Responder"),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)		
  } else {
    #		legend("topright",legend=c(newpvs,"Non-responder"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=j)
    legend(2.25,4.75,legend=c(newpvs,"Non-responder","Responder"),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)
  }
}	

#vs TG
k=virvars[12]
for(j in js[c(7,11)])
{
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
  
  if(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))<0.001)
  {
    pints=sprintf("%.1E",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
  } else {
    pints=round(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))),3)
  }
  if(pchisq(aa1,aa2, lower.tail = F)<0.001)
  {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))				
  } else {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints))		
  }
  
  plot(temp$logdays,temp$gate3,xlab="Days post-admission",yaxt='n',xaxt='n',pch=21,bg=adjustcolor("azure4",alpha.f = 0.5),ylab=paste0("%"),type='n',main=j)
  points(temp$logdays,temp$gate3,bg=temp$col,pch=21,cex=1)
  axis(2,at=log2(c(0.001,0.01,0.10,1.00,10.00,100.00)),labels=c(0.001,0.01,0.1,1,10,100),las=2)
  if(j%in%js[c(10,12)])
  {
    legend("topright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)		
  } else {
    legend("bottomright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)
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

plot(cytos$x,cytos$y,pch=20,cex=0.5,col=grps.cytos.cols[match(cytos$cytokine2,grps.cytos)],xlab="",ylab="",xaxt='n',yaxt='n',type='n',main="Viral rpM < median")
points(cytos$x[cytos$vrld==F],cytos$y[cytos$vrld==F],pch=20,cex=0.5,col=grps.cytos.cols[match(cytos$cytokine2[cytos$vrld==F],grps.cytos)])
lines(c(-10,12),c(-25,25),lty=3)
text(c(7,16),c(25,25),labels=c("CD4","CD8"))
arrows(-31,-26.5,37,-26.5,length=0.1)
arrows(-31,-26.5,-31,28,length=0.1)
legend("topright",grps.cytos[1:7],pch=20,col=grps.cytos.cols[1:7],bty='n',pt.cex=1.5)
#legend("bottomright",grps.cytos[5:7],pch=20,col=grps.cytos.cols[5:7],bty='n',pt.cex=1.5)
mtext("tSNE-1",side=1,line=0.5,cex=0.75)
mtext("tSNE-2",side=2,line=0,cex=0.75)

plot(cytos$x,cytos$y,pch=20,cex=0.5,col=grps.cytos.cols[match(cytos$cytokine2,grps.cytos)],xlab="",ylab="",xaxt='n',yaxt='n',type='n',main="Viral rpM > median")
points(cytos$x[cytos$vrld==T],cytos$y[cytos$vrld==T],pch=20,cex=0.5,col=grps.cytos.cols[match(cytos$cytokine2[cytos$vrld==T],grps.cytos)])
lines(c(-10,12),c(-25,25),lty=3)
text(c(7,16),c(25,25),labels=c("CD4","CD8"))
arrows(-31,-26.5,37,-26.5,length=0.1)
arrows(-31,-26.5,-31,28,length=0.1)
#legend("topright",grps.cytos[1:4],pch=20,col=grps.cytos.cols[1:4],bty='n',pt.cex=1.5)
#legend("bottomright",grps.cytos[5:7],pch=20,col=grps.cytos.cols[5:7],bty='n',pt.cex=1.5)
#mtext("tSNE-1",side=1,line=0.5,cex=0.75)
#mtext("tSNE-2",side=2,line=0,cex=0.75)

plot(cytos$x,cytos$y,pch=20,cex=0.5,col=grps.cytos.cols[match(cytos$cytokine2,grps.cytos)],xlab="",ylab="",xaxt='n',yaxt='n',type='n',main="TG123")
points(cytos$x[cytos$tg45==F],cytos$y[cytos$tg45==F],pch=20,cex=0.5,col=grps.cytos.cols[match(cytos$cytokine2[cytos$tg45==F],grps.cytos)])
lines(c(-10,12),c(-25,25),lty=3)
text(c(7,16),c(25,25),labels=c("CD4","CD8"))
arrows(-31,-26.5,37,-26.5,length=0.1)
arrows(-31,-26.5,-31,28,length=0.1)
legend("topright",grps.cytos[1:7],pch=20,col=grps.cytos.cols[1:7],bty='n',pt.cex=1.5)
#legend("bottomright",grps.cytos[5:7],pch=20,col=grps.cytos.cols[5:7],bty='n',pt.cex=1.5)
#mtext("tSNE-1",side=1,line=0.5,cex=0.75)
#mtext("tSNE-2",side=2,line=0,cex=0.75)

plot(cytos$x,cytos$y,pch=20,cex=0.5,col=grps.cytos.cols[match(cytos$cytokine2,grps.cytos)],xlab="",ylab="",xaxt='n',yaxt='n',type='n',main="TG45")
points(cytos$x[cytos$tg45==T],cytos$y[cytos$tg45==T],pch=20,cex=0.5,col=grps.cytos.cols[match(cytos$cytokine2[cytos$tg45==T],grps.cytos)])
lines(c(-10,12),c(-25,25),lty=3)
text(c(7,16),c(25,25),labels=c("CD4","CD8"))
arrows(-31,-26.5,37,-26.5,length=0.1)
arrows(-31,-26.5,-31,28,length=0.1)
#legend("topright",grps.cytos[1:4],pch=20,col=grps.cytos.cols[1:4],bty='n',pt.cex=1.5)
#legend("bottomright",grps.cytos[5:7],pch=20,col=grps.cytos.cols[5:7],bty='n',pt.cex=1.5)
#mtext("tSNE-1",side=1,line=0.5,cex=0.75)
#mtext("tSNE-2",side=2,line=0,cex=0.75)

dev.off()

##Supp Figure 9
svg(paste0("output/new_supp_fig9.svg"),width=10.4,height=14.9) 
layout(matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8,9,9,9,9,7,7,7,7,8,8,8,8,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,10,10,10,10,11,11,11,11,12,12,12,12),nrow=8,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))

for(j in names(cytos)[c(22,23,16,8,9)])
{
  colors=rep("azure3",length(unlist(cytos[which(cytos.names==j)])))
  colors[which(unlist(cytos[which(cytos.names==j)])>=cytos.cuts[j])]="coral3"
  colors[is.na(cytos$group)]="white"
  plot(cytos$x[!is.na(cytos[["type2"]])],cytos$y[!is.na(cytos[["type2"]])],pch=20,cex=0.5,col=colors[!is.na(cytos[["type2"]])],main=j,xlab="",ylab="",xaxt='n',yaxt='n')
  lines(c(-10,12),c(-25,25),lty=3)
  text(c(7,16),c(25,25),labels=c("CD4","CD8"))
  arrows(-31,-25.5,37,-25.5,length=0.1)
  arrows(-31,-25.5,-31,28,length=0.1)
  if(j=="CD40L")
  {
    mtext("tSNE-1",side=1,line=0.5,cex=0.75)
    mtext("tSNE-2",side=2,line=0,cex=0.75)
  }
}

##plot by AIM positivity
aimcol=c("azure4","orange","magenta")[match(cytos[["aim_pos_cat"]],levels(cytos[["aim_pos_cat"]]))]
plot(cytos$x[!is.na(cytos[["type2"]])],cytos$y[!is.na(cytos[["type2"]])],pch=20,cex=0.5,col=aimcol[!is.na(cytos[["type2"]])],main="AIM positivity",xlab="",ylab="",xaxt='n',yaxt='n')
legend("topright",legend=c("None","1",">1"),pch=21,pt.bg=c("azure4","orange","magenta"),bty='n',pt.cex=1.5)
lines(c(-10,12),c(-25,25),lty=3)
text(c(7,16),c(25,25),labels=c("CD4","CD8"))
arrows(-31,-25.5,37,-25.5,length=0.1)
arrows(-31,-25.5,-31,28,length=0.1)
dev.off()

##Supp Figure 7
svg(paste0("output/fig4_supp_omics_ics.svg"),width=10.4,height=14.9) 
layout(matrix(c(1,1,1,2,2,2,3,3,3,4,4,4,5,5,5,6,6,6,7,7,7,8,8,8,9,9,9,10,10,10,11,11,11,12,12,12,13,13,13,14,14,14,15,15,15,16,16,16,17,17,17,18,18,18,19,19,19,20,20,20,21,21,21,22,22,22,23,23,23,24,24,24),nrow=6,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))

i=is[1]
for(j in js[c(7:12)])
{
  for(k in omics.vars[c(51,48,15,29)])
  {
    temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tcell.final.acute[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin)) #get data
    if(omics.type[k]=="cytof")
    {
      tmp=cytof[gsub(".*[-]","",cytof$event_id)=="1",]
      tmpp=cytof.pct[gsub(".*[-]","",cytof$event_id)=="1",]
      temp$bvar=tmpp[match(temp$id,tmp$participant_id),k]
    } else {
      tmp=olink[gsub(".*[-]","",olink$event_id)=="1",]
      temp$bvar=olink[match(temp$id,olink$participant_id),k]
    }
    temp$gate=as.numeric(temp$gate)
    temp$bvar2=log10(temp$bvar+min(temp$bvar[which(temp$bvar>0)]))	
    if(omics.type[k]=="cytof")
    {
      temp$var=cytof.pct[match(temp$event_id,cytof$event_id),k]
    } else {
      temp$var=olink[match(temp$event_id,olink$event_id),k]
    }	
    temp$var2=log10(temp$var+min(temp$var[which(temp$var>0)]))
    temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
    temp$livecd3=livecd3.final[[1]][match(temp$event_id,names(livecd3.final[[1]]))] #get stim live count
    temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
    temp$responder=factor(temp$gate>0)
    temp$logdays=log(as.numeric(temp$days))	
    
    yrng=range(c(temp$var+min(temp$var[which(temp$var>0)]),temp$bvar+min(temp$bvar[which(temp$bvar>0)])),na.rm=T)
    
    if(yrng[1]>0)
    {
      boxplot(var2~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab="",ylab="",type='n',border="white",pch="",col="white",xaxt='n',log='y',las=2)
      par(new=T)
      boxplot(bvar2~responder,temp,xlim=c(0.5,5.5),ylim=log10(yrng),xlab=j,ylab=omics.nams2[k],type='n',border="white",pch="",col="white",yaxt='n',xaxt='n')
      stripchart(bvar2~responder,temp[which(temp$responder==F),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("azure4",alpha.f=0.5),cex=0.5,lwd=0.2)
      stripchart(bvar2~responder,temp[which(temp$responder==T),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("coral2",alpha.f=0.5),cex=0.5,lwd=0.2)
      par(new=T)
      boxplot(bvar2~responder,temp,xlim=c(0.5,5.5),ylim=log10(yrng),xlab="",ylab="",type='n',medcol="black",pch="",col=adjustcolor(c("azure4","coral2"),alpha.f=0.1),xaxt='n',yaxt='n')
      par(new=T)
      boxplot(var2~responder,temp,xlim=c(0.5,5.5),ylim=log10(yrng),xlab=j,ylab="",type='n',border="white",pch="",col="white",yaxt='n',xaxt='n',at=4:5)
      stripchart(var2~responder,temp[which(temp$responder==F),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("azure4",alpha.f=0.5),cex=0.5,at=4:5,lwd=0.2)
      stripchart(var2~responder,temp[which(temp$responder==T),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("coral2",alpha.f=0.5),cex=0.5,at=4:5,lwd=0.2)
      par(new=T)
      boxplot(var2~responder,temp,xlim=c(0.5,5.5),ylim=log10(yrng),xlab="",ylab="",type='n',medcol="black",pch="",col=adjustcolor(c("azure4","coral2"),alpha.f=0.25),xaxt='n',yaxt='n',at=4:5)
      axis(1,at=c(1:2,4:5),labels=rep(c("Non-R","R"),2),cex.axis=0.5)
      newpvs=sprintf("%.1e",c(tcell.bomics.or.pv[[i]][k,j],tcell.omics.or.pv[[i]][k,j]))
    } else {
      yrng=range(c(temp$var,temp$bvar),na.rm=T)
      boxplot(bvar~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab=j,ylab=omics.nams2[k],type='n',border="white",pch="",col="white",xaxt='n',las=2)
      stripchart(bvar~responder,temp[which(temp$responder==F),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("azure4",alpha.f=0.5),cex=0.25,lwd=0.2)
      stripchart(bvar~responder,temp[which(temp$responder==T),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("coral2",alpha.f=0.5),cex=0.25,lwd=0.2)
      par(new=T)
      boxplot(bvar~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab="",ylab="",type='n',medcol="black",pch="",col=adjustcolor(c("azure4","coral2"),alpha.f=0.1),xaxt='n',yaxt='n')
      par(new=T)
      boxplot(var~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab=j,ylab="",type='n',border="white",pch="",col="white",yaxt='n',xaxt='n',at=4:5)
      stripchart(var~responder,temp[which(temp$responder==F),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("azure4",alpha.f=0.5),cex=0.25,at=4:5,lwd=0.2)
      stripchart(var~responder,temp[which(temp$responder==T),],xlim=c(0.5,5.5),add=T,method="jitter",vertical=T,pch=21,bg=adjustcolor("coral2",alpha.f=0.5),cex=0.25,at=4:5,lwd=0.2)
      par(new=T)
      boxplot(var~responder,temp,xlim=c(0.5,5.5),ylim=yrng,xlab="",ylab="",type='n',medcol="black",pch="",col=adjustcolor(c("azure4","coral2"),alpha.f=0.25),xaxt='n',yaxt='n',at=4:5)
      axis(1,at=c(1:2,4:5),labels=rep(c("Non-R","R"),2),cex.axis=0.5)
      newpvs=sprintf("%.1e",c(tcell.bomics.or.pv[[i]][k,j],tcell.omics.or.pv[[i]][k,j]))
    }
    abline(v=3,lty=3)
    axis(3,at=c(1.5,4.5),tick=F,labels=paste0(c("Baseline-P=","P="),newpvs),line=-0.5)
    
  }
}	

dev.off()

newnewcols=c("black","azure3","orange","firebrick3")

##Figure 6
svg(paste0("output/new_fig6.svg"),width=10.4,height=14.9) 
layout(matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,1,1,1,1,2,2,2,2,4,4,4,4,5,5,5,5,8,8,8,8,11,11,11,11,5,5,5,5,8,8,8,8,11,11,11,11,6,6,6,6,9,9,9,9,12,12,12,12,6,6,6,6,9,9,9,9,12,12,12,12,7,7,7,7,10,10,10,10,13,13,13,13,7,7,7,7,10,10,10,10,13,13,13,13),nrow=8,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))

plot(aimics$x,aimics$y,pch=20,cex=0.5,col=grps.cols2[match(aimics$group,grps)],xaxt='n',yaxt='n',xlab="",ylab="")
lines(c(-30,1),c(-30,18))
legend("bottomleft",c("Naive","CM","EM","TEMRA"),pch=20,col=grps.cols2[1:4],bty='n',pt.cex=1.5)
text(c(-4,6),c(18,18),labels=c("CD8","CD4"))
arrows(-32,-31.5,29,-31.5,length=0.1)
arrows(-32,-31.5,-32,18,length=0.1)
mtext("tSNE-1",side=1,line=0.25,cex=0.75)
mtext("tSNE-2",side=2,line=0.25,cex=0.75)

plot(aimics$x,aimics$y,pch=20,cex=0.5,col=grps.aimics.cols[match(aimics$aimics_group,grps.aimics)],xlab="",ylab="",xaxt='n',yaxt='n')
legend("bottomleft",legend=grps.aimics,pch=20,col=grps.aimics.cols,bty='n',pt.cex=1.5)
lines(c(-30,1),c(-30,18))
text(c(-4,6),c(18,18),labels=c("CD8","CD4"))
arrows(-32,-31.5,29,-31.5,length=0.1)
arrows(-32,-31.5,-32,18,length=0.1)

##barplot
for(k in aim.vars[2:3])
{
  tmp=table(paste0(aimics$type2,"-",aimics$aimics_group),aimics[[k]])[1:6,]
  tt=barplot(t(tmp)/colSums(tmp)*100,beside=T,las=2,ylab="",xlab="",col=c("azure4","azure2"),xaxt='n')
  axis(1,colMeans(tt),labels=rep(grps.aimics,2),cex.axis=0.75,pos=6,tick=F)
  legend("topright",horiz=F,legend=aims.vars.labs[[k]],pch=22,pt.bg=c("azure4","azure2"),bty='n',pt.cex=1.5)
  mtext("%",side=2,line=2.5,cex=0.75)
  axis(1,c(mean(tt[1:6]),mean(tt[7:12])),labels=c("CD4+ T cells","CD8+ T cells"),cex.axis=1,pos=-2,tick=F)
  
  if(chisq.test(t(tmp)/colSums(tmp)*100)$p.value<0.001)
  {
    newpv=sprintf("%.1E",chisq.test(t(tmp)/colSums(tmp)*100)$p.value)
  } else {
    newpv=round(chisq.test(t(tmp)/colSums(tmp)*100)$p.value,3)
  }
  legend("right",legend=paste0("P = ",newpv),bty='n')
}

k=virvars[1]
l=js[5]

for(j in js[c(8,10,12)])
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
  
  temp$responder=as.factor(temp$gate>0) #make responder
  
  temp$days_binned=tcell.info.acute$event_date_binned[match(temp$event_id,tcell.info.acute$event_id)]
  temp$days_symptom_binned=tcell.info.acute$event_date_symptom_binned[match(temp$event_id,tcell.info.acute$event_id)]
  
  temp$logdays=log(temp$days) #needed for lme
  temp$logdays_symptom=log(temp$days_symptom) #needed for lme
  temp=temp[!is.na(temp$var),]	#remove missing
  
  lgth=length(levels(temp$responder))
  
  temp$aim=(tcell.final.acute[[i]][match(temp$event_id,tcell.info.acute$event_id),l]>0)+0
  temp$newvar=factor(paste0(temp$responder,"_",temp$aim))
  
  lgth=length(levels(temp$newvar))
  temp$col=newnewcols[match(temp$newvar,levels(temp$newvar))]
  
  
  m=gam(var~newvar+s(logdays,bs="cr",by=newvar),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$newvar),],random =list(id=~1),family="gaussian")
  a1=anova(m)
  aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
  aa2=sum(a1$s.table[-c(1), 2])
  
  pints=as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))
  pints2=as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))
  for(z in 1:length(pints))
  {
    if(pints[z]<0.001)
    {
      pints2[z]=sprintf("%.1e",pints[z])
    } else {
      pints2[z]=round(pints[z],3)
    }
  }
  if(pchisq(aa1,aa2, lower.tail = F)<0.001)
  {
    newpvs=paste0(c("Slope",paste0(c("R-","NR+","R+"),j)),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints2))				
  } else {
    newpvs=paste0(c("Slope",paste0(c("R-","NR+","R+"),j)),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints2))		
  }
  
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
  
  temp$aim=(tcell.final.acute[[i]][match(temp$event_id,tcell.info.acute$event_id),l]>0)+0
  temp$newvar=factor(paste0(temp$var,"_",temp$aim))
  temp$newvar=relevel(temp$newvar,"FALSE_1")
  
  lgth=length(levels(temp$newvar))
  temp$col=newnewcols[match(temp$newvar,levels(temp$newvar))]
  m=gam(gate3~newvar+s(logdays,bs="cr",by=newvar),data=temp[!is.na(temp$logdays)&!is.na(temp$newvar)&!is.na(temp$gate3),],random =list(id=~1))
  a1=anova(m)
  aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
  aa2=sum(a1$s.table[-c(1), 2])
  
  pints=as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))
  pints2=as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))
  for(z in 1:length(pints))
  {
    if(pints[z]<0.001)
    {
      pints2[z]=sprintf("%.1e",pints[z])
    } else {
      pints2[z]=round(pints[z],3)
    }
  }
  if(pchisq(aa1,aa2, lower.tail = F)<0.001)
  {
    newpvs=paste0(c("Slope","TG123+NR","TG45+NR","TG45+R"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints2))				
  } else {
    newpvs=paste0(c("Slope","TG123+NR","TG45+NR","TG45+R"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints2))		
  }
  
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
  temp$ics=tcell.final.acute[[i]][match(temp$event_id,tcell.info.acute$event_id),l]
  
  temp=temp[which(temp$ics>0),]
  temp$ics=temp$ics>median(temp$ics)
  temp$newvar=factor(paste0(temp$responder,"_",temp$ics))
  tmp=table(temp$var,temp$newvar)
  tttt=colSums(tmp)
  tmp2=tmp[2,]/colSums(tmp)*100
  m=glm(var~relevel(temp$newvar,"TRUE_FALSE"),temp,family="binomial")
  pints=as.numeric(unlist(summary(m)$coefficients[-1,4]))
  pints2=as.numeric(unlist(summary(m)$coefficients[-1,4]))
  for(z in 1:length(pints))
  {
    if(pints[z]<0.001)
    {
      pints2[z]=sprintf("%.1e",pints[z])
    } else {
      pints2[z]=round(pints[z],3)
    }
  } 
  ymax=100
  tmp=barplot(tmp2,main="",ylab="TG45 (%)",xaxt='n',col=newnewcols[c(1,3,2,4)],las=2,xlab=paste0(j," & ",l),ylim=c(0,ymax))
  axis(1,tmp,c("NR+ICSlow","NR+ICShigh","R+ICSlow","R+ICShigh"),tick=F,las=1,cex.axis=0.8)
  text(tmp,rep(5,length(tmp)),labels=unlist(tttt),col="white")
  text(tmp[-3],rep(70,length(tmp[-3])),labels=paste0("P=",pints2))
}

dev.off()

#Supp Figure 12
svg(paste0("output/new_supp_fig12.svg"),width=10.4,height=14.9) 
#layout(matrix(c(1,1,1,1,2,2,3,3,4,4,5,5,1,1,1,1,2,2,3,3,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8,7,7,7,7,7,7,8,8,8,8,8,8,9,9,10,10,10,10,11,11,12,12,12,12,9,9,10,10,10,10,11,11,12,12,12,12,13,13,13,13,13,13,14,14,14,15,15,15,13,13,13,13,13,13,16,16,16,17,17,17),nrow=8,byrow=T))
layout(matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8,9,9,9,9,7,7,7,7,8,8,8,8,10,10,10,10,rep(11,24)),nrow=8,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))

tmpj=c("CD137","CD69","CD25","CD40L","OX40","IFNg","IL-2","TNFa")
names(tmpj)=c("CD137","CD69","CD25","CD40L","OX40","IFN-g","IL-2","TNF-a")
for(j in c("CD137","CD69","CD25","CD40L","OX40","IFN-g","IL-2","TNF-a"))
{
  colors=rep("azure3",length(unlist(aimics[which(aimics.names==j)])))
  colors[which(unlist(aimics[which(aimics.names==j)])>=aimics.cuts[j])]="coral3"
  colors[is.na(aimics$group)]="white"
  plot(aimics$x,aimics$y,pch=20,cex=0.5,col=colors,main=tmpj[j],xlab="",ylab="",xaxt='n',yaxt='n')
  lines(c(-30,1),c(-30,18))
  text(c(-4,6),c(18,18),labels=c("CD8","CD4"))
  arrows(-32,-31.5,29,-31.5,length=0.1)
  arrows(-32,-31.5,-32,18,length=0.1)
}

tmp=table(aimics$group,aimics$aimics_group)
tt=barplot(t(tmp[c(3,1,2,4),]/rowSums(tmp[c(3,1,2,4),])*100),beside=T,col=grps.aimics.cols,xaxt='n',yaxt='n')
axis(2,seq(0,100,20),labels=seq(0,100,20),cex.axis=0.75,pos=0,las=2)
axis(1,colMeans(tt),labels=c("Naive","CM","EM","TEMRA"),cex.axis=0.9,pos=7.5,tick=F)
mtext("CD4+ T cells",side=1,line=1.5,cex=0.75)
mtext("%",side=2,line=2.5,cex=0.75)

tt=barplot(t(tmp[c(3,1,2,4)+4,]/rowSums(tmp[c(3,1,2,4)+4,])*100),beside=T,col=grps.aimics.cols,xaxt='n',yaxt='n')
axis(2,seq(0,100,20),labels=seq(0,100,20),cex.axis=0.75,pos=0,las=2)
axis(1,colMeans(tt),labels=c("Naive","CM","EM","TEMRA"),cex.axis=0.9,pos=7.5,tick=F)
mtext("CD8+ T cells",side=1,line=1.5,cex=0.75)
mtext("%",side=2,line=2.5,cex=0.75)
legend("topright",horiz=F,legend=grps.aimics,pch=22,pt.bg=grps.aimics.cols,bty='n',pt.cex=1.5)

dev.off()

#Supp Figure 13
svg(paste0("output/new_supp_fig13.svg"),width=10.4,height=14.9) 
layout(matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8,9,9,9,9,7,7,7,7,8,8,8,8,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,10,10,10,10,11,11,11,11,12,12,12,12),nrow=8,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))

k=virvars[1]
for(j in js[c(9,11)])
{
  for(l in js[c(2)])
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
    
    temp$responder=as.factor(temp$gate>0) #make responder
    
    temp$days_binned=tcell.info.acute$event_date_binned[match(temp$event_id,tcell.info.acute$event_id)]
    temp$days_symptom_binned=tcell.info.acute$event_date_symptom_binned[match(temp$event_id,tcell.info.acute$event_id)]
    
    temp$logdays=log(temp$days) #needed for lme
    temp$logdays_symptom=log(temp$days_symptom) #needed for lme
    temp=temp[!is.na(temp$var),]	#remove missing
    
    lgth=length(levels(temp$responder))
    
    temp$aim=(tcell.final.acute[[i]][match(temp$event_id,tcell.info.acute$event_id),l]>0)+0
    temp$newvar=factor(paste0(temp$responder,"_",temp$aim))
    
    lgth=length(levels(temp$newvar))
    temp$col=newnewcols[match(temp$newvar,levels(temp$newvar))]
    
    
    m=gam(var~newvar+s(logdays,bs="cr",by=newvar),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$newvar),],random =list(id=~1),family="gaussian")
    a1=anova(m)
    aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
    aa2=sum(a1$s.table[-c(1), 2])
    
    pints=as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))
    pints2=as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))
    for(z in 1:length(pints))
    {
      if(pints[z]<0.001)
      {
        pints2[z]=sprintf("%.1e",pints[z])
      } else {
        pints2[z]=round(pints[z],3)
      }
    }
    if(pchisq(aa1,aa2, lower.tail = F)<0.001)
    {
      newpvs=paste0(c("Slope",paste0(c("R-","NR+","R+"),j)),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints2))				
    } else {
      newpvs=paste0(c("Slope",paste0(c("R-","NR+","R+"),j)),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints2))		
    }		
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
}

l=js[5]
for(j in js[c(8:11)])
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
  
  temp$responder=as.factor(temp$gate>0) #make responder
  
  temp$days_binned=tcell.info.acute$event_date_binned[match(temp$event_id,tcell.info.acute$event_id)]
  temp$days_symptom_binned=tcell.info.acute$event_date_symptom_binned[match(temp$event_id,tcell.info.acute$event_id)]
  
  temp$logdays=log(temp$days) #needed for lme
  temp$logdays_symptom=log(temp$days_symptom) #needed for lme
  temp=temp[!is.na(temp$var),]	#remove missing
  
  lgth=length(levels(temp$responder))
  
  temp$aim=(tcell.final.acute[[i]][match(temp$event_id,tcell.info.acute$event_id),l]>0)+0
  temp$newvar=factor(paste0(temp$responder,"_",temp$aim))
  
  lgth=length(levels(temp$newvar))
  temp$col=newnewcols[match(temp$newvar,levels(temp$newvar))]
  
  
  m=gam(var~newvar+s(logdays,bs="cr",by=newvar),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$newvar),],random =list(id=~1),family="gaussian")
  a1=anova(m)
  aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
  aa2=sum(a1$s.table[-c(1), 2])
  
  pints=as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))
  pints2=as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))
  for(z in 1:length(pints))
  {
    if(pints[z]<0.001)
    {
      pints2[z]=sprintf("%.1e",pints[z])
    } else {
      pints2[z]=round(pints[z],3)
    }
  }
  if(pchisq(aa1,aa2, lower.tail = F)<0.001)
  {
    newpvs=paste0(c("Slope",paste0(c("R-","NR+","R+"),j)),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints2))				
  } else {
    newpvs=paste0(c("Slope",paste0(c("R-","NR+","R+"),j)),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints2))		
  }	
  
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

l=js[6]
for(j in js[c(8,11,12)])
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
  
  temp$responder=as.factor(temp$gate>0) #make responder
  
  temp$days_binned=tcell.info.acute$event_date_binned[match(temp$event_id,tcell.info.acute$event_id)]
  temp$days_symptom_binned=tcell.info.acute$event_date_symptom_binned[match(temp$event_id,tcell.info.acute$event_id)]
  
  temp$logdays=log(temp$days) #needed for lme
  temp$logdays_symptom=log(temp$days_symptom) #needed for lme
  temp=temp[!is.na(temp$var),]	#remove missing
  
  lgth=length(levels(temp$responder))
  
  temp$aim=(tcell.final.acute[[i]][match(temp$event_id,tcell.info.acute$event_id),l]>0)+0
  temp$newvar=factor(paste0(temp$responder,"_",temp$aim))
  
  lgth=length(levels(temp$newvar))
  temp$col=newnewcols[match(temp$newvar,levels(temp$newvar))]
  
  
  m=gam(var~newvar+s(logdays,bs="cr",by=newvar),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$newvar),],random =list(id=~1),family="gaussian")
  a1=anova(m)
  aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
  aa2=sum(a1$s.table[-c(1), 2])
  
  pints=as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))
  pints2=as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))
  for(z in 1:length(pints))
  {
    if(pints[z]<0.001)
    {
      pints2[z]=sprintf("%.1e",pints[z])
    } else {
      pints2[z]=round(pints[z],3)
    }
  }
  if(pchisq(aa1,aa2, lower.tail = F)<0.001)
  {
    newpvs=paste0(c("Slope",paste0(c("R-","NR+","R+"),j)),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints2))				
  } else {
    newpvs=paste0(c("Slope",paste0(c("R-","NR+","R+"),j)),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints2))		
  }		
  
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

dev.off()

#Supp Figure 14
svg(paste0("output/new_supp_fig14.svg"),width=10.4,height=14.9) 
layout(matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8,9,9,9,9,7,7,7,7,8,8,8,8,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,10,10,10,10,11,11,11,11,12,12,12,12),nrow=8,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))

k=virvars[12]
for(l in js[c(2,5)])
{
  for(j in js[c(8,10,12)])
  {
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
    
    temp$aim=(tcell.final.acute[[i]][match(temp$event_id,tcell.info.acute$event_id),l]>0)+0
    temp$newvar=factor(paste0(temp$var,"_",temp$aim))
    temp$newvar=relevel(temp$newvar,"FALSE_1")
    
    lgth=length(levels(temp$newvar))
    temp$col=newnewcols[match(temp$newvar,levels(temp$newvar))]
    m=gam(gate3~newvar+s(logdays,bs="cr",by=newvar),data=temp[!is.na(temp$logdays)&!is.na(temp$newvar)&!is.na(temp$gate3),],random =list(id=~1))
    a1=anova(m)
    aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
    aa2=sum(a1$s.table[-c(1), 2])
    
    pints=as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))
    pints2=as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))
    for(z in 1:length(pints))
    {
      if(pints[z]<0.001)
      {
        pints2[z]=sprintf("%.1e",pints[z])
      } else {
        pints2[z]=round(pints[z],3)
      }
    }
    if(pchisq(aa1,aa2, lower.tail = F)<0.001)
    {
      newpvs=paste0(c("Slope","TG123+NR","TG45+NR","TG45+R"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints2))				
    } else {
      newpvs=paste0(c("Slope","TG123+NR","TG45+NR","TG45+R"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints2))		
    }
    
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
}

k=virvars[12]
j=js[9]

for(l in js[c(1,2,5)])
{
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
  
  temp$aim=(tcell.final.acute[[i]][match(temp$event_id,tcell.info.acute$event_id),l]>0)+0
  temp$newvar=factor(paste0(temp$var,"_",temp$aim))
  temp$newvar=relevel(temp$newvar,"FALSE_1")
  
  lgth=length(levels(temp$newvar))
  temp$col=newnewcols[match(temp$newvar,levels(temp$newvar))]
  m=gam(gate3~newvar+s(logdays,bs="cr",by=newvar),data=temp[!is.na(temp$logdays)&!is.na(temp$newvar)&!is.na(temp$gate3),],random =list(id=~1))
  a1=anova(m)
  aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
  aa2=sum(a1$s.table[-c(1), 2])
  
  pints=as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))
  pints2=as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))
  for(z in 1:length(pints))
  {
    if(pints[z]<0.001)
    {
      pints2[z]=sprintf("%.1e",pints[z])
    } else {
      pints2[z]=round(pints[z],3)
    }
  }
  if(pchisq(aa1,aa2, lower.tail = F)<0.001)
  {
    newpvs=paste0(c("Slope","TG123+NR","TG45+NR","TG45+R"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints2))				
  } else {
    newpvs=paste0(c("Slope","TG123+NR","TG45+NR","TG45+R"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints2))		
  }
  
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

dev.off()

#Supp Figure 15
svg(paste0("output/new_supp_fig15.svg"),width=10.4,height=14.9) 
layout(matrix(c(1,1,1,1,2,2,2,2,3,3,3,3,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,8,9,9,9,9,7,7,7,7,8,8,8,8,9,9,9,9,10,10,10,10,11,11,11,11,12,12,12,12,10,10,10,10,11,11,11,11,12,12,12,12),nrow=8,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))

k=virvars[12]
j=js[1]

for(l in js[c(7,9,11)])
{
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
  temp$ics=tcell.final.acute[[i]][match(temp$event_id,tcell.info.acute$event_id),l]
  
  temp=temp[which(temp$ics>0),]
  temp$ics=temp$ics>median(temp$ics)
  temp$newvar=factor(paste0(temp$responder,"_",temp$ics))
  tmp=table(temp$var,temp$newvar)
  tttt=colSums(tmp)
  tmp2=tmp[2,]/colSums(tmp)*100
  m=glm(var~relevel(temp$newvar,"TRUE_FALSE"),temp,family="binomial")
  pints=as.numeric(unlist(summary(m)$coefficients[-1,4]))
  pints2=as.numeric(unlist(summary(m)$coefficients[-1,4]))
  for(z in 1:length(pints))
  {
    if(pints[z]<0.001)
    {
      pints2[z]=sprintf("%.1e",pints[z])
    } else {
      pints2[z]=round(pints[z],3)
    }
  }  
  ymax=100
  tmp=barplot(tmp2,main="",ylab="TG45 (%)",xaxt='n',col=newnewcols[c(1,3,2,4)],las=2,xlab=paste0(j," & ",l),ylim=c(0,ymax))
  axis(1,tmp,c("NR+ICSlow","NR+ICShigh","R+ICSlow","R+ICShigh"),tick=F,las=1,cex.axis=0.8)
  text(tmp,rep(5,length(tmp)),labels=unlist(tttt),col="white")
  text(tmp[-3],rep(70,length(tmp[-3])),labels=paste0("P=",pints2))
}

k=virvars[12]
j=js[2]

for(l in js[c(8,10,12)])
{
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
  temp$ics=tcell.final.acute[[i]][match(temp$event_id,tcell.info.acute$event_id),l]
  
  temp=temp[which(temp$ics>0),]
  temp$ics=temp$ics>median(temp$ics)
  temp$newvar=factor(paste0(temp$responder,"_",temp$ics))
  tmp=table(temp$var,temp$newvar)
  tttt=colSums(tmp)
  tmp2=tmp[2,]/colSums(tmp)*100
  m=glm(var~relevel(temp$newvar,"TRUE_FALSE"),temp,family="binomial")
  pints=as.numeric(unlist(summary(m)$coefficients[-1,4]))
  pints2=as.numeric(unlist(summary(m)$coefficients[-1,4]))
  for(z in 1:length(pints))
  {
    if(pints[z]<0.001)
    {
      pints2[z]=sprintf("%.1e",pints[z])
    } else {
      pints2[z]=round(pints[z],3)
    }
  } 
  ymax=100
  tmp=barplot(tmp2,main="",ylab="TG45 (%)",xaxt='n',col=newnewcols[c(1,3,2,4)],las=2,xlab=paste0(j," & ",l),ylim=c(0,ymax))
  axis(1,tmp,c("NR+ICSlow","NR+ICShigh","R+ICSlow","R+ICShigh"),tick=F,las=1,cex.axis=0.8)
  text(tmp,rep(5,length(tmp)),labels=unlist(tttt),col="white")
  text(tmp[-3],rep(70,length(tmp[-3])),labels=paste0("P=",pints2))
}

j=js[2]

for(l in js[c(7,9,11)])
{
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
  temp$ics=tcell.final.acute[[i]][match(temp$event_id,tcell.info.acute$event_id),l]
  
  temp=temp[which(temp$ics>0),]
  temp$ics=temp$ics>median(temp$ics)
  temp$newvar=factor(paste0(temp$responder,"_",temp$ics))
  tmp=table(temp$var,temp$newvar)
  tttt=colSums(tmp)
  tmp2=tmp[2,]/colSums(tmp)*100
  m=glm(var~relevel(temp$newvar,"TRUE_FALSE"),temp,family="binomial")
  pints=as.numeric(unlist(summary(m)$coefficients[-1,4]))
  pints2=as.numeric(unlist(summary(m)$coefficients[-1,4]))
  for(z in 1:length(pints))
  {
    if(pints[z]<0.001)
    {
      pints2[z]=sprintf("%.1e",pints[z])
    } else {
      pints2[z]=round(pints[z],3)
    }
  } 
  ymax=100
  tmp=barplot(tmp2,main="",ylab="TG45 (%)",xaxt='n',col=newnewcols[c(1,3,2,4)],las=2,xlab=paste0(j," & ",l),ylim=c(0,ymax))
  axis(1,tmp,c("NR+ICSlow","NR+ICShigh","R+ICSlow","R+ICShigh"),tick=F,las=1,cex.axis=0.8)
  text(tmp,rep(5,length(tmp)),labels=unlist(tttt),col="white")
  text(tmp[-3],rep(70,length(tmp[-3])),labels=paste0("P=",pints2))
}

j=js[5]

for(l in js[c(8,10,12)])
{
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
  temp$ics=tcell.final.acute[[i]][match(temp$event_id,tcell.info.acute$event_id),l]
  
  temp=temp[which(temp$ics>0),]
  temp$ics=temp$ics>median(temp$ics)
  temp$newvar=factor(paste0(temp$responder,"_",temp$ics))
  tmp=table(temp$var,temp$newvar)
  tttt=colSums(tmp)
  tmp2=tmp[2,]/colSums(tmp)*100
  m=glm(var~relevel(temp$newvar,"TRUE_FALSE"),temp,family="binomial")
  pints=as.numeric(unlist(summary(m)$coefficients[-1,4]))
  pints2=as.numeric(unlist(summary(m)$coefficients[-1,4]))
  for(z in 1:length(pints))
  {
    if(pints[z]<0.001)
    {
      pints2[z]=sprintf("%.1e",pints[z])
    } else {
      pints2[z]=round(pints[z],3)
    }
  } 
  ymax=100
  tmp=barplot(tmp2,main="",ylab="TG45 (%)",xaxt='n',col=newnewcols[c(1,3,2,4)],las=2,xlab=paste0(j," & ",l),ylim=c(0,ymax))
  axis(1,tmp,c("NR+ICSlow","NR+ICShigh","R+ICSlow","R+ICShigh"),tick=F,las=1,cex.axis=0.8)
  text(tmp,rep(5,length(tmp)),labels=unlist(tttt),col="white")
  text(tmp[-3],rep(70,length(tmp[-3])),labels=paste0("P=",pints2))
}

dev.off()

##DMSO figure 7
svg(paste0("output/new_fig7.svg"),width=10.4,height=14.9) 
layout(matrix(c(1,1,1,1,1,1,2,2,2,2,2,2,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,4,4,5,5,6,6,7,7,3,3,3,3,4,4,5,5,8,8,8,8,9,9,9,9,9,9,10,10,10,10,10,10,9,9,9,9,9,9,10,10,10,10,10,10,11,11,12,12,12,12,13,13,14,14,14,14,11,11,12,12,12,12,13,13,14,14,14,14),nrow=8,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))


#Spike, CD8 TNFa
i=is[1]
k=virvars[12]
j=js[12]
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

#DMSO
temp$gate=ncs[match(temp$event_id,rownames(ncs)),j]	
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

if(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))<0.001)
{
  pints=sprintf("%.1E",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
} else {
  pints=round(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))),3)
}
#newpvs=paste0(c("Slope",levels(temp$var)[-1]),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))		
if(pchisq(aa1,aa2, lower.tail = F)<0.001)
{
  newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))				
} else {
  newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints))		
}		

plot(temp$logdays,temp$gate3,xlab="Days post-admission",yaxt='n',xaxt='n',pch=21,bg=adjustcolor("azure4",alpha.f = 0.5),ylab=paste0(j," (%)"),type='n',main="DMSO")
points(temp$logdays,temp$gate3,bg=temp$col,pch=21,cex=1)
axis(2,at=log2(c(0.001,0.01,0.10,1.00,10.00,100.00)),labels=c(0.001,0.01,0.1,1,10,100),las=2)
if(j%in%js[c(10,12)])
{
  #legend("topright",legend=c(newpvs,"FALSE"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=virnams[k])		
  legend("topright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)		
} else {
  #legend("bottomright",legend=c(newpvs,"FALSE"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=virnams[k])
  legend("bottomright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)
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

#raw data
temp$gate=tcell.rawdata[[i]][match(temp$event_id,rownames(tcell.rawdata[[i]])),j]
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

if(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))<0.001)
{
  pints=sprintf("%.1E",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
} else {
  pints=round(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))),3)
}
#newpvs=paste0(c("Slope",levels(temp$var)[-1]),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))		
if(pchisq(aa1,aa2, lower.tail = F)<0.001)
{
  newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))				
} else {
  newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints))		
}		

plot(temp$logdays,temp$gate3,xlab="Days post-admission",yaxt='n',xaxt='n',pch=21,bg=adjustcolor("azure4",alpha.f = 0.5),ylab=paste0(j," (%)"),type='n',main="Spike-uncorrected")
points(temp$logdays,temp$gate3,bg=temp$col,pch=21,cex=1)
axis(2,at=log2(c(0.001,0.01,0.10,1.00,10.00,100.00)),labels=c(0.001,0.01,0.1,1,10,100),las=2)
if(j%in%js[c(10,12)])
{
  #legend("topright",legend=c(newpvs,"FALSE"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=virnams[k])		
  legend("topright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)		
} else {
  #legend("bottomright",legend=c(newpvs,"FALSE"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=virnams[k])
  legend("bottomright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)
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

#PCA and clusters
i=is[1]
#tmp=ncs.final.acute[apply(ncs.final.acute,1,function(x){sum(is.na(x))/length(x)*100})==0,]
#tmpp=prcomp(tmp[,js],scale=T)
#tmpc=rowSums(tmp[,js[1:12]]>0)

#tmpp2=prcomp(tmp[tmpp$x[,1]<8,js],scale=T)

ii=cut(tmpc,breaks=seq(min(tmpc),max(tmpc),len=100),include.lowest=T)
tcolors=colorRampPalette(spctrl[rev(c(2:5,8:10))])(99)[ii]
plot(tmpp2$x[,1],tmpp2$x[,2],xlab=paste0("PC1 (",round(summary(tmpp2)[[6]][2,1]*100,1),"%)"),ylab=paste0("PC2 (",round(summary(tmpp2)[[6]][2,2]*100,1),"%)"),xaxt='n',yaxt='n',pch=21,bg="azure3",cex=1)
axis(1,at=seq(-5,15,5),labels=seq(-5,15,5))
axis(2,at=seq(-15,10,5),labels=seq(-15,10,5),las=2)
abline(h=0,v=0,lty=3)
#legend("bottomleft",title="No. +",legend=seq(0,12,2),pch=21,pt.bg=tcolors[match(seq(0,12,2),tmpc)],bty='n')

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
plot(tmpp2$x[ncsk2$cluster==z,1],tmpp2$x[ncsk2$cluster==z,2],xlab="PC1",ylab="PC2",xaxt='n',yaxt='n',pch=21,bg="orange",cex=0.8,xlim=c(-2,10),ylim=c(-15,6),main=paste0("Low-responders"))
axis(1,at=seq(-5,15,5),labels=seq(-5,15,5),pos=-15)
axis(2,at=seq(-15,10,5),labels=seq(-15,10,5),las=2,pos=-2)
z=2
plot(tmpp2$x[ncsk2$cluster==z,1],tmpp2$x[ncsk2$cluster==z,2],xlab="PC1",ylab="PC2",xaxt='n',yaxt='n',pch=21,bg="magenta",cex=0.8,xlim=c(-2,10),ylim=c(-15,6),main=paste0("High-responders"))
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
legend("topright",legend=c("Low-responder","High-responder"),pch=22,pt.bg=c("orange","magenta"),pt.cex=1.5,bty='n')

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
  
  if(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))<0.001)
  {
    pints=sprintf("%.1E",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
  } else {
    pints=round(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))),3)
  }
  #newpvs=paste0(c("Slope",levels(temp$var)[-1]),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))		
  if(pchisq(aa1,aa2, lower.tail = F)<0.001)
  {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))				
  } else {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints))		
  }
  
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
    legend("bottomright",legend=c(newpvs,"Low-responder","High-responder"),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",c("orange","magenta")),pt.cex=1.5)		
  } else {
    legend("right",legend=c(newpvs,"Low-responder","High-responder"),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",c("orange","magenta")),pt.cex=1.5)
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
  if(summary(m)$coefficients[2,4]<0.001)
  {
    newpvs=sprintf("%.1E",summary(m)$coefficients[2,4])
  } else {
    newpvs=round(summary(m)$coefficients[2,4],3)
  }
  ymax=100
  tmp=barplot(tmp2,main="",ylab="High-responder (%)",xaxt='n',col=unlist(lapply(vircols[[k]],function(x){rep(x,1)})),las=2,xlab=virnams[k],ylim=c(0,ymax))
  axis(1,tmp,virlvls[[k]],tick=F,las=1)
  legend("topright",legend=paste0("P=",newpvs),bty='n')
  text(tmp,rep(2.5,length(tmp)),labels=unlist(tttt))
  
  
  m=gam(ncs_responder_c2~var+s(logdays,bs="cr",by=var),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$ncs_responder_c2),],random =list(id=~1),family="binomial")
  a1=anova(m)
  aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
  aa2=sum(a1$s.table[-c(1), 2])
  
  if(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))<0.001)
  {
    pints=sprintf("%.1E",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
  } else {
    pints=round(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))),3)
  }
  #newpvs=paste0(c("Slope",levels(temp$var)[-1]),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))		
  if(pchisq(aa1,aa2, lower.tail = F)<0.001)
  {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))				
  } else {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints))		
  }
  
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
  tmp=barplot(unlist(ttt),main="",ylab="High-responder (%)",xaxt='n',col=unlist(lapply(vircols[[k]],function(x){rep(x,3)})),las=2,xlab="Days post-admission",ylim=c(0,ymax))
  axis(1,tmp,rep(tnams,lgth),tick=F,las=1)
  #legend("topright",legend=c(newpvs,"TRUE","FALSE"),bty='n',pch=c(NA,rep(22,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=virnams[k])
  legend("topright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(22,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)		
  text(tmp,rep(2.5,length(tmp)),labels=unlist(tttt))
}

dev.off()

##DMSO supp figure 16
svg(paste0("output/new_supp_fig16.svg"),width=10.4,height=14.9) 
layout(matrix(c(1,1,1,1,1,1,2,2,2,2,2,2,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5,5,5,6,6,6,6,6,6,5,5,5,5,5,5,6,6,6,6,6,6,7,7,7,7,7,7,8,8,8,8,8,8,7,7,7,7,7,7,8,8,8,8,8,8),nrow=8,byrow=T))
par(bty='n',mai=c(0.62,0.62,0.22,0.22))


for(j in js[c(9,10,8,11)])
{
  i=is[1]
  k=virvars[12]
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
  
  #DMSO
  temp$gate=ncs[match(temp$event_id,rownames(ncs)),j]	
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
  
  if(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))<0.001)
  {
    pints=sprintf("%.1E",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
  } else {
    pints=round(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))),3)
  }
  #newpvs=paste0(c("Slope",levels(temp$var)[-1]),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))		
  if(pchisq(aa1,aa2, lower.tail = F)<0.001)
  {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))				
  } else {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints))		
  }
  
  plot(temp$logdays,temp$gate3,xlab="Days post-admission",yaxt='n',xaxt='n',pch=21,bg=adjustcolor("azure4",alpha.f = 0.5),ylab=paste0(j," (%)"),type='n',main="DMSO")
  points(temp$logdays,temp$gate3,bg=temp$col,pch=21,cex=1)
  axis(2,at=log2(c(0.001,0.01,0.10,1.00,10.00,100.00)),labels=c(0.001,0.01,0.1,1,10,100),las=2)
  if(j%in%js[c(10,12)])
  {
    #legend("topright",legend=c(newpvs,"FALSE"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=virnams[k])		
    legend("topright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)		
  } else {
    #legend("bottomright",legend=c(newpvs,"FALSE"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=virnams[k])
    legend("bottomright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)
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
  
  #raw data
  temp$gate=tcell.rawdata[[i]][match(temp$event_id,rownames(tcell.rawdata[[i]])),j]
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
  
  if(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))<0.001)
  {
    pints=sprintf("%.1E",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
  } else {
    pints=round(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))),3)
  }
  #newpvs=paste0(c("Slope",levels(temp$var)[-1]),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))		
  if(pchisq(aa1,aa2, lower.tail = F)<0.001)
  {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))				
  } else {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints))		
  }
  
  plot(temp$logdays,temp$gate3,xlab="Days post-admission",yaxt='n',xaxt='n',pch=21,bg=adjustcolor("azure4",alpha.f = 0.5),ylab=paste0(j," (%)"),type='n',main="Spike-uncorrected")
  points(temp$logdays,temp$gate3,bg=temp$col,pch=21,cex=1)
  axis(2,at=log2(c(0.001,0.01,0.10,1.00,10.00,100.00)),labels=c(0.001,0.01,0.1,1,10,100),las=2)
  if(j%in%js[c(10,12)])
  {
    #legend("topright",legend=c(newpvs,"FALSE"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=virnams[k])		
    legend("topright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)		
  } else {
    #legend("bottomright",legend=c(newpvs,"FALSE"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=virnams[k])
    legend("bottomright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)
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

dev.off()

##Supp Figure 1
svg(paste0("output/new_supp_fig1a.svg"),width=10.4,height=14.9) 
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
legend("topright",legend=c("Low-responder","High-responder"),pch=22,pt.bg=c("orange","magenta"),pt.cex=1.5,bty='n')

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
  
  if(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))<0.001)
  {
    pints=sprintf("%.1E",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
  } else {
    pints=round(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))),3)
  }
  #newpvs=paste0(c("Slope",levels(temp$var)[-1]),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))		
  if(pchisq(aa1,aa2, lower.tail = F)<0.001)
  {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))				
  } else {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints))		
  }
  
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
    legend("bottomright",legend=c(newpvs,"Low-responder","High-responder"),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",c("orange","magenta")),pt.cex=1.5)		
  } else {
    legend("right",legend=c(newpvs,"Low-responder","High-responder"),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",c("orange","magenta")),pt.cex=1.5)		
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
  if(summary(m)$coefficients[2,4]<0.001)
  {
    newpvs=sprintf("%.1e",summary(m)$coefficients[2,4])
  } else {
    newpvs=round(summary(m)$coefficients[2,4],3)
  }
  ymax=100
  tmp=barplot(tmp2,main="",ylab="High-responder (%)",xaxt='n',col=unlist(lapply(vircols[[k]],function(x){rep(x,1)})),las=2,xlab=virnams[k],ylim=c(0,ymax))
  axis(1,tmp,virlvls[[k]],tick=F,las=1)
  legend("topright",legend=paste0("P=",newpvs),bty='n')
  text(tmp,rep(5,length(tmp)),labels=unlist(tttt))
  
  
  m=gam(responder_c2~var+s(logdays,bs="cr",by=var),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$responder_c2),],random =list(id=~1),family="binomial")
  a1=anova(m)
  aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
  aa2=sum(a1$s.table[-c(1), 2])
  
  if(pchisq(aa1,aa2, lower.tail = F)<0.001)
  {
    newpvs=paste0(c("Slope"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F))))				
  } else {
    newpvs=paste0(c("Slope"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3)))		
  }
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
  tmp=barplot(unlist(ttt),main="",ylab="High-responder (%)",xaxt='n',col=unlist(lapply(vircols[[k]],function(x){rep(x,3)})),las=2,xlab="Days post-admission",ylim=c(0,ymax))
  axis(1,tmp,rep(tnams,lgth),tick=F,las=1)
  legend("topright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,rep(22,lgth)),pt.bg=c("white",vircols[[k]]),pt.cex=1.5)
  text(tmp,rep(5,length(tmp)),labels=unlist(tttt))
}

dev.off()

##Suppp Figure 2
svg(paste0("output/new_supp_fig1b.svg"),width=10.4,height=14.9) 
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
legend("topright",legend=c("Low-responder","High-responder"),pch=22,pt.bg=c("orange","magenta"),pt.cex=1.5,bty='n')

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
  
  if(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))<0.001)
  {
    pints=sprintf("%.1E",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
  } else {
    pints=round(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))),3)
  }
  #newpvs=paste0(c("Slope",levels(temp$var)[-1]),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))		
  if(pchisq(aa1,aa2, lower.tail = F)<0.001)
  {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))				
  } else {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints))		
  }
  
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
    legend("bottomright",legend=c(newpvs,"Low-responder","High-responder"),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",c("orange","magenta")),pt.cex=1.5)		
  } else {
    legend("right",legend=c(newpvs,"Low-responder","High-responder"),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",c("orange","magenta")),pt.cex=1.5)		
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
  if(summary(m)$coefficients[2,4]<0.001)
  {
    newpvs=sprintf("%.1e",summary(m)$coefficients[2,4])
  } else {
    newpvs=round(summary(m)$coefficients[2,4],3)
  }
  ymax=100
  tmp=barplot(tmp2,main="",ylab="High-responder (%)",xaxt='n',col=unlist(lapply(vircols[[k]],function(x){rep(x,1)})),las=2,xlab=virnams[k],ylim=c(0,ymax))
  axis(1,tmp,virlvls[[k]],tick=F,las=1)
  legend("topright",legend=paste0("P=",newpvs),bty='n')
  text(tmp,rep(5,length(tmp)),labels=unlist(tttt))
  
  
  m=gam(responder_c2~var+s(logdays,bs="cr",by=var),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$responder_c2),],random =list(id=~1),family="binomial")
  a1=anova(m)
  aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
  aa2=sum(a1$s.table[-c(1), 2])
  
  if(pchisq(aa1,aa2, lower.tail = F)<0.001)
  {
    newpvs=paste0(c("Slope"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F))))				
  } else {
    newpvs=paste0(c("Slope"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3)))		
  }
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
  tmp=barplot(unlist(ttt),main="",ylab="High-responder (%)",xaxt='n',col=unlist(lapply(vircols[[k]],function(x){rep(x,3)})),las=2,xlab="Days post-admission",ylim=c(0,ymax))
  axis(1,tmp,rep(tnams,lgth),tick=F,las=1)
  legend("topright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,rep(22,lgth)),pt.bg=c("white",vircols[[k]]),pt.cex=1.5)
  text(tmp,rep(5,length(tmp)),labels=unlist(tttt))
}

dev.off()

###CD4/CD8
##Supp figure 3
svg(paste0("output/new_supp_fig3.svg"),width=14,height=14.9)  
layout(matrix(c(1,1,1,5,5,5,9,9,9,13,13,13,1,1,1,5,5,5,9,9,9,13,13,13,2,2,2,6,6,6,10,10,10,14,14,14,2,2,2,6,6,6,10,10,10,14,14,14,3,3,3,7,7,7,11,11,11,15,15,15,3,3,3,7,7,7,11,11,11,15,15,15,4,4,4,8,8,8,12,12,12,16,16,16,4,4,4,8,8,8,12,12,12,16,16,16),nrow=8,byrow=T))
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
    
    if(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))<0.001)
    {
      pints=sprintf("%.1E",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
    } else {
      pints=round(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))),3)
    }
    if(pchisq(aa1,aa2, lower.tail = F)<0.001)
    {
      newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))				
    } else {
      newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints))		
    }
    
    plot(temp$logdays,temp$var,xlab="Days post-admission",pch=21,bg=adjustcolor("orange",alpha.f = 0.5),ylab=virnams[k],type='n',xaxt='n',las=2,main=j)
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
      legend(1.9,1,legend=c(newpvs,"Non-responder","Responder"),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)		
    } else {
      #		legend("topright",legend=c(newpvs,"Non-responder"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=j)
      legend(1.8,4.75,legend=c(newpvs,"Non-responder","Responder"),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)
    }
  }
}

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
    if(summary(m)$coefficients[2,4]<0.001)
    {
      oldnewpvs=sprintf("%.1e",summary(m)$coefficients[2,4])
    } else {
      oldnewpvs=round(summary(m)$coefficients[2,4],3)
    }
    ymax=100
    #tmp=barplot(tmp2,main="",ylab="Responder (%)",xaxt='n',col=unlist(lapply(vircols[[k]],function(x){rep(x,1)})),las=2,xlab=virnams[k],ylim=c(0,ymax))
    #axis(1,tmp,c("FALSE","TRUE"),tick=F,las=1)
    #legend("topright",legend=paste0("P=",newpvs),bty='n')
    #text(tmp,rep(5,length(tmp)),labels=unlist(tttt))
    
    m=gam(responder~var+s(logdays,bs="cr",by=var),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$responder),],random =list(id=~1),family="binomial")
    a1=anova(m)
    aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
    aa2=sum(a1$s.table[-c(1), 2])
    
    if(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))<0.001)
    {
      pints=sprintf("%.1E",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
    } else {
      pints=round(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))),3)
    }
    if(pchisq(aa1,aa2, lower.tail = F)<0.001)
    {
      newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),oldnewpvs))				
    } else {
      newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),oldnewpvs))		
    }
    
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
    tmp=barplot(unlist(ttt),main=j,ylab=paste0("Responder (%)"),xaxt='n',col=unlist(lapply(vircols[[k]],function(x){rep(x,3)})),las=2,xlab="Days post-admission",ylim=c(0,ymax))
    axis(1,tmp,rep(tnams,lgth),tick=F,las=1)
    legend("topright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(22,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)
    text(tmp,rep(5,length(tmp)),labels=unlist(tttt))
  }
}

dev.off()

##Supp figure 7
svg(paste0("output/new_supp_fig7.svg"),width=10.4,height=14.9) 
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
  
  if(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))<0.001)
  {
    pints=sprintf("%.1E",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
  } else {
    pints=round(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))),3)
  }
  if(pchisq(aa1,aa2, lower.tail = F)<0.001)
  {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))				
  } else {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints))		
  }
  
  plot(temp$logdays,temp$var,xlab="Days post-admission",pch=21,bg=adjustcolor("orange",alpha.f = 0.5),ylab=virnams[k],main=j,type='n',xaxt='n',las=2)
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
    legend(1.9,1,legend=c(newpvs,"Non-responder","Responder"),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)		
  } else {
    #		legend("topright",legend=c(newpvs,"Non-responder"),bty='n',pch=c(NA,rep(21,lgth)),pt.bg=c("white",rev(vircols[[k]])),pt.cex=1.5,title=j)
    legend(1.9,4.75,legend=c(newpvs,"Non-responder","Responder"),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)
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
  
  if(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/")))<0.001)
  {
    pints=sprintf("%.1E",as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))))
  } else {
    pints=round(as.numeric(unlist(strsplit(paste0(summary(m)$p.table[-1,4],collapse="/"),"/"))),3)
  }
  if(pchisq(aa1,aa2, lower.tail = F)<0.001)
  {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(sprintf("%.1E",pchisq(aa1,aa2, lower.tail = F)),pints))				
  } else {
    newpvs=paste0(c("Slope","Intercept"),"-P=",c(round(pchisq(aa1,aa2, lower.tail = F),3),pints))		
  }
  
  plot(temp$logdays,temp$gate3,xlab="Days post-admission",yaxt='n',xaxt='n',pch=21,bg=adjustcolor("azure4",alpha.f = 0.5),ylab=paste0("%"),type='n',main=j)
  points(temp$logdays,temp$gate3,bg=temp$col,pch=21,cex=1)
  axis(2,at=log2(c(0.001,0.01,0.10,1.00,10.00,100.00)),labels=c(0.001,0.01,0.1,1,10,100),las=2)
  if(j%in%js[c(10,12)])
  {
    legend("topright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)		
  } else {
    legend("bottomright",legend=c(newpvs,virlvls[[k]]),bty='n',pch=c(NA,NA,rep(21,lgth)),pt.bg=c("white","white",vircols[[k]]),pt.cex=1.5)
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


dev.off()
