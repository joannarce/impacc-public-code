library(nnet)
library(foreign)
library(plyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(gamm4)
library(RColorBrewer)
library(nlme)
library(pvca)
library(MASS)
library(flowCore)
library(reshape2)
library(pheatmap)

setwd("~/type2_tcell/")

# Data import

## Clinical metadata


# Patient metadata
clin.individ <- read.csv(
  "/data/clinical/current/2023-01-01-locked/2023-01-01-impacc-clin-individ-locked.csv",
  header = TRUE)

# Sample
clin.sample <- read.csv(
  "/data/clinical/current/2023-01-01-locked/2023-01-01-impacc-clin-sample-locked.csv",
  header = TRUE)

# Event
clin.event <- read.csv(
  "/data/clinical/current/2023-01-01-locked/2023-01-01-impacc-clin-event-locked.csv",
  header = TRUE)

#Convalescent
clin.conv <- read.csv(
  "/scratch/clinical/current/impacc-convalescent.csv",
  header=TRUE)

#Vaccines
clin.vacc <- read.csv(
  "/scratch/data-fullcohort/vaccination/June2023_vaccine_dataset.csv",
  header=TRUE)
clin.vacc=clin.vacc[!is.na(clin.vacc$vaccine),] #remove patients without vaccine info

#supplementary clinical data
clin.supp <- read.csv(
  "/data/clinical/current/2023-01-01-locked/2023-01-01-impacc-clin-individ-supp-locked.csv",
  header=TRUE)



Merge metadata


# Combine sample and event metadata
clin.sample.event <- merge(
  clin.sample, clin.event,
  by = "event_id",
  all.x = TRUE, all.y = FALSE,
  suffixes = c("",".y"),
  sort = FALSE
)

# Add patient metadata
clin.sample.event.individ <- merge(
  clin.sample.event, clin.individ,
  by = "participant_id", all = FALSE, sort = FALSE
)


Manually correct some typos.


# Correct the symptom_date of patient 008-0103: change from -376 to -13 (confirmed with data management team)
clin.sample.event.individ[clin.sample.event.individ$participant_id=="008-0103","symptom_date"] <- -13

# Correct event_date of 004-0012-8: change from 2 to 32 (it's 32 in 2021-11-11-frozen data)
clin.sample.event.individ[clin.sample.event.individ$event_id=="004-0012-8","event_date"] <- 32


## tcells data and metadata 


tcell.data <- readRDS(
  "2024.01.10_data_zerosubtracted.Robj")
tcell.rawdata <- readRDS("2024.01.10_data_raw.obj")

#get matched clinical info
keep_events=rownames(tcell.data[[1]])[rownames(tcell.data[[1]])%in%clin.sample.event.individ$event_id]
tcell.info=clin.sample.event.individ[match(keep_events,clin.sample.event.individ$event_id),]
tcell.final=lapply(tcell.data,function(x){x[keep_events,]})
tcell.raw.final=lapply(tcell.rawdata,function(x){x[keep_events,]})

# live/CD3+ counts for data and the negative control
livecd3=readRDS("2024.01.10_data_livecd3_count.Robj")
names(livecd3)=names(tcell.final)
nc.livecd3=readRDS("2024.01.10_nc_livecd3_count.Robj")

livecd3.final=lapply(livecd3,function(x){x[match(keep_events,rownames(tcell.data[[1]]))]})
for(i in 1:4)
{
  names(livecd3.final[[i]])=tcell.info$event_id
}
nc.livecd3.final=nc.livecd3[match(keep_events,rownames(tcell.data[[1]]))]
names(nc.livecd3.final)=tcell.info$event_id

## Import antibody data

#antibody data
abd <- read.csv("/data/serum-rbd-abtiters/current/serum-rbd-abtiters-Counts.csv",
                header = TRUE, row.names = 1
)
#metadata
abd.metadata <- read.csv("/data/serum-rbd-abtiters/current/serum-rbd-abtiters-Metadata.csv",
                         header = TRUE, row.names = 1
)

# Only keep samples that pass QC
abd=abd[abd.metadata$sample_status=="IMPACC sample assayed and passed QC",]
abd.metadata <- subset(abd.metadata, sample_status=="IMPACC sample assayed and passed QC")

abd.metadata$event_id=clin.sample.event.individ$event_id[match(rownames(abd),clin.sample.event.individ$sample_id)]


## Import viral load data


# mNGS metadata
mngs.metadata <- read.csv(
  # "/data/nasal-metagenomics/current/nasal-metagenomics-Metadata.csv",
  "/data/nasal-metagenomics/legacy/2022-09-03/nasal-metagenomics-Metadata.csv",
  header = TRUE, row.names = NULL,
  check.names = FALSE
)

mngs.metadata <- mngs.metadata %>%
  filter(!is.na(sc2_rpm))

# mNGS QC
mngs.metadata <- mngs.metadata %>%
  filter(sample_status == "IMPACC sample assayed and passed QC")


# Reformat visit_num: decapitalize "VISIT"
mngs.metadata$visit_num <- gsub(
  "^VISIT", "Visit", mngs.metadata$visit_num
)

# Merge with sample metadata
clin.sample.event.individ <- merge(
  clin.sample.event.individ, mngs.metadata[,c("participant_id","visit_num","sc2_rpm")],
  by.x=c("participant_id","event_type"), by.y=c("participant_id","visit_num"),
  all.x=TRUE, all.y=FALSE # keep samples that don't pass SARS-CoV-2 metatranscriptomics QC
  # all=FALSE # discard samples that don't pass SARS-CoV-2 metatranscriptomics QC <<< do NOT do this
)

# log-transform viral load
clin.sample.event.individ$sc2_rpm_log10 <- log10(as.numeric(clin.sample.event.individ$sc2_rpm)+1)


##import omics data

#serum Olink
olink=read.delim("/data/serum-olink/current/Olink-Counts.csv",header=T,sep=",")
olink=cbind(olink,clin.sample[match(olink$X,clin.sample$sample_id),])

#plasma proteomics
prt=read.delim("/data/proteomics/plasma-proteomics-targeted/current/ClassicalPlasmaProteomics-Counts.csv",header=T,sep=",")
prt.metadata=read.delim("/data/proteomics/plasma-proteomics-targeted/current/ClassicalPlasmaProteomics-Metadata.csv",header=T,sep=",")
prt$sample_id=prt.metadata$sample_id[match(prt$core_specific_ID,prt.metadata$core_specific_ID)]
prt$event_id=clin.sample$event_id[match(prt$sample_id,clin.sample$sample_id)]

#blood CyTOF
cytof=read.delim("/data/bld-cytof/current/bld-cytof-Counts.csv",header=T,sep=",")
cytof=cbind(cytof,clin.sample[match(cytof$sample_id,clin.sample$sample_id),])
cytof[,1]=NULL

cytof.metadata=read.delim("/data/bld-cytof/current/bld-cytof-Metadata.csv",header=T,sep=",")
cytof.metadata=cytof.metadata[match(cytof$sample_id,cytof.metadata$sample_id),]

#columns 1:66 are cell populations, optionally convert to %
cytof.pct=cytof
cytof.pct[,1:66]=(cytof[,1:66])/(cytof.metadata$total_cell_counts)*100

# Analysis: longitudinal

##keep samples up to visit-6 (n=356)
tcell.info.acute=tcell.info[tcell.info$event_type%in%paste0("Visit ",1:6),]
tcell.final.acute=lapply(tcell.final,function(x){x[tcell.info$event_type%in%paste0("Visit ",1:6),]})

is=names(tcell.final)[2:4]
js=colnames(tcell.final.acute[[1]])[c(3,8,9,10,15,21,48:53)]

for(i in 1:4)
{
  tmp1=rowSums(tcell.final.acute[[i]][,c(31,32,33)],na.rm=T)
  tmp2=rowSums(tcell.final.acute[[i]][,c(33,34,31)],na.rm=T)
  tcell.final.acute[[i]][,"CD4+IFNg+"]=apply(cbind(tmp1,tmp2),1,max,na.rm=T)
  tmp1=rowSums(tcell.final.acute[[i]][,c(40,41,42)],na.rm=T)
  tmp2=rowSums(tcell.final.acute[[i]][,c(40,42,43)],na.rm=T)
  tcell.final.acute[[i]][,"CD8+IFNg+"]=apply(cbind(tmp1,tmp2),1,max,na.rm=T)
  tcell.final.acute[[i]][,"CD4+IL2+"]=rowSums(tcell.final.acute[[i]][,c(31,35)],na.rm=T)
  tcell.final.acute[[i]][,"CD8+IL2+"]=rowSums(tcell.final.acute[[i]][,c(40,44)],na.rm=T)
  tcell.final.acute[[i]][,"CD4+TNFa+"]=rowSums(tcell.final.acute[[i]][,c(33,37)],na.rm=T)
  tcell.final.acute[[i]][,"CD8+TNFa+"]=rowSums(tcell.final.acute[[i]][,c(42,46)],na.rm=T)
}

#bin dates
tcell.info.acute$event_date_binned=NA
tcell.info.acute$event_date_binned[tcell.info.acute$event_date<4]="1-3"
tcell.info.acute$event_date_binned[tcell.info.acute$event_date>3&tcell.info.acute$event_date<11]="4-10"
tcell.info.acute$event_date_binned[tcell.info.acute$event_date>10]="11+"
tcell.info.acute$event_date_binned=factor(tcell.info.acute$event_date_binned,c("1-3","4-10","11+"))

tcell.info.acute$event_date_symptom_binned=NA
tcell.info.acute$event_date_symptom_binned[tcell.info.acute$event_date_symptom<11]="2-10"
tcell.info.acute$event_date_symptom_binned[tcell.info.acute$event_date_symptom>10&tcell.info.acute$event_date_symptom<23]="11-22"
tcell.info.acute$event_date_symptom_binned[tcell.info.acute$event_date_symptom>22]="23+"
tcell.info.acute$event_date_symptom_binned=factor(tcell.info.acute$event_date_symptom_binned,c("2-10","11-22","23+"))

#setup variables to analyse
#get viral load data
tcell.info.acute$sc2_rpm_log10=clin.sample.event.individ$sc2_rpm_log10[match(tcell.info.acute$event_id,clin.sample.event.individ$event_id)]
tcell.info.acute$rbd_igg=abd$AUC.RBD.IgG[match(tcell.info.acute$event_id,abd.metadata$event_id)]
tcell.info.acute$spike_igg=abd$AUC.Spike.IgG[match(tcell.info.acute$event_id,abd.metadata$event_id)]
#get max per patient
tcell.info.acute$sc2_rpm_log10_max=NA
tcell.info.acute$rbd_igg_max=NA
tcell.info.acute$spike_igg_max=NA
for(i in 1:nrow(tcell.info.acute))
{
  tcell.info.acute$sc2_rpm_log10_max[i]=max(clin.sample.event.individ$sc2_rpm_log10[clin.sample.event.individ$participant_id%in%tcell.info.acute$participant_id[i]],na.rm=T)
  tcell.info.acute$rbd_igg_max[i]=max(abd$AUC.RBD.IgG[abd.metadata$participant_id%in%tcell.info.acute$participant_id[i]],na.rm=T)
  tcell.info.acute$spike_igg_max[i]=max(abd$AUC.Spike.IgG[abd.metadata$participant_id%in%tcell.info.acute$participant_id[i]],na.rm=T)
}
tcell.info.acute$sc2_rpm_log10_max[is.infinite(tcell.info.acute$sc2_rpm_log10_max)]=NA
#first, dichotomise viral load (max and current)
quantile(tcell.info.acute$sc2_rpm_log10,na.rm=T)
quantile(tcell.info.acute$sc2_rpm_log10_max,na.rm=T)
table(tcell.info.acute$sc2_rpm_log10>1)/sum(table(tcell.info.acute$sc2_rpm_log10>1))
table(tcell.info.acute$sc2_rpm_log10_max>1)/sum(table(tcell.info.acute$sc2_rpm_log10_max>1))
#both medians around 1, 53% and 59% >1
tcell.info.acute$vrld=as.factor(tcell.info.acute$sc2_rpm_log10>1)
tcell.info.acute$vrld_max=as.factor(tcell.info.acute$sc2_rpm_log10_max>1)
#for antibodies, >1 as seropositive
tcell.info.acute$rbd_igg_pos=as.factor(tcell.info.acute$rbd_igg>1)
tcell.info.acute$spike_igg_pos=as.factor(tcell.info.acute$spike_igg>1)

for(i in c(7,11,15,19,20,21,22,23,27))
{
  tcell.info.acute[,i]=as.numeric(tcell.info.acute[,i])
}

##setup variables
tcell.info.acute$rbd_igg_log10=log10(tcell.info.acute$rbd_igg)
tcell.info.acute$spike_igg_log10=log10(tcell.info.acute$spike_igg)
tcell.info.acute$tg45=tcell.info.acute$trajectory_group>3
tcell.info.acute$died_ever=tcell.indivs.supp$diedever[match(tcell.info.acute$participant_id,tcell.indivs.supp$participant_id)]
tcell.info.acute$tg123_4=NA
tcell.info.acute$tg123_4[tcell.info.acute$trajectory_group<4]=F
tcell.info.acute$tg123_4[tcell.info.acute$trajectory_group==4]=T
tcell.info.acute$tg123=factor(tcell.info.acute$trajectory_group,levels=c(1,2,3))
tcell.info.acute$tg12_3=tcell.info.acute$tg123==3
tcell.info.acute$admit_icu=tcell.indivs.supp$hosp_admitlevel[match(tcell.info.acute$participant_id,tcell.indivs.supp$participant_id)]<3
tcell.info.acute$respiratory_status_peak=NA
for(i in 1:nrow(tcell.info.acute))
{
  tcell.info.acute$respiratory_status_peak[i]=max(clin.event[clin.event$participant_id%in%tcell.info.acute$participant_id[i],"respiratory_status"],na.rm=T)
}
tcell.info.acute$respiratory_status_peak_567=tcell.info.acute$respiratory_status_peak>4
tcell.info.acute$respiratory_status_567=tcell.info.acute$respiratory_status>4
tcell.info.acute$died_d28=tcell.info.acute$died_ever
tcell.info.acute$died_d28[tcell.info.acute$participant_id%in%tcell.indivs.supp$participant_id[which(tcell.indivs.supp$diedever==1 & tcell.indivs.supp$timetodeath>28)]]=0

##add in ABD to rpm ratio
tcell.info.acute$ratio_rpm=log10(tcell.info.acute$rbd_igg/tcell.info.acute$sc2_rpm_log10)
tcell.info.acute$ratio_rpm[is.infinite(tcell.info.acute$ratio_rpm)]=NA
tcell.info.acute$ratio_rpm_min=log10(tcell.info.acute$rbd_igg/(tcell.info.acute$sc2_rpm_log10+0.008))
tcell.info.acute$ratio2=tcell.info.acute$ratio_rpm_min>2

virvars=c("sc2_rpm_log10","vrld","rbd_igg_log10","rbd_igg_pos","spike_igg_log10","spike_igg_pos","respiratory_status","respiratory_status_567","respiratory_status_peak","respiratory_status_peak_567","trajectory_group","tg45","tg123_4","tg123","died_ever","died_d28","admit_icu","tg12_3","ratio_rpm_min","ratio2","sc2_rpm_log10_max","vrld_max")
virnams=c("Viral rpM","Viral rpM >median","RBD IgG","RBD sero+","Spike IgG","Spike sero+","Respiratory score","Respiratory score >4","Peak respiratory score","Peak respirartory score >4","Trajectory group","Trajectory group >3","TG1-3 vs 4","TG123","Mortality","D28 mortality","ICU upon admission","TG1-2 vs 3","RBD:viral rpM","RBD:viral rpM >2","Peak viral rpM","Peak viral rpM >median")
virmodel=c("lm","glm","lm","glm","lm","glm","polr","glm","polr","glm","polr","glm","glm","polr","glm","glm","glm","glm","lm","glm","lm","glm")
names(virmodel)=virvars
names(virnams)=virvars

#place to store results
acute.pvalues=vector("list",length(is))
names(acute.pvalues)=is
for(i in 1:length(acute.pvalues))
{
  acute.pvalues[[i]]=vector("list",length(comps))
  names(acute.pvalues[[i]])=comps
  for(j in 1:length(acute.pvalues[[i]]))
  {
    acute.pvalues[[i]][[j]]=data.frame(matrix(nrow=length(js),ncol=length(virvars)))
    colnames(acute.pvalues[[i]][[j]])=virvars
    rownames(acute.pvalues[[i]][[j]])=js
  }
}

#also store coefficent (coef oe exp(coef))
acute.coefs=acute.pvalues
acute.pvalues.int=acute.pvalues

#run 
##test for interaction of each clinical variable with time
for(i in is) #only run on the 3 SARS stimulations, not the positive control
{
  
  for(j in js)  #run on a few key gates
  {
    for(k in virvars)
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
      
      temp$logdays=log(temp$days) #needed for lme
      temp$logdays_symptom=log(temp$days_symptom) #needed for lme
      if(virmodel[k]!="lm")
      {
        temp$var=factor(temp$var)
        tryCatch({
          m=gam(gate2~var+s(days,bs="cr",by=var),data=temp[!is.na(temp$days)&!is.na(temp$var)&!is.na(temp$gate2),],random =list(id=~1),family="gaussian")
          a1=anova(m)
          aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
          aa2=sum(a1$s.table[-c(1), 2])
          acute.pvalues[[i]][[1]][j,k]=pchisq(aa1,aa2, lower.tail = F)
          acute.coefs[[i]][[1]][j,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
          acute.pvalues.int[[i]][[1]][j,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
        },error=function(e){})
        
        tryCatch({
          m=gam(gate2~var+s(logdays,bs="cr",by=var),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$gate2),],random =list(id=~1),family="gaussian")
          a1=anova(m)
          aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
          aa2=sum(a1$s.table[-c(1), 2])
          acute.pvalues[[i]][[2]][j,k]=pchisq(aa1,aa2, lower.tail = F)
          acute.coefs[[i]][[2]][j,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
          acute.pvalues.int[[i]][[2]][j,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
        },error=function(e){})
        
        tryCatch({
          m=gam(gate3~var+s(days,bs="cr",by=var),data=temp[!is.na(temp$days)&!is.na(temp$var)&!is.na(temp$gate3),],random =list(id=~1),family="gaussian")
          a1=anova(m)
          aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
          aa2=sum(a1$s.table[-c(1), 2])
          acute.pvalues[[i]][[3]][j,k]=pchisq(aa1,aa2, lower.tail = F)
          acute.coefs[[i]][[3]][j,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
          acute.pvalues.int[[i]][[3]][j,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
        },error=function(e){})
        
        tryCatch({
          m=gam(gate3~var+s(logdays,bs="cr",by=var),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$gate3),],random =list(id=~1),family="gaussian")
          a1=anova(m)
          aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
          aa2=sum(a1$s.table[-c(1), 2])
          acute.pvalues[[i]][[4]][j,k]=pchisq(aa1,aa2, lower.tail = F)
          acute.coefs[[i]][[4]][j,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
          acute.pvalues.int[[i]][[4]][j,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
        },error=function(e){})
        
        tryCatch({
          m=gam(gate2~var+s(days_symptom,bs="cr",by=var),data=temp[!is.na(temp$days_symptom)&!is.na(temp$var)&!is.na(temp$gate2),],random =list(id=~1),family="gaussian")
          a1=anova(m)
          aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
          aa2=sum(a1$s.table[-c(1), 2])
          acute.pvalues[[i]][[5]][j,k]=pchisq(aa1,aa2, lower.tail = F)
          acute.coefs[[i]][[5]][j,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
          acute.pvalues.int[[i]][[5]][j,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
        },error=function(e){})
        
        tryCatch({
          m=gam(gate2~var+s(logdays_symptom,bs="cr",by=var),data=temp[!is.na(temp$logdays_symptom)&!is.na(temp$var)&!is.na(temp$gate2),],random =list(id=~1),family="gaussian")
          a1=anova(m)
          aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
          aa2=sum(a1$s.table[-c(1), 2])
          acute.pvalues[[i]][[6]][j,k]=pchisq(aa1,aa2, lower.tail = F)
          acute.coefs[[i]][[6]][j,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
          acute.pvalues.int[[i]][[6]][j,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
        },error=function(e){})
        
        tryCatch({
          m=gam(gate3~var+s(days_symptom,bs="cr",by=var),data=temp[!is.na(temp$days_symptom)&!is.na(temp$var)&!is.na(temp$gate3),],random =list(id=~1),family="gaussian")
          a1=anova(m)
          aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
          aa2=sum(a1$s.table[-c(1), 2])
          acute.pvalues[[i]][[7]][j,k]=pchisq(aa1,aa2, lower.tail = F)
          acute.coefs[[i]][[7]][j,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
          acute.pvalues.int[[i]][[7]][j,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
        },error=function(e){})
        
        tryCatch({
          m=gam(gate3~var+s(logdays_symptom,bs="cr",by=var),data=temp[!is.na(temp$logdays_symptom)&!is.na(temp$var)&!is.na(temp$gate3),],random =list(id=~1),family="gaussian")
          a1=anova(m)
          aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
          aa2=sum(a1$s.table[-c(1), 2])
          acute.pvalues[[i]][[8]][j,k]=pchisq(aa1,aa2, lower.tail = F)
          acute.coefs[[i]][[8]][j,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
          acute.pvalues.int[[i]][[8]][j,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
        },error=function(e){})
        
        tryCatch({
          m=gam(responder~var+s(days,bs="cr",by=var),data=temp[!is.na(temp$days)&!is.na(temp$var)&!is.na(temp$responder),],random =list(id=~1),family="binomial")
          a1=anova(m)
          aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
          aa2=sum(a1$s.table[-c(1), 2])
          acute.pvalues[[i]][[9]][j,k]=pchisq(aa1,aa2, lower.tail = F)
          acute.coefs[[i]][[9]][j,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
          acute.pvalues.int[[i]][[9]][j,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
        },error=function(e){})
        
        tryCatch({
          m=gam(responder~var+s(logdays,bs="cr",by=var),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$responder),],random =list(id=~1),family="binomial")
          a1=anova(m)
          aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
          aa2=sum(a1$s.table[-c(1), 2])
          acute.pvalues[[i]][[10]][j,k]=pchisq(aa1,aa2, lower.tail = F)
          acute.coefs[[i]][[10]][j,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
          acute.pvalues.int[[i]][[10]][j,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
        },error=function(e){})
        
        tryCatch({
          m=gam(responder~var+s(days_symptom,bs="cr",by=var),data=temp[!is.na(temp$days_symptom)&!is.na(temp$var)&!is.na(temp$responder),],random =list(id=~1),family="binomial")
          a1=anova(m)
          aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
          aa2=sum(a1$s.table[-c(1), 2])
          acute.pvalues[[i]][[11]][j,k]=pchisq(aa1,aa2, lower.tail = F)
          acute.coefs[[i]][[11]][j,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
          acute.pvalues.int[[i]][[11]][j,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
        },error=function(e){})
        
        tryCatch({
          m=gam(responder~var+s(logdays_symptom,bs="cr",by=var),data=temp[!is.na(temp$logdays_symptom)&!is.na(temp$var)&!is.na(temp$responder),],random =list(id=~1),family="binomial")
          a1=anova(m)
          aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
          aa2=sum(a1$s.table[-c(1), 2])
          acute.pvalues[[i]][[12]][j,k]=pchisq(aa1,aa2, lower.tail = F)
          acute.coefs[[i]][[12]][j,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
          acute.pvalues.int[[i]][[12]][j,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
        },error=function(e){})
      } else {
        
        temp$responder=factor(temp$responder)
        
        tryCatch({
          m=lme(gate2~var+days+days*var,data=temp[!is.na(temp$days)&!is.na(temp$var)&!is.na(temp$gate2),],random =~1|id)
          m1=lme(gate2~var+days,data=temp[!is.na(temp$days)&!is.na(temp$var)&!is.na(temp$gate2),],random =~1|id)
          acute.pvalues[[i]][[1]][j,k]=summary(m)$tTable[4,5]
          acute.coefs[[i]][[1]][j,k]=summary(m)$tTable[4,1]
          acute.pvalues.int[[i]][[1]][j,k]=summary(m)$tTable[2,5]
        },error=function(e){})
        
        tryCatch({
          m=lme(gate2~var+logdays+logdays*var,data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$gate2),],random =~1|id)
          m1=lme(gate2~var+logdays,data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$gate2),],random =~1|id)
          acute.pvalues[[i]][[2]][j,k]=summary(m)$tTable[4,5]
          acute.coefs[[i]][[2]][j,k]=summary(m)$tTable[4,1]
          acute.pvalues.int[[i]][[2]][j,k]=summary(m)$tTable[2,5]
        },error=function(e){})
        
        tryCatch({
          m=lme(gate3~var+days+days*var,data=temp[!is.na(temp$days)&!is.na(temp$var)&!is.na(temp$gate3),],random =~1|id)
          m1=lme(gate3~var+days,data=temp[!is.na(temp$days)&!is.na(temp$var)&!is.na(temp$gate3),],random =~1|id)
          acute.pvalues[[i]][[3]][j,k]=summary(m)$tTable[4,5]
          acute.coefs[[i]][[3]][j,k]=summary(m)$tTable[4,1]
          acute.pvalues.int[[i]][[3]][j,k]=summary(m)$tTable[2,5]
        },error=function(e){})
        
        tryCatch({
          m=lme(gate3~var+logdays+logdays*var,data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$gate3),],random =~1|id)
          m1=lme(gate3~var+logdays,data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$gate3),],random =~1|id)
          acute.pvalues[[i]][[4]][j,k]=summary(m)$tTable[4,5]
          acute.coefs[[i]][[4]][j,k]=summary(m)$tTable[4,1]
          acute.pvalues.int[[i]][[4]][j,k]=summary(m)$tTable[2,5]
        },error=function(e){})
        
        tryCatch({
          m=lme(gate2~var+days_symptom+days_symptom*var,data=temp[!is.na(temp$days_symptom)&!is.na(temp$var)&!is.na(temp$gate2),],random =~1|id)
          m1=lme(gate2~var+days_symptom,data=temp[!is.na(temp$days_symptom)&!is.na(temp$var)&!is.na(temp$gate2),],random =~1|id)
          acute.pvalues[[i]][[5]][j,k]=summary(m)$tTable[4,5]
          acute.coefs[[i]][[5]][j,k]=summary(m)$tTable[4,1]
          acute.pvalues.int[[i]][[5]][j,k]=summary(m)$tTable[2,5]
        },error=function(e){})
        
        tryCatch({
          m=lme(gate2~var+logdays_symptom+logdays_symptom*var,data=temp[!is.na(temp$logdays_symptom)&!is.na(temp$var)&!is.na(temp$gate2),],random =~1|id)
          m1=lme(gate2~var+logdays_symptom,data=temp[!is.na(temp$logdays_symptom)&!is.na(temp$var)&!is.na(temp$gate2),],random =~1|id)
          acute.pvalues[[i]][[6]][j,k]=summary(m)$tTable[4,5]
          acute.coefs[[i]][[6]][j,k]=summary(m)$tTable[4,1]
          acute.pvalues.int[[i]][[6]][j,k]=summary(m)$tTable[2,5]
        },error=function(e){})
        
        tryCatch({
          m=lme(gate3~var+days_symptom+days_symptom*var,data=temp[!is.na(temp$days_symptom)&!is.na(temp$var)&!is.na(temp$gate3),],random =~1|id)
          m1=lme(gate3~var+days_symptom,data=temp[!is.na(temp$days_symptom)&!is.na(temp$var)&!is.na(temp$gate3),],random =~1|id)
          acute.pvalues[[i]][[7]][j,k]=summary(m)$tTable[4,5]
          acute.coefs[[i]][[7]][j,k]=summary(m)$tTable[4,1]
          acute.pvalues.int[[i]][[7]][j,k]=summary(m)$tTable[2,5]
        },error=function(e){})
        
        tryCatch({
          m=lme(gate3~var+logdays_symptom+logdays_symptom*var,data=temp[!is.na(temp$logdays_symptom)&!is.na(temp$var)&!is.na(temp$gate3),],random =~1|id)
          m1=lme(gate3~var+logdays_symptom,data=temp[!is.na(temp$logdays_symptom)&!is.na(temp$var)&!is.na(temp$gate3),],random =~1|id)
          acute.pvalues[[i]][[8]][j,k]=summary(m)$tTable[4,5]
          acute.coefs[[i]][[8]][j,k]=summary(m)$tTable[4,1]
          acute.pvalues.int[[i]][[8]][j,k]=summary(m)$tTable[2,5]
        },error=function(e){})
        
        tryCatch({
          m=gam(var~responder+s(days,bs="cr",by=responder),data=temp[!is.na(temp$days)&!is.na(temp$var)&!is.na(temp$responder),],random =list(id=~1),family="gaussian")
          a1=anova(m)
          aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
          aa2=sum(a1$s.table[-c(1), 2])
          acute.pvalues[[i]][[9]][j,k]=pchisq(aa1,aa2, lower.tail = F)
          acute.coefs[[i]][[9]][j,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
          acute.pvalues.int[[i]][[9]][j,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
        },error=function(e){})
        
        tryCatch({
          m=gam(var~responder+s(logdays,bs="cr",by=responder),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$responder),],random =list(id=~1),family="gaussian")
          a1=anova(m)
          aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
          aa2=sum(a1$s.table[-c(1), 2])
          acute.pvalues[[i]][[10]][j,k]=pchisq(aa1,aa2, lower.tail = F)
          acute.coefs[[i]][[10]][j,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
          acute.pvalues.int[[i]][[10]][j,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
        },error=function(e){})
        
        tryCatch({
          m=gam(var~responder+s(days_symptom,bs="cr",by=responder),data=temp[!is.na(temp$days_symptom)&!is.na(temp$var)&!is.na(temp$responder),],random =list(id=~1),family="gaussian")
          a1=anova(m)
          aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
          aa2=sum(a1$s.table[-c(1), 2])
          acute.pvalues[[i]][[11]][j,k]=pchisq(aa1,aa2, lower.tail = F)
          acute.coefs[[i]][[11]][j,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
          acute.pvalues.int[[i]][[11]][j,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
        },error=function(e){})
        
        tryCatch({
          m=gam(var~responder+s(logdays_symptom,bs="cr",by=responder),data=temp[!is.na(temp$logdays_symptom)&!is.na(temp$var)&!is.na(temp$responder),],random =list(id=~1),family="gaussian")
          a1=anova(m)
          aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
          aa2=sum(a1$s.table[-c(1), 2])
          acute.pvalues[[i]][[12]][j,k]=pchisq(aa1,aa2, lower.tail = F)
          acute.coefs[[i]][[12]][j,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
          acute.pvalues.int[[i]][[12]][j,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
        },error=function(e){})
        
      }
    }
    print(j)
  }
  print(i)
}	


##merged responder status
i=is[1]
#remove few with NAs
tmp=tcell.final.acute[[i]][apply(tcell.final.acute[[i]],1,function(x){sum(is.na(x))/length(x)*100})==0,]
#PCA to combine
tmpp=prcomp((tmp[,js]>0)+0,scale=F)
#get summaries per sample
tmpc=rowSums(tmp[,js[7:12]]>0)
tmpa=rowSums(tmp[,js[1:6]]>0)

#add 2- and 4-cluster assignment to sample info
tcell.info.acute$responder_c2=NA
tcell.info.acute$responder_c2[tcell.info.acute$event_id%in%rownames(tmp)[tmpk2$cluster==1]]="high"
tcell.info.acute$responder_c2[tcell.info.acute$event_id%in%rownames(tmp)[tmpk2$cluster==2]]="low"
tcell.info.acute$responder_c2=as.factor(tcell.info.acute$responder_c2)
tcell.info.acute$responder_c4=NA
tcell.info.acute$responder_c4[tcell.info.acute$event_id%in%rownames(tmp)[tmpk$cluster==1]]="AIM"
tcell.info.acute$responder_c4[tcell.info.acute$event_id%in%rownames(tmp)[tmpk$cluster==2]]="ICS"
tcell.info.acute$responder_c4[tcell.info.acute$event_id%in%rownames(tmp)[tmpk$cluster==3]]="low"
tcell.info.acute$responder_c4[tcell.info.acute$event_id%in%rownames(tmp)[tmpk$cluster==4]]="high"
tcell.info.acute$responder_c4=factor(tcell.info.acute$responder_c4,c("low","AIM","ICS","high"))

#responder vs outcome regressions
resp.acute.pvalues=vector("list",length(is))
names(resp.acute.pvalues)=is
for(i in 1:length(resp.acute.pvalues))
{
  resp.acute.pvalues[[i]]=vector("list",length(comps[9:12]))
  names(resp.acute.pvalues[[i]])=comps[9:12]
  for(j in 1:length(resp.acute.pvalues[[i]]))
  {
    resp.acute.pvalues[[i]][[j]]=data.frame(matrix(nrow=1,ncol=length(virvars)))
    colnames(resp.acute.pvalues[[i]][[j]])=virvars
  }
}

#also store coefficent (coef oe exp(coef))
resp.acute.coefs=resp.acute.pvalues
resp.acute.pvalues.int=resp.acute.pvalues

i=is[1]

for(k in virvars)
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
  
  temp$logdays=log(temp$days) #needed for lme
  temp$logdays_symptom=log(temp$days_symptom) #needed for lme
  temp$responder_c2=factor(tcell.info.acute$responder_c2[match(temp$event_id,tcell.info.acute$event_id)],c("low","high"))
  
  if(virmodel[k]!="lm")
  {
    temp$var=factor(temp$var)
    
    tryCatch({
      m=gam(responder_c2~var+s(days,bs="cr",by=var),data=temp[!is.na(temp$days)&!is.na(temp$var)&!is.na(temp$responder_c2),],random =list(id=~1),family="binomial")
      a1=anova(m)
      aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
      aa2=sum(a1$s.table[-c(1), 2])
      resp.acute.pvalues[[i]][[1]][1,k]=pchisq(aa1,aa2, lower.tail = F)
      resp.acute.coefs[[i]][[1]][1,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
      resp.acute.pvalues.int[[i]][[1]][1,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
    },error=function(e){})
    
    tryCatch({
      m=gam(responder_c2~var+s(logdays,bs="cr",by=var),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$responder_c2),],random =list(id=~1),family="binomial")
      a1=anova(m)
      aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
      aa2=sum(a1$s.table[-c(1), 2])
      resp.acute.pvalues[[i]][[2]][1,k]=pchisq(aa1,aa2, lower.tail = F)
      resp.acute.coefs[[i]][[2]][1,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
      resp.acute.pvalues.int[[i]][[2]][1,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
    },error=function(e){})
    
    tryCatch({
      m=gam(responder_c2~var+s(days_symptom,bs="cr",by=var),data=temp[!is.na(temp$days_symptom)&!is.na(temp$var)&!is.na(temp$responder_c2),],random =list(id=~1),family="binomial")
      a1=anova(m)
      aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
      aa2=sum(a1$s.table[-c(1), 2])
      resp.acute.pvalues[[i]][[3]][1,k]=pchisq(aa1,aa2, lower.tail = F)
      resp.acute.coefs[[i]][[3]][1,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
      resp.acute.pvalues.int[[i]][[3]][1,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
    },error=function(e){})
    
    tryCatch({
      m=gam(responder_c2~var+s(logdays_symptom,bs="cr",by=var),data=temp[!is.na(temp$logdays_symptom)&!is.na(temp$var)&!is.na(temp$responder_c2),],random =list(id=~1),family="binomial")
      a1=anova(m)
      aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
      aa2=sum(a1$s.table[-c(1), 2])
      resp.acute.pvalues[[i]][[4]][1,k]=pchisq(aa1,aa2, lower.tail = F)
      resp.acute.coefs[[i]][[4]][1,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
      resp.acute.pvalues.int[[i]][[4]][1,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
    },error=function(e){})
  } else {
    
    tryCatch({
      m=gam(var~responder_c2+s(days,bs="cr",by=responder_c2),data=temp[!is.na(temp$days)&!is.na(temp$var)&!is.na(temp$responder_c2),],random =list(id=~1),family="gaussian")
      a1=anova(m)
      aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
      aa2=sum(a1$s.table[-c(1), 2])
      resp.acute.pvalues[[i]][[2]][1,k]=pchisq(aa1,aa2, lower.tail = F)
      resp.acute.coefs[[i]][[2]][1,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
      resp.acute.pvalues.int[[i]][[2]][1,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
    },error=function(e){})
    
    tryCatch({
      m=gam(var~responder_c2+s(logdays,bs="cr",by=responder_c2),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$responder_c2),],random =list(id=~1),family="gaussian")
      a1=anova(m)
      aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
      aa2=sum(a1$s.table[-c(1), 2])
      resp.acute.pvalues[[i]][[2]][1,k]=pchisq(aa1,aa2, lower.tail = F)
      resp.acute.coefs[[i]][[2]][1,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
      resp.acute.pvalues.int[[i]][[2]][1,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
    },error=function(e){})
    
    tryCatch({
      m=gam(var~responder_c2+s(days_symptom,bs="cr",by=responder_c2),data=temp[!is.na(temp$days_symptom)&!is.na(temp$var)&!is.na(temp$responder_c2),],random =list(id=~1),family="gaussian")
      a1=anova(m)
      aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
      aa2=sum(a1$s.table[-c(1), 2])
      resp.acute.pvalues[[i]][[3]][1,k]=pchisq(aa1,aa2, lower.tail = F)
      resp.acute.coefs[[i]][[3]][1,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
      resp.acute.pvalues.int[[i]][[3]][1,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
    },error=function(e){})
    
    tryCatch({
      m=gam(var~responder_c2+s(logdays_symptom,bs="cr",by=responder_c2),data=temp[!is.na(temp$logdays_symptom)&!is.na(temp$var)&!is.na(temp$responder_c2),],random =list(id=~1),family="gaussian")
      a1=anova(m)
      aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
      aa2=sum(a1$s.table[-c(1), 2])
      resp.acute.pvalues[[i]][[4]][1,k]=pchisq(aa1,aa2, lower.tail = F)
      resp.acute.coefs[[i]][[4]][1,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
      resp.acute.pvalues.int[[i]][[4]][1,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
    },error=function(e){})
    
  }
  print(k)
}

#four group	
for(k in virvars)
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
  
  temp$logdays=log(temp$days) #needed for lme
  temp$logdays_symptom=log(temp$days_symptom) #needed for lme
  temp$responder_c2=factor(tcell.info.acute$responder_c2[match(temp$event_id,tcell.info.acute$event_id)],c("low","high"))
  temp$responder_c4=factor(tcell.info.acute$responder_c4[match(temp$event_id,tcell.info.acute$event_id)],c("low","AIM","ICS","high"))
  
  if(virmodel[k]!="lm")
  {
    temp$var=factor(temp$var)
    
    tryCatch({
      m=gam(responder_c4~var+s(days,bs="cr",by=var),data=temp[!is.na(temp$days)&!is.na(temp$var)&!is.na(temp$responder_c4),],random =list(id=~1),family="binomial")
      a1=anova(m)
      aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
      aa2=sum(a1$s.table[-c(1), 2])
      resp.acute.pvalues[[i]][[1]][2,k]=pchisq(aa1,aa2, lower.tail = F)
      resp.acute.coefs[[i]][[1]][2,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
      resp.acute.pvalues.int[[i]][[1]][2,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
    },error=function(e){})
    
    tryCatch({
      m=gam(responder_c4~var+s(logdays,bs="cr",by=var),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$responder_c4),],random =list(id=~1),family="binomial")
      a1=anova(m)
      aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
      aa2=sum(a1$s.table[-c(1), 2])
      resp.acute.pvalues[[i]][[2]][2,k]=pchisq(aa1,aa2, lower.tail = F)
      resp.acute.coefs[[i]][[2]][2,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
      resp.acute.pvalues.int[[i]][[2]][2,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
    },error=function(e){})
    
    tryCatch({
      m=gam(responder_c4~var+s(days_symptom,bs="cr",by=var),data=temp[!is.na(temp$days_symptom)&!is.na(temp$var)&!is.na(temp$responder_c4),],random =list(id=~1),family="binomial")
      a1=anova(m)
      aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
      aa2=sum(a1$s.table[-c(1), 2])
      resp.acute.pvalues[[i]][[3]][2,k]=pchisq(aa1,aa2, lower.tail = F)
      resp.acute.coefs[[i]][[3]][2,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
      resp.acute.pvalues.int[[i]][[3]][2,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
    },error=function(e){})
    
    tryCatch({
      m=gam(responder_c4~var+s(logdays_symptom,bs="cr",by=var),data=temp[!is.na(temp$logdays_symptom)&!is.na(temp$var)&!is.na(temp$responder_c4),],random =list(id=~1),family="binomial")
      a1=anova(m)
      aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
      aa2=sum(a1$s.table[-c(1), 2])
      resp.acute.pvalues[[i]][[4]][2,k]=pchisq(aa1,aa2, lower.tail = F)
      resp.acute.coefs[[i]][[4]][2,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
      resp.acute.pvalues.int[[i]][[4]][2,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
    },error=function(e){})
  } else {
    
    tryCatch({
      m=gam(var~responder_c4+s(days,bs="cr",by=responder_c4),data=temp[!is.na(temp$days)&!is.na(temp$var)&!is.na(temp$responder_c4),],random =list(id=~1),family="gaussian")
      a1=anova(m)
      aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
      aa2=sum(a1$s.table[-c(1), 2])
      resp.acute.pvalues[[i]][[2]][2,k]=pchisq(aa1,aa2, lower.tail = F)
      resp.acute.coefs[[i]][[2]][2,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
      resp.acute.pvalues.int[[i]][[2]][2,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
    },error=function(e){})
    
    tryCatch({
      m=gam(var~responder_c4+s(logdays,bs="cr",by=responder_c4),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$responder_c4),],random =list(id=~1),family="gaussian")
      a1=anova(m)
      aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
      aa2=sum(a1$s.table[-c(1), 2])
      resp.acute.pvalues[[i]][[2]][2,k]=pchisq(aa1,aa2, lower.tail = F)
      resp.acute.coefs[[i]][[2]][2,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
      resp.acute.pvalues.int[[i]][[2]][2,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
    },error=function(e){})
    
    tryCatch({
      m=gam(var~responder_c4+s(days_symptom,bs="cr",by=responder_c4),data=temp[!is.na(temp$days_symptom)&!is.na(temp$var)&!is.na(temp$responder_c4),],random =list(id=~1),family="gaussian")
      a1=anova(m)
      aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
      aa2=sum(a1$s.table[-c(1), 2])
      resp.acute.pvalues[[i]][[3]][2,k]=pchisq(aa1,aa2, lower.tail = F)
      resp.acute.coefs[[i]][[3]][2,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
      resp.acute.pvalues.int[[i]][[3]][2,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
    },error=function(e){})
    
    tryCatch({
      m=gam(var~responder_c4+s(logdays_symptom,bs="cr",by=responder_c4),data=temp[!is.na(temp$logdays_symptom)&!is.na(temp$var)&!is.na(temp$responder_c4),],random =list(id=~1),family="gaussian")
      a1=anova(m)
      aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
      aa2=sum(a1$s.table[-c(1), 2])
      resp.acute.pvalues[[i]][[4]][2,k]=pchisq(aa1,aa2, lower.tail = F)
      resp.acute.coefs[[i]][[4]][2,k]=paste0(summary(m)$p.table[-1,1],collapse="/")
      resp.acute.pvalues.int[[i]][[4]][2,k]=paste0(summary(m)$p.table[-1,4],collapse="/")
    },error=function(e){})
    
  }
  print(k)
}

tcell.info.acute$PC1=unlist(tmpp$x[match(tcell.info.acute$event_id,rownames(tmpp$x)),1])
tcell.info.acute$PC2=unlist(tmpp$x[match(tcell.info.acute$event_id,rownames(tmpp$x)),2])
tcell.info.acute$event_date_log=log(tcell.info.acute$event_date)

tcell.info.acute$age60=factor(as.numeric(tcell.info.acute$admit_age)>=60)
tcell.info.acute$age60=factor(tcell.info.acute$age60)
tcell.info.acute$sex=factor(tcell.info.acute$sex)
tcell.info.acute$bmi30=factor(as.numeric(tcell.info.acute$bmi)>=30)

# Analysis: Olink and CyTOF

#select variables
tmp=colnames(cytof)[unique(c(grep("B.Cell",colnames(cytof)),grep("T.Cell",colnames(cytof)),grep("^Mono",colnames(cytof)),grep("^Hemato",colnames(cytof)),grep("^Plasmacy",colnames(cytof))))]
x=read.xlsx("writeup/IMPACC_Tcells_AssaySurvey.xlsx",sheet=2)
tmp2=x[which(x[,1]>1),2]

omics.vars=c(tmp,tmp2)
omics.type=c(rep("cytof",length(tmp)),rep("olink",length(tmp2)))
names(omics.type)=omics.vars

#store correlation with gate and odds ratio for responder status
omics.res=vector("list",length(is))
names(omics.res)=is

omics.type.df=data.frame("Dataset"=omics.type)
rownames(omics.type.df)=omics.nams

#setup tables for results
for(i in 1:length(is))
{
  omics.res[[i]]=vector("list",length(js))
  names(omics.res[[i]])=js
  for(j in 1:length(js))
  {
    omics.res[[i]][[j]]=data.frame(matrix(ncol=10,nrow=length(omics.vars)))
    rownames(omics.res[[i]][[j]])=omics.vars
    colnames(omics.res[[i]][[j]])=c("rho","rho_pv","rho_nzero","rho_nzero_pv","lme","lme_pv","lme_nzero","lme_nzero_pv","OR","OR_pv")
  }
}

omics.nams=omics.vars
omics.nams[omics.type=="cytof"]=paste0(omics.vars[omics.type=="cytof"]," (%)")
names(omics.nams)=omics.vars
omics.nams[28]="Monocytes..CD14posCD16pos (%)"
omics.nams[29]="Monocytes..CD14posCD16neg (%)"
omics.nams[30]="Monocytes..CD14negCD16pos (%)"

#run it
for(i in is)
{
  for(j in js)
  {
    for(k in omics.vars)
    {
      temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tcell.final.acute[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin)) #get data
      if(omics.type[k]=="cytof")
      {
        temp$var=cytof.pct[match(temp$event_id,cytof$event_id),k]
      } else {
        temp$var=olink[match(temp$event_id,olink$event_id),k]
      }
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
      
      omics.res[[i]][[j]][k,1]=cor.test(temp$var,temp$gate2,method="spearman")$estimate
      omics.res[[i]][[j]][k,2]=cor.test(temp$var,temp$gate2,method="spearman")$p.value
      omics.res[[i]][[j]][k,3]=cor.test(temp$var,temp$gate3,method="spearman")$estimate
      omics.res[[i]][[j]][k,4]=cor.test(temp$var,temp$gate3,method="spearman")$p.value
      m=lme(gate2~var,temp[!is.na(temp$var)&!is.na(temp$gate2),],random =~1|id)
      omics.res[[i]][[j]][k,5]=summary(m)$tTable[2,1]
      omics.res[[i]][[j]][k,6]=summary(m)$tTable[2,5]
      m=lme(gate3~var,temp[!is.na(temp$var)&!is.na(temp$gate3),],random =~1|id)
      omics.res[[i]][[j]][k,7]=summary(m)$tTable[2,1]
      omics.res[[i]][[j]][k,8]=summary(m)$tTable[2,5]
      m=glmer(responder~var+(1|id),temp[!is.na(temp$var)&!is.na(temp$gate2),],family="binomial")
      omics.res[[i]][[j]][k,9]=exp(summary(m)$coefficients[2,1])
      omics.res[[i]][[j]][k,10]=summary(m)$coefficients[2,4]
    }
  }
  print(i)
}

#do log-transformed version for cytof
omics.res.log=omics.res

for(i in is)
{
  for(j in js)
  {
    for(k in omics.vars)
    {
      temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tcell.final.acute[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin)) #get data
      if(omics.type[k]=="cytof")
      {
        temp$var=cytof.pct[match(temp$event_id,cytof$event_id),k]
        
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
        
        temp$var2=log10(temp$var+min(temp$var[which(temp$var>0)]))
        
        omics.res.log[[i]][[j]][k,1]=cor.test(temp$var2,temp$gate2,method="spearman")$estimate
        omics.res.log[[i]][[j]][k,2]=cor.test(temp$var2,temp$gate2,method="spearman")$p.value
        omics.res.log[[i]][[j]][k,3]=cor.test(temp$var2,temp$gate3,method="spearman")$estimate
        omics.res.log[[i]][[j]][k,4]=cor.test(temp$var2,temp$gate3,method="spearman")$p.value
        m=lme(gate2~var2,temp[!is.na(temp$var2)&!is.na(temp$gate2),],random =~1|id)
        omics.res.log[[i]][[j]][k,5]=summary(m)$tTable[2,1]
        omics.res.log[[i]][[j]][k,6]=summary(m)$tTable[2,5]
        m=lme(gate3~var2,temp[!is.na(temp$var2)&!is.na(temp$gate3),],random =~1|id)
        omics.res.log[[i]][[j]][k,7]=summary(m)$tTable[2,1]
        omics.res.log[[i]][[j]][k,8]=summary(m)$tTable[2,5]
        m=glmer(responder~var2+(1|id),temp[!is.na(temp$var2)&!is.na(temp$gate2),],family="binomial")
        omics.res.log[[i]][[j]][k,9]=exp(summary(m)$coefficients[2,1])
        omics.res.log[[i]][[j]][k,10]=summary(m)$coefficients[2,4]
      }
    }
  }
  print(i)
}


#repeat but use baseline omics data per patient
omics.res.baseline=vector("list",length(is))
names(omics.res.baseline)=is
#setup tables for results
for(i in 1:length(is))
{
  omics.res.baseline[[i]]=vector("list",length(js))
  names(omics.res.baseline[[i]])=js
  for(j in 1:length(js))
  {
    omics.res.baseline[[i]][[j]]=data.frame(matrix(ncol=10,nrow=length(omics.vars)))
    rownames(omics.res.baseline[[i]][[j]])=omics.vars
    colnames(omics.res.baseline[[i]][[j]])=c("rho","rho_pv","rho_nzero","rho_nzero_pv","lme","lme_pv","lme_nzero","lme_nzero_pv","OR","OR_pv")
  }
}

#run it
for(i in is)
{
  for(j in js)
  {
    for(k in omics.vars)
    {
      temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tcell.final.acute[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin)) #get data
      if(omics.type[k]=="cytof")
      {
        tmp=cytof[gsub(".*[-]","",cytof$event_id)=="1",]
        tmpp=cytof.pct[gsub(".*[-]","",cytof$event_id)=="1",]
        temp$var=tmpp[match(temp$id,tmp$participant_id),k]
      } else {
        tmp=olink[gsub(".*[-]","",olink$event_id)=="1",]
        temp$var=olink[match(temp$id,olink$participant_id),k]
      }
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
      
      omics.res.baseline[[i]][[j]][k,1]=cor.test(temp$var,temp$gate2,method="spearman")$estimate
      omics.res.baseline[[i]][[j]][k,2]=cor.test(temp$var,temp$gate2,method="spearman")$p.value
      omics.res.baseline[[i]][[j]][k,3]=cor.test(temp$var,temp$gate3,method="spearman")$estimate
      omics.res.baseline[[i]][[j]][k,4]=cor.test(temp$var,temp$gate3,method="spearman")$p.value
      m=lme(gate2~var,temp[!is.na(temp$var)&!is.na(temp$gate2),],random =~1|id)
      omics.res.baseline[[i]][[j]][k,5]=summary(m)$tTable[2,1]
      omics.res.baseline[[i]][[j]][k,6]=summary(m)$tTable[2,5]
      m=lme(gate3~var,temp[!is.na(temp$var)&!is.na(temp$gate3),],random =~1|id)
      omics.res.baseline[[i]][[j]][k,7]=summary(m)$tTable[2,1]
      omics.res.baseline[[i]][[j]][k,8]=summary(m)$tTable[2,5]
      m=glmer(responder~var+(1|id),temp[!is.na(temp$var)&!is.na(temp$gate2),],family="binomial")
      omics.res.baseline[[i]][[j]][k,9]=exp(summary(m)$coefficients[2,1])
      omics.res.baseline[[i]][[j]][k,10]=summary(m)$coefficients[2,4]
    }
  }
  print(i)
}

#do log-transformed version for cytof
omics.res.baseline.log=omics.res.baseline

for(i in is)
{
  for(j in js)
  {
    for(k in omics.vars)
    {
      temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tcell.final.acute[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin)) #get data
      if(omics.type[k]=="cytof")
      {
        tmp=cytof[gsub(".*[-]","",cytof$event_id)=="1",]
        tmpp=cytof.pct[gsub(".*[-]","",cytof$event_id)=="1",]
        temp$var=tmpp[match(temp$id,tmp$participant_id),k]
        
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
        
        temp$var2=log10(temp$var+min(temp$var[which(temp$var>0)]))
        
        omics.res.baseline.log[[i]][[j]][k,1]=cor.test(temp$var2,temp$gate2,method="spearman")$estimate
        omics.res.baseline.log[[i]][[j]][k,2]=cor.test(temp$var2,temp$gate2,method="spearman")$p.value
        omics.res.baseline.log[[i]][[j]][k,3]=cor.test(temp$var2,temp$gate3,method="spearman")$estimate
        omics.res.baseline.log[[i]][[j]][k,4]=cor.test(temp$var2,temp$gate3,method="spearman")$p.value
        m=lme(gate2~var2,temp[!is.na(temp$var2)&!is.na(temp$gate2),],random =~1|id)
        omics.res.baseline.log[[i]][[j]][k,5]=summary(m)$tTable[2,1]
        omics.res.baseline.log[[i]][[j]][k,6]=summary(m)$tTable[2,5]
        m=lme(gate3~var2,temp[!is.na(temp$var2)&!is.na(temp$gate3),],random =~1|id)
        omics.res.baseline.log[[i]][[j]][k,7]=summary(m)$tTable[2,1]
        omics.res.baseline.log[[i]][[j]][k,8]=summary(m)$tTable[2,5]
        m=glmer(responder~var2+(1|id),temp[!is.na(temp$var2)&!is.na(temp$gate2),],family="binomial")
        omics.res.baseline.log[[i]][[j]][k,9]=exp(summary(m)$coefficients[2,1])
        omics.res.baseline.log[[i]][[j]][k,10]=summary(m)$coefficients[2,4]
      }
    }
  }
  print(i)
}

#repeat for total responder, 2 clusters
omics.res.responder=data.frame(matrix(ncol=6,nrow=length(omics.vars)))
rownames(omics.res.responder)=omics.vars
omics.res.responder.log=data.frame(matrix(ncol=6,nrow=length(omics.vars)))
rownames(omics.res.responder.log)=omics.vars

for(k in omics.vars)
{
  temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tcell.final.acute[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin,"responder"=tcell.info.acute$responder_c2,"PC1"=tcell.info.acute$PC1)) #get data
  if(omics.type[k]=="cytof")
  {
    temp$var=cytof.pct[match(temp$event_id,cytof$event_id),k]
  } else {
    temp$var=olink[match(temp$event_id,olink$event_id),k]
  }
  temp$PC1=as.numeric(temp$PC1)
  temp$var2=log10(temp$var+min(temp$var[which(temp$var>0)]))
  temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
  temp$livecd3=livecd3.final[[1]][match(temp$event_id,names(livecd3.final[[1]]))] #get stim live count
  temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
  temp$responder=factor(temp$responder=="1")
  temp$logdays=log(as.numeric(temp$days))
  
  tryCatch({	
    m=glmer(responder~var+(1|id),temp[!is.na(temp$var)&!is.na(temp$responder),],family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    omics.res.responder[k,1]=summary(m)$coefficients[2,4]
    omics.res.responder[k,2]=exp(fixef(m))[2]
    omics.res.responder[k,3:4]=unlist(exp(confint(m,devtol=Inf))[3,])
  },error=function(e){})	
  tryCatch({	
    m=glmer(responder~var2+(1|id),temp[!is.na(temp$var2)&!is.na(temp$responder),],family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    omics.res.responder.log[k,1]=summary(m)$coefficients[2,4]
    omics.res.responder.log[k,2]=exp(fixef(m))[2]
    omics.res.responder.log[k,3:4]=unlist(exp(confint(m,devtol=Inf))[3,])	
  },error=function(e){})
  
  omics.res.responder[k,5]=cor.test(temp$var,temp$PC1,method="spearman")$estimate
  omics.res.responder[k,6]=cor.test(temp$var,temp$PC1,method="spearman")$p.value
  omics.res.responder.log[k,5]=cor.test(temp$var2,temp$PC1,method="spearman")$estimate
  omics.res.responder.log[k,6]=cor.test(temp$var2,temp$PC1,method="spearman")$p.value
  
  tryCatch({	
    m=gam(var~responder+s(logdays,bs="cr",by=responder),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$responder),],random =list(id=~1),family="gaussian")
    a1=anova(m)
    aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
    aa2=sum(a1$s.table[-c(1), 2])
    omics.res.responder[k,7]=pchisq(aa1,aa2, lower.tail = F)
    omics.res.responder[k,8]=summary(m)$p.table[2,4]	
  },error=function(e){})	
  tryCatch({	
    m=gam(var2~responder+s(logdays,bs="cr",by=responder),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$responder),],random =list(id=~1),family="gaussian")
    a1=anova(m)
    aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
    aa2=sum(a1$s.table[-c(1), 2])
    omics.res.responder.log[k,7]=pchisq(aa1,aa2, lower.tail = F)
    omics.res.responder.log[k,8]=summary(m)$p.table[2,4]	
  },error=function(e){})
  
  print(k)
}

#use baseline omics measures
omics.baseline.res.responder=data.frame(matrix(ncol=6,nrow=length(omics.vars)))
rownames(omics.baseline.res.responder)=omics.vars
omics.baseline.res.responder.log=data.frame(matrix(ncol=6,nrow=length(omics.vars)))
rownames(omics.baseline.res.responder.log)=omics.vars

for(k in omics.vars)
{
  temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tcell.final.acute[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin,"responder"=tcell.info.acute$responder_c2,"PC1"=tcell.info.acute$PC1)) #get data
  if(omics.type[k]=="cytof")
  {
    tmp=cytof[gsub(".*[-]","",cytof$event_id)=="1",]
    tmpp=cytof.pct[gsub(".*[-]","",cytof$event_id)=="1",]
    temp$var=tmpp[match(temp$id,tmp$participant_id),k]
  } else {
    tmp=olink[gsub(".*[-]","",olink$event_id)=="1",]
    temp$var=olink[match(temp$id,olink$participant_id),k]
  }
  temp$PC1=as.numeric(temp$PC1)
  temp$var2=log10(temp$var+min(temp$var[which(temp$var>0)]))
  temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
  temp$livecd3=livecd3.final[[1]][match(temp$event_id,names(livecd3.final[[1]]))] #get stim live count
  temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
  temp$responder=factor(temp$responder=="1")
  temp$logdays=log(as.numeric(temp$days))
  
  tryCatch({	
    m=glmer(responder~var+(1|id),temp[!is.na(temp$var)&!is.na(temp$responder),],family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    omics.baseline.res.responder[k,1]=summary(m)$coefficients[2,4]
    omics.baseline.res.responder[k,2]=exp(fixef(m))[2]
    omics.baseline.res.responder[k,3:4]=unlist(exp(confint(m,devtol=Inf))[3,])
  },error=function(e){})	
  tryCatch({	
    m=glmer(responder~var2+(1|id),temp[!is.na(temp$var2)&!is.na(temp$responder),],family="binomial",control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=100000)))
    omics.baseline.res.responder.log[k,1]=summary(m)$coefficients[2,4]
    omics.baseline.res.responder.log[k,2]=exp(fixef(m))[2]
    omics.baseline.res.responder.log[k,3:4]=unlist(exp(confint(m,devtol=Inf))[3,])	
  },error=function(e){})
  
  omics.baseline.res.responder[k,5]=cor.test(temp$var,temp$PC1,method="spearman")$estimate
  omics.baseline.res.responder[k,6]=cor.test(temp$var,temp$PC1,method="spearman")$p.value
  omics.baseline.res.responder.log[k,5]=cor.test(temp$var2,temp$PC1,method="spearman")$estimate
  omics.baseline.res.responder.log[k,6]=cor.test(temp$var2,temp$PC1,method="spearman")$p.value
  
  tryCatch({	
    m=gam(var~responder+s(logdays,bs="cr",by=responder),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$responder),],random =list(id=~1),family="gaussian")
    a1=anova(m)
    aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
    aa2=sum(a1$s.table[-c(1), 2])
    omics.baseline.res.responder[k,7]=pchisq(aa1,aa2, lower.tail = F)
    omics.baseline.res.responder[k,8]=summary(m)$p.table[2,4]	
  },error=function(e){})	
  tryCatch({	
    m=gam(var2~responder+s(logdays,bs="cr",by=responder),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$responder),],random =list(id=~1),family="gaussian")
    a1=anova(m)
    aa1=sum(a1$s.table[-c(1), 3]*a1$s.table[-c(1), 2])
    aa2=sum(a1$s.table[-c(1), 2])
    omics.baseline.res.responder.log[k,7]=pchisq(aa1,aa2, lower.tail = F)
    omics.baseline.res.responder.log[k,8]=summary(m)$p.table[2,4]	
  },error=function(e){})
  
  print(k)
}	

#get medians per responder status to make plot
omics.responder.aves=data.frame(matrix(ncol=length(omics.vars),nrow=4))
colnames(omics.responder.aves)=omics.vars
rownames(omics.responder.aves)=c("low_baseline","high_baseline","low","high")	

for(k in omics.vars)
{
  temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tcell.final.acute[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin,"responder"=tcell.info.acute$responder_c2,"PC1"=tcell.info.acute$PC1)) #get data
  if(omics.type[k]=="cytof")
  {
    tmp=cytof[gsub(".*[-]","",cytof$event_id)=="1",]
    tmpp=cytof.pct[gsub(".*[-]","",cytof$event_id)=="1",]
    temp$var=tmpp[match(temp$id,tmp$participant_id),k]
  } else {
    tmp=olink[gsub(".*[-]","",olink$event_id)=="1",]
    temp$var=olink[match(temp$id,olink$participant_id),k]
  }
  temp$PC1=as.numeric(temp$PC1)
  temp$var2=log10(temp$var+min(temp$var[which(temp$var>0)]))
  temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
  temp$livecd3=livecd3.final[[1]][match(temp$event_id,names(livecd3.final[[1]]))] #get stim live count
  temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
  temp$responder=factor(temp$responder=="1")
  temp$logdays=log(as.numeric(temp$days))
  omics.responder.aves[1:2,k]=aggregate(temp$var,median,by=list(temp$responder),na.rm=T)[,2]
  
  temp=data.frame(cbind("days"=tcell.info.acute$event_date,"days_symptom"=tcell.info.acute$event_date_symptom,"gate"=unlist(tcell.final.acute[[i]][,j]),"id"=tcell.info.acute$participant_id,"event_id"=tcell.info.acute$event_id,"symptom_bin"=tcell.info.acute$symptom_date_bin,"responder"=tcell.info.acute$responder_c2,"PC1"=tcell.info.acute$PC1)) #get data
  if(omics.type[k]=="cytof")
  {
    temp$var=cytof.pct[match(temp$event_id,cytof$event_id),k]
  } else {
    temp$var=olink[match(temp$event_id,olink$event_id),k]
  }
  temp$PC1=as.numeric(temp$PC1)
  temp$var2=log10(temp$var+min(temp$var[which(temp$var>0)]))
  temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
  temp$livecd3=livecd3.final[[1]][match(temp$event_id,names(livecd3.final[[1]]))] #get stim live count
  temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
  temp$responder=factor(temp$responder=="1")
  temp$logdays=log(as.numeric(temp$days))
  omics.responder.aves[3:4,k]=aggregate(temp$var,median,by=list(temp$responder),na.rm=T)[,2]
  
  print(k)
}

#get pvalues to add to plot
omics.responder.pvals=data.frame(matrix(ncol=length(omics.vars),nrow=4))
omics.responder.pvals[1,]=omics.baseline.res.responder.log[,1]
omics.responder.pvals[2,]=omics.baseline.res.responder.log[,1]
omics.responder.pvals[3,]=omics.res.responder.log[,1]
omics.responder.pvals[4,]=omics.res.responder.log[,1]
omics.responder.ors=data.frame(matrix(ncol=length(omics.vars),nrow=4))
omics.responder.ors[1,]=omics.baseline.res.responder.log[,2]>1
omics.responder.ors[2,]=omics.baseline.res.responder.log[,2]>1
omics.responder.ors[3,]=omics.res.responder.log[,2]>1
omics.responder.ors[4,]=omics.res.responder.log[,2]>1

omics.responder.sigs=data.frame(matrix(ncol=length(omics.vars),nrow=4))
colnames(omics.responder.sigs)=omics.vars
rownames(omics.responder.sigs)=c("low_baseline","high_baseline","low","high")
omics.responder.sigs[1,omics.responder.pvals[1,]<=0.05 & omics.responder.ors[1,]==F]="*"
omics.responder.sigs[1,omics.responder.pvals[1,]<=0.05 & omics.responder.ors[1,]==T|omics.responder.pvals[1,]>0.05]=""
omics.responder.sigs[2,omics.responder.pvals[1,]<=0.05 & omics.responder.ors[1,]==T]="*"
omics.responder.sigs[2,omics.responder.pvals[1,]<=0.05 & omics.responder.ors[1,]==F|omics.responder.pvals[1,]>0.05]=""
omics.responder.sigs[3,omics.responder.pvals[3,]<=0.05 & omics.responder.ors[3,]==F]="*"
omics.responder.sigs[3,omics.responder.pvals[3,]<=0.05 & omics.responder.ors[3,]==T|omics.responder.pvals[3,]>0.05]=""
omics.responder.sigs[4,omics.responder.pvals[3,]<=0.05 & omics.responder.ors[3,]==T]="*"
omics.responder.sigs[4,omics.responder.pvals[3,]<=0.05 & omics.responder.ors[3,]==F|omics.responder.pvals[3,]>0.05]=""

omics.nams2=str_replace(gsub("[..]"," ",gsub("[ ].*","",omics.nams))," +"," ")
names(omics.nams2)=omics.vars
#plot

##redo with Tfh format
#correlate tcell populations with omics
#first with longitudinal data
tcell.omics.cor=vector("list",length(is))
names(tcell.omics.cor)=is
for(i in is)
{	
  tcell.omics.cor[[i]]=data.frame(matrix(nrow=length(omics.vars),ncol=length(js)))
  colnames(tcell.omics.cor[[i]])=js
  rownames(tcell.omics.cor[[i]])=omics.vars
}	
tcell.omics.cor.pv=tcell.omics.cor

for(i in is) #only run on the 3 SARS stimulations, not the positive control
{
  
  for(j in js)  #run on a few key gates
  {
    
    for(k in omics.vars)
    {
      tempp=tcell.info.acute[!is.na(tcell.final.acute[[i]][,j]),]
      temp=tcell.final.acute[[i]][!is.na(tcell.final.acute[[i]][,j]),]
      if(omics.type[k]=="cytof")
      {
        tempp$var=cytof.pct[match(tempp$event_id,cytof$event_id),k]
        tempp$var2=log2(tempp$var+min(tempp$var[which(tempp$var>0)]))
        tcell.omics.cor[[i]][k,j]=cor.test(temp[,j],tempp$var2,method="spearman")$estimate
        tcell.omics.cor.pv[[i]][k,j]=cor.test(temp[,j],tempp$var2,method="spearman")$p.value
      } else {
        tempp$var=olink[match(tempp$event_id,olink$event_id),k]
        tcell.omics.cor[[i]][k,j]=cor.test(temp[,j],tempp$var,method="spearman")$estimate
        tcell.omics.cor.pv[[i]][k,j]=cor.test(temp[,j],tempp$var,method="spearman")$p.value		
      }
    }
  }
  print(i)
}	

#then with baseline data
tcell.bomics.cor=vector("list",length(is))
names(tcell.bomics.cor)=is
for(i in is)
{	
  tcell.bomics.cor[[i]]=data.frame(matrix(nrow=length(omics.vars),ncol=length(js)))
  colnames(tcell.bomics.cor[[i]])=js
  rownames(tcell.bomics.cor[[i]])=omics.vars
}	
tcell.bomics.cor.pv=tcell.bomics.cor

for(i in is) #only run on the 3 SARS stimulations, not the positive control
{
  
  for(j in js)  #run on a few key gates
  {
    
    for(k in omics.vars)
    {
      tempp=tcell.info.acute[!is.na(tcell.final.acute[[i]][,j]),]
      temp=tcell.final.acute[[i]][!is.na(tcell.final.acute[[i]][,j]),]
      if(omics.type[k]=="cytof")
      {
        tempo=cytof[substr(cytof$event_id,10,10)=="1",]
        tempoc=cytof.pct[substr(cytof$event_id,10,10)=="1",]
        tempp$var=tempoc[match(tempp$participant_id,tempo$participant_id),k]
        tempp$var2=log2(tempp$var+min(tempp$var[which(tempp$var>0)]))
        tcell.bomics.cor[[i]][k,j]=cor.test(temp[,j],tempp$var2,method="spearman")$estimate
        tcell.bomics.cor.pv[[i]][k,j]=cor.test(temp[,j],tempp$var2,method="spearman")$p.value
      } else {
        tempo=olink[substr(olink$event_id,10,10)=="1",]
        tempp$var=tempo[match(tempp$participant_id,tempo$participant_id),k]	
        tcell.bomics.cor[[i]][k,j]=cor.test(temp[,j],tempp$var,method="spearman")$estimate
        tcell.bomics.cor.pv[[i]][k,j]=cor.test(temp[,j],tempp$var,method="spearman")$p.value		
      }
    }
  }
  print(i)
}	

#plot correlations with significance stars
#get pvalues to add to plot
tcell.omics.sigs=vector("list",length(is))
names(tcell.omics.sigs)=is
for(j in is)
{
  tcell.omics.sigs[[j]]=tcell.omics.cor.pv[[j]]
  for(i in 1:ncol(tcell.omics.cor.pv[[j]]))
  {
    tcell.omics.sigs[[j]][which(tcell.omics.cor.pv[[j]][,i]>(0.05/59)),i]=""
    tcell.omics.sigs[[j]][which(tcell.omics.cor.pv[[j]][,i]<=(0.05/59)),i]="*"
    tcell.omics.sigs[[j]][which(tcell.omics.cor.pv[[j]][,i]<=(0.01/59)),i]="**"
    tcell.omics.sigs[[j]][which(tcell.omics.cor.pv[[j]][,i]<=(0.001/59)),i]="***"
  }
}

tcell.bomics.sigs=vector("list",length(is))
names(tcell.bomics.sigs)=is
for(j in is)
{
  tcell.bomics.sigs[[j]]=tcell.bomics.cor.pv[[j]]
  for(i in 1:ncol(tcell.bomics.cor.pv[[j]]))
  {
    tcell.bomics.sigs[[j]][which(tcell.bomics.cor.pv[[j]][,i]>(0.05/59)),i]=""
    tcell.bomics.sigs[[j]][which(tcell.bomics.cor.pv[[j]][,i]<=(0.05/59)),i]="*"
    tcell.bomics.sigs[[j]][which(tcell.bomics.cor.pv[[j]][,i]<=(0.01/59)),i]="**"
    tcell.bomics.sigs[[j]][which(tcell.bomics.cor.pv[[j]][,i]<=(0.001/59)),i]="***"
  }
}

tcell.omics.sigs=vector("list",length(is))
names(tcell.omics.sigs)=is
for(j in is)
{
  tcell.omics.sigs[[j]]=tcell.omics.cor.pv[[j]]
  for(i in 1:ncol(tcell.omics.cor.pv[[j]]))
  {
    tcell.omics.sigs[[j]][which(tcell.omics.cor.pv[[j]][,i]>(0.05/1)),i]=""
    tcell.omics.sigs[[j]][which(tcell.omics.cor.pv[[j]][,i]<=(0.05/1)),i]="*"
    tcell.omics.sigs[[j]][which(tcell.omics.cor.pv[[j]][,i]<=(0.01/1)),i]="**"
    tcell.omics.sigs[[j]][which(tcell.omics.cor.pv[[j]][,i]<=(0.001/1)),i]="***"
  }
}

tcell.bomics.sigs=vector("list",length(is))
names(tcell.bomics.sigs)=is
for(j in is)
{
  tcell.bomics.sigs[[j]]=tcell.bomics.cor.pv[[j]]
  for(i in 1:ncol(tcell.bomics.cor.pv[[j]]))
  {
    tcell.bomics.sigs[[j]][which(tcell.bomics.cor.pv[[j]][,i]>(0.05/1)),i]=""
    tcell.bomics.sigs[[j]][which(tcell.bomics.cor.pv[[j]][,i]<=(0.05/1)),i]="*"
    tcell.bomics.sigs[[j]][which(tcell.bomics.cor.pv[[j]][,i]<=(0.01/1)),i]="**"
    tcell.bomics.sigs[[j]][which(tcell.bomics.cor.pv[[j]][,i]<=(0.001/1)),i]="***"
  }
}

##repeat but tcell>0 vs omics
#first with longitudinal data
tcell.omics.or=vector("list",length(is))
names(tcell.omics.or)=is
for(i in is)
{	
  tcell.omics.or[[i]]=data.frame(matrix(nrow=length(omics.vars),ncol=length(js)))
  colnames(tcell.omics.or[[i]])=js
  rownames(tcell.omics.or[[i]])=omics.vars
}	
tcell.omics.or.pv=tcell.omics.or

for(i in is) #only run on the 3 SARS stimulations, not the positive control
{
  
  for(j in js)  #run on a few key gates
  {
    
    for(k in omics.vars)
    {
      tempp=tcell.info.acute[!is.na(tcell.final.acute[[i]][,j]),]
      temp=tcell.final.acute[[i]][!is.na(tcell.final.acute[[i]][,j]),]
      tryCatch({
        if(omics.type[k]=="cytof")
        {
          tempp$var=cytof.pct[match(tempp$event_id,cytof$event_id),k]
          tempp$var2=log2(tempp$var+min(tempp$var[which(tempp$var>0)]))
          tempp$gate=temp[,j]>0
          m=glmer(tempp[!is.na(tempp$var2)&!is.na(tempp$gate),"gate"]>0~tempp[!is.na(tempp$var2)&!is.na(tempp$gate),"var2"]+(1|tempp[!is.na(tempp$var2)&!is.na(tempp$gate),"participant_id"]),family="binomial",control=glmerControl(optimizer = "nloptwrap"),nAGQ=20)
          tcell.omics.or.pv[[i]][k,j]=summary(m)$coefficients[2,4]
          tcell.omics.or[[i]][k,j]=exp(summary(m)$coefficients[2,1])
        } else {
          tempp$var=olink[match(tempp$event_id,olink$event_id),k]
          tempp$gate=temp[,j]>0
          m=glmer(tempp[!is.na(tempp$var)&!is.na(tempp$gate),"gate"]>0~tempp[!is.na(tempp$var)&!is.na(tempp$gate),"var"]+(1|tempp[!is.na(tempp$var)&!is.na(tempp$gate),"participant_id"]),family="binomial",control=glmerControl(optimizer = "nloptwrap"),nAGQ=20)
          tcell.omics.or.pv[[i]][k,j]=summary(m)$coefficients[2,4]
          tcell.omics.or[[i]][k,j]=exp(summary(m)$coefficients[2,1])	
        }
      },error=function(e){})
    }
  }
  print(i)
}	

#then with baseline data
tcell.bomics.or=vector("list",length(is))
names(tcell.bomics.or)=is
for(i in is)
{	
  tcell.bomics.or[[i]]=data.frame(matrix(nrow=length(omics.vars),ncol=length(js)))
  colnames(tcell.bomics.or[[i]])=js
  rownames(tcell.bomics.or[[i]])=omics.vars
}	
tcell.bomics.or.pv=tcell.bomics.or

for(i in is) #only run on the 3 SARS stimulations, not the positive control
{
  
  for(j in js)  #run on a few key gates
  {
    
    for(k in omics.vars)
    {
      tempp=tcell.info.acute[!is.na(tcell.final.acute[[i]][,j]),]
      temp=tcell.final.acute[[i]][!is.na(tcell.final.acute[[i]][,j]),]
      tryCatch({
        if(omics.type[k]=="cytof")
        {
          tempo=cytof[substr(cytof$event_id,10,10)=="1",]
          tempoc=cytof.pct[substr(cytof$event_id,10,10)=="1",]
          tempp$var=tempoc[match(tempp$participant_id,tempo$participant_id),k]
          tempp$var2=log2(tempp$var+min(tempp$var[which(tempp$var>0)]))
          tempp$gate=temp[,j]>0
          m=glmer(tempp[!is.na(tempp$var2)&!is.na(tempp$gate),"gate"]>0~tempp[!is.na(tempp$var2)&!is.na(tempp$gate),"var2"]+(1|tempp[!is.na(tempp$var2)&!is.na(tempp$gate),"participant_id"]),family="binomial",control=glmerControl(optimizer = "nloptwrap"),nAGQ=20)
          tcell.bomics.or.pv[[i]][k,j]=summary(m)$coefficients[2,4]
          tcell.bomics.or[[i]][k,j]=exp(summary(m)$coefficients[2,1])
        } else {
          tempo=olink[substr(olink$event_id,10,10)=="1",]
          tempp$var=tempo[match(tempp$participant_id,tempo$participant_id),k]	
          tempp$gate=temp[,j]>0
          m=glmer(tempp[!is.na(tempp$var)&!is.na(tempp$gate),"gate"]>0~tempp[!is.na(tempp$var)&!is.na(tempp$gate),"var"]+(1|tempp[!is.na(tempp$var)&!is.na(tempp$gate),"participant_id"]),family="binomial",control=glmerControl(optimizer = "nloptwrap"),nAGQ=20)
          tcell.bomics.or.pv[[i]][k,j]=summary(m)$coefficients[2,4]
          tcell.bomics.or[[i]][k,j]=exp(summary(m)$coefficients[2,1])	
        }
      },error=function(e){})
    }
  }
  print(i)
}		

#plot ORs with significance stars
#get pvalues to add to plot
tcell.omics.sigs=vector("list",length(is))
names(tcell.omics.sigs)=is
for(j in is)
{
  tcell.omics.sigs[[j]]=tcell.omics.or.pv[[j]]
  for(i in 1:ncol(tcell.omics.or.pv[[j]]))
  {
    tcell.omics.sigs[[j]][which(tcell.omics.or.pv[[j]][,i]>(0.05/59)),i]=""
    tcell.omics.sigs[[j]][which(tcell.omics.or.pv[[j]][,i]<=(0.05/59)),i]="*"
    tcell.omics.sigs[[j]][which(tcell.omics.or.pv[[j]][,i]<=(0.01/59)),i]="**"
    tcell.omics.sigs[[j]][which(tcell.omics.or.pv[[j]][,i]<=(0.001/59)),i]="***"
  }
}

tcell.bomics.sigs=vector("list",length(is))
names(tcell.bomics.sigs)=is
for(j in is)
{
  tcell.bomics.sigs[[j]]=tcell.bomics.or.pv[[j]]
  for(i in 1:ncol(tcell.bomics.or.pv[[j]]))
  {
    tcell.bomics.sigs[[j]][which(tcell.bomics.or.pv[[j]][,i]>(0.05/59)),i]=""
    tcell.bomics.sigs[[j]][which(tcell.bomics.or.pv[[j]][,i]<=(0.05/59)),i]="*"
    tcell.bomics.sigs[[j]][which(tcell.bomics.or.pv[[j]][,i]<=(0.01/59)),i]="**"
    tcell.bomics.sigs[[j]][which(tcell.bomics.or.pv[[j]][,i]<=(0.001/59)),i]="***"
  }
}

tcell.omics.sigs=vector("list",length(is))
names(tcell.omics.sigs)=is
for(j in is)
{
  tcell.omics.sigs[[j]]=tcell.omics.or.pv[[j]]
  for(i in 1:ncol(tcell.omics.or.pv[[j]]))
  {
    tcell.omics.sigs[[j]][which(tcell.omics.or.pv[[j]][,i]>(0.05/1)),i]=""
    tcell.omics.sigs[[j]][which(tcell.omics.or.pv[[j]][,i]<=(0.05/1)),i]="*"
    tcell.omics.sigs[[j]][which(tcell.omics.or.pv[[j]][,i]<=(0.01/1)),i]="**"
    tcell.omics.sigs[[j]][which(tcell.omics.or.pv[[j]][,i]<=(0.001/1)),i]="***"
  }
}

tcell.bomics.sigs=vector("list",length(is))
names(tcell.bomics.sigs)=is
for(j in is)
{
  tcell.bomics.sigs[[j]]=tcell.bomics.or.pv[[j]]
  for(i in 1:ncol(tcell.bomics.or.pv[[j]]))
  {
    tcell.bomics.sigs[[j]][which(tcell.bomics.or.pv[[j]][,i]>(0.05/1)),i]=""
    tcell.bomics.sigs[[j]][which(tcell.bomics.or.pv[[j]][,i]<=(0.05/1)),i]="*"
    tcell.bomics.sigs[[j]][which(tcell.bomics.or.pv[[j]][,i]<=(0.01/1)),i]="**"
    tcell.bomics.sigs[[j]][which(tcell.bomics.or.pv[[j]][,i]<=(0.001/1)),i]="***"
  }
}

###TFH

##read-in, process and integrate
tfh=read.delim("output/2024.06.04_merged_tfh.txt",header=T,sep="\t")
tfh$stim[tfh$stim=="SPIKE"]="SPK"
tfh$stim[tfh$stim=="CD4D"]="CD4"
tfh$stim[tfh$stim=="CD8D"]="CD8"

tfh.info=tfh[tfh$stim=="NC",1:2]
tfh.final=vector("list",5)
names(tfh.final)=c(names(tcell.final),"NC")
for(i in 1:5)
{
  tfh.final[[i]]=tfh[tfh$stim==names(tfh.final)[i],-c(1,4)]
  tfh.final[[i]]=tfh.final[[i]][tfh.info$event_id%in%tcell.info$event_id,]
}
tfh.info=tfh.info[tfh.info$event_id%in%tcell.info$event_id,]
tfh.info=tcell.info[match(tfh.info$event_id,tcell.info$event_id),]

#calculate % of each AIM gate that is Tfh
for(i in 1:4)
{
  tmp=tcell.raw.final[[i]][match(tfh.info$event_id,tcell.info$event_id),]
  tt=tfh.final[[i]][,17]/tmp[,3]*100
  tt[is.infinite(tt)]=NA
  tfh.final[[i]][,paste0("Tfh_of_",gsub("[_].*","",colnames(tfh.final[[i]])[17]))]=tt
  tt=tfh.final[[i]][,18]/tmp[,8]*100
  tt[is.infinite(tt)]=NA
  tfh.final[[i]][,paste0("Tfh_of_",gsub("[_].*","",colnames(tfh.final[[i]])[18]))]=tt
  tt=tfh.final[[i]][,19]/tmp[,10]*100
  tt[is.infinite(tt)]=NA
  tfh.final[[i]][,paste0("Tfh_of_",gsub("[_].*","",colnames(tfh.final[[i]])[19]))]=tt
  tt=tfh.final[[i]][,20]/tmp[,9]*100
  tt[is.infinite(tt)]=NA
  tfh.final[[i]][,paste0("Tfh_of_",gsub("[_].*","",colnames(tfh.final[[i]])[20]))]=tt
  tt=tfh.final[[i]][,21]/tmp[,15]*100
  tt[is.infinite(tt)]=NA
  tfh.final[[i]][,paste0("Tfh_of_",gsub("[_].*","",colnames(tfh.final[[i]])[21]))]=tt
}

#subset to early and acute data
tfh.final.acute=tfh.final
for(i in 1:5)
{
  tfh.final.acute[[i]]=tfh.final[[i]][tfh.info$event_id%in%tcell.info.acute$event_id,]
}
tfh.info.acute=tfh.info[tfh.info$event_id%in%tcell.info.acute$event_id,]
tfh.info.acute=tcell.info.acute[match(tfh.info.acute$event_id,tcell.info.acute$event_id),]

##compare Tfh in stimulations vs NC, see if actually changes
tfh.js=colnames(tfh.final[[1]])[c(4,16:24,15)]
tfh.js2=c("Tfh","PD1Tfh","CD25OX40Tfh","CD69CD40LTfh","CD137CD40LTfh","CD137OX40Tfh","OX40CD40LTfh","CD40LTfh","IL2Tfh","OX40Tfh","CD8Tfh")
names(tfh.js2)=tfh.js
tfh.js3=c("Tfh","PD1+Tfh","CD25+OX40+Tfh","CD69+CD40L+Tfh","CD137+CD40L+Tfh","CD137+OX40+Tfh","OX40+CD40L+Tfh","CD40L+Tfh","IL2+Tfh","OX40+Tfh","CD8+Tfh")
names(tfh.js3)=tfh.js

##change over time
comps=c("all","log(all)","symptom","log(symptom)")
comps.nams=c("Freq vs time","Freq vs log(time)","Freq vs symptom","Freq vs log(symptom)")

tfh.acute.time.pvalues=vector("list",length(is))
names(tfh.acute.time.pvalues)=is
for(i in 1:length(tfh.acute.time.pvalues))
{
  tfh.acute.time.pvalues[[i]]=data.frame(matrix(nrow=length(tfh.js),ncol=length(comps)))
  colnames(tfh.acute.time.pvalues[[i]])=comps
  rownames(tfh.acute.time.pvalues[[i]])=tfh.js
}

#also store coefficent (coef oe exp(coef))
tfh.acute.time.coefs=tfh.acute.time.pvalues

for(i in is) #only run on the 3 SARS stimulations, not the positive control
{
  
  for(j in tfh.js)  #run on a few key gates
  {
    
    temp=data.frame(cbind("days"=tfh.info.acute$event_date,"days_symptom"=tfh.info.acute$event_date_symptom,"gate"=unlist(tfh.final.acute[[i]][,j]),"id"=tfh.info.acute$participant_id,"event_id"=tfh.info.acute$event_id,"symptom_bin"=tfh.info.acute$symptom_date_bin)) #get data
    temp$nc_livecd3=nc.livecd3.final[match(temp$event_id,names(nc.livecd3.final))]  #get NC live count
    temp$livecd3=livecd3.final[[i]][match(temp$event_id,names(livecd3.final[[i]]))] #get stim live count
    temp=temp[which(temp$nc_livecd3>=100 & temp$livecd3>=100),] #remove samples with <100 NC/stim live cells
    temp$gate=as.numeric(temp$gate)
    temp$days=as.numeric(temp$days)
    temp$days_symptom=as.numeric(temp$days_symptom)
    
    temp$gate3=log2(temp$gate)  #as above but don't add minimum value
    temp$gate3[is.infinite(temp$gate3)]=NA
    
    temp$responder=temp$gate>0 #make responder
    
    temp$logdays=log(temp$days) #needed for lme
    temp$logdays_symptom=log(temp$days_symptom) #needed for lme
    
    m=gamm4(gate3~s(days,bs="cr"),data=temp[!is.na(temp$days),],random = ~(1|id),family="gaussian")
    tfh.acute.time.coefs[[i]][j,1]=summary(m$mer)$coefficients[2,1]
    tfh.acute.time.pvalues[[i]][j,1]=summary(m$gam)$s.pv
    m=gamm4(gate3~s(logdays,bs="cr"),data=temp[!is.na(temp$logdays),],random = ~(1|id),family="gaussian")
    tfh.acute.time.coefs[[i]][j,2]=summary(m$mer)$coefficients[2,1]
    tfh.acute.time.pvalues[[i]][j,2]=summary(m$gam)$s.pv
    m=gamm4(gate3~s(days_symptom,bs="cr"),data=temp[!is.na(temp$days_symptom),],random = ~(1|id),family="gaussian")
    tfh.acute.time.coefs[[i]][j,3]=summary(m$mer)$coefficients[2,1]
    tfh.acute.time.pvalues[[i]][j,3]=summary(m$gam)$s.pv
    m=gamm4(gate3~s(logdays_symptom,bs="cr"),data=temp[!is.na(temp$logdays_symptom),],random = ~(1|id),family="gaussian")
    tfh.acute.time.coefs[[i]][j,4]=summary(m$mer)$coefficients[2,1]
    tfh.acute.time.pvalues[[i]][j,4]=summary(m$gam)$s.pv					
  }
  
  print(i)
}

##change by outcome/time
#place to store results
tfh.acute.pvalues=vector("list",length(is))
names(tfh.acute.pvalues)=is
for(i in 1:length(tfh.acute.pvalues))
{
  tfh.acute.pvalues[[i]]=vector("list",length(comps))
  names(tfh.acute.pvalues[[i]])=comps
  for(j in 1:length(tfh.acute.pvalues[[i]]))
  {
    tfh.acute.pvalues[[i]][[j]]=data.frame(matrix(nrow=length(tfh.js),ncol=length(virvars)))
    colnames(tfh.acute.pvalues[[i]][[j]])=virvars
    rownames(tfh.acute.pvalues[[i]][[j]])=tfh.js
  }
}

#also store coefficent (coef oe exp(coef))
tfh.acute.coefs=tfh.acute.pvalues
tfh.acute.pvalues.int=tfh.acute.pvalues

#run 
##test for interaction of each clinical variable with time
for(i in is) #only run on the 3 SARS stimulations, not the positive control
{
  
  for(j in tfh.js)  #run on a few key gates
  {
    for(k in virvars)
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
      if(virmodel[k]!="lm")
      {
        temp$var=factor(temp$var)
        
        tryCatch({
          m=gamm(gate3~var+s(days,bs="cr",by=var),data=temp[!is.na(temp$days)&!is.na(temp$var)&!is.na(temp$gate3),],random =list(id=~1),family="gaussian")
          m1=gamm(gate3~var+s(days,bs="cr"),data=temp[!is.na(temp$days)&!is.na(temp$var)&!is.na(temp$gate3),],random =list(id=~1),family="gaussian")
          tfh.acute.pvalues[[i]][[1]][j,k]=anova(m1$gam,m$gam)$s.pv
          tfh.acute.coefs[[i]][[1]][j,k]=paste0(summary(m$gam)$p.table[-1,1],collapse="/")
          tfh.acute.pvalues.int[[i]][[1]][j,k]=paste0(summary(m$gam)$p.table[-1,4],collapse="/")
        },error=function(e){})
        
        tryCatch({
          m=gamm(gate3~var+s(logdays,bs="cr",by=var),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$gate3),],random =list(id=~1),family="gaussian")
          m1=gamm(gate3~var+s(logdays,bs="cr"),data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$gate3),],random =list(id=~1),family="gaussian")
          tfh.acute.pvalues[[i]][[2]][j,k]=anova(m1$gam,m$gam)$s.pv
          tfh.acute.coefs[[i]][[2]][j,k]=paste0(summary(m$gam)$p.table[-1,1],collapse="/")
          tfh.acute.pvalues.int[[i]][[2]][j,k]=paste0(summary(m$gam)$p.table[-1,4],collapse="/")
        },error=function(e){})
        
        tryCatch({
          m=gamm(gate3~var+s(days_symptom,bs="cr",by=var),data=temp[!is.na(temp$days_symptom)&!is.na(temp$var)&!is.na(temp$gate3),],random =list(id=~1),family="gaussian")
          m1=gamm(gate3~var+s(days_symptom,bs="cr"),data=temp[!is.na(temp$days_symptom)&!is.na(temp$var)&!is.na(temp$gate3),],random =list(id=~1),family="gaussian")
          tfh.acute.pvalues[[i]][[3]][j,k]=anova(m1$gam,m$gam)$s.pv
          tfh.acute.coefs[[i]][[3]][j,k]=paste0(summary(m$gam)$p.table[-1,1],collapse="/")
          tfh.acute.pvalues.int[[i]][[3]][j,k]=paste0(summary(m$gam)$p.table[-1,4],collapse="/")
        },error=function(e){})
        
        tryCatch({
          m=gamm(gate3~var+s(logdays_symptom,bs="cr",by=var),data=temp[!is.na(temp$logdays_symptom)&!is.na(temp$var)&!is.na(temp$gate3),],random =list(id=~1),family="gaussian")
          m1=gamm(gate3~var+s(logdays_symptom,bs="cr"),data=temp[!is.na(temp$logdays_symptom)&!is.na(temp$var)&!is.na(temp$gate3),],random =list(id=~1),family="gaussian")
          tfh.acute.pvalues[[i]][[4]][j,k]=anova(m1$gam,m$gam)$s.pv
          tfh.acute.coefs[[i]][[4]][j,k]=paste0(summary(m$gam)$p.table[-1,1],collapse="/")
          tfh.acute.pvalues.int[[i]][[4]][j,k]=paste0(summary(m$gam)$p.table[-1,4],collapse="/")
        },error=function(e){})
        
      } else {
        
        
        tryCatch({
          m=lme(gate3~var+days+days*var,data=temp[!is.na(temp$days)&!is.na(temp$var)&!is.na(temp$gate3),],random =~1|id)
          m1=lme(gate3~var+days,data=temp[!is.na(temp$days)&!is.na(temp$var)&!is.na(temp$gate3),],random =~1|id)
          tfh.acute.pvalues[[i]][[1]][j,k]=summary(m)$tTable[4,5]
          tfh.acute.coefs[[i]][[1]][j,k]=summary(m)$tTable[4,1]
          tfh.acute.pvalues.int[[i]][[1]][j,k]=summary(m)$tTable[2,5]
        },error=function(e){})
        
        tryCatch({
          m=lme(gate3~var+logdays+logdays*var,data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$gate3),],random =~1|id)
          m1=lme(gate3~var+logdays,data=temp[!is.na(temp$logdays)&!is.na(temp$var)&!is.na(temp$gate3),],random =~1|id)
          tfh.acute.pvalues[[i]][[2]][j,k]=summary(m)$tTable[4,5]
          tfh.acute.coefs[[i]][[2]][j,k]=summary(m)$tTable[4,1]
          tfh.acute.pvalues.int[[i]][[2]][j,k]=summary(m)$tTable[2,5]
        },error=function(e){})
        
        tryCatch({
          m=lme(gate3~var+days_symptom+days_symptom*var,data=temp[!is.na(temp$days_symptom)&!is.na(temp$var)&!is.na(temp$gate3),],random =~1|id)
          m1=lme(gate3~var+days_symptom,data=temp[!is.na(temp$days_symptom)&!is.na(temp$var)&!is.na(temp$gate3),],random =~1|id)
          tfh.acute.pvalues[[i]][[3]][j,k]=summary(m)$tTable[4,5]
          tfh.acute.coefs[[i]][[3]][j,k]=summary(m)$tTable[4,1]
          tfh.acute.pvalues.int[[i]][[3]][j,k]=summary(m)$tTable[2,5]
        },error=function(e){})
        
        tryCatch({
          m=lme(gate3~var+logdays_symptom+logdays_symptom*var,data=temp[!is.na(temp$logdays_symptom)&!is.na(temp$var)&!is.na(temp$gate3),],random =~1|id)
          m1=lme(gate3~var+logdays_symptom,data=temp[!is.na(temp$logdays_symptom)&!is.na(temp$var)&!is.na(temp$gate3),],random =~1|id)
          tfh.acute.pvalues[[i]][[4]][j,k]=summary(m)$tTable[4,5]
          tfh.acute.coefs[[i]][[4]][j,k]=summary(m)$tTable[4,1]
          tfh.acute.pvalues.int[[i]][[4]][j,k]=summary(m)$tTable[2,5]
        },error=function(e){})
        
      }
    }
    print(j)
  }
  print(i)
}	

##correlate Tfh populations with 1) total T cell response and 2) individual gates
tfh.final.acute.merge=lapply(tfh.final.acute,function(x){x[match(tcell.info.acute$event_id,tfh.info.acute$event_id),]})

#first with responder PC1
tfh.responder.cor=data.frame(matrix(nrow=length(is),ncol=length(tfh.js)))
rownames(tfh.responder.cor)=is
colnames(tfh.responder.cor)=tfh.js
tfh.responder.cor.pv=tfh.responder.cor

for(i in is) #only run on the 3 SARS stimulations, not the positive control
{
  
  for(j in tfh.js)  #run on a few key gates
  {
    tfh.responder.cor[i,j]=cor.test(tfh.final.acute.merge[[i]][,j],-tcell.info.acute$PC1,method="spearman")$estimate
    tfh.responder.cor.pv[i,j]=cor.test(tfh.final.acute.merge[[i]][,j],-tcell.info.acute$PC1,method="spearman")$p.value
  }
  print(i)
}	

#then with responder C2
tfh.responder.or=data.frame(matrix(nrow=length(is),ncol=length(tfh.js)))
rownames(tfh.responder.or)=is
colnames(tfh.responder.or)=tfh.js
tfh.responder.or.pv=tfh.responder.or

for(i in is) #only run on the 3 SARS stimulations, not the positive control
{
  
  for(j in tfh.js)  #run on a few key gates
  {
    m=glmer(relevel(tcell.info.acute$responder_c2[!is.na(tcell.info.acute$responder_c2)&!is.na(tfh.final.acute.merge[[i]][,j])],"low")~tfh.final.acute.merge[[i]][!is.na(tcell.info.acute$responder_c2)&!is.na(tfh.final.acute.merge[[i]][,j]),j]+(1|tcell.info.acute$participant_id[!is.na(tcell.info.acute$responder_c2)&!is.na(tfh.final.acute.merge[[i]][,j])]),family="binomial")
    tfh.responder.or.pv[i,j]=summary(m)$coefficients[2,4]
    tfh.responder.or[i,j]=exp(summary(m)$coefficients[2,1])
  }
  print(i)
}	##none significant

#then with individual gates
tfh.individual.cor=vector("list",length(is))
names(tfh.individual.cor)=is
for(i in 1:length(is))
{	
  tfh.individual.cor[[i]]=data.frame(matrix(nrow=length(js),ncol=length(tfh.js)))
  colnames(tfh.individual.cor[[i]])=tfh.js
  rownames(tfh.individual.cor[[i]])=js
}	
tfh.individual.cor.pv=tfh.individual.cor

for(i in is) #only run on the 3 SARS stimulations, not the positive control
{
  
  for(j in tfh.js)  #run on a few key gates
  {
    for(k in js)
    {
      tfh.individual.cor[[i]][k,j]=cor.test(tfh.final.acute.merge[[i]][,j],log2(tcell.final.acute[[i]][,k]),method="spearman")$estimate
      tfh.individual.cor.pv[[i]][k,j]=cor.test(tfh.final.acute.merge[[i]][,j],log2(tcell.final.acute[[i]][,k]),method="spearman")$p.value
    }
  }
  print(i)
}	

tfh.individual.or=vector("list",length(is))
names(tfh.individual.or)=is
for(i in 1:length(is))
{	
  tfh.individual.or[[i]]=data.frame(matrix(nrow=length(js),ncol=length(tfh.js)))
  colnames(tfh.individual.or[[i]])=tfh.js
  rownames(tfh.individual.or[[i]])=js
}	
tfh.individual.or.pv=tfh.individual.or

for(i in is) #only run on the 3 SARS stimulations, not the positive control
{
  
  for(j in tfh.js)  #run on a few key gates
  {
    for(k in js)
    {
      m=glmer(tcell.final.acute[[i]][!is.na(tcell.final.acute[[i]][,k])&!is.na(tfh.final.acute.merge[[i]][,j]),k]>0~tfh.final.acute.merge[[i]][!is.na(tcell.final.acute[[i]][,k])&!is.na(tfh.final.acute.merge[[i]][,j]),j]+(1|tcell.info.acute$participant_id[!is.na(tcell.final.acute[[i]][,k])&!is.na(tfh.final.acute.merge[[i]][,j])]),family="binomial")
      tfh.individual.or.pv[[i]][k,j]=summary(m)$coefficients[2,4]
      tfh.individual.or[[i]][k,j]=exp(summary(m)$coefficients[2,1])
    }
  }
  print(i)
}	

#plot correlations with significance stars
tfh.responder.pvals=rbind(tfh.responder.cor.pv[1,tfh.js[-1]],tfh.individual.cor.pv[[1]][,tfh.js[-1]])
tfh.responder.ors=rbind(tfh.responder.cor[1,tfh.js[-1]],tfh.individual.cor[[1]][,tfh.js[-1]])

tfh.responder.sigs=tfh.responder.pvals
for(i in 1:ncol(tfh.responder.pvals))
{
  tfh.responder.sigs[which(tfh.responder.pvals[,i]>(0.05/130)),i]=""
  tfh.responder.sigs[which(tfh.responder.pvals[,i]<=(0.05/130)),i]="*"
  tfh.responder.sigs[which(tfh.responder.pvals[,i]<=(0.01/130)),i]="**"
  tfh.responder.sigs[which(tfh.responder.pvals[,i]<=(0.001/130)),i]="***"
}

ryb=brewer.pal(11,"RdYlBu")

#correlate Tfh populations with omics
#first with longitudinal data
tfh.omics.cor=vector("list",length(is))
names(tfh.omics.cor)=is
for(i in is)
{	
  tfh.omics.cor[[i]]=data.frame(matrix(nrow=length(omics.vars),ncol=length(tfh.js)))
  colnames(tfh.omics.cor[[i]])=tfh.js
  rownames(tfh.omics.cor[[i]])=omics.vars
}	
tfh.omics.cor.pv=tfh.omics.cor

for(i in is) #only run on the 3 SARS stimulations, not the positive control
{
  
  for(j in tfh.js)  #run on a few key gates
  {
    
    for(k in omics.vars)
    {
      tempp=tcell.info.acute[!is.na(tfh.final.acute.merge[[i]][,j]),]
      temp=tfh.final.acute.merge[[i]][!is.na(tfh.final.acute.merge[[i]][,j]),]
      if(omics.type[k]=="cytof")
      {
        tempp$var=cytof.pct[match(tempp$event_id,cytof$event_id),k]
        tempp$var2=log2(tempp$var+min(tempp$var[which(tempp$var>0)]))
        tfh.omics.cor[[i]][k,j]=cor.test(temp[,j],tempp$var2,method="spearman")$estimate
        tfh.omics.cor.pv[[i]][k,j]=cor.test(temp[,j],tempp$var2,method="spearman")$p.value
      } else {
        tempp$var=olink[match(tempp$event_id,olink$event_id),k]
        tfh.omics.cor[[i]][k,j]=cor.test(temp[,j],tempp$var,method="spearman")$estimate
        tfh.omics.cor.pv[[i]][k,j]=cor.test(temp[,j],tempp$var,method="spearman")$p.value		
      }
    }
  }
  print(i)
}	

#then with baseline data
tfh.bomics.cor=vector("list",length(is))
names(tfh.bomics.cor)=is
for(i in is)
{	
  tfh.bomics.cor[[i]]=data.frame(matrix(nrow=length(omics.vars),ncol=length(tfh.js)))
  colnames(tfh.bomics.cor[[i]])=tfh.js
  rownames(tfh.bomics.cor[[i]])=omics.vars
}	
tfh.bomics.cor.pv=tfh.bomics.cor

for(i in is) #only run on the 3 SARS stimulations, not the positive control
{
  
  for(j in tfh.js)  #run on a few key gates
  {
    
    for(k in omics.vars)
    {
      tempp=tcell.info.acute[!is.na(tfh.final.acute.merge[[i]][,j]),]
      temp=tfh.final.acute.merge[[i]][!is.na(tfh.final.acute.merge[[i]][,j]),]
      if(omics.type[k]=="cytof")
      {
        tempo=cytof[substr(cytof$event_id,10,10)=="1",]
        tempoc=cytof.pct[substr(cytof$event_id,10,10)=="1",]
        tempp$var=tempoc[match(tempp$participant_id,tempo$participant_id),k]
        tempp$var2=log2(tempp$var+min(tempp$var[which(tempp$var>0)]))
        tfh.bomics.cor[[i]][k,j]=cor.test(temp[,j],tempp$var2,method="spearman")$estimate
        tfh.bomics.cor.pv[[i]][k,j]=cor.test(temp[,j],tempp$var2,method="spearman")$p.value
      } else {
        tempo=olink[substr(olink$event_id,10,10)=="1",]
        tempp$var=tempo[match(tempp$participant_id,tempo$participant_id),k]	
        tfh.bomics.cor[[i]][k,j]=cor.test(temp[,j],tempp$var,method="spearman")$estimate
        tfh.bomics.cor.pv[[i]][k,j]=cor.test(temp[,j],tempp$var,method="spearman")$p.value		
      }
    }
  }
  print(i)
}	


##repeat but tfh>median vs omics
#first with longitudinal data
tfh.omics.or=vector("list",length(is))
names(tfh.omics.or)=is
for(i in is)
{	
  tfh.omics.or[[i]]=data.frame(matrix(nrow=length(omics.vars),ncol=length(tfh.js)))
  colnames(tfh.omics.or[[i]])=tfh.js
  rownames(tfh.omics.or[[i]])=omics.vars
}	
tfh.omics.or.pv=tfh.omics.or

for(i in is) #only run on the 3 SARS stimulations, not the positive control
{
  
  for(j in tfh.js)  #run on a few key gates
  {
    
    for(k in omics.vars)
    {
      tempp=tcell.info.acute[!is.na(tfh.final.acute.merge[[i]][,j]),]
      temp=tfh.final.acute.merge[[i]][!is.na(tfh.final.acute.merge[[i]][,j]),]
      tryCatch({
        if(omics.type[k]=="cytof")
        {
          tempp$var=cytof.pct[match(tempp$event_id,cytof$event_id),k]
          tempp$var2=log2(tempp$var+min(tempp$var[which(tempp$var>0)]))
          tempp$gate=temp[,j]>median(temp[,j],na.rm=T)
          m=glmer(tempp[!is.na(tempp$var2)&!is.na(tempp$gate),"gate"]>0~tempp[!is.na(tempp$var2)&!is.na(tempp$gate),"var2"]+(1|tempp[!is.na(tempp$var2)&!is.na(tempp$gate),"participant_id"]),family="binomial",control=glmerControl(optimizer = "nloptwrap"),nAGQ=20)
          tfh.omics.or.pv[[i]][k,j]=summary(m)$coefficients[2,4]
          tfh.omics.or[[i]][k,j]=exp(summary(m)$coefficients[2,1])
        } else {
          tempp$var=olink[match(tempp$event_id,olink$event_id),k]
          tempp$gate=temp[,j]>median(temp[,j],na.rm=T)
          m=glmer(tempp[!is.na(tempp$var)&!is.na(tempp$gate),"gate"]>0~tempp[!is.na(tempp$var)&!is.na(tempp$gate),"var"]+(1|tempp[!is.na(tempp$var)&!is.na(tempp$gate),"participant_id"]),family="binomial",control=glmerControl(optimizer = "nloptwrap"),nAGQ=20)
          tfh.omics.or.pv[[i]][k,j]=summary(m)$coefficients[2,4]
          tfh.omics.or[[i]][k,j]=exp(summary(m)$coefficients[2,1])	
        }
      },error=function(e){})
    }
  }
  print(i)
}	

#then with baseline data
tfh.bomics.or=vector("list",length(is))
names(tfh.bomics.or)=is
for(i in is)
{	
  tfh.bomics.or[[i]]=data.frame(matrix(nrow=length(omics.vars),ncol=length(tfh.js)))
  colnames(tfh.bomics.or[[i]])=tfh.js
  rownames(tfh.bomics.or[[i]])=omics.vars
}	
tfh.bomics.or.pv=tfh.bomics.or

for(i in is) #only run on the 3 SARS stimulations, not the positive control
{
  
  for(j in tfh.js)  #run on a few key gates
  {
    
    for(k in omics.vars)
    {
      tempp=tcell.info.acute[!is.na(tfh.final.acute.merge[[i]][,j]),]
      temp=tfh.final.acute.merge[[i]][!is.na(tfh.final.acute.merge[[i]][,j]),]
      tryCatch({
        if(omics.type[k]=="cytof")
        {
          tempo=cytof[substr(cytof$event_id,10,10)=="1",]
          tempoc=cytof.pct[substr(cytof$event_id,10,10)=="1",]
          tempp$var=tempoc[match(tempp$participant_id,tempo$participant_id),k]
          tempp$var2=log2(tempp$var+min(tempp$var[which(tempp$var>0)]))
          tempp$gate=temp[,j]>median(temp[,j],na.rm=T)
          m=glmer(tempp[!is.na(tempp$var2)&!is.na(tempp$gate),"gate"]>0~tempp[!is.na(tempp$var2)&!is.na(tempp$gate),"var2"]+(1|tempp[!is.na(tempp$var2)&!is.na(tempp$gate),"participant_id"]),family="binomial",control=glmerControl(optimizer = "nloptwrap"),nAGQ=20)
          tfh.bomics.or.pv[[i]][k,j]=summary(m)$coefficients[2,4]
          tfh.bomics.or[[i]][k,j]=exp(summary(m)$coefficients[2,1])
        } else {
          tempo=olink[substr(olink$event_id,10,10)=="1",]
          tempp$var=tempo[match(tempp$participant_id,tempo$participant_id),k]	
          tempp$gate=temp[,j]>median(temp[,j],na.rm=T)
          m=glmer(tempp[!is.na(tempp$var)&!is.na(tempp$gate),"gate"]>0~tempp[!is.na(tempp$var)&!is.na(tempp$gate),"var"]+(1|tempp[!is.na(tempp$var)&!is.na(tempp$gate),"participant_id"]),family="binomial",,control=glmerControl(optimizer = "nloptwrap"),nAGQ=20)
          tfh.bomics.or.pv[[i]][k,j]=summary(m)$coefficients[2,4]
          tfh.bomics.or[[i]][k,j]=exp(summary(m)$coefficients[2,1])	
        }
      },error=function(e){})
    }
  }
  print(i)
}		

#plot ORs with significance stars
#get pvalues to add to plot
tfh.omics.sigs=vector("list",length(is))
names(tfh.omics.sigs)=is
for(j in is)
{
  tfh.omics.sigs[[j]]=tfh.omics.or.pv[[j]]
  for(i in 1:ncol(tfh.omics.or.pv[[j]]))
  {
    tfh.omics.sigs[[j]][which(tfh.omics.or.pv[[j]][,i]>(0.05/1)),i]=""
    tfh.omics.sigs[[j]][which(tfh.omics.or.pv[[j]][,i]<=(0.05/1)),i]="*"
    tfh.omics.sigs[[j]][which(tfh.omics.or.pv[[j]][,i]<=(0.01/1)),i]="**"
    tfh.omics.sigs[[j]][which(tfh.omics.or.pv[[j]][,i]<=(0.001/1)),i]="***"
  }
}

tfh.bomics.sigs=vector("list",length(is))
names(tfh.bomics.sigs)=is
for(j in is)
{
  tfh.bomics.sigs[[j]]=tfh.bomics.or.pv[[j]]
  for(i in 1:ncol(tfh.bomics.or.pv[[j]]))
  {
    tfh.bomics.sigs[[j]][which(tfh.bomics.or.pv[[j]][,i]>(0.05/1)),i]=""
    tfh.bomics.sigs[[j]][which(tfh.bomics.or.pv[[j]][,i]<=(0.05/1)),i]="*"
    tfh.bomics.sigs[[j]][which(tfh.bomics.or.pv[[j]][,i]<=(0.01/1)),i]="**"
    tfh.bomics.sigs[[j]][which(tfh.bomics.or.pv[[j]][,i]<=(0.001/1)),i]="***"
  }
}

##Analysis: unsupervised data
cytos[["sid"]]=gsub("D0","D1",paste0(cytos[["pid"]],"_",gsub("[-]","",cytos[["visit"]])))
xinfo=readRDS("2024.03.25_rajeshids_eventids.obj")
cytos[["event_id"]]=xinfo[[1]]$event_id[match(cytos[["sid"]],gsub("Day","D",xinfo[[1]]$sid))]

cytos[["days"]]=tcell.info$event_date[match(cytos[["event_id"]],tcell.info$event_id)]
cytos[["days_group"]]=NA
cytos[["days_group"]][which(cytos[["days"]]<4)]="1-3"
cytos[["days_group"]][which(cytos[["days"]]>=4)]="4-10"
cytos[["days_group"]][which(cytos[["days"]]>=11)]="11+"
cytos[["days_group"]]=factor(cytos[["days_group"]],c("1-3","4-10","11+"))