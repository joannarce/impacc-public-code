#lib
library(stringr)
library(tidyr)
library(dplyr)
library(plyr)
library(ggplot2)
library(knitr)
library(readxl)
library(reshape2)
library(tables)
library(doBy)
library(lcmm)
#library(xlsx)
library(readxl)
library(CluMix)
library(cluster)
library(fpc)
library(gplots)
library(tsne)


#trajectory data
dd1<-read_excel(path ="PROMIS_membership_06072022.xlsx", sheet="membership")
#remove symptoms
#dd1<-select(dd1,-c("class_organ","class_symptom"))
#mean and standard deviation for impact and mental
#impact
impact<-read_excel(path="PROMIS_scores_05272022.xlsx", sheet="Impact")
impact<-impact[,c("PID","PROMIS_Impact")]
impact_sum<-impact %>% group_by(PID) %>% summarise_each(funs(mean, sd))
#if only one dp, sd=0
#impact_sum[is.na(impact_sum)]<-0
colnames(impact_sum)[2:3]<-c("impact_mean","impact_sd")

#mental
mental<-read_excel(path="LCMM/PROMIS_scores_05272022.xlsx", sheet="Mental")
mental<-mental[,c("PID","PROMIS_Mental")]
mental_sum<-mental %>% group_by(PID) %>% summarise_each(funs(mean, sd))
#if only one dp, sd=0
#mental_sum[is.na(mental_sum)]<-0
colnames(mental_sum)[2:3]<-c("mental_mean","mental_sd")

#PROMIS physical
physical<-read_excel(path="PROMIS_scores_05272022.xlsx", sheet="Physical")
physical<-physical[,c("PID","PROMIS_Physical")]
physical_sum<-physical %>% group_by(PID) %>% summarise_each(funs(mean, sd))
colnames(physical_sum)[2:3]<-c("physical_mean","physical_sd")


#merge with dd1
dd1<-join_all(list(dd1,mental_sum,impact_sum,physical_sum), by = 'PID', type = 'full')
#remove class_physical
dd1<-subset(dd1, select = -class_physical)
dd1_noid<-subset(dd1, select = -PID) 
#try recode missing as another category?

dd1_noid$class_health<-ifelse(is.na(dd1_noid$class_health),2,dd1_noid$class_health)
#Gower
ds <- dist.subjects(dd1_noid) #do we need to weight the trajectory variables?

cl <- tibble(kk = 2:10,
             ds = list(ds))


library(purrr)
cl <- cl %>%
  mutate(
    pam=map2(kk, ds, ~ pam(x=.y, k=.x, diss=T, cluster.only=T)),
    mcquitty=map2(kk, ds,~ disthclustCBI(.y, .x, cut="number", method='mcquitty', scaling=F)$partition),
    complete= map2(kk, ds,
                   ~ disthclustCBI(.y, .x, cut="number", method='complete', scaling=F)$partition),
    average=map2(kk, ds,
                 ~ disthclustCBI(.y, .x, cut="number", method='average', scaling=F)$partition),
    ward=map2(kk, ds,
              ~ disthclustCBI(.y, .x, cut="number", method='ward.D2', scaling=F)$partition)
  )


cl1 <- cl %>%
  gather(cluster, partition, -c(1:2)) %>% 
  mutate(stats = map2(ds, partition, ~ cluster.stats(.x, .y))) %>% 
  mutate(within.cluster.ss = map_dbl(stats,"within.cluster.ss"),avg.silwidth = map_dbl(stats,"avg.silwidth"), dunn= map_dbl(stats,"dunn"),dunn2= map_dbl(stats,"dunn2"), wb.ratio= map_dbl(stats, 'wb.ratio'))
cl1_to_show<-cl1[,c(1,3,6:10)]

#final plots

cl1_to_show$sil<-"avg.silwidth"
p1<-ggplot(data = cl1_to_show, aes(x =kk, y = avg.silwidth,group=cluster))+geom_line(size=1)+geom_point(size=2)+
  scale_x_discrete(limits=c(2,4,6,8,10))+ theme(panel.spacing= unit(0.5, "in"),text = element_text(size = 20))+
  xlab("")+ylab("")+theme_bw(base_rect_size =1,base_size=15,base_line_size =1)+facet_wrap(~sil)

cl1_to_show$d<-"dunn"
p2<-ggplot(data = cl1_to_show, aes(x =kk, y = dunn2,group=cluster))+geom_line(size=1)+geom_point(size=2)+
  scale_x_discrete(limits=c(2,4,6,8,10))+ theme(panel.spacing= unit(0.5, "in"),text = element_text(size = 20))+
  xlab("")+ylab("")+theme_bw(base_rect_size =1,base_size=15,base_line_size =1)+facet_wrap(~d)

cl1_to_show$wb<-"wb.ratio"
p3<-ggplot(data = cl1_to_show, aes(x =kk, y = wb.ratio,group=cluster))+geom_line(size=1)+geom_point(size=2)+
  scale_x_discrete(limits=c(2,4,6,8,10))+ theme(panel.spacing= unit(0.5, "in"),text = element_text(size = 20))+
  xlab("")+ylab("")+theme_bw(base_rect_size =1,base_size=15,base_line_size =1)+facet_wrap(~wb)


cl1_to_show$wcss<-"within.cluster.ss"
p4<-ggplot(data = cl1_to_show, aes(x =kk, y = within.cluster.ss,group=cluster))+geom_line(size=1)+geom_point(size=2)+
  scale_x_discrete(limits=c(2,4,6,8,10))+ theme(panel.spacing= unit(0.5, "in"),text = element_text(size = 20))+
  xlab("")+ylab("")+theme_bw(base_rect_size =1,base_size=15,base_line_size =1)+facet_wrap(~wcss)

#panel of 4
library(patchwork)
nested <- (p1|p2|p3|p4)
nested

