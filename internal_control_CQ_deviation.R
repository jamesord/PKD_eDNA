# Calculating deviation in sample mean internal control Cq from that of the mean of the standards
# James Ord 25/08/23

getwd()
rm(list=ls())

Sys.setenv(LANG = "en")
library("readxl")
library(data.table)
library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)
`%not_in%`<-purrr::negate(`%in%`)
`%not_like%`<-purrr::negate(`%like%`)

# read in file including Millex, Steri, standards.
samples <- read_excel("data/environmental/qPCR_data_final_changed.xlsx")

# get the mean internal control Cq for each run
# we want columns 13 and 15
standard_means<-subset(samples,Filter_Type!="Sterivex"|Filter_Type!="Millex")[c(13,15)] %>%
  group_by(Run_ID)%>%summarise(runmean_std_ICCQ=mean(`IC Cq`))

# subset only environmental samples
samples<-subset(samples,Filter_Type=="Sterivex"|Filter_Type=="Millex")
# add ID of the reference site
samples$River<-ifelse(is.na(samples$River),"REF_FIWI",samples$River)
# create new sample ID
samples$Sample_ID<-paste(samples$River,samples$Date,samples$`Field Rep`,samples$Filter_Type,sep="_")
# get only desired columns
samples<-samples[c(17,15,13,3,6)]

# merge samples with standard means, by run ID
samples<-merge(samples,standard_means,by="Run_ID")
# get the difference between the IC Cq of the sample and the mean IC Cq of the run
samples$IC_CqDiff<-samples$`IC Cq`-samples$runmean_std_ICCQ
# then finally get the mean of the difference for each sample
sample_means<-samples%>%group_by(Sample_ID,River,Filter_Type)%>%summarise(meanIC_CqDiff=mean(IC_CqDiff))
# max, min (i.e. 'maximum negative'), and mean deviation in Cq
max(sample_means$meanIC_CqDiff) # most devation of sample in positive direction = 2.82 (AA, Millex)
min(sample_means$meanIC_CqDiff) # most devation of sample in negative direction = -0.88 (AA, Millex)
mean(sample_means$meanIC_CqDiff) # mean deviation in sample average Cq = -0.07

# make plot
mean_deviations_plot<-ggplot(sample_means,aes(x=Filter_Type,y=meanIC_CqDiff))+
  geom_point(position=position_jitter(0.1))+
  geom_hline(yintercept=c(3,1),lty="dashed")+
  facet_wrap(.~River,ncol=4)+
  labs(x="Filter type",y="Internal control CQ (sample mean) \u2013\nInternal control CQ (standard mean)")

# print plot
png(filename = "data_processing/IC_CQ_deviations_bysample_250823.png", width = 5000, height = 2500,res=600)
mean_deviations_plot
dev.off()
