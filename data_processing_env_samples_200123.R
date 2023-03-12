########################################################
# DATA PROCESSING OF eDNA SAMPLES | JAMES ORD | 200123 #
########################################################
# reformatting / processing qPCR and ddPCR data (environmental samples)
# plotting estimated copy numbers
# prepping and outputting binary results dataframe for input into occupancy models
   
getwd()
rm(list=ls())
 
Sys.setenv(LANG = "en")
library("readxl")
library(data.table)
library(readr)
library(tidyr)
library(ggpubr)
library(gridExtra)
library(dplyr)     
`%not_in%`<-purrr::negate(`%in%`)
 
# Read in and organise data qPCR
   
# millex samples qpcr
samples_q <- read.csv("data/environmental/qPCR_data_merged_200123.csv")
samples_q$Cq<-as.numeric(samples_q$Cq)

head(samples_q)

qconcs<-NULL
for(i in levels(as.factor(samples_q$Run_ID))){
    qrun<-subset(samples_q,Run_ID==i) # need separate standard curve for each qPCR run
    standards<-subset(qrun,Sample_ID %like% "cop") # get the standards
    standards$copies<-parse_number(standards$Sample_ID)
    cmod<-lm(Cq~log(copies),data=standards) # model for standard curve
    coefs<-coef(cmod)
    qrun$modeled_conc<-exp((qrun$Cq-coefs[1])/coefs[2]) # model the unknown concentration from the standard curve model
    qconcs<-rbind(qconcs,qrun)
}

samples_q<-subset(qconcs,(Filter_Type=="Sterivex"|Filter_Type=="Millex"))[c(1,2,3,4,6,8)]

# see number of millex samples
nrow(subset(samples_q,Filter_Type=="Millex")[!duplicated(subset(samples_q,Filter_Type=="Millex")$Sample_ID),])
nrow(subset(samples_q,Filter_Type=="Sterivex")[!duplicated(subset(samples_q,Filter_Type=="Sterivex")$Sample_ID),])

# class as qPCR   
samples_q$pcr<-"qPCR"
# add total copies (concentration * 4ul input)
samples_q$tot_copies_mod<-samples_q$modeled_conc*4
   
# Now we have MODELED quantities (total copies) for all reps of all millex and sterivex samples
     
LOD_qPCR<-3.569*4
# visualise the Cq values vs modeled quant
thresholds_plot1<-ggplot(samples_q,aes(x=Cq,y=tot_copies_mod,col=Filter_Type,shape=River))+ # for total copies in 4ul input
    theme_bw()+
    geom_point()+
    geom_hline(yintercept=0.1,lty="dashed")+
    geom_hline(yintercept=LOD_qPCR,lty="dashed")+
    geom_hline(yintercept=1,lty="dashed")+
    labs(x="Cq",y="Estimated total copies",title="A")+
    scale_shape_manual(values=c(0,1,2,8,5,6))
   
# a closer look at the higher Cq replicates
thresholds_plot2<-ggplot(subset(samples_q,tot_copies_mod<10),aes(x=Cq,y=tot_copies_mod,col=Filter_Type,shape=River))+ # for total copies in 4ul input
    theme_bw()+
    geom_point()+
    geom_hline(yintercept=0.1,lty="dashed")+
    geom_hline(yintercept=LOD_qPCR,lty="dashed")+
    geom_hline(yintercept=1,lty="dashed")+
    labs(x="Cq",y="Estimated total copies",title="B")+
    scale_shape_manual(values=c(0,1,2,8,5,6))
   
# subset the 'questionable' samples - those with extremely low modeled quantities
questionable<-subset(samples_q,tot_copies_mod<0.1)
   
# Read in and organise data ddPCR
     
ddpcr_data<-read.csv("data/environmental/ddpcr_counts_final.csv",sep=";")
ddpcr_data<-subset(ddpcr_data,Sample %in% samples_q$Sample_ID)

# see number of sterivex (ddPCR) samples
nrow(ddpcr_data[!duplicated(ddpcr_data$Sample),])
   
# COMPARING SAMPLE AVERAGE ESTIMATED COPY NUMBERS

# get average counts per sample (ddPCR)
ddpcr_data_avg<-ddpcr_data %>% group_by(Sample) %>% 
  summarise(mean_count=mean(count))
   
colnames(ddpcr_data_avg)[1]<-"Sample_ID" # make column name match that for the qPCR dats
   
# get average counts per sample (qPCR)
samples_q$modeled_conc<-ifelse(is.na(samples_q$modeled_conc),0,samples_q$modeled_conc)
qconcs_avg<-samples_q%>%group_by(Sample_ID)%>%summarise(mean_q=mean(tot_copies_mod))
   
# merge ddPCR with qPCR average counts
merged_data<-merge(ddpcr_data_avg,qconcs_avg,by="Sample_ID")
   
# for the comparison we will only consider samples with at least 'some' quantity estimate from each of qPCR and ddPCR
# (i.e. at least one rep has some estimated copy number)
merged_data<-subset(merged_data,mean_q>0&mean_count>0)
   
# add the river variable
river_index<-samples_q[c(1,5)]
river_index<-river_index[!duplicated(river_index),]
merged_data <-merge(merged_data,river_index,by="Sample_ID")
   
# plot of correlation between log qPCR and ddPCR copy number estimates
cor_nol<-ggplot(merged_data,aes(x=log(mean_q),y=log(mean_count)))+
    theme_bw()+
    geom_point(aes(shape=River))+
    geom_smooth(method="lm")+
    stat_cor(label.y = 5)+
    stat_regline_equation(label.y = 4.5)+
    labs(x="log of estimated total copies from qPCR",
             y="log of estimated total copies from ddPCR",
             title="A")+
    xlim(-3.5,4.5)+ylim(-3.5,5.5)+
    scale_shape_manual(values=c(0,1,2,8,5,6))
   
# plotting differences in concentration estimate
     
# put all replicate-level copy number estimates into one dataframe
quantscomp_df<-data.frame(Sample_ID=ddpcr_data$Sample,tot_copies_mod=ddpcr_data$count,pcr="ddPCR")

quantscomp_df<-rbind(quantscomp_df,samples_q[c(5,7:8)])

# again get only samples with some copy number estimate in both qPCR and ddPCR
quantscomp_df<-subset(quantscomp_df,Sample_ID %in% merged_data$Sample_ID)
# re-add River
river_index<-samples_q[c(1,5)]
river_index<-river_index[!duplicated(river_index),]
quantscomp_df <-merge(quantscomp_df,river_index,by="Sample_ID")
quantscomp_df<-subset(quantscomp_df,Sample_ID%in%merged_data$Sample_ID)  
# then rank by average copy number
merged_data$mean_count_q<-(merged_data$mean_count+merged_data$mean_q)/2
sorted_data<-merged_data[order(merged_data$mean_count_q),]
sorted_data$rank<-seq(1:nrow(sorted_data))
rank_index<-sorted_data[c(1,6)]
quantscomp_df<-merge(quantscomp_df,rank_index,by="Sample_ID")
   
# Then plotting...
     
quantscomp_df$pcr<-factor(quantscomp_df$pcr,levels=c("qPCR","ddPCR"))
   
quantplot_ranks<-ggplot(quantscomp_df,aes(x=rank,y=tot_copies_mod,col=pcr,shape=River))+
    theme_bw()+
    geom_pointrange(stat="summary",fun.min = min,
                        fun.max = max,
                        fun = mean,position=position_dodge(0.5))+
    scale_shape_manual(values=c(0,1,2,8,5,6))+
    labs(x="Rank of mean total copy estimate",y="Estimated total copies\n(sample mean, minimum and maximum)",
             title="B")+
    scale_color_manual(name="PCR method",values=c("black","darkred"))+
    geom_hline(yintercept=4.27,lty="dashed",col="black")+
    geom_hline(yintercept=3.66,lty="dashed",col="darkred")
   
# PREPARING DATA FOR OCCUPANCY MODEL INPUT
# prepare binary results dataframe
     
# ddPCR data
samples_dd <- data.frame(Sample_ID=ddpcr_data$Sample,Positives=ddpcr_data$count)
# add sample metadata. these are the same for the sterivex samples:
dd_sample_meta<-subset(samples_q,Filter_Type=="Sterivex")[c(1:3,5)]
dd_sample_meta<-dd_sample_meta[!duplicated(dd_sample_meta),]
samples_dd<-merge(samples_dd,dd_sample_meta,by="Sample_ID")
samples_dd$Result<-ifelse(samples_dd$Positives>=1,1,0)
samples_dd$Positives<-NULL
samples_dd$pcr<-"ddPCR"
   
# combine with the qPCR data
samples_q$Result0<-ifelse(is.na(samples_q$Cq),0,1)
# set Cq cutoff of 40
samples_q$Result <-ifelse(samples_q$Cq>40&samples_q$Result0==1,0,samples_q$Result0)
samples_q$Result0<-NULL
samples_q2<-samples_q[c(1,2,3,5,7,9)]
samples_all<-rbind(samples_q2,samples_dd)
   
# a little bit of pre-sanity-checking before we write out the final thing
# proportion of positive replicates and likewise for only positive samples
     
# get only samples with positive reps
pos_samps<-samples_all%>%group_by(Sample_ID)%>%summarise(sum_pos=sum(Result))%>%subset(sum_pos>0)
   
# proportion of positive replicates
pos_reps_summary1<-samples_all%>%group_by(Filter_Type,pcr)%>%summarise(sum_pos=sum(Result),n=n())
pos_reps_summary1$frac_pos<-pos_reps_summary1$sum_pos/pos_reps_summary1$n
   
# proportion of positive replicates considering only positive samples
pos_reps_summary2<-subset(samples_all,Sample_ID%in%pos_samps$Sample_ID)%>%group_by(Filter_Type,pcr)%>%summarise(sum_pos=sum(Result),n=n())
pos_reps_summary2$frac_pos<-pos_reps_summary2$sum_pos/pos_reps_summary2$n
# write the above as a table
# write.table(pos_reps_summary2,file="data_processing/prop_posreps_of_possamps_200123.txt",sep="\t",row.names = F,quote=F)

# Additional (16/01/23): proportions of samples with positive reps
pos_samps_summary1<-samples_all%>%group_by(Sample_ID,Filter_Type,pcr)%>%summarise(sum_pos=sum(Result))
pos_samps_summary1$pos_samp<-ifelse(pos_samps_summary1$sum_pos>0,1,0)
pos_samps_summary2<-pos_samps_summary1%>%group_by(Filter_Type,pcr)%>%summarise(sum_pos=sum(pos_samp),nsamps=n())
pos_samps_summary2$prop_pos<-pos_samps_summary2$sum_pos/pos_samps_summary2$nsamps
# write the above as a table
# write.table(pos_samps_summary2,file="data_processing/prop_possamps_of_totsamps_200123.txt",sep="\t",row.names = F,quote=F)
   
# Write out the final binary results dataframe (input for occ. models)
# write.table(samples_all,file="data/environmental/detections_qpcr_ddpcr_filtered_200123.txt",sep="\t",row.names = F,quote=F)

# print figures
# png(filename = "data_processing/cor_plus_quantplot_ranks_200123.png", width = 10000, height = 3500,res=600)
# grid.arrange(cor_nol,quantplot_ranks,ncol=2)
# dev.off()
# 
# png(filename = "data_processing/threshold_plots_200123.png", width = 10000, height = 3500,res=600)
# grid.arrange(thresholds_plot1,thresholds_plot2,ncol=2)
# dev.off()