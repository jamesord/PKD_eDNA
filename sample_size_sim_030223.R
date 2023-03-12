############################################################
# Sample size plot from occ model | James Ord | 03/02/2023 #
############################################################

# Here we run sample size simulations based on a multilevel occupancy model, calculating the number of samples required to detect PKD eDNA.

Sys.setLanguage("en")
library(msocc)
library(boot)
library("readxl")
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(scales)
`%not_in%`<-purrr::negate(`%in%`)

# load in and preprocess the data as before
alldata<-read.table("data/environmental/detections_qpcr_ddpcr_filtered_200123.txt",sep="\t",header=T)

head(alldata)

# add replicate IDs
repnums<-NULL
for (i in levels(as.factor(alldata$Sample_ID))){
  sample<-subset(alldata,Sample_ID==i)
  sample$rep.num<-seq(1:nrow(sample))
  repnums<-rbind(repnums,sample)
}
alldata<-repnums
# add sample IDs
sampnums<-NULL
for (i in levels(as.factor(alldata$River))){
  river<-subset(alldata,River==i)[c(1,4)]
  river<-river[!duplicated(river),]
  river$samp.num<-seq(1:nrow(river))
  sampnums<-rbind(sampnums,river)
}
alldata<-merge(alldata,sampnums[c(2,3)],by="Sample_ID")

# convert data to wide format
wide_data<-alldata[-c(3,5)] %>% pivot_wider( 
  names_from = rep.num,
  values_from = Result)
colnames(wide_data)[c(2,4,3)]<-c("site","sample","volume")

# prepare detections dataframe and covariate tables
detections<-wide_data[-c(1,3)]
samp_covs<-wide_data[c(2,4,3)]
rep_covs<-alldata[c(2,8,7,3,5)];colnames(rep_covs)<-c("site","sample","rep","filter","pcr_method")
site_covs<-wide_data[c(2,4)]

#############
# FIT MODEL #
#############

# the only model we care about this time (model 2 from before)
mod2 <- msocc_mod(wide_data = detections,
                  site = list(model = ~ 1, cov_tbl = site_covs),
                  sample = list(model = ~ volume, cov_tbl = samp_covs),
                  rep = list(model = ~ filter, cov_tbl = rep_covs),
                  progress = T,num.mcmc=11000,seed = 1598)

# As we deemed there not to be significant effects of pcr method or filter type, we omit this covariate.
# Rather, we will run our simulations based on the effect of sample volume.

###########################
# SAMPLE SIZE SIMULATIONS #
###########################

# extract posterior samples of delta parameters from the model (last 10k samples = 1k burn-in), and convert to p estimates for millex and sterivex
alphas<-tail(as.data.frame(mod2$alpha),10000)
deltas<-tail(as.data.frame(mod2$delta),10000)

samples_est_all<-data.frame(theta=c(inv.logit(alphas$`(Intercept)`+alphas$volume*300),
                                    inv.logit(alphas$`(Intercept)`+alphas$volume*600)),
                            volume=c(rep("300",10000),rep("600",10000)),
                            p=rep(inv.logit(deltas$`(Intercept)`),2))

# For scenarios of collecting 1 to 10 samples, we calculate cumulative probabilities of obtaining at least one 'positive sample', given N samples.
# 1-((1-prob_positive)^N) where prob_positive is the probability of any given sample being positive and N is number of samples
# prob_positive is the capture probability (theta) multiplied by the detection probability, p.
# However, we have multiple PCR replicates and we may require multiple reps to be positive before we deem the sample positive.
# Therefore, the cumulative probability is calculated for three 'positive sample' criteria:
# at least 1/3 PCR reps, at least 2/3 PCR reps, and 3/3 PCR reps,
sampNsim<-NULL
for (i in seq(1:10)){
  samples_est_N<-samples_est_all  
  # detection criterion
  samples_est_N_criteria<-NULL
  # positive criterion: at least 1/3 PCR reps (additive probability of 1/3, 2/3, and 3/3 reps)
  samples_est_N$cumul_prob1<- 1-((1-(samples_est_N$theta * (dbinom(1,3,samples_est_N$p)+dbinom(2,3,samples_est_N$p)+dbinom(3,3,samples_est_N$p)) ))^i)
  # positive criterion: at least 2/3 PCR reps (additive probability of 2/3, and 3/3 reps)
  samples_est_N$cumul_prob2<- 1-((1-(samples_est_N$theta * (dbinom(2,3,samples_est_N$p)+dbinom(3,3,samples_est_N$p)) ))^i)
  # positive criterion: 3/3 PCR reps (only the probability of 3/3 reps)
  samples_est_N$cumul_prob3<- 1-((1-(samples_est_N$theta * (dbinom(3,3,samples_est_N$p)) ))^i)
  samples_est_N$N_samples<-i
  sampNsim<-rbind(sampNsim,samples_est_N)}

cpr1<-sampNsim[-c(5,6)];colnames(cpr1)[4]<-"cumul_prob";cpr1$criterion="1 of 3 PCR replicates"
cpr2<-sampNsim[-c(4,6)];colnames(cpr2)[4]<-"cumul_prob";cpr2$criterion="2 of 3 PCR replicates"
cpr3<-sampNsim[-c(4,5)];colnames(cpr3)[4]<-"cumul_prob";cpr3$criterion="3 of 3 PCR replicates"

sampNsim<-rbind(cpr1,cpr2,cpr3)

# The above calculations were performed for all samples in the posterior distribution of the model.
# We can therefore represent the estimates as medians and credible intervals.

probs_toplotl<-sampNsim %>%
  group_by(volume,N_samples,criterion) %>%
  summarise(med.prob = median(cumul_prob, na.rm = TRUE),
            lower=quantile(cumul_prob,probs=c(.025,.975))[1],
            upper=quantile(cumul_prob,probs=c(.025,.975))[2])

# plot code
sample_size_plotl<-ggplot(probs_toplotl, aes(x=N_samples, y=med.prob,group=volume)) +
  theme_bw()+
  geom_errorbar(aes(ymin=lower, ymax=upper),
                width=0,position=position_dodge(0.5),colour="black")+
  geom_line(aes(lty=volume),position=position_dodge(0.5))+
  geom_point(size=3, position=position_dodge(width=0.5),aes(fill=volume),shape=21)+
  labs(x="Number of samples",
       y = "Cumulative probability of positive sample\n(median +/- 95% credible intervals)")+
  geom_hline(yintercept=0.95,lty="dashed")+
  theme(text = element_text(size = 15))+
  scale_x_continuous(breaks= pretty_breaks())+
  #theme(legend.position = "top")+
  facet_wrap(.~criterion,ncol=3)+
  theme(legend.position = "right")+
  scale_fill_manual(name = "Filtered water\nvolume (ml):",values=c("white","black"))+
  scale_linetype_manual(name = "Filtered water\nvolume (ml):",values=c("dashed","solid"))

png(filename = "sample_size_plot_120323.png", width = 7000, height = 2500,res=600)
sample_size_plotl
dev.off()

# APPENDIX

# EXAMPLE OF PROBABILITY CALCULATION WITH REP-LEVEL DETECTION PROBABILITY = 0.92
# ASSUMING SAMPLE IS POSITIVE

# Probability of exactly 1 positive rep, given 3 reps:
dbinom(1,3,0.92) # = 0.017
# Probability of exactly 2 positive reps, given 3 reps:
dbinom(2,3,0.92) # = 0.2
# Probability of exactly 3 positive reps, given 3 reps:
dbinom(3,3,0.92) # = 0.78

# To get the probability of 'at least N', we sum up the probability of each permissible outcome...

# Probability of at least 1 success = sum of three permissible outcomes:
dbinom(1,3,0.92)+dbinom(2,3,0.92)+dbinom(3,3,0.92)# = 0.99
# Probability of at least 2 successes = sum of two permissible outcomes:
dbinom(2,3,0.92)+dbinom(3,3,0.92)# = 0.98
# Probability of all three successes = probability of only one permissible outcome:
dbinom(3,3,0.92)# = 0.78
