##########################################################
# Occupancy models with 'msocc' | James Ord | 20/01/2023 #
##########################################################

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

###########################
# define custom functions #
###########################

# modified lppd function to include burn-in option
compute_lppd_burnin <- function(msocc_mod, burnin){
  #pull model info
  num.mcmc <- msocc_mod$model.info$num.mcmc-burnin
  J <- msocc_mod$model.info$J
  K <- msocc_mod$model.info$K
  z <- tail(msocc_mod$model.info$z.vec,num.mcmc)[, rep(1:ncol(msocc_mod$model.info$z.vec), K)]
  A <- tail(msocc_mod$model.info$A,num.mcmc)[, rep(1:ncol(msocc_mod$model.info$A), K)]
  y <- matrix(rep(msocc_mod$model.info$y, num.mcmc), nrow = num.mcmc, byrow = T)
  #convert MCMC samples to probability
  psi <- tail(psi_mcmc(msocc_mod),num.mcmc)
  theta <- tail(theta_mcmc(msocc_mod),num.mcmc)
  p <- tail(p_mcmc(msocc_mod),num.mcmc)
  site <- psi^z * (1 - psi)^(1 - z)
  sample <- (z*theta)^A * (1 - z*theta)^(1 - A)
  rep <- (A*p)^y * (1 - A*p)^(1 - y)
  tmp <- (site * sample * rep)[2:num.mcmc,]
  return(sum(log(colMeans(tmp))))
}

# modified pwaic2 function to include burn-in option
compute_pwaic2_burnin <- function(msocc_mod, burnin){
  #pull model info
  num.mcmc <- msocc_mod$model.info$num.mcmc-burnin
  J <- msocc_mod$model.info$J
  K <- msocc_mod$model.info$K
  z <- tail(msocc_mod$model.info$z.vec,num.mcmc)[, rep(1:ncol(msocc_mod$model.info$z.vec), K)]
  A <- tail(msocc_mod$model.info$A,num.mcmc)[, rep(1:ncol(msocc_mod$model.info$A), K)]
  y <- matrix(rep(msocc_mod$model.info$y, num.mcmc), nrow = num.mcmc, byrow = T)
  #first piece of the difference
  #convert MCMC samples to probability
  psi <- tail(psi_mcmc(msocc_mod),num.mcmc)
  theta <- tail(theta_mcmc(msocc_mod),num.mcmc)
  p <- tail(p_mcmc(msocc_mod),num.mcmc)
  site <- psi^z * (1 - psi)^(1 - z)
  sample <- (z*theta)^A * (1 - z*theta)^(1 - A)
  rep <- (A*p)^y * (1 - A*p)^(1 - y)
  tmp <- log(site * sample * rep)[2:num.mcmc,]
  pwaic2 <- sum(apply(tmp, 2, var))
  return(pwaic2)
}

# modified waic function (with type 2 penalty score) to include burn-in option
waic2_burnin <- function(msocc_mod, burnin){
  waic <- -2 * (compute_lppd_burnin(msocc_mod,burnin) - compute_pwaic2_burnin(msocc_mod,burnin))
  return(waic)
}

# function to get a table of waic and associated parameters from model fit
get_waic_table<-function(msocc_mod,burnin){
  formula<-paste("psi(",msocc_mod$model.info$site_mod[2],"),theta(",msocc_mod$model.info$samp_mod[2],"),p(",msocc_mod$model.info$rep_mod[2],")",sep="")
  waic<-waic2_burnin(msocc_mod,burnin)
  pwaic2<-compute_pwaic2_burnin(msocc_mod,burnin)
  lppd<-compute_lppd_burnin(msocc_mod,burnin)
  num.mcmc <- msocc_mod$model.info$num.mcmc
  return(data.frame(formula=formula,waic=waic,pwaic2=pwaic2,lppd=lppd,iter=num.mcmc,burnin=burnin))
}

########################
# LOAD AND FORMAT DATA #
########################

# Here we load in the detections data in long format, with one row per replicate.
# We need to reformat this into wide format for input into msocc.
# We also need to reassign Sample and Replicate IDs to be ordered numbers.

# load in the data
alldata<-read.table("data/environmental/detections_qpcr_ddpcr_filtered_200123.txt",sep="\t",header=T)

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

#################
# MODEL FITTING #
#################

# Here I will fit a series of models, beginning with a null model (no covariates)
# Each model will have 11k iterations, as we will 'burn' the first 1k before obtaining parameter estimates for further analyses.
# More information on MCMC sampling can be found here: https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo
# The specific sampling algorithm used by msocc is Gibbs sampling.

# NULL MODEL

# code for null model. Note that covariate tables need to be included even if there are no covariates in the model formula.
mod1 <- msocc_mod(wide_data = detections,
                  site = list(model = ~ 1, cov_tbl = site_covs),
                  sample = list(model = ~ 1, cov_tbl = samp_covs),
                  rep = list(model = ~ 1, cov_tbl = rep_covs),
                  progress = T,num.mcmc=11000,seed = 1598) # random seed is important for reproducibility
# I also obtain a one-row dataframe with the model formula, WAIC score and associated metrics (used to compare all subsequent models)
mod1_score<-get_waic_table(mod1,burnin=1000)

# Trace plots
# If no covariates are specified in the model, the posterior is already converted to the probability scale.
# For the trace plots I want to see the sampling distribution on the original logit scale, so I use the logit() function
par(mfrow=c(2,2))
plot(logit(mod1$psi),type='l',main="beta",ylab="sampled value")
plot(logit(mod1$theta),type='l',main="alpha",ylab="sampled value")
plot(logit(mod1$p),type='l',main="delta",ylab="sampled value")

# they look OK. The big jump at the start of the delta sampling shows why a burn-in is sometimes desired.

# COVARIATE MODELS

# I now fit a series of models with various combinations of three covariates (volume, filter type, and pcr type) at the sample and / or replicate level.
# Inclusion of covariates should be carefully justified.
# Sample volume is included as a  sample level covariate:
# - Sample occupancy / capture probability is likely to be affected because more sample volume gives more opportunity to capture target eDNA in the sample
# Filter type is included only at the replicate level
# - Composition of the filter could affect the efficacy of DNA extraction
# - As the two filter types, Millex and Sterivex have the same pore size, there is no reason to suppose that filter type would affect sample occupancy
# PCR method is also included only at the replicate level

# The following code fits seven further models.

mod2 <- msocc_mod(wide_data = detections,
                  site = list(model = ~ 1, cov_tbl = site_covs),
                  sample = list(model = ~ volume, cov_tbl = samp_covs),
                  rep = list(model = ~ filter, cov_tbl = rep_covs),
                  progress = T,num.mcmc=11000,seed = 1598)
mod2_score<-get_waic_table(mod2,burnin=1000)

mod3 <- msocc_mod(wide_data = detections,
                  site = list(model = ~ 1, cov_tbl = site_covs),
                  sample = list(model = ~ volume, cov_tbl = samp_covs),
                  rep = list(model = ~ 1, cov_tbl = rep_covs),
                  progress = T,num.mcmc=11000,seed = 1598)
mod3_score<-get_waic_table(mod3,burnin=1000)

mod4 <- msocc_mod(wide_data = detections,
                  site = list(model = ~ 1, cov_tbl = site_covs),
                  sample = list(model = ~ 1, cov_tbl = samp_covs),
                  rep = list(model = ~ filter, cov_tbl = rep_covs),
                  progress = T,num.mcmc=11000,seed = 1598)
mod4_score<-get_waic_table(mod4,burnin=1000)

# incorporating PCR method

mod5 <- msocc_mod(wide_data = detections,
                  site = list(model = ~ 1, cov_tbl = site_covs),
                  sample = list(model = ~ volume, cov_tbl = samp_covs),
                  rep = list(model = ~ filter+pcr_method, cov_tbl = rep_covs),
                  progress = T,num.mcmc=11000,seed = 1598)
mod5_score<-get_waic_table(mod5,burnin=1000)

mod6 <- msocc_mod(wide_data = detections,
                  site = list(model = ~ 1, cov_tbl = site_covs),
                  sample = list(model = ~ 1, cov_tbl = samp_covs),
                  rep = list(model = ~ filter+pcr_method, cov_tbl = rep_covs),
                  progress = T,num.mcmc=11000,seed = 1598)
mod6_score<-get_waic_table(mod6,burnin=1000)

mod7 <- msocc_mod(wide_data = detections,
                  site = list(model = ~ 1, cov_tbl = site_covs),
                  sample = list(model = ~ volume, cov_tbl = samp_covs),
                  rep = list(model = ~ pcr_method, cov_tbl = rep_covs),
                  progress = T,num.mcmc=11000,seed = 1598)
mod7_score<-get_waic_table(mod7,burnin=1000)

mod8 <- msocc_mod(wide_data = detections,
                  site = list(model = ~ 1, cov_tbl = site_covs),
                  sample = list(model = ~ 1, cov_tbl = samp_covs),
                  rep = list(model = ~ pcr_method, cov_tbl = rep_covs),
                  progress = T,num.mcmc=11000,seed = 1598)
mod8_score<-get_waic_table(mod8,burnin=1000)

#####################
# MODEL COMPARISONS #
#####################

# Here we use WAIC scores to compare models.
# The WAIC is an 'information criterion' reflecting a compromise between how well the model fits the data and how generalisable it is.
# Too far in either direction leads to underfitting or overfitting
# (underfitting = mdoel does not fit data very well; overfitting: model fits data too well and is less likely to reflect reality).
# Lower WAIC scores reflect a better compromise between model fit and 'generalisability', and therefore models with lower scores are considered 'better'.
# Note that informaiton criteria such as the WAIC are only of use to compare models: a single value is meaningless.
# More inforation can be found here:https://en.wikipedia.org/wiki/Model_selection
# (see specifically the section on 'Criteria')

# We already obtained the WAIC scores after fitting each model.
# Here I compile all model formulae and WAIC scores into a dataframe:
model_scores<-rbind(mod1_score,mod2_score,mod3_score,mod4_score,mod5_score,mod6_score,mod7_score,mod8_score)
model_scores<-model_scores[order(model_scores$waic),] # sort by WAIC score in ascending order (lower = better)
model_scores$model<-row.names(model_scores)

model_scores # view the table

write.table(model_scores,file="Occ_models/qPCR_ddPCR/msocc_model_scores_200123.txt",quote=F,col.names =T,row.names = F,sep="\t")

# According to Spiegelhalter et al (2002), 'rules of thumb' are applicable to differences in AIC and DIC scores:
# https://rss.onlinelibrary.wiley.com/doi/10.1111/1467-9868.00353
# Models with scores of 1-2 unit difference of the 'best' model deserve consideration, while those with 3-7 units of difference have less support.
# We could use this rule in the other direction to say that models with scores at least 3 units lower than the null model deserve consideration.

################################################
# VISUALISING PARAMETER ESTIMATES FROM MODEL 5 #
################################################

# TRACE PLOTS

# We should examine the traceplots of the parameter estimates before going further with estimates of theta and p
par(mfrow=c(2,3))
plot(logit(mod5$psi),type='l',main="beta",ylab="sampled value")
plot(as.data.frame(mod5$alpha)$`(Intercept)`,type='l',main="alpha.intercept",ylab="sampled value")
plot(as.data.frame(mod5$alpha)$volume,type='l',main="alpha.volume",ylab="sampled value")
plot(as.data.frame(mod5$delta)$`(Intercept)`,type='l',main="delta.intercept",ylab="sampled value")
plot(as.data.frame(mod5$delta)$filter,type='l',main="delta.filter(Sterivex)",ylab="sampled value")
plot(as.data.frame(mod5$delta)$pcr_method,type='l',main="delta.pcr_method(ddPCR)",ylab="sampled value")
# they all look hunky dory

# extract posterior samples of delta parameters from the model (last 10k samples = 1k burn-in), and convert to p estimates for millex and sterivex
alphas<-tail(as.data.frame(mod5$alpha),10000)
deltas<-tail(as.data.frame(mod5$delta),10000)

thetas<-data.frame(theta=c(inv.logit(alphas$`(Intercept)`+alphas$volume*300),
                           inv.logit(alphas$`(Intercept)`+alphas$volume*600)),
                   volume=c(rep("300",10000),rep("600",10000)))

ps<-data.frame(p=c(inv.logit(deltas$`(Intercept)`+deltas$pcr_methodqPCR), # Millex_qPCR
                   inv.logit(deltas$`(Intercept)`+deltas$pcr_methodqPCR+deltas$filterSterivex), # Sterivex_qPCR
                   inv.logit(deltas$`(Intercept)`+deltas$filterSterivex)), # Sterivex_ddPCR
               filter=c(rep("Millex",10000),rep("Sterivex",20000)),
               pcr_method=c(rep("qPCR",20000),rep("ddPCR",10000)))

# obtain median and credible intervals for theta
theta_summarised<-thetas %>% group_by(volume) %>%
  summarise(median=median(theta),
            lower=quantile(theta,probs=c(.025,.975))[1],
            upper=quantile(theta,probs=c(.025,.975))[2])

# obtain median and credible intervals for p
p_summarised<-ps %>% group_by(filter,pcr_method) %>%
  summarise(median=median(p),
            lower=quantile(p,probs=c(.025,.975))[1],
            upper=quantile(p,probs=c(.025,.975))[2])

# PLOTTING
# plot of p as a function of filter type + pcr
thetaplot1<-ggplot(theta_summarised,aes(x=volume,y=median))+
  theme_bw()+
  geom_point()+
  geom_errorbar(aes(ymin=lower,ymax=upper),width=0)+
  labs(x="Sample volume (ml)",y="\U03B8 (median +/- 95% credible intervals)",title="A")+
  ylim(0,1)

p_summarised$method_combo<-paste(p_summarised$filter,p_summarised$pcr_method,sep="_")
p_summarised$method_combo<-factor(p_summarised$method_combo,levels=c("Millex_qPCR","Sterivex_qPCR","Sterivex_ddPCR"))
# plot of p as a function of filter type + pcr
pplot1<-ggplot(p_summarised,aes(x=method_combo,y=median,col=method_combo))+
  theme_bw()+
  geom_point()+
  geom_errorbar(aes(ymin=lower,ymax=upper),width=0)+
  labs(x="Filter type / PCR method combination",y="p (median +/- 95% credible intervals)",title="B")+
  ylim(0,1)+
  scale_x_discrete(labels=c("Millex / qPCR", "Sterivex / qPCR", "Sterivex / ddPCR"))+
  scale_color_manual(values=c("black","black","darkred"))+
  theme(legend.position = "none")

#####################################
# EXAMINING SITE-SPECIFIC VARIATION #
#####################################

# An additional model fit to attempt to explain the slightly lower detection probability estimates for ddPCR.
# The following model is the same as model 5 but includes an interaction term of site with pcr_method.
# This results in an additional slope term (pcr effect) being estimated for each site

mod5site <- msocc_mod(wide_data = detections,
                      site = list(model = ~ 1, cov_tbl = site_covs),
                      sample = list(model = ~ volume, cov_tbl = samp_covs),
                      rep = list(model = ~ filter+pcr_method*site, cov_tbl = rep_covs),
                      progress = T,num.mcmc=11000,seed = 1598)
mod5site_score<-get_waic_table(mod5site,burnin=1000)

# Horrible code to extract the posterior estimates of theta and p, accounting for site specific interactions
deltas2<-tail(as.data.frame(mod5site$delta),10000)
deltas2$siteAA<-0 # add the slope parameter for siteAA, which is 0 as it is the 'base' condition, however we add this for simplicity
deltas2$`pcr_methodqPCR:siteAA`<-0 # add the interaction parameter for siteAA, which is 0 as it is the 'base' condition, however we add this for simplicity
siteslopes<-deltas2[c(4:8,14)]
siteinteracts<-deltas2[c(9:13,15)]

# loop through each of the site-specific slopes and interaction parameters (including the zero ones of AA) and calculate p
ps2<-NULL
for (i in c(1:6)){
  site_slopes<-siteslopes[i]
  site_interacts<-siteinteracts[i]
  colnames(site_slopes)[1]<-"site_slope"
  colnames(site_interacts)[1]<-"site_interaction"
  # get p for each filter type / pcr method combo
  ps2_site<-data.frame(p=c(inv.logit(deltas2$`(Intercept)`+deltas2$pcr_methodqPCR+site_slopes$site_slope+site_interacts$site_interaction), # Millex_qPCR
                           inv.logit(deltas2$`(Intercept)`+deltas2$pcr_methodqPCR+deltas2$filterSterivex+site_slopes$site_slope+site_interacts$site_interaction), # Sterivex_qPCR
                           inv.logit(deltas2$`(Intercept)`+deltas2$filterSterivex+site_slopes$site_slope)), # Sterivex_ddPCR
                       filter=c(rep("Millex",10000),rep("Sterivex",20000)),
                       pcr_method=c(rep("qPCR",20000),rep("ddPCR",10000)),
                       sitenum=i)
  ps2<-rbind(ps2,ps2_site)}
# re-add sites names
sitenames<-data.frame(sitenum=c(1:6),site=c("Fu","Ly","U","Wi RON","Wi ROT","AA"))
ps2=merge(ps2,sitenames,by="sitenum");ps2$sitenum<-NULL

# obtain median and credible intervals for p
p_summarised<-ps2 %>% group_by(filter,pcr_method,site) %>%
  summarise(median=median(p),
            lower=quantile(p,probs=c(.025,.975))[1],
            upper=quantile(p,probs=c(.025,.975))[2])

# PLOTTING (supplementary figure)
p_summarised$method_combo<-paste(p_summarised$filter,p_summarised$pcr_method,sep="_")
p_summarised$method_combo<-factor(p_summarised$method_combo,levels=c("Millex_qPCR","Sterivex_qPCR","Sterivex_ddPCR"))
# plot of p as a function of filter type + pcr * site
pplot2<-ggplot(p_summarised,aes(x=method_combo,y=median,col=method_combo))+
  theme_bw()+
  geom_point()+
  geom_errorbar(aes(ymin=lower,ymax=upper),width=0)+
  labs(x="Filter type / PCR method combination",y="p (median +/- 95% credible intervals)")+
  facet_wrap(~site)+
  ylim(0,1)+
  scale_x_discrete(labels=c("Millex / qPCR", "Sterivex / qPCR", "Sterivex / ddPCR"))+
  scale_color_manual(values=c("black","black","darkred"))+
  theme(legend.position = "none")

# print plots
lay <- rbind(c(1,1,2,2,2))
grobz <- lapply(list(thetaplot1,pplot1), ggplotGrob)

# png(filename = "Occ_models/qPCR_ddPCR/theta_p_msocc_200123.png", width = 4000, height = 2500,res=600)
# grid.arrange(grobs = grobz, layout_matrix = lay)
# dev.off()

# supplementary figure
# png(filename = "Occ_models/qPCR_ddPCR/pplot_rivers_170223.png", width = 6000, height = 4000,res=600)
# pplot2
# dev.off()
