####################################################################################
# LIMITS OF DETECTION AND QUANTIFICATION FOR qPCR AND ddPCR | JAMES ORD | 13/01/22 #
####################################################################################

Sys.setenv(LANG = "en")

library(ggplot2)
library(stringr)
library(dplyr)
library(gridExtra)
library(matrixStats)
library("drc")

getwd()
rm(list=ls())

# Here I predict limits of detection for Bettge and Hutchins PKD assays when using qPCR and ddPCR.
# I then calculate coefficients of variation from predicted concentration values to illustrate some limit of quantification.
# In case of LOD, detection is treated as a binary variable (in case of both qPCR and ddPCR).
# Here, I use concentration to refer to copies per ul input (where input was 4ul).
# This is derived by getting the copies per 20ul reaction (copies per ul reaction volume * 20) and dividing it by 4.

###################
# 1) READ IN DATA #
###################

# 1.1 qPCR

# read in the qPCR data
qPCRcor<-read.csv("data/LOD_LOQ/qPCR_for_LOD_LOQ_corrected_assaysrenamed.csv",sep=";")
qPCRcor$pcr<-"qPCR"

# I'll omit the blanks (SQ of 0)
qPCRcor<-subset(qPCRcor,SQ>0)
# qPCRcor$SQ<-qPCRcor$SQ*20 # CONVERT TO TOTAL EXPECTED NUMBER OF COPIES (copies/ul * 20ul reaction)
qPCRcor$SQ<-(qPCRcor$SQ*20)/4 # CONVERT TO EXPECTED CONCENTRATION IN 1ul INITIAL INPUT

# add binary response variable
qPCRcor$Pos_binary0<-ifelse(is.na(qPCRcor$Cq),0,1)
qPCRcor$Pos_binary <-ifelse(qPCRcor$Cq>40&qPCRcor$Pos_binary0==1,0,qPCRcor$Pos_binary0)
qPCRcor$Pos_binary0<-NULL

# MODELING QUANTITY FROM CQ

# because qPCR does not directly output a concentration estimate, we need to calculate this from the CQ values
# Here, I fit linear models for this purpose, modelling CQ as a function of SQ
# I do this separately for Bettge and Hutchins assays

m1Bettge<-lm(Cq~log(SQ),data=subset(qPCRcor,Target=="Bettge"))
m1Hutchins<-lm(Cq~log(SQ),data=subset(qPCRcor,Target=="Hutchins"))

# the models give coefficients from which a CQ value can be predicted for a given SQ.
ccBettge <- coef(m1Bettge)
ccHutchins <- coef(m1Hutchins)

# the same goes in reverse: using the model coefficients (intercept and slope estimates), we can 'back calculate' a concentration from a given CQ.

# proof of concept
# # regular forward calculation from equation
# ccBettge[1]+(ccBettge[2])*log(20000)
#
# # reverse calculation
# exp((25-ccBettge[1])/ccBettge[2])
# exp((25-ccHutchins[1])/ccHutchins[2])

# Here I calculate back from the equation and add these 'modeled' concentrations to the original data frame
qPCRcor$modeled_Q<-ifelse(qPCRcor$Target=="Bettge",
                          exp((qPCRcor$Cq-ccBettge[1])/ccBettge[2]),
                          exp((qPCRcor$Cq-ccHutchins[1])/ccHutchins[2]))

# Theoretically, one could do the above with just one model but I am too dumb to figure out how to do the reverse calculation with the interaction term.

# Some plots to show that the above didn't go horribly wrong

# the following plot shows the Cq values against the initial expected concentrations:
Cq_plot<-ggplot(qPCRcor,aes(x=log(SQ),y=Cq,color=Target,fill=Target))+
  theme_bw()+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE)+
  scale_color_manual(values=c("black","darkgrey"))+
  scale_fill_manual(values=c("black","darkgrey"))+
  labs(x="Log of expected concentration\n(copies / \U00B5l input)",y="Cq",title="A")+
  theme(legend.position = "top")

# the following plot shows the modeled concentrations against the initial expected ones:
modeled_Q_plot<-ggplot(qPCRcor,aes(x=log(SQ),y=log(modeled_Q),color=Target,fill=Target))+
  theme_bw()+
  geom_point()+
  geom_smooth(method = "lm", se = TRUE)+
  scale_color_manual(values=c("black","darkgrey"))+
  scale_fill_manual(values=c("black","darkgrey"))+
  labs(x="Log of expected concentration\n(copies / \U00B5l input)",y="Log of modeled concentration\n(copies / \U00B5l input)",title="B")+
  theme(legend.position = "top")+
  ylim(c(-2.5,12))

# Note that the above plot of modeled concentration does not include replicates that did not amplify, as without a Cq value, the concentration cannot be modeled.

# 1.2 ddPCR

# read in the ddPCR data
ddPCR<-read.csv("data/LOD_LOQ/ddPCR_for_LOD_LOQ_small_assaysrenamed.csv",sep=";")[c(2,3,5,6)]
ddPCR$pcr<-"ddPCR"
ddPCR$SQ<-(ddPCR$SQ*20)/4;ddPCR$modeled_Q<-(ddPCR$modeled_Q*20)/4 # CONVERT TO EXPECTED CONCENTRATION IN 1ul INITIAL INPUT

# As the concentrations already derive from a statistical model, I won't bother with any outlier identification and instead trust the estimates,
# except in the case of the arbitrary values for unquantifiable high concentrations, which we remove later prior to LOQ calculation.

# add binary response variable
ddPCR$Pos_binary<-ifelse(ddPCR$Positives>0,1,0)

head(qPCRcor)
head(ddPCR)

# combining
all_PCR<-rbind(qPCRcor[c(1,3:6)],ddPCR[c(1,3:6)])
all_PCR$assay_method<-paste(all_PCR$pcr,all_PCR$Target,sep="_")

##########
# 2) LOD #
##########

# Method adapted from klymus script - reproduces same LODs for qPCR assays

# define stuff first
LOD.FCTS <- list(LL.2(),LL.3(),LL.3u(),LL.4(),LL.5(),W1.2(),W1.3(),W1.4(),W2.2(),W2.3(),
                 W2.4(),AR.2(),AR.3(),MM.2(),MM.3())
LODs<-NULL
newdata<-NULL
for (i in levels(as.factor(all_PCR$assay_method))){
# Fit initial model
subdata<-subset(all_PCR,assay_method==i)
LOD.mod <- drm(Pos_binary~SQ,data=subdata,fct=W2.4())
## Test all available models and select the best one:
LOD.FCT2 <- row.names(mselect(LOD.mod,LOD.FCTS))[1]
LOD.FCT3 <- getMeanFunctions(fname=LOD.FCT2)[1]
# Fit the second model with the 'best' function
LOD.mod2 <- drm(Pos_binary~SQ,data=subset(all_PCR,assay_method==i),fct=LOD.FCT3[[1]])

# Put LOD and SE of LOD into dataframe
LOD_method<-as.data.frame(ED(LOD.mod2,0.95,type="absolute"))
LOD_method$assay_method<-i
LOD_method$mod_func<-LOD.FCT2
# append to existing df
LODs<-rbind(LODs,LOD_method)

# get predictions for plotting purposes
newdata_method <- data.frame(logSQ=seq(min(log(all_PCR$SQ)),max(log(all_PCR$SQ)),0.05))
newdata_method$SQ<-exp(newdata_method$logSQ)
preds<-as.data.frame(predict(LOD.mod2,newdata_method,se.fit =T))
newdata_method$pred<-preds$Prediction
newdata_method$pred.SE<-preds$SE
newdata_method$assay_method<-i
# append to existing df
newdata<-rbind(newdata,newdata_method)
}

# view the LOD predictions
LODs

write.csv(LODs, file="LOD_LOQ/LODs_130123.csv",row.names =FALSE,quote=FALSE)

LODs_forplot<-LODs
LODs_forplot$assay_method<-str_replace(LODs_forplot$assay_method,"_Bettge","\n(Bettge\net al.)")
LODs_forplot$assay_method<-str_replace(LODs_forplot$assay_method,"_Hutchins","\n(Hutchins\net al.)")
LODs_forplot$assay_method<-factor(LODs_forplot$assay_method,levels=c("qPCR\n(Bettge\net al.)","qPCR\n(Hutchins\net al.)","ddPCR\n(Bettge\net al.)","ddPCR\n(Hutchins\net al.)"))
LOD.plot1<-ggplot(LODs_forplot,aes(x=assay_method,y=Estimate,col=assay_method))+
  theme_bw()+
  geom_point()+
  geom_errorbar(aes(ymin=Estimate-`Std. Error`,ymax=Estimate+`Std. Error`),width=0)+
  labs(y="LOD estimate\n(copies / ul input +/-SE)",x="Assay method",title="")+
  scale_color_manual(values=c("black","darkgrey","darkred","lightpink"))+
  theme(legend.position = "none")+
  theme(panel.grid.major.x = element_blank())+
  theme(panel.grid.minor.x = element_blank())+
  theme(plot.background = element_rect(colour = "black", fill="white", linewidth=0.5))

newdata$assay_method<-str_replace(newdata$assay_method,"_Bettge"," (Bettge et al.)")
newdata$assay_method<-str_replace(newdata$assay_method,"_Hutchins"," (Hutchins et al.)")
newdata$assay_method<-factor(newdata$assay_method,levels=c("qPCR (Bettge et al.)","qPCR (Hutchins et al.)","ddPCR (Bettge et al.)","ddPCR (Hutchins et al.)"))

head(newdata)
# plot code
LOD.plot2<-ggplot(subset(newdata,SQ<=30),aes(x=SQ,y=pred))+
  # theme preset
  theme_bw()+
  # main trendlines
  geom_line(aes(col=assay_method))+
  geom_ribbon(aes(fill=assay_method,ymin=pred-pred.SE,ymax=pred+pred.SE),alpha=0.1)+
  # put a line at the 0.95 mark
  geom_hline(yintercept=0.95,lty="dashed")+
  # set manual color / fill scale
  scale_color_manual(values=c("black","darkgrey","darkred","lightpink"),name="Assay method")+
  scale_fill_manual(values=c("black","darkgrey","darkred","lightpink"),name="Assay method")+
  # set axis labels and plot title
  labs(x="Copies / \U00B5l input",y="Detection probability (+/- SE)",title="A")+
  # final legend aesthetics
  theme(legend.position = c(0.8,0.25))+
  theme(legend.background = element_rect(colour ="black"))+
  annotate("text", x = 27.5, y = 0.9, label = "Pr(positive) = 0.95")+
  # dashed lines to show LOD for each assay combination
  # qPCR-Bettge
  geom_vline(xintercept = LODs[3,1],lty="dashed",col="black")+
  # qPCR-Hutchins
  geom_vline(xintercept = LODs[4,1],lty="dashed",col="darkgrey")+
  # ddPCR-Bettge
  geom_vline(xintercept = LODs[1,1],lty="dashed",col="darkred")+
  # ddPCR-Hutchins
  geom_vline(xintercept = LODs[2,1],lty="dashed",col="lightpink")+
  annotation_custom(ggplotGrob(LOD.plot1), xmin = 7.5, xmax = 20, 
                    ymin = 0, ymax = 0.75)+
  scale_x_continuous(breaks = seq(0,30, by = 5))

################
#     1.2) LOQ #
################

# combining
all_PCR2<-all_PCR

# COEFFICIENT OF VARIATION

# for this, we include negatives as '0' values
all_PCR2$modeled_Q<-ifelse(all_PCR2$Pos_binary==0,0,all_PCR2$modeled_Q)

# Now we get the mean and SD of the modeled concentration value for each expected concentration and assay (Target)
# We consider a lowest possible SQ of 0.5 (= 2 copies per reaction). This I think is the lowest standard that is relevant to the LOQ.

CVs <- subset(all_PCR2,SQ>=0.5&SQ<=30) %>% group_by(SQ,Target,pcr) %>% summarise(mean_mQ=mean(modeled_Q),SD=sd(modeled_Q))
# Because the CV is just the SD divided by the mean
CVs$CV<-CVs$SD/CVs$mean_mQ
CVs<-CVs[complete.cases(CVs),] # remove those for which CV could not be computed (mean and SD of 0)

# POLYNOMIAL MODEL OF CV
# Here I use a model to approximate LOQ values as the value at which the CV is 0.35.
# The polynomial model-based method is the same as that used by Klymus et al.

# Make a load of models for the CV with different polynomial terms
# The klymus script fits more models, I believe, but 5 seems like a reasonable max.
LOQmod1 <- lm(CV ~ poly(SQ, 1)*Target*pcr, data = CVs)
LOQmod2 <- lm(CV ~ poly(SQ, 2)*Target*pcr, data = CVs)
LOQmod3 <- lm(CV ~ poly(SQ, 3)*Target*pcr, data = CVs)
LOQmod4 <- lm(CV ~ poly(SQ, 4)*Target*pcr, data = CVs)
LOQmod5 <- lm(CV ~ poly(SQ, 5)*Target*pcr, data = CVs)
AIC(LOQmod1,LOQmod2,LOQmod3,LOQmod4,LOQmod5) # model 5 seems best here

# now we have fit a model, we can predict from it:
# first establish a sequence of simulated concentration values for each assay
newdata_CV <- data.frame(SQ=rep(seq(min(CVs$SQ),max(CVs$SQ),0.05),4))

newdata_CV$Target=c(rep("Bettge",nrow(newdata_CV)/2),rep("Hutchins",nrow(newdata_CV)/2))
newdata_CV$pcr=c(rep("qPCR",nrow(newdata_CV)/4),rep("ddPCR",nrow(newdata_CV)/4),rep("qPCR",nrow(newdata_CV)/4),rep("ddPCR",nrow(newdata_CV)/4))

# obtain predicted CV values along the simulated concentration values
preds_CV<-predict(LOQmod4, newdata_CV)

# put predicted values in the dataframe with the simulated concentration values
# together, these will allow us to plot smooth curves
newdata_CV$fit<-preds_CV

# loop through all combinations and derive LOQ for each one
LOQ_all<-NULL
for (j in c("qPCR","ddPCR")){
  LOQ_pcr<-NULL
  for (i in c("Bettge","Hutchins")){
    LOQ_target<-data.frame(LOQ=approx(x = subset(newdata_CV,Target==i&pcr==j)$fit, y = subset(newdata_CV,Target==i&pcr==j)$SQ, xout=0.35)$y,
                           Target=i,
                           pcr=j)
    LOQ_pcr<-rbind(LOQ_pcr,LOQ_target)
  }
  LOQ_all<-rbind(LOQ_all,LOQ_pcr)
}

#write.csv(LOQ_all, file="LOD_LOQ/LOQs_130123.csv",row.names =FALSE,quote=FALSE)

newdata_CV$pcr_target<-paste(newdata_CV$pcr," (",newdata_CV$Target," et al.)",sep="")
CVs$pcr_target<-paste(CVs$pcr," (",CVs$Target," et al.)",sep="")
newdata_CV$pcr_target<-factor(newdata_CV$pcr_target,levels=c("qPCR (Bettge et al.)","qPCR (Hutchins et al.)","ddPCR (Bettge et al.)","ddPCR (Hutchins et al.)"))

#head(CVs)

LOQs_forplot<-LOQ_all
LOQs_forplot$assay_method<-paste(LOQs_forplot$pcr,"\n(",LOQs_forplot$Target,"\net al.)",sep="")
LOQs_forplot$assay_method<-factor(LOQs_forplot$assay_method,levels=c("qPCR\n(Bettge\net al.)","qPCR\n(Hutchins\net al.)","ddPCR\n(Bettge\net al.)","ddPCR\n(Hutchins\net al.)"))
LOQ.plot1<-ggplot(LOQs_forplot,aes(x=assay_method,y=LOQ,fill=assay_method))+
  theme_bw()+
  geom_bar(stat="identity",width=0.5)+
  labs(y="LOQ estimate\n(copies / ul input)",x="Assay method",title="")+
  scale_fill_manual(values=c("black","darkgrey","darkred","lightpink"))+
  theme(legend.position = "none")+
  theme(panel.grid.major.x = element_blank())+
  theme(panel.grid.minor.x = element_blank())+
  theme(plot.background = element_rect(colour = "black", fill="white", linewidth=0.5))

# plot code
LOQ.plot2<-ggplot(newdata_CV,aes(x=SQ,y=fit))+
  theme_bw()+
  geom_line(aes(col=pcr_target))+
  geom_point(data=CVs,aes(x=SQ,y=CV,color=pcr_target))+
  geom_hline(yintercept=0.35,lty="dashed")+
  # dashed lines to show LOQ for each assay combination
  # qPCR-Bettge
  geom_vline(xintercept = LOQ_all[1,1],lty="dashed",col="black")+
  # ddPCR-Bettge
  geom_vline(xintercept = LOQ_all[2,1],lty="dashed",col="darkgrey")+
  # qPCR-Hutchins
  geom_vline(xintercept = LOQ_all[3,1],lty="dashed",col="darkred")+
  # ddPCR-Hutchins
  geom_vline(xintercept = LOQ_all[4,1],lty="dashed",col="lightpink")+
  # set manual color / fill scale
  scale_color_manual(values=c("black","darkgrey","darkred","lightpink"),name="Assay method")+
  scale_fill_manual(values=c("black","darkgrey","darkred","lightpink"),name="Assay method")+
  # set axis labels and plot title
  labs(x="Copies / \U00B5l input",y="Coefficient of variation",title="B")+
  # final legend aesthetics
  theme(legend.position = c(0.8,0.75))+
  theme(legend.background = element_rect(colour ="black"))+
  annotate("text", x = 27.5, y = 0.5, label = "CV = 0.35")+
  annotation_custom(ggplotGrob(LOQ.plot1), xmin = 7.5, xmax = 20, 
                    ymin = 0.6, ymax = 2.1)+
  scale_y_continuous(labels = function(x) format(x, nsmall = 2))+
  scale_x_continuous(breaks = seq(0,30, by = 5))

LOQ.plot2

#################
# 3) FINAL PLOT #
#################

lay <- rbind(c(1,1),
             c(2,2))
grobz <- lapply(list(LOD.plot2,LOQ.plot2), ggplotGrob)

# png(filename = "LOD_LOQ/LOD_LOQ_qPCR_ddPCR_130123.png", width = 5000, height = 5000,res=600)
# grid.arrange(grobs = grobz, layout_matrix = lay)
# dev.off()
# 
# png(filename = "LOD_LOQ/qPCR_modeled_concs_130123.png", width = 5000, height = 2500,res=600)
# grid.arrange(Cq_plot, modeled_Q_plot,ncol=2)
# dev.off()
