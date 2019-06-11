#Libraries--------------------------------------
library(unmarked)
library(AICcmodavg)
library(ggplot2)
library(ggpubr)
library(plotrix)
library(corrplot)

data=read.csv("DataS1_unmarkedAnalyses.csv")

#Calculated correlation among explanatory habitat variables
sl16=data$SL16
em16=data$EM16
fl16=data$FL16
aq16=data$AQ16
cattle16=data$Cattle16
fish16=data$Fish16
sl16=data$SL16
LogArea=data$LogArea

sl17=data$SL17
em17=data$EM17
fl17=data$FL17
aq17=data$AQ17
cattle17=data$Cattle17
fish17=data$Fish17
sl17=data$SL17
ph=data$pH

correlations2016=cor(cbind(sl16,em16,fl16,aq16,cattle16,fish16,LogArea))
correlations2017=cor(cbind(sl17,em17,fl17,aq17,cattle17,fish17,LogArea,ph))
corrplot(correlations2016, type="upper", order="hclust")
corrplot(correlations2017, type="upper", order="hclust")


options(scipen=999)  # turn off scientific notation like 1e+06


#All four frog species are analyzed separately using a single species, single-season occupancy model for each year of the study (2016 + 2017).
#Species names are abbreviated as follows:
  #HYLA: Hyla chrysoscelis/versicolor (Cope's/Eastern gray treefrog)
  #ACBL: Acris blanchardi (Blanchard's cricket frog)
  #LIBL: Lithobates blairi (Plains leopard frog)
  #LICA: Lithobates catesbeianus (American bullfrog)
  #In all cases, the year being analyzed is denoted as either "16" or "17", and the number of sites included is denoted with "51"

#HYLA16_51_ssom--------------------------------------------------------
hyla16y<-data[,6:9]
n<-nrow(data)
#site level individual covariates
hyla16.site<-data[,55:64]
#put everything together in unmarked data frame

hyla16 <- unmarkedFrameOccu(y = hyla16y, siteCovs = data[55:64],obsCovs=list(Julian=data[35:38],Temp=data[39:42],Time=data[51:54]))
summary(hyla16) 
#Constant Detection and occupancy Model -->  follows format: model<-occu(~detection_formula ~occupancy_formula, dataframe)
#_______________________
#DETECTION
hylad61=occu(~1~1,hyla16)
hylad62=occu(~Time~1,hyla16)
hylad63=occu(~Temp~1,hyla16)
hylad64=occu(~Julian~1,hyla16)
hylad65=occu(~Julian+Time~1,hyla16)
hylad66=occu(~Temp*Time~1,hyla16)
hylad67=occu(~Temp+Time~1,hyla16)
hylad68=occu(~Temp+Julian~1,hyla16)
hylad69=occu(~Temp*Julian~1,hyla16)
#Goodness of Fit Test Bootstrap
#(gof_boot <- mb.gof.test(hylad61, nsim = 500, plot.hist=FALSE))
#AICc Model Selection Table
hyla16.Det.Cand.models=list(hylad61,hylad62,hylad63,hylad64,hylad66,hylad67,hylad68,hylad69)
hyla16.Det.Modnames=c("ConstantDet","TimeDet","TempDet","JulianDet","Temp.by.Time","Temp.Time","Temp.Julian","Temp.by.Julian")
print(aictab(cand.set = hyla16.Det.Cand.models,modnames = hyla16.Det.Modnames,second.ord = TRUE), digits = 4)
#Generate Detection-correct Estimates
backTransform(hylad61,'det') # for detection
backTransform(hylad63,'state') # for occupancy

#___________________________________
#Stage 1 - Univariate models
#NULL
hyla16OccuNull=occu(~Temp~1,hyla16)
#Habitat Covariate Models
hyla16OccuSL=occu(~Temp~SL16,hyla16)
hyla16OccuLogArea=occu(~Temp~LogArea,hyla16)
hyla16OccuCattle=occu(~Temp~Cattle16,hyla16)
hyla16OccuEM=occu(~Temp~EM16,hyla16)
hyla16OccuAQ=occu(~Temp~AQ16,hyla16)
hyla16OccuFL=occu(~Temp~FL16,hyla16)
hyla16OccuFish=occu(~Temp~Fish16,hyla16)
#Full Candidate Set
hyla16.Occu.Cand.models=list(hyla16OccuSL,hyla16OccuLogArea,hyla16OccuCattle,hyla16OccuEM,hyla16OccuAQ,hyla16OccuFL,hyla16OccuFish,hyla16OccuNull)
hyla16.Occu.Modnames=c("SL","LogArea","Cattle","EM","AQ","FL","Fish","Null")
print(aictab(cand.set = hyla16.Occu.Cand.models,modnames = hyla16.Occu.Modnames,second.ord = TRUE), digits = 4)

#        K     AICc Delta_AICc AICcWt Cum.Wt       LL
#EM      4 120.5026     0.0000 0.9072 0.9072 -55.8165
#LogArea 4 126.1859     5.6832 0.0529 0.9601 -58.6582
#Fish    4 128.3372     7.8346 0.0180 0.9781 -59.7338
#FL      4 130.0137     9.5111 0.0078 0.9860 -60.5721
#AQ      4 130.5104    10.0078 0.0061 0.9920 -60.8204
#Null    3 131.2899    10.7873 0.0041 0.9962 -62.3896
#SL      4 132.2410    11.7384 0.0026 0.9987 -61.6857
#Cattle  4 133.6435    13.1409 0.0013 1.0000 -62.3870

#_______________________________
#Stage 2 - Additive Model: for cases where >1 model comprises the 95% confidence set, a single additive model was generated comprising all those variables. 
hyla16OccuEM.Area=occu(~Temp~EM16+LogArea,hyla16)
#FINAL CANDIDATE SET
hyla16.FINAL.Occu.Cand.models=list(hyla16OccuEM.Area,hyla16OccuLogArea,hyla16OccuEM,hyla16OccuNull)
hyla16.FINAL.Occu.Modnames=c("EM.Area","LogArea","EM","Null")
#FINAL AICc Table
print(aictab(cand.set = hyla16.FINAL.Occu.Cand.models,modnames = hyla16.FINAL.Occu.Modnames,second.ord = TRUE), digits = 4)

#..Final Model Set
#        K     AICc Delta_AICc AICcWt Cum.Wt       LL
#EM.Area 5 117.1245     0.0000 0.8359 0.8359 -52.8956
#EM      4 120.5026     3.3782 0.1544 0.9903 -55.8165
#LogArea 4 126.1859     9.0614 0.0090 0.9993 -58.6582
#Null    3 131.2899    14.1654 0.0007 1.0000 -62.3896


#MODEL AVERAGING _____________________________
#Model averaging for each variable   
modavg(cand.set = hyla16.FINAL.Occu.Cand.models, modnames = hyla16.FINAL.Occu.Modnames, parm = "EM16",parm.type = "psi",conf.level = 0.85)
modavg(cand.set = hyla16.FINAL.Occu.Cand.models, modnames = hyla16.FINAL.Occu.Modnames, parm = "LogArea",parm.type = "psi",conf.level = 0.85)

#Create New data frames to store model averaged Predictions in 
dat.predEM <- data.frame(EM16 = seq(from = min(data$EM16),to = max(data$EM16), by = 5),LogArea = mean(data$LogArea))
dat.predArea <- data.frame(LogArea = seq(from = min(data$LogArea),to = max(data$LogArea), by = 0.1),EM16 = mean(data$EM16))
#Model-averaged predictions of psi across range of values
##of sitevar1 and entire model set
predictionsHyla16EM=modavgPred(cand.set = hyla16.FINAL.Occu.Cand.models, modnames = hyla16.FINAL.Occu.Modnames,newdata = dat.predEM, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionsHyla16EM,'predictionsHyla16EM.csv')
predictionsHyla16LogArea=modavgPred(cand.set = hyla16.FINAL.Occu.Cand.models, modnames = hyla16.FINAL.Occu.Modnames,newdata = dat.predArea, parm.type = "psi",conf.level = 0.85)

#write.csv(predictionsHyla16LogArea,'predictionsHyla16LogArea.csv')
#PLOTTING IN GGPLOT_____________________________
#Convert ModAvgPred Output to class Data.Frame to be used in ggplot
    Hyla16EMpreds=as.data.frame(predictionsHyla16EM)
    Hyla16LogAreapreds=as.data.frame(predictionsHyla16LogArea)

    #Write CSVs
#    write.csv(Hyla16EMpreds, file ="Hyla16EMpreds.csv")
#    write.csv(Hyla16LogAreapreds, file ="Hyla16LogAreapreds.csv")
#...HYLA 16 plots -------------------------------
(hyla16em <- ggplot(Hyla16EMpreds, aes(x=seq(from = 0, to = 65, by = 5), y=mod.avg.pred))+
   geom_line(data=Hyla16EMpreds)+
   geom_ribbon(data=Hyla16EMpreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
   theme_pubr()+
   scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
   scale_x_continuous(name="Emergent Vegetation Cover (%)")+
      ggtitle("Gray treefrog spp. (2016)"))
    
(hyla16logarea <- ggplot(Hyla16LogAreapreds, aes(x=seq(from = 2.3, to = 4.3, by = 0.1), y=mod.avg.pred))+
   geom_line(data=Hyla16LogAreapreds)+
   geom_ribbon(data=Hyla16LogAreapreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
   theme_pubr()+
   scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
   scale_x_continuous(name="Log Area")+
      ggtitle("Gray treefrog spp. (2016)"))
  
    
    
# ACBL16_51_ssom------------------------------------------------------
acbl16y<-data[,2:5]
n<-nrow(data)
#site level individual covariates
acbl16.site<-data[,55:64]
#put everything together in unmarked data frame

acbl16 <- unmarkedFrameOccu(y = acbl16y, siteCovs = data[55:64],obsCovs=list(Julian=data[35:38],Temp=data[39:42],Time=data[51:54]))
summary(acbl16) 
#Constant Detection and occupancy Model -->  follows format: model<-occu(~detection_formula ~occupancy_formula, dataframe)
#_______________________
#DETECTION
acbld61=occu(~1~1,acbl16)
acbld62=occu(~Time~1,acbl16)
acbld63=occu(~Temp~1,acbl16)
acbld64=occu(~Julian~1,acbl16)
acbld65=occu(~Julian+Time~1,acbl16)
acbld66=occu(~Temp*Time~1,acbl16)
acbld67=occu(~Temp+Time~1,acbl16)
acbld68=occu(~Temp+Julian~1,acbl16)
acbld69=occu(~Temp*Julian~1,acbl16)
#Goodness of Fit Test Bootstrap
#(gof_boot <- mb.gof.test(acbl16.global, nsim = 500, plot.hist=TRUE))
acbl16.global=occu(~Julian+Time+Temp~EM16+AQ16+FL16+LogArea+Fish16+SL16,acbl16)
#AICc Model Selection Table
acbl16.Det.Cand.models=list(acbld61,acbld62,acbld63,acbld64,acbld66,acbld67,acbld68,acbld69)
acbl16.Det.Modnames=c("ConstantDet","TimeDet","TempDet","JulianDet","Temp.by.Time","Temp.Time","Temp.Julian","Temp.by.Julian")
print(aictab(cand.set = acbl16.Det.Cand.models,modnames = acbl16.Det.Modnames,second.ord = TRUE), digits = 4)
#Generate Detection-corrected Estimates
backTransform(acbld61,'det') # for detection
backTransform(acbld68,'state') # for occupancy
#___________________________________
#Stage 1 - Abiotic and Biotic
#NULL
acbl16OccuNull=occu(~Temp+Julian~1,acbl16)
#Abiotic
acbl16OccuSL=occu(~Temp+Julian~SL16,acbl16)
acbl16OccuLogArea=occu(~Temp+Julian~LogArea,acbl16)
acbl16OccuCattle=occu(~Temp+Julian~Cattle16,acbl16)
#Biotic
acbl16OccuEM=occu(~Temp+Julian~EM16,acbl16)
acbl16OccuAQ=occu(~Temp+Julian~AQ16,acbl16)
acbl16OccuFL=occu(~Temp+Julian~FL16,acbl16)
acbl16OccuFish=occu(~Temp+Julian~Fish16,acbl16)
acbl16.Occu.Cand.models=list(acbl16OccuSL,acbl16OccuLogArea,acbl16OccuCattle,acbl16OccuEM,acbl16OccuAQ,acbl16OccuFL,acbl16OccuFish,acbl16OccuNull)
acbl16.Occu.Modnames=c("SL","LogArea","Cattle","EM","AQ","FL","Fish","Null")
print(aictab(cand.set = acbl16.Occu.Cand.models,modnames = acbl16.Occu.Modnames,second.ord = TRUE), digits = 4)

#        K     AICc Delta_AICc AICcWt Cum.Wt       LL
#AQ      5 126.7395     0.0000 0.6108 0.6108 -57.7031
#EM      5 128.9939     2.2544 0.1979 0.8086 -58.8303
#Null    4 131.4609     4.7214 0.0576 0.8663 -61.2957
#FL      5 132.4239     5.6844 0.0356 0.9019 -60.5453
#LogArea 5 132.6716     5.9321 0.0315 0.9333 -60.6691
#Cattle  5 132.8611     6.1216 0.0286 0.9620 -60.7639
#SL      5 133.4770     6.7375 0.0210 0.9830 -61.0718
#Fish    5 133.9017     7.1622 0.0170 1.0000 -61.2842

#_______________________________
#Stage 2 - Interactions
acbl16OccuEM.AQ=occu(~Temp+Julian~EM16+AQ16,acbl16)
#FINAL CANDIDATE SET
acbl16.FINAL.Occu.Cand.models=list(acbl16OccuEM.AQ,acbl16OccuEM,acbl16OccuAQ,acbl16OccuNull)
acbl16.FINAL.Occu.Modnames=c("EM.AQ","EM","AQ","Null")
#FINAL AICc Table
print(aictab(cand.set = acbl16.FINAL.Occu.Cand.models,modnames = acbl16.FINAL.Occu.Modnames,second.ord = TRUE), digits = 4)

#..Final Model Set
#      K     AICc Delta_AICc AICcWt Cum.Wt       LL
#EM.AQ 6 119.8780     0.0000 0.9561 0.9561 -52.9845
#AQ    5 126.7395     6.8615 0.0309 0.9871 -57.7031
#EM    5 128.9939     9.1159 0.0100 0.9971 -58.8303
#Null  4 131.4609    11.5829 0.0029 1.0000 -61.2957

#MODEL AVERAGING ______________________
#Model averaging for each variable  
modavg(cand.set = acbl16.FINAL.Occu.Cand.models, modnames = acbl16.FINAL.Occu.Modnames, parm = "EM16",parm.type = "psi",conf.level = 0.85)
modavg(cand.set = acbl16.FINAL.Occu.Cand.models, modnames = acbl16.FINAL.Occu.Modnames, parm = "AQ16",parm.type = "psi",conf.level = 0.85)
#Create New data frames to store model averaged Predictions in 
dat.predEM <- data.frame(EM16 = seq(from = min(data$EM16),to = max(data$EM16), by = 5),AQ16 = mean(data$AQ16))
dat.predAQ <- data.frame(AQ16 = seq(from = min(data$AQ16),to = max(data$AQ16), by = 5),EM16 = mean(data$EM16))
#Model-averaged predictions of psi across range of values
##of sitevar1 and entire model set
predictionsacbl16EM=modavgPred(cand.set = acbl16.FINAL.Occu.Cand.models, modnames = acbl16.FINAL.Occu.Modnames,newdata = dat.predEM, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionsacbl16EM,'predictionsacbl16EM.csv')
predictionsacbl16AQ=modavgPred(cand.set = acbl16.FINAL.Occu.Cand.models, modnames = acbl16.FINAL.Occu.Modnames,newdata = dat.predAQ, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionsacbl16LogArea,'predictionsacbl16LogArea.csv')
summary(data$AQ16)
#PLOTTING IN GGPLOT________________________________
#Convert ModAvgPred Output to class Data.Frame to be used in GGPLOT
acbl16EMpreds=as.data.frame(predictionsacbl16EM)
acbl16AQpreds=as.data.frame(predictionsacbl16AQ)
#write to csv
#write.csv(acbl16EMpreds, file ="acbl16EMpreds.csv")
#write.csv(acbl16AQpreds, file ="acbl16AQpreds.csv")

#...ACBL 16 plots -------------------------------
(acbl16em <- ggplot(acbl16EMpreds, aes(x=seq(from = 0, to = 65, by = 5), y=mod.avg.pred))+
    geom_line(data=acbl16EMpreds)+
    geom_ribbon(data=acbl16EMpreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
    theme_pubr()+
    scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
    scale_x_continuous(name="Emergent Vegetation Cover (%)")+
  ggtitle("Blanchard's cricket frog (2016)"))
(acbl16aq <- ggplot(acbl16AQpreds, aes(x=seq(from = 0, to = 100, by = 5), y=mod.avg.pred))+
    geom_line(data=acbl16AQpreds)+
    geom_ribbon(data=acbl16AQpreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
    theme_pubr()+
    scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
    scale_x_continuous(name="Submerged Vegetation Cover (%)")+
  ggtitle("Blanchard's cricket frog (2016)"))

#LIBL16_51_ssom ------------------------------------
libl16y<-data[,10:13]
n<-nrow(data)
#site level individual covariates
libl16.site<-data[,55:64]
#put everything together in unmarked data frame

libl16 <- unmarkedFrameOccu(y = libl16y, siteCovs = data[55:64],obsCovs=list(Julian=data[35:38],Temp=data[39:42],Time=data[51:54]))
summary(libl16) 
#Constant Detection and occupancy Model -->  follows format: model<-occu(~detection_formula ~occupancy_formula, dataframe)
#_______________________
#DETECTION
libld61=occu(~1~1,libl16)
libld62=occu(~Time~1,libl16)
libld63=occu(~Temp~1,libl16)
libld64=occu(~Julian~1,libl16)
libld65=occu(~Julian+Time~1,libl16)
libld66=occu(~Temp*Time~1,libl16)
libld67=occu(~Temp+Time~1,libl16)
libld68=occu(~Temp+Julian~1,libl16)
libld69=occu(~Temp*Julian~1,libl16)
#Goodness of Fit Test Bootstrap

#(gof_boot <- mb.gof.test(libld61, nsim = 500, plot.hist=FALSE))
#AICc Model Selection Table
libl16.Det.Cand.models=list(libld61,libld62,libld63,libld64,libld66,libld67,libld68,libld69)
libl16.Det.Modnames=c("ConstantDet","TimeDet","TempDet","JulianDet","Temp.by.Time","Temp.Time","Temp.Julian","Temp.by.Julian")
print(aictab(cand.set = libl16.Det.Cand.models,modnames = libl16.Det.Modnames,second.ord = TRUE), digits = 4)
#Generate Detection-correct Estimates
backTransform(libld61,'det') # for detection
backTransform(libld61,'state') # for occupancy
#___________________________________
#Stage 1 - Habitat covariates associated with occupancy
#NULL
libl16OccuNull=occu(~1~1,libl16)
#Abiotic
libl16OccuSL=occu(~1~SL16,libl16)
libl16OccuLogArea=occu(~1~LogArea,libl16)
libl16OccuCattle=occu(~1~Cattle16,libl16)
#Biotic
libl16OccuEM=occu(~1~EM16,libl16)
libl16OccuAQ=occu(~1~AQ16,libl16)
libl16OccuFL=occu(~1~FL16,libl16)
libl16OccuFish=occu(~1~Fish16,libl16)
libl16.Occu.Cand.models=list(libl16OccuSL,libl16OccuLogArea,libl16OccuCattle,libl16OccuEM,libl16OccuAQ,libl16OccuFL,libl16OccuFish,libl16OccuNull)
libl16.Occu.Modnames=c("SL","LogArea","Cattle","EM","AQ","FL","Fish","Null")
print(aictab(cand.set = libl16.Occu.Cand.models,modnames = libl16.Occu.Modnames,second.ord = TRUE), digits = 4)

# K     AICc Delta_AICc AICcWt Cum.Wt       LL
#Cattle  3 106.4201     0.0000 0.6940 0.6940 -49.9547
#SL      3 110.5189     4.0989 0.0894 0.7834 -52.0041
#EM      3 111.4816     5.0615 0.0552 0.8386 -52.4855
#Fish    3 111.4966     5.0765 0.0548 0.8935 -52.4930
#AQ      3 111.8005     5.3805 0.0471 0.9406 -52.6450
#Null    2 112.4227     6.0027 0.0345 0.9751 -54.0864
#LogArea 3 114.2658     7.8458 0.0137 0.9888 -53.8776
#FL      3 114.6755     8.2555 0.0112 1.0000 -54.0825

#_______________________________
#Stage 2 - Interactions
libl16OccuSL.Cattle.EM.Fish.AQ=occu(~1~SL16+Cattle16+EM16+Fish16+AQ16,libl16)

libl16OccuSL.Cattle.EM.Fish   =occu(~1~SL16+Cattle16+EM16+Fish16,libl16)
libl16OccuSL.Cattle.EM.AQ     =occu(~1~SL16+Cattle16+EM16+AQ16,libl16)
libl16OccuSL.Cattle.Fish.AQ   =occu(~1~SL16+Cattle16+Fish16+AQ16,libl16)
libl16OccuSL.EM.Fish.AQ       =occu(~1~SL16+EM16+Fish16+AQ16,libl16)
libl16OccuCattle.EM.Fish.AQ   =occu(~1~Cattle16+EM16+Fish16+AQ16,libl16)

libl16OccuSL.Cattle.EM        =occu(~1~SL16+Cattle16+EM16,libl16)
libl16OccuSL.Cattle.Fish      =occu(~1~SL16+Cattle16+Fish16,libl16)
libl16OccuSL.Cattle.AQ        =occu(~1~SL16+Cattle16+AQ16,libl16)
libl16OccuSL.EM.AQ            =occu(~1~SL16+EM16+AQ16,libl16)
libl16OccuSL.EM.Fish          =occu(~1~SL16+EM16+Fish16,libl16)
libl16OccuSL.Fish.AQ          =occu(~1~SL16+Fish16+AQ16,libl16)
libl16OccuCattle.EM.Fish      =occu(~1~Cattle16+EM16+Fish16,libl16)
libl16OccuCattle.Fish.AQ      =occu(~1~Cattle16+Fish16+AQ16,libl16)
libl16OccuCattle.EM.AQ        =occu(~1~Cattle16+EM16+AQ16,libl16)
libl16OccuEM.Fish.AQ          =occu(~1~EM16+Fish16+AQ16,libl16)

libl16OccuSL.Cattle           =occu(~1~SL16+Cattle16,libl16)
libl16OccuSL.EM               =occu(~1~SL16+EM16,libl16)
libl16OccuSL.Fish             =occu(~1~SL16+Fish16,libl16)
libl16OccuSL.AQ               =occu(~1~SL16+AQ16,libl16)
libl16OccuCattle.EM           =occu(~1~Cattle16+EM16,libl16)
libl16OccuCattle.Fish         =occu(~1~Cattle16+Fish16,libl16)
libl16OccuCattle.AQ           =occu(~1~Cattle16+AQ16,libl16)
libl16OccuEM.Fish             =occu(~1~EM16+Fish16,libl16)
libl16OccuEM.AQ               =occu(~1~EM16+AQ16,libl16)
libl16OccuFish.AQ             =occu(~1~Fish16+AQ16,libl16)

libl16OccuSL                  =occu(~1~SL16,libl16)
libl16OccuCattle              =occu(~1~Cattle16,libl16)
libl16OccuEM                  =occu(~1~EM16,libl16)
libl16OccuFish                =occu(~1~Fish16,libl16)
libl16OccuAQ                  =occu(~1~AQ16,libl16)



#                     K     AICc Delta_AICc AICcWt Cum.Wt       LL
#Cattle.EM.Fish       5 104.2627     0.0000 0.1647 0.1647 -46.4647
#Cattle.EM            4 105.5308     1.2681 0.0874 0.2521 -48.3306
#Cattle.EM.Fish.AQ    6 105.9024     1.6397 0.0725 0.3246 -45.9966
#SL.Cattle.EM         5 105.9438     1.6812 0.0711 0.3957 -47.3052
#SL.Cattle.EM.Fish    6 106.0798     1.8171 0.0664 0.4621 -46.0853
#Cattle.Fish          4 106.2076     1.9450 0.0623 0.5243 -48.6690
#Cattle.EM.AQ         5 106.3886     2.1259 0.0569 0.5812 -47.5276
#Cattle               3 106.4201     2.1574 0.0560 0.6372 -49.9547
#SL.Cattle.EM.AQ      6 107.1842     2.9215 0.0382 0.6755 -46.6375
#SL.EM                4 107.4661     3.2034 0.0332 0.7087 -49.2982
#SL.EM.AQ             5 107.5024     3.2397 0.0326 0.7413 -48.0845
#SL.Cattle            4 107.7709     3.5082 0.0285 0.7698 -49.4507
#SL.Cattle.EM.Fish.AQ 7 107.8186     3.5559 0.0278 0.7976 -45.6070
#Cattle.AQ            4 107.8275     3.5648 0.0277 0.8253 -49.4789
#Cattle.Fish.AQ       5 108.0770     3.8143 0.0245 0.8498 -48.3718
#SL.EM.Fish           5 108.1521     3.8894 0.0236 0.8733 -48.4094
#SL.Cattle.Fish       5 108.3443     4.0816 0.0214 0.8947 -48.5055
#SL.EM.Fish.AQ        6 108.8037     4.5410 0.0170 0.9117 -47.4473
#SL.Cattle.AQ         5 109.5226     5.2600 0.0119 0.9236 -49.0947
#EM.AQ                4 109.8565     5.5938 0.0100 0.9336 -50.4935
#EM.Fish.AQ           5 109.9179     5.6552 0.0097 0.9434 -49.2923
#EM.Fish              4 110.0055     5.7429 0.0093 0.9527 -50.5680
#SL.Cattle.Fish.AQ    6 110.3820     6.1193 0.0077 0.9604 -48.2365
#SL                   3 110.5189     6.2563 0.0072 0.9676 -52.0041
#SL.Fish              4 110.9923     6.7296 0.0057 0.9733 -51.0614
#SL.AQ                4 111.2731     7.0104 0.0049 0.9783 -51.2018
#EM                   3 111.4816     7.2189 0.0045 0.9827 -52.4855
#Fish                 3 111.4966     7.2339 0.0044 0.9872 -52.4930
#AQ                   3 111.8005     7.5379 0.0038 0.9910 -52.6450
#Fish.AQ              4 112.0957     7.8330 0.0033 0.9942 -51.6131
#SL.Fish.AQ           5 112.2957     8.0330 0.0030 0.9972 -50.4812
#Null                 2 112.4227     8.1600 0.0028 1.0000 -54.0864

#FINAL CANDIDATE SET
libl16.FINAL.Occu.Cand.models=list(libl16OccuSL.Cattle.EM.Fish.AQ,libl16OccuSL.Cattle.EM.Fish,libl16OccuSL.Cattle.EM.AQ,libl16OccuSL.Cattle.Fish.AQ,libl16OccuSL.EM.Fish.AQ,libl16OccuCattle.EM.Fish.AQ,libl16OccuSL.Cattle.EM,libl16OccuSL.Cattle.Fish,libl16OccuSL.Cattle.AQ,libl16OccuSL.EM.AQ,libl16OccuSL.EM.Fish,libl16OccuSL.Fish.AQ,libl16OccuCattle.EM.Fish,libl16OccuCattle.Fish.AQ,libl16OccuCattle.EM.AQ,libl16OccuEM.Fish.AQ,libl16OccuSL.Cattle,libl16OccuSL.EM,libl16OccuSL.Fish,libl16OccuSL.AQ,libl16OccuCattle.EM,libl16OccuCattle.Fish,libl16OccuCattle.AQ,libl16OccuEM.Fish,libl16OccuEM.AQ,libl16OccuFish.AQ,libl16OccuSL,libl16OccuCattle,libl16OccuEM,libl16OccuFish,libl16OccuAQ,libl16OccuNull)
libl16.FINAL.Occu.Modnames=c("SL.Cattle.EM.Fish.AQ","SL.Cattle.EM.Fish","SL.Cattle.EM.AQ","SL.Cattle.Fish.AQ","SL.EM.Fish.AQ","Cattle.EM.Fish.AQ","SL.Cattle.EM","SL.Cattle.Fish","SL.Cattle.AQ","SL.EM.AQ","SL.EM.Fish","SL.Fish.AQ","Cattle.EM.Fish","Cattle.Fish.AQ","Cattle.EM.AQ","EM.Fish.AQ","SL.Cattle","SL.EM","SL.Fish","SL.AQ","Cattle.EM","Cattle.Fish","Cattle.AQ","EM.Fish","EM.AQ","Fish.AQ","SL","Cattle","EM","Fish","AQ","Null")
#FINAL AICc Table
print(aictab(cand.set = libl16.FINAL.Occu.Cand.models,modnames=libl16.FINAL.Occu.Modnames,second.ord = TRUE), digits = 4)

#MODEL AVERAGING_____________________________
#Model averaging for each variable   
modavg(cand.set = libl16.FINAL.Occu.Cand.models, modnames = libl16.FINAL.Occu.Modnames, parm = "SL16",parm.type = "psi",conf.level = 0.85)
modavg(cand.set = libl16.FINAL.Occu.Cand.models, modnames = libl16.FINAL.Occu.Modnames, parm = "Cattle16",parm.type = "psi",conf.level = 0.85)
modavg(cand.set = libl16.FINAL.Occu.Cand.models, modnames = libl16.FINAL.Occu.Modnames, parm = "Fish16",parm.type = "psi",conf.level = 0.85)
modavg(cand.set = libl16.FINAL.Occu.Cand.models, modnames = libl16.FINAL.Occu.Modnames, parm = "AQ16",parm.type = "psi",conf.level = 0.85)
modavg(cand.set = libl16.FINAL.Occu.Cand.models, modnames = libl16.FINAL.Occu.Modnames, parm = "EM16",parm.type = "psi",conf.level = 0.85)

#Create New data frames to store model averaged Predictions in 
dat.predSL <- data.frame(SL16 = seq(from = min(data$SL16),to = max(data$SL16), by = 0.01),Cattle16=mean(data$Cattle16),EM16=mean(data$EM16),Fish16=mean(data$Fish16),AQ16=mean(data$AQ16))
dat.predCattle <- data.frame(Cattle16 = factor(c("0", "1")),SL16=mean(data$SL16),EM16=mean(data$EM16),Fish16=mean(data$Fish16),AQ16=mean(data$AQ16))
dat.predFish <- data.frame(Fish16 = factor(c("0", "1")),SL16=mean(data$SL16),EM16=mean(data$EM16),Cattle16=mean(data$Cattle16),AQ16=mean(data$AQ16))
dat.predAQ <- data.frame(AQ16 = seq(from = min(data$AQ16),to = max(data$AQ16), by = 5),Cattle16=mean(data$Cattle16),SL16=mean(data$SL16),Fish16=mean(data$Fish16),EM16=mean(data$EM16))
dat.predEM <- data.frame(EM16 = seq(from = min(data$EM16),to = max(data$EM16), by = 5),Cattle16=mean(data$Cattle16),SL16=mean(data$SL16),Fish16=mean(data$Fish16),AQ16=mean(data$AQ16))

#Model-averaged predictions of psi across range of values
predictionslibl16SL=modavgPred(cand.set = libl16.FINAL.Occu.Cand.models, modnames = libl16.FINAL.Occu.Modnames,newdata = dat.predSL, parm.type = "psi",conf.level = 0.85)
predictionslibl16Cattle=modavgPred(cand.set = libl16.FINAL.Occu.Cand.models, modnames = libl16.FINAL.Occu.Modnames,newdata = dat.predCattle, parm.type = "psi",conf.level = 0.85)
predictionslibl16Fish=modavgPred(cand.set = libl16.FINAL.Occu.Cand.models, modnames = libl16.FINAL.Occu.Modnames,newdata = dat.predFish, parm.type = "psi",conf.level = 0.85)
predictionslibl16AQ=modavgPred(cand.set = libl16.FINAL.Occu.Cand.models, modnames = libl16.FINAL.Occu.Modnames,newdata = dat.predAQ, parm.type = "psi",conf.level = 0.85)
predictionslibl16EM=modavgPred(cand.set = libl16.FINAL.Occu.Cand.models, modnames = libl16.FINAL.Occu.Modnames,newdata = dat.predEM, parm.type = "psi",conf.level = 0.85)

#PLOTTING IN GGPLOT_____________________________
#Convert ModAvgPred Output to class Data.Frame to be used in GGPLOT
libl16SLpreds=as.data.frame(predictionslibl16SL)
libl16Cattlepreds=as.data.frame(predictionslibl16Cattle)
libl16Fishpreds=as.data.frame(predictionslibl16Fish)
libl16AQpreds=as.data.frame(predictionslibl16AQ)
libl16EMpreds=as.data.frame(predictionslibl16EM)

#...LIBL 16 plots -------------------------------
(libl16sl <- ggplot(libl16SLpreds, aes(x=seq(from = 0.02, to = 0.35, by = 0.01), y=mod.avg.pred))+
    geom_line(data=libl16SLpreds)+
    geom_ribbon(data=libl16SLpreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
    theme_pubr()+
    scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
    scale_x_continuous(name="Pond Slope (rise:run)")+
  ggtitle("Plains leopard frog (2016)"))
(libl16Cattle <-ggplot(libl16Cattlepreds, aes(x=c(0,1), y=mod.avg.pred)) + 
    geom_bar(position=position_dodge(), stat="identity",alpha=0.3) +
    geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL),width=.2,# Width of the error bars
                  position=position_dodge(.9))
    +theme_pubr()+
    scale_x_discrete(name="Cattle")+
    scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
  ggtitle("Plains leopard frog (2016)"))

(libl16Fish <-ggplot(libl16Fishpreds, aes(x=c(0,1), y=mod.avg.pred)) + 
    geom_bar(position=position_dodge(), stat="identity",alpha=0.3) +
    geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL),width=.2,# Width of the error bars
                  position=position_dodge(.9))
    +theme_pubr()+
    scale_x_discrete(name="Fish")+
    scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
  ggtitle("Plains leopard frog (2016)"))

(libl16em <- ggplot(libl16EMpreds, aes(x=seq(from = 0, to = 65, by = 5), y=mod.avg.pred))+
    geom_line(data=libl16EMpreds)+
    geom_ribbon(data=libl16EMpreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
    theme_pubr()+
    scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
    scale_x_continuous(name="Emergent Vegetation Cover (%)")+
  ggtitle("Plains leopard frog (2016)"))

(libl16aq <- ggplot(libl16AQpreds, aes(x=seq(from = 0, to = 100, by = 5), y=mod.avg.pred))+
    geom_line(data=libl16AQpreds)+
    geom_ribbon(data=libl16AQpreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
    theme_pubr()+
    scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
    scale_x_continuous(name="Submerged Vegetation Cover (%)")+
  ggtitle("Plains leopard frog (2016)"))

#LICA16_51_ssom------------------------------------------------
lica16y<-data[,14:17]
n<-nrow(data)
#site level individual covariates
lica16.site<-data[,55:64]
#put everything together in unmarked data frame

lica16 <- unmarkedFrameOccu(y = lica16y, siteCovs = data[55:64],obsCovs=list(Julian=data[35:38],Temp=data[39:42],Time=data[51:54]))
summary(lica16) 
#Constant Detection and occupancy Model -->  follows format: model<-occu(~detection_formula ~occupancy_formula, dataframe)
#_______________________
#DETECTION
licad61=occu(~1~1,lica16)
licad62=occu(~Time~1,lica16)
licad63=occu(~Temp~1,lica16)
licad64=occu(~Julian~1,lica16)
licad65=occu(~Julian+Time~1,lica16)
licad66=occu(~Temp*Time~1,lica16)
licad67=occu(~Temp+Time~1,lica16)
licad68=occu(~Temp+Julian~1,lica16)
licad69=occu(~Temp*Julian~1,lica16)
#Goodness of Fit Test Bootstrap
#(gof_boot <- mb.gof.test(licad61, nsim = 500, plot.hist=FALSE))
#AICc Model Selection Table
lica16.Det.Cand.models=list(licad61,licad62,licad63,licad64,licad66,licad67,licad68,licad69)
lica16.Det.Modnames=c("ConstantDet","TimeDet","TempDet","JulianDet","Temp.by.Time","Temp.Time","Temp.Julian","Temp.by.Julian")
print(aictab(cand.set = lica16.Det.Cand.models,modnames = lica16.Det.Modnames,second.ord = TRUE), digits = 4)
#Generate Detection-correct Estimates
backTransform(licad61,'det') # for detection
backTransform(licad61,'state') # for occupancy
#___________________________________
#Stage 1 - Abiotic and Biotic
#NULL
lica16OccuNull=occu(~1~1,lica16)
#Abiotic
lica16OccuSL=occu(~1~SL16,lica16)
lica16OccuLogArea=occu(~1~LogArea,lica16)
lica16OccuCattle=occu(~1~Cattle16,lica16)
#Biotic
lica16OccuEM=occu(~1~EM16,lica16)
lica16OccuAQ=occu(~1~AQ16,lica16)
lica16OccuFL=occu(~1~FL16,lica16)
lica16OccuFish=occu(~1~Fish16,lica16)
lica16.Occu.Cand.models=list(lica16OccuSL,lica16OccuLogArea,lica16OccuCattle,lica16OccuEM,lica16OccuAQ,lica16OccuFL,lica16OccuFish,lica16OccuNull)
lica16.Occu.Modnames=c("SL","LogArea","Cattle","EM","AQ","FL","Fish","Null")
print(aictab(cand.set = lica16.Occu.Cand.models,modnames = lica16.Occu.Modnames,second.ord = TRUE), digits = 4)

# K     AICc Delta_AICc AICcWt Cum.Wt        LL
#LogArea 3 240.7361     0.0000 0.4515 0.4515 -117.1127
#FL      3 242.7546     2.0185 0.1646 0.6161 -118.1220
#Null    2 242.8527     2.1166 0.1567 0.7727 -119.3013
#SL      3 244.7916     4.0555 0.0594 0.8322 -119.1405
#EM      3 244.8414     4.1054 0.0580 0.8901 -119.1654
#Cattle  3 244.9825     4.2465 0.0540 0.9442 -119.2359
#Fish    3 245.0562     4.3201 0.0521 0.9962 -119.2728
#AQ      3 250.3048     9.5687 0.0038 1.0000 -121.8971

#_______________________________
#Stage 2 - Interactions
lica16OccuFL.Area=occu(~1~FL16+LogArea,lica16)
#FINAL CANDIDATE SET
lica16.FINAL.Occu.Cand.models=list(lica16OccuFL.Area,lica16OccuLogArea,lica16OccuFL,lica16OccuNull)
lica16.FINAL.Occu.Modnames=c("FL.Area","LogArea","FL","Null")
#FINAL AICc Table
print(aictab(cand.set = lica16.FINAL.Occu.Cand.models,modnames = lica16.FINAL.Occu.Modnames,second.ord = TRUE), digits = 4)

#..Final Model Set
#        K     AICc Delta_AICc AICcWt Cum.Wt        LL
#LogArea 3 240.7361     0.0000 0.4176 0.4176 -117.1127
#FL.Area 4 241.4984     0.7623 0.2853 0.7029 -116.3144
#FL      3 242.7546     2.0185 0.1522 0.8551 -118.1220
#Null    2 242.8527     2.1166 0.1449 1.0000 -119.3013


#MODEL AVERAGING __________________________
#Model averaging for each variable   
modavg(cand.set = lica16.FINAL.Occu.Cand.models, modnames = lica16.FINAL.Occu.Modnames, parm = "FL16",parm.type = "psi",conf.level = 0.85)
modavg(cand.set = lica16.FINAL.Occu.Cand.models, modnames = lica16.FINAL.Occu.Modnames, parm = "LogArea",parm.type = "psi",conf.level = 0.85)
#Create New data frames to store model averaged Predictions in 
dat.predFL<- data.frame(FL16 = seq(from = min(data$FL16),to = max(data$FL16), by = 5),LogArea = mean(data$LogArea))
dat.predArea <- data.frame(LogArea = seq(from = min(data$LogArea),to = max(data$LogArea), by = 0.1),FL16 = mean(data$FL16))
#Model-averaged predictions of psi across range of values
##of sitevar1 and entire model set
predictionslica16FL=modavgPred(cand.set = lica16.FINAL.Occu.Cand.models, modnames = lica16.FINAL.Occu.Modnames,newdata = dat.predFL, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionslica16EM,'predictionslica16EM.csv')
predictionslica16Area=modavgPred(cand.set = lica16.FINAL.Occu.Cand.models, modnames = lica16.FINAL.Occu.Modnames,newdata = dat.predArea, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionslica16LogArea,'predictionslica16LogArea.csv')

#PLOTTING IN GGPLOT___________________________
#Convert ModAvgPred Output to class Data.Frame to be used in GGPLOT
lica16FLpreds=as.data.frame(predictionslica16FL)
lica16LogAreapreds=as.data.frame(predictionslica16Area)
summary(data$FL16)
#write to csv
#write.csv(lica16FLpreds, file ="lica16FLpreds.csv")
#write.csv(lica16LogAreapreds, file ="lica16LogAreapreds.csv")

#...LICA 16 plots -------------------------------
(lica16fl <- ggplot(lica16FLpreds, aes(x=seq(from = 0, to = 100, by = 5), y=mod.avg.pred))+
    geom_line(data=lica16FLpreds)+
    geom_ribbon(data=lica16FLpreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
    theme_pubr()+
    scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
    scale_x_continuous(name="Floating Vegetation Cover (%)")+
    ggtitle("American bullfrog (2016)"))
(lica16logarea <- ggplot(lica16LogAreapreds, aes(x=seq(from = 2.3, to = 4.3, by = 0.1), y=mod.avg.pred))+
    geom_line(data=lica16LogAreapreds)+
    geom_ribbon(data=lica16LogAreapreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
    theme_pubr()+
    scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
    scale_x_continuous(name="Log Area")+
  ggtitle("American bullfrog (2016)"))


#HYLA17_51_SSOM---------------------------------------

hyla17y<-data[,18:21]
n<-nrow(data)
#site level individual covariates
hyla17.site<-data[,61:70]
#put everything together in unmarked data frame

hyla17 <- unmarkedFrameOccu(y = hyla17y, siteCovs = data[61:70],obsCovs=list(Julian=data[43:46],Temp=data[47:50],Time=data[51:54]))
summary(hyla17) 
#Constant Detection and occupancy Model -->  follows format: model<-occu(~detection_formula ~occupancy_formula, dataframe)
#_______________________
#DETECTION
hylad71=occu(~1~1,hyla17)
hylad72=occu(~Time~1,hyla17)
hylad73=occu(~Temp~1,hyla17)
hylad74=occu(~Julian~1,hyla17)
hylad75=occu(~Julian+Time~1,hyla17)
hylad76=occu(~Temp*Time~1,hyla17)
hylad77=occu(~Temp+Time~1,hyla17)
hylad78=occu(~Temp+Julian~1,hyla17)
hylad79=occu(~Temp*Julian~1,hyla17)
#Goodness of Fit Test Bootstrap
#(gof_boot <- mb.gof.test(hylad61, nsim = 500, plot.hist=FALSE))
#AICc Model Selection Table
hyla17.Det.Cand.models=list(hylad71,hylad72,hylad73,hylad74,hylad76,hylad77,hylad78,hylad79)
hyla17.Det.Modnames=c("ConstantDet","TimeDet","TempDet","JulianDet","Temp.by.Time","Temp.Time","Temp.Julian","Temp.by.Julian")
print(aictab(cand.set = hyla17.Det.Cand.models,modnames = hyla17.Det.Modnames,second.ord = TRUE), digits = 4)
#Generate Detection-correct Estimates
backTransform(hylad71,'det') # for detection
backTransform(hylad71,'state') # for occupancy

#___________________________________
#Stage 1 - Abiotic and Biotic
#NULL
hyla17OccuNull=occu(~1~1,hyla17)
#Abiotic
hyla17OccuSL=occu(~1~SL17,hyla17)
hyla17OccuLogArea=occu(~1~LogArea,hyla17)
hyla17OccupH=occu(~1~pH,hyla17)
hyla17OccuCattle=occu(~1~Cattle17,hyla17)
#Biotic
hyla17OccuEM=occu(~1~EM17,hyla17)
hyla17OccuAQ=occu(~1~AQ17,hyla17)
hyla17OccuFL=occu(~1~FL17,hyla17)
hyla17OccuFish=occu(~1~Fish17,hyla17)
hyla17.Occu.Cand.models=list(hyla17OccuSL,hyla17OccuLogArea,hyla17OccupH,hyla17OccuCattle,hyla17OccuEM,hyla17OccuAQ,hyla17OccuFL,hyla17OccuFish,hyla17OccuNull)
hyla17.Occu.Modnames=c("SL","LogArea","pH","Cattle","EM","AQ","FL","Fish","Null")
print(aictab(cand.set = hyla17.Occu.Cand.models,modnames = hyla17.Occu.Modnames,second.ord = TRUE), digits = 4)

#        K     AICc Delta_AICc AICcWt Cum.Wt       LL
#pH      3 122.4327     0.0000 0.9999 0.9999 -57.9611
#EM      3 141.0913    18.6585 0.0001 1.0000 -67.2903
#LogArea 3 145.3764    22.9436 0.0000 1.0000 -69.4329
#Fish    3 146.6883    24.2555 0.0000 1.0000 -70.0888
#SL      3 148.3873    25.9545 0.0000 1.0000 -70.9383
#Null    2 154.4907    32.0580 0.0000 1.0000 -75.1204
#FL      3 154.9380    32.5052 0.0000 1.0000 -74.2137
#AQ      3 156.5631    34.1304 0.0000 1.0000 -75.0262
#Cattle  3 156.6778    34.2450 0.0000 1.0000 -75.0836

#_______________________________
#Stage 2 - Interactions #Not needed, as pH model carries 0.9999 weight.

#FINAL CANDIDATE SET
hyla17.FINAL.Occu.Cand.models=list(hyla17OccupH,hyla17OccuNull)
hyla17.FINAL.Occu.Modnames=c("pH","Null")
#FINAL AICc Table
print(aictab(cand.set = hyla17.FINAL.Occu.Cand.models,modnames = hyla17.FINAL.Occu.Modnames,second.ord = TRUE), digits = 4)

#..Final Model Set
#      K     AICc Delta_AICc AICcWt Cum.Wt       LL
#pH   3 122.4327      0.000      1      1 -57.9611
#Null 2 154.4907     32.058      0      1 -75.1204

#MODEL AVERAGING __________________________
#Model averaging for each variable   
modavg(cand.set = hyla17.FINAL.Occu.Cand.models, modnames = hyla17.FINAL.Occu.Modnames, parm = "pH",parm.type = "psi",conf.level = 0.85)
#Create New data frames to store model averaged Predictions in 
dat.predEM<- data.frame(EM17 = seq(from = min(data$EM17),to = max(data$EM17), by = 5),pH = mean(data$pH))
dat.predpH <- data.frame(pH = seq(from = min(data$pH),to = max(data$pH), by = 0.1),EM17 = mean(data$EM17))
#Model-averaged predictions of psi across range of values
##of sitevar1 and entire model set
predictionshyla17EM=modavgPred(cand.set = hyla17.FINAL.Occu.Cand.models, modnames = hyla17.FINAL.Occu.Modnames,newdata = dat.predEM, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionshyla17EM,'predictionshyla17EM.csv')
predictionshyla17pH=modavgPred(cand.set = hyla17.FINAL.Occu.Cand.models, modnames = hyla17.FINAL.Occu.Modnames,newdata = dat.predpH, parm.type = "psi",conf.level = 0.85)
predictionshyla17pH95=modavgPred(cand.set = hyla17.FINAL.Occu.Cand.models, modnames = hyla17.FINAL.Occu.Modnames,newdata = dat.predpH, parm.type = "psi",conf.level = 0.95)

#write.csv(predictionshyla17LogArea,'predictionshyla17LogArea.csv')

#PLOTTING IN GGPLOT___________________________
#Convert ModAvgPred Output to class Data.Frame to be used in GGPLOT
hyla17EMpreds=as.data.frame(predictionshyla17EM)
hyla17pHpreds=as.data.frame(predictionshyla17pH)
hyla17pHpreds95=as.data.frame(predictionshyla17pH95)

summary(data$EM17)
#write to csv
#write.csv(hyla17EMpreds, file ="hyla17EMpreds.csv")
#write.csv(hyla17pHpreds, file ="hyla17pHpreds.csv")

#...HYLA 17 plots -------------------------------
(hyla17pH <- ggplot(hyla17pHpreds, aes(x=seq(from = 6.94, to = 10.95, by = 0.1), y=mod.avg.pred))+
    geom_line(data=hyla17pHpreds)+
    geom_ribbon(data=hyla17pHpreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
    theme_pubr()+
    scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
    scale_x_continuous(name="pH")+
  ggtitle("Gray treefrog spp. (2017)"))


#ACBL17_51_SSOM---------------------------------------

acbl17y<-data[,22:25]
n<-nrow(data)
#site level individual covariates
acbl17.site<-data[,61:70]
#put everything together in unmarked data frame

acbl17 <- unmarkedFrameOccu(y = acbl17y, siteCovs = data[61:70],obsCovs=list(Julian=data[43:46],Temp=data[47:50],Time=data[51:54]))
summary(acbl17) 
#Constant Detection and occupancy Model -->  follows format: model<-occu(~detection_formula ~occupancy_formula, dataframe)
#_______________________
#DETECTION
acbld71=occu(~1~1,acbl17)
acbld72=occu(~Time~1,acbl17)
acbld73=occu(~Temp~1,acbl17)
acbld74=occu(~Julian~1,acbl17)
acbld75=occu(~Julian+Time~1,acbl17)
acbld76=occu(~Temp*Time~1,acbl17)
acbld77=occu(~Temp+Time~1,acbl17)
acbld78=occu(~Temp+Julian~1,acbl17)
acbld79=occu(~Temp*Julian~1,acbl17)
#Goodness of Fit Test Bootstrap
#(gof_boot <- mb.gof.test(acbld61, nsim = 500, plot.hist=FALSE))
#AICc Model Selection Table
acbl17.Det.Cand.models=list(acbld71,acbld72,acbld73,acbld74,acbld76,acbld77,acbld78,acbld79)
acbl17.Det.Modnames=c("ConstantDet","TimeDet","TempDet","JulianDet","Temp.by.Time","Temp.Time","Temp.Julian","Temp.by.Julian")
print(aictab(cand.set = acbl17.Det.Cand.models,modnames = acbl17.Det.Modnames,second.ord = TRUE), digits = 4)
#Generate Detection-correct Estimates
backTransform(acbld71,'det') # for detection
backTransform(acbld74,'state') # for occupancy
#___________________________________
#Stage 1 - Abiotic and Biotic
#NULL
acbl17OccuNull=occu(~Julian~1,acbl17)
#Abiotic
acbl17OccuSL=occu(~Julian~SL17,acbl17)
acbl17OccuLogArea=occu(~Julian~LogArea,acbl17)
acbl17OccupH=occu(~Julian~pH,acbl17)
acbl17OccuCattle=occu(~Julian~Cattle17,acbl17)
#Biotic
acbl17OccuEM=occu(~Julian~EM17,acbl17)
acbl17OccuAQ=occu(~Julian~AQ17,acbl17)
acbl17OccuFL=occu(~Julian~FL17,acbl17)
acbl17OccuFish=occu(~Julian~Fish17,acbl17)
acbl17.Occu.Cand.models=list(acbl17OccuSL,acbl17OccuLogArea,acbl17OccupH,acbl17OccuCattle,acbl17OccuEM,acbl17OccuAQ,acbl17OccuFL,acbl17OccuFish,acbl17OccuNull)
acbl17.Occu.Modnames=c("SL","LogArea","pH","Cattle","EM","AQ","FL","Fish","Null")
print(aictab(cand.set = acbl17.Occu.Cand.models,modnames = acbl17.Occu.Modnames,second.ord = TRUE), digits = 4)

#        K     AICc Delta_AICc AICcWt Cum.Wt       LL
#EM      4 117.2362     0.0000 0.7931 0.7931 -54.1833
#pH      4 122.6483     5.4121 0.0530 0.8461 -56.8894
#FL      4 122.7421     5.5059 0.0506 0.8967 -56.9362
#AQ      4 123.8245     6.5883 0.0294 0.9261 -57.4775
#Null    3 123.8367     6.6006 0.0292 0.9553 -58.6630
#SL      4 125.3883     8.1521 0.0135 0.9688 -58.2593
#Fish    4 125.7237     8.4876 0.0114 0.9802 -58.4271
#Cattle  4 125.9355     8.6993 0.0102 0.9904 -58.5330
#LogArea 4 126.0674     8.8312 0.0096 1.0000 -58.5989

#_______________________________
#Stage 2 - Additive, all-subsets models
acbl17OccuEM.pH.FL.AQ=occu(~Julian~EM17+pH+FL17+AQ17,acbl17)

acbl17OccuEM.pH.FL=occu(~Julian~EM17+pH+FL17,acbl17)
acbl17OccupH.FL.AQ=occu(~Julian~pH+FL17+AQ17,acbl17)
acbl17OccuEM.FL.AQ=occu(~Julian~EM17+FL17+AQ17,acbl17)
acbl17OccuEM.pH.AQ=occu(~Julian~EM17+pH+AQ17,acbl17)

acbl17OccuEM.pH   =occu(~Julian~EM17+pH,acbl17)
acbl17OccuEM.FL   =occu(~Julian~EM17+FL17,acbl17)
acbl17OccuEM.AQ   =occu(~Julian~EM17+AQ17,acbl17)

acbl17OccupH.AQ   =occu(~Julian~pH+AQ17,acbl17)
acbl17OccupH.FL   =occu(~Julian~pH+FL17,acbl17)

acbl17OccuFL.AQ   =occu(~Julian~FL17+AQ17,acbl17)

acbl17OccuEM      =occu(~Julian~EM17,acbl17)
acbl17OccuAQ      =occu(~Julian~AQ17,acbl17)
acbl17OccuFL      =occu(~Julian~FL17,acbl17)
acbl17OccupH      =occu(~Julian~pH,acbl17)


#FINAL CANDIDATE SET
acbl17.FINAL.Occu.Cand.models=list(acbl17OccuEM.pH.FL.AQ,acbl17OccuEM.pH.FL, acbl17OccupH.FL.AQ, acbl17OccuEM.FL.AQ, acbl17OccuEM.pH.AQ, acbl17OccuEM.pH, acbl17OccuEM.FL,acbl17OccuEM.AQ,acbl17OccupH.AQ,acbl17OccupH.FL,acbl17OccuFL.AQ,acbl17OccuEM,acbl17OccuAQ,acbl17OccuFL,acbl17OccupH,acbl17OccuNull)
acbl17.FINAL.Occu.Modnames=c("EM.pH.FL.AQ","EM.pH.FL", "pH.FL.AQ","EM.FL.AQ","EM.pH.AQ","EM.pH","EM.FL","EM.AQ","pH.AQ","pH.FL","FL.AQ","EM","AQ","FL","pH","Null")
#FINAL AICc Table
print(aictab(cand.set = acbl17.FINAL.Occu.Cand.models,modnames = acbl17.FINAL.Occu.Modnames,second.ord = TRUE), digits = 4)

#..Final Model Set
#            K     AICc Delta_AICc AICcWt Cum.Wt       LL
#EM.FL       5 116.6566     0.0000 0.2330 0.2330 -52.6616
#EM          4 117.2362     0.5796 0.1744 0.4074 -54.1833
#EM.AQ       5 117.7970     1.1404 0.1317 0.5391 -53.2318
#EM.pH.AQ    6 118.4617     1.8051 0.0945 0.6336 -52.2763
#EM.FL.AQ    6 118.7660     2.1094 0.0812 0.7148 -52.4285
#EM.pH.FL    6 118.8851     2.2285 0.0765 0.7912 -52.4880
#EM.pH       5 118.9835     2.3269 0.0728 0.8640 -53.8251
#pH.AQ       5 120.4086     3.7521 0.0357 0.8997 -54.5376
#EM.pH.FL.AQ 7 120.4517     3.7951 0.0349 0.9346 -51.9235
#pH.FL       5 122.4111     5.7546 0.0131 0.9478 -55.5389
#pH.FL.AQ    6 122.5715     5.9149 0.0121 0.9599 -54.3312
#pH          4 122.6483     5.9917 0.0116 0.9715 -56.8894
#FL          4 122.7421     6.0855 0.0111 0.9826 -56.9362
#AQ          4 123.8245     7.1679 0.0065 0.9891 -57.4775
#Null        3 123.8367     7.1802 0.0064 0.9955 -58.6630
#FL.AQ       5 124.5600     7.9034 0.0045 1.0000 -56.6133


#MODEL AVERAGING __________________________
#Model averaging for each variable   
modavg(cand.set = acbl17.FINAL.Occu.Cand.models, modnames = acbl17.FINAL.Occu.Modnames, parm = "EM17",parm.type = "psi",conf.level = 0.85)
modavg(cand.set = acbl17.FINAL.Occu.Cand.models, modnames = acbl17.FINAL.Occu.Modnames, parm = "pH",parm.type = "psi",conf.level = 0.85)
modavg(cand.set = acbl17.FINAL.Occu.Cand.models, modnames = acbl17.FINAL.Occu.Modnames, parm = "FL17",parm.type = "psi",conf.level = 0.85)
modavg(cand.set = acbl17.FINAL.Occu.Cand.models, modnames = acbl17.FINAL.Occu.Modnames, parm = "AQ17",parm.type = "psi",conf.level = 0.85)

#Create New data frames to store model averaged Predictions in 
dat.predEM<- data.frame(EM17 = seq(from = min(data$EM17),to = max(data$EM17), by = 5),pH = mean(data$pH), AQ17=mean(data$AQ17),FL17=mean(data$FL17))
dat.predpH <- data.frame(pH = seq(from = min(data$pH),to = max(data$pH), by = 0.1),EM17 = mean(data$EM17), AQ17=mean(data$AQ17),FL17=mean(data$FL17))
dat.predaq <- data.frame(AQ17 = seq(from = min(data$AQ17),to = max(data$AQ17), by = 5),EM17 = mean(data$EM17),pH = mean(data$pH), FL17=mean(data$FL17))
dat.predfl <- data.frame(FL17 = seq(from = min(data$FL17),to = max(data$FL17), by = 5),EM17 = mean(data$EM17),pH = mean(data$pH), AQ17=mean(data$AQ17))

#Model-averaged predictions of psi across range of values
##of sitevar1 and entire model set
predictionsacbl17EM=modavgPred(cand.set = acbl17.FINAL.Occu.Cand.models, modnames = acbl17.FINAL.Occu.Modnames,newdata = dat.predEM, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionsacbl17EM,'predictionsacbl17EM.csv')
predictionsacbl17pH=modavgPred(cand.set = acbl17.FINAL.Occu.Cand.models, modnames = acbl17.FINAL.Occu.Modnames,newdata = dat.predpH, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionsacbl17LogArea,'predictionsacbl17LogArea.csv')
predictionsacbl17aq=modavgPred(cand.set = acbl17.FINAL.Occu.Cand.models, modnames = acbl17.FINAL.Occu.Modnames,newdata = dat.predaq, parm.type = "psi",conf.level = 0.85)
predictionsacbl17fl=modavgPred(cand.set = acbl17.FINAL.Occu.Cand.models, modnames = acbl17.FINAL.Occu.Modnames,newdata = dat.predfl, parm.type = "psi",conf.level = 0.85)

#PLOTTING IN GGPLOT___________________________
#Convert ModAvgPred Output to class Data.Frame to be used in GGPLOT
acbl17EMpreds=as.data.frame(predictionsacbl17EM)
acbl17pHpreds=as.data.frame(predictionsacbl17pH)
acbl17aqpreds=as.data.frame(predictionsacbl17aq)
acbl17flpreds=as.data.frame(predictionsacbl17fl)

summary(data$EM17)
summary(data$AQ17)
summary(data$FL17)

#write to csv
#write.csv(acbl17EMpreds, file ="acbl17EMpreds.csv")
#write.csv(acbl17pHpreds, file ="acbl17pHpreds.csv")

#...ACBL 17 plots -------------------------------
(acbl17em <- ggplot(acbl17EMpreds, aes(x=seq(from = 0, to = 70, by = 5), y=mod.avg.pred))+
   geom_line(data=acbl17EMpreds)+
   geom_ribbon(data=acbl17EMpreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
   theme_pubr()+
   scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
   scale_x_continuous(name="Emergent Vegetation Cover (%)")+
  ggtitle("Blanchard's cricket frog (2017)"))

(acbl17pH <- ggplot(acbl17pHpreds, aes(x=seq(from = 6.94, to = 10.95, by = 0.1), y=mod.avg.pred))+
    geom_line(data=acbl17pHpreds)+
    geom_ribbon(data=acbl17pHpreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
    theme_pubr()+
    scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
    scale_x_continuous(name="pH")+
  ggtitle("Blanchard's cricket frog (2017)"))

(acbl17aq <- ggplot(acbl17aqpreds, aes(x=seq(from = 0, to = 70, by = 5), y=mod.avg.pred))+
    geom_line(data=acbl17aqpreds)+
    geom_ribbon(data=acbl17aqpreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
    theme_pubr()+
    scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
    scale_x_continuous(name="Submerged Vegetation Cover (%)")+
  ggtitle("Blanchard's cricket frog (2017)"))

(acbl17fl <- ggplot(acbl17flpreds, aes(x=seq(from = 0, to = 95, by = 5), y=mod.avg.pred))+
    geom_line(data=acbl17flpreds)+
    geom_ribbon(data=acbl17flpreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
    theme_pubr()+
    scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
    scale_x_continuous(name="Floating Vegetation Cover (%)")+
ggtitle("Blanchard's cricket frog (2017)"))

#LIBL17_51_SSOM---------------------------------------

libl17y<-data[,26:29]
n<-nrow(data)
#site level individual covariates
libl17.site<-data[,61:70]
#put everything together in unmarked data frame

libl17 <- unmarkedFrameOccu(y = libl17y, siteCovs = data[61:70],obsCovs=list(Julian=data[43:46],Temp=data[47:50],Time=data[51:54]))
summary(libl17) 
#Constant Detection and occupancy Model -->  follows format: model<-occu(~detection_formula ~occupancy_formula, dataframe)
#_______________________
#DETECTION
libld71=occu(~1~1,libl17)
libld72=occu(~Time~1,libl17)
libld73=occu(~Temp~1,libl17)
libld74=occu(~Julian~1,libl17)
libld75=occu(~Julian+Time~1,libl17)
libld76=occu(~Temp*Time~1,libl17)
libld77=occu(~Temp+Time~1,libl17)
libld78=occu(~Temp+Julian~1,libl17)
libld79=occu(~Temp*Julian~1,libl17)
#Goodness of Fit Test Bootstrap
#(gof_boot <- mb.gof.test(libld61, nsim = 500, plot.hist=FALSE))
#AICc Model Selection Table
libl17.Det.Cand.models=list(libld71,libld72,libld73,libld74,libld76,libld77,libld78,libld79)
libl17.Det.Modnames=c("ConstantDet","TimeDet","TempDet","JulianDet","Temp.by.Time","Temp.Time","Temp.Julian","Temp.by.Julian")
print(aictab(cand.set = libl17.Det.Cand.models,modnames = libl17.Det.Modnames,second.ord = TRUE), digits = 4)
#Generate Detection-correct Estimates
backTransform(libld71,'det') # for detection
backTransform(libld78,'state') # for occupancy
#___________________________________
#Stage 1 - Abiotic and Biotic
#NULL
libl17OccuNull=occu(~Temp+Julian~1,libl17)
#Abiotic
libl17OccuSL=occu(~Temp+Julian~SL17,libl17)
libl17OccuLogArea=occu(~Temp+Julian~LogArea,libl17)
libl17OccupH=occu(~Temp+Julian~pH,libl17)
libl17OccuCattle=occu(~Temp+Julian~Cattle17,libl17)
#Biotic
libl17OccuEM=occu(~Temp+Julian~EM17,libl17)
libl17OccuAQ=occu(~Temp+Julian~AQ17,libl17)
libl17OccuFL=occu(~Temp+Julian~FL17,libl17)
libl17OccuFish=occu(~Temp+Julian~Fish17,libl17)
libl17.Occu.Cand.models=list(libl17OccuSL,libl17OccuLogArea,libl17OccupH,libl17OccuCattle,libl17OccuEM,libl17OccuAQ,libl17OccuFL,libl17OccuFish,libl17OccuNull)
libl17.Occu.Modnames=c("SL","LogArea","pH","Cattle","EM","AQ","FL","Fish","Null")
print(aictab(cand.set = libl17.Occu.Cand.models,modnames = libl17.Occu.Modnames,second.ord = TRUE), digits = 4)

#K     AICc Delta_AICc AICcWt Cum.Wt       LL
#AQ      5 132.3933     0.0000 0.7442 0.7442 -60.5300
#Fish    5 137.1436     4.7503 0.0692 0.8134 -62.9051
#Null    4 137.6138     5.2205 0.0547 0.8681 -64.3721
#EM      5 138.7118     6.3185 0.0316 0.8997 -63.6892
#LogArea 5 138.9886     6.5953 0.0275 0.9272 -63.8276
#FL      5 139.5081     7.1148 0.0212 0.9484 -64.0874
#Cattle  5 139.8002     7.4069 0.0183 0.9667 -64.2334
#SL      5 139.9195     7.5262 0.0173 0.9840 -64.2931
#pH      5 140.0735     7.6802 0.0160 1.0000 -64.3701

#_______________________________
#Stage 2 - Interactions
libl17OccuAQ.Fish=occu(~Temp+Julian~AQ17+Fish17,libl17)
#FINAL CANDIDATE SET
libl17.FINAL.Occu.Cand.models=list(libl17OccuAQ.Fish,libl17OccuAQ,libl17OccuFish,libl17OccuNull)
libl17.FINAL.Occu.Modnames=c("AQ.Fish","AQ","Fish","Null")
#FINAL AICc Table
print(aictab(cand.set = libl17.FINAL.Occu.Cand.models,modnames = libl17.FINAL.Occu.Modnames,second.ord = TRUE), digits = 4)

#..Final model set
#        K     AICc Delta_AICc AICcWt Cum.Wt       LL
#AQ      5 132.3933     0.0000 0.5900 0.5900 -60.5300
#AQ.Fish 6 133.6692     1.2759 0.3117 0.9018 -59.8800
#Fish    5 137.1436     4.7503 0.0549 0.9566 -62.9051
#Null    4 137.6138     5.2205 0.0434 1.0000 -64.3721


#MODEL AVERAGING __________________________
#Model averaging for each variable   
modavg(cand.set = libl17.FINAL.Occu.Cand.models, modnames = libl17.FINAL.Occu.Modnames, parm = "AQ17",parm.type = "psi",conf.level = 0.85)
modavg(cand.set = libl17.FINAL.Occu.Cand.models, modnames = libl17.FINAL.Occu.Modnames, parm = "Fish17",parm.type = "psi",conf.level = 0.85)
#Create New data frames to store model averaged Predictions in 
dat.predAQ<- data.frame(AQ17 = seq(from = min(data$AQ17),to = max(data$AQ17), by = 5),Fish17 = mean(data$Fish17))
#dat.predFish <- data.frame(Fish17 = seq(from = min(data$Fish17),to = max(data$Fish17), by = 0.1),AQ17 = mean(data$AQ17))
dat.predFish<-data.frame(Fish17 = factor(c("0", "1")),AQ17 = mean(data$SL17))
#Model-averaged predictions of psi across range of values
##of sitevar1 and entire model set
predictionslibl17AQ=modavgPred(cand.set = libl17.FINAL.Occu.Cand.models, modnames = libl17.FINAL.Occu.Modnames,newdata = dat.predAQ, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionslibl17EM,'predictionslibl17EM.csv')
predictionslibl17Fish=modavgPred(cand.set = libl17.FINAL.Occu.Cand.models, modnames = libl17.FINAL.Occu.Modnames,newdata = dat.predFish, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionslibl17LogArea,'predictionslibl17LogArea.csv')

#PLOTTING IN GGPLOT___________________________
#Convert ModAvgPred Output to class Data.Frame to be used in GGPLOT
libl17AQpreds=as.data.frame(predictionslibl17AQ)
libl17Fishpreds=as.data.frame(predictionslibl17Fish)
summary(data$AQ17)
#write to csv
#write.csv(libl17AQpreds, file ="libl17AQpreds.csv")
#write.csv(libl17Fishpreds, file ="libl17Fishpreds.csv")

#...LIBL 17 plots -------------------------------
(libl17aq <- ggplot(libl17AQpreds, aes(x=seq(from = 0, to = 70, by = 5), y=mod.avg.pred))+
   geom_line(data=libl17AQpreds)+
   geom_ribbon(data=libl17AQpreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
   theme_pubr()+
   scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
   scale_x_continuous(name="Submerged Vegetation Cover (%)")+
  ggtitle("Plains leopard frog (2017)"))
(libl17Fish <-ggplot(libl17Fishpreds, aes(x=c(0,1), y=mod.avg.pred)) + 
    geom_bar(position=position_dodge(), stat="identity",alpha=0.3) +
    geom_errorbar(aes(ymin=lower.CL, ymax=upper.CL),width=.2,# Width of the error bars
                  position=position_dodge(.9))
    +theme_pubr()+
    scale_x_discrete(name="Fish")+
    scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
  ggtitle("Plains leopard frog (2017)"))
  
#LICA17_51_SSOM---------------------------------------

lica17y<-data[,30:33]
n<-nrow(data)
#site level individual covariates
lica17.site<-data[,61:70]
#put everything together in unmarked data frame

lica17 <- unmarkedFrameOccu(y = lica17y, siteCovs = data[61:70],obsCovs=list(Julian=data[43:46],Temp=data[47:50],Time=data[51:54]))
summary(lica17) 
#Constant Detection and occupancy Model -->  follows format: model<-occu(~detection_formula ~occupancy_formula, dataframe)
#_______________________
#DETECTION
licad71=occu(~1~1,lica17)
licad72=occu(~Time~1,lica17)
licad73=occu(~Temp~1,lica17)
licad74=occu(~Julian~1,lica17)
licad75=occu(~Julian+Time~1,lica17)
licad76=occu(~Temp*Time~1,lica17)
licad77=occu(~Temp+Time~1,lica17)
licad78=occu(~Temp+Julian~1,lica17)
licad79=occu(~Temp*Julian~1,lica17)
#Goodness of Fit Test Bootstrap
#(gof_boot <- mb.gof.test(licad61, nsim = 500, plot.hist=FALSE))
#AICc Model Selection Table
lica17.Det.Cand.models=list(licad71,licad72,licad73,licad74,licad76,licad77,licad78,licad79)
lica17.Det.Modnames=c("ConstantDet","TimeDet","TempDet","JulianDet","Temp.by.Time","Temp.Time","Temp.Julian","Temp.by.Julian")
print(aictab(cand.set = lica17.Det.Cand.models,modnames = lica17.Det.Modnames,second.ord = TRUE), digits = 4)
#Generate Detection-correct Estimates
backTransform(licad71,'det') # for detection
backTransform(licad71,'state') # for occupancy
#___________________________________
#Stage 1 - Abiotic and Biotic
#NULL
lica17OccuNull=occu(~1~1,lica17)
#Abiotic
lica17OccuSL=occu(~1~SL17,lica17)
lica17OccuLogArea=occu(~1~LogArea,lica17)
lica17OccupH=occu(~1~pH,lica17)
lica17OccuCattle=occu(~1~Cattle17,lica17)
#Biotic
lica17OccuEM=occu(~1~EM17,lica17)
lica17OccuAQ=occu(~1~AQ17,lica17)
lica17OccuFL=occu(~1~FL17,lica17)
lica17OccuFish=occu(~1~Fish17,lica17)
lica17.Occu.Cand.models=list(lica17OccuSL,lica17OccuLogArea,lica17OccupH,lica17OccuCattle,lica17OccuEM,lica17OccuAQ,lica17OccuFL,lica17OccuFish,lica17OccuNull)
lica17.Occu.Modnames=c("SL","LogArea","pH","Cattle","EM","AQ","FL","Fish","Null")
print(aictab(cand.set = lica17.Occu.Cand.models,modnames = lica17.Occu.Modnames,second.ord = TRUE), digits = 4)

#        K     AICc Delta_AICc AICcWt Cum.Wt        LL
#LogArea 3 241.9236     0.0000 0.8216 0.8216 -117.7065
#FL      3 247.6046     5.6810 0.0480 0.8696 -120.5470
#Null    2 248.0657     6.1421 0.0381 0.9077 -121.9079
#pH      3 249.3876     7.4639 0.0197 0.9274 -121.4385
#Fish    3 249.6430     7.7194 0.0173 0.9447 -121.5662
#EM      3 249.7400     7.8164 0.0165 0.9612 -121.6147
#Cattle  3 250.1384     8.2148 0.0135 0.9747 -121.8139
#AQ      3 250.2392     8.3155 0.0129 0.9876 -121.8643
#SL      3 250.3044     8.3808 0.0124 1.0000 -121.8969

#_______________________________
#Stage 2 - Interactions
lica17OccuFL.Area=occu(~1~LogArea+FL17,lica17)
#FINAL CANDIDATE SET
lica17.FINAL.Occu.Cand.models=list(lica17OccuFL.Area,lica17OccuLogArea,lica17OccuFL,lica17OccuNull)
lica17.FINAL.Occu.Modnames=c("FL.Area","LogArea","FL","Null")
#FINAL AICc Table
print(aictab(cand.set = lica17.FINAL.Occu.Cand.models,modnames = lica17.FINAL.Occu.Modnames,second.ord = TRUE), digits = 4)

#..Final Model Set
#        K     AICc Delta_AICc AICcWt Cum.Wt        LL
#FL.Area 4 241.5213     0.0000 0.5254 0.5254 -116.3259
#LogArea 3 241.9236     0.4023 0.4296 0.9550 -117.7065
#FL      3 247.6046     6.0833 0.0251 0.9801 -120.5470
#Null    2 248.0657     6.5444 0.0199 1.0000 -121.9079

#MODEL AVERAGING __________________________
#Model averaging for each variable  
modavg(cand.set = lica17.FINAL.Occu.Cand.models, modnames = lica17.FINAL.Occu.Modnames, parm = "FL17",parm.type = "psi",conf.level = 0.85)
modavg(cand.set = lica17.FINAL.Occu.Cand.models, modnames = lica17.FINAL.Occu.Modnames, parm = "LogArea",parm.type = "psi",conf.level = 0.85)
#Create New data frames to store model averaged Predictions in 
dat.predFL<- data.frame(FL17 = seq(from = min(data$FL17),to = max(data$FL17), by = 5),LogArea = mean(data$LogArea))
dat.predLogArea <- data.frame(LogArea = seq(from = min(data$LogArea),to = max(data$LogArea), by = 0.1),FL17 = mean(data$FL17))
#Model-averaged predictions of psi across range of values
##of sitevar1 and entire model set
predictionslica17Area=modavgPred(cand.set = lica17.FINAL.Occu.Cand.models, modnames = lica17.FINAL.Occu.Modnames,newdata = dat.predLogArea, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionslica17EM,'predictionslica17EM.csv')
predictionslica17FL=modavgPred(cand.set = lica17.FINAL.Occu.Cand.models, modnames = lica17.FINAL.Occu.Modnames,newdata = dat.predFL, parm.type = "psi",conf.level = 0.85)
#write.csv(predictionslica17LogArea,'predictionslica17LogArea.csv')

#PLOTTING IN GGPLOT___________________________
#Convert ModAvgPred Output to class Data.Frame to be used in GGPLOT
lica17Areapreds=as.data.frame(predictionslica17Area)
lica17FLpreds=as.data.frame(predictionslica17FL)
summary(data$LogArea)
summary(data$FL17)
#write to csv
#write.csv(lica17Areapreds, file ="lica17Areapreds.csv")
#write.csv(lica17FLpreds, file ="lica17FLpreds.csv")

#...LICA 17 plots -------------------------------
(lica17fl <- ggplot(lica17FLpreds, aes(x=seq(from = 0, to = 95, by = 5), y=mod.avg.pred))+
   geom_line(data=lica17FLpreds)+
   geom_ribbon(data=lica17FLpreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
   theme_pubr()+
   scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
   scale_x_continuous(name="Floating Vegetation Cover (%)")+
  ggtitle("American Bullfrog (2017)"))

(lica17area <- ggplot(lica17Areapreds, aes(x=seq(from = 2.3, to = 4.3, by = 0.1), y=mod.avg.pred))+
    geom_line(data=lica17Areapreds)+
    geom_ribbon(data=lica17Areapreds,aes(ymin=lower.CL,ymax=upper.CL),alpha=0.3)+
    theme_pubr()+
    scale_y_continuous(name="Breeding Occupancy Probability",limits = c(0,1))+
    scale_x_continuous(name="Log Area")+
    ggtitle("American Bullfrog (2017)"))


#Final Figures--------------------------------------
ggarrange( acbl16aq,acbl16em,acbl17em,
           hyla16em,hyla16logarea,hyla17pH,
           libl16Cattle,libl16em,libl17aq,
           lica16logarea,lica17area,ncol = 3, nrow = 4)




