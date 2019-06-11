data=read.csv("DataS2_Chronosequence_Analysis.csv")

#Packages to load---------------------------------
library(ggplot2)
library(ggpubr)
library(dplyr)
library(gridExtra)
library(multcomp)
library(tweedie)
library(statmod)

#For each variable, the year of the measurement was taken is noted as either a "1617" to indicate that 
#measurements from 2016 and 2017 were averaged for this comparison. pH was only measured in 2017, so no year is listed.
#The following abbreviations are used throughout:
  #fl= Floating vegetation percent cover
  #em= Emergent vegetation percent cover
  #aq= submerged/aquatic vegetation percent cover
  #cattle= presence (1) or absence (0) of cattle around the pond

#Dependent Variables (Habitat)
sl1617=data$SL1617
em1617=data$EM1617
fl1617=data$FL1617
aq1617=data$AQ1617
ph=data$pH

#Independent Variables
Age=data$AgeCategory
Cattle1617=data$Cattle1617
treesedge=data$TreesEdge


#Below, glms are generated for each response variable. for variables that are not normally distributed, we use Tweedie distributions https://cran.r-project.org/web/packages/tweedie/tweedie.pdf#

#estimate index parameter p, a number between 1 and 2 to determine shape of tweedie distribution
xi.vec <- seq(1.01, 1.99, by=0.1) #vector containing possible values for p

#GLM- Emergent Vegetation####
out.em <- tweedie.profile(em1617~1, xi.vec=xi.vec, do.plot=TRUE, verbose=TRUE)
#Emergent Vegetation glm
m1.em=(glm(em1617 ~ Age, family=tweedie(var.power=out.em$xi.max, link.power=0) ))
m2.em<- glm(em1617~Age)
m3.em=(glm(em1617 ~ Age+Cattle1617+treesedge, family=tweedie(var.power=out.em$xi.max, link.power=0) ))
m4.em<- glm(em1617~Age+Cattle1617+treesedge)

AICtweedie(m1.em) #395.1952
AIC(m2.em) #432.2528
AICtweedie(m3.em) #395.1208
AIC(m4.em) #429.9347

summary(m3.em) #effects for Cattle and Edge Tree cover
#             Estimate Std. Error t value Pr(>|t|)  
#Cattle1617   -0.5270     0.3235  -1.629   0.1104    
#treesedge    -0.8961     0.4251  -2.108   0.0408 *  


tukEm=glht(m3.em,mcp(Age="Tukey"))
tuk.cld.em<- cld(tukEm)
plot(tuk.cld.em) #derive letter assignments from plot for significance of differences between categories (alpha = 0.05)
#<30 30-39 40-49 50-59 60+
#a   a     a     a     a
summary(tukEm)
#Multiple Comparisons of Means: Tukey Contrasts
#                    Estimate Std. Error z value Pr(>|z|)
#30-39 - <30 == 0   -0.08200    0.50058  -0.164    1.000
#40-49 - <30 == 0    0.36525    0.43085   0.848    0.914
#50-59 - <30 == 0    0.45289    0.37723   1.201    0.748
#60+ - <30 == 0      0.95269    0.41522   2.294    0.144
#40-49 - 30-39 == 0  0.44725    0.51195   0.874    0.905
#50-59 - 30-39 == 0  0.53489    0.46757   1.144    0.780
#60+ - 30-39 == 0    1.03469    0.49905   2.073    0.228
#50-59 - 40-49 == 0  0.08764    0.39132   0.224    0.999
#60+ - 40-49 == 0    0.58744    0.42873   1.370    0.643
#60+ - 50-59 == 0    0.49980    0.37468   1.334    0.666
#(Adjusted p values reported -- single-step method)


#GLM- Floating Vegetation Cover (%)####
out.fl <- tweedie.profile(fl1617~1, xi.vec=xi.vec, do.plot=TRUE, verbose=TRUE)
# Fit the glm
m1.fl=(glm(fl1617 ~ Age, family=tweedie(var.power=out.fl$xi.max, link.power=0) ))
m2.fl<- glm(fl1617~Age)
m3.fl=(glm(fl1617 ~ Age+Cattle1617+treesedge, family=tweedie(var.power=out.fl$xi.max, link.power=0) ))
m4.fl<- glm(fl1617~Age+Cattle1617+treesedge)

AICtweedie(m1.fl) #313.9577
AIC(m2.fl) #473.6372
AICtweedie(m3.fl) #316.7752
AIC(m4.fl) #476.3082

summary(m3.fl) #effects for Cattle and Edge Tree cover
#             Estimate Std. Error t value Pr(>|t|)  
#Cattle1617   -0.5724     0.6111  -0.937  0.35403   
#treesedge     0.2207     0.7210   0.306  0.76094 

tukfl=(glht(m3.fl,mcp(Age="Tukey"))) #using normal dist model (m2) since normal and AIC not diff
summary(tukfl)
tuk.cld.fl<- cld(tukfl)
plot(tuk.cld.fl) #derive letter assignments from plot for significance of differences between categories (alpha = 0.05)
#<30 30-39 40-49 50-59 60+
#a   a     a     a     a

#Multiple Comparisons of Means: Tukey Contrasts
#                    Estimate Std. Error z value Pr(>|z|)
#30-39 - <30 == 0    0.768924   0.870106   0.884    0.902
#40-49 - <30 == 0    0.521700   0.801463   0.651    0.966
#50-59 - <30 == 0    0.774503   0.698205   1.109    0.799
#60+ - <30 == 0      0.331478   0.817126   0.406    0.994
#40-49 - 30-39 == 0 -0.247224   0.869103  -0.284    0.999
#50-59 - 30-39 == 0  0.005579   0.773704   0.007    1.000
#60+ - 30-39 == 0   -0.437446   0.883588  -0.495    0.988
#50-59 - 40-49 == 0  0.252803   0.694500   0.364    0.996
#60+ - 40-49 == 0   -0.190222   0.815601  -0.233    0.999
#60+ - 50-59 == 0   -0.443025   0.712580  -0.622    0.971
#(Adjusted p values reported -- single-step method)


#GLM- Submerged Vegetation Cover (%)####
out.aq <- tweedie.profile( aq1617~1, xi.vec=xi.vec, do.plot=TRUE, verbose=TRUE)
# Fit the glm
m1.aq=(glm(aq1617 ~ Age, family=tweedie(var.power=out.aq$xi.max, link.power=0) ))
m2.aq<- glm(aq1617~Age)
m3.aq=(glm(aq1617 ~ Age+Cattle1617+treesedge, family=tweedie(var.power=out.aq$xi.max, link.power=0) ))
m4.aq=glm(aq1617~Age+Cattle1617+treesedge)

AICtweedie(m1.aq) #288.3357
AIC(m2.aq) #427.7133
AICtweedie(m3.aq) #285.7538
AIC(m4.aq) #425.6009

summary(m3.aq) #effects of cattle and edge tree cover 
#             Estimate Std. Error t value Pr(>|t|)  
#Cattle1617  -1.47243    0.56842  -2.590   0.0130 *  
#treesedge   -1.64848    0.74559  -2.211   0.0323 * 

tukaq=(glht(m3.aq,mcp(Age="Tukey"))) #using normal dist model (m2) since normal and AIC not diff
summary(tukaq)
tuk.cld.aq<- cld(tukaq)
plot(tuk.cld.aq) #derive letter assignments from plot for significance of differences between categories (alpha = 0.05)
#<30 30-39 40-49 50-59 60+
#a   a     a     a     a

#Multiple Comparisons of Means: Tukey Contrasts
#                   Estimate Std. Error z value Pr(>|z|)
#30-39 - <30 == 0   -0.49342    0.85036  -0.580    0.977
#40-49 - <30 == 0    0.30884    0.67735   0.456    0.991
#50-59 - <30 == 0    0.06062    0.60853   0.100    1.000
#60+ - <30 == 0     -0.22941    0.72095  -0.318    0.998
#40-49 - 30-39 == 0  0.80225    0.86608   0.926    0.884
#50-59 - 30-39 == 0  0.55404    0.81357   0.681    0.960
#60+ - 30-39 == 0    0.26401    0.89951   0.294    0.998
#50-59 - 40-49 == 0 -0.24821    0.62883  -0.395    0.995
#60+ - 40-49 == 0   -0.53824    0.73863  -0.729    0.949
#60+ - 50-59 == 0   -0.29003    0.67651  -0.429    0.993
#(Adjusted p values reported -- single-step method)

#GLM- ph#### 
out.ph <- tweedie.profile( ph~1, xi.vec=xi.vec, do.plot=TRUE, verbose=TRUE)
# Fit the glm
m1.ph=(glm( ph ~ Age, family=tweedie(var.power=out.ph$xi.max, link.power=0) ))
m2.ph<- glm(ph~Age)
m3.ph=(glm( ph ~ Age+Cattle1617+treesedge, family=tweedie(var.power=out.ph$xi.max, link.power=0) ))
m4.ph=(glm( ph ~ Age+Cattle1617+treesedge))

AICtweedie(m1.ph) #136.6735
AIC(m2.ph) #137.8237
AICtweedie(m3.ph) #135.2255
AIC(m4.ph) #136.5233

summary(m4.ph) #effects for Cattle and Edge Tree cover
#             Estimate Std. Error t value Pr(>|t|)  
#Cattle1617   0.08839    0.29483   0.300 0.765741    
#treesedge   -0.66842    0.38307  -1.745 0.087987 .

tukph=(glht(m4.ph,mcp(Age="Tukey"))) #using normal dist model (m2) since normal and AIC not diff
summary(tukph)
tuk.cld.ph<- cld(tukph)
plot(tuk.cld.ph) #derive letter assignments from plot for significance of differences between categories (alpha = 0.05)
#<30 30-39 40-49 50-59 60+
#b   ab    b     a     a

#Multiple Comparisons of Means: Tukey Contrasts
#                   Estimate Std. Error z value Pr(>|z|)    
#30-39 - <30 == 0    -0.3237     0.4317  -0.750  0.94381    
#40-49 - <30 == 0     0.0377     0.3826   0.099  0.99998    
#50-59 - <30 == 0    -1.1478     0.3348  -3.428  0.00531 ** 
#60+ - <30 == 0      -1.4837     0.3825  -3.879  < 0.001 ***
#40-49 - 30-39 == 0   0.3614     0.4475   0.807  0.92741    
#50-59 - 30-39 == 0  -0.8241     0.4069  -2.025  0.25076    
#60+ - 30-39 == 0    -1.1601     0.4475  -2.592  0.07061 .  
#50-59 - 40-49 == 0  -1.1855     0.3542  -3.347  0.00722 ** 
#60+ - 40-49 == 0    -1.5214     0.4003  -3.801  0.00131 ** 
#60+ - 50-59 == 0    -0.3360     0.3543  -0.948  0.87632    

#GLM- Slope (rise:run)####
out.sl <- tweedie.profile( sl1617~1, xi.vec=xi.vec, do.plot=TRUE, verbose=TRUE)
# Fit the glm
m1.sl=glm(sl1617 ~ Age, family=tweedie(var.power=out.sl$xi.max, link.power=0))
m2.sl=glm(sl1617 ~ Age)
m3.sl=glm(sl1617 ~ Age+Cattle1617+treesedge, family=tweedie(var.power=out.sl$xi.max, link.power=0))
m4.sl=glm(sl1617 ~ Age+Cattle1617+treesedge)

AICtweedie(m1.sl) #-154.0144
AIC(m2.sl) #-154.0668
AICtweedie(m3.sl) #-166.7661
AIC(m4.sl) #-167.6125

#AIC values suggest tweedie unnecessary, response variable sufficiently normal
summary(m4.sl) #effects for Cattle and Edge Tree cover
#             Estimate Std. Error t value Pr(>|t|)  
#Cattle1617  -0.033538   0.014949  -2.243 0.029950 *  
#treesedge    0.039565   0.019424   2.037 0.047701 * 

tuksl=(glht(m4.sl,mcp(Age="Tukey"))) 
summary(tuksl)
tuk.sl.cld<- cld(tuksl)
par(mar=c(5,5,5,5))
plot(tuk.sl.cld) #derive letter assignments from plot for significance of differences between categories (alpha = 0.05)
#<30 30-39 40-49 50-59 60+
#b   b     b     ab    a


#Multiple Comparisons of Means: Tukey Contrasts
#                    Estimate Std. Error z value Pr(>|z|)   
#30-39 - <30 == 0    0.004220   0.021890   0.193  0.99969   
#40-49 - <30 == 0    0.002317   0.019398   0.119  0.99995   
#50-59 - <30 == 0   -0.027015   0.016976  -1.591  0.49921   
#60+ - <30 == 0     -0.069460   0.019393  -3.582  0.00307 **
#40-49 - 30-39 == 0 -0.001903   0.022693  -0.084  0.99999   
#50-59 - 30-39 == 0 -0.031235   0.020631  -1.514  0.54966   
#60+ - 30-39 == 0   -0.073680   0.022693  -3.247  0.01009 * 
#50-59 - 40-49 == 0 -0.029332   0.017960  -1.633  0.47227   
#60+ - 40-49 == 0   -0.071778   0.020297  -3.536  0.00366 **
#60+ - 50-59 == 0   -0.042446   0.017965  -2.363  0.12378   


#Plot- Emergent Vegetation Cover (%) -----------------------------------------------------

Emsummary=group_by(data, AgeCategory) %>%
  summarise(
    count = n(),
    mean = mean(EM1617, na.rm = TRUE),
    sd = sd(EM1617, na.rm = TRUE)
  )
(Emplot=ggplot(Emsummary,aes(AgeCategory,mean))+
    geom_line(aes(group=1))+
    geom_errorbar(aes(ymin = mean - (sd/sqrt(count)), ymax = mean + (sd/sqrt(count))), width=0.2)+
    geom_point(data=data, aes(AgeCategory,EM1617, size="1"),position=position_jitterdodge(jitter.width = 0.4, jitter.height = 0,
                                                                                                           dodge.width = 0.4))+
    geom_line(aes(x=AgeCategory,y=mean)))
(Emnice= Emplot + labs(y="Emergent Vegetation Cover (%)", x = "Pond Age (years)") + theme_pubr()+geom_point()+ theme(legend.position="none"))



#Plot- Floating Vegetation Cover (%)--------------------------------
FLsummary=group_by(data, AgeCategory) %>%
  summarise(
    count = n(),
    mean = mean(FL1617, na.rm = TRUE),
    sd = sd(FL1617, na.rm = TRUE)
  )
(FLplot=ggplot(FLsummary,aes(AgeCategory,mean))+
    geom_line(aes(group=1))+
    geom_errorbar(aes(ymin = mean - (sd/sqrt(count)), ymax = mean + (sd/sqrt(count))), width=0.2)+
    geom_point(data=data, aes(AgeCategory,FL1617, size="1"),position=position_jitterdodge(jitter.width = 0.4, jitter.height = 0,dodge.width = 0.4))+
    geom_line(aes(x=AgeCategory,y=mean)))
(FLnice=  FLplot + labs(y="Floating Vegetation Cover (%)", x = "Pond Age (years)") + theme_pubr()+geom_point()+ theme(legend.position="none"))

#Plot- Submerged Vegetation Cover (%)------------------------------------------------

Subsummary=group_by(data, AgeCategory) %>%
  summarise(
    count = n(),
    mean = mean(AQ1617, na.rm = TRUE),
    sd = sd(AQ1617, na.rm = TRUE)
  )
(Subplot=ggplot(Subsummary,aes(AgeCategory,mean))+
    geom_line(aes(group=1))+
    geom_errorbar(aes(ymin = mean - (sd/sqrt(count)), ymax = mean + (sd/sqrt(count))), width=0.2)+
    geom_point(data=data, aes(AgeCategory,AQ1617, size="1"),position=position_jitterdodge(jitter.width = 0.4, jitter.height = 0,
                                                                                                            dodge.width = 0.4))+
    geom_line(aes(x=AgeCategory,y=mean)))
(Subnice= Subplot + labs(y="Submerged Vegetation Cover (%)", x = "Pond Age (years)") + theme_pubr() +geom_point()+ theme(legend.position="none"))


#Plot- pH ------------------------------------------------

pHsummary=group_by(data, AgeCategory) %>%
  summarise(
    count = n(),
    mean = mean(pH, na.rm = TRUE),
    sd = sd(pH, na.rm = TRUE)
  )
(pHplot=ggplot(pHsummary,aes(AgeCategory,mean))+
    geom_line(aes(group=1))+
    geom_errorbar(aes(ymin = mean - (sd/sqrt(count)), ymax = mean + (sd/sqrt(count))), width=0.2)+
    geom_point(data=data, aes(AgeCategory,pH, size="1"),position=position_jitterdodge(jitter.width = 0.4, jitter.height = 0,
                                                                                                        dodge.width = 0.4))+
    geom_line(aes(x=AgeCategory,y=mean)))
(pHnice= pHplot + labs(y="pH", x = "Pond Age (years)") + theme_pubr() +geom_point() + theme(legend.position="none"))


#Plot- Slope (rise:run)--------------------------------------------

SLsummary=group_by(data, AgeCategory) %>%
  summarise(
    count = n(),
    mean = mean(SL1617, na.rm = TRUE),
    sd = sd(SL1617, na.rm = TRUE)
  )
(SLplot=ggplot(SLsummary,aes(AgeCategory,mean))+
    geom_line(aes(group=1))+
    geom_errorbar(aes(ymin = mean - (sd/sqrt(count)), ymax = mean + (sd/sqrt(count))), width=0.2)+
    geom_point(data=data, aes(AgeCategory,SL1617, size="1"),position=position_jitterdodge(jitter.width = 0.4, jitter.height = 0,
                                                                                                            dodge.width = 0.4))+
    geom_line(aes(x=AgeCategory,y=mean)))
(SLnice= SLplot + labs(y="Slope", x = "Pond Age (years)") + theme_pubr() +geom_point()+theme(legend.position="none"))


# Generate final multi-paneled plot of all 5 variables -------------------
grid.arrange(Emnice, FLnice, Subnice, pHnice,SLnice, nrow = 3)
