#title: October 2014 Experiment R Analysis Script
#author: Sheila Saia
#date: January 12, 2017

# ----
## Introduction
# ----

#This is the R code used to run various analyses for the paper titled 'Evidence for polyphosphate accumulating organism (PAO)-mediated phosphorus cycling in stream biofilms under alternating aerobic/anaerobic conditions' by Saia et al. (2017) in Freshwater Science. See the readme.txt file in the associated GitHub repository for additional information and metadata about this project (https://github.com/sheilasaia/paper-p-cycling-in-stream-biofilms).

#If you have questions, comments, or would like to report an error in the code, please contact Sheila Saia at sms493 at cornell dot edu.

# ----
## Load Packages & Set File Directory
# ----

# load packages
library(Cairo) #for plotting
library(mgcv) #for gam modeling

# set directory (edit this for personal use)
setwd("C:\\Users\\Sheila\\Documents\\Cornell\\Research\\PhD!\\biofilms\\Oct2014_experiment\\forPaper\\forGitHub")

# ----
## Import Data
# ----

# import data
data=read.table("allTUBdata_oct2014_forPaper.txt",header=T,sep="\t")
# T1 = aerobic/anaerobic
# T2 = aerobic

# ----
## Reformat Data
# ----

# wet biofilm mass (g) per treatment
T1BFmass=8.949
T2BFmass=10.223

# rock SA (m^2) per treatment
T1BFsa=0.554+0.545
T2BFsa=0.406+0.544

# total SA and mass per treatment
rockdata=data.frame(cbind(rbind(T1BFmass,T2BFmass),rbind(T1BFsa,T2BFsa)))
colnames(rockdata)=c("wetBFg","SAm2")
rownames(rockdata)=c("T1","T2")

# wet biofilm mass per m2 for each treatment
rockdata$gBFperm2=c(T1BFmass/T1BFsa,T2BFmass/T2BFsa)

# calc adjusted AvgSRPppm based on mass of biofilm in tub
data$AvgSRPppmAdj[data$TubID=='T1']=data$AvgSRPppm[data$TubID=='T1']/T1BFmass
data$AvgSRPppmAdj[data$TubID=='T2']=data$AvgSRPppm[data$TubID=='T2']/T2BFmass

# calc adjusted StdSRPppm based on mass of biofilm in tub
data$StdSRPppmAdj[data$TubID=='T1']=data$StdSRPppm[data$TubID=='T1']/T1BFmass
data$StdSRPppmAdj[data$TubID=='T2']=data$StdSRPppm[data$TubID=='T2']/T2BFmass

# calc adjusted TotalPppm based on mass of biofilm in tub
data$TotalPppmAdj[data$TubID=='T1']=data$TotalPppm[data$TubID=='T1']/T1BFmass
data$TotalPppmAdj[data$TubID=='T2']=data$TotalPppm[data$TubID=='T2']/T2BFmass

# calc adjusted Fe2ppm based on mass of biofilm in tub
data$Fe2ppmAdj[data$TubID=='T1']=data$Fe2ppm[data$TubID=='T1']/T1BFmass
data$Fe2ppmAdj[data$TubID=='T2']=data$Fe2ppm[data$TubID=='T2']/T2BFmass

# calc adjusted TotalFeppm based on mass of biofilm in tub
data$TotalFeppmAdj[data$TubID=='T1']=data$TotalFeppm[data$TubID=='T1']/T1BFmass
data$TotalFeppmAdj[data$TubID=='T2']=data$TotalFeppm[data$TubID=='T2']/T2BFmass

# calc adjusted TotalCappm based on mass of biofilm in tub
data$TotalCappmAdj[data$TubID=='T1']=data$TotalCappm[data$TubID=='T1']/T1BFmass
data$TotalCappmAdj[data$TubID=='T2']=data$TotalCappm[data$TubID=='T2']/T2BFmass

# calc adjusted TotalSppm based on mass of biofilm in tub
data$TotalSppmAdj[data$TubID=='T1']=data$TotalSppm[data$TubID=='T1']/T1BFmass
data$TotalSppmAdj[data$TubID=='T2']=data$TotalSppm[data$TubID=='T2']/T2BFmass

# calc adjusted TotalKppm based on mass of biofilm in tub
data$TotalKppmAdj[data$TubID=='T1']=data$TotalKppm[data$TubID=='T1']/T1BFmass
data$TotalKppmAdj[data$TubID=='T2']=data$TotalKppm[data$TubID=='T2']/T2BFmass

# calc adjusted TotalMgppm based on mass of biofilm in tub
data$TotalMgppmAdj[data$TubID=='T1']=data$TotalMgppm[data$TubID=='T1']/T1BFmass
data$TotalMgppmAdj[data$TubID=='T2']=data$TotalMgppm[data$TubID=='T2']/T2BFmass

# calc adjusted TotalMnppm based on mass of biofilm in tub
data$TotalMnppmAdj[data$TubID=='T1']=data$TotalMnppm[data$TubID=='T1']/T1BFmass
data$TotalMnppmAdj[data$TubID=='T2']=data$TotalMnppm[data$TubID=='T2']/T2BFmass

# rearrange data
T1data=data[data$TubID=='T1',]
T2data=data[data$TubID=='T2',]

# ----
## Phosphate GAM (Phosphate vs Hours from Start)
# ----

# gam model of both treatments together
my.gamp=gam(AvgSRPppmAdj~s(HoursFromStart,by=as.numeric(TubID=="T1"))+s(HoursFromStart,by=as.numeric(TubID=="T2")),data=data,method="REML")
summary(my.gamp)
#plot(resid(my.gamp))
#qqnorm(resid(my.gamp))
#qqline(resid(my.gamp)) #looks normal

# with pH
my.gamp2=gam(AvgSRPppmAdj~pH+s(HoursFromStart,by=as.numeric(TubID=="T1"))+s(HoursFromStart,by=as.numeric(TubID=="T2")),data=data,method="REML")
#summary(my.gamp2)

# pH with smooth
my.gamp3=gam(AvgSRPppmAdj~s(pH)+s(HoursFromStart,by=as.numeric(TubID=="T1"))+s(HoursFromStart,by=as.numeric(TubID=="T2")),data=data,method="REML")
#summary(my.gamp3)
#plot(my.gamp3)

# pH with smooth 2
my.gamp4=gam(AvgSRPppmAdj~s(pH,by=as.numeric(TubID=="T1"))+s(pH,by=as.numeric(TubID=="T2"))+s(HoursFromStart,by=as.numeric(TubID=="T1"))+s(HoursFromStart,by=as.numeric(TubID=="T2")),data=data,method="REML")
#summary(my.gamp4)
#plot(my.gamp4)

# lowest AIC is best fit 
# must be 2 AIC units diferent to be statistically different
AIC(my.gamp,my.gamp2,my.gamp3,my.gamp4)

# predictions using my.gamp
newdatagamT1p=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T1",length(seq(0,47,0.5))))
newdatagamT2p=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T2",length(seq(0,47,0.5))))
gam.predict.T1p=predict(my.gamp,newdata=newdatagamT1p,se.fit=T)
gam.predict.T2p=predict(my.gamp,newdata=newdatagamT2p,se.fit=T)

# t critical calculation
alpha=0.05
gamp.n=length(newdatagamT1p$HoursFromStart)
gamp.p=length(coef(my.gamp))
tcrit95.gamp=qt(1-(alpha/2),gamp.n-gamp.p)

# calculate upper and lower bound y's using se's model fit
T1gamp.UP=gam.predict.T1p$fit+tcrit95.gamp*gam.predict.T1p$se.fit
T1gamp.LO=gam.predict.T1p$fit-tcrit95.gamp*gam.predict.T1p$se.fit
T2gamp.UP=gam.predict.T2p$fit+tcrit95.gamp*gam.predict.T2p$se.fit
T2gamp.LO=gam.predict.T2p$fit-tcrit95.gamp*gam.predict.T2p$se.fit

# ----
## Ca GAM (Ca vs Hours from Start)
# ----

# gam model of both treatments together
my.gamca=gam(TotalCappmAdj~s(HoursFromStart,by=as.numeric(TubID=="T1"))+s(HoursFromStart,by=as.numeric(TubID=="T2")),data=data,method="REML")
summary(my.gamca)
#plot(my.gamca)
#qqnorm(resid(my.gamca))
#qqline(resid(my.gamca)) #looks normal

# predictions
newdataT1ca=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T1",length(seq(0,47,0.5))))
newdataT2ca=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T2",length(seq(0,47,0.5))))
gam.predict.T1ca=predict(my.gamca,newdata=newdataT1ca,se.fit=T)
gam.predict.T2ca=predict(my.gamca,newdata=newdataT2ca,se.fit=T)

# t critical calculation
alpha=0.05
gamca.n=length(newdataT1ca$HoursFromStart)
gamca.p=length(coef(my.gamca))
tcrit95.gamca=qt(1-(alpha/2),gamca.n-gamca.p)

# calculate upper and lower bound y's using se's model fit
T1gamca.UP=gam.predict.T1ca$fit+tcrit95.gamca*gam.predict.T1ca$se.fit
T1gamca.LO=gam.predict.T1ca$fit-tcrit95.gamca*gam.predict.T1ca$se.fit
T2gamca.UP=gam.predict.T2ca$fit+tcrit95.gamca*gam.predict.T2ca$se.fit
T2gamca.LO=gam.predict.T2ca$fit-tcrit95.gamca*gam.predict.T2ca$se.fit

# ----
## K GAM (K vs Hours from Start)
# ----

# gam model of both treatments together
my.gamk=gam(TotalKppmAdj~s(HoursFromStart,by=as.numeric(TubID=="T1"))+s(HoursFromStart,by=as.numeric(TubID=="T2")),data=data,method="REML")
summary(my.gamk)
#plot(my.gamk)
#qqnorm(resid(my.gamk))
#qqline(resid(my.gamk)) #looks normal

# predictions
newdataT1k=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T1",length(seq(0,47,0.5))))
newdataT2k=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T2",length(seq(0,47,0.5))))
gam.predict.T1k=predict(my.gamk,newdata=newdataT1k,se.fit=T)
gam.predict.T2k=predict(my.gamk,newdata=newdataT2k,se.fit=T)

# t critical calculation
alpha=0.05
gamk.n=length(newdataT1k$HoursFromStart)
gamk.p=length(coef(my.gamk))
tcrit95.gamk=qt(1-(alpha/2),gamk.n-gamk.p)

# calculate upper and lower bound y's using se's model fit
T1gamk.UP=gam.predict.T1k$fit+tcrit95.gamk*gam.predict.T1k$se.fit
T1gamk.LO=gam.predict.T1k$fit-tcrit95.gamk*gam.predict.T1k$se.fit
T2gamk.UP=gam.predict.T2k$fit+tcrit95.gamk*gam.predict.T2k$se.fit
T2gamk.LO=gam.predict.T2k$fit-tcrit95.gamk*gam.predict.T2k$se.fit

# ----
## Mg GAM (Mg vs Hours from Start)
# ----

# gam model of both treatments together
my.gammg=gam(TotalMgppmAdj~s(HoursFromStart,by=as.numeric(TubID=="T1"))+s(HoursFromStart,by=as.numeric(TubID=="T2")),data=data,method="REML")
summary(my.gammg)
#plot(my.gammg)
#qqnorm(resid(my.gammg))
#qqline(resid(my.gammg)) #looks normal

# predictions
newdataT1mg=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T1",length(seq(0,47,0.5))))
newdataT2mg=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T2",length(seq(0,47,0.5))))
gam.predict.T1mg=predict(my.gammg,newdata=newdataT1mg,se.fit=T)
gam.predict.T2mg=predict(my.gammg,newdata=newdataT2mg,se.fit=T)

# t critical calculation
alpha=0.05
gammg.n=length(newdataT1mg$HoursFromStart)
gammg.p=length(coef(my.gammg))
tcrit95.gammg=qt(1-(alpha/2),gammg.n-gammg.p)

# calculate upper and lower bound y's using se's model fit
T1gammg.UP=gam.predict.T1mg$fit+tcrit95.gammg*gam.predict.T1mg$se.fit
T1gammg.LO=gam.predict.T1mg$fit-tcrit95.gammg*gam.predict.T1mg$se.fit
T2gammg.UP=gam.predict.T2mg$fit+tcrit95.gammg*gam.predict.T2mg$se.fit
T2gammg.LO=gam.predict.T2mg$fit-tcrit95.gammg*gam.predict.T2mg$se.fit

# ----
## Mn GAM (Mn vs Hours from Start)
# ----

# gam model of both treatments together
my.gammn=gam(TotalMnppmAdj~s(HoursFromStart,by=as.numeric(TubID=="T1"))+s(HoursFromStart,by=as.numeric(TubID=="T2")),data=data,method="REML")
summary(my.gammn)
#plot(resid(my.gammn))
#qqnorm(resid(my.gammn))
#qqline(resid(my.gammn)) #looks normal

# predictions
newdatagamT1mn=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T1",length(seq(0,47,0.5))))
newdatagamT2mn=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T2",length(seq(0,47,0.5))))
gam.predict.T1mn=predict(my.gammn,newdata=newdatagamT1mn,se.fit=T)
gam.predict.T2mn=predict(my.gammn,newdata=newdatagamT2mn,se.fit=T)

# t critical calculation
alpha=0.05
gammn.n=length(newdatagamT1mn$HoursFromStart)
gammn.p=length(coef(my.gammn))
tcrit95.gammn=qt(1-(alpha/2),gammn.n-gammn.p)

# calculate upper and lower bound y's using se's model fit
T1gammn.UP=gam.predict.T1mn$fit+tcrit95.gammn*gam.predict.T1mn$se.fit
T1gammn.LO=gam.predict.T1mn$fit-tcrit95.gammn*gam.predict.T1mn$se.fit
T2gammn.UP=gam.predict.T2mn$fit+tcrit95.gammn*gam.predict.T2mn$se.fit
T2gammn.LO=gam.predict.T2mn$fit-tcrit95.gammn*gam.predict.T2mn$se.fit

# ----
## FeII GAM (FeII vs Hours from Start)
# ----

# gam model of both treatments together
noNAFe2data=data[which(data$Fe2ppmAdj!="NA"),] #take out NA's but gives the same result if you just use 'data'
my.gamfe=gam(Fe2ppmAdj~s(HoursFromStart,by=as.numeric(TubID=="T1"))+s(HoursFromStart,by=as.numeric(TubID=="T2")),data=data,method="REML")
summary(my.gamfe)
#plot(my.gamfe)
#qqnorm(resid(my.gamfe))
#qqline(resid(my.gamfe)) #looks normal

# gam model with no trend (for comparison of T2)
my.lmfenotrend=lm(Fe2ppmAdj~1,data=noNAFe2data)
summary(my.lmfenotrend)
#
# compare AICs
AIC(my.gamfe,my.lmfenotrend)

# predictions
newdataT1fe=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T1",length(seq(0,47,0.5))))
newdataT2fe=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T2",length(seq(0,47,0.5))))
gam.predict.T1fe=predict(my.gamfe,newdata=newdataT1fe,se.fit=T)
gam.predict.T2fe=predict(my.gamfe,newdata=newdataT2fe,se.fit=T)

# t critical calculation
alpha=0.05
gamfe.n=length(newdataT1fe$HoursFromStart)
gamfe.p=length(coef(my.gamfe))
tcrit95.gamfe=qt(1-(alpha/2),gamfe.n-gamfe.p)

# calculate upper and lower bound y's using se's model fit
T1gamfe.UP=gam.predict.T1fe$fit+tcrit95.gamfe*gam.predict.T1fe$se.fit
T1gamfe.LO=gam.predict.T1fe$fit-tcrit95.gamfe*gam.predict.T1fe$se.fit
T2gamfe.UP=gam.predict.T2fe$fit+tcrit95.gamfe*gam.predict.T2fe$se.fit
T2gamfe.LO=gam.predict.T2fe$fit-tcrit95.gamfe*gam.predict.T2fe$se.fit

# ----
## Total S GAM (Total S vs Hours from Start)
# ----

# gam model of both treatments together
my.gams=gam(TotalSppmAdj~s(HoursFromStart,by=as.numeric(TubID=="T1"))+s(HoursFromStart,by=as.numeric(TubID=="T2")),data=data,method="REML")
summary(my.gams)
#plot(my.gams)
#qqnorm(resid(my.gams))
#qqline(resid(my.gams)) #looks normal

# predictions
newdataT1s=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T1",length(seq(0,47,0.5))))
newdataT2s=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T2",length(seq(0,47,0.5))))
gam.predict.T1s=predict(my.gams,newdata=newdataT1s,se.fit=T)
gam.predict.T2s=predict(my.gams,newdata=newdataT2s,se.fit=T)

# t critical calculation
alpha=0.05
gams.n=length(newdataT1s$HoursFromStart)
gams.p=length(coef(my.gams))
tcrit95.gams=qt(1-(alpha/2),gams.n-gams.p)

# calculate upper and lower bound y's using se's model fit
T1gams.UP=gam.predict.T1s$fit+tcrit95.gams*gam.predict.T1s$se.fit
T1gams.LO=gam.predict.T1s$fit-tcrit95.gams*gam.predict.T1s$se.fit
T2gams.UP=gam.predict.T2s$fit+tcrit95.gams*gam.predict.T2s$se.fit
T2gams.LO=gam.predict.T2s$fit-tcrit95.gams*gam.predict.T2s$se.fit

# ----
## Figure of Phophate vs Hours from Start (Normalized)
# ----

#CairoPDF("JustPvsTime.pdf", width= 10, height=10,pointsize=14)
plot(data$AvgSRPppmAdj[data$TubID=="T1"]~data$HoursFromStart[data$TubID=="T1"],pch=16,ylim=c(0,0.01),xlab="Hours from start",ylab="Phosphate-P (mg/L/g wet biofilm)")
lines(newdatagamT1p$HoursFromStart,gam.predict.T1p$fit,lwd=3,col="black")
lines(newdatagamT1p$HoursFromStart,T1gamp.UP,lwd=0.5,col="black")
lines(newdatagamT1p$HoursFromStart,T1gamp.LO,lwd=0.5,col="black")
points(data$AvgSRPppmAdj[data$TubID=="T2"]~data$HoursFromStart[data$TubID=="T2"],pch=1)
lines(newdatagamT2p$HoursFromStart,gam.predict.T2p$fit,lwd=3,lty=2,col="black")
lines(newdatagamT2p$HoursFromStart,T2gamp.UP,lwd=0.5,col="black")
lines(newdatagamT2p$HoursFromStart,T2gamp.LO,lwd=0.5,col="black")
rect(11.083,0,22.583,0.002/10,col='black')
rect(35.083,0,46.583,0.002/10,col='black')
legend("topright",c("T1 (aerobic/anaerobic)","T2 (aerobic)","T1 GAM","T2 GAM","95% confidence intervals"),pch=c(16,1,NA,NA,NA),lwd=c(NA,NA,3,3,1.5),lty=c(NA,NA,1,2,1))
#dev.off() # writes file to working directory

# ----
## Figure of All Cations vs Hours from Start (Normalized)
# ----

#CairoPDF("AllCationsvsTime.pdf", width= 11, height=14,pointsize=14)
par(mfrow=c(3,2))
# Fe2+ figure
plot(data$Fe2ppmAdj[data$TubID=="T1"]~data$HoursFromStart[data$TubID=="T1"],pch=16,ylim=c(-0.3,0.5),xlab="",ylab="FeII (mg/L/g wet biofilm)")
lines(newdataT1fe$HoursFromStart,gam.predict.T1fe$fit,lwd=3,col="black")
lines(newdataT1fe$HoursFromStart,T1gamfe.UP,lwd=0.5,col="black")
lines(newdataT1fe$HoursFromStart,T1gamfe.LO,lwd=0.5,col="black")
points(data$Fe2ppmAdj[data$TubID=="T2"]~data$HoursFromStart[data$TubID=="T2"],pch=1)
lines(newdataT2fe$HoursFromStart,gam.predict.T2fe$fit,lwd=3,lty=2,col="black")
lines(newdataT2fe$HoursFromStart,T2gamfe.UP,lwd=0.5,col="black")
lines(newdataT2fe$HoursFromStart,T2gamfe.LO,lwd=0.5,col="black")
rect(11.083,-0.3,22.583,-0.3+0.2/10,col='black')
rect(35.083,-0.3,46.583,-0.3+0.2/10,col='black')
legend("topleft",c("T1 (aerobic/anaerobic)","T2 (aerobic)","T1 GAM","T2 GAM","95% confidence intervals"),pch=c(16,1,NA,NA,NA),lwd=c(NA,NA,3,3,1.5),lty=c(NA,NA,1,2,1))

# TS figure
plot(data$TotalSppmAdj[data$TubID=="T1"]~data$HoursFromStart[data$TubID=="T1"],pch=16,ylim=c(0,40),xlab="",ylab="Total S (mg/L/g wet biofilm)")
lines(newdataT1s$HoursFromStart,gam.predict.T1s$fit,lwd=3,col="black")
lines(newdataT1s$HoursFromStart,T1gams.UP,lwd=0.5,col="black")
lines(newdataT1s$HoursFromStart,T1gams.LO,lwd=0.5,col="black")
points(data$TotalSppmAdj[data$TubID=="T2"]~data$HoursFromStart[data$TubID=="T2"],pch=1)
lines(newdataT2s$HoursFromStart,gam.predict.T2s$fit,lwd=3,lty=2,col="black")
lines(newdataT2s$HoursFromStart,T2gams.UP,lwd=0.5,col="black")
lines(newdataT2s$HoursFromStart,T2gams.LO,lwd=0.5,col="black")
rect(11.083,0,22.583,10/10,col='black')
rect(35.083,0,46.583,10/10,col='black')

# Mn figure
plot(data$TotalMnppmAdj[data$TubID=="T1"]~data$HoursFromStart[data$TubID=="T1"],pch=16,ylim=c(-0.02,0.18),xlab="",ylab="Total Mn (mg/L/g wet biofilm)")
lines(newdatagamT1mn$HoursFromStart,gam.predict.T1mn$fit,lwd=3,col="black")
lines(newdatagamT1mn$HoursFromStart,T1gammn.UP,lwd=0.5,col="black")
lines(newdatagamT1mn$HoursFromStart,T1gammn.LO,lwd=0.5,col="black")
points(data$TotalMnppmAdj[data$TubID=="T2"]~data$HoursFromStart[data$TubID=="T2"],pch=1)
lines(newdatagamT2mn$HoursFromStart,gam.predict.T2mn$fit,lwd=3,lty=2,col="black")
lines(newdatagamT2mn$HoursFromStart,T2gammn.UP,lwd=0.5,col="black")
lines(newdatagamT2mn$HoursFromStart,T2gammn.LO,lwd=0.5,col="black")
rect(11.083,-0.02,22.583,-0.018,col='black')
rect(35.083,-0.02,46.583,-0.018,col='black')

# Ca figure
plot(data$TotalCappmAdj[data$TubID=="T1"]~data$HoursFromStart[data$TubID=="T1"],pch=16,ylim=c(0,13),xlab="",ylab="Ca (mg/L/g wet biofilm)")
lines(newdataT1ca$HoursFromStart,gam.predict.T1ca$fit,lwd=3,col="black")
lines(newdataT1ca$HoursFromStart,T1gamca.UP,lwd=0.5,col="black")
lines(newdataT1ca$HoursFromStart,T1gamca.LO,lwd=0.5,col="black")
points(data$TotalCappmAdj[data$TubID=="T2"]~data$HoursFromStart[data$TubID=="T2"],pch=1)
lines(newdataT2ca$HoursFromStart,gam.predict.T2ca$fit,lwd=3,lty=2,col="black")
lines(newdataT2ca$HoursFromStart,T2gamca.UP,lwd=0.5,col="black")
lines(newdataT2ca$HoursFromStart,T2gamca.LO,lwd=0.5,col="black")
rect(11.083,0,22.583,2/10,col='black')
rect(35.083,0,46.583,2/10,col='black')

# K figure
plot(data$TotalKppmAdj[data$TubID=="T1"]~data$HoursFromStart[data$TubID=="T1"],pch=16,ylim=c(0.2,0.7),xlab="Hours from start",ylab="K (mg/L/g wet biofilm)")
lines(newdataT1k$HoursFromStart,gam.predict.T1k$fit,lwd=3,col="black")
lines(newdataT1k$HoursFromStart,T1gamk.UP,lwd=0.5,col="black")
lines(newdataT1k$HoursFromStart,T1gamk.LO,lwd=0.5,col="black")
points(data$TotalKppmAdj[data$TubID=="T2"]~data$HoursFromStart[data$TubID=="T2"],pch=1)
lines(newdataT2k$HoursFromStart,gam.predict.T2k$fit,lwd=3,lty=2,col="black")
lines(newdataT2k$HoursFromStart,T2gamk.UP,lwd=0.5,col="black")
lines(newdataT2k$HoursFromStart,T2gamk.LO,lwd=0.5,col="black")
rect(11.083,0.2,22.583,0.2+0.1/10,col='black')
rect(35.083,0.2,46.583,0.2+0.1/10,col='black')

# Mg figure
plot(data$TotalMgppmAdj[data$TubID=="T1"]~data$HoursFromStart[data$TubID=="T1"],pch=16,ylim=c(0.6,1.8),xlab="Hours from start",ylab="Mg (mg/L/g wet biofilm)")
lines(newdataT1mg$HoursFromStart,gam.predict.T1mg$fit,lwd=3,col="black")
lines(newdataT1mg$HoursFromStart,T1gammg.UP,lwd=0.5,col="black")
lines(newdataT1mg$HoursFromStart,T1gammg.LO,lwd=0.5,col="black")
points(data$TotalMgppmAdj[data$TubID=="T2"]~data$HoursFromStart[data$TubID=="T2"],pch=1)
lines(newdataT2mg$HoursFromStart,gam.predict.T2mg$fit,lwd=3,lty=2,col="black")
lines(newdataT2mg$HoursFromStart,T2gammg.UP,lwd=0.5,col="black")
lines(newdataT2mg$HoursFromStart,T2gammg.LO,lwd=0.5,col="black")
rect(11.083,0.6,22.583,0.6+0.2/10,col='black')
rect(35.083,0.6,46.583,0.6+0.2/10,col='black')
#dev.off() # writes file to working directory

# ----
## GAM and Figure of Phophate vs Hours from Start (Not Normalized, For Suppliment)
# ----

# P gam model of both treatments together
my.gampNN=gam(AvgSRPppm~s(HoursFromStart,by=as.numeric(TubID=="T1"))+s(HoursFromStart,by=as.numeric(TubID=="T2")),data=data,method="REML")
#summary(my.gampNN)
#plot(resid(my.gampNN))
#qqnorm(resid(my.gampNN))
#qqline(resid(my.gampNN)) #looks ok

# predictions
newdatagamT1pNN=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T1",length(seq(0,47,0.5))))
newdatagamT2pNN=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T2",length(seq(0,47,0.5))))
gam.predict.T1pNN=predict(my.gampNN,newdata=newdatagamT1pNN,se.fit=T)
gam.predict.T2pNN=predict(my.gampNN,newdata=newdatagamT2pNN,se.fit=T)

# t critical calculation
alpha=0.05
gampNN.n=length(newdatagamT1pNN$HoursFromStart)
gampNN.p=length(coef(my.gampNN))
tcrit95.gampNN=qt(1-(alpha/2),gampNN.n-gampNN.p)

# calculate upper and lower bound y's using se's model fit
T1gampNN.UP=gam.predict.T1pNN$fit+tcrit95.gampNN*gam.predict.T1pNN$se.fit
T1gampNN.LO=gam.predict.T1pNN$fit-tcrit95.gampNN*gam.predict.T1pNN$se.fit
T2gampNN.UP=gam.predict.T2pNN$fit+tcrit95.gampNN*gam.predict.T2pNN$se.fit
T2gampNN.LO=gam.predict.T2pNN$fit-tcrit95.gampNN*gam.predict.T2pNN$se.fit

# P figure
#CairoPDF("JustPNotNormalizedvsTime.pdf", width= 10, height=10,pointsize=14)
par(mfrow=c(1,1))
plot(data$AvgSRPppm[data$TubID=="T1"]~data$HoursFromStart[data$TubID=="T1"],pch=16,ylim=c(0,0.1),xlab="Hours from start",ylab="Phosphate-P (mg/L)")
lines(newdatagamT1pNN$HoursFromStart,gam.predict.T1pNN$fit,lwd=3,col="black")
lines(newdatagamT1pNN$HoursFromStart,T1gampNN.UP,lwd=0.5,col="black")
lines(newdatagamT1pNN$HoursFromStart,T1gampNN.LO,lwd=0.5,col="black")
points(data$AvgSRPppm[data$TubID=="T2"]~data$HoursFromStart[data$TubID=="T2"],pch=1)
lines(newdatagamT2pNN$HoursFromStart,gam.predict.T2pNN$fit,lwd=3,lty=2,col="black")
lines(newdatagamT2pNN$HoursFromStart,T2gampNN.UP,lwd=0.5,col="black")
lines(newdatagamT2pNN$HoursFromStart,T2gampNN.LO,lwd=0.5,col="black")
rect(11.083,0,22.583,0.002,col='black')
rect(35.083,0,46.583,0.002,col='black')
legend("topright",c("T1 (aerobic/anaerobic)","T2 (aerobic)","T1 GAM","T2 GAM","95% confidence intervals"),pch=c(16,1,NA,NA,NA),lwd=c(NA,NA,3,3,1.5),lty=c(NA,NA,1,2,1))
#dev.off() # writes file to working directory

# ----
## GAMs and Figures of All Cations vs Hours from Start (Not Normalized - NN, For Suppliment)
# ----
# FeII gam model of both treatments together
noNAFe2dataNN=data[which(data$Fe2ppm!="NA"),] #take out NA's but gives the same result if you just use 'data'
my.gamfeNN=gam(Fe2ppm~s(HoursFromStart,by=as.numeric(TubID=="T1"))+s(HoursFromStart,by=as.numeric(TubID=="T2")),data=noNAFe2dataNN,method="REML")
#summary(my.gamfeNN)
#qqnorm(resid(my.gamfeNN))
#qqline(resid(my.gamfeNN)) #looks ok

# predictions
newdataT1feNN=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T1",length(seq(0,47,0.5))))
newdataT2feNN=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T2",length(seq(0,47,0.5))))
gam.predict.T1feNN=predict(my.gamfeNN,newdata=newdataT1feNN,se.fit=T)
gam.predict.T2feNN=predict(my.gamfeNN,newdata=newdataT2feNN,se.fit=T)

# t critical calculation
alpha=0.05
gamfeNN.n=length(newdataT1feNN$HoursFromStart)
gamfeNN.p=length(coef(my.gamfeNN))
tcrit95.gamfeNN=qt(1-(alpha/2),gamfeNN.n-gamfeNN.p)

# calculate upper and lower bound y's using se's model fit
T1gamfeNN.UP=gam.predict.T1feNN$fit+tcrit95.gamfeNN*gam.predict.T1feNN$se.fit
T1gamfeNN.LO=gam.predict.T1feNN$fit-tcrit95.gamfeNN*gam.predict.T1feNN$se.fit
T2gamfeNN.UP=gam.predict.T2feNN$fit+tcrit95.gamfeNN*gam.predict.T2feNN$se.fit
T2gamfeNN.LO=gam.predict.T2feNN$fit-tcrit95.gamfeNN*gam.predict.T2feNN$se.fit

# FeII figure
#CairoPDF("AllCationsNotNormalizedvsTime.pdf", width= 11, height=14,pointsize=14)
par(mfrow=c(3,2))
plot(data$Fe2ppm[data$TubID=="T1"]~data$HoursFromStart[data$TubID=="T1"],pch=16,ylim=c(-1,3),xlab="",ylab="FeII (mg/L)")
lines(newdataT1feNN$HoursFromStart,gam.predict.T1feNN$fit,lwd=3,col="black")
lines(newdataT1feNN$HoursFromStart,T1gamfeNN.UP,lwd=0.5,col="black")
lines(newdataT1feNN$HoursFromStart,T1gamfeNN.LO,lwd=0.5,col="black")
points(data$Fe2ppm[data$TubID=="T2"]~data$HoursFromStart[data$TubID=="T2"],pch=1)
lines(newdataT2feNN$HoursFromStart,gam.predict.T2feNN$fit,lwd=3,lty=2,col="black")
lines(newdataT2feNN$HoursFromStart,T2gamfeNN.UP,lwd=0.5,col="black")
lines(newdataT2feNN$HoursFromStart,T2gamfeNN.LO,lwd=0.5,col="black")
rect(11.083,-1,22.583,-1+0.2/10,col='black')
rect(35.083,-1,46.583,-1+0.2/10,col='black')
legend("topleft",c("T1 (aerobic/anaerobic)","T2 (arobic)","T1 GAM","T2 GAM","95% confidence intervals"),pch=c(16,1,NA,NA,NA),lwd=c(NA,NA,3,3,1.5),lty=c(NA,NA,1,2,1))

# Total S gam model of both treatments together
my.gamsNN=gam(TotalSppm~s(HoursFromStart,by=as.numeric(TubID=="T1"))+s(HoursFromStart,by=as.numeric(TubID=="T2")),data=data,method="REML")
#summary(my.gamsNN)
#qqnorm(resid(my.gamsNN))
#qqline(resid(my.gamsNN)) #looks ok

# predictions
newdataT1sNN=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T1",length(seq(0,47,0.5))))
newdataT2sNN=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T2",length(seq(0,47,0.5))))
gam.predict.T1sNN=predict(my.gamsNN,newdata=newdataT1sNN,se.fit=T)
gam.predict.T2sNN=predict(my.gamsNN,newdata=newdataT2sNN,se.fit=T)

# t critical calculation
alpha=0.05
gamsNN.n=length(newdataT1sNN$HoursFromStart)
gamsNN.p=length(coef(my.gamsNN))
tcrit95.gamsNN=qt(1-(alpha/2),gamsNN.n-gamsNN.p)

# calculate upper and lower bound y's using se's model fit
T1gamsNN.UP=gam.predict.T1sNN$fit+tcrit95.gamsNN*gam.predict.T1sNN$se.fit
T1gamsNN.LO=gam.predict.T1sNN$fit-tcrit95.gamsNN*gam.predict.T1sNN$se.fit
T2gamsNN.UP=gam.predict.T2sNN$fit+tcrit95.gamsNN*gam.predict.T2sNN$se.fit
T2gamsNN.LO=gam.predict.T2sNN$fit-tcrit95.gamsNN*gam.predict.T2sNN$se.fit

# Total S figure
plot(data$TotalSppm[data$TubID=="T1"]~data$HoursFromStart[data$TubID=="T1"],pch=16,ylim=c(100,300),xlab="",ylab="Total S (mg/L)")
lines(newdataT1sNN$HoursFromStart,gam.predict.T1sNN$fit,lwd=3,col="black")
lines(newdataT1sNN$HoursFromStart,T1gamsNN.UP,lwd=0.5,col="black")
lines(newdataT1sNN$HoursFromStart,T1gamsNN.LO,lwd=0.5,col="black")
points(data$TotalSppm[data$TubID=="T2"]~data$HoursFromStart[data$TubID=="T2"],pch=1)
lines(newdataT2sNN$HoursFromStart,gam.predict.T2sNN$fit,lwd=3,lty=2,col="black")
lines(newdataT2sNN$HoursFromStart,T2gamsNN.UP,lwd=0.5,col="black")
lines(newdataT2sNN$HoursFromStart,T2gamsNN.LO,lwd=0.5,col="black")
rect(11.083,100,22.583,101,col='black')
rect(35.083,100,46.583,101,col='black')

# Mn gam model of both treatments together
my.gammnNN=gam(TotalMnppm~s(HoursFromStart,by=as.numeric(TubID=="T1"))+s(HoursFromStart,by=as.numeric(TubID=="T2")),data=data,method="REML")
#summary(my.gammnNN)
#qqnorm(resid(my.gammnNN))
#qqline(resid(my.gammnNN)) #looks ok

# predictions
newdataT1mnNN=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T1",length(seq(0,47,0.5))))
newdataT2mnNN=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T2",length(seq(0,47,0.5))))
gam.predict.T1mnNN=predict(my.gammnNN,newdata=newdataT1mnNN,se.fit=T)
gam.predict.T2mnNN=predict(my.gammnNN,newdata=newdataT2mnNN,se.fit=T)

# t critical calculation
alpha=0.05
gammnNN.n=length(newdataT1mnNN$HoursFromStart)
gammnNN.p=length(coef(my.gammnNN))
tcrit95.gammnNN=qt(1-(alpha/2),gammnNN.n-gammnNN.p)

# calculate upper and lower bound y's using se's model fit
T1gammnNN.UP=gam.predict.T1mnNN$fit+tcrit95.gammnNN*gam.predict.T1mnNN$se.fit
T1gammnNN.LO=gam.predict.T1mnNN$fit-tcrit95.gammnNN*gam.predict.T1mnNN$se.fit
T2gammnNN.UP=gam.predict.T2mnNN$fit+tcrit95.gammnNN*gam.predict.T2mnNN$se.fit
T2gammnNN.LO=gam.predict.T2mnNN$fit-tcrit95.gammnNN*gam.predict.T2mnNN$se.fit

# Mn figure
plot(data$TotalMnppm[data$TubID=="T1"]~data$HoursFromStart[data$TubID=="T1"],pch=16,ylim=c(-0.2,1.5),xlab="",ylab="Total Mn (mg/L)")
lines(newdataT1mnNN$HoursFromStart,gam.predict.T1mnNN$fit,lwd=3,col="black")
lines(newdataT1mnNN$HoursFromStart,T1gammnNN.UP,lwd=0.5,col="black")
lines(newdataT1mnNN$HoursFromStart,T1gammnNN.LO,lwd=0.5,col="black")
points(data$TotalMnppm[data$TubID=="T2"]~data$HoursFromStart[data$TubID=="T2"],pch=1)
lines(newdataT2mnNN$HoursFromStart,gam.predict.T2mnNN$fit,lwd=3,lty=2,col="black")
lines(newdataT2mnNN$HoursFromStart,T2gammnNN.UP,lwd=0.5,col="black")
lines(newdataT2mnNN$HoursFromStart,T2gammnNN.LO,lwd=0.5,col="black")
rect(11.083,-0.2,22.583,-0.22,col='black')
rect(35.083,-0.2,46.583,-0.22,col='black')

# Ca gam model of both treatments together
my.gamcaNN=gam(TotalCappm~s(HoursFromStart,by=as.numeric(TubID=="T1"))+s(HoursFromStart,by=as.numeric(TubID=="T2")),data=data,method="REML")
#summary(my.gamcaNN)
#qqnorm(resid(my.gamcaNN))
#qqline(resid(my.gamcaNN)) #looks ok

# predictions
newdataT1caNN=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T1",length(seq(0,47,0.5))))
newdataT2caNN=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T2",length(seq(0,47,0.5))))
gam.predict.T1caNN=predict(my.gamcaNN,newdata=newdataT1caNN,se.fit=T)
gam.predict.T2caNN=predict(my.gamcaNN,newdata=newdataT2caNN,se.fit=T)

# t critical calculation
alpha=0.05
gamcaNN.n=length(newdataT1caNN$HoursFromStart)
gamcaNN.p=length(coef(my.gamcaNN))
tcrit95.gamcaNN=qt(1-(alpha/2),gamcaNN.n-gamcaNN.p)

# calculate upper and lower bound y's using se's model fit
T1gamcaNN.UP=gam.predict.T1caNN$fit+tcrit95.gamcaNN*gam.predict.T1caNN$se.fit
T1gamcaNN.LO=gam.predict.T1caNN$fit-tcrit95.gamcaNN*gam.predict.T1caNN$se.fit
T2gamcaNN.UP=gam.predict.T2caNN$fit+tcrit95.gamcaNN*gam.predict.T2caNN$se.fit
T2gamcaNN.LO=gam.predict.T2caNN$fit-tcrit95.gamcaNN*gam.predict.T2caNN$se.fit

# Ca figure
plot(data$TotalCappm[data$TubID=="T1"]~data$HoursFromStart[data$TubID=="T1"],pch=16,ylim=c(0,150),xlab="",ylab="Ca (mg/L)")
lines(newdataT1caNN$HoursFromStart,gam.predict.T1caNN$fit,lwd=3,col="black")
lines(newdataT1caNN$HoursFromStart,T1gamcaNN.UP,lwd=0.5,col="black")
lines(newdataT1caNN$HoursFromStart,T1gamcaNN.LO,lwd=0.5,col="black")
points(data$TotalCappm[data$TubID=="T2"]~data$HoursFromStart[data$TubID=="T2"],pch=1)
lines(newdataT2caNN$HoursFromStart,gam.predict.T2caNN$fit,lwd=3,lty=2,col="black")
lines(newdataT2caNN$HoursFromStart,T2gamcaNN.UP,lwd=0.5,col="black")
lines(newdataT2caNN$HoursFromStart,T2gamcaNN.LO,lwd=0.5,col="black")
rect(11.083,0,22.583,2/10,col='black')
rect(35.083,0,46.583,2/10,col='black')

# K gam model of both treatments together
my.gamkNN=gam(TotalKppm~s(HoursFromStart,by=as.numeric(TubID=="T1"))+s(HoursFromStart,by=as.numeric(TubID=="T2")),data=data,method="REML")
#summary(my.gamkNN)
#qqnorm(resid(my.gamkNN))
#qqline(resid(my.gamkNN)) #looks ok

# predictions
newdataT1kNN=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T1",length(seq(0,47,0.5))))
newdataT2kNN=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T2",length(seq(0,47,0.5))))
gam.predict.T1kNN=predict(my.gamkNN,newdata=newdataT1kNN,se.fit=T)
gam.predict.T2kNN=predict(my.gamkNN,newdata=newdataT2kNN,se.fit=T)

# t critical calculation
alpha=0.05
gamkNN.n=length(newdataT1kNN$HoursFromStart)
gamkNN.p=length(coef(my.gamkNN))
tcrit95.gamkNN=qt(1-(alpha/2),gamkNN.n-gamkNN.p)

# calculate upper and lower bound y's using se's model fit
T1gamkNN.UP=gam.predict.T1kNN$fit+tcrit95.gamkNN*gam.predict.T1kNN$se.fit
T1gamkNN.LO=gam.predict.T1kNN$fit-tcrit95.gamkNN*gam.predict.T1kNN$se.fit
T2gamkNN.UP=gam.predict.T2kNN$fit+tcrit95.gamkNN*gam.predict.T2kNN$se.fit
T2gamkNN.LO=gam.predict.T2kNN$fit-tcrit95.gamkNN*gam.predict.T2kNN$se.fit

# K figure
plot(data$TotalKppm[data$TubID=="T1"]~data$HoursFromStart[data$TubID=="T1"],pch=16,ylim=c(1,6),xlab="Hours from start",ylab="K (mg/L)")
lines(newdataT1kNN$HoursFromStart,gam.predict.T1kNN$fit,lwd=3,col="black")
lines(newdataT1kNN$HoursFromStart,T1gamkNN.UP,lwd=0.5,col="black")
lines(newdataT1kNN$HoursFromStart,T1gamkNN.LO,lwd=0.5,col="black")
points(data$TotalKppm[data$TubID=="T2"]~data$HoursFromStart[data$TubID=="T2"],pch=1)
lines(newdataT2kNN$HoursFromStart,gam.predict.T2kNN$fit,lwd=3,lty=2,col="black")
lines(newdataT2kNN$HoursFromStart,T2gamkNN.UP,lwd=0.5,col="black")
lines(newdataT2kNN$HoursFromStart,T2gamkNN.LO,lwd=0.5,col="black")
rect(11.083,1,22.583,1.01,col='black')
rect(35.083,1,46.583,1.01,col='black')

# Mg gam model of both treatments together
my.gammgNN=gam(TotalMgppm~s(HoursFromStart,by=as.numeric(TubID=="T1"))+s(HoursFromStart,by=as.numeric(TubID=="T2")),data=data,method="REML")
#summary(my.gammgNN)
#qqnorm(resid(my.gammgNN))
#qqline(resid(my.gammgNN)) #looks ok

# predictions
newdataT1mgNN=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T1",length(seq(0,47,0.5))))
newdataT2mgNN=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T2",length(seq(0,47,0.5))))
gam.predict.T1mgNN=predict(my.gammgNN,newdata=newdataT1mgNN,se.fit=T)
gam.predict.T2mgNN=predict(my.gammgNN,newdata=newdataT2mgNN,se.fit=T)

# t critical calculation
alpha=0.05
gammgNN.n=length(newdataT1mgNN$HoursFromStart)
gammgNN.p=length(coef(my.gammgNN))
tcrit95.gammgNN=qt(1-(alpha/2),gammgNN.n-gammgNN.p)

# calculate upper and lower bound y's using se's model fit
T1gammgNN.UP=gam.predict.T1mgNN$fit+tcrit95.gammgNN*gam.predict.T1mgNN$se.fit
T1gammgNN.LO=gam.predict.T1mgNN$fit-tcrit95.gammgNN*gam.predict.T1mgNN$se.fit
T2gammgNN.UP=gam.predict.T2mgNN$fit+tcrit95.gammgNN*gam.predict.T2mgNN$se.fit
T2gammgNN.LO=gam.predict.T2mgNN$fit-tcrit95.gammgNN*gam.predict.T2mgNN$se.fit

# Mg figure
plot(data$TotalMgppm[data$TubID=="T1"]~data$HoursFromStart[data$TubID=="T1"],pch=16,ylim=c(5,15),xlab="Hours from start",ylab="Mg (mg/L)")
lines(newdataT1mgNN$HoursFromStart,gam.predict.T1mgNN$fit,lwd=3,col="black")
lines(newdataT1mgNN$HoursFromStart,T1gammgNN.UP,lwd=0.5,col="black")
lines(newdataT1mgNN$HoursFromStart,T1gammgNN.LO,lwd=0.5,col="black")
points(data$TotalMgppm[data$TubID=="T2"]~data$HoursFromStart[data$TubID=="T2"],pch=1)
lines(newdataT2mgNN$HoursFromStart,gam.predict.T2mgNN$fit,lwd=3,lty=2,col="black")
lines(newdataT2mgNN$HoursFromStart,T2gammgNN.UP,lwd=0.5,col="black")
lines(newdataT2mgNN$HoursFromStart,T2gammgNN.LO,lwd=0.5,col="black")
rect(11.083,5,22.583,5+0.2/10,col='black')
rect(35.083,5,46.583,5+0.2/10,col='black')
#dev.off() # writes file to working directory

# ----
## pH vs Dissolved Oxygen Analysis
# ----

# plot and linear model
DOvspHlm=lm(DOmgL~pH,data=data)
summary(DOvspHlm)

# ----
## pH GAM (pH vs Hours from Start, For Suppliment)
# ----

# gam model of both treatments together
my.gamph=gam(pH~s(HoursFromStart,by=as.numeric(TubID=="T1"))+s(HoursFromStart,by=as.numeric(TubID=="T2")),data=data)
summary(my.gamph)
#plot(my.gamph)
#qqnorm(resid(my.gamph))
#qqline(resid(my.gamph)) #looks normal

# predictions
newdataT1ph=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T1",length(seq(0,47,0.5))))
newdataT2ph=data.frame(HoursFromStart=seq(0,47,0.5),TubID=rep("T2",length(seq(0,47,0.5))))
gam.predict.T1ph=predict(my.gamph,newdata=newdataT1ph,se.fit=T)
gam.predict.T2ph=predict(my.gamph,newdata=newdataT2ph,se.fit=T)

# t critical calculation
alpha=0.05
gamph.n=length(newdataT1ph$HoursFromStart)
gamph.p=length(coef(my.gamph))
tcrit95=qt(1-(alpha/2),gamph.n-gamph.p)

# calculate upper and lower bound y's using se's model fit
T1gamph.UP=gam.predict.T1ph$fit+tcrit95*gam.predict.T1ph$se.fit
T1gamph.LO=gam.predict.T1ph$fit-tcrit95*gam.predict.T1ph$se.fit
T2gamph.UP=gam.predict.T2ph$fit+tcrit95*gam.predict.T2ph$se.fit
T2gamph.LO=gam.predict.T2ph$fit-tcrit95*gam.predict.T2ph$se.fit

# ----
## Figure of pH GAM (pH vs Hours from Start, For Suppliment)
# ----

#CairoPDF("pHvsTimeGAM.pdf", width= 10, height=10,pointsize=14)
par(mfrow=c(1,1))
plot(data$pH[data$TubID=="T1"]~data$HoursFromStart[data$TubID=="T1"],pch=16,ylim=c(5,10),xlab="Hours from start",ylab="pH")
lines(newdataT1ph$HoursFromStart,gam.predict.T1ph$fit,lwd=3,col="black")
lines(newdataT1ph$HoursFromStart,T1gamph.UP,lwd=0.5,col="black")
lines(newdataT1ph$HoursFromStart,T1gamph.LO,lwd=0.5,col="black")
points(data$pH[data$TubID=="T2"]~data$HoursFromStart[data$TubID=="T2"],pch=1)
lines(newdataT2ph$HoursFromStart,gam.predict.T2ph$fit,lwd=3,lty=2,col="black")
lines(newdataT2ph$HoursFromStart,T2gamph.UP,lwd=0.5,col="black")
lines(newdataT2ph$HoursFromStart,T2gamph.LO,lwd=0.5,col="black")
legend("topright",c("T1 (aerobic/anaerobic)","T2 (aerobic)","T1 GAM","T2 GAM","95% confidence intervals"),pch=c(16,1,NA,NA,NA),lwd=c(NA,NA,3,3,1.5),lty=c(NA,NA,1,2,1))
#dev.off() # writes file to working directory

# ----
# AIC Scores for Analytes vs Hours from Start (Normalized)
# ----

AIC(my.gamp,my.gamca,my.gamk,my.gammg,my.gammn,my.gamfe,my.gams,my.gamph)

# ----
## Phosphate vs Ca Linear Regression and Plot (T1 Only)
# ----

precipP.lm1=lm(AvgSRPppmAdj~TotalCappmAdj,data=T1data)
summary(precipP.lm1)

# ----
## Phosphate vs K Linear Regression (T1 Only)
# ----

precipP.lm2=lm(AvgSRPppmAdj~TotalKppmAdj,data=T1data)
summary(precipP.lm2)

# ----
## Phosphate vs Mg Linear Regression (T1 Only)
# ----

precipP.lm3=lm(AvgSRPppmAdj~TotalMgppmAdj,data=T1data)
summary(precipP.lm3)

# ----
## Phosphate vs Mn Linear Regression (T1 Only)
# ----

precipP.lm4=lm(AvgSRPppmAdj~TotalMnppmAdj,data=T1data)
summary(precipP.lm4)

# ----
## Phosphate vs FeII Linear Regression (T1 Only)
# ----

precipP.lm5=lm(AvgSRPppmAdj~Fe2ppmAdj,data=T1data)
summary(precipP.lm5)

# ----
## Phosphate vs Total S Linear Regression (T1 Only)
# ----

precipP.lm6=lm(AvgSRPppmAdj~TotalSppmAdj,data=T1data)
summary(precipP.lm6)

# ----
## AIC Scores for Phosphate vs Analytes Linear Regression (T1 Only)
# ----

AIC(precipP.lm1,precipP.lm2,precipP.lm3,precipP.lm4,precipP.lm5,precipP.lm6)

# ----
## Polyphosphate Biofilm Extract Analysis
# ----

# import data
PPdata=read.table("PPextAll_oct2014_forPaper.txt",header=T)

# adjust [P] for water, acid, and oxidizer added per 1g wet BF (units = mg P/g wet BF)
# *assume ppm=mg/L
#PPdata$Pmgperg=(PPdata$myPppmFix*(PPdata$VolAddedmL/1000))/PPdata$wetBFg
PPdata$Pmgperg=PPdata$myPppmFix/PPdata$wetBFg

# normalize to SA (mg P/m2)
PPdata$PmgpergadjSAm2=PPdata$Pmgperg/PPdata$AvgSAm2

# Anova and Tukey Test
PPbf.aov=aov(PmgpergadjSAm2~factor(SampleID),data=PPdata)
summary(PPbf.aov)
TukeyHSD(PPbf.aov)

# averages start and end
PPavgTable=as.list(tapply(PPdata$PmgpergadjSAm2,PPdata$SampleID,mean))
PPsdTable=as.list(tapply(PPdata$PmgpergadjSAm2,PPdata$SampleID,sd))

# ----
## Total P Biofilm Extract Analysis
# ----

# import data
TPdata=read.table("TPextAll_oct2014_forPaper.txt",header=T)

# adjust [P] for water, acid, and oxidizer added per 1g wet BF (units = mg P/g wet BF)
# *assume ppm=mg/L
#TPdata$Pmgperg=(TPdata$myPppmFix*(TPdata$VolAddedmL/1000))/TPdata$wetBFg
TPdata$Pmgperg=TPdata$myPppmFix/TPdata$wetBFg #Hunter suggests not multiplying by 26.5mL

# normalize to SA (mg P/m2)
TPdata$PmgpergadjSAm2=TPdata$Pmgperg/TPdata$AvgSAm2

# Anova and Tukey Test
TPbf.aov=aov(PmgpergadjSAm2~factor(SampleID),data=TPdata)
summary(TPbf.aov)
TukeyHSD(TPbf.aov)

# averages start and end
TPavgTable=tapply(TPdata$PmgpergadjSAm2,TPdata$SampleID,mean)
TPsdTable=tapply(TPdata$PmgpergadjSAm2,TPdata$SampleID,sd)

# ----
## Percent PAOs Cell Count Analysis
# ----

# import data
celldata=read.table("cellCounts_oct2014_forPaper.txt",header=T)

# define variables
filterdiam=16 #diameter of filtered solution on filter (mm)
filterareamm2=pi*(filterdiam/2)^2 #area of filter (mm^2)
filterareaum2=filterareamm2*1000^2 #area of filter (um*2)
viewht=83.78 #field of view height (um)
viewwd=111.3 #field of view width (um)
viewareaum2=viewht*viewwd #field of view area (um^2)
amtfiltered=1 #amount of solution filtered (ml)

# calculate the number of DNA and polyP cells per ml
celldata$DNAcellsml=(celldata$DNACounts*(filterareaum2/viewareaum2)*celldata$Dillution)/amtfiltered
celldata$PolyPcellsml=(celldata$PolyPCounts*(filterareaum2/viewareaum2)*celldata$Dillution)/amtfiltered
celldata$PerPolyP=(celldata$PolyPcellsml/celldata$DNAcellsml)*100

# t-test for differences in percent polyP cell counts
# check normality of T1
#hist((celldata$PerPolyP[celldata$SampleID=="T1"]))
#qqnorm((celldata$PerPolyP[celldata$SampleID=="T1"]))
#qqline((celldata$PerPolyP[celldata$SampleID=="T1"])) # looks normal

# check normality of T2
#hist((celldata$PerPolyP[celldata$SampleID=="T2"]))
#qqnorm((celldata$PerPolyP[celldata$SampleID=="T2"]))
#qqline((celldata$PerPolyP[celldata$SampleID=="T2"])) # looks normal

# check variance (here we put greater variance on top)
cellsAIRvar=var((celldata$PerPolyP[celldata$SampleID=="T2"]),na.rm=T) #larger (thus on top)
cellsAIRlth=length(na.omit(celldata$PerPolyP[celldata$SampleID=="T2"]))
cellsN2var=var((celldata$PerPolyP[celldata$SampleID=="T1"]),na.rm=T)
cellsN2lth=length(na.omit(celldata$PerPolyP[celldata$SampleID=="T1"]))

# 2 * (1 - pf(var(data2)/var(data1), length(data1)-1, length(data2)-1))
2*(1-pf(cellsAIRvar/cellsN2var,cellsN2lth-1,cellsAIRlth-1)) #=4.91e-6 
# reject that variances are equal if p < 0.05

# run two-sided t test with unequal var (because pval=4.91e-6)
t.test((celldata$PerPolyP[celldata$SampleID=="T2"]),(celldata$PerPolyP[celldata$SampleID=="T1"]),alternative="two.sided",var.equal=F)
# pval=0.0036 mean of air is sign different than mean of n2 (when 'two.sided')

# ----
## Figure Percent PAOs Cell Count
# ----

#CairoPDF("cellCounts.pdf", width= 10, height=10,pointsize=14)
par(mfrow=c(1,1))
boxplot(PerPolyP~SampleID,data=celldata,col=c('grey','white'),cex.axis=1.5,cex.lab=1.5,ylab=c("Percent of biofilm cells with polyP at end (%)"),ylim=c(0,80),names=c("T1 (anaerobic)","T2 (aerobic)"))
#dev.off() # writes file to working directory



