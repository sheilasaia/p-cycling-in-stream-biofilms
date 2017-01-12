#
# set working directory
setwd("C:\\Users\\Sheila\\Documents\\Cornell\\Research\\PhD!\\biofilms\\Oct2014_experiment\\forPaper\\forGitHub\\srp_anslysis")
#
# read in data
pdata1=read.table("pData1.txt",header=T) #DP sample periods 1 to 12
pdata2=read.table("pData2.txt",header=T) #DP sample periods 13 to 25
pdata3=read.table("pData3.txt",header=T) #hot water polyP biofilm ext (10x diluted)
pdata4=read.table("pData4.txt",header=T) #TP biofilm ext (10x diluted)
#
# read in calibrations
pcalib1=read.table("pCalib1.txt",header=T)
pcalib2=read.table("pCalib2.txt",header=T)
pcalib3=read.table("pCalib3.txt",header=T) # for data 3 and 4
#
# read in check calibrations
pcheck1=read.table("pCheck1.txt",header=T) 
pcheck2=read.table("pCheck2.txt",header=T) 
pcheck3=read.table("pCheck3.txt",header=T) # for data 3 and 4
#
# read in filter method blanks: di and di filtered (dif)
pDI1=read.table("pDI1.txt",header=T)
pDI2=read.table("pDI2.txt",header=T)
pDI3=read.table("pDI3.txt",header=T) #hot water ext polyP method blank
pDI4=read.table("pDI4.txt",header=T) #TP method blank
#
#
# --------------------------------------------------------------
# linear model data 1 (samples 0 to 12, all treatments)
p.lm1=lm(KnownPppm~HeightAu,data=pcalib1)
summary(p.lm1)
#
# save coeff's
inter1=p.lm1$coeff[1]
slope1=p.lm1$coeff[2]
R2.1=summary(p.lm1)$r.squared #check this should be close to 1.00
#
# make prediction
p.predict1=predict(p.lm1)
#
# plot calibration points and model
plot(pcalib1$KnownPppm~pcalib1$HeightAu,pch=16) 
lines(p.predict1~pcalib1$HeightAu,lwd=3,col="red")
#
# add check points
points(pcheck1$KnownPppm~pcheck1$HeightAu,pch=1)
legend("topleft",c("calibration","checks"),pch=c(16,1))
#
# recalc di conc's
pDI1$myPppm=pDI1$HeightAu*slope1+inter1
concPinFilter1=mean(pDI1$myPppm[pDI1$SampleID=="dif"]-pDI1$myPppm[pDI1$SampleID=="di"])
# if negative set to zero
if(concPinFilter1[1]<0) {
  concPinFilter1=0
}
#
# recalc p conc's in data1
pdata1$myPppm=(pdata1$HeightAu*slope1+inter1)-concPinFilter1
#
#
# --------------------------------------------------------------
# linear model data 2  (samples 13 to 25, all treatments)
p.lm2=lm(KnownPppm~HeightAu,data=pcalib2)
summary(p.lm2)
#
# save coeff's
inter2=p.lm2$coeff[1]
slope2=p.lm2$coeff[2]
R2.2=summary(p.lm2)$r.squared #check this should be close to 1.00
#
# make prediction
p.predict2=predict(p.lm2)
#
# plot calibration points and model
plot(pcalib2$KnownPppm~pcalib2$HeightAu,pch=16) 
lines(p.predict2~pcalib2$HeightAu,lwd=3,col="red")
#
# add check points
points(pcheck2$KnownPppm~pcheck2$HeightAu,pch=1)
legend("topleft",c("calibration","checks"),pch=c(16,1))
#
# recalc di conc's
pDI2$myPppm=pDI2$HeightAu*slope2+inter2
concPinFilter2=mean(pDI2$myPppm[pDI2$SampleID=="dif"]-pDI2$myPppm[pDI2$SampleID=="di"])
# if negative set to zero
if(concPinFilter2[1]<0) {
  concPinFilter2=0
}
#
# recalc p conc's in data2
pdata2$myPppm=(pdata2$HeightAu*slope2+inter2)-concPinFilter2
#
#
# --------------------------------------------------------------
# linear model data 3 (hot water ext polyP)
p.lm3=lm(KnownPppm~HeightAu,data=pcalib3)
summary(p.lm3)
#
# save coeff's
inter3=p.lm3$coeff[1]
slope3=p.lm3$coeff[2]
R2.3=summary(p.lm3)$r.squared #check this should be close to 1.00
#
# make prediction
p.predict3=predict(p.lm3)
#
# plot calibration points and model
plot(pcalib3$KnownPppm~pcalib3$HeightAu,pch=16) 
lines(p.predict3~pcalib3$HeightAu,lwd=3,col="red")
#
# add check points
points(pcheck3$KnownPppm~pcheck3$HeightAu,pch=1)
legend("topleft",c("calibration","checks"),pch=c(16,1))
#
# recalc di conc's
pDI3$myPppm=pDI3$HeightAu*slope3+inter3
concPinFilter3=mean(pDI3$myPppm[pDI3$SampleID=="diPP"]-pDI3$myPppm[pDI3$SampleID=="di"])
# if negative set to zero
if(concPinFilter3[1]<0) {
  concPinFilter3=0
}
#
# recalc p conc's in data3
pdata3$myPppm=(pdata3$HeightAu*slope3+inter3)-concPinFilter3
pdata3$myPppmFix=pdata3$myPppm #space holder for dilution section
#
# dilution
dilutionFactor=10 #10x diluted
for (i in 1:dim(pdata3)[1]){
  if(pdata3$myPppm[i]!="NA") {
    pdata3$myPppmFix[i]=as.numeric(pdata3$myPppm[i])*dilutionFactor
  }
}
#
#
# --------------------------------------------------------------
# linear model data 4 (TP ext)
# use lm3 b/c was analyzed at same time
#
# recalc di conc's
pDI4$myPppm=pDI4$HeightAu*slope3+inter3
concPinFilter4=mean(pDI4$myPppm[pDI4$SampleID=="diTP"]-pDI4$myPppm[pDI4$SampleID=="di"])
# if negative set to zero
if(concPinFilter4[1]<0) {
  concPinFilter4=0
}
#
# recalc p conc's in data3
pdata4$myPppm=(pdata4$HeightAu*slope3+inter3)-concPinFilter4
pdata4$myPppmFix=pdata4$myPppm #space holder for dilution section
#
# dilution
dilutionFactor=10 #10x diluted
for (i in 1:dim(pdata4)[1]){
  if(pdata4$myPppm[i]!="NA") {
    pdata4$myPppmFix[i]=as.numeric(pdata4$myPppm[i])*dilutionFactor
  }
}
#
#
# -----------------------------------------------------
#
# export newly calculated data
write.table(pdata1, "pNewData1.txt", sep="\t")
write.table(pdata2, "pNewData2.txt", sep="\t")
write.table(pdata3, "pNewData3.txt", sep="\t")
write.table(pdata4, "pNewData4.txt", sep="\t")
