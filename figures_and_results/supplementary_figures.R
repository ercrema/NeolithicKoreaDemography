# Load Libraries, Data, and Results ####
library(raster)
library(maptools)
library(rworldmap)
library(rcarbon)
library(Bchron)
library(truncnorm)
library(coda)
library(latex2exp)
library(here)
library(nimbleCarbon)

# Load Results & Data
load(here('R_image_files','koreanC14.RData'))
load(here('R_image_files','spd_test_results.RData'))
load(here('R_image_files','agedepthmodels.RData'))
load(here('R_image_files','mcmcdiagnostic_postpredcheck.RData'))
load(here('R_image_files','mcmc_samples_coastal.RData'))
load(here('R_image_files','mcmc_samples_inland.RData'))
load(here('R_image_files','mcmc_samples_all.RData'))


## Figure S1 - Millet Dates ####
pdf(file = here('figures_and_results','figureS1.pdf'),width = 5,height = 8)
millet.dates = subset(koreaC14,koreaC14$milletAsso==TRUE)
millet.dates$col=ifelse(millet.dates$region=='inland',"#FC8D62","#66C2A5")
millet.caldates = calibrate(millet.dates$c14age, millet.dates$c14error,verbose=F,ids=millet.dates$labcode)
multiplot(millet.caldates,decreasing = TRUE, label = TRUE,col.fill=millet.dates$col,gapFactor = 0.2,rescale=TRUE,cex.id=0.3)
legend('bottomright',legend=c('Coastal Dates','Inland Dates'),fill=c("#66C2A5","#FC8D62"),bty='n')
dev.off()
## Figure S2 - Millet SPD ####
pdf(file = here('figures_and_results','figureS2.pdf'),width = 6,height = 5)
millet.spd = stackspd(millet.caldates,group=millet.dates$region,timeRange = c(7000,3000),runm=100)
plot(millet.spd)
med.dates =medCal(millet.caldates)
barCodes(med.dates,yrng=c(0,0.002),col=rgb(0,0,0,0.5),width =10)
dev.off()

## Figure S3 - Millet Permutation Test ####
pdf(file = here('figures_and_results','figureS3.pdf'),width = 6,height = 5)
koreaC14$millets='no'
koreaC14$millets[koreaC14$milletAsso==TRUE]='yes'
millet.permtest=permTest(caldates,marks=koreaC14$millets,timeRange=c(7000,3000),runm=100,nsim=1000)
plot(millet.permtest,focalm='yes')
legend('topright',bty='n', legend=c('Millet SPD','Null SPD','Positive Deviation','Negative Deviation'), lwd=c(1,5,5,5),col=c(1,'lightgrey',rgb(0.7,0,0,0.2),rgb(0,0,0.7,0.2)), cex=0.8,bg='white')
legend('topleft',legend=paste0('Global P-value=',round(millet.permtest$pValueList[1],5)),bty='n')
dev.off()


## Figure S4 - SST Median Time Series ####
pdf(file = here('figures_and_results','figureS4.pdf'),width = 8,height = 8)
par(mfrow=c(3,1),mar=c(2.5,4,1,1))
plot(0,xlim=c(7000,3000),type='n',ylim=c(range(SSDP102.temp$T2L_SSDP_102_uk37_SST_from_uk37)+c(-0.1,0.1)),xlab='',ylab='temperature (deg C)',axes=F,main='a')
rect(xleft=med.SSDP102.model[80],xright=med.SSDP102.model[81],ybottom=10,ytop=30,border=NA,col=rgb(0.67,0.84,0.9,0.5))
lines(med.SSDP102.model,SSDP102.temp$T2L_SSDP_102_uk37_SST_from_uk37,type='l',lty=1,col='darkgrey')
points(med.SSDP102.model,SSDP102.temp$T2L_SSDP_102_uk37_SST_from_uk37,pch=20,col='darkgrey')
points(med.SSDP102.model[80:81],SSDP102.temp$T2L_SSDP_102_uk37_SST_from_uk37[80:81],pch=20,col=1,cex=1.2)
text(4949.801,22,TeX('$1050cm \\, (d_{1})$'),cex=0.8)
text(4700.733,20.33,TeX('$1035cm \\, (d_{2})$'),cex=0.8)
axis(2)
par(mar=c(2.5,4,1,1))
plot(0,xlim=c(7000,3000),type='n',ylim=c(0.75,0.95),xlab='cal BP',ylab='AP/TP',axes=F,main='b')
rect(xleft=med.pomaeho.model[72],xright=med.pomaeho.model[73],ybottom=0,ytop=1,border=NA,col=rgb(0.67,0.84,0.9,0.5))
lines(med.pomaeho.model,pomaeho.apt$AP_T_Ratio,type='l',lty=1,col='darkgrey')
points(med.pomaeho.model,pomaeho.apt$AP_T_Ratio,pch=20,col='darkgrey')
points(med.pomaeho.model[72:73],pomaeho.apt$AP_T_Ratio[72:73],pch=20,col=1,cex=1.2)
text(4348.611,0.95,TeX('$1020cm \\, (f_{1})$'),cex=0.8)
text(4262,0.76,TeX('$1017cm \\, (f_{2})$'),cex=0.8)
axis(2)
par(mar=c(2.5,4,1,1))
plot(0,xlim=c(7000,3000),type='n',ylim=c(0.75,0.95),xlab='cal BP',ylab='AP/TP',axes=F,main='b')
rect(xleft=med.gy.model[78],xright=med.gy.model[79],ybottom=0.76,ytop=1,border=NA,col=rgb(0.67,0.84,0.9,0.5))
lines(med.gy.model,gy.apt$AP_TP_Ratio,type='l',lty=1,col='darkgrey')
points(med.gy.model,gy.apt$AP_TP_Ratio,pch=20,col='darkgrey')
points(med.gy.model[78:79],gy.apt$AP_TP_Ratio[78:79],pch=20,col=1,cex=1.2)
text(4698,0.9,TeX('$890cm \\, (g_{1})$'),cex=0.8)
text(4500,0.82,TeX('$887cm \\, (g_{2})$'),cex=0.8)
axis(2)
axis(1,at=seq(7000,3000,-1000),line=-1)
axis(1,at=seq(7000,3000,-500),labels=NA,line=-1)
axis(1,at=seq(7000,3000,-100),tck=-0.01,labels=NA,line=-1)
mtext('Cal BP',side=1,line=1.5,cex=0.8)
dev.off()
## Figure S5 Age-Dept Model Kim et al ####
pdf(file = here('figures_and_results','figureS5.pdf'),width = 9.5,height = 4.5)
plot(SSDP102.model)
dev.off()

## Figure S6 Age-Dept Model Pomaeho ####
pdf(file = here('figures_and_results','figureS6.pdf'),width = 9.5,height = 4.5)
plot(pomaeho.model)
dev.off()

## Figure S7 Age-Dept Model GY ####
pdf(file = here('figures_and_results','figureS7.pdf'),width = 9.5,height = 4.5)
plot(gy.model)
dev.off()

## Figure S8 - Posterior Events  ####
pdf(file = here('figures_and_results','figureS8.pdf'),width = 7,height = 10)

plotEventPosterior = function(x,hpd,legloc='topright',legsize=1,...)
{
  hpdi.interval=HPDinterval(mcmc(x),prob = hpd)
  d.event=density(x)
  
  plot(d.event$x,d.event$y,type='n',xlab='cal BP',ylab='Probability Density',...)
  hpdi.x = d.event$x[which(d.event$x>=hpdi.interval[1]&d.event$x<=hpdi.interval[2])]
  hpdi.y = d.event$y[which(d.event$x>=hpdi.interval[1]&d.event$x<=hpdi.interval[2])]
  polygon(x=c(hpdi.x,rev(hpdi.x)),y=c(hpdi.y,rep(0,length(hpdi.y))),border=NA,col='lightblue')
  lines(d.event)
  abline(v=median(x),lty=2)
  legend(legloc,legend=c(paste0(hpd*100,'%HPDI:\n',paste(rev(hpdi.interval),collapse='-'),' cal BP')),bty='n',cex=legsize)
}

# Extract Posteriors
point_a1=SSDP102.model$thetaPredict[,81]
point_a2=SSDP102.model$thetaPredict[,80]
point_b1=pomaeho.model$thetaPredict[,73]
point_b2=pomaeho.model$thetaPredict[,72]
point_c1=gy.model$thetaPredict[,79]
point_c2=gy.model$thetaPredict[,78]

par(mfrow=c(3,2))
plotEventPosterior(x=point_a1,hpd=0.90,xlim=c(5500,3500),main=TeX('$Event \\, d_{1}$'),legsize = 0.8)
plotEventPosterior(x=point_a2,hpd=0.90,xlim=c(5500,3500),main=TeX('$Event \\, d_{2}$'),legsize = 0.8)
plotEventPosterior(x=point_b1,hpd=0.90,xlim=c(5500,3500),main=TeX('$Event \\, f_{1}$'),legloc = 'topleft',legsize = 0.8)
plotEventPosterior(x=point_b2,hpd=0.90,xlim=c(5500,3500),main=TeX('$Event \\, f_{2}$'),legloc = 'topleft',legsize = 0.8)
plotEventPosterior(x=point_c1,hpd=0.90,xlim=c(5500,3500),main=TeX('$Event \\, g_{1}$'),legloc = 'topleft',legsize = 0.8)
plotEventPosterior(x=point_c2,hpd=0.90,xlim=c(5500,3500),main=TeX('$Event \\, g_{2}$'),legloc = 'topleft',legsize = 0.8)
dev.off()


## Figure S9 - Prior Predictive Check ####
pdf(file = here('figures_and_results','figureS9.pdf'),width = 5,height = 5)
prior.params=list(r1=rnorm(500,mean=0,sd=0.0004),r2=rnorm(500,mean=0,sd=0.0004),mu=rtruncnorm(500,a=3000,b=7000,mean=5000,sd=1000))
modelPlot(model=dDoubleExponentialGrowth,params = prior.params,a = 7000,b=3000, nsample = 500,alpha=0.1,main='Prior Predictive Check')
dev.off()

## Figure S10 - Trace Plot ####
pdf(file = here('figures_and_results','figureS10.pdf'),width = 7.5,height = 7.5)
par(mfcol=c(3,3),mar=c(4.5,5.5,2,1))
plot(as.numeric(mcmc.samples.coastal$chain1[,'r1']),type='l',ylim=range(params.coastal$r1),lwd=0.3,col='lightblue',xlab='MCMC Samples',ylab='',las=1,main=TeX('$r_1\\, Coastal$'))
lines(as.numeric(mcmc.samples.coastal$chain2[,'r1']),col='pink',lwd=0.3)
lines(as.numeric(mcmc.samples.coastal$chain3[,'r1']),col='grey',lwd=0.3)
mtext(TeX('$r_1$'),2,line=4,cex=0.7,las=2)
plot(as.numeric(mcmc.samples.coastal$chain1[,'r2']),type='l',ylim=range(params.coastal$r2),lwd=0.3,col='lightblue',xlab='MCMC Samples',ylab='',las=1,main=TeX('$r_2\\, Coastal$'))
lines(as.numeric(mcmc.samples.coastal$chain2[,'r2']),col='pink',lwd=0.3)
lines(as.numeric(mcmc.samples.coastal$chain3[,'r2']),col='grey',lwd=0.3)
mtext(TeX('$r_2$'),2,line=4,cex=0.7,las=2)
plot(as.numeric(mcmc.samples.coastal$chain1[,'chp']),type='l',ylim=range(params.coastal$mu),lwd=0.3,col='lightblue',xlab='MCMC Samples',ylab='',las=1,main=TeX('$c\\, Coastal$'))
lines(as.numeric(mcmc.samples.coastal$chain2[,'chp']),col='pink',lwd=0.3)
lines(as.numeric(mcmc.samples.coastal$chain3[,'chp']),col='grey',lwd=0.3)
mtext(TeX('$c$'),2,line=4,cex=0.7,las=2)
legend('topleft',legend=c('Chain 1','Chain 2','Chain 3'),col=c('lightblue','pink','grey'),lty=1,cex=0.7,bty='n')

plot(as.numeric(mcmc.samples.inland$chain1[,'r1']),type='l',ylim=range(params.inland$r1),lwd=0.3,col='lightblue',xlab='MCMC Samples',ylab='',las=1,main=TeX('$r_1\\, Inland$'))
lines(as.numeric(mcmc.samples.inland$chain2[,'r1']),col='pink',lwd=0.3)
lines(as.numeric(mcmc.samples.inland$chain3[,'r1']),col='grey',lwd=0.3)
mtext(TeX('$r_1$'),2,line=4,cex=0.7,las=2)
plot(as.numeric(mcmc.samples.inland$chain1[,'r2']),type='l',ylim=range(params.inland$r2),lwd=0.3,col='lightblue',xlab='MCMC Samples',ylab='',las=1,main=TeX('$r_2\\, Inland$'))
lines(as.numeric(mcmc.samples.inland$chain2[,'r2']),col='pink',lwd=0.3)
lines(as.numeric(mcmc.samples.inland$chain3[,'r2']),col='grey',lwd=0.3) 
mtext(TeX('$r_2$'),2,line=4,cex=0.7,las=2)
plot(as.numeric(mcmc.samples.inland$chain1[,'chp']),type='l',ylim=range(params.inland$mu),lwd=0.3,col='lightblue',xlab='MCMC Samples',ylab='',las=1,main=TeX('$c\\, Inland$'))
lines(as.numeric(mcmc.samples.inland$chain2[,'chp']),col='pink',lwd=0.3)
lines(as.numeric(mcmc.samples.inland$chain3[,'chp']),col='grey',lwd=0.3)
mtext(TeX('$c$'),2,line=4,cex=0.7,las=2)

plot(as.numeric(mcmc.samples.all$chain1[,'r1']),type='l',ylim=range(params.all$r1),lwd=0.3,col='lightblue',xlab='MCMC Samples',ylab='',las=1,main=TeX('$r_1\\, All$'))
lines(as.numeric(mcmc.samples.all$chain2[,'r1']),col='pink',lwd=0.3)
lines(as.numeric(mcmc.samples.all$chain3[,'r1']),col='grey',lwd=0.3)
mtext(TeX('$r_1$'),2,line=4,cex=0.7,las=2)
plot(as.numeric(mcmc.samples.all$chain1[,'r2']),type='l',ylim=range(params.all$r2),lwd=0.3,col='lightblue',xlab='MCMC Samples',ylab='',las=1,main=TeX('$r_2\\, All$'))
lines(as.numeric(mcmc.samples.all$chain2[,'r2']),col='pink',lwd=0.3)
lines(as.numeric(mcmc.samples.all$chain3[,'r2']),col='grey',lwd=0.3) 
mtext(TeX('$r_2$'),2,line=4,cex=0.7,las=2)
plot(as.numeric(mcmc.samples.all$chain1[,'chp']),type='l',ylim=range(params.all$mu),lwd=0.3,col='lightblue',xlab='MCMC Samples',ylab='',las=1,main=TeX('$c\\, All$'))
lines(as.numeric(mcmc.samples.all$chain2[,'chp']),col='pink',lwd=0.3)
lines(as.numeric(mcmc.samples.all$chain3[,'chp']),col='grey',lwd=0.3)
mtext(TeX('$c$'),2,line=4,cex=0.7,las=2)

dev.off()
## Figure S11 - Fitted Model ####
pdf(file = here('figures_and_results','figureS11.pdf'),width = 7,height = 3)
par(mfrow=c(1,3))
modelPlot(dDoubleExponentialGrowth,params =params.coastal,nsample=500,alpha=0.05,a=7000,b=3000,main='a')
modelPlot(dDoubleExponentialGrowth,params =params.inland,nsample=500,alpha=0.05,a=7000,b=3000,main='b')
modelPlot(dDoubleExponentialGrowth,params =params.all,nsample=500,alpha=0.05,a=7000,b=3000,main='c')
dev.off()

## Figure S12 - Posterior Marginal Distributions ####
options(scipen=9999)
pdf(file = here('figures_and_results','figureS12.pdf'),width = 7.5,height = 7.5)
par(mfrow=c(3,3),mar=c(5,5,2,1))
postHPDplot(params.coastal$r1*100,xlab='Annual Growth Rate (%)',ylab='Probability',main=TeX('$r_1\\,Coastal$'))
postHPDplot(params.coastal$r2*100,xlab='Annual Growth Rate (%)',ylab='Probability',main=TeX('$r_2\\,Coastal$'))
postHPDplot(params.coastal$mu,xlab='cal BP',ylab='Probability',main=TeX('$c\\,Coastal$'),xlim=c(6500,4000))
postHPDplot(params.inland$r1*100,xlab='Annual Growth Rate (%)',ylab='Probability',main=TeX('$r_1\\,Inland$'))
postHPDplot(params.inland$r2*100,xlab='Annual Growth Rate (%)',ylab='Probability',main=TeX('$r_2\\,Inland$'))
postHPDplot(params.inland$mu,xlab='cal BP',ylab='Probability',main=TeX('$c\\,Inland$'),xlim=c(5400,4700))
postHPDplot(params.all$r1*100,xlab='Annual Growth Rate (%)',ylab='Probability',main=TeX('$r_1\\,All$'))
postHPDplot(params.all$r2*100,xlab='Annual Growth Rate (%)',ylab='Probability',main=TeX('$r_2\\,All$'))
postHPDplot(params.all$mu,xlab='cal BP',ylab='Probability',main=TeX('$c\\,All$'),xlim=c(6500,4000))
dev.off()




## Figure S13 - Posterior Predictive Checks ####
pdf(file = here('figures_and_results','figureS13.pdf'),width = 5,height = 6)
par(mfrow=c(3,1),mar=c(5,4,2,0.5))
plot(params.coastal.uncalsample,interval = 0.95,main='Coastal')
plot(params.inland.uncalsample,interval = 0.95,main='Inland')
plot(params.all.uncalsample,interval = 0.95,main='All')
dev.off()

