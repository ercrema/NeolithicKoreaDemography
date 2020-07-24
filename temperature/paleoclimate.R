library(Bchron)
# Read Kim 2004
kim2004.dates = read.csv('SSDP102_Kim/SSDP_102.Kim.2004-chron.csv',skip=1)
kim2004.temp  = read.csv('SSDP102_Kim/SSDP_102.Kim.2004.csv',skip=1)

# Fit Compound Poisson-Gamma chronology model
kim2004.model = Bchronology(ages=kim2004.dates$T2L_SSDP_102_c14_date - kim2004.dates$T2L_SSDP_102_delta_r,ageSds=sqrt(kim2004.dates$T2L_SSDP_102_c14_1s_err^2+kim2004.dates$T2L_SSDP_102_delta_r_1s_error^2),calCurves = rep('marine13',nrow(kim2004.dates)),ids=kim2004.dates$T2L_SSDP_102_labcode,positions=kim2004.dates$T2L_SSDP_102_depth_top,predictPositions=kim2004.temp$T2L_SSDP_102_depth)

# Median Date at each level
medPred = apply(kim2004.model$thetaPredict,2,median)
plot(medPred,kim2004.temp$T2L_SSDP_102_age_medianMedianBacon)
abline(a=0,b=1)

# Interpolate across 7000:3000
interp.temp = matrix(NA,nrow=length(7000:3000),ncol=1000)
for (i in 1:1000)
{
  interp.temp[,i]=approx(x=kim2004.model$thetaPredict[i,],y=kim2004.temp$T2L_SSDP_102_uk37_SST_from_uk37,xout=3000:7000)$y
}

med=apply(interp.temp,1,median,na.rm=TRUE)
lo=apply(interp.temp,1,min,na.rm=TRUE)
hi=apply(interp.temp,1,max,na.rm=TRUE)

# Plot Timing of Rapid Temperature Decrease
pdf(file = "./figureClimate1.pdf",width = 10,height = 5)
plot(3000:7000,med,xlim=c(7000,3000),type='n',ylim=c(range(c(lo,hi))),xlab='cal BP',ylab='temperature (deg C)')
polygon(c(3000:7000,7000:3000),c(lo,rev(hi)),border=NA,col='lightblue')
lines(3000:7000,med,col='darkblue')
lines(kim2004.temp$T2L_SSDP_102_age_medianMedianBacon,kim2004.temp$T2L_SSDP_102_uk37_SST_from_uk37,type='l',lty=2,col=1)
points(kim2004.temp$T2L_SSDP_102_age_medianMedianBacon[80:81],kim2004.temp$T2L_SSDP_102_uk37_SST_from_uk37[80:81],pch=20,col='red')
lines(kim2004.temp$T2L_SSDP_102_age_medianMedianBacon[80:81],kim2004.temp$T2L_SSDP_102_uk37_SST_from_uk37[80:81],lwd=1,col='red')
dev.off()

# Timing Cooling Event
library(coda)
coolEvent=kim2004.model$thetaPredict[,81]
dcool.hpdi=HPDinterval(mcmc(coolEvent),prob = 0.90)
d.coolEvent=density(coolEvent)

pdf(file = "./figureClimate2.pdf",width = 6,height = 5)

plot(d.coolEvent$x,d.coolEvent$y,type='n',xlab='cal BP',ylab='Probability Density',axes=FALSE,xlim=c(5500,4500))
title('Timing of Cooling Event')
axis(1)
axis(2)
hpdi.x = d.coolEvent$x[which(d.coolEvent$x>=dcool.hpdi[1]&d.coolEvent$x<=dcool.hpdi[2])]
hpdi.y = d.coolEvent$y[which(d.coolEvent$x>=dcool.hpdi[1]&d.coolEvent$x<=dcool.hpdi[2])]
polygon(x=c(hpdi.x,rev(hpdi.x)),y=c(hpdi.y,rep(0,length(hpdi.y))),border=NA,col='lightblue')
polygon(x=c(d.coolEvent$x,rev(d.coolEvent$x)),y=c(d.coolEvent$y,rep(0,length(d.coolEvent$y))))
abline(v=median(coolEvent),lty=2)
dev.off()

# Timing Difference Cooling Change in Growth Rate
load('../results_images/resABC_laplace.RData')
tol=0.05
library(coda)
post = res[order(res$euc.uncal)[1:(nrow(res)*tol)],]
diff = sample(coolEvent,size=1000,replace=TRUE)-sample(post$c,size=1000,replace=TRUE)
diff.hpdi=HPDinterval(mcmc(diff),prob = 0.90)
d.diff=density(diff)

pdf(file = "./figureClimate3.pdf",width = 6,height = 5)
plot(d.diff$x,d.diff$y,type='n',xlab='Years',ylab='Probability Density',axes=FALSE)
title('Number of Years after Cooling Event')
axis(1)
axis(2)
hpdi.x = d.diff$x[which(d.diff$x>=diff.hpdi[1]&d.diff$x<=diff.hpdi[2])]
hpdi.y = d.diff$y[which(d.diff$x>=diff.hpdi[1]&d.diff$x<=diff.hpdi[2])]
polygon(x=c(hpdi.x,rev(hpdi.x)),y=c(hpdi.y,rep(0,length(hpdi.y))),border=NA,col='lightblue')
polygon(x=c(d.diff$x,rev(d.diff$x)),y=c(d.diff$y,rep(0,length(d.diff$y))))
abline(v=median(diff),lty=2)
dev.off()



