# Load Libraries, Data, and Results ####
library(raster)
library(maptools)
library(rworldmap)
library(rcarbon)
library(Bchron)
library(coda)
library(latex2exp)

# Load Results & Data
load('../R_image_files/koreanC14.RData')
load('../R_image_files/spd_test_results.RData')
load('../R_image_files/agedepthmodels.RData')
load('../R_image_files/resABC_laplace_general.RData')
load('../R_image_files/resABC_laplace_coastal.RData')
load('../R_image_files/resABC_laplace_inland.RData')
load('../R_image_files/prior_predictive_check.RData')
load('../R_image_files/predcheck_results_general.RData')
load('../R_image_files/predcheck_results_coastal.RData')
load('../R_image_files/predcheck_results_inland.RData')

## Figure S1 - NHST Logistic and Exponential ####
pdf(file = "./figureS1.pdf",width = 9,height = 6)
par(mfcol=c(2,3))
plot(exp.general.test)
lines(exp.general.test$fit,lty=2,lwd=1,col='grey22')
title('Exponential (Combined)')
legend('topleft',bg = 'white',legend=c('Observed SPD','Fitted Model','Simulation Envelope','Positive Deviation','Negative Deviation'),lwd=c(1,1,5,5,5),col=c(1,1,'lightgrey',rgb(0.7,0,0,0.2),rgb(0,0,0.7,0.2)),cex=0.65,bty='o',box.lwd=0,lty=c(1,2,1,1,1))
legend('topright',legend=c(paste0('Global Pvalue=',round(exp.general.test$pval,6))),cex=0.9,bty='n')

plot(logistic.general.test)
lines(logistic.general.test$fit,lty=2,lwd=1,col='grey22')
legend('topright',legend=c(paste0('Global Pvalue=',
                                  round(logistic.general.test$pval,6))),cex=0.9,bty='n')
title('Logistic (Combined)')

plot(exp.coastal.test)
lines(exp.coastal.test$fit,lty=2,lwd=1,col='grey22')
title('Exponential (Coastal)')
legend('topright',legend=c(paste0('Global Pvalue=',round(exp.coastal.test$pval,6))),cex=0.9,bty='n')

plot(logistic.coastal.test)
lines(logistic.coastal.test$fit,lty=2,lwd=1,col='grey22')
title('Logistic (Coastal)')
legend('topright',legend=c(paste0('Global Pvalue=',round(logistic.coastal.test$pval,6))),cex=0.9,bty='n')

plot(exp.inland.test)
lines(exp.inland.test$fit,lty=2,lwd=1,col='grey22')
title('Exponential (inland)')
legend('topright',legend=c(paste0('Global Pvalue=',round(exp.inland.test$pval,6))),cex=0.9,bty='n')

plot(logistic.inland.test)
lines(logistic.inland.test$fit,lty=2,lwd=1,col='grey22')
title('Logistic (inland)')
legend('topright',legend=c(paste0('Global Pvalue=', round(logistic.inland.test$pval,6))),cex=0.9,bty='n')
dev.off()



## Figure S2 - Millet Dates ####
pdf(file = "./figureS2.pdf",width = 5,height = 8)
millet.dates = subset(koreaC14,koreaC14$milletAsso==TRUE)
millet.dates$col=ifelse(millet.dates$region=='inland',"#FC8D62","#66C2A5")
millet.caldates = calibrate(millet.dates$c14age, millet.dates$c14error,verbose=F,ids=millet.dates$labcode)
multiplot(millet.caldates,decreasing = TRUE, label = TRUE,col.fill=millet.dates$col,gapFactor = 0.2,rescale=TRUE,cex.id=0.3)
legend('bottomright',legend=c('Coastal Dates','Inland Dates'),fill=c("#66C2A5","#FC8D62"),bty='n')
dev.off()
## Figure S3 - Millet SPD ####
pdf(file = "./figureS3.pdf",width = 6,height = 5)
millet.spd = stackspd(millet.caldates,group=millet.dates$region,timeRange = c(7000,3000),runm=100)
plot(millet.spd)
med.dates =medCal(millet.caldates)
barCodes(med.dates,yrng=c(0,0.002),col=rgb(0,0,0,0.5),width =10)
dev.off()

## Figure S4 - Millet Permutation Test ####
pdf(file = "./figureS4.pdf",width = 6,height = 5)
koreaC14$millets='no'
koreaC14$millets[koreaC14$milletAsso==TRUE]='yes'
millet.permtest=permTest(caldates,marks=koreaC14$millets,timeRange=c(7000,3000),runm=100,nsim=1000)
plot(millet.permtest,focalm='yes')
legend('topright',bty='n', legend=c('Millet SPD','Null SPD','Positive Deviation','Negative Deviation'), lwd=c(1,5,5,5),col=c(1,'lightgrey',rgb(0.7,0,0,0.2),rgb(0,0,0.7,0.2)), cex=0.8,bg='white')
legend('topleft',legend=paste0('Global P-value=',round(millet.permtest$pValueList[1],5)),bty='n')
dev.off()


## Figure S5 - SST Median Time Series ####
pdf(file = "./figureS5.pdf",width = 8,height = 8)
par(mfrow=c(3,1),mar=c(2.5,4,1,1))
plot(0,xlim=c(7000,3000),type='n',ylim=c(range(SSDP102.temp$T2L_SSDP_102_uk37_SST_from_uk37)+c(-0.1,0.1)),xlab='',ylab='temperature (deg C)',axes=F,main='a')
rect(xleft=med.SSDP102.model[80],xright=med.SSDP102.model[81],ybottom=10,ytop=30,border=NA,col=rgb(0.67,0.84,0.9,0.5))
lines(med.SSDP102.model,SSDP102.temp$T2L_SSDP_102_uk37_SST_from_uk37,type='l',lty=1,col='darkgrey')
points(med.SSDP102.model,SSDP102.temp$T2L_SSDP_102_uk37_SST_from_uk37,pch=20,col='darkgrey')
points(med.SSDP102.model[80:81],SSDP102.temp$T2L_SSDP_102_uk37_SST_from_uk37[80:81],pch=20,col=1,cex=1.2)
text(4949.801,22,TeX('$1050cm \\, (a_{1})$'),cex=0.8)
text(4700.733,20.33,TeX('$1035cm \\, (a_{2})$'),cex=0.8)
axis(2)
par(mar=c(2.5,4,1,1))
plot(0,xlim=c(7000,3000),type='n',ylim=c(0.75,0.95),xlab='cal BP',ylab='AP/TP',axes=F,main='b')
rect(xleft=med.pomaeho.model[72],xright=med.pomaeho.model[73],ybottom=0,ytop=1,border=NA,col=rgb(0.67,0.84,0.9,0.5))
lines(med.pomaeho.model,pomaeho.apt$AP_T_Ratio,type='l',lty=1,col='darkgrey')
points(med.pomaeho.model,pomaeho.apt$AP_T_Ratio,pch=20,col='darkgrey')
points(med.pomaeho.model[72:73],pomaeho.apt$AP_T_Ratio[72:73],pch=20,col=1,cex=1.2)
text(4348.611,0.95,TeX('$1020cm \\, (b_{1})$'),cex=0.8)
text(4262,0.76,TeX('$1017cm \\, (b_{2})$'),cex=0.8)
axis(2)
par(mar=c(2.5,4,1,1))
plot(0,xlim=c(7000,3000),type='n',ylim=c(0.75,0.95),xlab='cal BP',ylab='AP/TP',axes=F,main='b')
rect(xleft=med.gy.model[78],xright=med.gy.model[79],ybottom=0.76,ytop=1,border=NA,col=rgb(0.67,0.84,0.9,0.5))
lines(med.gy.model,gy.apt$AP_TP_Ratio,type='l',lty=1,col='darkgrey')
points(med.gy.model,gy.apt$AP_TP_Ratio,pch=20,col='darkgrey')
points(med.gy.model[78:79],gy.apt$AP_TP_Ratio[78:79],pch=20,col=1,cex=1.2)
text(4698,0.9,TeX('$890cm \\, (c_{1})$'),cex=0.8)
text(4500,0.82,TeX('$887cm \\, (c_{2})$'),cex=0.8)
axis(2)
axis(1,at=seq(7000,3000,-1000),line=-1)
axis(1,at=seq(7000,3000,-500),labels=NA,line=-1)
axis(1,at=seq(7000,3000,-100),tck=-0.01,labels=NA,line=-1)
mtext('Cal BP',side=1,line=1.5,cex=0.8)
dev.off()
## FIgure S6 - Posterior Events  ####
pdf(file = "./figureS6.pdf",width = 7,height = 10)

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
plotEventPosterior(x=point_a1,hpd=0.90,xlim=c(5500,3500),main=TeX('$Event \\, a_{1}$'),legsize = 0.8)
plotEventPosterior(x=point_a2,hpd=0.90,xlim=c(5500,3500),main=TeX('$Event \\, a_{2}$'),legsize = 0.8)
plotEventPosterior(x=point_b1,hpd=0.90,xlim=c(5500,3500),main=TeX('$Event \\, b_{1}$'),legloc = 'topleft',legsize = 0.8)
plotEventPosterior(x=point_b2,hpd=0.90,xlim=c(5500,3500),main=TeX('$Event \\, b_{2}$'),legloc = 'topleft',legsize = 0.8)
plotEventPosterior(x=point_c1,hpd=0.90,xlim=c(5500,3500),main=TeX('$Event \\, c_{1}$'),legloc = 'topleft',legsize = 0.8)
plotEventPosterior(x=point_c2,hpd=0.90,xlim=c(5500,3500),main=TeX('$Event \\, c_{2}$'),legloc = 'topleft',legsize = 0.8)
dev.off()


## Figure S7 - ABC Concept ####
pdf(file = "./figureS7.pdf",width = 7,height = 10)
layout(mat=matrix(c(1,4,5,6,2,4,5,6,3,4,5,6),nrow=4,ncol=3),widths = c(0.3,0.3,0.3))

par(mar=c(5,4,2,1))
set.seed(123)
gr = seq(0,0.02,length.out = 500)
dens = dexp(gr,rate=500)
plot(gr,dens,type='n',xlab=TeX('$r_1$'),ylab='Probability',axes=FALSE)
axis(1)
polygon(x=c(gr,rev(gr)),y=c(dens,rep(0,500)),border=NA,col='lightgrey')
r1_candidate=rexp(1,500)
lines(x=c(r1_candidate,r1_candidate),y=c(0,250))
text(x=r1_candidate,y=300,labels=round(r1_candidate,8),cex=0.8,adj=0)
title('a')

gr = seq(0,0.02,length.out = 500)
dens = dexp(gr,rate=500)
plot(gr,dens,type='n',xlab=TeX('$r_2$'),ylab='Probability',axes=FALSE)
axis(1,at=pretty(gr),labels=-pretty(gr))
polygon(x=c(gr,rev(gr)),y=c(dens,rep(0,500)),border=NA,col='lightgrey')
r2_candidate=rexp(1,500)
lines(x=c(r2_candidate,r2_candidate),y=c(0,250))
text(x=r2_candidate,y=270,labels=paste0('-',round(r2_candidate,8)),cex=0.8,adj=0)
title('b')

plot(1,1,type='n',xlab=TeX('$c$'),ylab='Probability',axes=FALSE,xlim=c(6500,3500),ylim=c(0,1.2))
rect(xleft=6000,xright=4000,ybottom=0,ytop=1,border=NA,col='lightgrey')
axis(1)
c_candidate=round(runif(1,4000,6000))
text(x=c_candidate,y=0.648,labels=c_candidate,cex=0.8,adj=0)
lines(x=c(c_candidate,c_candidate),y=c(0,0.6))
title('c')


CalBP = 7000:3000
tpoint = max(CalBP)-round(c_candidate)
enpoint = round(c_candidate)-min(CalBP)
increase = 1*(1 + r1_candidate)^(1:tpoint)
decrease = increase[tpoint]*(1-r2_candidate)^(1:(enpoint+1))
PrDens=c(increase,decrease)/sum(c(increase,decrease))
d = data.frame(CalBP=CalBP,PrDens=PrDens/sum(PrDens))
plot(1,1,type='n',xlab='Cal BP',ylab='Probability',axes=FALSE,xlim=c(7000,3000),ylim=range(d$PrDens))
polygon(x=c(d$CalBP,rev(d$CalBP)),y=c(d$PrDens,rep(0,nrow(d))),border=NA,col='lightgrey')
axis(1)
title('d')

class(d)='CalGrid'
d = uncalibrate(d,verbose=F)

#Uncalsample
plot(1,1,type='n',xlab='C14 Age',ylab='Probability',axes=FALSE,xlim=rev(range(d$CRA)),ylim=range(d$PrDens))
polygon(x=c(d$CRA,rev(d$CRA)),y=c(d$PrDens,rep(0,nrow(d))),border=NA,col='lightgrey')
axis(1)
bins.coastal=bins[which(koreaC14$region=='coastal')]
nbins.coastal=length(unique(bins.coastal))
sampled.c14.dates.uncal=sample(d$CRA,size=nbins.coastal,prob=d$PrDens)
sampled.c14errors = sample(subset(koreaC14,region=='coastal')$c14error,size=nbins.coastal,replace=T)

title('e')
barCodes(sampled.c14.dates.uncal,yrng = c(0,0.0003))

# #CalSample
# plot(1,1,type='n',xlab='C14 Age',ylab='Probability',axes=FALSE,xlim=rev(range(d$CRA)),ylim=range(d$Raw))
# polygon(x=c(d$CRA,rev(d$CRA)),y=c(d$Raw,rep(0,nrow(d))),border=NA,col='lightgrey')
# axis(1)
# sampled.c14.dates.cal=sample(d$CRA,size=nbins.coastal,prob=d$Raw)
# barCodes(sampled.c14.dates.cal,yrng = c(0,0.0001))
# 

#Uncalsample SPD
calibrated.uncal=calibrate(sampled.c14.dates.uncal,sampled.c14errors)
uncal.spd = spd(calibrated.uncal,timeRange = c(7000,3000),spdnormalised = TRUE)
plot(1,1,type='n',xlab='Cal BP',ylab='Summed Probability',axes=FALSE,xlim=rev(range(uncal.spd$grid$calBP)),ylim=range(uncal.spd$grid$PrDens))
lines(uncal.spd$grid$calBP,uncal.spd$grid$PrDens/(sum(uncal.spd$grid$PrDens)))
axis(1)

load('../R_image_files/koreanC14.RData')
coastal.koreaC14 = subset(koreaC14,region=='coastal')
caldates.coastal = caldates[which(koreaC14$region=='coastal')]
bins.coastal = bins[which(koreaC14$region=='coastal')]
index=thinDates(coastal.koreaC14$c14age,coastal.koreaC14$c14error,bins.coastal,size=1,thresh=0,method='random')
target.spd =spd(caldates.coastal[index],timeRange = c(7000,3000),datenormalised=TRUE,spdnormalised = TRUE,verbose=FALSE)
lines(target.spd$grid$calBP,target.spd$grid$PrDens,lty=2,col='red')

title('f')
legend('topleft',legend=c('Candidate SPD','Target SPD'),lty=c(1,2),col=c(1,2))
dev.off()

## Figure S8 - Prior Predictive Check ####
load('../R_image_files/koreanC14.RData')

koreaC14.inland = subset(koreaC14,region=='inland')
caldates.inland = caldates[which(koreaC14$region=='inland')]
bins.inland = bins[which(koreaC14$region=='inland')]
thinindex.inland = thinDates(koreaC14.inland$c14age,koreaC14.inland$c14error,bins=bins.inland,size=1,thresh=0,method='random')
target.spd.inland =spd(caldates.inland[thinindex.inland],timeRange = c(7000,3000),datenormalised=TRUE,spdnormalised = TRUE,verbose=FALSE)

koreaC14.coastal = subset(koreaC14,region=='coastal')
caldates.coastal = caldates[which(koreaC14$region=='coastal')]
bins.coastal = bins[which(koreaC14$region=='coastal')]
thinindex.coastal = thinDates(koreaC14.coastal$c14age,koreaC14.coastal$c14error,bins=bins.coastal,size=1,thresh=0,method='random')
target.spd.coastal =spd(caldates.coastal[thinindex.coastal],timeRange = c(7000,3000),datenormalised=TRUE,spdnormalised = TRUE,verbose=FALSE)

thinindex.general = thinDates(koreaC14$c14age,koreaC14$c14error,bins=bins,size=1,thresh=0,method='random')
target.spd.general = spd(caldates[thinindex.general],timeRange = c(7000,3000),datenormalised=TRUE,spdnormalised = TRUE,verbose=FALSE)


pdf(file = "./figureS8.pdf",width = 4,height = 7)
par(mfrow=c(2,1))
plot(0,0,type='n',xlim=c(7000,3000),ylim=c(0,max(prior.check.uncal.coastal)),xlab='Cal BP',ylab='Summed Probability')
apply(prior.check.uncal.coastal,2,lines,x=c(7000:3000),col=rgb(0,0,0,0.05))
lines(target.spd.coastal$grid$calBP,target.spd.coastal$grid$PrDens,col=2)
title('Prior Predictive Check - Coastal')

plot(0,0,type='n',xlim=c(7000,3000),ylim=c(0,max(prior.check.uncal.inland)),xlab='Cal BP',ylab='Summed Probability')
apply(prior.check.uncal.inland,2,lines,x=c(7000:3000),col=rgb(0,0,0,0.05))
lines(target.spd.inland$grid$calBP,target.spd.inland$grid$PrDens,col=2)
title('Prior Predictive Check - Inland')
dev.off()

## Figure S10 and S11 - Posterior Marginal Distributions ####
tol=0.01 #tolerance level
options(scipen = 9999)

# Extract posterior
post.general.uncal = abc.general[order(abc.general$euc.uncal)[1:(nrow(abc.general)*tol)],]
post.general.cal = abc.general[order(abc.general$euc.cal)[1:(nrow(abc.general)*tol)],]
post.coastal.uncal = abc.coastal[order(abc.coastal$euc.uncal)[1:(nrow(abc.coastal)*tol)],]
post.coastal.cal = abc.coastal[order(abc.coastal$euc.cal)[1:(nrow(abc.coastal)*tol)],]
post.inland.uncal = abc.inland[order(abc.inland$euc.uncal)[1:(nrow(abc.inland)*tol)],]
post.inland.cal = abc.inland[order(abc.inland$euc.cal)[1:(nrow(abc.inland)*tol)],]

# Compute HPD Intervals and KDE
r1.hpdi.general.uncal=HPDinterval(mcmc(post.general.uncal$bl),prob = 0.90)
d.r1.general.uncal= density(post.general.uncal$bl)
r1.hpdi.general.cal=HPDinterval(mcmc(post.general.cal$bl),prob = 0.90)
d.r1.general.cal= density(post.general.cal$bl)
r2.hpdi.general.uncal=HPDinterval(mcmc(post.general.uncal$br),prob = 0.90)
d.r2.general.uncal= density(post.general.uncal$br)
r2.hpdi.general.cal=HPDinterval(mcmc(post.general.cal$br),prob = 0.90)
d.r2.general.cal= density(post.general.cal$br)
c.hpdi.general.uncal=HPDinterval(mcmc(post.general.uncal$c),prob = 0.90)
d.c.general.uncal= density(post.general.uncal$c)
c.hpdi.general.cal=HPDinterval(mcmc(post.general.cal$c),prob = 0.90)
d.c.general.cal= density(post.general.cal$c)

r1.hpdi.coastal.uncal=HPDinterval(mcmc(post.coastal.uncal$bl),prob = 0.90)
d.r1.coastal.uncal= density(post.coastal.uncal$bl)
r1.hpdi.coastal.cal=HPDinterval(mcmc(post.coastal.cal$bl),prob = 0.90)
d.r1.coastal.cal= density(post.coastal.cal$bl)
r2.hpdi.coastal.uncal=HPDinterval(mcmc(post.coastal.uncal$br),prob = 0.90)
d.r2.coastal.uncal= density(post.coastal.uncal$br)
r2.hpdi.coastal.cal=HPDinterval(mcmc(post.coastal.cal$br),prob = 0.90)
d.r2.coastal.cal= density(post.coastal.cal$br)
c.hpdi.coastal.uncal=HPDinterval(mcmc(post.coastal.uncal$c),prob = 0.90)
d.c.coastal.uncal= density(post.coastal.uncal$c)
c.hpdi.coastal.cal=HPDinterval(mcmc(post.coastal.cal$c),prob = 0.90)
d.c.coastal.cal= density(post.coastal.cal$c)

r1.hpdi.inland.uncal=HPDinterval(mcmc(post.inland.uncal$bl),prob = 0.90)
d.r1.inland.uncal= density(post.inland.uncal$bl)
r1.hpdi.inland.cal=HPDinterval(mcmc(post.inland.cal$bl),prob = 0.90)
d.r1.inland.cal= density(post.inland.cal$bl)
r2.hpdi.inland.uncal=HPDinterval(mcmc(post.inland.uncal$br),prob = 0.90)
d.r2.inland.uncal= density(post.inland.uncal$br)
r2.hpdi.inland.cal=HPDinterval(mcmc(post.inland.cal$br),prob = 0.90)
d.r2.inland.cal= density(post.inland.cal$br)
c.hpdi.inland.uncal=HPDinterval(mcmc(post.inland.uncal$c),prob = 0.90)
d.c.inland.uncal= density(post.inland.uncal$c)
c.hpdi.inland.cal=HPDinterval(mcmc(post.inland.cal$c),prob = 0.90)
d.c.inland.cal= density(post.inland.cal$c)

#plot
r1_r2_plot = function(d,hpdi,med)
{
  plot(d$x,d$y,type='n',xlab='% Annual Growth Rate',ylab='Probability Density',axes=FALSE,xlim=c(0,0.006))
  axis(1,at=axTicks(1),labels=axTicks(1)*100)
  axis(2)
  hpdi.x = d$x[which(d$x>=hpdi[1]&d$x<=hpdi[2])]
  hpdi.y = d$y[which(d$x>=hpdi[1]&d$x<=hpdi[2])]
  polygon(x=c(hpdi.x,rev(hpdi.x)),y=c(hpdi.y,rep(0,length(hpdi.y))),border=NA,col='lightblue')
  polygon(x=c(d$x,rev(d$x)),y=c(d$y,rep(0,length(d$y))))
  abline(v=med,lty=2)
}

c_plot = function(d,hpdi,med)
{
  plot(d$x,d$y,type='n',xlab='Cal BP',ylab='Probability Density',axes=FALSE,xlim=c(6000,4000))
  axis(1)
  axis(2)
  hpdi.x = d$x[which(d$x>=hpdi[1]&d$x<=hpdi[2])]
  hpdi.y = d$y[which(d$x>=hpdi[1]&d$x<=hpdi[2])]
  polygon(x=c(hpdi.x,rev(hpdi.x)),y=c(hpdi.y,rep(0,length(hpdi.y))),border=NA,col='lightblue')
  polygon(x=c(d$x,rev(d$x)),y=c(d$y,rep(0,length(d$y))))
  abline(v=med,lty=2)
}

pdf(file = "./figureS9.pdf",width = 8,height = 8)
par(mfrow=c(3,3))
leg.cex=0.85
r1_r2_plot(d.r1.general.uncal,r1.hpdi.general.uncal,median(post.general.uncal$bl))
title(TeX('$r_1$ Posterior Panregional'))
legend('topright',legend=c('90% HPDI Interval',paste(round(r1.hpdi.general.uncal*100,4),collapse='-')),bty='n',cex=leg.cex)
r1_r2_plot(d.r2.general.uncal,r2.hpdi.general.uncal,median(post.general.uncal$br))
title(TeX('$r_2$ Posterior Panregional'))
legend('topright',legend=c('90% HPDI Interval',paste(round(r2.hpdi.general.uncal*100,4),collapse='-')),bty='n',cex=leg.cex)
c_plot(d.c.general.uncal,c.hpdi.general.uncal,median(post.general.uncal$c))
title(TeX('$c$ Posterior Panregional'))
legend('topright',legend=c('90% HPDI Interval',paste(rev(round(c.hpdi.general.uncal)),collapse='-')),bty='n',cex=leg.cex)

r1_r2_plot(d.r1.coastal.uncal,r1.hpdi.coastal.uncal,median(post.coastal.uncal$bl))
title(TeX('$r_1$ Posterior Coastal'))
legend('topright',legend=c('90% HPDI Interval',paste(round(r1.hpdi.coastal.uncal*100,4),collapse='-')),bty='n',cex=leg.cex)
r1_r2_plot(d.r2.coastal.uncal,r2.hpdi.coastal.uncal,median(post.coastal.uncal$br))
title(TeX('$r_2$ Posterior Coastal'))
legend('topright',legend=c('90% HPDI Interval',paste(round(r2.hpdi.coastal.uncal*100,4),collapse='-')),bty='n',cex=leg.cex)
c_plot(d.c.coastal.uncal,c.hpdi.coastal.uncal,median(post.coastal.uncal$c))
title(TeX('$c$ Posterior Coastal'))
legend('topright',legend=c('90% HPDI Interval',paste(rev(round(c.hpdi.coastal.uncal)),collapse='-')),bty='n',cex=leg.cex)

r1_r2_plot(d.r1.inland.uncal,r1.hpdi.inland.uncal,median(post.inland.uncal$bl))
title(TeX('$r_1$ Posterior Inland'))
legend('topright',legend=c('90% HPDI Interval',paste(round(r1.hpdi.inland.uncal*100,4),collapse='-')),bty='n',cex=leg.cex)
r1_r2_plot(d.r2.inland.uncal,r2.hpdi.inland.uncal,median(post.inland.uncal$br))
title(TeX('$r_2$ Posterior Inland'))
legend('topright',legend=c('90% HPDI Interval',paste(round(r2.hpdi.inland.uncal*100,4),collapse='-')),bty='n',cex=leg.cex)
c_plot(d.c.inland.uncal,c.hpdi.inland.uncal,median(post.inland.uncal$c))
title(TeX('$c$ Posterior Inland'))
legend('topright',legend=c('90% HPDI Interval',paste(rev(round(c.hpdi.inland.uncal)),collapse='-')),bty='n',cex=leg.cex)
dev.off()


pdf(file = "./figureS10.pdf",width = 8,height = 8)
par(mfrow=c(3,3))
leg.cex=0.85
r1_r2_plot(d.r1.general.cal,r1.hpdi.general.cal,median(post.general.cal$bl))
title(TeX('$r_1$ Posterior Panregional'))
legend('topright',legend=c('90% HPDI Interval',paste(round(r1.hpdi.general.cal*100,4),collapse='-')),bty='n',cex=leg.cex)
r1_r2_plot(d.r2.general.cal,r2.hpdi.general.cal,median(post.general.cal$br))
title(TeX('$r_2$ Posterior Panregional'))
legend('topright',legend=c('90% HPDI Interval',paste(round(r2.hpdi.general.cal*100,4),collapse='-')),bty='n',cex=leg.cex)
c_plot(d.c.general.cal,c.hpdi.general.cal,median(post.general.cal$c))
title(TeX('$c$ Posterior Panregional'))
legend('topright',legend=c('90% HPDI Interval',paste(rev(round(c.hpdi.general.cal)),collapse='-')),bty='n',cex=leg.cex)

r1_r2_plot(d.r1.coastal.cal,r1.hpdi.coastal.cal,median(post.coastal.cal$bl))
title(TeX('$r_1$ Posterior Coastal'))
legend('topright',legend=c('90% HPDI Interval',paste(round(r1.hpdi.coastal.cal*100,4),collapse='-')),bty='n',cex=leg.cex)
r1_r2_plot(d.r2.coastal.cal,r2.hpdi.coastal.cal,median(post.coastal.cal$br))
title(TeX('$r_2$ Posterior Coastal'))
legend('topright',legend=c('90% HPDI Interval',paste(round(r2.hpdi.coastal.cal*100,4),collapse='-')),bty='n',cex=leg.cex)
c_plot(d.c.coastal.cal,c.hpdi.coastal.cal,median(post.coastal.cal$c))
title(TeX('$c$ Posterior Coastal'))
legend('topright',legend=c('90% HPDI Interval',paste(rev(round(c.hpdi.coastal.cal)),collapse='-')),bty='n',cex=leg.cex)

r1_r2_plot(d.r1.inland.cal,r1.hpdi.inland.cal,median(post.inland.cal$bl))
title(TeX('$r_1$ Posterior Inland'))
legend('topright',legend=c('90% HPDI Interval',paste(round(r1.hpdi.inland.cal*100,4),collapse='-')),bty='n',cex=leg.cex)
r1_r2_plot(d.r2.inland.cal,r2.hpdi.inland.cal,median(post.inland.cal$br))
title(TeX('$r_2$ Posterior Inland'))
legend('topright',legend=c('90% HPDI Interval',paste(round(r2.hpdi.inland.cal*100,4),collapse='-')),bty='n',cex=leg.cex)
c_plot(d.c.inland.cal,c.hpdi.inland.cal,median(post.inland.cal$c))
title(TeX('$c$ Posterior Inland'))
legend('topright',legend=c('90% HPDI Interval',paste(rev(round(c.hpdi.inland.cal)),collapse='-')),bty='n',cex=leg.cex)
dev.off()




## Figure S11 - Posterior Predictive Checks ####

ppmedian.general=apply(ppcheck.uncal.general,1,median)
pplo.general=apply(ppcheck.uncal.general,1,quantile,0.025)
pphi.general=apply(ppcheck.uncal.general,1,quantile,0.975)

ppmedian.coastal=apply(ppcheck.uncal.coastal,1,median)
pplo.coastal=apply(ppcheck.uncal.coastal,1,quantile,0.025)
pphi.coastal=apply(ppcheck.uncal.coastal,1,quantile,0.975)

ppmedian.inland=apply(ppcheck.uncal.inland,1,median)
pplo.inland=apply(ppcheck.uncal.inland,1,quantile,0.025)
pphi.inland=apply(ppcheck.uncal.inland,1,quantile,0.975)

pdf(file = "./figureS11.pdf",width = 10,height = 3.5)
par(mfrow=c(1,3))
plot(target.spd.general$grid$calBP,target.spd.general$grid$PrDens,type='n',xlim=c(7000,3000),ylim=c(0,max(c(target.spd.general$grid$PrDens,pphi.general))),xlab='cal BP',ylab='Summed Probability')
polygon(c(7000:3000,rev(7000:3000)),c(pplo.general,rev(pphi.general)),border=NA,col='lightgrey')
lines(7000:3000,ppmedian.general,col=2,lty=2)
lines(target.spd.general$grid$calBP,target.spd.general$grid$PrDens,lwd=1)
legend('topleft',legend=c('Observed (thinned)','Median Posterior Predictive Check','95% Posterior Predictive Interval'),col=c(1,2,'lightgrey'),lwd=c(1,1,5),lty=c(1,2,1),cex=0.6,bg='white')
title('Panregional SPD Posterior Predictive Check')

plot(target.spd.coastal$grid$calBP,target.spd.coastal$grid$PrDens,type='n',xlim=c(7000,3000),ylim=c(0,max(c(target.spd.coastal$grid$PrDens,pphi.coastal))),xlab='cal BP',ylab='Summed Probability')
polygon(c(7000:3000,rev(7000:3000)),c(pplo.coastal,rev(pphi.coastal)),border=NA,col='lightgrey')
lines(7000:3000,ppmedian.coastal,col=2,lty=2)
lines(target.spd.coastal$grid$calBP,target.spd.coastal$grid$PrDens,lwd=1)
title('Coastal SPD Posterior Predictive Check')

plot(target.spd.inland$grid$calBP,target.spd.inland$grid$PrDens,type='n',xlim=c(7000,3000),ylim=c(0,max(c(target.spd.inland$grid$PrDens,pphi.inland))),xlab='cal BP',ylab='Summed Probability')
polygon(c(7000:3000,rev(7000:3000)),c(pplo.inland,rev(pphi.inland)),border=NA,col='lightgrey')
lines(7000:3000,ppmedian.inland,col=2,lty=2)
lines(target.spd.inland$grid$calBP,target.spd.inland$grid$PrDens,lwd=1)
title('Inland SPD Posterior Predictive Check')
dev.off()

