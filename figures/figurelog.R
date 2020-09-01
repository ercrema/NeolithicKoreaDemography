# Load Libraries, Data, and Results ####
library(raster)
library(maptools)
library(rworldmap)
library(rcarbon)
library(Bchron)
library(coda)

# Load Results
load('../results_images/test_results.RData')
load('../results_images/resABC_laplace_general.RData')
load('../results_images/resABC_laplace_coastal.RData')
load('../results_images/resABC_laplace_inland.RData')
load('../results_images/predcheck_results_general.RData')
load('../results_images/predcheck_results_coastal.RData')
load('../results_images/predcheck_results_inland.RData')


# Site Distribution Figure ####
load('../data/koreanC14.RData')
sites.sp <- unique(data.frame(SiteID=koreaC14$site_id,latitude=koreaC14$latitude,longitude=koreaC14$longitude))
coast<-getMap(resolution = "high")
proj4string(coast)

coordinates(sites.sp)<-c("longitude","latitude")
proj4string(sites.sp)<-proj4string(coast)
dem=getData('alt', country='KOR', mask=TRUE)
projection(dem) <- proj4string(coast)

slope <- terrain(dem, opt='slope')
aspect <- terrain(dem, opt='aspect')
hs <- hillShade(slope, aspect, 40, 270)
#plot(hs, col=grey(0:100/100), legend=FALSE)

pdf(file = "./figure_site_map.pdf",width = 5,height = 5)
plot(coast,col="grey68",xlim=extent(hs)[1:2],ylim=extent(hs)[3:4],border=NA,xlab="",ylab="")
abline(v=122:135,lwd=0.5,col='grey88')
abline(h=33:38,lwd=0.5,col='grey88')
plot(hs, col=grey(0:100/100), legend=FALSE,add=TRUE)
points(sites.sp,pch=20,cex=0.7)
box()
axis(side=1,at=122:135,cex=0.7,las=2,cex.axis=0.5,hadj=0.1,tck=-0.01)
axis(side=2,at=33:38,cex=0.7,las=2,cex.axis=0.5,hadj=-0.5,tck=-0.01)
mtext('Longitude',1,1,cex=0.7)
mtext('Latitude',2,1,cex=0.7)
dev.off()

# Stacked SPD ####
runm = 100 #smoothing window
timeRange = c(7000,3000)
combined.spd=stackspd(x=caldates,timeRange=timeRange,bins=bins,runm = runm,group = koreaC14$region)
pdf(file = "./figure_stacked_spd.pdf",width = 8,height = 4.5)
plot(combined.spd)
dev.off()

# CKDE ####
bw=100
s.coastal = sampleDates(caldates[which(koreaC14$region=='coastal')],bins=bins[which(koreaC14$region=='coastal')],nsim=1000,boot=FALSE)
s.inland = sampleDates(caldates[which(koreaC14$region=='inland')],bins=bins[which(koreaC14$region=='inland')],nsim=1000,boot=FALSE)
ckde.coastal = ckde(s.coastal,timeRange=c(7000,3000),bw=bw)
ckde.inland = ckde(s.inland,timeRange=c(7000,3000),bw=bw)

pdf(file = "./figure_ckde.pdf",width = 7,height = 6)
par(mfrow=c(2,1),mar=c(4,4,3,1))
plot(ckde.inland)
title('Bootstrapped CKDE (Inland Sites)',cex.main=1,line=1)
legend('topright',bty='n',legend=c(paste0('n(bins)=',length(unique(bins[which(koreaC14$region=='inland')]))),paste0('n(dates)=',length(caldates[which(koreaC14$region=='inland')])),paste0('bw=',bw)),cex=0.7)
legend('topleft',legend=c('Average CKDE','95% C.I.'),lwd=c(2,5),lty=c(2,1),col=c(1,'lightgrey'),cex=0.7,bty='n')
plot(ckde.coastal)
title('Bootstrapped CKDE (Coastal Sites)',cex.main=1,line=1)
legend('topright',bty='n',legend=c(paste0('n(bins)=',length(unique(bins[which(koreaC14$region=='coastal')]))),paste0('n(dates)=',length(caldates[which(koreaC14$region=='coastal')])),paste0('bw=',bw)),cex=0.7)
dev.off()

# Permutation Tests ####
pdf(file = "./figure_permtest.pdf",width = 7,height = 6)
par(mfrow=c(2,1))
plot(coastal.inland.permtest,focalm='coastal')
title('Coastal SPD',cex.main=1,line=1)
plot(coastal.inland.permtest,focalm='inland')
title('Inland SPD',cex.main=1,line=1)
legend('topright',bty='n',legend=c('Observed SPD','Null SPD','Positive Deviation','Negative Deviation'),lwd=c(1,5,5,5),col=c(1,'lightgrey',rgb(0.7,0,0,0.2),rgb(0,0,0.7,0.2)),cex=0.7,bg='white')
dev.off()


# Model Tests ####
pdf(file = "./figure_modelTests.pdf",width = 9,height = 6)
par(mfcol=c(2,3))
plot(exp.general.test)
lines(exp.general.test$fit,lty=2,lwd=1,col='grey22')
title('Exponential (Combined)')
legend('topleft',bg = 'white',legend=c('Observed SPD','Fitted Model','Simulation Envelope','Positive Deviation','Negative Deviation'),lwd=c(1,1,5,5,5),col=c(1,1,'lightgrey',rgb(0.7,0,0,0.2),rgb(0,0,0.7,0.2)),cex=0.68,bty='o',box.lwd=0,lty=c(1,2,1,1,1))
plot(logistic.general.test)
lines(logistic.general.test$fit,lty=2,lwd=1,col='grey22')
title('Logistic (Combined)')

plot(exp.coastal.test)
lines(exp.coastal.test$fit,lty=2,lwd=1,col='grey22')
title('Exponential (Coastal)')
plot(logistic.coastal.test)
lines(logistic.coastal.test$fit,lty=2,lwd=1,col='grey22')
title('Logistic (Coastal)')

plot(exp.inland.test)
lines(exp.inland.test$fit,lty=2,lwd=1,col='grey22')
title('Exponential (inland)')
plot(logistic.inland.test)
lines(logistic.inland.test$fit,lty=2,lwd=1,col='grey22')
title('Logistic (inland)')
dev.off()


# Posterior Distributions (General) ####
load('../results_images/resABC_laplace_general.RData')
tol=0.05
library(coda)
post = abc.general[order(abc.general$euc.uncal)[1:(nrow(abc.general)*tol)],]
options(scipen = 9999)
pdf(file = "./figure_posterior_general.pdf",width = 10,height = 3.5)
par(mfrow=c(1,3))
#bl
bl.hpdi=HPDinterval(mcmc(post$bl),prob = 0.90)
d.bl=density(post$bl)
plot(d.bl$x,d.bl$y,type='n',xlab='% Annual Growth Rate',ylab='Probability Density',axes=FALSE)
title('Growing Phase Growth Rate Posterior')
axis(1,at=axTicks(1),labels=axTicks(1)*100)
axis(2)
hpdi.x = d.bl$x[which(d.bl$x>=bl.hpdi[1]&d.bl$x<=bl.hpdi[2])]
hpdi.y = d.bl$y[which(d.bl$x>=bl.hpdi[1]&d.bl$x<=bl.hpdi[2])]
polygon(x=c(hpdi.x,rev(hpdi.x)),y=c(hpdi.y,rep(0,length(hpdi.y))),border=NA,col='lightblue')
polygon(x=c(d.bl$x,rev(d.bl$x)),y=c(d.bl$y,rep(0,length(d.bl$y))))
abline(v=median(post$bl),lty=2)

#br
br.hpdi=HPDinterval(mcmc(post$br),prob = 0.90)
d.br=density(post$br)
plot(d.br$x,d.br$y,type='n',xlab='% Growth Rate',ylab='Probability Density',axes=FALSE)
title('Declining Phase Growth Rate Posterior')
axis(1,at=axTicks(1),labels=-axTicks(1)*100)
axis(2)
hpdi.x = d.br$x[which(d.br$x>=br.hpdi[1]&d.br$x<=br.hpdi[2])]
hpdi.y = d.br$y[which(d.br$x>=br.hpdi[1]&d.br$x<=br.hpdi[2])]
polygon(x=c(hpdi.x,rev(hpdi.x)),y=c(hpdi.y,rep(0,length(hpdi.y))),border=NA,col='lightblue')
polygon(x=c(d.br$x,rev(d.br$x)),y=c(d.br$y,rep(0,length(d.br$y))))
abline(v=median(post$br),lty=2)

#c
c.hpdi=HPDinterval(mcmc(post$c),prob = 0.90)
d.c=density(post$c)
plot(d.c$x,d.c$y,type='n',xlab='Cal BP',ylab='Probability Density',axes=FALSE,xlim=rev(range(d.c$x)))
title('Change Point Posterior')
axis(1)
axis(2)
hpdi.x = d.c$x[which(d.c$x>=c.hpdi[1]&d.c$x<=c.hpdi[2])]
hpdi.y = d.c$y[which(d.c$x>=c.hpdi[1]&d.c$x<=c.hpdi[2])]
polygon(x=c(hpdi.x,rev(hpdi.x)),y=c(hpdi.y,rep(0,length(hpdi.y))),border=NA,col='lightblue')
polygon(x=c(d.c$x,rev(d.c$x)),y=c(d.c$y,rep(0,length(d.c$y))))
abline(v=median(post$c),lty=2)
dev.off()

# Posterior Distributions (Coastal vs Inland ) ####
load('../results_images/resABC_laplace_coastal.RData')
load('../results_images/resABC_laplace_inland.RData')

tol=0.05
library(coda)
post.coastal = abc.coastal[order(abc.coastal$euc.uncal)[1:(nrow(abc.coastal)*tol)],]
post.inland = abc.inland[order(abc.inland$euc.uncal)[1:(nrow(abc.inland)*tol)],]
options(scipen = 9999)
pdf(file = "./figure_posterior_coastal_vs_inland.pdf",width = 10,height = 3.5)
par(mfrow=c(1,3))

#bl
bl.hpdi.coastal=HPDinterval(mcmc(post.coastal$bl),prob = 0.90)
bl.hpdi.inland=HPDinterval(mcmc(post.inland$bl),prob = 0.90)

d.bl.coastal=density(post.coastal$bl)
d.bl.inland=density(post.inland$bl,bw=d.bl.coastal$bw)

plot(0,0,type='n',xlab='% Annual Growth Rate',ylab='Probability Density',axes=FALSE,xlim=range(c(d.bl.coastal$x,d.bl.inland$x)),ylim=range(c(d.bl.coastal$y,d.bl.inland$y)))

title('Growing Phase Growth Rate Posterior')
axis(1,at=axTicks(1),labels=axTicks(1)*100)
axis(2)
hpdi.x.coastal = d.bl.coastal$x[which(d.bl.coastal$x>=bl.hpdi.coastal[1]&d.bl.coastal$x<=bl.hpdi.coastal[2])]
hpdi.x.inland = d.bl.inland$x[which(d.bl.inland$x>=bl.hpdi.inland[1]&d.bl.inland$x<=bl.hpdi.inland[2])]
hpdi.y.coastal = d.bl.coastal$y[which(d.bl.coastal$x>=bl.hpdi.coastal[1]&d.bl.coastal$x<=bl.hpdi.coastal[2])]
hpdi.y.inland = d.bl.inland$y[which(d.bl.inland$x>=bl.hpdi.inland[1]&d.bl.inland$x<=bl.hpdi.inland[2])]

polygon(x=c(hpdi.x.coastal,rev(hpdi.x.coastal)),y=c(hpdi.y.coastal,rep(0,length(hpdi.y.coastal))),border=NA,col=rgb(0.6784314,0.8470588,0.9019608,0.5))
polygon(x=c(d.bl.coastal$x,rev(d.bl.coastal$x)),y=c(d.bl.coastal$y,rep(0,length(d.bl.coastal$y))),border='lightblue')

polygon(x=c(hpdi.x.inland,rev(hpdi.x.inland)),y=c(hpdi.y.inland,rep(0,length(hpdi.y.inland))),border=NA,col=rgb(1.0000000,0.7529412,0.7960784,0.5))
polygon(x=c(d.bl.inland$x,rev(d.bl.inland$x)),y=c(d.bl.inland$y,rep(0,length(d.bl.inland$y))),border='lightpink')

abline(v=median(post.coastal$bl),lty=2,col='blue')
abline(v=median(post.inland$bl),lty=2,col='red')
legend('topright',legend=c('Coastal','Inland'),fill=c('lightblue','lightpink'))


#br
br.hpdi.coastal=HPDinterval(mcmc(post.coastal$br),prob = 0.90)
br.hpdi.inland=HPDinterval(mcmc(post.inland$br),prob = 0.90)

d.br.coastal=density(post.coastal$br)
d.br.inland=density(post.inland$br,bw=d.br.coastal$bw)

plot(0,0,type='n',xlab='% Annual Growth Rate',ylab='Probability Density',axes=FALSE,xlim=range(c(d.br.coastal$x,d.br.inland$x)),ylim=range(c(d.br.coastal$y,d.br.inland$y)))

title('Declining Phase Growth Rate Posterior')
axis(1,at=axTicks(1),labels=-axTicks(1)*100)
axis(2)
hpdi.x.coastal = d.br.coastal$x[which(d.br.coastal$x>=br.hpdi.coastal[1]&d.br.coastal$x<=br.hpdi.coastal[2])]
hpdi.x.inland = d.br.inland$x[which(d.br.inland$x>=br.hpdi.inland[1]&d.br.inland$x<=br.hpdi.inland[2])]
hpdi.y.coastal = d.br.coastal$y[which(d.br.coastal$x>=br.hpdi.coastal[1]&d.br.coastal$x<=br.hpdi.coastal[2])]
hpdi.y.inland = d.br.inland$y[which(d.br.inland$x>=br.hpdi.inland[1]&d.br.inland$x<=br.hpdi.inland[2])]

polygon(x=c(hpdi.x.coastal,rev(hpdi.x.coastal)),y=c(hpdi.y.coastal,rep(0,length(hpdi.y.coastal))),border=NA,col=rgb(0.6784314,0.8470588,0.9019608,0.5))
polygon(x=c(d.br.coastal$x,rev(d.br.coastal$x)),y=c(d.br.coastal$y,rep(0,length(d.br.coastal$y))),border='lightblue')

polygon(x=c(hpdi.x.inland,rev(hpdi.x.inland)),y=c(hpdi.y.inland,rep(0,length(hpdi.y.inland))),border=NA,col=rgb(1.0000000,0.7529412,0.7960784,0.5))
polygon(x=c(d.br.inland$x,rev(d.br.inland$x)),y=c(d.br.inland$y,rep(0,length(d.br.inland$y))),border='lightpink')

abline(v=median(post.coastal$br),lty=2,col='blue')
abline(v=median(post.inland$br),lty=2,col='red')


#c
c.hpdi.coastal=HPDinterval(mcmc(post.coastal$c),prob = 0.90)
c.hpdi.inland=HPDinterval(mcmc(post.inland$c),prob = 0.90)

d.c.coastal=density(post.coastal$c)
d.c.inland=density(post.inland$c,bw=d.c.coastal$bw)


plot(0,0,type='n',xlab='Cal BP',ylab='Probability Density',axes=FALSE,xlim=rev(range(c(d.c.coastal$x,d.c.inland$x))),ylim=range(c(d.c.coastal$y,d.c.inland$y)))



title('Change Point Posterior')
axis(1)
axis(2)

hpdi.x.coastal = d.c.coastal$x[which(d.c.coastal$x>=c.hpdi.coastal[1]&d.c.coastal$x<=c.hpdi.coastal[2])]
hpdi.x.inland = d.c.inland$x[which(d.c.inland$x>=c.hpdi.inland[1]&d.c.inland$x<=c.hpdi.inland[2])]
hpdi.y.coastal = d.c.coastal$y[which(d.c.coastal$x>=c.hpdi.coastal[1]&d.c.coastal$x<=c.hpdi.coastal[2])]
hpdi.y.inland = d.c.inland$y[which(d.c.inland$x>=c.hpdi.inland[1]&d.c.inland$x<=c.hpdi.inland[2])]

polygon(x=c(hpdi.x.coastal,rev(hpdi.x.coastal)),y=c(hpdi.y.coastal,rep(0,length(hpdi.y.coastal))),border=NA,col=rgb(0.6784314,0.8470588,0.9019608,0.5))
polygon(x=c(d.c.coastal$x,rev(d.c.coastal$x)),y=c(d.c.coastal$y,rep(0,length(d.c.coastal$y))),border='lightblue')

polygon(x=c(hpdi.x.inland,rev(hpdi.x.inland)),y=c(hpdi.y.inland,rep(0,length(hpdi.y.inland))),border=NA,col=rgb(1.0000000,0.7529412,0.7960784,0.5))
polygon(x=c(d.c.inland$x,rev(d.c.inland$x)),y=c(d.c.inland$y,rep(0,length(d.c.inland$y))),border='lightpink')

abline(v=median(post.coastal$c),lty=2,col='blue')
abline(v=median(post.inland$c),lty=2,col='red')
dev.off()


# Posterior Predictive Check (General) ####
load('../results_images/predcheck_results_general.RData')
load('../data/koreanC14.RData')
library(rcarbon)
observed = spd(caldates,bins,timeRange=c(7000,3000),spdnormalised = TRUE)
ppmedian=apply(ppcheck.cal.general,1,median)
pplo=apply(ppcheck.cal.general,1,quantile,0.025)
pphi=apply(ppcheck.cal.general,1,quantile,0.975)

pdf(file = "./figure_ppcheck_general.pdf",width = 4,height = 4)
plot(observed$grid$calBP,observed$grid$PrDens,type='n',xlim=c(7000,3000),ylim=c(0,max(c(observed$grid$PrDens,pphi))),xlab='cal BP',ylab='Summed Probability')
polygon(c(7000:3000,rev(7000:3000)),c(pplo,rev(pphi)),border=NA,col='lightgrey')
lines(7000:3000,ppmedian,col=2,lty=2)
lines(observed$grid$calBP,observed$grid$PrDens,lwd=1)
legend('topleft',legend=c('Observed','Median Posterior Predictive Check','95% Posterior Predictive Interval'),col=c(1,2,'lightgrey'),lwd=c(1,1,5),lty=c(1,2,1),cex=0.6,bg='white')
dev.off()

# Posterior Predictive Check (Coastal vs Inland) ####
load('../results_images/predcheck_results_coastal.RData')
load('../results_images/predcheck_results_inland.RData')
load('../data/koreanC14.RData')
coastal.index = which(koreaC14$region=='coastal')
inland.index = which(koreaC14$region=='inland')

observed.coastal = spd(caldates[coastal.index],bins[coastal.index],timeRange=c(7000,3000),spdnormalised = TRUE)
observed.inland = spd(caldates[inland.index],bins[inland.index],timeRange=c(7000,3000),spdnormalised = TRUE)

ppmedian.coastal=apply(ppcheck.cal.coastal,1,median)
pplo.coastal=apply(ppcheck.cal.coastal,1,quantile,0.025)
pphi.coastal=apply(ppcheck.cal.coastal,1,quantile,0.975)
ppmedian.inland=apply(ppcheck.cal.inland,1,median)
pplo.inland=apply(ppcheck.cal.inland,1,quantile,0.025)
pphi.inland=apply(ppcheck.cal.inland,1,quantile,0.975)

pdf(file = "./figure_ppcheck_coastal_vs_inland.pdf",width = 4,height = 8)
par(mfrow=c(2,1))
plot(observed.coastal$grid$calBP,observed.coastal$grid$PrDens,type='n',xlim=c(7000,3000),ylim=c(0,max(c(observed.coastal$grid$PrDens,pphi.coastal))),xlab='cal BP',ylab='Summed Probability')
polygon(c(7000:3000,rev(7000:3000)),c(pplo.coastal,rev(pphi.coastal)),border=NA,col='lightgrey')
lines(7000:3000,ppmedian.coastal,col=2,lty=2)
lines(observed.coastal$grid$calBP,observed.coastal$grid$PrDens,lwd=1)
legend('topleft',legend=c('Observed','Median Posterior Predictive Check','95% Posterior Predictive Interval'),col=c(1,2,'lightgrey'),lwd=c(1,1,5),lty=c(1,2,1),cex=0.6,bg='white')

plot(observed.inland$grid$calBP,observed.inland$grid$PrDens,type='n',xlim=c(7000,3000),ylim=c(0,max(c(observed.inland$grid$PrDens,pphi.inland))),xlab='cal BP',ylab='Summed Probability')
polygon(c(7000:3000,rev(7000:3000)),c(pplo.inland,rev(pphi.inland)),border=NA,col='lightgrey')
lines(7000:3000,ppmedian.inland,col=2,lty=2)
lines(observed.inland$grid$calBP,observed.inland$grid$PrDens,lwd=1)
legend('topleft',legend=c('Observed','Median Posterior Predictive Check','95% Posterior Predictive Interval'),col=c(1,2,'lightgrey'),lwd=c(1,1,5),lty=c(1,2,1),cex=0.6,bg='white')

dev.off()


# Age-Depth Model of Kim 2004 ####
# Read Kim 2004
kim2004.dates = read.csv('../data/SSDP_102.Kim.2004-chron.csv',skip=1)
kim2004.temp  = read.csv('../data/SSDP_102.Kim.2004.csv',skip=1)

# Read marine20 calibration curve
marine20 = read.csv('http://intcal.org/curves/marine20.14c', encoding="UTF-8",skip=11,header=F)
createCalCurve(name='marine20',calAges=marine20[,1],uncalAges=marine20[,2],oneSigma=marine20[,3])
file.copy(from = 'marine20.rda',to = system.file('data',package='Bchron'))

kim2004.model.marine20 = Bchronology(ages=round(kim2004.dates$T2L_SSDP_102_c14_date - kim2004.dates$T2L_SSDP_102_delta_r),ageSds=round(sqrt(kim2004.dates$T2L_SSDP_102_c14_1s_err^2+kim2004.dates$T2L_SSDP_102_delta_r_1s_error^2)),calCurves = rep('marine20',nrow(kim2004.dates)),ids=kim2004.dates$T2L_SSDP_102_labcode,positions=kim2004.dates$T2L_SSDP_102_depth_top,predictPositions=kim2004.temp$T2L_SSDP_102_depth)

# Extract Median Dates
med.bchron.marine20 = apply(kim2004.model.marine20$thetaPredict,2,median)


# Plot Median Predicted Temperature Change
pdf(file = "./figure_kim2004_reanalysis.pdf",width = 6,height = 5.5)
layout(matrix(c(1,2,1,3),2,2))
par(mar=c(5,4,1,1))
plot(0,xlim=c(7000,3000),type='n',ylim=c(range(kim2004.temp$T2L_SSDP_102_uk37_SST_from_uk37)+c(-0.1,0.1)),xlab='cal BP',ylab='temperature (deg C)',axes=F)
rect(xleft=med.bchron.marine20[80],xright=med.bchron.marine20[81],ybottom=10,ytop=30,border=NA,col=rgb(0.67,0.84,0.9,0.5))
lines(med.bchron.marine20,kim2004.temp$T2L_SSDP_102_uk37_SST_from_uk37,type='l',lty=1,col='darkgrey')
points(med.bchron.marine20,kim2004.temp$T2L_SSDP_102_uk37_SST_from_uk37,pch=20,col='darkgrey')
points(med.bchron.marine20[80:81],kim2004.temp$T2L_SSDP_102_uk37_SST_from_uk37[80:81],pch=20,col=1,cex=1.2)
text(4949.801,21.89769,'a')
text(4700.733,20.45735,'b')
#lines(med.bchron.marine20[80:81],kim2004.temp$T2L_SSDP_102_uk37_SST_from_uk37[80:81],lwd=2,lty=2,col=1)
axis(1,at=seq(7000,3000,-1000),tck=-0.1)
axis(1,at=seq(7000,3000,-500),tck=-0.07,labels=NA)
axis(1,at=seq(7000,3000,-100),tck=-0.05,labels=NA)
axis(2)
point_a=kim2004.model.marine20$thetaPredict[,81]
point_b=kim2004.model.marine20$thetaPredict[,80]

dcool.hpdi.a=HPDinterval(mcmc(point_a),prob = 0.90)
d.event.a=density(point_a)
dcool.hpdi.b=HPDinterval(mcmc(point_b),prob = 0.90)
d.event.b=density(point_b)

plot(d.event.a$x,d.event.a$y,type='n',xlab='cal BP',ylab='Probability Density',axes=FALSE,xlim=c(5500,4000))
legend('topright',legend='Timing of a',bty='n',cex=0.8)
axis(1)
axis(1,at=seq(7000,3000,-100),tck=-0.03,labels=NA)

axis(2)
hpdi.x = d.event.a$x[which(d.event.a$x>=dcool.hpdi.a[1]&d.event.a$x<=dcool.hpdi.a[2])]
hpdi.y = d.event.a$y[which(d.event.a$x>=dcool.hpdi.a[1]&d.event.a$x<=dcool.hpdi.a[2])]
polygon(x=c(hpdi.x,rev(hpdi.x)),y=c(hpdi.y,rep(0,length(hpdi.y))),border=NA,col='lightblue')
polygon(x=c(d.event.a$x,rev(d.event.a$x)),y=c(d.event.a$y,rep(0,length(d.event.a$y))))
abline(v=median(point_a),lty=2)

plot(d.event.b$x,d.event.b$y,type='n',xlab='cal BP',ylab='Probability Density',axes=FALSE,xlim=c(5500,4000))
legend('topright',legend='Timing of b',bty='n',cex=0.8)
axis(1)
axis(2)
axis(1,at=seq(7000,3000,-100),tck=-0.03,labels=NA)

hpdi.x = d.event.b$x[which(d.event.b$x>=dcool.hpdi.b[1]&d.event.b$x<=dcool.hpdi.b[2])]
hpdi.y = d.event.b$y[which(d.event.b$x>=dcool.hpdi.b[1]&d.event.b$x<=dcool.hpdi.b[2])]
polygon(x=c(hpdi.x,rev(hpdi.x)),y=c(hpdi.y,rep(0,length(hpdi.y))),border=NA,col='lightblue')
polygon(x=c(d.event.b$x,rev(d.event.b$x)),y=c(d.event.b$y,rep(0,length(d.event.b$y))))
abline(v=median(point_b),lty=2)
dev.off()

# Millet SPD ####
millet.dates = subset(koreaC14,koreaC14$milletAsso==TRUE)
millet.caldates = calibrate(millet.dates$c14age,millet.dates$c14error,ids=millet.dates$labcode)
millet.spd = stackspd(millet.caldates,group=millet.dates$region,timeRange = c(7000,3000),runm=100)
pdf(file = "./figure_millet_spd.pdf",width = 6,height = 5)
plot(millet.spd)
med.dates =medCal(millet.caldates)
barCodes(med.dates,yrng=c(0,0.002),col=rgb(0,0,0,0.5),width =10)
dev.off()

pdf(file = "./figure_millet_multiplot.pdf",width = 5,height = 7)
col=ifelse(millet.dates$region=='inland',rgb(0.99,0.55,0.38),rgb(0.40,0.76,0.65))
multiplot(millet.caldates,decreasing = TRUE,credMass = TRUE,label = TRUE,col.fill=col,gapFactor = 0.2,rescale=TRUE,cex.id=0.3)
legend('bottomright',legend=c('Coastal Dates','Inland Dates'),fill=c(rgb(0.40,0.76,0.65),rgb(0.99,0.55,0.38)),bty='n')
dev.off()


# Event Comparison Plot ####
# Millet SPDs
millet.dates = subset(koreaC14,koreaC14$milletAsso==TRUE)
millet.caldates = calibrate(millet.dates$c14age,millet.dates$c14error)
millet.spd = stackspd(millet.caldates,group=millet.dates$region,timeRange = c(7000,3000),runm=100)
coastal.millet = millet.spd$spds$coastal$grid$PrDens
inland.millet = millet.spd$spds$inland$grid$PrDens
total.millet = coastal.millet + inland.millet

pdf(file = "./figure_event_comparisons.pdf",width = 6,height = 9)
par(mfrow=c(3,1),mar=c(0,4,1,1))
plot(0,0,xlim=c(7000,3000),ylim=range(total.millet),axes=FALSE,xlab='',ylab='Summed Probability')
polygon(c(3000:7000,7000:3000),c(inland.millet,rep(0,length(inland.millet))),col=rgb(0.99,0.55,0.38),lwd=0.5,border=rgb(0.99,0.55,0.38))
polygon(c(3000:7000,7000:3000),c(total.millet,rev(inland.millet)),col=rgb(0.40,0.76,0.65),border=rgb(0.40,0.76,0.65),lwd=0.5)
axis(2)
legend('topleft',legend='a',bty='n',cex=2)

# Cooling Event
plot(0,0,xlim=c(7000,3000),ylim=c(0,max(d.event.a$y,d.event.b$y)),axes=FALSE,xlab='',ylab='Probability')
hpdi.x = d.event.a$x[which(d.event.a$x>=dcool.hpdi.a[1]&d.event.a$x<=dcool.hpdi.a[2])]
hpdi.y = d.event.a$y[which(d.event.a$x>=dcool.hpdi.a[1]&d.event.a$x<=dcool.hpdi.a[2])]
polygon(x=c(d.event.a$x,rev(d.event.a$x)),y=c(d.event.a$y,rep(0,length(d.event.a$y))),col=rgb(0.67,0.84,0.9,0.5),border=NA)

hpdi.x = d.event.b$x[which(d.event.b$x>=dcool.hpdi.b[1]&d.event.b$x<=dcool.hpdi.b[2])]
hpdi.y = d.event.b$y[which(d.event.b$x>=dcool.hpdi.b[1]&d.event.b$x<=dcool.hpdi.b[2])]
polygon(x=c(d.event.b$x,rev(d.event.b$x)),y=c(d.event.b$y,rep(0,length(d.event.b$y))),col=rgb(1,0.71,0.76,0.5),border=NA)

abline(v=median(point_a),lty=2)
text(median(point_a)+50,y=0.0016,label='Cooling Start (Median)',srt=90,cex=0.8)
abline(v=median(point_b),lty=2)
text(median(point_b)+50,y=0.0016,label='Cooling End (Median)',srt=90,cex=0.8)
axis(2)
legend('topleft',legend='b',bty='n',cex=2)


# Change Point
par(mar=c(4,4,1,1))
plot(0,0,xlim=c(7000,3000),ylim=c(0,max(d.c.inland$y,d.c.coastal$y)),axes=FALSE,xlab='',ylab='Probability')
polygon(x=c(d.c.coastal$x,rev(d.c.coastal$x)),y=c(d.c.coastal$y,rep(0,length(d.c.coastal$y))),col=rgb(0.99,0.55,0.38,0.5))
polygon(x=c(d.c.inland$x,rev(d.c.inland$x)),y=c(d.c.inland$y,rep(0,length(d.c.inland$y))),col=rgb(0.40,0.76,0.65,0.5))
axis(2)
axis(1,at=seq(7000,3000,-1000),tck=-0.05)
axis(1,at=seq(7000,3000,-500),tck=-0.02,labels=NA)
axis(1,at=seq(7000,3000,-100),tck=-0.01,labels=NA)
mtext('Cal BP',side=1,line=3,cex=0.7)
legend('topleft',legend='c',bty='n',cex=2)
legend('topright',legend=c('Inland','Coastal'),fill=c(rgb(0.40,0.76,0.65),rgb(0.99,0.55,0.38)),bty='n')
dev.off()



# Temporal Distance Plots ####
a_coast=sample(point_a)-sample(post.coastal$c,size=length(point_a))
a_coast.hpdi_left=c(HPDinterval(mcmc(a_coast),prob = 0.90)[1],0)
a_coast.hpdi_right=c(0,HPDinterval(mcmc(a_coast),prob = 0.90)[2])
a_coast.dens=density(a_coast)
a_coast.hpdi.left.x = a_coast.dens$x[which(a_coast.dens$x>=a_coast.hpdi_left[1]&a_coast.dens$x<=a_coast.hpdi_left[2])]
a_coast.hpdi.left.y = a_coast.dens$y[which(a_coast.dens$x>=a_coast.hpdi_left[1]&a_coast.dens$x<=a_coast.hpdi_left[2])]
a_coast.hpdi.right.x = a_coast.dens$x[which(a_coast.dens$x>=a_coast.hpdi_right[1]&a_coast.dens$x<=a_coast.hpdi_right[2])]
a_coast.hpdi.right.y = a_coast.dens$y[which(a_coast.dens$x>=a_coast.hpdi_right[1]&a_coast.dens$x<=a_coast.hpdi_right[2])]


b_coast=sample(point_b)-sample(post.coastal$c,size=length(point_a))
b_coast.hpdi_left=c(HPDinterval(mcmc(b_coast),prob = 0.90)[1],0)
b_coast.hpdi_right=c(0,HPDinterval(mcmc(b_coast),prob = 0.90)[2])
b_coast.dens=density(b_coast)
b_coast.hpdi.left.x = b_coast.dens$x[which(b_coast.dens$x>=b_coast.hpdi_left[1]&b_coast.dens$x<=b_coast.hpdi_left[2])]
b_coast.hpdi.left.y = b_coast.dens$y[which(b_coast.dens$x>=b_coast.hpdi_left[1]&b_coast.dens$x<=b_coast.hpdi_left[2])]
b_coast.hpdi.right.x = b_coast.dens$x[which(b_coast.dens$x>=b_coast.hpdi_right[1]&b_coast.dens$x<=b_coast.hpdi_right[2])]
b_coast.hpdi.right.y = b_coast.dens$y[which(b_coast.dens$x>=b_coast.hpdi_right[1]&b_coast.dens$x<=b_coast.hpdi_right[2])]

a_inland=sample(point_a)-sample(post.inland$c,size=length(point_a))
a_inland.hpdi_left=c(HPDinterval(mcmc(a_inland),prob = 0.90)[1],0)
a_inland.hpdi_right=c(0,HPDinterval(mcmc(a_inland),prob = 0.90)[2])
a_inland.dens=density(a_inland)
a_inland.hpdi.left.x = a_inland.dens$x[which(a_inland.dens$x>=a_inland.hpdi_left[1]&a_inland.dens$x<=a_inland.hpdi_left[2])]
a_inland.hpdi.left.y = a_inland.dens$y[which(a_inland.dens$x>=a_inland.hpdi_left[1]&a_inland.dens$x<=a_inland.hpdi_left[2])]
a_inland.hpdi.right.x = a_inland.dens$x[which(a_inland.dens$x>=a_inland.hpdi_right[1]&a_inland.dens$x<=a_inland.hpdi_right[2])]
a_inland.hpdi.right.y = a_inland.dens$y[which(a_inland.dens$x>=a_inland.hpdi_right[1]&a_inland.dens$x<=a_inland.hpdi_right[2])]


b_inland=sample(point_b)-sample(post.inland$c,size=length(point_a))
b_inland.hpdi_left=c(HPDinterval(mcmc(b_inland),prob = 0.90)[1],0)
b_inland.hpdi_right=c(0,HPDinterval(mcmc(b_inland),prob = 0.90)[2])
b_inland.dens=density(b_inland)
b_inland.hpdi.left.x = b_inland.dens$x[which(b_inland.dens$x>=b_inland.hpdi_left[1]&b_inland.dens$x<=b_inland.hpdi_left[2])]
b_inland.hpdi.left.y = b_inland.dens$y[which(b_inland.dens$x>=b_inland.hpdi_left[1]&b_inland.dens$x<=b_inland.hpdi_left[2])]
b_inland.hpdi.right.x = b_inland.dens$x[which(b_inland.dens$x>=b_inland.hpdi_right[1]&b_inland.dens$x<=b_inland.hpdi_right[2])]
b_inland.hpdi.right.y = b_inland.dens$y[which(b_inland.dens$x>=b_inland.hpdi_right[1]&b_inland.dens$x<=b_inland.hpdi_right[2])]



pdf(file = "./figure_climate_vs_changepoint.pdf",width = 7,height = 7)
par(mfrow=c(2,2),mar=c(5,4,2,1))
plot(a_coast.dens$x,a_coast.dens$y,type='n',xlab='Years',ylab='Probability Density',axes=FALSE,xlim=c(-1000,1000),main='Coastal Changepoint vs Event A')
polygon(x=c(a_coast.hpdi.left.x,rev(a_coast.hpdi.left.x)),y=c(a_coast.hpdi.left.y,rep(0,length(a_coast.hpdi.left.y))),border=NA,col='lightblue')
polygon(x=c(a_coast.hpdi.right.x,rev(a_coast.hpdi.right.x)),y=c(a_coast.hpdi.right.y,rep(0,length(a_coast.hpdi.right.y))),border=NA,col='lightpink')

polygon(x=c(a_coast.dens$x,rev(a_coast.dens$x)),y=c(a_coast.dens$y,rep(0,length(a_coast.dens$y))),border='lightgrey')
axis(1,at=seq(-1000,1000,200),labels=abs(seq(-1000,1000,200)))
axis(2)
box()
text(x=500,y=median(par('usr')[3:4]),label=paste('Changepoint after\n P=',round(sum(a_coast>0)/1000,2)),cex=0.8)
text(x=-500,y=median(par('usr')[3:4]),label=paste('Changepoint before\n P=',round(sum(a_coast<0)/1000,2)),cex=0.8)
abline(v=0,lty=2,lwd=2)


plot(b_coast.dens$x,b_coast.dens$y,type='n',xlab='Years',ylab='Probability Density',axes=FALSE,xlim=c(-1000,1000),main='Coastal Changepoint vs Event B')
polygon(x=c(b_coast.hpdi.left.x,rev(b_coast.hpdi.left.x)),y=c(b_coast.hpdi.left.y,rep(0,length(b_coast.hpdi.left.y))),border=NA,col='lightblue')
polygon(x=c(b_coast.hpdi.right.x,rev(b_coast.hpdi.right.x)),y=c(b_coast.hpdi.right.y,rep(0,length(b_coast.hpdi.right.y))),border=NA,col='lightpink')

polygon(x=c(b_coast.dens$x,rev(b_coast.dens$x)),y=c(b_coast.dens$y,rep(0,length(b_coast.dens$y))),border='lightgrey')
axis(1,at=seq(-1000,1000,200),labels=abs(seq(-1000,1000,200)))
axis(2)
box()
text(x=500,y=median(par('usr')[3:4]),label=paste('Changepoint after\n P=',round(sum(b_coast>0)/1000,2)),cex=0.8)
text(x=-500,y=median(par('usr')[3:4]),label=paste('Changepoint before\n P=',round(sum(b_coast<0)/1000,2)),cex=0.8)
abline(v=0,lty=2,lwd=2)

plot(a_inland.dens$x,a_inland.dens$y,type='n',xlab='Years',ylab='Probability Density',axes=FALSE,xlim=c(-1000,1000),main='Inland Changepoint vs Event A')
polygon(x=c(a_inland.hpdi.left.x,rev(a_inland.hpdi.left.x)),y=c(a_inland.hpdi.left.y,rep(0,length(a_inland.hpdi.left.y))),border=NA,col='lightblue')
polygon(x=c(a_inland.hpdi.right.x,rev(a_inland.hpdi.right.x)),y=c(a_inland.hpdi.right.y,rep(0,length(a_inland.hpdi.right.y))),border=NA,col='lightpink')


polygon(x=c(a_inland.dens$x,rev(a_inland.dens$x)),y=c(a_inland.dens$y,rep(0,length(a_inland.dens$y))),border='lightgrey')
axis(1,at=seq(-1000,1000,200),labels=abs(seq(-1000,1000,200)))
axis(2)
box()
text(x=500,y=median(par('usr')[3:4]),label=paste('Changepoint after\n P=',round(sum(a_inland>0)/1000,2)),cex=0.8)
text(x=-500,y=median(par('usr')[3:4]),label=paste('Changepoint before\n P=',round(sum(a_inland<0)/1000,2)),cex=0.8)
abline(v=0,lty=2,lwd=2)


plot(b_inland.dens$x,b_inland.dens$y,type='n',xlab='Years',ylab='Probability Density',axes=FALSE,xlim=c(-1000,1000),main='Inland Changepoint vs Event B')
polygon(x=c(b_inland.hpdi.left.x,rev(b_inland.hpdi.left.x)),y=c(b_inland.hpdi.left.y,rep(0,length(b_inland.hpdi.left.y))),border=NA,col='lightblue')
polygon(x=c(b_inland.hpdi.right.x,rev(b_inland.hpdi.right.x)),y=c(b_inland.hpdi.right.y,rep(0,length(b_inland.hpdi.right.y))),border=NA,col='lightpink')
polygon(x=c(b_inland.dens$x,rev(b_inland.dens$x)),y=c(b_inland.dens$y,rep(0,length(b_inland.dens$y))),border='lightgrey')
axis(1,at=seq(-1000,1000,200),labels=abs(seq(-1000,1000,200)))
axis(2)
box()
text(x=500,y=median(par('usr')[3:4]),label=paste('Changepoint after\n P=',round(sum(b_inland>0)/1000,2)),cex=0.8)
text(x=-500,y=median(par('usr')[3:4]),label=paste('Changepoint before\n P=',round(sum(b_inland<0)/1000,2)),cex=0.8)
abline(v=0,lty=2,lwd=2)
dev.off()

