# Load Libraries, Data, and Results ####
library(raster)
library(maptools)
library(rworldmap)
library(rcarbon)
library(Bchron)
library(coda)


load('../data/koreanC14.RData')
load('../results_images/test_results.RData')

# Site Distribution Figure ####
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

# Posterior Distributions ####
load('../results_images/resABC_laplace.RData')
tol=0.05
library(coda)
post = res[order(res$euc.uncal)[1:(nrow(res)*tol)],]
options(scipen = 9999)
pdf(file = "./figure4.pdf",width = 10,height = 3.5)
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



# Posterior Predictive Check ####
load('../results_images/predcheck_results.RData')
load('../data/koreanC14.RData')
library(rcarbon)
observed = spd(caldates,bins,timeRange=c(7000,3000),spdnormalised = TRUE)
ppmedian=apply(ppcheck.uncal,1,median)
pplo=apply(ppcheck.uncal,1,quantile,0.025)
pphi=apply(ppcheck.uncal,1,quantile,0.975)



pdf(file = "./figure5.pdf",width = 4,height = 4)
plot(observed$grid$calBP,observed$grid$PrDens,type='n',xlim=c(7000,3000),ylim=c(0,max(c(observed$grid$PrDens,pphi))),xlab='cal BP',ylab='Summed Probability')
polygon(c(7000:3000,rev(7000:3000)),c(pplo,rev(pphi)),border=NA,col='lightgrey')
lines(7000:3000,ppmedian,col=2,lty=2)
lines(observed$grid$calBP,observed$grid$PrDens,lwd=1)
legend('topleft',legend=c('Observed','Median Posterior Predictive Check','95% Posterior Predictive Interval'),col=c(1,2,'lightgrey'),lwd=c(1,1,5),lty=c(1,2,1),cex=0.6,bg='white')
dev.off()



#  



# Age-Depth Model of Kim 2004 ####

# Read Kim 2004
kim2004.dates = read.csv('../data/SSDP_102.Kim.2004-chron.csv',skip=1)
kim2004.temp  = read.csv('../data/SSDP_102.Kim.2004.csv',skip=1)

# Read marine20 calibration curve
marine20 = read.csv('http://intcal.org/curves/marine20.14c', encoding="UTF-8",skip=11,header=F)
createCalCurve(name='marine20',calAges=marine20[,1],uncalAges=marine20[,2],oneSigma=marine20[,3])
file.copy(from = 'marine20.rda',to = system.file('data',package='Bchron'))

# Fit Compound Poisson-Gamma chronology model
# kim2004.model.marine13 = Bchronology(ages=round(kim2004.dates$T2L_SSDP_102_c14_date - kim2004.dates$T2L_SSDP_102_delta_r),ageSds=round(sqrt(kim2004.dates$T2L_SSDP_102_c14_1s_err^2+kim2004.dates$T2L_SSDP_102_delta_r_1s_error^2)),calCurves = rep('marine13',nrow(kim2004.dates)),ids=kim2004.dates$T2L_SSDP_102_labcode,positions=kim2004.dates$T2L_SSDP_102_depth_top,predictPositions=kim2004.temp$T2L_SSDP_102_depth)

kim2004.model.marine20 = Bchronology(ages=round(kim2004.dates$T2L_SSDP_102_c14_date - kim2004.dates$T2L_SSDP_102_delta_r),ageSds=round(sqrt(kim2004.dates$T2L_SSDP_102_c14_1s_err^2+kim2004.dates$T2L_SSDP_102_delta_r_1s_error^2)),calCurves = rep('marine20',nrow(kim2004.dates)),ids=kim2004.dates$T2L_SSDP_102_labcode,positions=kim2004.dates$T2L_SSDP_102_depth_top,predictPositions=kim2004.temp$T2L_SSDP_102_depth)

# Compare predicted ages
# med.bchron.marine13 = apply(kim2004.model.marine13$thetaPredict,2,median)
# med.bchron.marine20 = apply(kim2004.model.marine20$thetaPredict,2,median)
# par(mfrow=c(2,2))
# plot(kim2004.temp$T2L_SSDP_102_age_medianMedianBacon,med.bchron.marine13,pch=20)
# abline(a=0,b=1,col='red')
# plot(kim2004.temp$T2L_SSDP_102_age_medianMedianBacon,med.bchron.marine20,pch=20)
# abline(a=0,b=1,col='red')

# Plot Median Predicted Temperature Change
pdf(file = "./figure_kim2004_reanalysis.pdf",width = 6,height = 5.5)
layout(matrix(c(1,2,1,3),2,2))
par(mar=c(5,4,1,1))
plot(3000:7000,med,xlim=c(7000,3000),type='n',ylim=c(range(kim2004.temp$T2L_SSDP_102_uk37_SST_from_uk37)),xlab='cal BP',ylab='temperature (deg C)')
rect(xleft=med.bchron.marine20[80],xright=med.bchron.marine20[81],ybottom=10,ytop=30,border=NA,col='lightblue')
lines(med.bchron.marine20,kim2004.temp$T2L_SSDP_102_uk37_SST_from_uk37,type='l',lty=1,col='darkgrey')
points(med.bchron.marine20[80:81],kim2004.temp$T2L_SSDP_102_uk37_SST_from_uk37[80:81],pch=20,col=1)
text(4949.801,21.89769,'a')
text(4716.733,20.45735,'b')
lines(med.bchron.marine20[80:81],kim2004.temp$T2L_SSDP_102_uk37_SST_from_uk37[80:81],lwd=2,lty=2,col=1)

point_a=kim2004.model.marine20$thetaPredict[,81]
point_b=kim2004.model.marine20$thetaPredict[,80]

dcool.hpdi.a=HPDinterval(mcmc(point_a),prob = 0.90)
d.event.a=density(point_a)
dcool.hpdi.b=HPDinterval(mcmc(point_b),prob = 0.90)
d.event.b=density(point_b)

plot(d.event.a$x,d.event.a$y,type='n',xlab='cal BP',ylab='Probability Density',axes=FALSE,xlim=c(5500,4000))
legend('topright',legend='Timing of a',bty='n',cex=0.8)
axis(1)
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
hpdi.x = d.event.b$x[which(d.event.b$x>=dcool.hpdi.b[1]&d.event.b$x<=dcool.hpdi.b[2])]
hpdi.y = d.event.b$y[which(d.event.b$x>=dcool.hpdi.b[1]&d.event.b$x<=dcool.hpdi.b[2])]
polygon(x=c(hpdi.x,rev(hpdi.x)),y=c(hpdi.y,rep(0,length(hpdi.y))),border=NA,col='lightblue')
polygon(x=c(d.event.b$x,rev(d.event.b$x)),y=c(d.event.b$y,rep(0,length(d.event.b$y))))
abline(v=median(point_b),lty=2)
dev.off()