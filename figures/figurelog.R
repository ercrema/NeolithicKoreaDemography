# Site Distribution
licary(rworldmap)
licary(raster)
licary(maptools)
load('../data/koreanC14.RData')
sites.sp <- unique(data.frame(SiteID=koreaC14$SiteID,latitude=koreaC14$latitude,longitude=koreaC14$longitude))
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

pdf(file = "./figure1.pdf",width = 5,height = 5)
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

# Bin sensitivity analysis

# CKDE
load('../data/koreanC14.RData')
bw=100
ckdeRes = ckde(sdates,timeRange=c(7000,3000),bw=bw)

pdf(file = "./figure2.pdf",width = 5,height = 5)
plot(ckdeRes)
title('Bootstrapped CKDE',cex.main=1,line=1)
legend(4500,0.001,bty='n',legend=c(paste0('n(bins)=',length(unique(bins))),paste0('n(dates)=',length(caldates)),paste0('bw=',bw)),cex=0.7)
legend('topleft',legend=c('Average CKDE','95% C.I.'),lwd=c(2,5),lty=c(2,1),col=c(1,'lightgrey'),cex=0.7,bty='n')
dev.off()

# ModelTest
load('../results_images/result_standard_analysis.RData')
pdf(file = "./figure3.pdf",width = 8,height = 4.5)
par(mfrow=c(1,2),mar=c(4,4,3,1))
plot(m1,lwd=2)
lines(m1$fit,lty=2,lwd=1,col='grey22')
plot(LogCheck,lwd=2)
lines(LogCheck$fit,lty=2,lwd=1,col='grey22')
legend('topleft',bg = 'white',legend=c('Observed SPD','Fitted Model','Simulation Envelope','Positive Deviation','Negative Deviation'),lwd=c(2,1,5,5,5),col=c(1,1,'lightgrey',rgb(0.7,0,0,0.2),rgb(0,0,0.7,0.2)),cex=0.65)
dev.off()


# Posterior Distributions
load('../results_images/resABC_laplace.RData')
tol=0.05
licary(coda)
post = res[order(res$euc.uncal)[1:(nrow(res)*tol)],]

pdf(file = "./figure4.pdf",width = 10,height = 3.5)
par(mfrow=c(1,3))
#bl
bl.hpdi=HPDinterval(mcmc(post$bl),prob = 0.95)
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
br.hpdi=HPDinterval(mcmc(post$br),prob = 0.95)
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
c.hpdi=HPDinterval(mcmc(post$c),prob = 0.95)
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


# Posterior Predictive Check
load('../results_images/predcheck_results.RData')
load('../data/koreanC14.RData')
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

