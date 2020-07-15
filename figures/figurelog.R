# Site Distribution
library(rworldmap)
library(raster)
library(maptools)
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





# Posterior Distributions

# Posterior Predictive Check