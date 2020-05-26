
library(rworldmap)
library(raster)
library(maptools)
tmp<-read.csv("~/Dropbox/PrehistoricPopulationKorea/R/data/koreaC14dates.csv")
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

plot(coast,col="lightgrey",xlim=extent(hs)[1:2],ylim=extent(hs)[3:4],border=NA,xlab="Longitude",ylab="Latitude")
plot(hs, col=grey(0:100/100), legend=FALSE,add=TRUE)
points(sites.sp,pch=20,col=rgb(1,0,0,0.6)))
box()
axis(side=1)
axis(side=2)

dev.print(device=pdf,"~/Dropbox/PrehistoricPopulationKorea/R/figures/figure1.pdf")

load("results.RData")
library(rcarbon)
library(extrafont)
loadfonts()

# Figure 2 (Observed SPD with BootStrap Interval)

spdKorea.norm <- spd(x=caldates,timeRange=c(7000,3000),bins=bins,runm=200,spdnormalised=T)
spdKoreaBase.norm <- spd(x=caldates,timeRange=c(7000,3000),bins=bins,spdnormalised=T)

par(family="Arial",bg="white")
plot(1,1,xlim=c(7000,3000),ylim=c(0,max(c(spdKoreaBase.norm$grid$PrDens,spdKorea.norm$grid$PrDens,bootSPD$bootmatrix))),type="n",ylab="",xlab="cal BP")
polygon(x=c(10000,0,0,10000),y=c(-100,-100,100,100),col="lightgrey",border=NA)

polygon(c(7000:3000,3000:7000),c(apply(bootBins$bootmatrix,1,quantile,0.025),rev(apply(bootBins$bootmatrix,1,quantile,0.975))),border=NA,col="deepskyblue")
lines(spdKoreaBase.norm$grid$calBP,spdKoreaBase.norm$grid$PrDens,xlim=c(7000,3000),lwd=0.5,col="blue")
lines(spdKorea.norm$grid$calBP,spdKorea.norm$grid$PrDens,col="darkblue",lwd=2)
abline(v=seq(7000,3000,-500),lty=3,col="white")
abline(h=c(0,0.0002,0.0004,0.0006,0.0008),lty=3,col="white")
axis(side=1,at=seq(7000,3000,-500),labels=NA,tck=-0.01)
legend("topleft",legend=c("SPD (raw)","SPD (smoothed)","Bootstrap CI"),col=c("blue","darkblue","deepskyblue"),lwd=c(0.5,2,8),lty=1,bg="white",cex=0.9, inset=0.01, box.lty=0)

dev.print(device=pdf,"~/Dropbox/PrehistoricPopulationKorea/R/figures/figure2.pdf")

# Figure 3 (Three Panel Plot with Model Test)
bbEXP=plot(ExpCheck,bbty="b")
bbLOG=plot(LogCheck,bbty="b")
bbLAP=plot(LapCheck,bbty="b")



par(family="Arial",bg="white",mfrow=c(1,3))

plot(1,1,xlim=c(7000,3000),ylim=c(0,0.2600765),type="n",ylab="",xlab="cal BP",main="")
polygon(x=c(10000,0,0,10000),y=c(-100,-100,100,100),col="lightgrey",border=NA)

polygon(c(7000:3000,3000:7000),c(ExpCheck$result$hi,rev(ExpCheck$result$lo)),border=NA,col="lightpink")
abline(v=seq(7000,3000,-500),lty=3,col="white")
abline(h=seq(0,0.2,0.05),lty=3,col="white")
lines(expFitDens$calBP,expFitDens$PrDens,lty=2,lwd=1.5,col="red")
lines(spdKorea$grid$calBP,spdKorea$grid$PrDens,col="darkblue",lwd=2)
#abline(v=seq(7000,3000,-500),lty=3,col="white")
#abline(h=seq(0,0.2,0.05),lty=3,col="white")
axis(side=1,at=seq(7000,3000,-500),labels=NA,tck=-0.01)
bbpolygons(bbEXP)
text(x=6000,y=0.25,"Exponential (p-value<0.001)",cex=1.4)



plot(1,1,xlim=c(7000,3000),ylim=c(0,0.2600765),type="n",ylab="",xlab="cal BP",main="")
polygon(x=c(10000,0,0,10000),y=c(-100,-100,100,100),col="lightgrey",border=NA)

polygon(c(7000:3000,3000:7000),c(LogCheck$result$hi,rev(LogCheck$result$lo)),border=NA,col="lightpink")
abline(v=seq(7000,3000,-500),lty=3,col="white")
abline(h=seq(0,0.2,0.05),lty=3,col="white")
lines(logFitDens$calBP,logFitDens$PrDens,lty=2,lwd=1.5,col="red")
lines(spdKorea$grid$calBP,spdKorea$grid$PrDens,col="darkblue",lwd=2)
#abline(v=seq(7000,3000,-500),lty=3,col="white")
#abline(h=seq(0,0.2,0.05),lty=3,col="white")
axis(side=1,at=seq(7000,3000,-500),labels=NA,tck=-0.01)
bbpolygons(bbLOG)
text(x=6000,y=0.25,"Logistic (p-value<0.001)",cex=1.4)


plot(1,1,xlim=c(7000,3000),ylim=c(0,0.2600765),type="n",ylab="",xlab="cal BP",main="")
polygon(x=c(10000,0,0,10000),y=c(-100,-100,100,100),col="lightgrey",border=NA)

polygon(c(7000:3000,3000:7000),c(LapCheck$result$hi,rev(LapCheck$result$lo)),border=NA,col="lightpink")
abline(v=seq(7000,3000,-500),lty=3,col="white")
abline(h=seq(0,0.2,0.05),lty=3,col="white")
lines(lapFitDens$calBP,lapFitDens$PrDens,lty=2,lwd=1.5,col="red")
lines(spdKorea$grid$calBP,spdKorea$grid$PrDens,col="darkblue",lwd=2)
#abline(v=seq(7000,3000,-500),lty=3,col="white")
#abline(h=seq(0,0.2,0.05),lty=3,col="white")
axis(side=1,at=seq(7000,3000,-500),labels=NA,tck=-0.01)
bbpolygons(bbLAP)
text(x=6000,y=0.25,"Laplace (p-value=0.2421)",cex=1.4)

legend("topright",legend=c("Observed SPD (smoothed)","Fitted Model", "Simulation Envelope", "Positive Deviation","Negative Deviation"),col=c("blue","red","lightpink",rgb(0.7,0,0,0.2),rgb(0,0,0.7,0.2)),lty=c(1,2,1,NA,NA),bg="white",lwd=c(2,1.5,6,NA,NA),pch=c(NA,NA,NA,15,15),pt.cex=2,cex=0.9, inset=0.01, box.lty=0)


dev.print(device=pdf,"~/Dropbox/PrehistoricPopulationKorea/R/figures/figure3.pdf")

