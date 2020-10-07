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
load('../R_image_files/kim2004_agedepthmodel.RData')
load('../R_image_files/resABC_laplace_general.RData')
load('../R_image_files/resABC_laplace_coastal.RData')
load('../R_image_files/resABC_laplace_inland.RData')
load('../R_image_files/predcheck_results_general.RData')
load('../R_image_files/predcheck_results_coastal.RData')
load('../R_image_files/predcheck_results_inland.RData')

## Figure 1 ####
pdf(file = "./figure1.pdf",width = 3.4,height = 2,pointsize=1.5)
par(mfrow=c(1,2),mar=c(4,3,3,1))

## Left Panel:

# Create SpatialPointsDataFrame
sites.sp <- unique(data.frame(SiteID=koreaC14$site_id,latitude=koreaC14$latitude,longitude=koreaC14$longitude,region=koreaC14$region))
coordinates(sites.sp) <- c("longitude","latitude")
proj4string(sites.sp) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

# Get Background
coast <- getMap(resolution = "high")

# Plot Sites
sites.sp$col = "#FC8D62"
sites.sp$col[sites.sp$region=='coastal'] = "#66C2A5"
plot(coast,col="grey78",xlim=extent(sites.sp)[1:2],
     ylim=extent(sites.sp)[3:4],border=NA,xlab="",ylab="")
abline(v=122:135,lwd=0.5,col='grey88')
abline(h=33:38,lwd=0.5,col='grey88')
points(sites.sp,pch=20,cex=0.7,col=sites.sp$col)
box()
axis(side=1,at=122:135,cex=0.7,las=2,cex.axis=0.6,hadj=0.3,tck=-0.01)
axis(side=2,at=33:38,cex=0.7,las=2,cex.axis=0.6,hadj=-0.3,tck=-0.01)
mtext('Longitude',1,2,cex=0.9)
mtext('Latitude',2,2,cex=0.9)
#legend('topright', legend=c('Coastal Sites','Inland Sites'), pch=20,col=c("#66C2A5","#FC8D62"),bg='white',cex=0.8)
title('a')
## Right Panel:

runm = 100 #smoothing window
timeRange = c(7000,3000)
combined.spd=stackspd(x=caldates,timeRange=timeRange,bins=bins,runm = runm,group = koreaC14$region)
plot(combined.spd,legend.arg=list(x='topright'))
title('b')
dev.off()


## Figure 2 ####
boomBust.inland=plot(coastal.inland.permtest,bbty='b',focalm='inland')
boomBust.coastal=plot(coastal.inland.permtest,bbty='b',focalm='coastal')

pdf(file = "./figure2.pdf",width = 3.4,height = 4,pointsize=1.5)
par(mfrow=c(2,1),mar=c(5,4,2,1.1))

plot(coastal.inland.permtest$observed$coastal$calBP,coastal.inland.permtest$observed$coastal$PrDens,type='n',xlab='Year cal BP',ylab='Summed Probability',xlim=c(7000,3000),ylim=c(0,0.18))

for (i in 1:length(boomBust.coastal$booms))
{
  indexCalBP = boomBust.coastal$booms[[i]][[2]]
  indexPrDens = which(coastal.inland.permtest$observed$coastal$calBP%in%indexCalBP)
  indexPrDens = c(indexPrDens,indexPrDens[length(indexPrDens)])
  polygon(c(indexCalBP,rev(indexCalBP)),c(coastal.inland.permtest$envelope$coastal[indexPrDens,2],rev(coastal.inland.permtest$observed$coastal$PrDens[indexPrDens])),border=NA,col=rgb(0.80,0.36,0.36,0.8))
}

for (i in 1:length(boomBust.coastal$busts))
{
  indexCalBP = boomBust.coastal$busts[[i]][[2]]
  indexPrDens = which(coastal.inland.permtest$observed$coastal$calBP%in%indexCalBP)
  indexPrDens = c(indexPrDens,indexPrDens[length(indexPrDens)])
  polygon(c(indexCalBP,rev(indexCalBP)),c(coastal.inland.permtest$envelope$coastal[indexPrDens,1],rev(coastal.inland.permtest$observed$coastal$PrDens[indexPrDens])),border=NA,col=rgb(0.25,0.41,0.88,0.8))
}
polygon(c(coastal.inland.permtest$observed$coastal$calBP,rev(coastal.inland.permtest$observed$coastal$calBP)),c(coastal.inland.permtest$envelope$coastal[,1],rev(coastal.inland.permtest$envelope$coastal[,2])),border=NA,col='lightgrey')
lines(coastal.inland.permtest$observed$coastal$calBP,coastal.inland.permtest$observed$coastal$PrDens,lwd=1.5)
text(x=6600,y=0.17,label='Coastal SPD')
text(x=6500,y=0.155,label=paste0('Global P-value<',round(coastal.inland.permtest$pValueList[1],5)),cex=0.8)



plot(coastal.inland.permtest$observed$inland$calBP,coastal.inland.permtest$observed$inland$PrDens,type='n',xlab='Year cal BP',ylab='Summed Probability',xlim=c(7000,3000),ylim=c(0,0.18))

for (i in 1:length(boomBust.coastal$booms))
{
  indexCalBP = boomBust.coastal$booms[[i]][[2]]
  indexPrDens = which(coastal.inland.permtest$observed$inland$calBP%in%indexCalBP)
  indexPrDens = c(indexPrDens,indexPrDens[length(indexPrDens)])
  polygon(c(indexCalBP,rev(indexCalBP)),c(coastal.inland.permtest$envelope$inland[indexPrDens,2],rev(coastal.inland.permtest$observed$inland$PrDens[indexPrDens])),border=NA,col=rgb(0.80,0.36,0.36,0.8))
}

for (i in 1:length(boomBust.coastal$busts))
{
  indexCalBP = boomBust.coastal$busts[[i]][[2]]
  indexPrDens = which(coastal.inland.permtest$observed$inland$calBP%in%indexCalBP)
  indexPrDens = c(indexPrDens,indexPrDens[length(indexPrDens)])
  polygon(c(indexCalBP,rev(indexCalBP)),c(coastal.inland.permtest$envelope$inland[indexPrDens,1],rev(coastal.inland.permtest$observed$inland$PrDens[indexPrDens])),border=NA,col=rgb(0.25,0.41,0.88,0.8))
}
polygon(c(coastal.inland.permtest$observed$inland$calBP,rev(coastal.inland.permtest$observed$inland$calBP)),c(coastal.inland.permtest$envelope$inland[,1],rev(coastal.inland.permtest$envelope$inland[,2])),border=NA,col='lightgrey')
lines(coastal.inland.permtest$observed$inland$calBP,coastal.inland.permtest$observed$inland$PrDens,lwd=1.5)
text(x=6600,y=0.17,label='Inland SPD')
text(x=6500,y=0.155,label=paste0('Global P-value<',round(coastal.inland.permtest$pValueList[2],5)),cex=0.8)
legend('topright',legend=c('Observed SPD','Null SPD','Positive Deviation','Negative Deviation'),lwd=c(1,5,5,5),col=c(1,'lightgrey',rgb(0.80,0.36,0.36,0.8),rgb(0.25,0.41,0.88,0.8)),cex=0.8,bty='n')
dev.off()












## Figure 3 ####
tol=0.01
library(coda)
post.coastal = abc.coastal[order(abc.coastal$euc.uncal)[1:(nrow(abc.coastal)*tol)],]
post.inland = abc.inland[order(abc.inland$euc.uncal)[1:(nrow(abc.inland)*tol)],]
coastal.col.alpha = rgb(0.4,0.76,0.65,0.5)
inland.col.alpha = rgb(0.98,0.55,0.38,0.5)
coastal.col="#66C2A5"
inland.col="#FC8D62"

options(scipen = 9999)
pdf(file = "./figure3.pdf",width = 3.4,height = 8)
par(mfrow=c(3,1))

#bl
bl.hpdi.coastal=HPDinterval(mcmc(post.coastal$bl),prob = 0.90)
bl.hpdi.inland=HPDinterval(mcmc(post.inland$bl),prob = 0.90)

d.bl.coastal=density(post.coastal$bl)
d.bl.inland=density(post.inland$bl,bw=d.bl.coastal$bw)

plot(0,0,type='n',xlab='% Annual Growth Rate',ylab='Probability Density',axes=FALSE,xlim=range(c(d.bl.coastal$x,d.bl.inland$x)),ylim=range(c(d.bl.coastal$y,d.bl.inland$y)))

title(TeX('$r_1$ Posterior'))
axis(1,at=axTicks(1),labels=axTicks(1)*100)
axis(2)
hpdi.x.coastal = d.bl.coastal$x[which(d.bl.coastal$x>=bl.hpdi.coastal[1]&d.bl.coastal$x<=bl.hpdi.coastal[2])]
hpdi.x.inland = d.bl.inland$x[which(d.bl.inland$x>=bl.hpdi.inland[1]&d.bl.inland$x<=bl.hpdi.inland[2])]
hpdi.y.coastal = d.bl.coastal$y[which(d.bl.coastal$x>=bl.hpdi.coastal[1]&d.bl.coastal$x<=bl.hpdi.coastal[2])]
hpdi.y.inland = d.bl.inland$y[which(d.bl.inland$x>=bl.hpdi.inland[1]&d.bl.inland$x<=bl.hpdi.inland[2])]

polygon(x=c(hpdi.x.coastal,rev(hpdi.x.coastal)),y=c(hpdi.y.coastal,rep(0,length(hpdi.y.coastal))),border=NA,col=coastal.col.alpha)
polygon(x=c(d.bl.coastal$x,rev(d.bl.coastal$x)),y=c(d.bl.coastal$y,rep(0,length(d.bl.coastal$y))),border=coastal.col)

polygon(x=c(hpdi.x.inland,rev(hpdi.x.inland)),y=c(hpdi.y.inland,rep(0,length(hpdi.y.inland))),border=NA,col=inland.col.alpha)
polygon(x=c(d.bl.inland$x,rev(d.bl.inland$x)),y=c(d.bl.inland$y,rep(0,length(d.bl.inland$y))),border=inland.col)

abline(v=median(post.coastal$bl),lty=2,col=coastal.col)
abline(v=median(post.inland$bl),lty=2,col=inland.col)
legend('topright',legend=c('Coastal','Inland'),fill=c(coastal.col,inland.col))


#br
br.hpdi.coastal=HPDinterval(mcmc(post.coastal$br),prob = 0.90)
br.hpdi.inland=HPDinterval(mcmc(post.inland$br),prob = 0.90)

d.br.coastal=density(post.coastal$br)
d.br.inland=density(post.inland$br,bw=d.br.coastal$bw)

plot(0,0,type='n',xlab='% Annual Growth Rate',ylab='Probability Density',axes=FALSE,xlim=range(c(d.br.coastal$x,d.br.inland$x)),ylim=range(c(d.br.coastal$y,d.br.inland$y)))

title(TeX('$r_2$ Posterior'))
axis(1,at=axTicks(1),labels=-axTicks(1)*100)
axis(2)
hpdi.x.coastal = d.br.coastal$x[which(d.br.coastal$x>=br.hpdi.coastal[1]&d.br.coastal$x<=br.hpdi.coastal[2])]
hpdi.x.inland = d.br.inland$x[which(d.br.inland$x>=br.hpdi.inland[1]&d.br.inland$x<=br.hpdi.inland[2])]
hpdi.y.coastal = d.br.coastal$y[which(d.br.coastal$x>=br.hpdi.coastal[1]&d.br.coastal$x<=br.hpdi.coastal[2])]
hpdi.y.inland = d.br.inland$y[which(d.br.inland$x>=br.hpdi.inland[1]&d.br.inland$x<=br.hpdi.inland[2])]

polygon(x=c(hpdi.x.coastal,rev(hpdi.x.coastal)),y=c(hpdi.y.coastal,rep(0,length(hpdi.y.coastal))),border=NA,col=coastal.col.alpha)
polygon(x=c(d.br.coastal$x,rev(d.br.coastal$x)),y=c(d.br.coastal$y,rep(0,length(d.br.coastal$y))),border=coastal.col)

polygon(x=c(hpdi.x.inland,rev(hpdi.x.inland)),y=c(hpdi.y.inland,rep(0,length(hpdi.y.inland))),border=NA,col=inland.col.alpha)
polygon(x=c(d.br.inland$x,rev(d.br.inland$x)),y=c(d.br.inland$y,rep(0,length(d.br.inland$y))),border=inland.col)

abline(v=median(post.coastal$br),lty=2,col=coastal.col)
abline(v=median(post.inland$br),lty=2,col=inland.col)


#c
c.hpdi.coastal=HPDinterval(mcmc(post.coastal$c),prob = 0.90)
c.hpdi.inland=HPDinterval(mcmc(post.inland$c),prob = 0.90)

d.c.coastal=density(post.coastal$c)
d.c.inland=density(post.inland$c,bw=d.c.coastal$bw)


plot(0,0,type='n',xlab='Cal BP',ylab='Probability Density',axes=FALSE,xlim=rev(range(c(d.c.coastal$x,d.c.inland$x))),ylim=range(c(d.c.coastal$y,d.c.inland$y)))



title(TeX('$c$ Posterior'))
axis(1)
axis(2)

hpdi.x.coastal = d.c.coastal$x[which(d.c.coastal$x>=c.hpdi.coastal[1]&d.c.coastal$x<=c.hpdi.coastal[2])]
hpdi.x.inland = d.c.inland$x[which(d.c.inland$x>=c.hpdi.inland[1]&d.c.inland$x<=c.hpdi.inland[2])]
hpdi.y.coastal = d.c.coastal$y[which(d.c.coastal$x>=c.hpdi.coastal[1]&d.c.coastal$x<=c.hpdi.coastal[2])]
hpdi.y.inland = d.c.inland$y[which(d.c.inland$x>=c.hpdi.inland[1]&d.c.inland$x<=c.hpdi.inland[2])]

polygon(x=c(hpdi.x.coastal,rev(hpdi.x.coastal)),y=c(hpdi.y.coastal,rep(0,length(hpdi.y.coastal))),border=NA,col=coastal.col.alpha)
polygon(x=c(d.c.coastal$x,rev(d.c.coastal$x)),y=c(d.c.coastal$y,rep(0,length(d.c.coastal$y))),border=coastal.col)

polygon(x=c(hpdi.x.inland,rev(hpdi.x.inland)),y=c(hpdi.y.inland,rep(0,length(hpdi.y.inland))),border=NA,col=inland.col.alpha)
polygon(x=c(d.c.inland$x,rev(d.c.inland$x)),y=c(d.c.inland$y,rep(0,length(d.c.inland$y))),border=inland.col)

abline(v=median(post.coastal$c),lty=2,col=coastal.col)
abline(v=median(post.inland$c),lty=2,col=inland.col)
dev.off()


## Figure 4 ####
point_a=kim2004.model.marine20$thetaPredict[,81]
point_b=kim2004.model.marine20$thetaPredict[,80]

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



pdf(file = "./figure4.pdf",width=3.4,height = 3.4,pointsize=1.5)
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
abline(v=0,lty=2,lwd=1.5)


plot(b_coast.dens$x,b_coast.dens$y,type='n',xlab='Years',ylab='Probability Density',axes=FALSE,xlim=c(-1000,1000),main='Coastal Changepoint vs Event B')
polygon(x=c(b_coast.hpdi.left.x,rev(b_coast.hpdi.left.x)),y=c(b_coast.hpdi.left.y,rep(0,length(b_coast.hpdi.left.y))),border=NA,col='lightblue')
polygon(x=c(b_coast.hpdi.right.x,rev(b_coast.hpdi.right.x)),y=c(b_coast.hpdi.right.y,rep(0,length(b_coast.hpdi.right.y))),border=NA,col='lightpink')

polygon(x=c(b_coast.dens$x,rev(b_coast.dens$x)),y=c(b_coast.dens$y,rep(0,length(b_coast.dens$y))),border='lightgrey')
axis(1,at=seq(-1000,1000,200),labels=abs(seq(-1000,1000,200)))
axis(2)
box()
text(x=500,y=median(par('usr')[3:4]),label=paste('Changepoint after\n P=',round(sum(b_coast>0)/1000,2)),cex=0.8)
text(x=-500,y=median(par('usr')[3:4]),label=paste('Changepoint before\n P=',round(sum(b_coast<0)/1000,2)),cex=0.8)
abline(v=0,lty=2,lwd=1.5)

plot(a_inland.dens$x,a_inland.dens$y,type='n',xlab='Years',ylab='Probability Density',axes=FALSE,xlim=c(-1000,1000),main='Inland Changepoint vs Event A')
polygon(x=c(a_inland.hpdi.left.x,rev(a_inland.hpdi.left.x)),y=c(a_inland.hpdi.left.y,rep(0,length(a_inland.hpdi.left.y))),border=NA,col='lightblue')
polygon(x=c(a_inland.hpdi.right.x,rev(a_inland.hpdi.right.x)),y=c(a_inland.hpdi.right.y,rep(0,length(a_inland.hpdi.right.y))),border=NA,col='lightpink')


polygon(x=c(a_inland.dens$x,rev(a_inland.dens$x)),y=c(a_inland.dens$y,rep(0,length(a_inland.dens$y))),border='lightgrey')
axis(1,at=seq(-1000,1000,200),labels=abs(seq(-1000,1000,200)))
axis(2)
box()
text(x=500,y=median(par('usr')[3:4]),label=paste('Changepoint after\n P=',round(sum(a_inland>0)/1000,2)),cex=0.8)
text(x=-500,y=median(par('usr')[3:4]),label=paste('Changepoint before\n P=',round(sum(a_inland<0)/1000,2)),cex=0.8)
abline(v=0,lty=2,lwd=1.5)


plot(b_inland.dens$x,b_inland.dens$y,type='n',xlab='Years',ylab='Probability Density',axes=FALSE,xlim=c(-1000,1000),main='Inland Changepoint vs Event B')
polygon(x=c(b_inland.hpdi.left.x,rev(b_inland.hpdi.left.x)),y=c(b_inland.hpdi.left.y,rep(0,length(b_inland.hpdi.left.y))),border=NA,col='lightblue')
polygon(x=c(b_inland.hpdi.right.x,rev(b_inland.hpdi.right.x)),y=c(b_inland.hpdi.right.y,rep(0,length(b_inland.hpdi.right.y))),border=NA,col='lightpink')
polygon(x=c(b_inland.dens$x,rev(b_inland.dens$x)),y=c(b_inland.dens$y,rep(0,length(b_inland.dens$y))),border='lightgrey')
axis(1,at=seq(-1000,1000,200),labels=abs(seq(-1000,1000,200)))
axis(2)
box()
text(x=500,y=median(par('usr')[3:4]),label=paste('Changepoint after\n P=',round(sum(b_inland>0)/1000,2)),cex=0.8)
text(x=-500,y=median(par('usr')[3:4]),label=paste('Changepoint before\n P=',round(sum(b_inland<0)/1000,2)),cex=0.8)
abline(v=0,lty=2,lwd=1.5)
dev.off()