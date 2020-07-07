library(rcarbon)
library(sp)
library(dbscan)
library(Bchron)

### Read 14C Data ####
koreaC14<-read.csv("./data/koreaC14dates.csv")

#### Use DBSCAN to Clusters site in proximity #### 
sites <- unique(data.frame(SiteID=koreaC14$SiteID,latitude=koreaC14$latitude,longitude=koreaC14$longitude))
sites.sp = sites
sp::coordinates(sites.sp) <- c("longitude","latitude")
sp::proj4string(sites.sp) <- sp::CRS("+proj=longlat +datum=WGS84")
D <- sp::spDists(sites.sp,sites.sp,longlat=TRUE)
eps=1 #eps set to 1 km
clust<-dbscan(as.dist(D),eps=eps,minPts=1)$clust
sites$cluster<-clust
clusterNames=unique(clust)
sites$clat=NA
sites$clon=NA

for (x in 1:length(clusterNames))
{
  i <- which(sites$cluster==clusterNames[x])
  sites$clat[i] = mean(sites$latitude[i])
  sites$clon[i] = mean(sites$longitude[i])
}
sites <- data.frame(SiteID=sites$SiteID,cluster=sites$cluster,clon=sites$clon,clat=sites$clat)

koreaC14 = unique(merge(x=koreaC14,y=sites,by.x="SiteID",by.y="SiteID",x.all=TRUE))

#### rcarbon analysis ####
caldates <- calibrate(koreaC14$c14age,koreaC14$c14error)
bins=binPrep(ages=caldates,h=100,sites=koreaC14$cluster)

save(caldates,bins,koreaC14,file='d')



# CKDE
s = sampleDates(caldates,bins=bins,nsim=100,boot=TRUE)
ckdeNorm = ckde(s,timeRange=c(7000,3000),bw=100)
plot(ckdeNorm,type='multiline')

# ABC Routine ####

# step 1: random thinning and creating of target SPD
index=thinDates(koreaC14$c14age,koreaC14$c14error,bins,size=1,thresh=0,method='random')
target.spd=spd(caldates[index],timeRange = c(7000,3000),spdnormalised = TRUE)

# step 2: prior predictive check
x <- target.spd$grid$calBP
y <- target.spd$grid$PrDens

laplace <- function(x, a, bl, br, c) {
  ifelse(x < c, a * exp((x-c) / br), a * exp(-(x-c) / bl))
}
nsim=100
a=runif(nsim,0.0001,0.001)
bl=1/rexp(nsim,500)
br=1/rexp(nsim,500)
c=rnorm(nsim,mean=5000,sd=200)

prior.predict=matrix(NA,nrow=length(x),ncol=nsim)

for (i in 1:nsim)
{
  prior.predict[,i]=laplace(x=x,a=a[i],bl=bl[i],br=br[i],c=c[i])
}
prior.predict = apply(prior.predict,2,function(x){x/sum(x)})
plot(x,y,xlim=c(7000,3000),ylim=range(prior.predict),type='l',lwd=2,col=2)
apply(prior.predict,2,lines,x=x,col=rgb(0,0,0,0.1))
lines(x,prior.predict[,sample(100,size=1)],col='green')


# step 3: construct ABC model
sim.model = function(x,caldates,bins,a,bl,br,c,timeRange=c(7000,3000))
{
  # Generate Random Target
  index=thinDates(x$c14age,x$c14error,bins,size=1,thresh=0,method='random')
  target.spd =spd(caldates[index],timeRange = timeRange,datenormalised=TRUE,spdnormalised = TRUE,verbose=FALSE)
  trange = target.spd$grid$calBP
  tprob = target.spd$grid$PrDens

  # Create Model
  CalBP = trange
  PrDens=ifelse(CalBP < c, a * exp((CalBP-c) / br), a * exp(-(CalBP-c) / bl))
  d = data.frame(CalBP=CalBP,PrDens=PrDens/sum(PrDens))
  # Collect Samples
  class(d)='CalGrid'
  d = uncalibrate(d,verbose=F)
  n = length(index)
  sampled.c14.dates.uncal=sample(d$CRA,size=n,replace=T,prob=d$PrDens)
  sampled.c14.dates.cal=sample(d$CRA,size=n,replace=T,prob=d$Raw)
  # Calibrate and generate theoretical spds
  sampled.calibrated.dates.uncal = BchronCalibrate(ages=sampled.c14.dates.uncal,ageSds = sample(x$c14error[index],replace=TRUE,size=n),calCurves=rep('intcal13',n))
  sampled.calibrated.dates.uncal = as.CalDates(sampled.calibrated.dates.uncal)
  sampled.calibrated.dates.uncal$metadata$StartBP=50000
  sampled.calibrated.dates.uncal$metadata$EndBP=0
  sampled.spd.uncal = spd(sampled.calibrated.dates.uncal,timeRange=timeRange,spdnormalised = TRUE,verbose=FALSE)
  candidate.uncal = sampled.spd.uncal$grid$PrDens
  sampled.calibrated.dates.cal = BchronCalibrate(ages=sampled.c14.dates.cal,ageSds = sample(x$c14error[index],replace=TRUE,size=n),calCurves=rep('intcal13',n))
  sampled.calibrated.dates.cal = as.CalDates(sampled.calibrated.dates.cal)
  sampled.calibrated.dates.cal$metadata$StartBP=50000
  sampled.calibrated.dates.cal$metadata$EndBP=0
  sampled.spd.cal = spd(sampled.calibrated.dates.cal,timeRange=timeRange,spdnormalised = TRUE,verbose=FALSE)
  candidate.cal = sampled.spd.cal$grid$PrDens
  
  # Compute dissimilarity between target and simulation
  euc_epsilon_uncalsample=sqrt(sum((candidate.uncal-tprob)^2))
  euc_epsilon_calsample=sqrt(sum((candidate.cal-tprob)^2))
  ks_epsilon_uncalsample=max(cumsum(candidate.uncal)-cumsum(tprob))
  ks_epsilon_calsample=max(cumsum(candidate.cal)-cumsum(tprob))
  
  return(list(a=a,bl=bl,br=br,c=c,euc_epsilon_uncalsample=euc_epsilon_uncalsample,euc_epsilon_calsample=euc_epsilon_calsample,ks_epsilon_uncalsample=ks_epsilon_uncalsample,ks_epsilon_calsample=ks_epsilon_calsample))
}

# test run:
nsim=5
a=runif(nsim,0.0001,0.001)
bl=1/rexp(nsim,500)
br=1/rexp(nsim,500)
c=rnorm(nsim,mean=5000,sd=200)
tmp.res = vector('list',length=nsim)
for (i in 1:nsim)
{
  print(i)
  tmp.res[[i]]=sim.model(x=koreaC14,caldates=caldates,bins=bins,a=a[i],bl=bl[i],br=br[i],c=c[i],timeRange=c(7000,3000))
}



