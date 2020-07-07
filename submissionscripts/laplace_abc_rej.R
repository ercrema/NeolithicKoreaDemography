library(rcarbon)
library(Bchron)
library(foreach)
library(doParallel)
library(doSNOW)

### Load Observed Data ###
load('../data/koreanC14.RData')

### Core Simulation Model ###
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

### Generate Parameters ####
nsim=100000
a=runif(nsim,0.0001,0.001)
bl=1/rexp(nsim,500)
br=1/rexp(nsim,500)
c=runif(nsim,4000,6000)

### ABC Settings ###
tol=0.01
ncores = 3
cl <- makeCluster(ncores)
registerDoSNOW(cl)
pb <- txtProgressBar(max = nsim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

reslist <- foreach (i=1:nsim,.packages=c('rcarbon','Bchron'),.options.snow = opts) %dopar% {
  res=sim.model(x=koreaC14,caldates=caldates,bins=bins,a=a[i],bl=bl[i],br=br[i],c=c[i],timeRange=c(7000,3000))
  return(res)
}


a = unlist(lapply(reslist,function(x){x[[1]]}))
bl = 1/unlist(lapply(reslist,function(x){x[[2]]}))
br = 1/unlist(lapply(reslist,function(x){x[[3]]}))
c = unlist(lapply(reslist,function(x){x[[4]]}))
euc.uncal = unlist(lapply(reslist,function(x){x[[5]]}))
euc.cal = unlist(lapply(reslist,function(x){x[[6]]}))
ks.uncal = unlist(lapply(reslist,function(x){x[[7]]}))
ks.cal = unlist(lapply(reslist,function(x){x[[8]]}))

res=data.frame(a=a,bl=bl,br=br,c=c,euc.uncal=euc.uncal,euc.cal=euc.cal,ks.uncal=ks.uncal,ks.cal=ks.cal)

save(res,file='resABC_laplace.RData')

