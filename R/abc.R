### Generate Sample Data with Known Growth Rate ####

# Load Required R packages
library(rcarbon)
library(Bchron)
library(foreach)
library(doParallel)
library(doSNOW)

# Set up population parameters
start = 7000
end = 3000
r = -0.002
CalBP = 7000:3000
PrDens = (1+r)^CalBP / sum((1+r)^CalBP)

# Generate hypothetical sample
n = 200 #sample size
set.seed(1)
sampled.calendar.dates = sample(CalBP,size=n,replace=TRUE,prob=PrDens)
sampled.c14.dates.backcalibrate = uncalibrate(sampled.calendar.dates)$rCRA
d = data.frame(CalBP=CalBP,PrDens=PrDens)
class(d)='CalGrid'
d = uncalibrate(d,verbose=F)
sampled.c14.dates.uncalsample=sample(d$CRA,size=n,replace=T,prob=d$PrDens) #Andy's assumption
sampled.c14.dates.calsample=sample(d$CRA,size=n,replace=T,prob=d$Raw) #traditional assumption


sampled.calibrated.dates.uncalsample = calibrate(sampled.c14.dates.uncalsample,errors=rep(20,n)) 
sampled.calibrated.dates.calsample = calibrate(sampled.c14.dates.calsample,errors=rep(20,n)) 

sampled.spd.uncalsample = spd(sampled.calibrated.dates.uncalsample,timeRange=c(7000,3000),spdnormalised = TRUE)
sampled.spd.calsample = spd(sampled.calibrated.dates.calsample,timeRange=c(7000,3000),spdnormalised = TRUE)

observed.uncalsample = sampled.spd.uncalsample$grid$PrDens
observed.calsample = sampled.spd.calsample$grid$PrDens


simExponential = function(x,method)
{
  require(Bchron)
  start = 7000
  end = 3000
  diff = abs(start-end)
  CalBP = 1:(diff+1)
  errors= c(20,20)
  n=200
  r = -x
  PrDens = (1+r)^CalBP / sum((1+r)^CalBP)
  d = data.frame(CalBP=start:end,PrDens=PrDens)
  class(d)='CalGrid'
  d = uncalibrate(d,verbose=F)
  if (method=='uncalsample'){sampled.c14.dates=sample(d$CRA,size=n,replace=T,prob=d$PrDens)}
  if (method=='calsample'){sampled.c14.dates=sample(d$CRA,size=n,replace=T,prob=d$Raw)}
  #sampled.calibrated.dates = calibrate(sampled.c14.dates,errors=sample(errors,replace=TRUE,size=n),verbose=FALSE) 
  sampled.calibrated.dates = BchronCalibrate(ages=sampled.c14.dates,ageSds = sample(errors,replace=TRUE,size=n),calCurves=rep('intcal13',n))
  sampled.calibrated.dates = as.CalDates(sampled.calibrated.dates)
  sampled.calibrated.dates$metadata$StartBP=50000
  sampled.calibrated.dates$metadata$EndBP=0
  sampled.spd = spd(sampled.calibrated.dates,timeRange=c(start,end),spdnormalised = TRUE,verbose=FALSE)
  observed = sampled.spd$grid$PrDens
  return(observed)
}


tol=0.01
ncores = 5
nsim = 50000
cl <- makeCluster(ncores)
registerDoSNOW(cl)
params=rexp(nsim,30)
observed=observed.uncalsample
pb <- txtProgressBar(max = nsim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

reslist <- foreach (i=1:nsim,.packages='rcarbon',.options.snow = opts) %dopar% {
  set.seed(i)
  res.uncalsample=rev(simExponential(params[i],method='uncalsample'))
  set.seed(i)
  res.calsample=rev(simExponential(params[i],method='calsample'))
  euc_epsilon_uncalsample=sqrt(sum((res.uncalsample-observed)^2))
  euc_epsilon_calsample=sqrt(sum((res.calsample-observed)^2))
  ks_epsilon_uncalsample=max(cumsum(res.uncalsample)-cumsum(observed))
  ks_epsilon_calsample=max(cumsum(res.calsample)-cumsum(observed))
  
  tmp.uncal = res.uncalsample/sum(res.uncalsample) - observed/sum(observed)
  r2_uncal = 1-sum(tmp.uncal^2)/sum((observed-mean(observed))^2)
  tmp.cal = res.calsample/sum(res.calsample) - observed/sum(observed)
  r2_cal = 1-sum(tmp.cal^2)/sum((observed-mean(observed))^2)
  
  r_used = params[i]
  return(list(euc_epsilon_uncalsample,euc_epsilon_calsample,ks_epsilon_uncalsample,ks_epsilon_calsample,r2_uncal=r2_uncal,r2_cal=r2_cal,r_used))
}
close(pb)
stopCluster(cl) 


euc_epsilon_uncalsample = unlist(lapply(reslist,function(x){x[[1]]}))
euc_epsilon_calsample = unlist(lapply(reslist,function(x){x[[2]]}))
ks_epsilon_uncalsample = unlist(lapply(reslist,function(x){x[[3]]}))
ks_epsilon_calsample = unlist(lapply(reslist,function(x){x[[4]]}))
r2_calsample = unlist(lapply(reslist,function(x){x[[5]]}))
r2_uncalsample = unlist(lapply(reslist,function(x){x[[6]]}))
r_used = unlist(lapply(reslist,function(x){x[[7]]}))

save.image('testABC.RData')





