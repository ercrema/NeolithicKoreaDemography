### Generate Sample Data with Known Growth Rate ###
library(rcarbon)
library(Bchron)
start = 7000
end = 3000
r = -0.0002
n = 200 #sample size
CalBP = 7000:3000
PrDens = (1+r)^CalBP / sum((1+r)^CalBP)
#plot(CalBP,PrDens,xlim=c(7000,3000))

## uncalibrate option 1 (individual dates):
set.seed(1)
sampled.calendar.dates = sample(CalBP,size=n,replace=TRUE,prob=PrDens)
sampled.c14.dates.uncalsample = uncalibrate(sampled.calendar.dates)$rCRA
d = data.frame(CalBP=CalBP,PrDens=PrDens)
class(d)='CalGrid'
d = uncalibrate(d,verbose=F)
sampled.c14.dates.calsample=sample(d$CRA,size=n,replace=T,prob=d$PrDens)

# par(mfrow=c(3,1))
# plot(as.numeric(names(table(sampled.c14.dates1))),table(sampled.c14.dates1),type='l')
# plot(as.numeric(names(table(sampled.c14.dates2))),table(sampled.c14.dates2),type='l')
# plot(as.numeric(names(table(sampled.c14.dates3))),table(sampled.c14.dates3),type='l')

sampled.calibrated.dates.uncalsample = calibrate(sampled.c14.dates.uncalsample,errors=rep(20,n)) 
sampled.calibrated.dates.calsample = calibrate(sampled.c14.dates.calsample,errors=rep(20,n)) 

sampled.spd.uncalsample = spd(sampled.calibrated.dates.uncalsample,timeRange=c(7000,3000),spdnormalised = TRUE)
sampled.spd.calsample = spd(sampled.calibrated.dates.calsample,timeRange=c(7000,3000),spdnormalised = TRUE)

observed = sampled.spd.uncalsample$grid$PrDens


simExponential = function(x)
{
  require(Bchron)
  start = 7000
  end = 3000
  diff = abs(start-end)
  CalBP = 1:(diff+1)
  errors= c(20,20)
  method='calsample'
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


set.seed(1)
nsim=10000
params=rexp(nsim,30)
euc_epsilon=ks_epsilon=numeric(length=nsim)
for (i in 1:nsim)
{
  print(i)
  set.seed(i)
  res=simExponential(params[i])
  euc_epsilon[i]=sqrt(sum((res-observed)^2))
  ks_epsilon[i]=max(cumsum(res)-cumsum(observed))
}

best01euc = order(euc_epsilon)[1:100]
best01kc = order(ks_epsilon)[1:100]

save.image('testABC.RData')



