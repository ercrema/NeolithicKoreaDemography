sim.model = function(x,caldates,bins,a,bl,br,c,timeRange=c(7000,3000),simonly=FALSE,method='both')
{
  if (!simonly)
  {
    # Generate Random Target
    index=thinDates(x$c14age,x$c14error,bins,size=1,thresh=0,method='random')
    target.spd =spd(caldates[index],timeRange = timeRange,datenormalised=TRUE,spdnormalised = TRUE,verbose=FALSE)
    tprob = target.spd$grid$PrDens
  }
  trange = timeRange[1]:timeRange[2]
  
  # Create Model
  CalBP = trange
  PrDens=ifelse(CalBP < c, a * exp((CalBP-c) / br), a * exp(-(CalBP-c) / bl))
  d = data.frame(CalBP=CalBP,PrDens=PrDens/sum(PrDens))
  # Collect Samples
  class(d)='CalGrid'
  d = uncalibrate(d,verbose=F)
  n = length(unique(bins))
  if (simonly) {errors=x$c14error}
  if (!simonly) {errors=x$c14error[index]}
    
    
  if (method=='uncalsample'|method=='both')
  {
    sampled.c14.dates.uncal=sample(d$CRA,size=n,replace=T,prob=d$PrDens)
    sampled.calibrated.dates.uncal = BchronCalibrate(ages=sampled.c14.dates.uncal,ageSds = sample(errors,replace=TRUE,size=n),calCurves=rep('intcal13',n))
    sampled.calibrated.dates.uncal = as.CalDates(sampled.calibrated.dates.uncal)
    sampled.calibrated.dates.uncal$metadata$StartBP=50000
    sampled.calibrated.dates.uncal$metadata$EndBP=0
    sampled.spd.uncal = spd(sampled.calibrated.dates.uncal,timeRange=timeRange,spdnormalised = TRUE,verbose=FALSE)
    candidate.uncal = sampled.spd.uncal$grid$PrDens
  }
  if (method=='calsample'|method=='both')
  {
    sampled.c14.dates.cal=sample(d$CRA,size=n,replace=T,prob=d$Raw)
    sampled.calibrated.dates.cal = BchronCalibrate(ages=sampled.c14.dates.cal,ageSds = sample(errors,replace=TRUE,size=n),calCurves=rep('intcal13',n))
    sampled.calibrated.dates.cal = as.CalDates(sampled.calibrated.dates.cal)
    sampled.calibrated.dates.cal$metadata$StartBP=50000
    sampled.calibrated.dates.cal$metadata$EndBP=0
    sampled.spd.cal = spd(sampled.calibrated.dates.cal,timeRange=timeRange,spdnormalised = TRUE,verbose=FALSE)
    candidate.cal = sampled.spd.cal$grid$PrDens
  }
  
  if (!simonly)
  {
    # Compute dissimilarity between target and simulation
    if (method=='uncalsample'|method=='both')
    {
    euc_epsilon_uncalsample=sqrt(sum((candidate.uncal-tprob)^2))
    }
    if (method=='calsample'|method=='both')
    {
      euc_epsilon_calsample=sqrt(sum((candidate.cal-tprob)^2))
    }
  }
    
  if (!simonly)
  {
    if (method=='both')
    {
    return(list(a=a,bl=bl,br=br,c=c,euc_epsilon_uncalsample=euc_epsilon_uncalsample,euc_epsilon_calsample=euc_epsilon_calsample))
    }
    if (method=='uncalsample')
    {
      return(list(a=a,bl=bl,br=br,c=c,euc_epsilon_uncalsample=euc_epsilon_uncalsample))
    }
    if (method=='calsample')
    {
      return(list(a=a,bl=bl,br=br,c=c,euc_epsilon_calsample=euc_epsilon_calsample))
    }
  }
  
  if (simonly)
  {
    if (method=='both')
    {
      return(list(sampled.spd.cal=sampled.spd.cal,sampled.spd.uncal=sampled.spd.uncal))
    }
    if (method=='uncalsample')
    {
      return(sampled.spd.uncal)
    }
    if (method=='calsample')
    {
      return(sampled.spd.cal)
    }
  }
}
