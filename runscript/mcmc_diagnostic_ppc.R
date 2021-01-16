# Load libraries ####
library(nimbleCarbon)
library(coda)
library(here)

# Load MCMC results  ####
load(here('R_image_files/','mcmc_samples_coastal.RData'))
load(here('R_image_files/','mcmc_samples_inland.RData'))
load(here('R_image_files/','mcmc_samples_all.RData'))

### Load and Pre-Process Observed Data ####
load(here('R_image_files','koreanC14.RData'))

# subset coastal/inland
coastal.koreaC14 = subset(koreaC14,region=='coastal')
caldates.coastal = caldates[which(koreaC14$region=='coastal')]
bins.coastal = bins[which(koreaC14$region=='coastal')]

inland.koreaC14 = subset(koreaC14,region=='inland')
caldates.inland = caldates[which(koreaC14$region=='inland')]
bins.inland = bins[which(koreaC14$region=='inland')]

# Thin Dates per Bin
thinned.dates.coastal = caldates.coastal[thinDates(ages=coastal.koreaC14$c14age,  errors=coastal.koreaC14$c14error, bins=bins.coastal, size=1, thresh=1,seed=123,method='splitsample')]
thinned.dates.inland = caldates.inland[thinDates(ages=inland.koreaC14$c14age,  errors=inland.koreaC14$c14error, bins=bins.inland, size=1, thresh=1,seed=123,method='splitsample')]
thinned.dates = caldates[thinDates(ages=koreaC14$c14age,  errors=koreaC14$c14error, bins=bins, size=1, thresh=1,seed=123,method='splitsample')]

# Consider samples within window of analysis
thinned.dates = subset(thinned.dates,BP<=7000&BP>=3000,p=0.5)
thinned.dates.inland = subset(thinned.dates.inland,BP<=7000&BP>=3000,p=0.5)
thinned.dates.coastal = subset(thinned.dates.coastal,BP<=7000&BP>=3000,p=0.5)

# Extract CRA and CRAErrors
obs.CRA.all = thinned.dates$metadata$CRA
obs.Errors.all = thinned.dates$metadata$Error
obs.CRA.inland = thinned.dates.inland$metadata$CRA
obs.Errors.inland = thinned.dates.inland$metadata$Error
obs.CRA.coastal = thinned.dates.coastal$metadata$CRA
obs.Errors.coastal = thinned.dates.coastal$metadata$Error



# Compute Gelman-Rubin Convergence Statistics and ESS ####
rhat.coastal=gelman.diag(mcmc.samples.coastal)
rhat.inland=gelman.diag(mcmc.samples.inland)
rhat.all=gelman.diag(mcmc.samples.all)

ess.coastal=effectiveSize(mcmc.samples.coastal)
ess.inland=effectiveSize(mcmc.samples.inland)
ess.all=effectiveSize(mcmc.samples.all)

# Posterior Predictive Check ####
params.coastal = list(r1 = c(mcmc.samples.coastal$chain1[,'r1'],mcmc.samples.coastal$chain2[,'r1'],mcmc.samples.coastal$chain3[,'r1']),
                 r2 = c(mcmc.samples.coastal$chain1[,'r2'],mcmc.samples.coastal$chain2[,'r2'],mcmc.samples.coastal$chain3[,'r2']),
                 mu = round(c(mcmc.samples.coastal$chain1[,'chp'],mcmc.samples.coastal$chain2[,'chp'],mcmc.samples.coastal$chain3[,'chp'])))

params.inland = list(r1 = c(mcmc.samples.inland$chain1[,'r1'],mcmc.samples.inland$chain2[,'r1'],mcmc.samples.inland$chain3[,'r1']),
                      r2 = c(mcmc.samples.inland$chain1[,'r2'],mcmc.samples.inland$chain2[,'r2'],mcmc.samples.inland$chain3[,'r2']),
                      mu = round(c(mcmc.samples.inland$chain1[,'chp'],mcmc.samples.inland$chain2[,'chp'],mcmc.samples.inland$chain3[,'chp'])))

params.all = list(r1 = c(mcmc.samples.all$chain1[,'r1'],mcmc.samples.all$chain2[,'r1'],mcmc.samples.all$chain3[,'r1']),
                     r2 = c(mcmc.samples.all$chain1[,'r2'],mcmc.samples.all$chain2[,'r2'],mcmc.samples.all$chain3[,'r2']),
                     mu = round(c(mcmc.samples.all$chain1[,'chp'],mcmc.samples.all$chain2[,'chp'],mcmc.samples.all$chain3[,'chp'])))


params.coastal.calsample=postPredSPD(obs.CRA.coastal,obs.Errors.coastal,calCurve = 'intcal20',model = dDoubleExponentialGrowth,a = 7000,b=3000,params=params.coastal,nsim = 500,ncores = 5,verbose=FALSE,method='calsample')
params.coastal.uncalsample=postPredSPD(obs.CRA.coastal,obs.Errors.coastal,calCurve = 'intcal20',model = dDoubleExponentialGrowth,a = 7000,b=3000,params=params.coastal,nsim = 500,ncores = 5,verbose=FALSE,method='uncalsample')

params.inland.calsample=postPredSPD(obs.CRA.inland,obs.Errors.inland,calCurve = 'intcal20',model = dDoubleExponentialGrowth,a = 7000,b=3000,params=params.inland,nsim = 500,ncores = 5,verbose=FALSE,method='calsample')
params.inland.uncalsample=postPredSPD(obs.CRA.inland,obs.Errors.inland,calCurve = 'intcal20',model = dDoubleExponentialGrowth,a = 7000,b=3000,params=params.inland,nsim = 500,ncores = 5,verbose=FALSE,method='uncalsample')

params.all.calsample=postPredSPD(obs.CRA.all,obs.Errors.all,calCurve = 'intcal20',model = dDoubleExponentialGrowth,a = 7000,b=3000,params=params.all,nsim = 500,ncores = 5,verbose=FALSE,method='calsample')
params.all.uncalsample=postPredSPD(obs.CRA.all,obs.Errors.all,calCurve = 'intcal20',model = dDoubleExponentialGrowth,a = 7000,b=3000,params=params.all,nsim = 500,ncores = 5,verbose=FALSE,method='uncalsample')


save(rhat.coastal,ess.coastal,params.coastal,params.coastal.calsample,params.coastal.uncalsample,rhat.inland,ess.inland,params.inland,params.inland.calsample,params.inland.uncalsample,rhat.all,ess.all,params.all,params.all.calsample,params.all.uncalsample,file=here('R_image_files','mcmcdiagnostic_postpredcheck.RData'))




