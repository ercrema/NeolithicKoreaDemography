### Load Required R packages ####
library(nimbleCarbon)
library(here)
library(truncnorm)

### Load and Pre-Process Observed Data ####
load(here('R_image_files','koreanC14.RData'))
inland.koreaC14 = subset(koreaC14,region=='inland')
caldates = caldates[which(koreaC14$region=='inland')]
bins = bins[which(koreaC14$region=='inland')]
# Thin Dates per Bin
thinned.dates = caldates[thinDates(ages=inland.koreaC14$c14age,  errors=inland.koreaC14$c14error, bins=bins, size=1, thresh=1,seed=123,method='splitsample')]
# Consider samples within window of analysis
thinned.dates = subset(thinned.dates,BP<=7000&BP>=3000,p=0.5)

# Extract CRA and CRAErrors
obs.CRA = thinned.dates$metadata$CRA
obs.Errors = thinned.dates$metadata$Error

### Bayesian Analysis ####

# Setup Data and Constants
data("intcal20")
constants <- list(N=length(obs.CRA),calBP=intcal20$CalBP,C14BP=intcal20$C14Age,C14err=intcal20$C14Age.sigma,start=7000,end=3000)
data <- list(X=obs.CRA,sigma=obs.Errors)

# Define Model
model <- nimbleCode({
  for (i in 1:N){
    theta[i] ~ dDoubleExponentialGrowth(a=7000,b=3000,r1=r1,r2=r2,mu=changept);
    mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
    sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
    sd[i] <- (sigma[i]^2+sigmaCurve[i]^2)^(1/2);
    X[i] ~ dnorm(mean=mu[i],sd=sd[i]);
  }
  r1 ~ dnorm(0,sd=0.0004); 
  r2 ~ dnorm(0,sd=0.0004);
  chp ~ T(dnorm(5000,sd=1000),3000,7000);
  changept <- round(chp);
})  

# Create Init Function
a = 7000
b = 3000
m.dates = medCal(thinned.dates)
if(any(m.dates>a|m.dates<b)){m.dates[m.dates>a]=a;m.dates[m.dates<b]=b}
initsFunction = function() list(r1=rnorm(1,sd=0.0004),r2=rnorm(1,sd=0.0004),chp=round(rtruncnorm(1,mean=5000,sd=1000,a=3000,b=7000)),theta=as.numeric(m.dates))

# Run MCMC
mcmc.samples.inland<- nimbleMCMC(code = model,constants = constants,data = data,niter = 100000, nchains = 3, thin=6, nburnin = 10000, summary = FALSE, monitors=c('r1','r2','chp'),WAIC=FALSE,samplesAsCodaMCMC=TRUE,inits=initsFunction,setSeed=c(1,2,3))

gelman.diag(mcmc.samples.inland)$psrf[1:3,]

# Store output
save(mcmc.samples.inland,file=here('R_image_files','mcmc_samples_inland.RData'))


