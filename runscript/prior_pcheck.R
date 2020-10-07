# Load ABC results and simulation model
load('../R_image_files/koreanC14.RData')
source('../src/sim_model.R')
source('../src/fastCalibrate.R')
library(rcarbon)
nsim = 100 
params = data.frame(bl=rexp(nsim,500),br=rexp(nsim,500),c=round(runif(nsim,4000,6000)))
prior.check.uncal.coastal= prior.check.uncal.inland =matrix(NA,nrow=length(7000:3000),ncol=nsim)
inland.index = which(koreaC14$region=='inland')
coastal.index = which(koreaC14$region=='coastal')


for (i in 1:nsim)
{
  set.seed(i)
  print(i)
  tmp.uncal=sim.model(x=koreaC14[coastal.index,],caldates=caldates[coastal.index],bins=bins[coastal.index],bl=params$bl[i],br=params$br[i],c=params$c[i],timeRange=c(7000,3000),simonly=TRUE,method='uncalsample')
  
  prior.check.uncal.coastal[,i] = tmp.uncal$grid$PrDens
  
  tmp.uncal=sim.model(x=koreaC14[inland.index,],caldates=caldates[inland.index],bins=bins[inland.index],bl=params$bl[i],br=params$br[i],c=params$c[i],timeRange=c(7000,3000),simonly=TRUE,method='uncalsample')
  
  prior.check.uncal.inland[,i] = tmp.uncal$grid$PrDens
  
  
}

save(prior.check.uncal.coastal,prior.check.uncal.inland,file='../R_image_files/prior_predictive_check.RData')
