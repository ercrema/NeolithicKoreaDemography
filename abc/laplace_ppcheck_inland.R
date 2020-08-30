# Load ABC results and simulation model
load('../results_images/resABC_laplace_inland.RData')
load('../data/koreanC14.RData')
source('./sim_model.R')
source('./fastCalibrate.R')
library(rcarbon)
# Retrieve parameters for the pred.check.size best fit to observed SPD
pred.check.size = 100
best_fit_uncal.param.general = abc.inland[order(abc.inland$euc.uncal)[1:pred.check.size],]
best_fit_cal.param.general = abc.inland[order(abc.inland$euc.cal)[1:pred.check.size],]
inland.index = which(koreaC14$region=='inland')
ppcheck.uncal.inland = ppcheck.cal.inland =matrix(NA,nrow=length(7000:3000),ncol=pred.check.size)
pb = txtProgressBar(max = pred.check.size, style = 3)

for (i in 1:pred.check.size)
{
  setTxtProgressBar(pb, i)
  tmp.uncal=sim.model(x=koreaC14[inland.index,],caldates=caldates[inland.index],bins=bins[inland.index],bl=best_fit_uncal.param.general$bl[i],br=best_fit_uncal.param.general$br[i],c=best_fit_uncal.param.general$c[i],timeRange=c(7000,3000),simonly=TRUE,method='uncalsample')
  
  ppcheck.uncal.inland[,i] = tmp.uncal$grid$PrDens
  
  tmp.cal=sim.model(x=koreaC14[inland.index,],caldates=caldates[inland.index],bins=bins[inland.index],bl=best_fit_cal.param.general$bl[i],br=best_fit_cal.param.general$br[i],c=best_fit_cal.param.general$c[i],timeRange=c(7000,3000),simonly=TRUE,method='calsample')
  
  ppcheck.cal.inland[,i] = tmp.cal$grid$PrDens
}
close(pb)

save(ppcheck.uncal.inland,ppcheck.cal.inland,file='../results_images/predcheck_results_inland.RData')
