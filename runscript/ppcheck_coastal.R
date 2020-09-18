# Load ABC results and simulation model
load('../R_image_files/resABC_laplace_coastal.RData')
load('../R_image_files/koreanC14.RData')
source('../src/sim_model.R')
source('../src/fastCalibrate.R')
library(rcarbon)
# Retrieve parameters for the pred.check.size best fit to observed SPD
pred.check.size = 100
best_fit_uncal.param.general = abc.coastal[order(abc.coastal$euc.uncal)[1:pred.check.size],]
best_fit_cal.param.general = abc.coastal[order(abc.coastal$euc.cal)[1:pred.check.size],]
coastal.index = which(koreaC14$region=='coastal')
ppcheck.uncal.coastal = ppcheck.cal.coastal =matrix(NA,nrow=length(7000:3000),ncol=pred.check.size)
pb = txtProgressBar(max = pred.check.size, style = 3)

for (i in 1:pred.check.size)
{
  setTxtProgressBar(pb, i)
  tmp.uncal=sim.model(x=koreaC14[coastal.index,],caldates=caldates[coastal.index],bins=bins[coastal.index],bl=best_fit_uncal.param.general$bl[i],br=best_fit_uncal.param.general$br[i],c=best_fit_uncal.param.general$c[i],timeRange=c(7000,3000),simonly=TRUE,method='uncalsample')
  
  ppcheck.uncal.coastal[,i] = tmp.uncal$grid$PrDens
  
  tmp.cal=sim.model(x=koreaC14[coastal.index,],caldates=caldates[coastal.index],bins=bins[coastal.index],bl=best_fit_cal.param.general$bl[i],br=best_fit_cal.param.general$br[i],c=best_fit_cal.param.general$c[i],timeRange=c(7000,3000),simonly=TRUE,method='calsample')
  
  ppcheck.cal.coastal[,i] = tmp.cal$grid$PrDens
}
close(pb)

save(ppcheck.uncal.coastal,ppcheck.cal.coastal,file='../R_image_files/predcheck_results_coastal.RData')
