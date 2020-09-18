# Load ABC results and simulation model
load('../R_image_files/resABC_laplace_general.RData')
load('../R_image_files/koreanC14.RData')
source('../src/sim_model.R')
source('../src/fastCalibrate.R')
library(rcarbon)
pred.check.size = 100 
best_fit_uncal.param.general = abc.general[order(abc.general$euc.uncal)[1:pred.check.size],]
best_fit_cal.param.general = abc.general[order(abc.general$euc.cal)[1:pred.check.size],]

ppcheck.uncal.general = ppcheck.cal.general =matrix(NA,nrow=length(7000:3000),ncol=pred.check.size)

for (i in 1:pred.check.size)
{
  set.seed(i)
  print(i)
  tmp.uncal=sim.model(x=koreaC14,caldates=caldates,bins=bins,bl=best_fit_uncal.param.general$bl[i],br=best_fit_uncal.param.general$br[i],c=best_fit_uncal.param.general$c[i],timeRange=c(7000,3000),simonly=TRUE,method='uncalsample')
  
  ppcheck.uncal.general[,i] = tmp.uncal$grid$PrDens
  
  set.seed(i)
  print(i)
  tmp.cal=sim.model(x=koreaC14,caldates=caldates,bins=bins,bl=best_fit_cal.param.general$bl[i],br=best_fit_cal.param.general$br[i],c=best_fit_cal.param.general$c[i],timeRange=c(7000,3000),simonly=TRUE,method='calsample')
  
  ppcheck.cal.general[,i] = tmp.cal$grid$PrDens
}

save(ppcheck.uncal.general,ppcheck.cal.general,file='../results_images/predcheck_results_general.RData')
