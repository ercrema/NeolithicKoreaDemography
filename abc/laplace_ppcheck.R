# Load ABC results and simulation model
load('../results_images/resABC_laplace.RData')
load('../data/koreanC14.RData')
source('./sim_model.R')
library(rcarbon)
library(Bchron)
pred.check.size = 100 
best_fit_uncal.param = res[order(res$euc.uncal)[1:pred.check.size],]
best_fit_cal.param = res[order(res$euc.cal)[1:pred.check.size],]

ppcheck.uncal = ppcheck.cal =matrix(NA,nrow=length(7000:3000),ncol=pred.check.size)

for (i in 1:pred.check.size)
{
  set.seed(i)
  print(i)
  tmp.uncal=sim.model(x=koreaC14,caldates=caldates,bins=bins,bl=best_fit_uncal.param$bl[i],br=best_fit_uncal.param$br[i],c=best_fit_uncal.param$c[i],timeRange=c(7000,3000),simonly=TRUE,method='uncalsample')
  
  ppcheck.uncal[,i] = tmp.uncal$grid$PrDens
  
  set.seed(i)
  print(i)
  tmp.cal=sim.model(x=koreaC14,caldates=caldates,bins=bins,bl=best_fit_cal.param$bl[i],br=best_fit_cal.param$br[i],c=best_fit_cal.param$c[i],timeRange=c(7000,3000),simonly=TRUE,method='calsample')
  
  ppcheck.cal[,i] = tmp.cal$grid$PrDens
}

save(ppcheck.uncal,ppcheck.cal,file='../results_images/predcheck_results.RData')
