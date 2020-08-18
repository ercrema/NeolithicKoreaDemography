library(rcarbon)
library(foreach)
library(doParallel)
library(doSNOW)

### Load Observed Data ###
load('../data/koreanC14.RData')
inland.koreaC14 = subset(koreaC14,region=='inland')
caldates = caldates[which(koreaC14$region=='inland')]
bins = bins[which(koreaC14$region=='inland')]
### Load Core Simulation Model and fastCalibrate() function ###
source('./sim_model.R')
source('./fastCalibrate.R')

### Generate Parameters ####
nsim=100000
bl=rexp(nsim,500)
br=rexp(nsim,500)
c=runif(nsim,4000,6000)
### ABC Settings ###
tol=0.01
ncores = 4
cl <- makeCluster(ncores)
registerDoSNOW(cl)
pb <- txtProgressBar(max = nsim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

reslist <- foreach (i=1:nsim,.packages=c('rcarbon'),.options.snow = opts) %dopar% {
  res=sim.model(x=inland.koreaC14,caldates=caldates,bins=bins,bl=bl[i],br=br[i],c=c[i],timeRange=c(7000,3000))
  return(res)
}


save(reslist,file='../results_images/resABC_laplace_inland.RData')

bl = unlist(lapply(reslist,function(x){x[[1]]}))
br = unlist(lapply(reslist,function(x){x[[2]]}))
c = unlist(lapply(reslist,function(x){x[[3]]}))
euc.uncal = unlist(lapply(reslist,function(x){x[[4]]}))
euc.cal = unlist(lapply(reslist,function(x){x[[5]]}))


abc.inland=data.frame(bl=bl,br=br,c=c,euc.uncal=euc.uncal,euc.cal=euc.cal)

save(abc.inland,file='../results_images/resABC_laplace_inland.RData')

