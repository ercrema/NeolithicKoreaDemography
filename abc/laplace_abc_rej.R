library(rcarbon)
library(Bchron)
library(foreach)
library(doParallel)
library(doSNOW)

### Load Observed Data ###
load('../data/koreanC14.RData')

### Load Core Simulation Model ###
source('./sim_model.R')

### Generate Parameters ####
nsim=100000
a=runif(nsim,0.0001,0.001)
bl=1/rexp(nsim,500)
br=1/rexp(nsim,500)
c=runif(nsim,4000,6000)

### ABC Settings ###
tol=0.01
ncores = 5
cl <- makeCluster(ncores)
registerDoSNOW(cl)
pb <- txtProgressBar(max = nsim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

reslist <- foreach (i=1:nsim,.packages=c('rcarbon','Bchron'),.options.snow = opts) %dopar% {
  res=sim.model(x=koreaC14,caldates=caldates,bins=bins,a=a[i],bl=bl[i],br=br[i],c=c[i],timeRange=c(7000,3000))
  return(res)
}


a = unlist(lapply(reslist,function(x){x[[1]]}))
bl = 1/unlist(lapply(reslist,function(x){x[[2]]}))
br = 1/unlist(lapply(reslist,function(x){x[[3]]}))
c = unlist(lapply(reslist,function(x){x[[4]]}))
euc.uncal = unlist(lapply(reslist,function(x){x[[5]]}))
euc.cal = unlist(lapply(reslist,function(x){x[[6]]}))
ks.uncal = unlist(lapply(reslist,function(x){x[[7]]}))
ks.cal = unlist(lapply(reslist,function(x){x[[8]]}))

res=data.frame(a=a,bl=bl,br=br,c=c,euc.uncal=euc.uncal,euc.cal=euc.cal,ks.uncal=ks.uncal,ks.cal=ks.cal)

save(res,file='../results_images/resABC_laplace.RData')

