library(rcarbon)
library(foreach)
library(doParallel)
library(doSNOW)

### Load Observed Data ###
load('../R_image_files/koreanC14.RData')
coastal.koreaC14 = subset(koreaC14,region=='coastal')
caldates = caldates[which(koreaC14$region=='coastal')]
bins = bins[which(koreaC14$region=='coastal')]
### Load Core Simulation Model and fastCalibrate() function ###
source('../src/sim_model.R')
source('../src/fastCalibrate.R')

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
  res=sim.model(x=coastal.koreaC14,caldates=caldates,bins=bins,bl=bl[i],br=br[i],c=c[i],timeRange=c(7000,3000))
  return(res)
}


save(reslist,file='../results_images/resABC_laplace_coastal.RData')

bl = unlist(lapply(reslist,function(x){x[[1]]}))
br = unlist(lapply(reslist,function(x){x[[2]]}))
c = unlist(lapply(reslist,function(x){x[[3]]}))
euc.uncal = unlist(lapply(reslist,function(x){x[[4]]}))
euc.cal = unlist(lapply(reslist,function(x){x[[5]]}))


abc.coastal=data.frame(bl=bl,br=br,c=c,euc.uncal=euc.uncal,euc.cal=euc.cal)

save(abc.coastal,file='../R_image_files/resABC_laplace_coastal.RData')

