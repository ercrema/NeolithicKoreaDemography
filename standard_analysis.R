# Load Files
load('./data/koreanC14.RData')
library(rcarbon)
# Core Parameters
bw = 100 #kernel bandwidth
runm = 100 #smoothing window
nsim = 5000 #number of MC iterations for modelTest()
timeRange = c(7000,3000) #time range of analysis
ncores = 5 #number of cores

# CKDE ####
ckdeNorm = ckde(sdates,timeRange=timeRange,bw=bw,normalised=TRUE)

# Inland/Coastal Comparison
costal.index = which(koreaC14$region=='coastal')
inland.index = which(koreaC14$region=='inland')

sdates.coastal = sampleDates(caldates[costal.index], bins = bins[costal.index], nsim=1000, boot = TRUE, verbose = TRUE)
sdates.inland = sampleDates(caldates[inland.index], bins = bins[inland.index], nsim=1000, boot = TRUE, verbose = TRUE)

ckdeNorm.coastal = ckde(sdates.coastal,timeRange=timeRange,bw=bw,normalised=TRUE)
ckdeNorm.inland = ckde(sdates.inland,timeRange=timeRange,bw=bw,normalised=TRUE)

# Permutation Test
CoastVsInland = permTest(caldates,marks=koreaC14$region,timeRange=timeRange,bins=bins,nsim=nsim,runm=runm)


# Against Theoretical Growth Models:

## Exponential ##
set.seed(123)
m1 <- modelTest(x=caldates,bins=bins,errors=koreaC14$c14error,model="exponential",runm=runm,nsim=nsim,ncores=ncores,timeRange=timeRange)

## Logistic ##
spdKorea <- spd(x=caldates,timeRange=timeRange,bins=bins,runm=runm,spdnormalised = TRUE)
obs.calBP <- spdKorea$grid$calBP
obs.PrDens <- spdKorea$grid$PrDens
logFit <- nls(obs.PrDens~SSlogis(obs.calBP, Asym, xmid, scale),control=nls.control(maxiter=200),start=list(Asym=0.2,xmid=5500,scale=-100))
set.seed(123)
logFitDens=data.frame(calBP=obs.calBP,PrDens=SSlogis(obs.calBP,coefficients(logFit)[1],coefficients(logFit)[2],coefficients(logFit)[3]))
LogCheck <- modelTest(caldates,bins=bins,errors=koreaC14$c14error,model="custom",runm=runm,nsim=nsim,ncores=ncores,timeRange=timeRange,predgrid=logFitDens,spdnormalised = TRUE)

save(m1,LogCheck,ckdeNorm.coastal,ckdeNorm.inland,CoastVsInland,file='./results_images/result_standard_analysis.RData')
