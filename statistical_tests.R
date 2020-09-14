# Load Files
load('./data/koreanC14.RData')
library(rcarbon)
# Core Parameters
runm = 100 #smoothing window
nsim = 1000 #number of MC iterations for modelTest()
timeRange = c(7000,3000) #time range of analysis
ncores = 6 #number of cores

# Inland/Coastal Index ####
coastal.index = which(koreaC14$region=='coastal')
inland.index = which(koreaC14$region=='inland')

# Permutation Test (Coastal vs Inland) ####
coastal.inland.permtest = permTest(caldates,marks=koreaC14$region,timeRange=timeRange,bins=bins,nsim=nsim,runm=runm)

# Permutation Test ####
koreaC14$milletAsso[which(koreaC14$milletAsso)]='yes'
koreaC14$milletAsso[is.na(koreaC14$milletAsso)]='no'
allDates=calibrate(koreaC14$c14age,koreaC14$c14error)
millet.permtest=permTest(allDates,marks=koreaC14$milletAsso,timeRange=timeRange,runm=runm,nsim=nsim)

# Theoretical Growth Model Tests ####

## Exponential ##
set.seed(123)
exp.general.test <- modelTest(x=caldates,bins=bins,errors=koreaC14$c14error,model="exponential",runm=runm,nsim=nsim,ncores=ncores,timeRange=timeRange)

exp.coastal.test <- modelTest(x=caldates[coastal.index],bins=bins[coastal.index],errors=koreaC14$c14error[coastal.index],model="exponential",runm=runm,nsim=nsim,ncores=ncores,timeRange=timeRange)

exp.inland.test <- modelTest(x=caldates[inland.index],bins=bins[inland.index],errors=koreaC14$c14error[inland.index],model="exponential",runm=runm,nsim=nsim,ncores=ncores,timeRange=timeRange)


## Logistic ##
spd.general <- spd(x=caldates,timeRange=timeRange,bins=bins,runm=runm,spdnormalised = TRUE)
spd.coastal <- spd(x=caldates[coastal.index],timeRange=timeRange,bins=bins[coastal.index],runm=runm,spdnormalised = TRUE)
spd.inland <- spd(x=caldates[inland.index],timeRange=timeRange,bins=bins[inland.index],runm=runm,spdnormalised = TRUE)

obs.calBP <- timeRange[1]:timeRange[2]
prdens.general <- spd.general$grid$PrDens
prdens.coastal <- spd.coastal$grid$PrDens
prdens.inland <- spd.inland$grid$PrDens

logFit.general <- nls(prdens.general~SSlogis(obs.calBP, Asym, xmid, scale),control=nls.control(maxiter=200),start=list(Asym=0.2,xmid=5500,scale=-100))
logFit.coastal <- nls(prdens.coastal~SSlogis(obs.calBP, Asym, xmid, scale),control=nls.control(maxiter=200),start=list(Asym=0.2,xmid=5500,scale=-100))
logFit.inland<- nls(prdens.inland~SSlogis(obs.calBP, Asym, xmid, scale),control=nls.control(maxiter=200),start=list(Asym=0.2,xmid=5500,scale=-100))

logistic.general=data.frame(calBP=obs.calBP,PrDens=SSlogis(obs.calBP,coefficients(logFit.general)[1],coefficients(logFit.general)[2],coefficients(logFit.general)[3]))
logistic.coastal=data.frame(calBP=obs.calBP,PrDens=SSlogis(obs.calBP,coefficients(logFit.coastal)[1],coefficients(logFit.coastal)[2],coefficients(logFit.coastal)[3]))
logistic.inland=data.frame(calBP=obs.calBP,PrDens=SSlogis(obs.calBP,coefficients(logFit.inland)[1],coefficients(logFit.inland)[2],coefficients(logFit.inland)[3]))

logistic.general.test <- modelTest(caldates,bins=bins,errors=koreaC14$c14error,model="custom",runm=runm,nsim=nsim,ncores=ncores,timeRange=timeRange,predgrid=logistic.general,spdnormalised = TRUE)
logistic.coastal.test <- modelTest(caldates[coastal.index],bins=bins[coastal.index],errors=koreaC14$c14error[coastal.index],model="custom",runm=runm,nsim=nsim,ncores=ncores,timeRange=timeRange,predgrid=logistic.general,spdnormalised = TRUE)
logistic.inland.test <- modelTest(caldates[inland.index],bins=bins[inland.index],errors=koreaC14$c14error[inland.index],model="custom",runm=runm,nsim=nsim,ncores=ncores,timeRange=timeRange,predgrid=logistic.general,spdnormalised = TRUE)

# Save Output Image #####

# assign NULL to fitobject to reduce image size
exp.general.test$fitobject=NULL
exp.coastal.test$fitobject=NULL
exp.inland.test$fitobject=NULL

save(coastal.inland.permtest,millet.permtest,exp.general.test,exp.coastal.test,exp.inland.test,logistic.general.test,logistic.coastal.test,logistic.inland.test,file='./results_images/test_results.RData')


