# Load Libraries
library(here)
library(rcarbon)
load(here('R_image_files','koreanC14.RData'))

# Core Parameters
runm = 100 #smoothing window
nsim = 1000 #number of MC iterations for modelTest()
timeRange = c(7000,3000) #time range of analysis
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

save(coastal.inland.permtest,millet.permtest,file=here('R_image_files','spd_test_results.RData'))


