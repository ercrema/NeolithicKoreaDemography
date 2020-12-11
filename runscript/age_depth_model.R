library(Bchron)
library(readxl)
library(dplyr)
## Setup ####

# Read Data
kim2004.dates = read.csv('../data/kim2004/SSDP_102.Kim.2004-chron.csv',skip=1)
kim2004.temp  = read.csv('../data/kim2004/SSDP_102.Kim.2004.csv',skip=1)
constantine2020.dates = read.csv('../data/Constantine2020/pomaeho_c14.csv')
constantine2020.apt = read_xlsx('../data/Constantine2020/pomaeho_pollen.xlsx',sheet = 'Sheet2')


## Constantine et al 2020, AP/T
DeltaR = -213
DeltaRError = 52
constantine2020.dates = subset(constantine2020.dates,CRA>1000)
constantine2020.dates$CRA[which(constantine2020.dates$Material2=='marine')]=constantine2020.dates$CRA[which(constantine2020.dates$Material2=='marine')]-DeltaR
constantine2020.dates$CRAError[which(constantine2020.dates$Material2=='marine')]=round(sqrt(constantine2020.dates$CRAError[which(constantine2020.dates$Material2=='marine')]^2+DeltaRError^2))
constantine2020.dates$calCurves = 'marine20'
constantine2020.dates$calCurves[which(constantine2020.dates$Material2=='terrestrial')] = 'intcal20'

outlierPresence = TRUE
constantine2020.dates$OutlierProb = 0
constantine2020.dates$Depth[nrow(constantine2020.dates)]=constantine2020.dates$Depth[nrow(constantine2020.dates)]+1


set.seed(12345)
while(outlierPresence)
{
  index = which(constantine2020.dates$OutlierProb<0.1)
  constantine2020.model = Bchronology(ages=constantine2020.dates$CRA[index],ageSds=constantine2020.dates$CRAError[index],calCurves = constantine2020.dates$calCurves[index],ids=constantine2020.dates$LabCode[index],positions=constantine2020.dates$Depth[index],predictPositions=constantine2020.apt$Depth,jitterPositions = FALSE,iterations=50000, burn = 10000)
  outliers=summary(constantine2020.model, type='outliers')
  if (any(outliers$OutlierProb>0.1))
  {
    constantine2020.dates$OutlierProb[match(outliers$Date,constantine2020.dates$LabCode)]=outliers$OutlierProb
  } else {outlierPresence=FALSE}
}

## Kim et al 2004, Alkenone


# Delta R (From http://calib.org/marine/index.html?npoints=1&clat=34.72845262222496&clon=128.2532094998296)
deltaRs = c(-94,-71)
deltaRSigmas = c(22,24)
mu = mean(deltaRs)
n = length(deltaRs)
deltaR.weightedMean = round(sum(deltaRs/deltaRSigmas)/sum(1/deltaRSigmas))
deltaRerr = round(sqrt(((1/(n-1)) * sum(((deltaRs-mean(deltaRs))/deltaRSigmas)^2))/((1/(n-1))*sum((1/deltaRSigmas)^2))))

outlierPresence = TRUE
kim2004.dates$OutlierProb = 0

set.seed(12345)
while(outlierPresence)
{
  index = which(kim2004.dates$OutlierProb<0.1)
  kim2004.model = Bchronology(ages=round(kim2004.dates$T2L_SSDP_102_c14_date - deltaR.weightedMean)[index],ageSds=round(sqrt(kim2004.dates$T2L_SSDP_102_c14_1s_err^2+deltaRerr^2))[index],calCurves = rep('marine20',nrow(kim2004.dates))[index],ids=kim2004.dates$T2L_SSDP_102_labcode[index],positions=kim2004.dates$T2L_SSDP_102_depth_top[index],predictPositions=kim2004.temp$T2L_SSDP_102_depth,iterations=50000, burn = 10000)
  outliers=summary(kim2004.model, type='outliers')
  if (any(outliers$OutlierProb>0.1))
  {
    kim2004.dates$OutlierProb[match(outliers$Date,kim2004.dates$T2L_SSDP_102_labcode)]=outliers$OutlierProb
  } else {outlierPresence=FALSE}
}


# Extract Median Dates
med.kim2004.model= apply(kim2004.model$thetaPredict,2,median)
med.constantine2020.model= apply(constantine2020.model$thetaPredict,2,median)

save(med.kim2004.model,med.constantine2020.model,kim2004.model,constantine2020.model,kim2004.dates,kim2004.temp,constantine2020.dates,constantine2020.apt,file='../R_image_files/agedepthmodels.RData')


