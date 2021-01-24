## Load Required R Libraries
library(Bchron)
library(readxl)
library(dplyr)
library(here)
## Setup ####

# Read Data
SSDP102.dates = read.csv(here('data','Kim_etal_2004','SSDP102_Kim_etlal_2004_c14.csv'),skip=1)
SSDP102.temp  = read.csv(here('data','Kim_etal_2004','SSDP102_Kim_etlal_2004_temp.csv'),skip=1)
pomaeho.dates = read.csv(here('data','Constantine_etal_2020','pomaeho_constantine_etal_2020_c14.csv'))
pomaeho.apt = read.csv(here('data','Constantine_etal_2020','pomaeho_constantine_etal_2020_aptp.csv'))
gy.dates = read.csv(here('data','Park_etal_2019','GY_Park_etal_2019_c14.csv'))
gy.apt = read.csv(here('data','Park_etal_2019','GY_Park_etal_2019_aptp.csv'))
## Constantine et al 2020, AP/T ####
## Delta R values calculated from the two dats at depth 1184 cm, using Reimar and Reimar 2016 calculation (DOI:10.1017/RDC.2016.117) and applet (http://calib.org/JS/JSdeltar20/)
DeltaR = -213
DeltaRError = 52
pomaeho.dates = subset(pomaeho.dates,CRA>1000)
pomaeho.dates$CRA[which(pomaeho.dates$Material2=='marine')]=pomaeho.dates$CRA[which(pomaeho.dates$Material2=='marine')]-DeltaR
pomaeho.dates$CRAError[which(pomaeho.dates$Material2=='marine')]=round(sqrt(pomaeho.dates$CRAError[which(pomaeho.dates$Material2=='marine')]^2+DeltaRError^2))
pomaeho.dates$calCurves = 'marine20'
pomaeho.dates$calCurves[which(pomaeho.dates$Material2=='terrestrial')] = 'intcal20'

outlierPresence = TRUE
pomaeho.dates$OutlierProb = 0
pomaeho.dates$Depth[nrow(pomaeho.dates)]=pomaeho.dates$Depth[nrow(pomaeho.dates)]+1 #added artificial offset as Bchron does not handle two dates on the same depth


set.seed(12345)
while(outlierPresence)
{
  index = which(pomaeho.dates$OutlierProb<0.1)
  pomaeho.model = Bchronology(ages=pomaeho.dates$CRA[index],ageSds=pomaeho.dates$CRAError[index],calCurves = pomaeho.dates$calCurves[index],ids=pomaeho.dates$LabCode[index],positions=pomaeho.dates$Depth[index],predictPositions=pomaeho.apt$Depth,jitterPositions = FALSE,iterations=50000, burn = 10000)
  outliers=summary(pomaeho.model, type='outliers')
  if (any(outliers$OutlierProb>0.1))
  {
    pomaeho.dates$OutlierProb[match(outliers$Date,pomaeho.dates$LabCode)]=outliers$OutlierProb
  } else {outlierPresence=FALSE}
}

## Kim et al 2004, Alkenone ####
# Calculated Weighted Mean by combining Delta R estimates from Lee & Kim XXXX and Kong et al 2005
deltaRs = c(-94,-71,-296,-253)
deltaRSigmas = c(22,24,35,45)
n = length(deltaRs)


deltaR.weightedMean= round(sum(deltaRs/deltaRSigmas^2)/sum(1/deltaRSigmas^2))
deltaRerr = round(sqrt(( (1/(n-1)) * sum(((deltaRs-deltaR.weightedMean)/deltaRSigmas)^2) )/((1/n)*sum(1/deltaRSigmas^2))))

outlierPresence = TRUE
SSDP102.dates$OutlierProb = 0

set.seed(12345)
while(outlierPresence)
{
  index = which(SSDP102.dates$OutlierProb<0.1)
  SSDP102.model = Bchronology(ages=round(SSDP102.dates$T2L_SSDP_102_c14_date - deltaR.weightedMean)[index],ageSds=round(sqrt(SSDP102.dates$T2L_SSDP_102_c14_1s_err^2+deltaRerr^2))[index],calCurves = rep('marine20',nrow(SSDP102.dates))[index],ids=SSDP102.dates$T2L_SSDP_102_labcode[index],positions=SSDP102.dates$T2L_SSDP_102_depth_top[index],predictPositions=SSDP102.temp$T2L_SSDP_102_depth,iterations=50000, burn = 10000)
  outliers=summary(SSDP102.model, type='outliers')
  if (any(outliers$OutlierProb>0.1))
  {
    SSDP102.dates$OutlierProb[match(outliers$Date,SSDP102.dates$T2L_SSDP_102_labcode)]=outliers$OutlierProb
  } else {outlierPresence=FALSE}
}

## Park et al 2019, AP/TP ####
outlierPresence = TRUE
gy.dates$OutlierProb = 0
while(outlierPresence)
{
  index = which(gy.dates$OutlierProb<0.1)
  gy.model = Bchronology(ages=gy.dates$CRA[index],ageSds=gy.dates$CRAError[index],calCurves = rep('intcal20',nrow(gy.dates))[index],ids=gy.dates$LabCode[index],positions=gy.dates$Depth[index],predictPositions=gy.apt$Depth,iterations=50000, burn = 10000)
  outliers=summary(gy.model, type='outliers')
  if (any(outliers$OutlierProb>0.1))
  {
    gy.dates$OutlierProb[match(outliers$Date,gy.dates$LabCode)]=outliers$OutlierProb
  } else {outlierPresence=FALSE}
}




# Extract Median Dates
med.pomaeho.model= apply(pomaeho.model$thetaPredict,2,median)
med.SSDP102.model= apply(SSDP102.model$thetaPredict,2,median)
med.gy.model= apply(gy.model$thetaPredict,2,median)



save(med.pomaeho.model,med.SSDP102.model,med.gy.model,pomaeho.apt,pomaeho.dates,pomaeho.model,SSDP102.dates,SSDP102.temp,SSDP102.model,gy.dates,gy.apt,gy.model,file=here('R_image_files/','agedepthmodels.RData'))


