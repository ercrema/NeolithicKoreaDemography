library(Bchron)
# Read Kim 2004
kim2004.dates = read.csv('../data/kim2004/SSDP_102.Kim.2004-chron.csv',skip=1)
kim2004.temp  = read.csv('../data/kim2004/SSDP_102.Kim.2004.csv',skip=1)

# Delta R (From http://calib.org/marine/index.html?npoints=1&clat=34.72845262222496&clon=128.2532094998296)
deltaR = -296
deltaRerr = 35


# Read marine20 calibration curve
marine20 = read.csv('http://intcal.org/curves/marine20.14c', encoding="UTF-8",skip=11,header=F)
createCalCurve(name='marine20',calAges=marine20[,1],uncalAges=marine20[,2],oneSigma=marine20[,3])
file.copy(from = 'marine20.rda',to = system.file('data',package='Bchron'))

set.seed(12345)
kim2004.model.marine20 = Bchronology(ages=round(kim2004.dates$T2L_SSDP_102_c14_date - deltaR),ageSds=round(sqrt(kim2004.dates$T2L_SSDP_102_c14_1s_err^2+deltaRerr^2)),calCurves = rep('marine20',nrow(kim2004.dates)),ids=kim2004.dates$T2L_SSDP_102_labcode,positions=kim2004.dates$T2L_SSDP_102_depth_top,predictPositions=kim2004.temp$T2L_SSDP_102_depth)

# Extract Median Dates
med.bchron.marine20 = apply(kim2004.model.marine20$thetaPredict,2,median)

save(med.bchron.marine20,kim2004.model.marine20,kim2004.dates,kim2004.temp,file='../R_image_files/kim2004_agedepthmodel.RData')
