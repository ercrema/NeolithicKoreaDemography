# Load Libraries, Data, and Results ####
library(rcarbon)
library(coda)
library(here)
load(here('R_image_files','spd_test_results.RData'))
load(here('R_image_files','agedepthmodels.RData'))
load(here('R_image_files','mcmcdiagnostic_postpredcheck.RData'))




# Age Depth Model Posterior ####
point_d1=SSDP102.model$thetaPredict[,81]
point_d2=SSDP102.model$thetaPredict[,80]
point_f1=pomaeho.model$thetaPredict[,73]
point_f2=pomaeho.model$thetaPredict[,72]
point_g1=gy.model$thetaPredict[,79]
point_g2=gy.model$thetaPredict[,78]

age.depth.res90HPD = data.frame(Event=c('d1','d2','f1','f2','g1','g2'),Start=NA,End=NA)
age.depth.res90HPD$Start[1] = HPDinterval(mcmc(point_d1),prob = 0.90)[2]
age.depth.res90HPD$End[1] = HPDinterval(mcmc(point_d1),prob = 0.90)[1]
age.depth.res90HPD$Start[2] = HPDinterval(mcmc(point_d2),prob = 0.90)[2]
age.depth.res90HPD$End[2] = HPDinterval(mcmc(point_d2),prob = 0.90)[1]
age.depth.res90HPD$Start[3] = HPDinterval(mcmc(point_f1),prob = 0.90)[2]
age.depth.res90HPD$End[3] = HPDinterval(mcmc(point_f1),prob = 0.90)[1]
age.depth.res90HPD$Start[4] = HPDinterval(mcmc(point_f2),prob = 0.90)[2]
age.depth.res90HPD$End[4] = HPDinterval(mcmc(point_f2),prob = 0.90)[1]
age.depth.res90HPD$Start[5] = HPDinterval(mcmc(point_g1),prob = 0.90)[2]
age.depth.res90HPD$End[5] = HPDinterval(mcmc(point_g1),prob = 0.90)[1]
age.depth.res90HPD$Start[6] = HPDinterval(mcmc(point_g2),prob = 0.90)[2]
age.depth.res90HPD$End[6] = HPDinterval(mcmc(point_g2),prob = 0.90)[1]

write.csv(age.depth.res90HPD,file=here('figures_and_results','age_depth_hpd90.csv'),row.names=FALSE)


# Bayesian Model Diagnostics and Posterior HPDs ####
options(scipen=9999)
mcmc_summary = data.frame(parameter=c(rep('r1',3),rep('r2',3),rep('c',3)),data=rep(c('Inland','Coastal','Combined'),3),HPDstart=NA,HPDEnd=NA,ESS=NA,Rhat=NA)
mcmc_summary$HPDstart[1] = HPDinterval(mcmc(params.inland$r1))[1]
mcmc_summary$HPDEnd[1] = HPDinterval(mcmc(params.inland$r1))[2]
mcmc_summary$HPDstart[2] = HPDinterval(mcmc(params.coastal$r1))[1]
mcmc_summary$HPDEnd[2] = HPDinterval(mcmc(params.coastal$r1))[2]
mcmc_summary$HPDstart[3] = HPDinterval(mcmc(params.all$r1))[1]
mcmc_summary$HPDEnd[3] = HPDinterval(mcmc(params.all$r1))[2]
mcmc_summary$HPDstart[4] = HPDinterval(mcmc(params.inland$r2))[1]
mcmc_summary$HPDEnd[4] = HPDinterval(mcmc(params.inland$r2))[2]
mcmc_summary$HPDstart[5] = HPDinterval(mcmc(params.coastal$r2))[1]
mcmc_summary$HPDEnd[5] = HPDinterval(mcmc(params.coastal$r2))[2]
mcmc_summary$HPDstart[6] = HPDinterval(mcmc(params.all$r2))[1]
mcmc_summary$HPDEnd[6] = HPDinterval(mcmc(params.all$r2))[2]
mcmc_summary$HPDstart[7] = round(HPDinterval(mcmc(params.inland$mu))[2])
mcmc_summary$HPDEnd[7] = round(HPDinterval(mcmc(params.inland$mu))[1])
mcmc_summary$HPDstart[8] = round(HPDinterval(mcmc(params.coastal$mu))[2])
mcmc_summary$HPDEnd[8] = round(HPDinterval(mcmc(params.coastal$mu))[1])
mcmc_summary$HPDstart[9] = round(HPDinterval(mcmc(params.all$mu))[2])
mcmc_summary$HPDEnd[9] = round(HPDinterval(mcmc(params.all$mu))[1])

mcmc_summary$ESS[1]=ess.inland[2]
mcmc_summary$ESS[2]=ess.coastal[2]
mcmc_summary$ESS[3]=ess.all[2]
mcmc_summary$ESS[4]=ess.inland[3]
mcmc_summary$ESS[5]=ess.coastal[3]
mcmc_summary$ESS[6]=ess.all[3]
mcmc_summary$ESS[7]=ess.inland[1]
mcmc_summary$ESS[8]=ess.coastal[1]
mcmc_summary$ESS[9]=ess.all[1]

mcmc_summary$Rhat[1]=rhat.inland$psrf[2,1]
mcmc_summary$Rhat[2]=rhat.coastal$psrf[2,1]
mcmc_summary$Rhat[3]=rhat.all$psrf[2,1]
mcmc_summary$Rhat[4]=rhat.inland$psrf[3,1]
mcmc_summary$Rhat[5]=rhat.coastal$psrf[3,1]
mcmc_summary$Rhat[6]=rhat.all$psrf[3,1]
mcmc_summary$Rhat[7]=rhat.inland$psrf[1,1]
mcmc_summary$Rhat[8]=rhat.coastal$psrf[1,1]
mcmc_summary$Rhat[9]=rhat.all$psrf[1,1]


write.csv(mcmc_summary,file=here('figures_and_results','mcmc_summary.csv'),row.names=FALSE)




