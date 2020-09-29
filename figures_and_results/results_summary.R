# Load Libraries, Data, and Results ####
library(rcarbon)
library(Bchron)
library(coda)

load('../R_image_files/spd_test_results.RData')
load('../R_image_files/kim2004_agedepthmodel.RData')
load('../R_image_files/resABC_laplace_general.RData')
load('../R_image_files/resABC_laplace_coastal.RData')
load('../R_image_files/resABC_laplace_inland.RData')


# modelTest and Permutation Test Results ####
nhst.res = data.frame(Test=c('Exponential (both regions)','Exponential (inland)','Exponential (coastal)','Logistic (both regions)','Logistic (inland)','Logistic (coastal)','Permutation (inland vs coastal)','Permutation (millet vs all dates)'),Pval=c(exp.general.test$pval,exp.inland.test$pval,exp.coastal.test$pval,logistic.general.test$pval,logistic.inland.test$pval,logistic.coastal.test$pval,coastal.inland.permtest$pValueList[1],millet.permtest$pValueList[1]))

write.csv(nhst.res,file='nhst_result.csv',row.names=FALSE)

# Age Depth Model Posterior
point_a=kim2004.model.marine20$thetaPredict[,81]
point_b=kim2004.model.marine20$thetaPredict[,80]
age.depth.res90HPD = data.frame(Event=c('a','b'),Start=NA,End=NA)
age.depth.res90HPD$Start[1] = HPDinterval(mcmc(point_a),prob = 0.90)[2]
age.depth.res90HPD$End[1] = HPDinterval(mcmc(point_a),prob = 0.90)[1]
age.depth.res90HPD$Start[2] = HPDinterval(mcmc(point_b),prob = 0.90)[2]
age.depth.res90HPD$End[2] = HPDinterval(mcmc(point_b),prob = 0.90)[1]

write.csv(age.depth.res90HPD,file='age_depth_hpd90.csv',row.names=FALSE)


# ABC Posterior
tol=0.01 #1%
post.general = abc.general[order(abc.general$euc.uncal)[1:(nrow(abc.general)*tol)],]
post.coastal = abc.coastal[order(abc.coastal$euc.uncal)[1:(nrow(abc.coastal)*tol)],]
post.inland = abc.inland[order(abc.inland$euc.uncal)[1:(nrow(abc.inland)*tol)],]

abc.res90HPD = data.frame(Region=c('Both','Coastal','Inland'),r1_start=NA,r1_end=NA,r2_start=NA,r2_end=NA,c_start=NA,c_end=NA)

abc.res90HPD$r1_start[1]=HPDinterval(mcmc(post.general$bl),prob = 0.90)[1] * 100
abc.res90HPD$r1_end[1]=HPDinterval(mcmc(post.general$bl),prob = 0.90)[2] * 100
abc.res90HPD$r1_start[2]=HPDinterval(mcmc(post.coastal$bl),prob = 0.90)[1] * 100
abc.res90HPD$r1_end[2]=HPDinterval(mcmc(post.coastal$bl),prob = 0.90)[2] * 100
abc.res90HPD$r1_start[3]=HPDinterval(mcmc(post.inland$bl),prob = 0.90)[1] * 100
abc.res90HPD$r1_end[3]=HPDinterval(mcmc(post.inland$bl),prob = 0.90)[2] * 100

abc.res90HPD$r2_start[1]=HPDinterval(mcmc(post.general$br),prob = 0.90)[2] * -100
abc.res90HPD$r2_end[1]=HPDinterval(mcmc(post.general$br),prob = 0.90)[1] * -100
abc.res90HPD$r2_start[2]=HPDinterval(mcmc(post.coastal$br),prob = 0.90)[2] * -100
abc.res90HPD$r2_end[2]=HPDinterval(mcmc(post.coastal$br),prob = 0.90)[1] * -100
abc.res90HPD$r2_start[3]=HPDinterval(mcmc(post.inland$br),prob = 0.90)[2] * -100
abc.res90HPD$r2_end[3]=HPDinterval(mcmc(post.inland$br),prob = 0.90)[1] * -100

abc.res90HPD$c_start[1]=round(HPDinterval(mcmc(post.general$c),prob = 0.90)[2])
abc.res90HPD$c_end[1]=round(HPDinterval(mcmc(post.general$c),prob = 0.90)[1])
abc.res90HPD$c_start[2]=round(HPDinterval(mcmc(post.coastal$c),prob = 0.90)[2])
abc.res90HPD$c_end[2]=round(HPDinterval(mcmc(post.coastal$c),prob = 0.90)[1])
abc.res90HPD$c_start[3]=round(HPDinterval(mcmc(post.inland$c),prob = 0.90)[2])
abc.res90HPD$c_end[3]=round(HPDinterval(mcmc(post.inland$c),prob = 0.90)[1])

write.csv(abc.res90HPD,file='abc_result_hpd90.csv',row.names=FALSE)




