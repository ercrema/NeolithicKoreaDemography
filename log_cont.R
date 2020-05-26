library(rcarbon)
load("temp1.RData")
#### Model Fitting and Hypothesis Testing #### 


## Exponential Model ##
set.seed(123)
expFit<-lm(log(PrDens)~calBP,data=spdKorea$grid)
expFitDens=data.frame(calBP=spdKorea$grid$calBP,PrDens=exp(predict(expFit)))
ExpCheck <- modelTest(x=caldates,bins=bins,errors=koreaC14$c14error,model="custom",runm=200,nsim=5000,ncores=3,timeRange=c(7000,3000),predgrid=expFitDens)
#plot(spdKorea)
#lines(7000:3000,exp(ExpCheck$coefficients[1]+c(7000:3000)*ExpCheck$coefficients[2]),col=2,lty=2)
save.image("temp2.RData")



## Logistic Model ##
x <- spdKorea$grid$calBP
y <- spdKorea$grid$PrDens
logFit <- nls(y~SSlogis(x, Asym, xmid, scale),control=nls.control(maxiter=200),start=list(Asym=0.2,xmid=5500,scale=-100))

#plot(spdKorea)
#lines(x,SSlogis(x,coefficients(logFit)[1],coefficients(logFit)[2],coefficients(logFit)[3]),col=2,lty=2)
set.seed(123)
logFitDens=data.frame(calBP=x,PrDens=SSlogis(x,coefficients(logFit)[1],coefficients(logFit)[2],coefficients(logFit)[3]))
LogCheck <- modelTest(bins=bins,errors=koreaC14$c14error,model="custom",runm=200,nsim=5000,ncores=3,x=caldates,timeRange=c(7000,3000),predgrid=logFitDens)
save.image("temp3.RData")


## Double Exponential Model ##
laplace <- function(x, a, bl, br, c) {
    ifelse(x < c, a * exp((x-c) / br), a * exp(-(x-c) / bl))
}
lapFit <- nls(y ~ laplace(x, a, bl, br, c), start=list(a=0.3, bl=500, br=500, c=5000))
#plot(spdKorea)
#lines(x,laplace(x, coefficients(lapFit)[1], coefficients(lapFit)[2], coefficients(lapFit)[3],coefficients(lapFit)[4]),col=2,lty=2)
#growth rates are computed from 1/coefficients(lapFit)[2] and 1/coefficients(lapFit)[3]

set.seed(123)
lapFitDens=data.frame(calBP=x,PrDens=laplace(x, coefficients(lapFit)[1], coefficients(lapFit)[2], coefficients(lapFit)[3],coefficients(lapFit)[4]))
LapCheck <- modelTest(bins=bins,errors=koreaC14$c14error,model="custom",runm=200,nsim=5000,ncores=3,x=caldates,timeRange=c(7000,3000),predgrid=lapFitDens)
save.image("temp4.RData")



#### CI of Laplace (bin bootstrap approach) ####
bls<-numeric()
brs<-numeric()
timeRange=c(7000,3000)
nsim=5000
midPT<-numeric()
obs.bins=bins
calBP=spdKorea$grid$calBP




for (x in 1:nsim)
{
     print(x)	
     boot.bins<-sample(unique(obs.bins),replace=TRUE)
     k=0
     obins <- numeric()
     newbins <- numeric()
     index <- numeric()
     for (j in 1:length(boot.bins))
	{
	  k=k+1
	  i<-which(obs.bins==boot.bins[j])
          index=c(index,i)
	  obins=c(obins,obs.bins[i])
	  newbins=c(newbins,rep(k,length(i)))
	}
     boot.caldates=caldates[index]
     boot.spd <- spd(x=boot.caldates,bins=newbins,timeRange=timeRange,spdnormalised=FALSE,verbose=F,runm=200)
     pr = boot.spd$grid$PrDens
	maxPr=max(pr)
     
     whichPeak=which(pr==maxPr)
     peak=calBP[whichPeak]
	up=lm(log(pr[1:whichPeak]+0.0001)~calBP[1:whichPeak])
	up=1/abs(coef(up)[2])
	down=lm(log(pr[-c(1:whichPeak)]+0.0001)~calBP[-c(1:whichPeak)])
	down=1/abs(coef(down)[2])
        lapFittmp <- tryCatch(nls(pr ~ laplace(x=calBP, a, bl, br, c), start=list(a=maxPr, bl=up, br=down, c=peak)),error=function(e){NA})
	bls[x]=NA
	brs[x]=NA
	if (length(lapFittmp)>1)
	{
	bls[x]=1/coefficients(lapFittmp)[2]
	brs[x]=1/coefficients(lapFittmp)[3]
	midPT[x]=coefficients(lapFittmp)[4]
	}
}



save.image("~/Dropbox/PrehistoricPopulationKorea/R/results.RData")





# 
# 
# #### CI of Laplace Model (smooth bootstrap approach) #### 
# 
# bls<-numeric()
# brs<-numeric()
# calBP=spdKorea$grid$calBP
# samplesize=length(unique(bins))
# errors= caldates$metadata$Error
# predgrid <- spdKorea$grid
# cragrid <- uncalibrate(as.CalGrid(predgrid), verbose=FALSE)
# obscras <- caldates$metadata$CRA
# calCurves <-  caldates$metadata$CalCurve
# cragrid$PrDens[cragrid$CRA > max(obscras) | cragrid$CRA < min(obscras)] <- 0
# timeRange=c(7000,3000)
# nsim=5000
# midPT<-numeric()
# set.seed(123)
# 
# for (x in 1:nsim)
# {
# 	print(x)
# 	randomDates <- sample(cragrid$CRA, replace=TRUE, size=samplesize, prob=cragrid$PrDens)
#         randomSDs <- sample(size=length(randomDates), errors, replace=TRUE)
#         tmp <- calibrate(ages=randomDates,errors=randomSDs, resOffsets=0 ,resErrors=0, timeRange=timeRange, calCurves='intcal13', ncores=3, verbose=FALSE, calMatrix=TRUE)
#         simDateMatrix <- tmp$calmatrix
#         prtemp<- apply(simDateMatrix,1,sum)
#         pr <- runMean(prtemp, 200, edge="fill")
# 	maxPr=max(pr)
# 	whichPeak=which(pr==maxPr)
# 	peak=calBP[whichPeak]
# 	up=lm(log(pr[1:whichPeak]+0.0001)~calBP[1:whichPeak])
# 	up=1/abs(coef(up)[2])
# 	down=lm(log(pr[-c(1:whichPeak)]+0.0001)~calBP[-c(1:whichPeak)])
# 	down=1/abs(coef(down)[2])
#         lapFittmp <- tryCatch(nls(pr ~ laplace(x=calBP, a, bl, br, c), start=list(a=maxPr, bl=up, br=down, c=peak)),error=function(e){NA})
# 	bls[x]=NA
# 	brs[x]=NA
# 	if (length(lapFittmp)>1)
# 	{
# 	bls[x]=1/coefficients(lapFittmp)[2]
# 	brs[x]=1/coefficients(lapFittmp)[3]
# 	midPT[x]=coefficients(lapFittmp)[4]
# 	}
# 
# }
# save.image("temp5.RData")





















save.image("~/Dropbox/PrehistoricPopulationKorea/R/results.RData")
