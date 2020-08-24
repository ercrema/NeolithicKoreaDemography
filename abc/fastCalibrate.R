fastCalibrate <- function(x, errors){
  require(rcarbon)  
  timeRange=c(55000,0)
  eps=1e-5
  ids <- as.character(1:length(x))
  
  calCurveFile <- paste(system.file("extdata", package="rcarbon"), "/intcal20.14c", sep="")
  options(warn=-1)
  calcurve <- readLines(calCurveFile, encoding="UTF-8")
  calcurve <- calcurve[!grepl("[#]",calcurve)]
  calcurve.con <- textConnection(calcurve)
  calcurve <- as.matrix(read.csv(calcurve.con, header=FALSE, stringsAsFactors=FALSE))[,1:3]
  close(calcurve.con)
  options(warn=0)
  colnames(calcurve) <- c("CALBP","C14BP","Error")
  calBP = calcurve[, 1]
  c14BP = calcurve[, 2]
  calSd = calcurve[, 3]
  calBPrange = timeRange[1]:timeRange[2]
  mu = stats::approx(calBP, c14BP, xout = calBPrange, rule = 2)$y
  tau1 = stats::approx(calBP, calSd, xout = calBPrange,rule = 2)$y
  
  tempList = vector('list',length=length(x))
  
  for (i in 1:length(x)) {
    tau = errors[i]^2 + tau1^2
    dens = dnorm(x[i], mean=mu, sd=sqrt(tau))
    dens[dens < eps] <- 0
    dens <- dens/sum(dens)
    dens[dens < eps] <- 0
    dens <- dens/sum(dens)
    tempList[[i]] = list(dens=dens[dens > eps],calBP = calBPrange[dens > eps])
    }

    reslist <- vector(mode="list", length=2)
    sublist <- vector(mode="list", length=length(x))
    names(sublist) <- names(x)
    names(reslist) <- c("metadata","grids")
    ## metadata
    df <- as.data.frame(matrix(ncol=11, nrow=length(x)), stringsAFactors=FALSE)
    names(df) <- c("DateID","CRA","Error","Details","CalCurve","ResOffsets","ResErrors","StartBP","EndBP","Normalised","CalEPS")
    df$DateID <- ids
    df$CRA <- x
    df$Error <- errors
    df$CalCurve <- 'intcal20'
    df$ResOffsets <- NA
    df$ResErrors <- NA
    df$StartBP <- 55000
    df$EndBP <- 0
    df$Normalised <- TRUE
    reslist[["metadata"]] <- df
    ## grids
    for (i in 1:length(x)){
      tmp <- x[[i]]
      res <- data.frame(calBP=tempList[[i]][[2]],PrDens=tempList[[i]][[1]])
      class(res) <- append(class(res),"calGrid")        
      sublist[[i]] <- res
    }
    reslist[["grids"]] <- sublist
    reslist[["calmatrix"]] <- NA
    class(reslist) <- c("CalDates",class(reslist))
    return(reslist)
  }


