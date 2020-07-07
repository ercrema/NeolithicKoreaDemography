library(rcarbon)
library(sp)
library(dbscan)

### Read 14C Data ####
koreaC14<-read.csv("./koreaC14dates.csv")

#### Use DBSCAN to Clusters site in proximity #### 
sites <- unique(data.frame(SiteID=koreaC14$SiteID,latitude=koreaC14$latitude,longitude=koreaC14$longitude))
sites.sp = sites
sp::coordinates(sites.sp) <- c("longitude","latitude")
sp::proj4string(sites.sp) <- sp::CRS("+proj=longlat +datum=WGS84")
D <- sp::spDists(sites.sp,sites.sp,longlat=TRUE)
eps=1 #eps set to 1 km
clust<-dbscan(as.dist(D),eps=eps,minPts=1)$clust
sites$cluster<-clust
clusterNames=unique(clust)
sites$clat=NA
sites$clon=NA

for (x in 1:length(clusterNames))
{
  i <- which(sites$cluster==clusterNames[x])
  sites$clat[i] = mean(sites$latitude[i])
  sites$clon[i] = mean(sites$longitude[i])
}
sites <- data.frame(SiteID=sites$SiteID,cluster=sites$cluster,clon=sites$clon,clat=sites$clat)

koreaC14 = unique(merge(x=koreaC14,y=sites,by.x="SiteID",by.y="SiteID",x.all=TRUE))

#### Calibration ####
caldates <- calibrate(koreaC14$c14age,koreaC14$c14error)
bins=binPrep(ages=caldates,h=100,sites=koreaC14$cluster)

#### Store in R image ####
save(caldates,bins,koreaC14,file='./koreanC14.RData')

