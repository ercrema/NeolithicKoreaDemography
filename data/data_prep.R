library(rcarbon)
library(sp)
library(dbscan)
library(rworldmap)
library(maptools)
library(rgdal)
library(rgeos)

### General Settings ###
coastalThreshold = 2000 #in m
clusteringThreshold = 1 #in km
temporalBinSize = 100 #in years

### Read 14C Data ####
koreaC14<-read.csv("./koreaC14dates.csv")

### Read Coastal Data from rworldmap
basemap <- getMap(resolution = "high")
koreanPeninsula=subset(basemap,SOVEREIGNT%in%c('South Korea','North Korea','Korea No Mans Area'))
koreanPeninsula <- unionSpatialPolygons(koreanPeninsula,IDs=c(1,1,1))
koreanPeninsulaCoast = as(koreanPeninsula, "SpatialLines")

#### Use DBSCAN to Clusters site in proximity #### 
sites <- unique(data.frame(SiteID=koreaC14$SiteID,latitude=koreaC14$latitude,longitude=koreaC14$longitude))
sites.sp = sites
sp::coordinates(sites.sp) <- c("longitude","latitude")
sp::proj4string(sites.sp) <- sp::CRS("+proj=longlat +datum=WGS84")
D <- sp::spDists(sites.sp,sites.sp,longlat=TRUE)
eps=clusteringThreshold
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

### Compute Distance to Coast from Cluster Centers
clusters = koreaC14
sp::coordinates(clusters) <- c("clon","clat")
sp::proj4string(clusters) <- sp::CRS("+proj=longlat +datum=WGS84")
clusters_utm <- spTransform(clusters,CRS("+proj=utm +zone=52 ellps=WGS84")) # Project to UTM
koreanPeninsulaCoast_utm <- spTransform(koreanPeninsulaCoast,CRS("+proj=utm +zone=52 ellps=WGS84")) # Project to UTM
koreanPeninsula_utm <- spTransform(koreanPeninsula,CRS("+proj=utm +zone=52 ellps=WGS84"))
distfromcoast <- gDistance(koreanPeninsulaCoast_utm,clusters_utm,byid=TRUE)
koreaC14$withinPoly=gContains(koreanPeninsula_utm,clusters_utm,byid=TRUE)
koreaC14$coastDist = as.numeric(distfromcoast)
koreaC14$region = 'inland'
koreaC14$region[which(koreaC14$coastDist<coastalThreshold)]='coastal'


## test plot:
# clusters_utm$region=koreaC14$region
# clusters_utm$coastDist=koreaC14$coastDist
# clusters_utm$withinPoly=koreaC14$withinPoly
# plot(koreanPeninsulaCoast_utm)
# points(clusters_utm,pch=20,cex=clusters_utm$coastDist/50000,col=rgb(1,0,0,0.5))
# points(clusters_utm,pch=20,cex=0.6,col=as.numeric(as.factor(clusters_utm$region)))
# points(clusters_utm,pch=20,cex=0.6,col=as.numeric(clusters_utm$withinPoly)+2)

#### Calibration ####
caldates <- calibrate(koreaC14$c14age,koreaC14$c14error)
bins=binPrep(ages=caldates,h=temporalBinSize,sites=koreaC14$cluster)

### CKDE Random Bootstrapped Samples ####
set.seed(123)
sdates = sampleDates(caldates, bins = bins, nsim=1000, boot = TRUE, verbose = TRUE)

#### Store in R image ####
save(sdates,caldates,bins,koreaC14,file='./koreanC14.RData')

