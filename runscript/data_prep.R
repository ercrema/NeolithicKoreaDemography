library(rcarbon)
library(sp)
library(dbscan)
library(rworldmap)
library(maptools)
library(rgdal)
library(rgeos)
library(raster)

### General Settings ###
coastalThreshold = 2000 #in m
clusteringThreshold = 1 #in km
temporalBinSize = 100 #in years

### Read 14C Data ####
koreaC14<-read.csv("../data/Neolithic_C14_dates_cleaned.csv")

### Read Coastal Data
koreanPeninsulaCoast = shapefile("../data/shp/polyline_korea.shp") 
koreanPeninsulaCoast = gLineMerge(koreanPeninsulaCoast)

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

### Compute Distance to Coast from Sites
sitepoints = koreaC14
sp::coordinates(sitepoints) <- c("longitude","latitude")
sp::proj4string(sitepoints) <- sp::CRS("+proj=longlat +datum=WGS84")
sitepoints_utm <- spTransform(sitepoints,CRS("+proj=utm +zone=52 ellps=WGS84")) # Project to UTM
koreanPeninsulaCoast_utm <- spTransform(koreanPeninsulaCoast,CRS("+proj=utm +zone=52 ellps=WGS84")) # Project to UTM
distfromcoast <- gDistance(koreanPeninsulaCoast_utm,sitepoints_utm,byid=TRUE)
koreaC14$coastDistSite = as.numeric(distfromcoast)
koreaC14$region_site = 'inland'
koreaC14$region_site[which(koreaC14$coastDistSite<coastalThreshold)]='coastal'
koreaC14$region_site[which(koreaC14$coastM == TRUE)]='coastal' #Manual adjustment for sites mislabelled as inland. Individual check of the site reports indicated that these sites were most likely located in the coastal area as they either yielded marine shells or wooden dugout boat (Bibongri). 

### Compute Distance to Coast from Cluster Centres
clusters = koreaC14
sp::coordinates(clusters) <- c("clon","clat")
sp::proj4string(clusters) <- sp::CRS("+proj=longlat +datum=WGS84")
clusters_utm <- spTransform(clusters,CRS("+proj=utm +zone=52 ellps=WGS84")) # Project to UTM
distfromcoast <- gDistance(koreanPeninsulaCoast_utm,clusters_utm,byid=TRUE)
koreaC14$coastDistClusters = as.numeric(distfromcoast)
koreaC14$region_clusters = 'inland'
koreaC14$region_clusters[which(koreaC14$coastDistClusters<coastalThreshold)]='coastal'

### Check Inconsistencies between cluster based region and site based regions
inconsistent_clusters=koreaC14[which(koreaC14$region_clusters!=koreaC14$region_site),]$cluster
inconsistent_cases =unique(dplyr::select(subset(koreaC14,cluster%in%c(inconsistent_clusters)),SiteID,sitename,latitude,longitude,clat,clon,coastDistSite,region_site,coastM,coastDistClusters,region_clusters,cluster))  
print(inconsistent_cases)
# All Inconsistencies (for coastalThreshold = 2000 meters and clusteringThreshold = 1km) are clusters with comprising a single site with manual adjustment (see line 55 above). 

#Repair inconsistencies
koreaC14$region_clusters=koreaC14$region_site

#Select Relevant fields for analyses
koreaC14 = dplyr::select(koreaC14,labcode,c14age,c14error,material,deltaC13,site_id=SiteID,site_kor=site,site_en=sitename,cluster_id=cluster,latitude,longitude,clat,clon,milletAsso=milletAsso,region=region_clusters)

#### Calibration ####
caldates <- calibrate(koreaC14$c14age,koreaC14$c14error)
bins=binPrep(ages=caldates,h=temporalBinSize,sites=koreaC14$cluster_id)

#### Store in R image ####
save(caldates,bins,koreaC14,file='../R_image_files/koreanC14.RData')

