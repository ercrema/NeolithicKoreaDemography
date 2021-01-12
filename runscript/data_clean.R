### Load Required R Packages
library(dplyr)
library(here)

### Read 14C Data ###
koreaC14<-read.csv(here('data','Neolithic_C14_dates.csv'),stringsAsFactors = FALSE)

### Clean Data ###
koreaC14 <- koreaC14 %>%mutate(milletAsso= case_when(
  stringr::str_detect(seed_association, "FM") | stringr::str_detect(seed_association, "BM") | 
    stringr::str_detect(Materail_cmt, "millet") | stringr::str_detect(Materail_cmt, "Millet")  ~ TRUE))

koreaC14 <- data.frame(EntryNo=koreaC14$Entry,
                       labcode=koreaC14$Labcode,
                       site=koreaC14$site_name_kor,
                       latitude=koreaC14$Lat,
                       longitude=koreaC14$Long,
                       deltaC13=koreaC14$X13C.0.00.,
                       c14age=koreaC14$uncal_bp,
                       c14error=koreaC14$uncal_range,
                       material=koreaC14$Material,
                       coastM=koreaC14$coast_manual, 
                       sitename=koreaC14$site_name_eng,
                       milletAsso=koreaC14$milletAsso)

# Data Check #
rownames(koreaC14)=koreaC14$EntryNo
koreaC14<-koreaC14[,-1]

# Remove cases without 14C Dates
koreaC14<-koreaC14[which(!is.na(koreaC14$c14age)),]
# Convert blank space in NA for labcodes:
koreaC14$labcode[which(koreaC14$labcode=="")]=NA

## Remove Marine, no LabCode, and 14C Age between start.date and end.date
start.date=6200
end.date=2800


koreaC14 <- subset(koreaC14,
                   !is.na(labcode)&
                     material%in%c("charcoal",
                                   "charcoal (annual)",
                                   "wood")&
                     c14age<start.date&c14age>=end.date)




# Convert Lat/Lon in numeric
koreaC14$latitude<-as.numeric(as.character(koreaC14$latitude))
koreaC14$longitude<-as.numeric(as.character(koreaC14$longitude))


#Assign Site ID to unique SiteName-Lat-Lon Combination
siteList = unique(data.frame(Site=koreaC14$sitename,Latitude=koreaC14$latitude,Longitude=koreaC14$longitude,stringsAsFactors = FALSE))

# Site Names with different coordinates
subset(siteList,Site%in%c(names(which(table(siteList$Site)>1))))

siteList$SiteID = 1:nrow(siteList)
koreaC14$SiteID = NA
for (i in 1:nrow(koreaC14))
{
  k=which(siteList$Site==koreaC14$sitename[i]&siteList$Latitude==koreaC14$latitude[i]&siteList$Longitude==koreaC14$longitude[i])
  koreaC14$SiteID[i]=siteList$SiteID[k]
}



# Remove Dates without Lat/Long
koreaC14 <- subset(koreaC14,!is.na(latitude)&!is.na(longitude))

write.csv(koreaC14,file=here('data','Neolithic_C14_dates_cleaned.csv'),row.names=F, fileEncoding="UTF-8")

          