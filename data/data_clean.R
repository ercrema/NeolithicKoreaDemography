### Read 14C Data ###
koreaC14<-read.csv("./2016_Neolithic_C14_dates_collection_v9.3.csv",stringsAsFactors = FALSE)

### Clean Data ###
koreaC14 <- data.frame(EntryNo=koreaC14$Entry,
                       labcode=koreaC14$Labcode,
                       site=koreaC14$site_name_kor,
                       latitude=koreaC14$Lat,
                       longitude=koreaC14$Long,
                       deltaC13=koreaC14$X13C.0.00.,
                       c14age=koreaC14$uncal_bp,
                       c14error=koreaC14$uncal_range,
                       material=koreaC14$Material,
                       coastM=koreaC14$coast_manual, #During the Chulmun period, the sites marked with 'T' were located in the coastal area. We know this because a) these sites yielded marine shells such as oyster, b) paleocoastal construction studies at few of these sites suggested the higher sea-level during Chulmun made them coastal site, and c) archaeologists found a wooden dugout boat from one of these sites called 'Bibongri.' 
                       sitename=koreaC14$site_name_eng)

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

write.csv(koreaC14,file="./koreaC14dates.csv",row.names=F, fileEncoding="UTF-8")
