### Read 14C Data ###
koreaC14<-read.csv("./2016_Neolithic_C14_dates_collection_v9.2.csv",stringsAsFactors = TRUE)

### Clean Data ###
koreaC14 <- data.frame(EntryNo=koreaC14$Entry,
		       labcode=koreaC14$Labcode,
		       site=koreaC14$Site.name,
		       latitude=koreaC14$Lat,
		       longitude=koreaC14$Long,
		       deltaC13=koreaC14$X13C.0.00.,
		       c14age=koreaC14$uncal_bp,
		       c14error=koreaC14$uncal_range,
		       material=koreaC14$Material)

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



koreaC14$SiteID <- as.numeric(koreaC14$site)
koreaC14$latitude<-as.numeric(as.character(koreaC14$latitude))
koreaC14$longitude<-as.numeric(as.character(koreaC14$longitude))


problems.id=c()

for (x in 1:length(unique(koreaC14$SiteID)))
{
ref=unique(koreaC14$SiteID)[x]	
tmp=subset(koreaC14,SiteID==ref)
if(length(unique(tmp$latitude))>1|length(unique(tmp$longitude))>1)
{
problems.id=c(problems.id,ref)
}	

}	



# Remove Dates without Lat/Long
koreaC14 <- subset(koreaC14,!is.na(latitude)&!is.na(longitude))

write.csv(koreaC14,file="./koreaC14dates.csv",row.names=F, fileEncoding="UTF-8")
