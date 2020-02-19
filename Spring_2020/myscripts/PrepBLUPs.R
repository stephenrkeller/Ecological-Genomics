setwd("~/OneDrive - University of Vermont/spruce/CASRI/TD_calc/")  #  Set your wd as appropriate

finalHT <- read.csv("height_data.csv", header=T)

require(dplyr)
require(tidyr)
require(raster)

### gather worldClim data

# Note, I already downloaded these data to my local dir, so I've set the "download=F" flag to avoid re-downloading each time.  If you try for the first time and get an error, you probably need to download the data once.

r <- getData("worldclim",var="bio",res=0.5,
             lon=-79,
             lat=39,
             download = FALSE,
             path = "./") 

require(tidyverse)
#require(tidyr)
finalHT$mLongitude = (finalHT$Longitude)*-1
finalHT$plant_ID <- as.character(finalHT$plant_ID)
coords <- data.frame(finalHT$plant_ID, finalHT$mLongitude,finalHT$Latitude)
coords <- coords[!is.na(coords[,3]),]
nrow(coords) #4405 - five rows dropped
colnames(coords) <- c("plantID","Longitude","Latitude")
str(coords)
detach("package:tidyr", unload = TRUE) # package tidyr has issues with running raster package

library(sp)

head(coords) # coordinates in decimal degrees
r@crs # to check the projection system # crs(r) - does the same function

points <- SpatialPoints(coords = coords[,2:3], proj4string = r@crs)

values <- extract(r,points)

df <- cbind.data.frame(coordinates(points),values)

head(df)

bioInfo <- cbind(coords, df)  # combining the plant ID with bio variables 
nrow(bioInfo)
head(bioInfo)
str(bioInfo)

source2 <- merge(finalHT,bioInfo,by.x="plant_ID", by.y="plantID")
nrow(source2)
summary(source2)

Ht_ALL_VTMDNC <- lmer(data = source2, Height_Fall2019 ~ Garden + initial_ht + (1|mBed) + (1|Family) )
Ht_ALL_MD <- lmer(data = source2[which(source2$Garden=="MD"),], Height_Fall2019 ~ initial_ht + (1|mBed) + (1|Family) )
Ht_ALL_NC <- lmer(data = source2[which(source2$Garden=="NC"),], Height_Fall2019 ~ initial_ht + (1|mBed) + (1|Family) )

ranef(Ht_ALL_VTMDNC)
ranef(Ht_ALL_MD)
ranef(Ht_ALL_NC)

Ht_ALL_VTMDNC_BLUPs <- ranef(Ht_ALL_VTMDNC)$Family
Ht_ALL_VTMDNC_BLUPs$Ind = rownames(Ht_ALL_VTMDNC_BLUPs)
Ht_ALL_MD_BLUPs <- ranef(Ht_ALL_MD)$Family
Ht_ALL_MD_BLUPs$Ind = rownames(Ht_ALL_MD_BLUPs)
Ht_ALL_NC_BLUPs <- ranef(Ht_ALL_NC)$Family
Ht_ALL_NC_BLUPs$Ind = rownames(Ht_ALL_NC_BLUPs)

EDGE <- read.table("EDGE.annot")
str(EDGE)

Ht_EDGE_VTMDNC_BLUPs <- Ht_ALL_VTMDNC_BLUPs[which(as.factor(Ht_ALL_VTMDNC_BLUPs$Ind) %in% EDGE$V1==TRUE),]
Ht_EDGE_MD_BLUPs <- Ht_ALL_MD_BLUPs[which(as.factor(Ht_ALL_MD_BLUPs$Ind) %in% EDGE$V1==TRUE),]
Ht_EDGE_NC_BLUPs <- Ht_ALL_NC_BLUPs[which(as.factor(Ht_ALL_NC_BLUPs$Ind) %in% EDGE$V1==TRUE),]

plot(Ht_EDGE_MD_BLUPs[,1],Ht_EDGE_NC_BLUPs[,1])
abline(0,1)

write.table(Ht_EDGE_VTMDNC_BLUPs[,1],"Ht_EDGE_VTMDNC.blups",row.names = F,col.names = F,quote = F)
write.table(Ht_EDGE_MD_BLUPs[,1],"Ht_EDGE_MD.blups",row.names = F,col.names = F,quote = F)
write.table(Ht_EDGE_NC_BLUPs[,1],"Ht_EDGE_NC.blups",row.names = F,col.names = F,quote = F)

EDGE2 = droplevels(source2[which(source2$Region=="E"),])

EDGE_df <- aggregate(EDGE2,by=list(EDGE2$Family),FUN="mean")

EDGE_df2 <- EDGE_df[which(EDGE_df$Group.1 %in% EDGE$V1),]

EDGE_df3 = EDGE_df2[,c(1,10,24,25,32,45)]

write.table(EDGE_df3,"Bioclim_EDGE.blups",row.names = F,col.names = T,quote = F)

