############ R CODE PART 1: ############
############ Dataset setup, filtering, and management, plus presence-absence matrix generation with letsR ############

################ Step 0.1 LOAD PACKAGES ################
# setwd("your working directory")
wd <- getwd() # main working directory that we can easily return to

packages <- c("sp","sf","dismo","raster","rgeos","rvertnet","rgdal","rworldmap","rworldxtra",
              "letsR","xlsx","dplyr","ecostructure","epm","lwgeom","tmap","terra","doParallel",
              "tidyverse","ggnetwork","vegan","ape","cluster","RColorBrewer","dendextend","NbClust",
              "factoextra","fpc","clValid","pvclust","colorspace","gclus","mapproj","viridis","cowplot")
lapply(packages, require, character.only=T)



################ Step 0.2 Set up plottable maps ################
pro <- paste("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
wgs84 <- CRS(pro)
worldmap <- getMap(resolution="high"); plot(worldmap)

# Tropical Asia Map
asia_tropic_bound <- as(extent(65, 132, -12, 37), "SpatialPolygons")
crs(worldmap)<-wgs84; crs(asia_tropic_bound) <- wgs84
asia_tropic_map <- raster::intersect(worldmap, asia_tropic_bound)
shapefile(asia_tropic_map, "Asia_Tropic_Map_2023.shp")


################ Step 0.3 Set up a plottable shapefile of Indomalaya ################
# The terrestrial ecoregions of the world can be downloaded from https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world (accessed 19 Apr 2023)
# On QGIS, we exported only ecoregions marked as part of the "Indomalayan Realm" and export it as a unique shapefile
im_ecoreg <- shapefile("Indomalaya_ecoregions_fixed.shp")
crs(im_ecoreg) <- wgs84

extent(im_ecoreg)

plot(asia_tropic_map); plot(im_ecoreg, col="#00443380", add=T)
indomalaya_map <- raster::intersect(im_ecoreg, asia_tropic_map) # trim off small mangrove features from the Arabian Peninsula
plot(indomalaya_map)

writeOGR(indomalaya_map, ".", "indomalaya_map_2023", 
         driver = "ESRI Shapefile")

# The individual features can be dissolved in QGIS (using the "Dissolve" function)
im_ecoreg_dissolved <- shapefile("IM_ecoregions_dissolved_fixed.shp")
im_map_simpl <- raster::intersect(im_ecoreg_dissolved, asia_tropic_bound)

plot(im_map_simpl, col="#BB00AA")
plot(as(extent(im_map_simpl), "SpatialPolygons"), add=T)


################ 01 Shapefile generation ################
# NOTE: Requires working in ArcGIS and QGIS in addition to R.
# On ArcGIS, we subsetted the original birdlife BOTW.gdb geodatabase into smaller batches (7000+ features per batch as opposed to 17500+).
# This was done based on scientific names that begin from A to G, then H to R, and finally S to Z.
### This subsetting was done though "select by attributes" with 
### SCINAME LIKE 'A%'
### OR SCINAME LIKE 'B%'
### OR SCINAME LIKE 'C%'
### OR SCINAME LIKE 'D%'
### (etc. etc.) for a subset of species names starting with A to D (or however far down the alphabet one chooses to go)
# Then with "select by location" and "intersect the source layer feature".
# We extracted only features (polygons) that intersect with the WWF Indomalaya shapefile.
# We then extracted the species names (SCINAME) and use this list to return to BOTW.gdb
# and extracted the species' full range polygons, using the following custom script (having listed all relevant species names in an R vector called "scinames")
# for(i in 1:length(scinames)){
#  cat("SCINAME = ")
#  cat(paste0("\'",(scinames[i] %>% str_replace("_"," ")),"\'"))
#  cat(" OR\n")
# }
# sink()
#
# We obtained 2827 unique species, which we have provided as "any_indomalaya.csv"

# These three large geodatabases can be split into individual shapefiles based on species names (or "SCINAME").
# We did this through the "Split By Attributes" function on ArcGIS, based on the field "SCINAME".

################ 02 Fix geometry (QGIS) ################
# These individual shapefiles may contain "invalid geometry" issues.
# We fixed them in QGIS, using the "fix geometries" tool, which will generate a new batch of 2827 shapefiles.


################ 03 Attributes Table ################
# The attributes tables of the dataset can be exported from ArcGIS/QGIS into Excel as "ALL_Attributes_raw.xlsx".
# Interpretations of the code are available in the Metadata document provided by BirdLife. 

# Since the geodatabase features were extracted in batches, so we combined all three spreadsheets into one,
# which can then be read into R.
ALL_attributes <- read.xlsx("ALL_Attributes_raw.xlsx", sheetName = "Sheet1")
names(ALL_attributes)
# Keep only 4 relevant columns
ALL_attributes <- ALL_attributes[,which(colnames(ALL_attributes) %in% c("SCINAME", "PRESENCE","ORIGIN", "SEASONAL"))]
View(ALL_attributes)

################ 04 Remove Non-Extant, Non-Native Polygons ################
# LIST OF SPECIES with non extant (PRESENCE != 1) features
spp.nonextant <- ALL_attributes$SCINAME[which(ALL_attributes$PRESENCE != 1)] %>% str_replace(" ","_") 
spp.nonextant <- unique(sort(spp.nonextant))

# LIST OF SPECIES with non native (ORIGIN != 1) features
spp.nonnative <- ALL_attributes$SCINAME[which(ALL_attributes$ORIGIN != 1)] %>% str_replace(" ","_")
spp.nonnative <- unique(sort(spp.nonnative))

# LIST OF SPECIES with non extant (PRESENCE != 1) features
spp.nonextant <- ALL_attributes$SCINAME[which(ALL_attributes$PRESENCE != 1)] %>% str_replace(" ","_")

# LIST OF SPECIES with non native (ORIGIN != 1) features
spp.nonnative <- ALL_attributes$SCINAME[which(ALL_attributes$ORIGIN != 1)] %>% str_replace(" ","_")

# Set shapefile directory for all extant and native species to be evaluated
path <- "/Shapefiles/ALL_Indomalaya_extnat"
filenames <- list.files(path, pattern = ".shp$") %>% str_replace(".shp","")

### REMOVING NON-EXTANT & NON-NATIVE FEATURES
sort(unique(c(spp.nonextant, spp.nonnative)))

setwd(path)

 for(i in 1:length(filenames)){
   print(paste0("Evaluating species ", i, ": ", filenames[i]))
   if(filenames[i] %in% spp.nonextant == T){
     print(paste0("Non-extant features of species ", i, " (", filenames[i], ") are being removed"))
     polygon <- shapefile(filenames[i])
     polygon <- polygon[polygon$PRESENCE==1,]
     shapefile(polygon, paste0(filenames[i],".shp"), overwrite=T)
   }
   else{
     print(paste0("Species ", i, " (", filenames[i], ") is good"))
   }
 }
 
 for(i in 1:length(filenames)){
   print(paste0("Evaluating species ", i, ": ", filenames[i]))
   if(filenames[i] %in% spp.nonnative == T){
     print(paste0("Non-native features of species ", i, " (", filenames[i], ") are being removed"))
     polygon <- shapefile(filenames[i])
     polygon <- polygon[polygon$ORIGIN==1,]
     shapefile(polygon, paste0(filenames[i],".shp"), overwrite=T)
   }
   else{
     print(paste0("Species ", i, " (", filenames[i], ") is good"))
   }
 }

setwd(wd) # Reset directory


################ 05 Remove Passage Polygons ################
spp.passage <- ALL_attributes$SCINAME[which(ALL_attributes$SEASONAL == 4)] %>% str_replace(" ","_") 

setwd(path)

for(i in 1:length(filenames)){
  print(paste0("Evaluating species ", i, ": ", filenames[i]))
  if(filenames[i] %in% spp.passage == T){
    print(paste0("Passage features of species ", i, " (", filenames[i], ") are being removed"))
    polygon <- shapefile(filenames[i])
    polygon <- polygon[polygon$SEASONAL != 4,]
    shapefile(polygon, paste0(filenames[i],".shp"), overwrite=T)
  }
  else{
    print(paste0("Species ", i, " (", filenames[i], ") is good"))
  }
}
setwd(wd)

################ 05 Remove "Uncertain" Polygons ################
spp.uncert <- ALL_attributes$SCINAME[which(ALL_attributes$SEASONAL == 5)] %>% str_replace(" ", "_")

setwd(path)

for(i in 1:length(filenames)){
  print(paste0("Evaluating species ", i, ": ", filenames[i]))
  if(filenames[i] %in% spp.uncert == T){
    print(paste0("Seasonality uncertain features of species ", i, " (", filenames[i], ") are being removed"))
    polygon <- shapefile(filenames[i])
    polygon <- polygon[polygon$SEASONAL != 5,]
    shapefile(polygon, paste0(filenames[i],".shp"), overwrite=T)
  }
  else{
    print(paste0("Species ", i, " (", filenames[i], ") is good"))
  }
}
setwd(wd)

################ 06 Area calculation ################
# At this point, all species' range map polygons should
# only contain features that are "native", "extant", 
# and with only seasonal codes "1", "2", and "3" (resident, breeding, nonbreeding)

# Create a new data frame to hold the calculations
percent_area <- data.frame(matrix(ncol = 2, nrow = length(filenames)))
colnames(percent_area) <- c("Species","Area_Presence")


for(i in 1:length(filenames)){
#for(i in1 1:12){ # small scale test
  print(paste0("Species ",i, ": ", filenames[i])) # track progress
  polygon <- shapefile(paste0(path,"/",filenames[i],".shp"))
  intersection <- gIntersection(polygon, im_map_simpl) # create an intersection polygon between species range and Indomalaya
  if (is.null(intersection) == TRUE){
    percent_area$Species[i] <- filenames[i]
    percent_area$Area_Presence[i] <- 0 # In case there are somehow non-intersection shapefiles in the directory (stops the loop from choking)
    print(paste0("Species ", i, " overlap: 0"))
  } else {
    percent_area$Species[i] <- filenames[i]
    percent_area$Area_Presence[i]<- sum(area(intersection))/sum(area(polygon)) # Calculates the percentage of the species range that is present in Indomalaya
    print(paste0("Species ", i, " overlap: ", percent_area$Area_Presence[i]))
  }
}
View(percent_area)

write.xlsx(percent_area, "percent_area.xlsx", sheetName = "percent_area", row.names = F)
write.table(percent_area, "percent_area.txt", row.names = F)

# Now we can view the number of species retained at a given minimum area presence threshold
nrow(percent_area[which(percent_area$Area_Presence >= 0.95),]) # If the minimum is 95%
nrow(percent_area[which(percent_area$Area_Presence >= 0.90),])
nrow(percent_area[which(percent_area$Area_Presence >= 0.80),])
# etc. etc.

################ 07 Visual assessment of appropriate filter threshold with maps ################
# Make a vector of species with between 50% and 95% of their total ranges present in Indomalaya
# These are the "gray area" species that needs visual inspection
forassess <- percent_area[which(percent_area$Area_Presence < 0.95 & percent_area$Area_Presence >= 0.50),]
spp.forassess <- forassess$Species
spp.forassess

setwd(path)
pdf("your_chosen_name.pdf")
for(i in 1:length(spp.forassess)){
  # for(i in 1:4){ # small scale trial
  print(paste0("species number ",i, " : ", spp.forassess[i])) # tracker
  polygon <- shapefile(paste0(spp.forassess[i],".shp"))
  scodes <- as.factor(polygon$SEASONAL) # change SEASONAL codes from characters to factors
  
  for(j in 1:length(scodes)){ # for loop to substitute SEASONAL code with real words
    levels(scodes)[j] <- seasonal[as.integer(levels(scodes)[j])]
  }
  seasoncols <- as.factor(polygon$SEASONAL)
  for(k in 1:length(seasoncols)){
    levels(seasoncols)[k] <- colorset1[as.integer(levels(seasoncols)[k])]
  }
  plot(asia_tropic_map, main=spp.forassess[i]) # plot base map
  plot(polygon, col=as.character(seasoncols), add=T) # plot species (color by factorized SEASONAL codes)
  legend("bottomleft", legend = scodes, fill = as.character(seasoncols), cex=0.9)
  mtext(paste0(forassess$Species[i], ": ", forassess$Area_Presence[i]), side = 3, cex = 0.85) # add a subheader with area presence percentage
  
}
dev.off()
setwd(wd)

################ 08 Create a unique directory to house only the 1544 selected species ################
# 60% is our decision; these species will be the main study species for all downstream analyses
above60 <- percent_area$Species[which(percent_area$Area_Presence >= 0.60)]
dir.60 <-"/Indomalaya_60percent"

setwd(path)
for(i in 1:length(above60)){
  print(paste0(i, ": ",above60[i]))
  #file.copy(paste0(above60[i],".cpg"),dir.60)
  file.copy(paste0(above60[i],".dbf"),dir.60)
  file.copy(paste0(above60[i],".prj"),dir.60)
  #file.copy(paste0(above60[i],".sbn"),dir.60)
  #file.copy(paste0(above60[i],".sbx"),dir.60)
  file.copy(paste0(above60[i],".shp"),dir.60)
  #file.copy(paste0(above60[i],".shp.xml"),dir.60)
  file.copy(paste0(above60[i],".shx"),dir.60)
}
setwd(wd)

################ 09 Import HBW checklist ################
# read in HBW checklist on which the taxonomy of the BirdLife data is based
checklist <- read.xlsx2("BirdLife_HBW_Taxonomic_Checklist_V4.xls", sheetName = "BirdLife_HBW_Taxonomic_Checklis")
View(checklist)
# create an Indomalayan-only checklist
checklist.60 <- checklist[which(checklist$Scientific_Name %in% (above60 %>% str_replace("_", " "))),]

# Check for possible pelagic species (there shouldn't be any)
orderchecks <- c("PELECANIFORMES", "SULIFORMES", 'CHARADRIIFORMES', 'ANSERIFORMES')

pdf("RangeMapPlots/pelagic_checks.pdf")
setwd(dir.60)
for(i in 1:nrow(to.check)){
  spnm <- to.check$Scientific_Name[i] %>% str_replace(" ","_")
  print(paste0("species number ",i, " : ", spnm)) # tracker
  polygon <- shapefile(paste0(spnm,".shp"))
  scodes <- as.factor(polygon$SEASONAL) # change SEASONAL codes from characters to factors
  for(j in 1:length(scodes)){ # for loop to substitute SEASONAL code with real words
    levels(scodes)[j] <- seasonal[as.integer(levels(scodes)[j])]
  }
  seasoncols <- as.factor(polygon$SEASONAL)
  for(k in 1:length(seasoncols)){
    levels(seasoncols)[k] <- colorset1[as.integer(levels(seasoncols)[k])]
  }
  plot(asia_tropic_map, main=to.check$Scientific_Name[i]) # plot base map
  plot(polygon, col=as.character(seasoncols), add=T) # plot species (color by factorized SEASONAL codes)
  legend("bottomleft", legend = scodes, fill = as.character(seasoncols), cex=0.9)
  mtext(paste0("Order: ", to.check$Order_[i], "; Family: ", to.check$Family_Name[i]), side = 3, cex = 0.85) # add a subheader
}
dev.off()
setwd(wd)

################ 10 Generate presence-absence matrix with letsR ################
# 1 degree resolution
PAM_60_1 <- lets.presab.birds(dir.60, xmn = 65, xmx = 132, ymn = -12, ymx = 37,
                              resol = 1, crs = wgs84, crs.grid = wgs84, count = T)
							  
# Visualize richness raster in a map
dev.new()
plot(PAM_60_1$Richness_Raster, main="Indomalaya species richness (1°, n=1544)",
     xlab="Longitude (°)",ylab="Latitude (°)"); plot(asia_tropic_map, border="gray68", add=T)
# plot(asia_tropic_map, border="gray", main="Indomalaya species richness (1°, n=1544)",
#      xlab="Longitude (°)",ylab="Latitude (°)"); plot(PAM_60_1$Richness_Raster, , add=T)
plot(im_ecoreg_dissolved, border="blue", add=T)
dev.off()

# 0.5 (half) degree resolution
PAM_60_half <- lets.presab.birds(dir.60, xmn = 65, xmx = 132, ymn = -12, ymx = 37,
                                 resol=0.5, crs=wgs84, crs.grid= wgs84, count=T)

# 2 degree resolution
PAM_60_2 <- lets.presab.birds(dir.60, xmn = 65, xmx = 132, ymn = -12, ymx = 37,
                                 resol=2, crs=wgs84, crs.grid= wgs84, count=T)

# 5 degree resolution
PAM_60_5 <- lets.presab.birds(dir.60, xmn = 65, xmx = 132, ymn = -12, ymx = 37,
                              resol=5, crs=wgs84, crs.grid= wgs84, count=T)


# Export presence-absence matrices as csv files 
write.csv(as.data.frame(PAM_60_1$P), "IM60_Presab_1deg.csv")
write.csv(as.data.frame(PAM_60_half$P), "IM60_Presab_halfdeg.csv")
write.csv(as.data.frame(PAM_60_2$P), "IM60_Presab_2deg.csv")
write.csv(as.data.frame(PAM_60_5$P), "IM60_Presab_5deg.csv")




