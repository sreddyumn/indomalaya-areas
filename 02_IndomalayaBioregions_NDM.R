############ R CODE PART 2: ############
############ NDM analyses ############
### NOTE: This script may involve R objects from "Part 1"

# setwd("your working directory")
wd <- getwd() # main working directory that we can easily return to

packages <- c("sp","sf","dismo","raster","rgeos","rvertnet","rgdal","rworldmap","rworldxtra",
              "letsR","xlsx","dplyr","ecostructure","epm","lwgeom","tmap","terra","doParallel",
              "tidyverse","ggnetwork","vegan","ape","cluster","RColorBrewer","dendextend","NbClust",
              "factoextra","fpc","clValid","pvclust","colorspace","gclus","mapproj","viridis","cowplot")
lapply(packages, require, character.only=T)


################ 01 NDM input file setup ################
getwd()

##### Note: this section's code needs to be repeated a total of four times, one for each analytical resolution
# Read in the presence-absence matrix csv file generated from letsR
letsR.PA <- read.csv(file.choose()) # manually select relevant file
str(letsR.PA)
summary(letsR.PA)
head(letsR.PA)

### Making a hash ###
# vectorize assign, get and exists functions for convenience
assign_hash <- Vectorize(assign, vectorize.args = c("x", "value"))
get_hash <- Vectorize(get, vectorize.args = "x")
exists_hash <- Vectorize(exists, vectorize.args = "x")

# Establish longitudes keys and associated integer values
longs_key <- as.character(sort(unique(letsR.PA$Longitude.x.), decreasing = F)) # Non-decreasing longitude goes west to east
longs_value <- c(0:(length(longs_key)-1))

# Establish latitudes keys and associated integer values
lats_key <- as.character(sort(unique(letsR.PA$Latitude.y.), decreasing = T)) # Decreasing latitude goes north to south
lats_value <- c(0:(length(lats_key)-1))

## Initialize a hash ##
hash <- new.env(hash = TRUE, parent = emptyenv(), size = 400L)

### Longitudes hash (column 2) ###
assign_hash(longs_key, longs_value, hash)
# Test to see if it works
get_hash(as.character(letsR.PA[5:10,2]), hash) # Extract some random longitude values from original pres-abs matrix

### Test Longs For Loop ###
# Make a for loop to see if I can extract integer values based on Longitude keys
testLongVector <- c()
for (i in 1:length(letsR.PA$X)){
  testLongVector[i]<-(get_hash(as.character(letsR.PA$Longitude.x.[i]), hash))
}
# Test to see if the value assignments worked
testLongVector[1:19]
get_hash(as.character(letsR.PA$Longitude.x.[1:19]),hash)

### Latitudes hash (column 3) ###
assign_hash(lats_key, lats_value, hash)
# Test to see if the value assignments worked
get_hash(as.character(letsR.PA[c(3,100:102,200:202,300),3]), hash) # Extract some random latitude values from original pres-abs matrix

### Test Longs For Loop ###
# Do the same for loop test for latitudes too
testLatVector<-c()
for (i in 1:length(letsR.PA$X)){
  testLatVector[i]<-(get_hash(as.character(letsR.PA$Latitude.y.[i]), hash))
}
testLatVector[90:145]
get_hash(as.character(letsR.PA$Latitude.y.[90:145]), hash)

# Crete a matrix of just presence-absences (a matrix of numbers for each species/column)
presab_matrix <- as.matrix(letsR.PA[1:length(letsR.PA$X),4:length(letsR.PA)]) 
head(presab_matrix)
as.numeric(presab_matrix[1,]) # Check to see it works

# Create an object resol to calculate resolution for the input file
sort(unique(letsR.PA$Longitude.x.), decreasing = T)
resol <- sort(unique(letsR.PA$Longitude.x.), decreasing = T)[1] - sort(unique(letsR.PA$Longitude.x.), decreasing = T)[2]; resol

# Create an object SppN (species number) for the input file
SppN<-length(names(letsR.PA[-c(1:3)])); SppN


### The below lines create the input file ###
resolstring <- XXX # manually input the resolution string as a character string

NDM.filename <- paste0("NDM/IM60_",resolstring,"deg.dat")
sink(NDM.filename)
cat(paste("spp",SppN, sep = " ")) # This writes the species number line
cat("\n")
cat(paste("rows",length(lats_key), sep = " ")) # This writes the number of rows in the matrix
cat("\n")
cat(paste("cols",length(longs_key), sep = " ")) # This writes the number of lines in the matrix
cat("\n")
cat(paste("gridx",min(letsR.PA$Longitude.x.), resol, sep=" ")) # This writes the longitudinal point of reference 
cat("\n")
cat(paste("gridy",max(letsR.PA$Latitude.y.), resol, sep=" ")) # This writes the latitudinal point of reference
cat("\ndata\n") # The following bit is the actual matrix
for(i in 1:length(letsR.PA$X)){
  cat(paste("A",testLatVector[i],"-",testLongVector[i], sep=''))
  cat("\t")
  
  cat(as.numeric(presab_matrix[i,]), sep="")
  
  cat("\n")
}
sink()

##### !! Rinse and repeat until files for all resolutions are written

################ 02 Run analyses on NDM/VNDM ################
# This needs to be done externally with the VNDM GUI
# After importing the .dat input file, one can run the analysis with the "Analyze with NDM" function,
# 	which prompts a dialog box where one can adjust settings.
# The results can be saved into a text file by the "create file for output" function (as a .out file)
# Area images can be saved as metafiles.
# All consensus areas' coordinates can be saved as a txt file and visualized in DIVAGIS with the "Save coordinates for GIS" function.
# 	The coordinates can be exported as shapefiles from the DIVAGIS software.
# The entire run can be saved as a log file (.ndm format) with the "Save to Log file" function.

################ 03 Run NDM Species ID ################
# Because NDM processes species names only as integers starting with 0, 
# we need a way to identify any given species of interest in the results based on its NDM species id
ndmwd <- "~/NDM"
### NDM-associated id lookup ###
n <- c(0:(length(above60)-1)) # A vector of integers starting from 0 to end of 1544 species list
species_names <- above60

# Create a dataframe where I can look up the species name based on the NDM species ID
NDM_species <- as.data.frame(cbind(n, species_names))
NDM_species$n <- as.numeric(NDM_species$n) # turn number column into numbers
class(NDM_species$n);class(NDM_species$species_names) # Check class
write.xlsx2(NDM_species, "NDM/NDM_species_60percent.xlsx", sheetName = "NDM_species", row.names = F)


################ 03 Manipulate .OUT files from NDM to extract information on consensus areas and endemics ################
library("stringr")
library("xlsx")

# Set working directory as the NDM directory
setwd(ndmwd)

# SET NAME of the specific NDM run for file extraction
NDMrun <- "IM60_halfdeg_2percent" # Manually change character string to subdirectory containing the specific analysis of interest
								  # in this case we're examining the results of the 0.5 degree resolution, 2 percent consensus areas

# build character string of filename 
outfile <- list.files()[grep(paste0(NDMrun,"_DATA.out"),list.files())]

# Extract all lines from the .out file
outfile.lines <- readLines(outfile)
head(outfile.lines, 20)

# Check to see if searching for "Consensus" will successfully give the relevant lines
# Detect (and give index of) strings with the substring "Consensus"
outfile.lines[str_detect(outfile.lines, "Consensus")] 
str_extract(outfile.lines[str_detect(outfile.lines, "Consensus")], regex("^Consensus area {1}(\\d+)"))
# regex translation: ^ anchors "Consensus area " as the start, {1} means there's only 1 space, 
# and \\d for any digits, + means one or more (digits)
# This yields substrings "Consensus area X", where XX can be 1- or 2-digits

# the line numbers where all Consensus Area's endemic species info starts (where it says "Consensus area XX (from XX areas), ... give score:")
Consensus_indices_start <- grep("Consensus", outfile.lines) 
# the line numbers where all Consensus Area's endemic species info ends (where it says "Areas included: ...")
Consensus_indices_end <- grep("Areas included", outfile.lines) 

# Create an empty list to contain dataframes of all consensus areas and its endemic species
Consensus_area_endemics <- list()

for(i in 1:length(Consensus_indices_start)){
  subset <- outfile.lines[Consensus_indices_start[i]:(Consensus_indices_end[i]+1)] # Limit search within individual consensus areas
  
  # Pull out strings with endemic species information (given the regex below)
  endemics <- subset[grep(regex("^ *(\\d+) +\\({1}(\\d+)\\){1}:{1} \\({1}.*", multiline=T), subset)]
  # Regex translation:
  # ^ * --> any number (even 0) of spaces at the beginning of the string
  # (\\d+) --> any digit, occurring once or more
  #  + --> one or more spaces
  # \\({1} --> exactly one open parenthesis symbol (double backslashes \\ makes it literal)
  # (\\d+) --> same thing: any digit, occurring once or more
  # \\){1} --> exactly one close parenthesis symbol (double backslashes \\ makes it literal)
  #  \\({1} --> exactly one space-open parenthesis combination
  # .* --> any character (period), for any number of times (asterisk)
  
  # pull out only the numbers at the start of the strings, before the first open parenthesis,
  # and make them numeric
  endemics <- as.numeric(str_extract(endemics, (regex("(\\d+){1}", multiline=T))))
  
  # From the premade NDM_species data frame,
  # subset based on the specific vector of numeric codes to pull out the correct species names
  Consensus_area_endemics[[i]] <- NDM_species[which(NDM_species$n %in% endemics),]
}

# Rename the list elements with "Consensus area X"
names(Consensus_area_endemics) <- str_extract(outfile.lines[str_detect(outfile.lines, "Consensus")], regex("^Consensus area {1}(\\d+)"))


Consensus_area_endemics

lapply(Consensus_area_endemics, function(x) write.table( data.frame(x), 
                                                         paste0(NDMrun, "_spList.csv"), 
                                                         append= T, sep=',' ))
write.xlsx(Consensus_area_endemics[[1]], paste0(NDMrun, "_spList.xlsx"), 
           sheetName=names(Consensus_area_endemics[1]), row.names = F)
for(i in 2:length(Consensus_area_endemics)){
  print(names(Consensus_area_endemics[i]))
  write.xlsx(Consensus_area_endemics[[i]], paste0(NDMrun, "_spList.xlsx"),
             sheetName=names(Consensus_area_endemics[i]), row.names = F, append = T)  
}

# Reset working directory before returning to the main script
setwd(wd)


######### 04 Group consensus areas geographically for ease of plotting #########
#### 1 degree
AOE1deg.dir <- ("NDM/IM60_1deg_AOE_2pct")
setwd(AOE1deg.dir)

AOEs <- list.files(pattern = ".shp$")
AOEs1deg <- list()
AOEs1deg[[1]] <- AOEs[c(30, 6, 21)] # Indian subcontinent
AOEs1deg[[2]] <- AOEs[c(27, 24, 1, 22, 3, 17)] # Philippines
AOEs1deg[[3]] <- AOEs[c(31, 25)] # Andaman-Nicobar
AOEs1deg[[4]] <- AOEs[c(13, 26, 2, 19, 7, 20, 10, 5, 4, 8)] # Sundaland
AOEs1deg[[5]] <- AOEs[c(23, 15, 12, 18, 11, 13)] # Mainland SEA
AOEs1deg[[6]] <- AOEs[c(28, 9, 14, 29)] # East Asia
names(AOEs1deg) <- c("Ind", "Phi", "And", "Sun", "MSE", "Eas")

setwd(wd)

#### half degree
AOEhalfdeg.dir <- ("NDM/IM60_halfdeg_AOE_2pct")
setwd(AOEhalfdeg.dir)

AOEs <- list.files(pattern = ".shp$"); AOEs
AOEshalfdeg <- list()

AOEshalfdeg[[1]] <- AOEs[c(23, 15)] # Indian subcontinent
AOEshalfdeg[[2]] <- AOEs[c(25, 14, 21, 16, 5, 9, 2, 10, 22, 31, 11, 12)] # Philippines
AOEshalfdeg[[3]] <- AOEs[c(27, 6)] # Andaman-Nicobar
AOEshalfdeg[[4]] <- AOEs[c(20, 1, 26, 24, 4, 32, 29, 17, 8, 19, 7)] # Sundaland
AOEshalfdeg[[5]] <- AOEs[c(18)] # Mainland SEA
AOEshalfdeg[[6]] <- AOEs[c(28, 33, 13, 3, 30)] # East Asia
names(AOEshalfdeg) <- c("Ind", "Phi", "And", "Sun", "MSE", "Eas")


setwd(wd)

#### 5 degrees 
AOE5deg.dir <- ("NDM/IM60_5deg_AOE_2pct")
setwd(AOE5deg.dir)

AOEs <- list.files(pattern = ".shp$"); AOEs
AOEs5deg <- list()

AOEs5deg[[1]] <- AOEs[3] # Indian subcontinent
AOEs5deg[[2]] <- AOEs[2] # Philippines
AOEs5deg[[3]] <- AOEs[10] # Andaman-Nicobar
AOEs5deg[[4]] <- AOEs[c(7, 8, 6)] # Sundaland
AOEs5deg[[5]] <- AOEs[4] # Mainland SEA
AOEs5deg[[6]] <- AOEs[c(1, 9, 5)] # East Asia
names(AOEs5deg) <- c("Ind", "Phi", "And", "Sun", "MSE", "Eas")


setwd(wd)


#### 2 degrees 
AOE2deg.dir <- ("NDM/IM60_2deg_AOE_2pct")
setwd(AOE2deg.dir)

AOEs <- list.files(pattern = ".shp$"); AOEs
AOEs2deg <- list()

AOEs2deg[[1]] <- AOEs[c(19, 7)] # Indian subcontinent
AOEs2deg[[2]] <- AOEs[c(15, 20, 12)] # Philippines
AOEs2deg[[3]] <- AOEs[14] # Andaman-Nicobar
AOEs2deg[[4]] <- AOEs[c(2, 9, 1, 4, 8, 3)] # Sundaland
AOEs2deg[[5]] <- AOEs[c(13, 6, 10, 5)] # Mainland SEA
AOEs2deg[[6]] <- AOEs[c(17, 16, 18, 11)] # East Asia
names(AOEs2deg) <- c("Ind", "Phi", "And", "Sun", "MSE", "Eas")

setwd(wd)

################ 05 AOE "by-resolution" figure ################
# Using 5° areas as base, plot the "byresol" figure 
names(AOEs5deg) <- c("IND", "PHI", "AND", "SUN", "MSEA", "EAS")
names(AOEs2deg) <- c("IND", "PHI", "AND", "SUN", "MSEA", "EAS")
names(AOEs1deg) <- c("IND", "PHI", "AND", "SUN", "MSEA", "EAS")
names(AOEshalfdeg) <- c("IND", "PHI", "AND", "SUN", "MSEA", "EAS")

pdf(paste0("NDM/NDM_byresol", "_v9", ".pdf"), height = 25, width = (25*3/2))

par(mfrow=c(4,6))
par(mar=c(1,2,5,1))
for(i in 1:length(AOEs5deg)){
  print(names(AOEs5deg)[i])
  group <- AOEs5deg[[i]]
  plot(asia_tropic_map, main=names(AOEs5deg)[i],
       cex.main=4
       #cex.main=2
       )
  for(j in 1:length(AOEs5deg[[i]])){
    AOE <- shapefile(paste0(AOE5deg.dir, "/", group[j]))
    crs(AOE) <- wgs84
    plot(AOE, col=paste0(map.col.pal[i],"80"), border=NA, add=T)
  }
}

for(i in 1:length(AOEs2deg)){
  print(names(AOEs2deg)[i])
  group <- AOEs2deg[[i]]
  plot(asia_tropic_map, main=names(AOEs5deg)[i],
       cex.main=4
       #cex.main=2
  )
  for(j in 1:length(AOEs2deg[[i]])){
    AOE <- shapefile(paste0(AOE2deg.dir, "/", group[j]))
    crs(AOE) <- wgs84
    plot(AOE, col=paste0(map.col.pal[i],"80"), border=NA, add=T)
  }
}

for(i in 1:length(AOEs1deg)){
  print(names(AOEs1deg)[i])
  group <- AOEs1deg[[i]]
  plot(asia_tropic_map, main=names(AOEs5deg)[i],
       cex.main=4
       #cex.main=2
  )
  for(j in 1:length(AOEs1deg[[i]])){
    AOE <- shapefile(paste0(AOE1deg.dir, "/", group[j]))
    crs(AOE) <- wgs84
    plot(AOE, col=paste0(map.col.pal[i],"80"), border=NA, add=T)
  }
}

for(i in 1:length(AOEshalfdeg)){
  print(names(AOEshalfdeg)[i])
  group <- AOEshalfdeg[[i]]
  plot(asia_tropic_map, main=names(AOEs5deg)[i],
       cex.main=4
       #cex.main=2
  )
  for(j in 1:length(AOEshalfdeg[[i]])){
    AOE <- shapefile(paste0(AOEhalfdeg.dir, "/", group[j]))
    crs(AOE) <- wgs84
    plot(AOE, col=paste0(map.col.pal[i],"80"), border=NA, add=T)
  }
}
dev.off()

################ 06 Figures for each resolution's consensus areas (or AOEs) ################
# Based on the results' .out files, we manually recorded the area names and endemicity scores of all AOEs across all four analytical resolutions.
# We saved these numbers in a spreadsheet named "NDMconsensus.xlsx".

# Set up data frames with each resolution's consensus areas, their scores, and shp filenames
consensus_5deg <- read.xlsx2("NDMconsensus.xlsx", 1)[c(1:10),]
consensus_2deg <- read.xlsx2("NDMconsensus.xlsx", 2)[c(1:20),]
consensus_1deg <- read.xlsx2("NDMconsensus.xlsx", 3)[c(1:31),]
consensus_halfdeg <- read.xlsx2("NDMconsensus.xlsx", 4)[c(1:33),]

consensus_5deg$files <- list.files(AOE5deg.dir,".shp")
consensus_2deg$files <- list.files(AOE2deg.dir,".shp")
consensus_1deg$files <- list.files(AOE1deg.dir,".shp")
consensus_halfdeg$files <- list.files(AOEhalfdeg.dir,".shp")


pdf(paste0("NDM/NDM_AOEs_5deg", "_v2", ".pdf"), height = 25, width = (25*2))
par(mar=c(1,1,10,1))
par(mfrow=c(2,5))

for(i in 1:length(unlist(AOEs5deg))){
  plot(asia_tropic_map, 
       main = paste0("Area ",i,": \n",consensus_5deg$Name[which(consensus_5deg$files == unlist(AOEs5deg)[i])]), 
       cex.main=5
       )
  print(unlist(AOEs5deg)[i])
  AOE <- shapefile(paste0(AOE5deg.dir, "/", unlist(AOEs5deg)[i]))
  crs(AOE) <- wgs84
  plot(AOE, col=paste0(map.col.pal[4],"60"), border = "grey50", add=T)
  text(73,c(-2,-5), cex=3, adj=0,
         labels=c(paste0("Min E score: ", round(as.numeric(consensus_5deg$Min[which(consensus_5deg$files == unlist(AOEs5deg)[i])]), 2)),
                  paste0("Max E score: ", round(as.numeric(consensus_5deg$Max[which(consensus_5deg$files == unlist(AOEs5deg)[i])]), 2))))
  
}

dev.off()


pdf(paste0("NDM/NDM_AOEs_2deg", "_v2", ".pdf"), height = 25, width = (25*3/2))
par(mar=c(1,1,7,1))
par(mfrow=c(4,5))

for(i in 1:length(unlist(AOEs2deg))){
  plot(asia_tropic_map, 
       main = paste0("Area ",i,": \n",consensus_2deg$Name[which(consensus_2deg$files == unlist(AOEs2deg)[i])]), 
       cex.main=3.7
  )
  print(unlist(AOEs2deg)[i])
  AOE <- shapefile(paste0(AOE2deg.dir, "/", unlist(AOEs2deg)[i]))
  crs(AOE) <- wgs84
  plot(AOE, col=paste0(map.col.pal[4],"60"), border = "grey50", add=T)
  text(73,c(-2,-5), cex=2.4, adj=0,
       labels=c(paste0("Min E score: ", round(as.numeric(consensus_2deg$Min[which(consensus_2deg$files == unlist(AOEs2deg)[i])]), 2)),
                paste0("Max E score: ", round(as.numeric(consensus_2deg$Max[which(consensus_2deg$files == unlist(AOEs2deg)[i])]), 2))))
  
}

dev.off()


pdf(paste0("NDM/NDM_AOEs_1deg", "_v3", ".pdf"), height = (10), width = (10*2))
par(mar=c(1,1,14,1))
par(mfrow=c(2,4))

for(i in 1:length(unlist(AOEs1deg))){
  plot(asia_tropic_map, 
       main = paste0("Area ",i,": \n",consensus_1deg$Name[which(consensus_1deg$files == unlist(AOEs1deg)[i])]), 
       cex.main=2
  )
  print(unlist(AOEs1deg)[i])
  AOE <- shapefile(paste0(AOE1deg.dir, "/", unlist(AOEs1deg)[i]))
  crs(AOE) <- wgs84
  plot(AOE, col=paste0(map.col.pal[4],"80"), border = "grey50", add=T)
  text(73,c(-2,-5), cex=1.4, adj=0,
       labels=c(paste0("Min E score: ", round(as.numeric(consensus_1deg$Min[which(consensus_1deg$files == unlist(AOEs1deg)[i])]), 2)),
                paste0("Max E score: ", round(as.numeric(consensus_1deg$Max[which(consensus_1deg$files == unlist(AOEs1deg)[i])]), 2))))
  
}

dev.off()


pdf(paste0("NDM/NDM_AOEs_halfdeg", "_v3", ".pdf"), height = (8), width = (8*3/2))
par(mar=c(1,1,10,1))
par(mfrow=c(2,3))

for(i in 1:length(unlist(AOEshalfdeg))){
  plot(asia_tropic_map, 
       main = paste0("Area ",i,": \n",consensus_halfdeg$Name[which(consensus_halfdeg$files == unlist(AOEshalfdeg)[i])]), 
       cex.main=2
  )
  print(unlist(AOEshalfdeg)[i])
  AOE <- shapefile(paste0(AOEhalfdeg.dir, "/", unlist(AOEshalfdeg)[i]))
  crs(AOE) <- wgs84
  plot(AOE, col=paste0(map.col.pal[4],"80"), border = "grey50", add=T)
  text(73,c(-2,-5), cex=1.4, adj=0,
       labels=c(paste0("Min E score: ", round(as.numeric(consensus_halfdeg$Min[which(consensus_halfdeg$files == unlist(AOEshalfdeg)[i])]), 2)),
                paste0("Max E score: ", round(as.numeric(consensus_halfdeg$Max[which(consensus_halfdeg$files == unlist(AOEshalfdeg)[i])]), 2))))
  
}

dev.off()
