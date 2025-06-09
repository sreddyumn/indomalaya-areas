############ R CODE PART 4: ############
############ Hierarchical clustering analysis and plotting ############
### NOTE: This script may involve R objects from "Part 1", "Part 2", and "Part 3"

# setwd("your working directory")
wd <- getwd() # main working directory that we can easily return to

packages <- c("sp","sf","dismo","raster","rgeos","rvertnet","rgdal","rworldmap","rworldxtra",
              "letsR","xlsx","dplyr","ecostructure","epm","lwgeom","tmap","terra","doParallel",
              "tidyverse","ggnetwork","vegan","ape","cluster","RColorBrewer","dendextend","NbClust",
              "factoextra","fpc","clValid","pvclust","colorspace","gclus","mapproj","viridis","cowplot")
lapply(packages, require, character.only=T)


################## 01 Hierarchical Clustering algorithm selection ##################
# References: 
# https://doi.org/10.1007/978-3-319-71404-2_4
# https://doi.org/10.1111/ddi.13535

# Load letsR presence-absence matrix 
matrix_1 <- as.data.frame(PAM_60_1$P) # This object was created in "Part 1" using letsR
View(matrix_1)

# Sort data first by descending latitude, and then by ascending longitude
library(dplyr)
matrix_1 <- arrange(matrix_1, desc(`Latitude(y)`), `Longitude(x)`)
View(matrix_1)

# Create a character vector with longitudes followed by latitudes connected with an "_"
longitudes <- matrix_1$`Longitude(x)` # Extract vector of longitude values
latitudes <- matrix_1$`Latitude(y)` # Extact vector of latitude values
species <- colnames(matrix_1[,-c(1:2)]) # Extract vector of species names
longlat <- paste(longitudes,latitudes,sep="_")
rownames(matrix_1) <- longlat # set rownames as long-lat values
matrix_1 <- matrix_1[,-c(1:2)] # remove original lat and long columns after assigning rownames

# Now REMOVE cells with less than 5 species, to simplify the dataset
matrix_2 <- matrix_1[which(apply(matrix_1, 1, sum) >= 5),]
# View(matrix_2)

# Calculate beta-sim index (which is independent of the effects of species richness)
betasim <- betadiver(matrix_2, method = 22); class(betasim)
# betadiver(help=T)
# "sim" = pmin(b,c)/(pmin(b,c)+a)

?hclust
# 8 hclust methods to assess (actually 7, but there are two versions of Ward's method
hclust.methods <- c("ward.D", "ward.D2", "single", "complete", "average", 
                    "mcquitty", # (= WPGMA), 
                    "median", #(= WPGMC)
                    "centroid") # (=UPGMC)

# Run hclust on betasim matrix for each of the 8 methods
for(i in 1:length(hclust.methods)){
  hc <- hclust(betasim, method=hclust.methods[i])
  assign(paste0("betasim.",hclust.methods[i]), hc)
}

# Assess cophenetic correlation coefficient to determine the best method to proceed with
for(i in 1:length(hclust.methods)){
  hc <- get(paste0("betasim.",hclust.methods[i]))
  coph <- cophenetic(hc) # cophenetic distance matrix
  print(paste0(hclust.methods[i], ": ", cor(betasim, coph)))
}

# "ward.D :  0.699971416410193"
# "ward.D2 :  0.725530452504486"
# "single :  0.114159298842664"
# "complete :  0.596504501358013"
# "average :  0.759539736415385"
# "mcquitty :  0.671040274577702"
# "median :  0.61828037212314"
# "centroid :  0.746679746666982"

# UPGMA is the best one

# Rename object for better clarity 
betasim.UPGMA <- betasim.average
plot(betasim.UPGMA)

################## 02a Assess optimal clusters: Fusion levels ##################
## Assessing fusion levels (p. 75 in Borcard et al. 2018) 
## "When read from right to left, long horizontal lines preceding steep increases suggest cutting levels"
dev.new()
plot(
  #betasim.UPGMA$height,
  #nrow(betasim):2,
  tail(betasim.UPGMA$height, 49),
  50:2,
  yaxt = "n",
  type = "S",
  main = "Fusion levels - UPGMA",
  ylab = "k (number of clusters)",
  xlab = "h (node height)",
  col = "grey"
)
axis(2, at = seq(0, 50, by = 5), las=2)
text(#betasim.UPGMA$height,
     tail(betasim.UPGMA$height, 49),
     #nrow(betasim):2,
     #nrow(betasim):2,
     50:2,
     50:2,
     col = "red",
     cex = 0.8)
# abline(h=2, col = "red)
### Note: Starting from 50 clusters seem biologically more realistic than 1000+ 
### By this method, it looks like 9, 7, and 5 are possible choices to assess

nth(betasim.UPGMA$height, nrow(betasim)-7)
tail(betasim.UPGMA$height, 7) # height at k=7 is 0.4781859

?pvclust

################## 02b Assess optimal clusters: Silhouette Index ##################
# "[The] larger the value is, the better the object is clustered. Negative values
# suggest that the corresponding objects may have been placed in the wrong cluster."

# Choose and rename the dendrogram ("hclust" object)
hc <- betasim.average
# Plot average silhouette widths for all partitions 
# except for the trivial partitions (k = 1 or k = n)
Si <- numeric(nrow(betasim)) # Note: vector of zeros
for (k in 2:(nrow(betasim) - 1)) {
  print(k)
  sil <- silhouette(cutree(hc, k = k), betasim)
  Si[k] <- summary(sil)$avg.width
}


k.best.si <- which.max(Si)
sort(Si, partial=(length(Si)))[length(Si)]
which(Si == sort(Si, partial=(length(Si)-1))[length(Si)-1]) # second best k
which(Si == sort(Si, partial=(length(Si)-1))[length(Si)-2]) # third best k
# k=6 is the 'most optimal' cluster by Silhouette index

dev.new()
plot(
  # 1:nrow(betasim),
  c(1:100),
  Si[c(1:100)],
  type = "h",
  main = "Average Silhouette Width (ASW) by number of clusters",
  xlab = "k (number of clusters)",
  ylab = "Average silhouette width"
)
axis(
  1,
  k.best.si,
  paste("optimum", k.best.si, sep = "\n"),
  col = "red",
  font = 2,
  col.axis = "red"
)
points(k.best.si,
       max(Si),
       pch = 16,
       col = "red",
       cex = 1)

# Make a sortable Average Silhouette Width table
Si.df <- data.frame(c(1:length(Si)),Si)
colnames(Si.df) <- c("k", "ASW")

head(Si.df[order(Si.df$ASW, decreasing=T),], 10)

# Top 3 partitions: k = 6, 7, and 5
dev.off()

################## 02c Assess optimal clusters: Matrix correlation ##################
# Comparison Between the Dissimilarity Matrix and Binary Matrices Representing Partitions
# Optimal number of clusters according to matrix correlation statistic (Pearson)

# Function to compute a binary dissimilarity matrix from clusters
grpdist <- function(X) {
  require(cluster)
  gr <- as.data.frame(as.factor(X))
  distgr <- daisy(gr, "gower")
  distgr
}

kt <- data.frame(k = 1:nrow(betasim), r = 0)
for (i in 2:(nrow(betasim) - 1)) {
  print(i)
  gr <- cutree(hc, i)
  distgr <- grpdist(gr)
  mt <- cor(betasim, distgr, method = "pearson")
  kt[i, 2] <- mt
}
k.best.gr <- which.max(kt$r)
sort(kt$r, partial=(length(kt$r)))[length(kt$r)] # highest corr value (for k=9)
# k=9 is the 'most optimal' cluster by matrix correlation
which(kt$r > 0.71) # k= 7 to 15 are all above 0.7; there is very little power to differentiate among k values

which(kt$r == sort(kt$r, partial=(length(kt$r)-1))[length(kt$r)-1])
sort(kt$r, partial=(length(kt$r)-1))[length(kt$r)-1] # 2nd highest corr (for k=10)
which(kt$r == sort(kt$r, partial=(length(kt$r)-2))[length(kt$r)-2])
sort(kt$r, partial=(length(kt$r)-2))[length(kt$r)-2] # 3rd highest corr (for k=8)
which(kt$r == sort(kt$r, partial=(length(kt$r)-3))[length(kt$r)-3])
sort(kt$r, partial=(length(kt$r)-3))[length(kt$r)-3] # 4th highest corr (for k=7)
which(kt$r == sort(kt$r, partial=(length(kt$r)-4))[length(kt$r)-4])
sort(kt$r, partial=(length(kt$r)-4))[length(kt$r)-4] # 5th highest corr (for k=11)

plot(
  kt$k[1:100],
  kt$r[1:100],
  type = "h",
  main = "Matrix correlation-optimal number of clusters",
  xlab = "k (number of clusters)",
  ylab = "Pearson's correlation"
)
axis(
  1,
  k.best.gr,
  paste("optimum", k.best.gr, sep = "\n"),
  col = "red",
  font = 2,
  col.axis = "red"
)
points(k.best.gr,
       max(kt$r),
       pch = 16,
       col = "red",
       cex = 1)



?cutree
betasim.UPGMA.k9 <- cutree(betasim.UPGMA, k=9)
betasim.UPGMA.k5 <- cutree(betasim.UPGMA, k=5)
betasim.UPGMA.k6 <- cutree(betasim.UPGMA, k=6)
betasim.UPGMA.k7 <- cutree(betasim.UPGMA, k=7)

class(betasim.UPGMA.k7) # cutree creates a vector of group memberships (as integers) for all cells


################## 03 Cluster Membership Visualization (k=7) ##################
# Reorder the betasim UPGMA hclust object to be as close as possible to the original betasim matrix arrangment
betasim.UPGMA.r <- reorder.hclust(betasim.UPGMA, betasim)
# View(as.matrix(betasim))

# Create a dendrogram from the hclust object
UPGMA.dend <- as.dendrogram(betasim.UPGMA.r)

# Color the dendrogram branches by the following colors
map.col.pal <- c("#88CCEE", "#DDCC77", "#CC6677", "#117733", "#AA4499", "#332288", "#888888")

UPGMA.dend <- UPGMA.dend %>% color_branches(k = 7, col=map.col.pal[c(1:7)]) # %>%
              # set("leaves_cex" ) # %>%
              # set("branches_lwd", c(2,1,2)) # %>%
              # set("branches_lty", c(1,2,1))
plot(UPGMA.dend)

# Create an attributes table to track cells (leaves), colors, and cluster id
dendro.attributes <- data.frame(get_leaves_attr(UPGMA.dend, attribute = "label"), get_leaves_branches_attr(UPGMA.dend))
colnames(dendro.attributes) <- c("cell", "color")
for(i in 1:nrow(dendro.attributes)){
  dendro.attributes$long[i] <- as.numeric(substring(dendro.attributes$cell[i], 1,regexpr("_",dendro.attributes$cell[i])-1))
  dendro.attributes$lat[i] <- as.numeric(substring(dendro.attributes$cell[i], regexpr("_",dendro.attributes$cell[i])+1))
  dendro.attributes$cluster[i] <- which(map.col.pal == dendro.attributes$color[i])
}
dendro.attributes <- arrange(dendro.attributes, desc(lat), long) # rearrange this table by lat, long
View(dendro.attributes)

# now the coordinates and the attribute table's cell values are aligned
identical(rownames(longlat_matrix2@coords), dendro.attributes$cell)


# designate branches to draw k=7 rectangles
rd7 <- rect.dendrogram(UPGMA.dend, k=7, lty = 5, lwd = 0, 
                       lower_rect = -0.08)
# Set up the label text
rect.text.7 <- c("IND", "HIM", "MLD", "TWN", "SUN", "PHI", "AND")

# Thicken the branches
UPGMA.dend.thick <- UPGMA.dend %>% set("branches_lwd", 2)

##### PLOT COLORED DENDROGRAM AT K=7 #####
pdf("HClust/hclust-k7_v6.pdf", width = 20, height = 10)
par(mfrow=c(1,2))

# Plot the dendrogram
# UPGMA.dend the original k=7-colored dendrogram
plot(UPGMA.dend.thick, leaflab =  "none", ylab="beta-sim",
     #edgePar=list(lwd=2, p.lwd=2)
     # nodePar = list(lab.cex = 0.6, pch = c(NA, 19), 
     #                cex = 0.7, col = "blue")
     )
?plot.dendrogram

tail(betasim.UPGMA$height, 9) # First 9 partition heights
betasim.UPGMA$height[length(betasim.UPGMA$height)] # k=2
# minus n from the subsetting variable to get k=2+n
# so to plot a line marking k=7, I need length(betasim.UPGMA$height)-5

abline(h=betasim.UPGMA$height[length(betasim.UPGMA$height)-5], # k=7
       col="red", lty=2)
text(x=-20, y=0.01+betasim.UPGMA$height[length(betasim.UPGMA$height)-5],
     labels="k=7", cex=0.8)

abline(h=betasim.UPGMA$height[length(betasim.UPGMA$height)-4], # k=6
       col="red", lty=2)
text(x=-20, y=0.01+betasim.UPGMA$height[length(betasim.UPGMA$height)-4],
     labels="k=6", cex=0.8)

abline(h=betasim.UPGMA$height[length(betasim.UPGMA$height)-3], # k=5
       col="red", lty=2)
text(x=-20, y=0.01+betasim.UPGMA$height[length(betasim.UPGMA$height)-3],
     labels="k=5", cex=0.8)

abline(h=betasim.UPGMA$height[length(betasim.UPGMA$height)-7], # k=9
       col="red", lty=2)
text(x=-20, y=0.01+betasim.UPGMA$height[length(betasim.UPGMA$height)-7],
     labels="k=9", cex=0.8)


# Draw rectangles around each cluster
for(i in 1:length(rd7)){
  print(i)
  rect.dendrogram(UPGMA.dend, k=7, lty = 5, lwd = 0, which = i,
                  col = paste0(map.col.pal[i], 
                              # "20"),
                               "30"),
                  lower_rect = -0.02,
                  upper_rect = 0.01,
                  border=0)
}

# Print text for clusters
text(x=(head(cumsum(c(1, lengths(rd7))), -1) + ((lengths(rd7)-1)/2)), 
     y=0.5, font=2, cex = 1.2,
     col=map.col.pal[1:7], 
     labels=rect.text.7)

dev.off()

##### PLOT CLUSTERS AT K=7 ON A MAP #####
## Assign group ID to each lat-long coordinate for visualizing group membership on map
## Create a dataframe (long, lat, and spp richness columns) of all grid cells...
longlat_matrix2 <- data.frame(long = as.numeric(substring(rownames(matrix_2), 1,regexpr("_",rownames(matrix_2))-1)),
                              lat = as.numeric(substring(rownames(matrix_2), regexpr("_",rownames(matrix_2))+1)),
                              richness = apply(matrix_2, 1, sum))
rownames(longlat_matrix2) <- rownames(matrix_2) # these are individual cells 
# ...and create a spatial object (points) from the dataframe
coordinates(longlat_matrix2) <- ~ long + lat
proj4string(longlat_matrix2) <- wgs84 # assign CRS

pdf("HClust/hclust-k7-MAP_v2.pdf", width = 8, height = 8)
plot(longlat_matrix2, pch=20, 
     col = map.col.pal[as.factor(dendro.attributes$cluster)],
     cex=2); plot(asia_tropic_map, add=T)
dev.off()



################## 04 Standardize plotting other k values, up to 36 ##################
col36 <- c("#88CCEE","#DDCC77","#CC6677","#117733","#AA4499","#332288","#888888",
           "#999933","#661100","#882255","#44AA99","dodgerblue2","#E31A1C","maroon",
           "#6A3D9A","#FF7F00", "brown","gold1","skyblue2","#FB9A99","palegreen2",
           "#CAB2D6","#FDBF6F","khaki2","green1","gray70","orchid1","deeppink1",
           "steelblue4","blue1","darkturquoise","green4","yellow4","yellow3",
           "darkorange4","gray12")

hclustmap <- function(X, k, colset=col36,labels_cex=0.5, matrix=longlat_matrix2,
                      plot.pch=15, basemap=asia_tropic_map) {
  # X should be betasim.UPGMA.r
  cut_tree <- as.dendrogram(X) %>% color_branches(k = k, col=colset[c(1:k)]) %>%
    set("labels_cex", 0.5)
  # assign(paste0("UPGMA.dend.k",k), cut_tree)
  attrib <- data.frame(get_leaves_attr(cut_tree, attribute = "label"), get_leaves_branches_attr(cut_tree))
  colnames(attrib) <- c("cell", "color")
  for(i in 1:nrow(attrib)){
    attrib$long[i] <- as.numeric(substring(attrib$cell[i], 1,regexpr("_",attrib$cell[i])-1))
    attrib$lat[i] <- as.numeric(substring(attrib$cell[i], regexpr("_",attrib$cell[i])+1))
    attrib$cluster[i] <- which(colset == attrib$color[i])
  }
  attrib <- arrange(attrib, desc(lat), long)
  plot(matrix, pch=plot.pch, main = paste0("UPGMA Clusters, k = ", k),
       col = colset[as.factor(attrib$cluster)],
       cex=1.5); plot(basemap, add=T)
  
}
hclustmap(betasim.UPGMA.r, k=7)

pdf("HClust/hclust-k8-36_MAPS_v1.pdf", width = 12, height = 8)
# pdf("HClust/hclust-k8-36_MAPS_v2.pdf", width = 12/1.5, height = 8/1.5)
for(i in 8:36){
  print(i)
  hclustmap(betasim.UPGMA.r, k=i)
}
dev.off()


