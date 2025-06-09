############ R CODE PART 3: ############
############ ECOSTRUCTURE analyses and plotting ############
### NOTE: This script may involve R objects from "Part 1" and "Part 2"

# setwd("your working directory")
wd <- getwd() # main working directory that we can easily return to

packages <- c("sp","sf","dismo","raster","rgeos","rvertnet","rgdal","rworldmap","rworldxtra",
              "letsR","xlsx","dplyr","ecostructure","epm","lwgeom","tmap","terra","doParallel",
              "tidyverse","ggnetwork","vegan","ape","cluster","RColorBrewer","dendextend","NbClust",
              "factoextra","fpc","clValid","pvclust","colorspace","gclus","mapproj","viridis","cowplot")
lapply(packages, require, character.only=T)


################ 01 Ecostructure installation ################
### installing dependencies ###
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

BiocManager::install(version = "3.17")
library(BiocManager)

BiocManager::install("Biobase")
install.packages("devtools")
library(devtools)
devtools::install_github("kkdey/methClust");  
devtools::install_github("kkdey/CountClust")
devtools::install_github("TaddyLab/maptpx")

### Loading dependencies ###
library(Biobase)
library(methClust)
library(CountClust)
library(maptpx)

### installing and loading ecostructure ###
devtools::install_github("kkdey/ecostructure")
library(ecostructure)

# help(package = "ecostructure") # for comprehensive list of functions

################ 02 Ecostructure modal fitting ################
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

k_max <- 20
Kvec <- c(2:k_max) # K values: prior settings of motif number

# function to run the model for parallel computing
ecos_fit_K <- function(n){
  # K = n, the number of the motif to model
  pres_ab_fit <- ecos_fit(matrix_1, K = n, tol = 0.1, num_trials = 100)
  saveRDS(pres_ab_fit, file=paste0("Ecostructure/IM60_gom_fit_K",n,".rds"))
}

print("start model fitting:")
# initiate cluster

?makeCluster
#cl <- makeCluster(length(Kvec), type="PSOCK")
cl <- makeCluster(detectCores(), type="PSOCK")
registerDoParallel(cl)
clusterEvalQ(cl, {
  library(ecostructure)
})
clusterExport(cl, "matrix_1")

parLapply(cl, Kvec, ecos_fit_K)
stopCluster(cl) # shutdown cluster
print("DONE")


################## 03 Extract motif contribution rasters from the results of model fitting ################## 
# References: https://rdrr.io/github/kkdey/ecostructure/src/R/ecos_fit.R 

# This function allows the creation of a raster brick from the output of each K run
ecos.rasterize <- function(x, res = c(1,1), crs="+proj=longlat +datum=WGS84 +no_defs"){
  runfile <- readRDS(paste0("Ecostructure/IM60_gom_fit_K",x,".rds"))
  omega <- as.data.frame(runfile$omega)
  omega$long <- substring(rownames(omega),1,regexpr("_",rownames(omega))-1)
  # substring the long_lat string, starting at character 1 until 1 before "_"
  omega$lat <- substring(rownames(omega),regexpr("_",rownames(omega))+1)
  # substring the long_lat string, always starting 1 character after "_"
  omega <- omega[,c(ncol(omega)-1,ncol(omega),1:(ncol(omega)-2))]
  # re-order the data frame with the long and lat columns coming first and second respectively
  omega.rast <- rasterFromXYZ(omega, res=res,crs=crs)
  return(omega.rast)
}

# Create a raster brick for each K run
for (i in 2:k_max) {
  assign(paste0("ecosomega", i, ".rast"), ecos.rasterize(i)) 
  # assign the "value" of the ecos.rasterize() function to the name stated within paste0()
}

# Set up an empty raster for raster::extract 
#ncol(ecosomega4.rast$X2) (checking a real output for reference)
#nrow(ecosomega4.rast$X2)
rasterbase <- raster(ncol=49, nrow=67, xmn=65, xmx=132, ymn=-12, ymx=37, res=1)
rasterbase[]<-0 # set all values to 0
rasterbase
# Set up an empty spatial polygon with the same extent
rasterbound <- as(extent(65, 132, -12, 37), "SpatialPolygons")
crs(rasterbound) <- wgs84

# For coloring the raster when plotting
pie_darkblue <- "#332288"

# Big for loop to generate a list (for each K) of motif dataframes
for(K in 2:k_max){ # loop from K=2 to K=20
  # for(K in 2:3){
  print(paste0("run ", K)) # tracker
  omegabricks <- ecos.rasterize(K) # raster brick for the given K run (a list of K raster layers)
  hold <- list() # generate an empty list to hold all K dfs (one for each motif)
  for(i in 1:length(omegabricks@data@names)){ # loop over all K motifs within the Kth run
    print(paste0("motif ", i)) # tracker
    hold[[i]] <- raster::extract(x=omegabricks[[i]], # extract values for ith motif using the motif raster
                                  y=rasterbound, # using the extent polygon as reference
                                  df=T, cellnumbers=T) %>%  # output as a dataframe and include cell numbers 
      arrange(cell) %>% rename(m_contrib = 3) # arrange the rows by cell number, and then rename3rd column as "m_contrib"
    rastercoords <- coordinates(omegabricks[[i]])[hold[[i]][,2],] # extract lat long values
    hold[[i]] <- cbind(hold[[i]], rastercoords) 
    for(j in 1:nrow(hold[[i]])){ # loop through all cells within a motif df
      formatted <- as.character(format(round(hold[[i]][,3], 2), nsmall = 2)) # round motif contribution (omega) to 2 decimal points
      if(is.na(hold[[i]][,3][j]) == TRUE){ # create a color column, then assign a hex code with transparency based on m_contrib
        # hold[[i]]$color[j] <- paste0(pie_darkblue,"00")
        hold[[i]]$color[j] <- paste0(pie_darkblue,"00")
      } else if(formatted[j] == "1.00") {
        # hold[[i]]$color[j] <- paste0(pie_darkblue,"99")
        hold[[i]]$color[j] <- paste0(pie_darkblue,"99")
      } else {
        # hold[[i]]$color[j] <- paste0(pie_darkblue,substring(formatted[j], 3,5))
        hold[[i]]$color[j] <- paste0(pie_darkblue,substring(formatted[j], 3,5))
      }
    }
  }
  assign(paste0("ecosomega.K",K), hold) # each assigned object will be a list (of lists) for each Kth run 
                                        # (containing K number of dfs (one for each motif))
}



####### 04 Calculate motif overlap (Jaccard index) #######
# Code by Li et al. (2021) New Phytologist, https://doi.org/10.1111/nph.17443
# from https://github.com/reelab/hengduan-ecostructure/blob/main/Rscript_1_floristic_structure_analysis.R
# and https://github.com/reelab/hengduan-ecostructure/blob/main/Rscript_0_functions.R

fitK <- lapply(paste0("Ecostructure/IM60_gom_fit_K",c(2:20),".rds"), readRDS)
# extract motif contribution (omega) matrix into a list and rename column names with "m_" prefix
motif_list <- lapply(fitK, function(x) x$omega %>% as.data.frame() %>% 
                       rename_all(list(~paste0("m_", .))))

motif_bi_list <- list() # a binary(1,0) list to hold dfs (for each K run) where we can assign a site with the dominant motif
for(i in 1:length(motif_list)){ # loop through all 19 omega matrices
  print(i) # tracker 
  df <- motif_list[[i]] 
  df_bi <- df # create a second df of the same dimension as the omega matrix
  df_bi[] <- 0 # make all values 0
  
  # assigning a site with the dominant motif
  df$group <- sapply(1:nrow(df), function(x) names(df)[which.max(df %>% slice(x))])
  # create a group column
  # "loop" (using sapply) through each row (using function slice() on the row number)
  # extract the name of the column where the omega value is the largest for that given cell
  # and output as a value in the df$group column for that given site
  
  for(g_id in 1:nrow(df)){ # loop through all rows within the given omega matrix
    row2gp <- df$group[g_id] # take the group id (motif name) of that row (site) 
    df_bi[g_id, which(names(df_bi) == row2gp)] <- 1 # and give the corresponding site in the binary df a value of 1
  }
  motif_bi_list[[i]] <- df_bi # assign this binary df to the ith element in the bigger binary list
}

# By having a list of the dominant motifs for each given site (map grid cell), 
# we can trace motif identity across K values, 
# and assess the overlap between all motifs in K and all motifs in K+1.
# From K to K+1, the maximum pairwise overlap for each motif at K, 
# when averaged across all motifs, can be used to assess motif stability (?)

# function to calculate overlap between motifs
cal.overlap <- function(motif_bi_list, K_no=c(2:3)){
  
  df1_K <- K_no[1] # K (just the number)
  df2_K <- K_no[2] # K+1 (just the number)
  
  df1 <- motif_bi_list[[c(df1_K-1)]] # binary matrix for K = K
  df2 <- motif_bi_list[[c(df2_K-1)]] # binary matrix for K = K+1
  
  m_comb <- expand.grid(df1_col=1:df1_K,df2_col=1:df2_K) # create a df from all combinations of the supplied vectors
  m_comb$df1_motif <- paste0("K",df1_K,"_m",m_comb$df1_col) # Using the above df, create two more columns...
  m_comb$df2_motif <- paste0("K",df2_K,"_m",m_comb$df2_col) # ...with the K and m values 
  
  for(j in 1:nrow(m_comb)){
    df_temp <- apply(cbind(df1[,m_comb$df1_col[j]], df2[,m_comb$df2_col[j]]), 1, sum)
    # TEST CODE: df_temp <- apply(cbind(motif_bi_list[[1]][,(expand.grid(df1_col=1:2, df2_col=1:3))$df1_col[1]], motif_bi_list[[2]][,(expand.grid(df1_col=1:2, df2_col=1:3))$df2_col[1]]),1,sum)
    # sum across all combinations of the motifs between K & K+1 to get 3 possible values: 
    # 0: this cell doesn't belong to either motif (no intersection)
    # 1: this cell belongs to one motif but not the other 
    # 2: this cell belongs to both motifs (intersection)
    
    df_count <- as.data.frame(table(df_temp))
    no_AB <- ifelse("1" %in% df_count$df_temp,df_count$Freq[df_count$df_temp =="1"],0) 
    # looks complicated but it's just the number of cells with "1"
    # which also translates to "number of elements one set but not the other"
    no_C <- ifelse("2" %in% df_count$df_temp,df_count$Freq[df_count$df_temp =="2"],0) 
    # number of cells with "2", or "number of elements in the intersection"
    no_D <- ifelse("0" %in% df_count$df_temp,df_count$Freq[df_count$df_temp =="0"],0)
    # number of cells with "0", or "number of elements in neither set"
    
    m_comb$overlap_1[j] <- (no_C) / (no_AB + no_C) # Jaccard index (used in main analysis)
    m_comb$overlap_2[j] <- (no_D + no_C) / (no_D + no_AB + no_C)
    m_comb$overlap_3[j] <- (sqrt(no_C * no_D) + no_C) / (sqrt(no_C * no_D) + no_AB + no_C)
    
  }
  return(m_comb)
}

# overlap between motifs across K by Jaccard index 
K_comb <- data.frame(K1=2:(k_max-1), K2=3:k_max)
overlap_list <- list()
for(i in 1:nrow(K_comb)){
  print(i)
  overlap_list[[i]] <- cal.overlap(motif_bi_list=motif_bi_list, K_no=as.numeric(K_comb[i,]))
}
# saveRDS(overlap_list, file="Ecostructure/comb_motif_overlap.rds")

View(overlap_list[[2]])

# Export results as excel spreadsheet
write.xlsx(overlap_list[[1]], "Ecostructure/Motif_Overlap_List_v2.xlsx", sheetName="K2-K3", row.names = F)
for(i in 2:length(overlap_list)){
  print(paste0("K", i+1, "-K", i+2))
  write.xlsx(overlap_list[[i]], "Ecostructure/Motif_Overlap_List_v2.xlsx", 
             sheetName=paste0("K", i+1, "-K", i+2), row.names = F, append = T)  
}


####### 05 Track motif identity & calculate motif stability #######
K_comb # paired df of K & K+1 numbers created earlier

max_overlap_list <- list() # a list to hold all maximum overlap values for all motifs across all K runs

for(i in 1:length(overlap_list)){ # Go through overlap_list [[1]] to [[18]] (K2-K3 to K19-K20)
  print(paste0(K_comb[i,1], " to ", K_comb[i,2])) # tracker
  max_overlap_list[[i]] <- vector("numeric", (i+1)) # establish vector for the max overlap value 
                                                    # for all (i+1) motifs within the Kth run
  for(j in 1:(i+1)){ # Inside each overlap_list[[i]], there are (i+1) motifs (for j in 1:(i+1))
    # so pick out all the motif at K = (i+1) (which will have all i+1 values of the K_(i+1)_m_j)
    # and find the maximum overlap value for each (i+1)th motif
    # 
    if(any(overlap_list[[i]]$overlap_1[which(overlap_list[[i]]$df1_motif == paste0("K",(i+1),"_m",j))] == 'NaN') == FALSE) {
      max_overlap_list[[i]][j] <- max(overlap_list[[i]]$overlap_1[which(overlap_list[[i]]$df1_motif == paste0("K",(i+1),"_m",j))])
    }
    else{
      culprit <- which(is.nan(overlap_list[[i]]$overlap_1[which(
        overlap_list[[i]]$df1_motif == paste0("K",(i+1),"_m",j))])==TRUE)
      overlap_list[[i]]$overlap_1[which(overlap_list[[i]]$df1_motif == paste0("K",(i+1),"_m",j))][culprit] <- 0
      max_overlap_list[[i]][j] <- max(overlap_list[[i]]$overlap_1[which(overlap_list[[i]]$df1_motif == paste0("K",(i+1),"_m",j))])
    }
      
  }
  names(max_overlap_list[[i]]) <- paste0("motif",c(1:(i+1)))
}
names(max_overlap_list) <- paste0("K=",c(2:19)) # assign number of K run to each slot's name
max_overlap_list

avg_overlap_list <- lapply(max_overlap_list, mean) # calculate average overlap 
names(avg_overlap_list) <- str_replace(names(avg_overlap_list), "K=", "") # leave out "K=" in names
avg_overlap <- unlist(avg_overlap_list)
avg_overlap <- as.data.frame(avg_overlap) %>% tibble::rownames_to_column(var="K") # Turn into df, with K values as another variable

stdev_overlap_list <- lapply(max_overlap_list, sd) # calculate st dev
sdev <- unlist(stdev_overlap_list) # unlist into a named vector

# Plot motif stability 
pdf("Ecostructure/motifstability_v2.pdf")
plot(avg_overlap, pch = 18,  ylim=range(0.2,1.1),
     xaxt = "n", xlab="Clustering scheme increments", ylab="avg max overlap", 
     main = "Motif Stability")
arrows(x0=as.numeric(avg_overlap$K), y0=avg_overlap$avg_overlap-sdev, 
       x1=as.numeric(avg_overlap$K), y1=avg_overlap$avg_overlap+sdev, 
       code=3, angle=90, length = 0.04)
axis(1, at = seq(2,19, by = 1),
     las=2, cex.axis=0.6,
     labels = paste0(K_comb[,1], " to ", K_comb[,2]))
abline(h=0.9, col="red") # arbitrary red line to mark 0.9 average similarity
dev.off()



####### 06 Similarity network graph #######
overlap_list_F <- do.call(rbind, overlap_list) # create a new *filtered* "overlap list" squeezed into a df
View(overlap_list_F)
overlap_list_F <- overlap_list_F %>% # pass this through dplyr::filter (keep/subset rows based on column values)
                  filter(df1_col!= 1 & df2_col!= 1) %>% # keep all rows that aren't motif 1
                  filter(!grepl("K2",df1_motif)) # and filter out K2 since K2 is pointless without motif 1


overlap_list_F$overlap <- overlap_list_F$overlap_1 # we use Jaccard overlap index
overlap_thr <- 0.1


### NOTE: Manual node rearrangement happens here ###
### We manually assigned each motif from each K run a name (based as closely as the most salient geographical location or feature possible)
### and also assigned each one a "graphing index", which dictates the position their nodes will be in the network graph.
### This information is provided in the file "motif_names_v4.csv"

# Node position data for plotting
pos_dat <- rbind(overlap_list_F %>% select(df1_motif) %>% rename(node=df1_motif),
                overlap_list_F %>% select(df2_motif) %>% rename(node=df2_motif)) %>% unique()

# add names for nodes for each motif of each K
motif_names <- read.csv("Ecostructure/motif_names_v4.csv", fileEncoding="UTF-8-BOM")
motif_node_label <- motif_names %>% mutate(from = paste0(dt_set,"_",motif)) %>% 
  mutate(from = sub("motif_","m",from)) %>% select(from, motif_name)
#View(motif_names)
#View(motif_node_label)

# motifs as x coordinate for graph 
pos_dat$x <- as.numeric(sub("m","",word(pos_dat$node,2,2,sep=fixed("_")))) 
pos_dat$x <- as.numeric(sub("index_","",word(motif_names$graph_index,2,2,sep=fixed("_")))) 
# K level as y coordinate for graph (modified so that K goes downwards as it gets bigger)
pos_dat$y <- as.numeric(sub("K","",word(pos_dat$node,1,1,sep=fixed("_")))) * (-1) 
#View(pos_dat)

# a new set of x coordinates that makes the plot symmetric
pos_dat$x_2 <- sapply(1:nrow(pos_dat), function(k) pos_dat$x[k]+(12+pos_dat$y[k])/2)

# Segment data frame
seg_dat <- overlap_list_F %>% select(overlap, df1_motif, df2_motif) %>% # pick out only the relevant 3 columns from overlap_list_F
  rename(from = df1_motif, to = df2_motif) # and rename the 1st motif set to "from" and the 2nd to "to"
# Add coordinate data for graphing purpose
seg_dat <- seg_dat %>% left_join(pos_dat %>% rename(from = node)) %>% # designate starting xy coords
  left_join(pos_dat %>% rename(to = node, xend = x, yend = y, xend_2 = x_2)) # designate ending xy coords
#View(seg_dat)

# Prepare final data frame from pos_dat and seg_dat for plotting
m_net <- pos_dat %>% rename(from=node) %>% # take pos_dat rename "node" as "from"
  mutate(to=from, xend=x, yend=y, xend_2=x_2, overlap=NA) %>% # change other column names
  rbind(seg_dat %>% arrange(overlap)) # add seg_data and arrange by increasing overlap value
m_net$node_label <- sub("m","",word(m_net$from,2,2,sep=fixed("_"))) # add node labels (motif numbers)
#View(m_net)


m_net <- m_net %>% left_join(motif_node_label) # attach node names
# m_net$node_label <- ifelse(is.na(m_net$motif_name), m_net$node_label, "") # replace motif name with motif label if no motif name available
# m_net$motif_name <- sub(" ","\n",m_net$motif_name) # replace spaces with line breaks in motif_name column(???)

m_net_plot = ggplot(m_net, aes(x = x_2, y = y, xend = xend_2, yend = yend)) +
  geom_edges(aes(size=overlap, colour=overlap), alpha=0.8, 
             data = function(x){x[which(x$overlap >= overlap_thr),]}) +
  scale_colour_gradient(low="gray90", high="gray2",
                        name="overlap\nindex",
                        limits=c(overlap_thr,1.0),
                        breaks=c(0.1,0.4,0.7,1.0)) +
  # geom_nodes(size=rel(7), shape = 21, color = "gray10", fill = "white") +
  geom_nodes(size=rel(8), shape = 21, color = "gray10", fill = "white") +
  # geom_nodetext(aes(label=node_label), fontface="bold") + 
  geom_nodetext(aes(label=motif_name), size = rel(2), fontface="bold") + 
  annotate("text",x=-3.5,y=c(-3:-20), label=paste("K =",3:20),fontface = "bold") + 
  scale_size_continuous(guide = "none") +
  theme_blank() +
  guides(colour = guide_colourbar(title.position="top",ticks = FALSE, title.vjust=1))
plot(m_net_plot)
ggsave(filename="Ecostructure/motif_stability_IM60_v9.pdf",m_net_plot,width=13.5,height=9)
ggsave(filename="Ecostructure/motif_stability_IM60_v8.pdf",m_net_plot, width=18,height=9)


####### 07 Examine species contribution #######
# At each K, one can examine the n top contributing species to a given motif

IM60_gom_fit_K6 <- readRDS(paste0("Ecostructure/IM60_gom_fit_K",6,".rds"))
IM60_gom_fit_K6$theta

plot(ecosomega6.rast)

# use the function ExtractTopFeatures from the R package CountClust 
# to compute and record the top contributing species to each motif
featurescount <- 20 # change the number of species to extract
K6features <- CountClust::ExtractTopFeatures(IM60_gom_fit_K6$theta, top_features = featurescount,
                                             method = "poisson", options = "max")
K6spp <- t(apply(K6features$indices, c(1,2), function(x) return(rownames(IM60_gom_fit_K6$theta)[x])))
K6spp <- K6spp[,-1] 
colnames(K6spp) <- c("India", "Mainland SEA", "SE China", "Sundaland", "Philippines")
rownames(K6spp) <- c(1:featurescount)

K6spp
write.xlsx2(K6spp, "Ecostructure/K6spp.xlsx")

# One can repeat this process for any K run, for any n number of "top contributing" species 

########## 08 Maps of all motifs ##########
# Write all motif maps into a pdf
motif_names$full_name
motif_names$full_name[which(motif_names$dt_set == "K3")]

motif_names$full_name[which(motif_names$dt_set == paste0("K",K) & 
                              motif_names$motif == paste0("motif_",i))]

motif_names$motif

identical(motif_names$motif_name[which(motif_names$dt_set == "K2" & 
                               motif_names$motif == "motif_1")], character(0))

pdf("Ecostructure/allrunsallmotifs_v5.pdf")
for(K in 2:k_max){
  testcolorraster <- rasterbase
  selected.df <- get(paste0("ecosomega.K",K))
  for(i in 1:K){
    print(paste0("K = ", K, ", area ", i))
    testcolorraster[selected.df[[i]]$cell] <- selected.df[[i]]$color
    if(identical(motif_names$motif_name[which(motif_names$dt_set == paste0("K",K) & 
                                    motif_names$motif == paste0("motif_",i))], character(0)) == T){
      plot(testcolorraster, main=paste0("K = ", K, ", motif ", i),
           xlab="longitude (°)", ylab="latitude (°)")
      plot(raster::intersect(worldmap,testcolorraster),add=T)
    } else{
      plot(testcolorraster, xlab="longitude (°)", ylab="latitude (°)",
           main=paste0("K = ", K, ", motif ", i, "\n", 
                       motif_names$motif_name[which(motif_names$dt_set == paste0("K",K) & 
                                                      motif_names$motif == paste0("motif_",i))], ": ",
                       motif_names$full_name[which(motif_names$dt_set == paste0("K",K) & 
                                                     motif_names$motif == paste0("motif_",i))]))
      plot(raster::intersect(worldmap,testcolorraster),add=T)  
    }
  }
}
dev.off()


########## 09 Inset maps for figure ##########
# inset maps
map.col.pal <- c("#88CCEE", "#DDCC77", "#CC6677", "#117733",
                 "#AA4499", "#332288", "#44AA99", "#882255",
                 "#661100", "#999933", "#888888")

rastcol.india <- sub("#332288", "#88CCEE", ecosomega.K6[[2]]$color)
rastcol.msea <- sub("#332288", "#CC6677", ecosomega.K6[[3]]$color)
rastcol.sunda <- sub("#332288", "#117733", ecosomega.K6[[5]]$color)
rastcol.phi <- sub("#332288", "#999933", ecosomega.K6[[6]]$color)

K6raster2.india <- rasterbase
K6raster2.india[ecosomega.K6[[2]]$cell] <- rastcol.india

K6raster3.msea <- rasterbase
K6raster3.msea[ecosomega.K6[[3]]$cell] <- rastcol.msea

K6raster4.sec <- rasterbase
K6raster4.sec[ecosomega.K6[[4]]$cell] <- ecosomega.K6[[4]]$color

K6raster5.sunda <- rasterbase
K6raster5.sunda[ecosomega.K6[[5]]$cell] <- rastcol.sunda

K6raster6.phi <- rasterbase
K6raster6.phi[ecosomega.K6[[6]]$cell] <- rastcol.phi

pdf(paste0("Ecostructure/ecostructure_figs_tables/panel_", "india", "_v02", ".pdf"), 
    height = 8, width = 8)
plot(K6raster2.india, axes= F); plot(asia_tropic_map,add=T)
dev.off()


pdf(paste0("Ecostructure/ecostructure_figs_tables/panel_", "MSEA", "_v02", ".pdf"), 
    height = 8, width = 8)
plot(K6raster3.msea, axes= F); plot(asia_tropic_map,add=T)
dev.off()

pdf(paste0("Ecostructure/ecostructure_figs_tables/panel_", "SEC", "_v02", ".pdf"), 
    height = 8, width = 8)
plot(K6raster4.sec, axes= F); plot(asia_tropic_map,add=T)
dev.off()

pdf(paste0("Ecostructure/ecostructure_figs_tables/panel_", "sundaland", "_v02", ".pdf"), 
    height = 8, width = 8)
plot(K6raster5.sunda, axes= F); plot(asia_tropic_map,add=T)
dev.off()

pdf(paste0("Ecostructure/ecostructure_figs_tables/panel_", "philippines", "_v02", ".pdf"), 
    height = 8, width = 8)
plot(K6raster6.phi, axes= F); plot(asia_tropic_map,add=T)
dev.off()

