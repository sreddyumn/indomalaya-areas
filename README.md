# indomalaya-areas
This repository contains a compilation of code, data, and visualizations associated with the paper "Biogeographic regions within Indomalaya: an integrated approach using bird distributions."
Authors: Zheng Oong, Sharon A. Jansa, Sushma Reddy

## R scripts ##
1. **01_IndomalayaBioregions_Dataset.R**: Dataset setup and filtering steps (including information about work done externally on ArcGIS and QGIS), plus presence-absence matrix generation with letsR.
2. **02_IndomalayaBioregions_NDM.R**: NDM input file generation and output management.
3. **03_IndomalayaBioregions_Ecostructure.R**: Ecostructure Grade of Membership model fitting, visualization of results as rasters, calculation of motif overlap, and network graph for motif stability.
4. **04_IndomalayaBioregions_hclust.R**: Hierarchical clustering algorithm selection, assessment of optimal cluster, visualization of cluster assignments on dendrogram and map.

## Data files ##
The 'data' directory contains the following files (listed alphabetically):
1. **any_indomalaya.csv**: A list of all 2827 species with any native and resident occurrence in Indomalaya. This list is the starting point for the area percentage filtering process in described in the script _01_IndomalayaBioregions_Dataset.R_.
2. **IM_Presab_XXdeg.csv**: Four presence-absence matrices generated at four spatial resolutions from the final dataset of 1544 study species. The end of the filenames indicate the spatial resolution of the matrix: 5deg, 2deg, 1deg, and halfdeg for 5°, 2°, 1°, and 0.5° respectively. Rows are individual cells, while columns are the 1544 study species plus two additional ones for longitude (x) and latitude (y).
3. **motif_names_v4.csv**: A list of all motif names (plus three-letter abbreviations) from Ecostructure, plus a column containing a series of "graph index" for use in plotting the network graph tracking motif stability.

In addition, the 'data' directory contains the following subdirectories (listed alphabetically):
1. **base_maps**: Contains the Tropical Asia and Indomalaya maps (shapefile) used as base maps for plotting purposes.
2. **ecostructure_model_fits**: Contains all 19 model fit results (K = 2 to 20) from Ecostructure in RDS format.
3. **hclust**: Contains tables (.xlsx) detailing the cophenetic correlation coefficients calculated for algorithm selection, as well as the silhouette widths and matrix correlation values for selecting optimal number(s) of clusters; also contains image files (.png) showing the fusion level, silhouette widths, and matrix correlation values for selecting optimal clusters (based on guidelines in [Borcard et al., 2018](https://doi.org/10.1007/978-3-319-71404-2_4)).

