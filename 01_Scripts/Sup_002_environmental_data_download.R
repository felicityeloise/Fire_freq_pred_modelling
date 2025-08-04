# Written by Felicity Charles
# Date:30/07/2025


##### Fire frequency analysis ----
# This script downloads yearly environmental data needed for predicting  fire frequency from online databases

# R version 4.3.1

# 1. Load required packages ----
library(terra) # terra_1.7-78 
library(dplyr) # dplyr_1.1.4

# Install gdalUtilities for working with large datasets
#library(devtools)
#install_github("JoshOBrien/gdalUtilities")
library(gdalUtilities) # gdalUtilities_1.2.5

# 1. Read in, crop and aggregate environmental data ----
# Download Interim Biogeographic Regionalisation for Australia (IBRA) Version 7.1 (Regions) from DCCEEW 

IBRA <- vect('./')

e <- ext(1902033, 2111776, -3257627, -2954985)
SEQ <- as.polygons(e, 'EPSG:3577')
writeVector(SEQ, './00_Data/SEQ_bound/SEQ.gpkg')

SEQ <- vect('./00_Data/SEQ_bound/SEQ.gpkg')


