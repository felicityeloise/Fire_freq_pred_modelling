# Caveat emptor
# Written by Felicity Charles 
# Date: 28/02/2022

##### Fire frequency analysis ----
# This script gathers together data needed for determining the fire frequency from readily available fire history products including the QLD spatial catalogue and TERN.
# R version 4.3.1

# Load required packages
library(fasterize) # fasterize_1.0.5
library(terra) # terra_1.7-78
library(ggplot2) # ggplot2_3.5.1
library(dplyr) # dplyr_1.1.4
library(sf) # sf_1.0-14


# 1. Read in QPWS fire history data ----
# This data can be downloaded from the QLD spatial catalogue

QPWS_fire_hist <- vect('./00_Data/Fire_data/QPWS_fire_history/Fire_history___QPWS.shp')



# 1.2 Subset QPWS fire history data  ----
# Keep only the data from 1987 onwards to match with TERN fire history data

QPWS_fire_hist_1987 <- subset(QPWS_fire_hist, QPWS_fire_hist$OUTYEAR >= 1987)
unique(QPWS_fire_hist_1987$OUTYEAR) # Check that this has worked


# We also need to remove any data that has an outyear as 2024 as we only want data from 1987-2023. 

QPWS_fire_hist_1987 <- subset(QPWS_fire_hist_1987, QPWS_fire_hist_1987$OUTYEAR <2024)
unique(QPWS_fire_hist_1987$OUTYEAR) # Check that this has worked


# Reproject this data, as we know the TERN data is in EPSG:3577

QPWS_fire_hist_1987 <- project(QPWS_fire_hist_1987,'EPSG:3577')

# Save this data as an output
writeVector(QPWS_fire_hist_1987, './00_Data/Fire_data/Outputs/QPWS_fire_hist_1987.gpkg')


QPWS_fire_1987 <- st_read('./00_Data/Fire_data/Outputs/QPWS_fire_hist_1987.gpkg')
QPWS_fire_1987 <- st_transform(QPWS_fire_1987, crs = 'EPSG:3577') # Have to run this transformation again after reading in the data

plot(QPWS_fire_1987)




# 1.3 Rasterize QPWS fire data using fasterize for SEQ only ----
# Get bounding box from SEQ
SEQ <- vect('./00_Data/SEQ_bound/SEQ.gpkg')
SEQ
rtemp <- raster::raster(xmn = 1902033, xmx = 2111776, ymn = -3257627, ymx = -2954985, res = 5, crs = 'EPSG:3577')
QPWS_SEQ_freq_rast <- fasterize(QPWS_fire_1987, rtemp, field = 'OUTYEAR', fun = 'count')

plot(QPWS_SEQ_freq_rast)
raster::writeRaster(QPWS_SEQ_freq_rast, './00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_freq_raster.tif')

QPWS <- rast('./00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_freq_raster.tif')
plot(QPWS)




# 2. Lets create the hydrogaphic features mask for generating random points and masking these areas on the fire frequency map
# 2.1 Read in hydrographic feature data
canal <- vect('./00_Data/Environmental_data/Hydrographic_features/Canal_areas.shp') %>% 
  project('EPSG:3577') %>% 
  crop(SEQ)
writeVector(canal, './00_Data/Environmental_data/Outputs/Hydrographic_features/Canal_SEQ.gpkg')

lake <- vect('./00_Data/Environmental_data/Hydrographic_features/Lakes.shp') %>% 
  project('EPSG:3577') %>% 
  crop(SEQ)
writeVector(lake, './00_Data/Environmental_data/Outputs/Hydrographic_features/Lakes_SEQ.gpkg')


pond <- vect('./00_Data/Environmental_data/Hydrographic_features/Pondage.shp') %>% 
  project('EPSG:3577') %>% 
  crop(SEQ)
writeVector(pond, './00_Data/Environmental_data/Outputs/Hydrographic_features/Ponds_SEQ.gpkg')


reservoir <- vect('./00_Data/Environmental_data/Hydrographic_features/Reservoirs.shp') %>% 
  project('EPSG:3577') %>% 
  crop(SEQ)
writeVector(reservoir, './00_Data/Environmental_data/Outputs/Hydrographic_features/Reservoirs_SEQ.gpkg')


watercourse <- vect('./00_Data/Environmental_data/Hydrographic_features/Watercourse_areas.shp') %>% 
  project('EPSG:3577') %>% 
  crop(SEQ)
writeVector(watercourse, './00_Data/Environmental_data/Outputs/Hydrographic_features/Watercourses_SEQ.gpkg')

# Mask the hydrological features 
QPWS_ff1 <- mask(QPWS_ff, canal, inverse = T)
plot(QPWS_ff1)
QPWS_ff2 <- mask(QPWS_ff1, lake, inverse = T)
plot(QPWS_ff2)
QPWS_ff3 <- mask(QPWS_ff2, pond, inverse = T)
plot(QPWS_ff3)
QPWS_ff4 <- mask(QPWS_ff3, reservoir, inverse = T)
plot(QPWS_ff4)
QPWS_ff5 <- mask(QPWS_ff4, watercourse, inverse = T)
plot(QPWS_ff5)
writeRaster(QPWS_ff5, './00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_freq_hydrographical_mask.tif')



# Attempted to run the following but had to use a remote desktop to do so and took weeks without completing so was abandoned. 
# Rasterize QPWS fire data using fasterize
# Create the raster template 
#rtemp <- raster::raster(xmn = 593219.8, xmx = 2121654, ymn = -3312368, ymx = -1010422, res = 5, crs = 'EPSG:3577') 
#QPWS_fire_1987_rast <-  fasterize(QPWS_fire_1987, rtemp, field = 'OUTYEAR', fun = 'count')  


# Date: 14/04/2023

# 2. Read in the TERN QLD fire history Sentinel data ----
# Landsat fire scar imagery records the earliest month in which a fire was detected where there are multiple fire scars for a given pixel. There is possibly some inaccuracy in the fire frequency calculation given this as we cannot determine if there was a fire that burnt a given pixel more than once in a year.

#  Sentinel-2 analysis does not provide a complete record of fire history. Fire scars may be missed or under-mapped due to: 1) Lack of visibility due to cloud, haze and smoke, and cloud shadow; 2) Misclassification as non-fire related change or cloud shadow; 3) Lack of detection due to size or patchiness. Fire scars smaller than 2 ha may not be included; 4) Lack of detection due to rapid regrowth of vegetation. This is particularly an issue when there have been multiple cloud-affected images in the time series; 5) Lack of detection for cool grass/understorey fires, obscured by unburnt vegetation

# 2.1 What is the spatial extent we are interested in? ----
# Create SpatExtent based on coordinates - we know we need to split this into tiles as the datasets are so large. This might be more splits than necessary but makes life easier to not crash R constantly.

extent1 <- ext(1902030, 2006910, -3030660, -2954970)# Order xmin, xmax, ymin, ymax (standard order for terra)
extent2 <- ext(2006895, 2111805, -3030660, -2954970)
extent3 <- ext(1902030, 2006910, -3106320, -3030630)
extent4 <- ext(2006895, 2111805, -3106320, -3030630)
extent5 <- ext(1902030, 2006910, -3181980, -3106290)
extent6 <- ext(2006895, 2111805, -3181980, -3106290)
extent7 <- ext(1902030, 2006910, -3257640, -3181950) 
extent8 <- ext(2006895, 2111805, -3257640, -3181950)




# 2.2 What years of Sentinel 1 and Sentinel 2 data are we interested in? ----
# Read in the data directly from the data repository, in this case TERN. As the data is read in to R, crop it to the area of interest and reclassify the data. 

## CHANGE CROPPING EXTENT FOR EACH EXTENT AS ABOVE, RUNNING THE BELOW CODE FOR EACH EXTENT

Y1987 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1987_dkaa2.tif") %>%
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254,0))
  

Y1988 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1988_dkaa2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))


Y1989 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1989_dkaa2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))


Y1990 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1990_dkaa2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))


Y1991 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1991_dkaa2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))


Y1992 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1992_dkaa2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))



Y1993 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1993_dkaa2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))


Y1994 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1994_dkaa2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))


Y1995 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1995_dkaa2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))
 


Y1996 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1996_dkaa2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))


Y1997 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1997_dkaa2.tif") %>% 
  crop(extent7) %>%  
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))


Y1998 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1998_dkaa2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))


Y1999 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1999_dkaa2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))


Y2000 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2000_dkaa2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))




Y2001 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2001_dkaa2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))



Y2002 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2002_dkaa2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))


Y2003 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2003_dkaa2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))


Y2004 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2004_dkaa2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))


Y2005 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2005_dkaa2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))


Y2006 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2006_dkaa2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))


Y2007 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2007_dkaa2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))


Y2008 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2008_dkaa2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))




Y2009 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2009_dkaa2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))


Y2010 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2010_dkaa2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))


Y2011 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2011_dkaa2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))

Y2012 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2012_dkaa2.tif") %>%
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254,0))
  


Y2013 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2013_dkda2.tif") %>%
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))


Y2014 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2014_dkda2.tif") %>%
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))


Y2015 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2015_dkda2.tif") %>%
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(254, 0))


Y2016 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2016_dkga2.tif") %>%
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) %>% 
  classify(cbind(34, 0))%>% 
  classify(cbind(255, 0))




Y2017 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/sentinel2/fire_scars/annual_fire_scars/cvmsre_qld_2017_afma2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) 


Y2018 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/sentinel2/fire_scars/annual_fire_scars/cvmsre_qld_2018_afma2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1))



Y2019 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/sentinel2/fire_scars/annual_fire_scars/cvmsre_qld_2019_afma2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1))



Y2020 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/sentinel2/fire_scars/annual_fire_scars/cvmsre_qld_2020_afma2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1)) 



Y2021 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/sentinel2/fire_scars/annual_fire_scars/cvmsre_qld_2021_afma2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1))


Y2022 <- rast("./vsicurl/https://data.tern.org.au/rs/public/data/sentinel2/fire_scars/annual_fire_scars/cvmsre_qld_2022_afma2.tif") %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1))

Y2023 <- rast('./vsicurl/https://data.tern.org.au/rs/public/data/sentinel2/fire_scars/annual_fire_scars/cvmsre_qld_2023_afma2.tif') %>% 
  crop(extent7) %>% 
  classify(cbind(NaN, 0)) %>% 
  classify(cbind(0,0)) %>% 
  classify(cbind(1,1)) %>% 
  classify(cbind(2,1)) %>% 
  classify(cbind(3,1)) %>% 
  classify(cbind(4,1)) %>% 
  classify(cbind(5,1)) %>% 
  classify(cbind(6,1)) %>% 
  classify(cbind(7,1)) %>% 
  classify(cbind(8,1)) %>% 
  classify(cbind(9,1)) %>% 
  classify(cbind(10,1)) %>% 
  classify(cbind(11,1)) %>% 
  classify(cbind(12,1))



# 3. Calculate fire frequency for Sentinel 1 ----
Sentinel1 <- sum(Y1987, Y1988, Y1989, Y1990, Y1991, Y1992, Y1993, Y1994, Y1995, Y1996, Y1997, Y1998, Y1999, Y2000, Y2001, Y2002, Y2003, Y2004, Y2005, Y2006, Y2007, Y2008, Y2009, Y2010, Y2011, Y2012, Y2013, Y2014, Y2015, Y2016)
writeRaster(Sentinel1, './00_Data/Fire_data/Outputs/Sentinel/Extent7/Sentinel1_ff_ext7.tif', overwrite = T)


# 4. Calculate the fire frequency of Sentinel 2 and aggregate the Sentinel 2 data so that it matches Sentinel 1 data resolution ----

Sentinel2 <- sum(Y2017, Y2018, Y2019, Y2020, Y2021, Y2022, Y2023) %>% 
  terra::aggregate(fact = 3)
writeRaster(Sentinel2, './00_Data/Fire_data/Outputs/Sentinel/Extent8/Sentinel2_ff_ext8.tif', overwrite = T)

# Clean the workspace before continuing
rm(list = ls())
gc()

# 5. Combine Sentinel 1 and 2 data fire frequency ----
# Read the data back in 
Sent1_ext1 <- rast('./00_Data/Fire_data/Outputs/Sentinel/Extent1/Sentinel1_ff_ext1.tif')
Sent2_ext1 <- rast('./00_Data/Fire_data/Outputs/Sentinel/Extent1/Sentinel2_ff_ext1.tif')

Sent1_ext2 <- rast('./00_Data/Fire_data/Outputs/Sentinel/Extent2/Sentinel1_ff_ext2.tif')
Sent2_ext2 <- rast('./00_Data/Fire_data/Outputs/Sentinel/Extent2/Sentinel2_ff_ext2.tif')


Sent1_ext3 <- rast('./00_Data/Fire_data/Outputs/Sentinel/Extent3/Sentinel1_ff_ext3.tif')
Sent2_ext3 <- rast('./00_Data/Fire_data/Outputs/Sentinel/Extent3/Sentinel2_ff_ext3.tif')

Sent1_ext4 <- rast('./00_Data/Fire_data/Outputs/Sentinel/Extent4/Sentinel1_ff_ext4.tif')
Sent2_ext4 <- rast('./00_Data/Fire_data/Outputs/Sentinel/Extent4/Sentinel2_ff_ext4.tif')

Sent1_ext5 <- rast('./00_Data/Fire_data/Outputs/Sentinel/Extent5/Sentinel1_ff_ext5.tif')
Sent2_ext5 <- rast('./00_Data/Fire_data/Outputs/Sentinel/Extent5/Sentinel2_ff_ext5.tif')


Sent1_ext6 <- rast('./00_Data/Fire_data/Outputs/Sentinel/Extent6/Sentinel1_ff_ext6.tif')
Sent2_ext6 <- rast('./00_Data/Fire_data/Outputs/Sentinel/Extent6/Sentinel2_ff_ext6.tif')

Sent1_ext7 <- rast('./00_Data/Fire_data/Outputs/Sentinel/Extent7/Sentinel1_ff_ext7.tif')
Sent2_ext7 <- rast('./00_Data/Fire_data/Outputs/Sentinel/Extent7/Sentinel2_ff_ext7.tif')
plot(Sent1_ext7)
plot(Sent2_ext7)

Sent1_ext8 <- rast('./00_Data/Fire_data/Outputs/Sentinel/Extent8/Sentinel1_ff_ext8.tif')
Sent2_ext8 <- rast('./00_Data/Fire_data/Outputs/Sentinel/Extent8/Sentinel2_ff_ext8.tif')





# Now when we put the two extents for each 'row' of the Sentinel image together
Sentinel_ff_ext1 <- mosaic(Sent1_ext1, Sent2_ext1, fun = "sum")
writeRaster(Sentinel_ff_ext1, './00_Data/Fire_data/Outputs/Sentinel/Extent1/Sentinelff_ext1.tif')

Sentinel_ff_ext2 <- mosaic(Sent1_ext2, Sent2_ext2, fun = "sum")
writeRaster(Sentinel_ff_ext2, './00_Data/Fire_data/Outputs/Sentinel/Extent2/Sentinelff_ext2.tif')

Sentinel_ff_ext3 <- mosaic(Sent1_ext3, Sent2_ext3, fun = "sum")
writeRaster(Sentinel_ff_ext3, './00_Data/Fire_data/Outputs/Sentinel/Extent3/Sentinelff_ext3.tif')

Sentinel_ff_ext4 <- mosaic(Sent1_ext4, Sent2_ext4, fun = "sum")
writeRaster(Sentinel_ff_ext4, './00_Data/Fire_data/Outputs/Sentinel/Extent4/Sentinelff_ext4.tif')

Sentinel_ff_ext5 <- mosaic(Sent1_ext5, Sent2_ext5, fun = "sum")
writeRaster(Sentinel_ff_ext5, './00_Data/Fire_data/Outputs/Sentinel/Extent5/Sentinelff_ext5.tif')

Sentinel_ff_ext6 <- mosaic(Sent1_ext6, Sent2_ext6, fun = "sum")
writeRaster(Sentinel_ff_ext6, './00_Data/Fire_data/Outputs/Sentinel/Extent6/Sentinelff_ext6.tif')

Sentinel_ff_ext7 <- mosaic(Sent1_ext7, Sent2_ext7, fun = "sum")

plot(Sentinel_ff_ext7)
writeRaster(Sentinel_ff_ext7, './00_Data/Fire_data/Outputs/Sentinel/Extent7/Sentinelff_ext7.tif')

Sentinel_ff_ext8 <- mosaic(Sent1_ext8, Sent2_ext8, fun = "sum")
writeRaster(Sentinel_ff_ext8, './00_Data/Fire_data/Outputs/Sentinel/Extent8/Sentinelff_ext8.tif')


# 6. Create QLD Sentinel fire frequency data layer ----
# Merge the extents into one layer, here we will use merge rather than mosaic as we keep getting a vertical line when we use mosiac on the right hand side extents that does not appear before we get to this step

Sentinel_ff <- merge(Sentinel_ff_ext1, Sentinel_ff_ext2, Sentinel_ff_ext3, Sentinel_ff_ext4, Sentinel_ff_ext5, Sentinel_ff_ext6, Sentinel_ff_ext7, Sentinel_ff_ext8)
plet(Sentinel_ff)
Sentinel_ff

writeRaster(Sentinel_ff, overwrite = T, './00_Data/Fire_data/Outputs/Sentinel/Sentinel_ff.tif')
# Clean the workspace before continuing
rm(list = ls())
gc()


# Mask Sentinel_ff by the hydrological features, masking by one hydrographic feature and then the next
Sentinel_ff <- rast('./00_Data/Fire_data/Outputs/Sentinel/Sentinel_ff.tif')
Sentinel_ff1 <- mask(Sentinel_ff, canal, inverse = T)
plot(Sentinel_ff1)
Sentinel_ff2 <- mask(Sentinel_ff1, lake, inverse = T)
plot(Sentinel_ff2)
Sentinel_ff3 <- mask(Sentinel_ff2, pond, inverse = T)
plot(Sentinel_ff3)
Sentinel_ff4 <- mask(Sentinel_ff3, reservoir, inverse = T)
plot(Sentinel_ff4)
Sentinel_ff5 <- mask(Sentinel_ff4, watercourse, inverse = T)
plot(Sentinel_ff5)


unique(Sentinel_ff5$lztmre_qld_1987_dkaa2)

writeRaster(Sentinel_ff5, './00_Data/Fire_data/Outputs/Sentinel/Sentinel_ff_hydrographical_mask_SEQ.tif')

# Use focal to fill in any areas of Sentinel data that have NA values

Sent_foc <- terra::focal(Sentinel_ff5, fun = "mean", na.policy = "only", na.rm = T)
Sent_foc # The same range of values
names(Sent_foc) <- "Sentinel_fire_freq"
Sent_foc <- round(Sent_foc)
writeRaster(Sent_foc, './00_Data/Fire_data/Outputs/Sentinel/Sentinel_ff_hydrographical_mask_SEQ_focal.tif')

plot(Sent_foc)


# 7. Compare the two fire frequency datasets ----
# Now that we have fire frequency calculate for both datasets and they are both raster files, we want to compare these datasets to see how well correlated they are.
Sentinel_ff <- rast("./00_Data/Fire_data/Outputs/Sentinel/Sentinel_ff_hydrographical_mask_SEQ_focal.tif")
QPWS_SEQ_ff <- rast("./00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_freq_hydrographical_mask.tif")

# 7.1 Limit QPWS data to QPWS estates
protected_land <- vect('./00_Data/Protected_areas/Protected_areas.shp')
protected_land <- project(protected_land, y = 'EPSG:3577')
protected_land <- crop(protected_land, QPWS_SEQ_ff)

QPWS_SEQ_ff_mask <- mask(QPWS_SEQ_ff, protected_land) # Restrict QPWS fire frequency data to only those areas that are managed by QPWS
QPWS_SEQ_ff_mask # Check how this looks 

plot(QPWS_SEQ_ff_mask)

# 7.2 Transect level fire frequency data ----
Transects <- vect('./00_Data/Transects/Transects.shp')
Transects <- project(Transects, y = 'EPSG:3577')

Fire_transects <- terra::extract(Sentinel_ff, Transects)
colnames(Fire_transects) <- c("ID", "Sentinel_fire_freq")

QPWS <- terra::extract(QPWS_SEQ_ff_mask, Transects)
Fire_transects$QPWS_SEQ_ff <- QPWS$QPWS_SEQ_freq_rast

# Add in transect information to this dataframe
Fire_transects$T_no <- Transects$Name
Fire_transects$Location <- Transects$Location
Fire_transects 
# Some QPWS values have NA values as they fall just outside the boundary of the polygon of the QPWS estate. In some cases this is because the point falls close to where a road is located on the ground with the boundary not quite all the way up to this point, so some data is lost here. 

st_write(Fire_transects, './00_Data/Fire_data/Outputs/Fire_frequency_transects.csv')


# 7.3 Convert to long data format for plotting ----
Fire_transects_long <- cbind(Fire_transects[4:5], stack(Fire_transects[2:3]))
colnames(Fire_transects_long) <- c("T_no", "Location", "Frequency", "Dataset")
Fire_transects_long$Frequency <- round(Fire_transects_long$Frequency, digits = 0)
st_write(Fire_transects_long, './00_Data/Fire_data/Outputs/Fire_frequency_transects_long.csv')
Fire_transects_long <- read.csv('./00_Data/Fire_data/Outputs/Fire_frequency_transects_long.csv')



# 7.4 Explore the data using plots ----
boxplot(Frequency ~ Dataset, data = Fire_transects_long) # Sentinel has lower values overall

ggplot(data = Fire_transects_long, aes(x = Dataset, y = Frequency, fill = Dataset, col = Dataset))+
  geom_violin()

ggplot(data = Fire_transects_long, 
       aes(x = Dataset, 
           y = Frequency, 
           fill = Dataset,
           col = Dataset))+
  geom_boxplot() +
  geom_jitter(col = "black", size = 0.7, alpha = 0.4)



# Test the correlation between the values

# Lets make a simple scatterplot
ggplot(Fire_transects_long) +
  aes(x = Frequency, y = Location, fill = Dataset, col = Dataset)+
  geom_point(position = position_dodge(.5)) +
  theme_minimal()


plot(Fire_transects_long$Frequency[Fire_transects_long$Dataset == "Sentinel_fire_freq"], Fire_transects_long$Frequency[Fire_transects_long$Dataset == "QPWS_SEQ_ff"], xlim = c(0,7), ylim = c(0,7), xlab = "Sentinel fire frequency", ylab = "QPWS fire frequency", pch = 20)


summary(Fire_transects_long$Frequency[Fire_transects_long$Dataset == "Sentinel_fire_freq"])

summary(Fire_transects_long$Frequency[Fire_transects_long$Dataset == "QPWS_SEQ_ff"])
# We know some of the QPWS values are supposed to be 0 such as those that are not National parks like WATER1, Gil, Bul, and Bar, HV, some LTLs as well but not all

Fire_transect_cor <- cor.test(Fire_transects$Sentinel_fire_freq, Fire_transects$QPWS_SEQ_ff)
Fire_transect_cor # Vary in the same direction, but a weak correlation between the datasets. 



# 8. Compare accuracy of Sentinel and QPWS fire frequency data based on a random sample of points ----
# 8.1 Create 10,000 random points ----
# We want to create the 10,000 random points using the QPWS data as vector data that is masked by protected areas and cropped to the Sentinel extent as we want to have 10,000 random points where we know that we should have QPWS data that is accurate. To do this we also need to make sure we only have the data from 1987 onwards

# Take the protected areas of SEQ and produce 10,000 random points in these areas. 
# Read in the data that will be used to crop and mask the protected areas. We want to crop it by SEQ land masses as points need to fall on land. 
SEQ <- vect('./00_Data/SEQ_bound/SEQ.gpkg')
Aus <- vect('./00_Data/Australia_shapefile/STE11aAust.shp') %>% 
  project("EPSG:3577")
canal <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Canal_SEQ.gpkg')
lake <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Lakes_SEQ.gpkg')
pond <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Ponds_SEQ.gpkg')
reservoir <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Reservoirs_SEQ.gpkg')
watercourse <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Watercourses_SEQ.gpkg')


protected_land <-  vect('./00_Data/Protected_areas/Protected_areas.shp') %>% 
  project('EPSG:3577') %>% 
  crop(Aus) %>% 
  crop(SEQ) %>% 
  erase(canal) %>%
  erase(lake) %>%
  erase(pond) %>%
  erase(reservoir) %>%
  erase(watercourse)

rm(Aus, canal, lake, pond, reservoir, SEQ, watercourse)
gc()

# Convert to a simple features object
protectedland <- sf::st_as_sf(protected_land, coords = c(protected_land$LONGITUDE, protected_land$LATITUDE), crs = 'EPSG:3577')


# For the purposes of predictive modelling we also need to consider that we do not know if areas that lack QPWS fire frequency data are true 0s (areas that have had 0 fires since 1987) or if these areas are missing data. We cannot confidently just reassign any NA values as 0 as this may be inaccurate. We also only want to produce presence points in areas that have a fire frequency from 1- 11.
QPWS_fire <- mask(QPWS_SEQ_ff, protected_land) # Limit the QPWS fire data to QPWS land
plet(QPWS_fire) # Check this looks right 

set.seed(33)
QPWSrand <- spatSample(QPWS_fire, 31197, na.rm = T, xy = T) # To actually produce 10,000 random points we need to request a larger number of points. This will return a warning message that fewer cells were returned than requested. Another issue we have is that when we then convert this to a spatvector we have more points than what this returns so we need to account for this as we go as well
unique(QPWSrand$QPWS_SEQ_freq_raster) # Check the fire frequencies that are being sampled, would like to capture most of the range from 1-13. With those at higher frequencies likely to be much less prevalent and probably small patches so less likely to be captured even by these 10,000 random points
QPWS_rand <- vect(QPWSrand, geom = c('x', 'y'), crs = "EPSG:3577") # Convert back to a terra spatVector data format
length(QPWS_rand) # Check the number of points 
plet(QPWS_rand) # Check this looks right
writeVector(QPWS_rand, './00_Data/Fire_data/Outputs/QPWS_random.gpkg', overwrite = T)

QPWS_rand <- vect('./00_Data/Fire_data/Outputs/QPWS_random.gpkg')



# 8.2 Extract values for 10,000 random points ----
Rand_fire <- terra::extract(Sentinel_ff, QPWS_rand)
QPWS_ff_rand <- terra::extract(QPWS_SEQ_ff, QPWS_rand)
unique(QPWS_ff_rand$QPWS_SEQ_freq_raster) # As we have limited the points to areas with known fire history (areas that had burnt) using the raster dataset we should have no NA values. 

range(unique(Rand_fire$Sentinel_fire_freq))
unique(Rand_fire$Sentinel_fire_freq)
unique(is.na(Rand_fire$Sentinel_fire_freq)) # We can see that our focal spatraster calculation step for NA values has worked as we are returned no NA values. 

# Combine the data into one dataset
Rand_fire$QPWS_ff_rand <- QPWS_ff_rand$QPWS_SEQ_freq_rast
colnames(Rand_fire) <- c("ID", "Sentinel_rand_firefreq", "QPWS_rand_firefreq")
head(Rand_fire)
str(Rand_fire) # The data looks as expected

# Add the latitude and longitude to the dataset
# Get the coords of QPWS_rand
coords <- crds(QPWS_rand)
head(coords)
Rand_fire$Lon <- coords[,1]
Rand_fire$Lat <- coords[,2]
head(Rand_fire)
str(Rand_fire)


write.csv(Rand_fire, './00_Data/Fire_data/Outputs/Random_points_data/Fire_frequency_random_points_QPWS_pres.csv')


# 8.4 Convert the data to long format ----
Random_fire <- cbind(Rand_fire[1], stack(Rand_fire[2:3]), Rand_fire[4:5])
tail(Random_fire)
colnames(Random_fire) <- c("ID", "Frequency", "Dataset", "Lon", "Lat")
head(Random_fire); tail(Random_fire); dim(Random_fire)
write.csv(Random_fire, './00_Data/Fire_data/Outputs/Random_points_data/Fire_frequency_randompoints_longformat_QPWS_presences.csv')


Rand_fire <- read.csv('./00_Data/Fire_data/Outputs/Random_points_data/Fire_frequency_random_points_QPWS_pres.csv', header = T)
head(Rand_fire); dim(Rand_fire)
Random_fire <- read.csv('./00_Data/Fire_data/Outputs/Random_points_data/Fire_frequency_randompoints_longformat_QPWS_presences.csv', header = T)
head(Rand_fire); dim(Rand_fire)


# Explore the data to test for correlations ----

# Lets make a simple scatterplot
ggplot(Random_fire) +
  aes(x = Frequency, y = ID, fill = Dataset, col = Dataset)+
  geom_point() +
  theme_minimal()

Rand <- cor.test(Rand_fire$Sentinel_rand_firefreq, Rand_fire$QPWS_rand_firefreq)
Rand # Vary in the same direction but weak correlation between these two variables. 

# As these two datasets are not very well correlated can we predict from QPWS for areas outside of QPWS estates?


head(Random_fire)
plot(Random_fire$Frequency[Random_fire$Dataset == "Sentinel_rand_ff"], Random_fire$Frequency[Random_fire$Dataset == "QPWS_rand_ff"], pch = 20, xlab = "Sentinel fire frequency", ylab = "QPWS fire frequency", ylim = c(0,15))

m1 <- lm(QPWS_rand_ff ~ Sentinel_rand_ff, data = Rand_fire)
summary(m1)
head(Rand_fire)

predict(m1)
p1 <- data.frame(Sentinel_rand_ff = seq(0,15, length.out = 100))
head(p1)
p1b <- predict(m1, newdata = p1, se.fit = T)
p1c <- p1
p1c$fit <- p1b$fit # Predicted QPWS 
p1c$se <- p1b$se.fit
head(p1c)
p1c$lci <- p1c$fit - (p1c$se*1.96)
p1c$uci <- p1c$fit + (p1c$se*1.96)

plot(p1c$Sentinel_rand_ff, p1c$fit, ylim = c(min(Rand_fire$QPWS_rand_ff), max(Rand_fire$QPWS_rand_ff)), type = "l", xlab = "Predicted Sentinel fire frequency", ylab = "Predicted QPWS fire frequency")
lines(p1c$Sentinel_rand_ff, p1c$lci, lty = 2)
lines(p1c$Sentinel_rand_ff, p1c$uci, lty = 2)

# This predictive model works well at low fire frequencies but accuracy is lost as fire frequency increases. 

ggplot(data = Random_fire, 
       aes(x = Dataset, 
           y = Frequency, 
           fill = Dataset,
           col = Dataset))+
  geom_boxplot() +
  geom_jitter(col = "black", size = 0.7, alpha = 0.4)


# Take a look at the spatial autocorrelation


m2 <- gls(QPWS_rand_ff ~ Sentinel_rand_ff, data = Rand_fire)
vario <- Variogram(m2, form = ~ Lon + Lat)
plot(vario, smooth = T, ylim = c(0,1.3)) # Semi-variance is decreasing with distance. 


# Our problem now is that we are working on public and private land but the fire data available for private land is not adequate, it tends to underestimate in areas we know that have high fire frequency in the QPWS estates. Can we improve this by adding in environmental data?
# Landsat fire scar imagery records the earliest month in which a fire was detected where there are multiple fire scars for a given pixel. There is possibly some inaccuracy in the fire frequency calculation given this as we cannot determine if there was a fire that burnt a given pixel more than once in a year.