# Caveat emptor
# Written by Felicity Charles 
# Date: 28/02/2022

##### Fire frequency analysis ----
# This script gathers together data needed for determining the fire frequency from readily available fire history products including the QLD spatial catalogue and TERN.
# R version 4.3.1 # update to 4.5.1

# Load required packages
library(terra) # terra_1.7-78 # updated to 1.8-60
library(dplyr) # dplyr_1.1.4
library(sf) # sf_1.0-14 # updated to 1.0-21
sessionInfo()

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

#plot(QPWS_fire_1987)




# 1.3 Rasterize QPWS fire data for SEQ only ----
# Download Interim Biogeographic Regionalisation for Australia (IBRA) Version 7.1 (Regions) from DCCEEW 
IBRA <- vect('./00_Data/Environmental_data/IBRA/Interim_Biogeographic_Regionalisation_for_Australia_(IBRA)_Version_7.1_(Regions).shp') %>% 
  project('EPSG:3577') # SEQ Bioregion extends into northern NSW, need to crop at Queensland border as TERN fire history mapping does not extend to northern NSW


Aus <- download.file("https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs-edition-3/jul2021-jun2026/access-and-downloads/digital-boundary-files/STE_2021_AUST_SHP_GDA2020.zip", destfile = './00_Data/Spatial_data/Australia.zip', mode = "wb", cacheOK = F)
unzip(zipfile = './00_Data/Spatial_data/Australia.zip', exdir = './00_Data/Spatial_data/Australia')
Aus <- vect('./00_Data/Spatial_data/Australia/STE_2021_AUST_GDA2020.shp') %>%
  project("EPSG:3577")

QLD <- Aus[Aus$STATE_NAME == "Queensland"]
SEQ <- IBRA[IBRA$REG_NAME_7 == "South Eastern Queensland"] %>% 
  crop(QLD)

writeVector(SEQ, './00_Data/SEQ_bound/SEQ_IBRA.gpkg')

SEQ <- vect('./00_Data/SEQ_bound/SEQ_IBRA.gpkg')

# Get bounding box from SEQ
SEQ
rtemp <- raster(xmn = 1881028, xmx = 2121594, ymn = -3442593, ymx = -2656210, res = 5, crs = 'EPSG:3577')

# Very slow to complete, took approximately 48 hours on MacBook Pro with M4 Pro chip and 24 GB memory
QPWS_SEQ_freq_rast <- rasterize(QPWS_fire_1987, rtemp, field = 'OUTYEAR', fun = 'count')
# Fails to run with new IBRA bioregion extent as it is too large QPWS_SEQ_freq_rast <- fastrize(QPWS_fire_1987, rtemp, field = 'OUTYEAR', fun = 'count')

plot(QPWS_SEQ_freq_rast)
QPWS_SEQ_freq_rast
raster::writeRaster(QPWS_SEQ_freq_rast, './00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_IBRA_freq_raster.tif')

QPWS_ff <- rast('./00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_IBRA_freq_raster.tif')
plot(QPWS_ff)




# 2. Lets create the hydrogaphic features mask for generating random points and masking these areas on the fire frequency map
# 2.1 Read in hydrographic feature data
canal <- vect('./00_Data/Environmental_data/Hydrographic_features/Canal_areas.shp') %>% 
  project('EPSG:3577') %>% 
  crop(SEQ)
writeVector(canal, './00_Data/Environmental_data/Outputs/Hydrographic_features/Canal_SEQ_IBRA.gpkg')

lake <- vect('./00_Data/Environmental_data/Hydrographic_features/Lakes.shp') %>% 
  project('EPSG:3577') %>% 
  crop(SEQ)
writeVector(lake, './00_Data/Environmental_data/Outputs/Hydrographic_features/Lakes_SEQ_IBRA.gpkg')


pond <- vect('./00_Data/Environmental_data/Hydrographic_features/Pondage.shp') %>% 
  project('EPSG:3577') %>% 
  crop(SEQ)
writeVector(pond, './00_Data/Environmental_data/Outputs/Hydrographic_features/Ponds_SEQ_IBRA.gpkg')


reservoir <- vect('./00_Data/Environmental_data/Hydrographic_features/Reservoirs.shp') %>% 
  project('EPSG:3577') %>% 
  crop(SEQ)
writeVector(reservoir, './00_Data/Environmental_data/Outputs/Hydrographic_features/Reservoirs_SEQ_IBRA.gpkg')


watercourse <- vect('./00_Data/Environmental_data/Hydrographic_features/Watercourse_areas.shp') %>% 
  project('EPSG:3577') %>% 
  crop(SEQ)
writeVector(watercourse, './00_Data/Environmental_data/Outputs/Hydrographic_features/Watercourses_SEQ_IBRA.gpkg')

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
writeRaster(QPWS_ff5, './00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_IBRA_freq_hydrographical_mask.tif')



# 2. Read in the TERN QLD fire history Sentinel data ----
# Landsat fire scar imagery records the earliest month in which a fire was detected where there are multiple fire scars for a given pixel. There is possibly some inaccuracy in the fire frequency calculation given this as we cannot determine if there was a fire that burnt a given pixel more than once in a year.

#  Sentinel-2 analysis does not provide a complete record of fire history. Fire scars may be missed or under-mapped due to: 1) Lack of visibility due to cloud, haze and smoke, and cloud shadow; 2) Misclassification as non-fire related change or cloud shadow; 3) Lack of detection due to size or patchiness. Fire scars smaller than 2 ha may not be included; 4) Lack of detection due to rapid regrowth of vegetation. This is particularly an issue when there have been multiple cloud-affected images in the time series; 5) Lack of detection for cool grass/understorey fires, obscured by unburnt vegetation

# 2.1 What is the spatial extent we are interested in? ----
SEQ
# In the example case, we are interested in the SEQ IBRA region of Australia with a spatial extent of 1881028, 2121562, -3247627, -2660998  
# In some instances it may be required to split the spatial extent of interest into smaller tiles which can be created using other GIS visualisation software or through terra::makeTiles(). 

# 2.2 What years of Sentinel 1 and Sentinel 2 data are we interested in? ----
# Read in the data directly from the data repository, in this case TERN. As the data is read in to R, crop it to the area of interest and reclassify the data. 
# NOTE: This process cannot be automated as there is a requirement to log on to the TERN data portal prior to running the code to download the datasets.

Y1987 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1987_dkaa2.tif") %>%
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)

Y1988 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1988_dkaa2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y1989 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1989_dkaa2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)

Y1990 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1990_dkaa2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y1991 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1991_dkaa2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y1992 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1992_dkaa2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)



Y1993 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1993_dkaa2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y1994 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1994_dkaa2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y1995 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1995_dkaa2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y1996 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1996_dkaa2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y1997 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1997_dkaa2.tif") %>% 
  crop(SEQ) %>%  
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)

Y1998 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1998_dkaa2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y1999 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_1999_dkaa2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y2000 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2000_dkaa2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y2001 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2001_dkaa2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y2002 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2002_dkaa2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y2003 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2003_dkaa2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y2004 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2004_dkaa2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y2005 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2005_dkaa2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)

Y2006 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2006_dkaa2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y2007 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2007_dkaa2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)

Y2008 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2008_dkaa2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y2009 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2009_dkaa2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y2010 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2010_dkaa2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)

Y2011 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2011_dkaa2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y2012 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2012_dkaa2.tif") %>%
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y2013 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2013_dkda2.tif") %>%
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y2014 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2014_dkda2.tif") %>%
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y2015 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2015_dkda2.tif") %>%
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    254,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y2016 <- rast("/vsicurl/https://data.tern.org.au/rs/public/data/landsat/burnt_area/qld_annual/lztmre_qld_2016_dkga2.tif") %>%
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1,
    34,0,
    255,0
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y2017 <- rast("/vsicurl/https://data.tern.org.au//rs/public/data/sentinel2/fire_scars/annual_fire_scars/cvmsre_qld_2017_afma2.tif") %>%  
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)

Y2018 <- rast("/vsicurl/https://data.tern.org.au//rs/public/data/sentinel2/fire_scars/annual_fire_scars/cvmsre_qld_2018_afma2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y2019 <- rast("/vsicurl/https://data.tern.org.au//rs/public/data/sentinel2/fire_scars/annual_fire_scars/cvmsre_qld_2019_afma2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)



Y2020 <- rast("/vsicurl/https://data.tern.org.au//rs/public/data/sentinel2/fire_scars/annual_fire_scars/cvmsre_qld_2020_afma2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y2021 <- rast("/vsicurl/https://data.tern.org.au//rs/public/data/sentinel2/fire_scars/annual_fire_scars/cvmsre_qld_2021_afma2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y2022 <- rast("/vsicurl/https://data.tern.org.au//rs/public/data/sentinel2/fire_scars/annual_fire_scars/cvmsre_qld_2022_afma2.tif") %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


Y2023 <- rast('/vsicurl/https://data.tern.org.au//rs/public/data/sentinel2/fire_scars/annual_fire_scars/cvmsre_qld_2023_afma2.tif') %>% 
  crop(SEQ) %>% 
  classify(matrix(c(
    0,0,
    1,1,
    2,1,
    3,1,
    4,1,
    5,1,
    6,1,
    7,1,
    8,1,
    9,1,
    10,1,
    11,1,
    12,1
  ), ncol = 2, byrow = T)) %>% 
  subst(NA,0)


# 3. Calculate fire frequency for Sentinel 1 ----
Sentinel1 <- sum(Y1987, Y1988, Y1989, Y1990, Y1991, Y1992, Y1993, Y1994, Y1995, Y1996, Y1997, Y1998, Y1999, Y2000, Y2001, Y2002, Y2003, Y2004, Y2005, Y2006, Y2007, Y2008, Y2009, Y2010, Y2011, Y2012, Y2013, Y2014, Y2015, Y2016)
writeRaster(Sentinel1, './00_Data/Fire_data/Outputs/Sentinel/SEQ_IBRA/Sentinel1_ff_SEQ_IBRA.tif', overwrite = T)


# 4. Calculate the fire frequency of Sentinel 2 and aggregate the Sentinel 2 data so that it matches Sentinel 1 data resolution ----

Sentinel2 <- sum(Y2017, Y2018, Y2019, Y2020, Y2021, Y2022, Y2023) %>% 
  terra::aggregate(fact = 3) # What type of aggregation here! Maybe max??
writeRaster(Sentinel2, './00_Data/Fire_data/Outputs/Sentinel/SEQ_IBRA/Sentinel2_ff_SEQ_IBRA.tif', overwrite = T)


# 5. Combine Sentinel 1 and 2 data fire frequency ----
Sentinel_ff <- mosaic(Sentinel1, Sentinel2, fun = "sum")
Sentinel_ff
plot(Sentinel_ff)
writeRaster(Sentinel_ff, './00_Data/Fire_data/Outputs/Sentinel/SEQ_IBRA/Sentinel_ff_SEQ_IBRA.tif')


# 6. Mask Sentinel_ff by the hydrological features ----
# Masking by one hydrographic feature and then the next
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


unique(Sentinel_ff5$layer)

writeRaster(Sentinel_ff5, './00_Data/Fire_data/Outputs/Sentinel/SEQ_IBRA/Sentinel_ff_hydrographical_mask_SEQ_IBRA.tif')



# Use focal to fill in any areas of Sentinel data that have NA values

Sent_foc <- terra::focal(Sentinel_ff5, fun = "max", na.policy = "only", na.rm = T)
Sent_foc # The same range of values
names(Sent_foc) <- "Sentinel_fire_freq"
plot(Sent_foc)

# Remove data for areas that are outside the SEQ IBRA polygon
Sent_foc <- mask(Sent_foc, SEQ)
writeRaster(Sent_foc, './00_Data/Fire_data/Outputs/Sentinel/SEQ_IBRA/Sentinel_ff_hydrographical_mask_SEQ_IBRA_focal.tif')

# 7. Compare the two fire frequency datasets ----
# Now that we have fire frequency calculate for both datasets and they are both raster files, we want to compare these datasets to see how well correlated they are.
QPWS_SEQ_IBRA_ff <- rast("./00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_IBRA_freq_hydrographical_mask.tif")


# 7.1 Limit QPWS data to QPWS estates
canal <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Canal_SEQ_IBRA.gpkg')
lake <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Lakes_SEQ_IBRA.gpkg')
pond <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Ponds_SEQ_IBRA.gpkg')
reservoir <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Reservoirs_SEQ_IBRA.gpkg')
watercourse <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Watercourses_SEQ_IBRA.gpkg')
SEQ <- vect('./00_Data/SEQ_bound/SEQ_IBRA.gpkg')

protected_land <- vect('./00_Data/Protected_areas/Protected_areas.shp') %>% 
  project('EPSG:3577') %>% 
  crop(SEQ) %>% 
  erase(canal) %>% 
  erase(lake) %>%
  erase(pond) %>%
  erase(reservoir) %>%
  erase(watercourse)
plot(protected_land)


QPWS_SEQ_ff_mask <- mask(QPWS_SEQ_IBRA_ff, protected_land) # Restrict QPWS fire frequency data to only those areas that are managed by QPWS
QPWS_SEQ_ff_mask # Check how this looks 
plot(QPWS_SEQ_ff_mask)

rm(canal, lake, pond, reservoir, watercourse, SEQ, protected_land, QPWS_SEQ_IBRA_ff)
gc()

# 8. Compare accuracy of Sentinel and QPWS fire frequency data based on a random sample of points ----
# 8.1 Create 10,000 random points ----
# We want to create the 10,000 random points using the QPWS data as vector data that is masked by protected areas and cropped to the Sentinel extent as we want to have 10,000 random points where we know that we should have QPWS data that is accurate. While stratified sampling would be more appropriate to gain coverage of all fire frequencies, the large nature of the raster file does not work well with this type of sampling. 
# Sample 1 million values and find the unique fire frequencies
set.seed(42)
QPWSrand <- spatSample(QPWS_SEQ_ff_mask$layer, size = 1e6, method = "random", na.rm = TRUE, as.points = T)
unique(QPWSrand$layer) # Fire frequencies over 10 are just so rare and small in the landscape that even 1 million random points cannot capture these frequencies 

# Now cut down to 10,000 points
set.seed(41)
QPWS_rand <- QPWSrand[sample(1:nrow(QPWSrand), 10000), ]
unique(QPWS_rand$layer) # Make sure with the set.seed() that when we subsample we still keep the range of fire frequencies.
length(QPWS_rand) # Check the number of points 
plet(QPWS_rand) # Check this looks right
writeVector(QPWS_rand, './00_Data/Fire_data/Outputs/QPWS_random_SEQ_IBRA.gpkg', overwrite = T)

QPWS_rand <- vect('./00_Data/Fire_data/Outputs/QPWS_random_SEQ_IBRA.gpkg')



# 8.2 Extract values for 10,000 random points ----
Sentinel_ff <- rast('./00_Data/Fire_data/Outputs/Sentinel/SEQ_IBRA/Sentinel_ff_hydrographical_mask_SEQ_IBRA_focal.tif')
Rand_fire <- terra::extract(Sentinel_ff, QPWS_rand)
QPWS_ff_rand <- terra::extract(QPWS_SEQ_ff_mask, QPWS_rand)
unique(QPWS_ff_rand$layer) # As we have limited the points to areas with known fire history (areas that had burnt) using the raster dataset we should have no NA values. 

range(unique(Rand_fire$Sentinel_fire_freq))
unique(Rand_fire$Sentinel_fire_freq)
unique(is.na(Rand_fire$Sentinel_fire_freq)) # We can see that our focal spatraster calculation step for NA values has worked as we are returned no NA values. 

# Combine the data into one dataset
Rand_fire$QPWS_ff_rand <- QPWS_ff_rand$layer
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
Rand_fire <- Rand_fire[, c(2:ncol(Rand_fire))]
head(Rand_fire)

write.csv(Rand_fire, './00_Data/Fire_data/Outputs/Random_points_data/Fire_frequency_random_points_QPWS_pres_IBRA.csv', row.names = F)


# 8.4 Convert the data to long format ----
Random_fire <- cbind(stack(Rand_fire[1:2]), Rand_fire[3:4])
tail(Random_fire)
colnames(Random_fire) <- c("Frequency", "Dataset", "Lon", "Lat")
head(Random_fire); tail(Random_fire); dim(Random_fire)
write.csv(Random_fire, './00_Data/Fire_data/Outputs/Random_points_data/Fire_frequency_randompoints_longformat_QPWS_presences_IBRA.csv', row.names = F)


Rand_fire <- read.csv('./00_Data/Fire_data/Outputs/Random_points_data/Fire_frequency_random_points_QPWS_pres.csv', header = T)
head(Rand_fire); dim(Rand_fire)
Random_fire <- read.csv('./00_Data/Fire_data/Outputs/Random_points_data/Fire_frequency_randompoints_longformat_QPWS_presences.csv', header = T)
head(Rand_fire); dim(Rand_fire)


# Explore the data to test for correlations ----
Rand <- cor.test(Rand_fire$Sentinel_rand_firefreq, Rand_fire$QPWS_rand_firefreq)
Rand # Vary in the same direction but weak correlation between these two variables. 

# As these two datasets are not very well correlated can we predict from QPWS for areas outside of QPWS estates?


head(Random_fire)
plot(Random_fire$Frequency[Random_fire$Dataset == "Sentinel_rand_firefreq"], Random_fire$Frequency[Random_fire$Dataset == "QPWS_rand_firefreq"], pch = 20, xlab = "Sentinel fire frequency", ylab = "QPWS fire frequency", ylim = c(0,15))

m1 <- lm(QPWS_rand_firefreq ~ Sentinel_rand_firefreq, data = Rand_fire)
summary(m1)
head(Rand_fire)

# Our problem now is that we are working on public and private land but the fire data available for private land is not adequate, it tends to underestimate in areas we know that have high fire frequency in the QPWS estates. Can we improve this by adding in environmental data?
# Landsat fire scar imagery records the earliest month in which a fire was detected where there are multiple fire scars for a given pixel. There is possibly some inaccuracy in the fire frequency calculation given this as we cannot determine if there was a fire that burnt a given pixel more than once in a year.