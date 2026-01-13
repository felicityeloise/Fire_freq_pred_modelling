# Written by Felicity Charles
# Date:1/08/2023
# Updated: 9/10/2025

##### Fire frequency analysis ----
# This script begins to produce the data needed for producing the predictive model for fire frequency

# R version 4.5.1

# 1. Load required packages ----
library(terra) # terra_1.7-78 # Updated to 1.8-60
library(dplyr) # dplyr_1.1.4
library(raster) # raster_3.6-23 # Updated to 3.6-32
library(tidyverse) # tidyverse_2.0.0
library(sf) # sf_1.0-14 # Updated t0 1.0-21
library(tmap) # tmap_3.3-3 # Updated to 4.1

# Other attached packages not called directly
# lubridate_1.9.2 
# forcats_1.0.0  
# stringr_1.5.1   
# purrr_1.0.2     
# readr_2.1.4     
# tidyr_1.3.1    
# tibble_3.2.1    
# ggplot2_3.5.1
# sp_2.0-0

# 2. Read in the data ----
SEQ <- vect('./00_Data/SEQ_bound/SEQ_IBRA.gpkg')
template <- rast(SEQ)
res(template) <- c(30,30)
ext(template) <- ext(SEQ)

QPWS_SEQ_ff <- rast('./00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_IBRA_freq_hydrographical_mask_cropped.tif')%>% 
  resample(template, method = 'bilinear')
names(QPWS_SEQ_ff) <- 'QPWS_SEQ_ff' 
plet(QPWS_SEQ_ff)
QPWS_rand <- vect('./00_Data/Fire_data/Outputs/QPWS_random_SEQ_IBRA.gpkg')

Sentinel_ff <- rast('./00_Data/Fire_data/Outputs/Sentinel/SEQ_IBRA/Sentinel_ff_hydrographical_mask_SEQ_IBRA_focal.tif')
plot(Sentinel_ff)
unique(Sentinel_ff$Sentinel_fire_freq)


TWI <- rast('./00_Data/Environmental_data/Outputs/TWI/SEQ_IBRA_TWI_cropped_focal_masked.tif')
tempseason <- rast('./00_Data/Environmental_data/Outputs/SILO_Temperature/Average_Temp_seasonality_SEQ_IBRA_reproj_cropped_focal_masked.tif')
precipseason <- rast('./00_Data/Environmental_data/Outputs/SILO_Rainfall/Average_precipseason_SEQ_IBRA_reproj_cropped_focal_masked.tif')
FPC <- rast('./00_Data/Environmental_data/Outputs/FPC/Average_FPC_SEQ_IBRA_cropped_focal_masked.tif')
NDVI <- rast('./00_Data/Environmental_data/Outputs/BOM_NDVI/Average_NDVI_cropped_focal_masked.tif')
soil_clay <- rast('./00_Data/Environmental_data/Outputs/Soil_clay/SEQ_IBRA_soilclay_cropped_focal_masked.tif')
slope <- rast('./00_Data/Environmental_data/Outputs/DEM/SEQ_IBRA_slope_cropped_focal_masked.tif')
aspect <- rast('./00_Data/Environmental_data/Outputs/DEM/SEQ_IBRA_aspect_cropped_focal_masked.tif')
topo_position <- rast('00_Data/Environmental_data/Outputs/DEM/SEQ_IBRA_TPI_cropped_focal_masked.tif')
elev <- rast('./00_Data/Environmental_data/Outputs/DEM/SEQ_IBRA_DEM_cropped_focal_masked.tif')
BVG <- rast('./00_Data/Environmental_data/Outputs/BVG/BVG_cropped_masked.tif')

# 2.2 Read in the 10,000 random points to act as presences ----
Rand_fire <- read.csv('./00_Data/Fire_data/Outputs/Random_points_data/Fire_frequency_random_points_QPWS_pres_IBRA.csv', header = T)
head(Rand_fire); dim(Rand_fire)

# 3. Extract environmental data for random points ----
TWI_rand <- terra::extract(TWI, QPWS_rand)
colnames(TWI_rand) <- c("ID", "TWI")


tempseason_rand <- terra::extract(tempseason, QPWS_rand) # Getting some NA values

precipseason_rand <- terra::extract(precipseason, QPWS_rand)# Getting some NA values


FPC_rand <- terra::extract(FPC, QPWS_rand) 

NDVI_rand <- terra::extract(NDVI, QPWS_rand)

soil_clay_rand <- terra::extract(soil_clay, QPWS_rand)
colnames(soil_clay_rand) <- c("ID", "soil_clay")


slope_rand <- terra::extract(slope, QPWS_rand)
colnames(slope_rand) <- c("ID", "slope")


aspect_rand <- terra::extract(aspect, QPWS_rand)
colnames(aspect_rand) <- c("ID", "aspect")


topo_rand <- terra::extract(topo_position, QPWS_rand)
colnames(topo_rand) <- c("ID", "topo_position")

elev_rand <- terra::extract(elev, QPWS_rand)
colnames(elev_rand) <- c("ID", "elevation")

BVG_rand <- terra::extract(BVG, QPWS_rand)
colnames(BVG_rand) <- c("ID", "BVG")


# 4.1 Add the environmental data to the fire dataframe ----
Rand_fire$TWI <- TWI_rand$TWI
Rand_fire$temp_season <- tempseason_rand$Avg_Temperature_seasonality
Rand_fire$precip_season <- precipseason_rand$Avg_Precipitation_seasonality
Rand_fire$FPC <- FPC_rand$Avg_FPC
Rand_fire$NDVI <- NDVI_rand$Avg_NDVI
Rand_fire$soil_clay <- soil_clay_rand$soil_clay
Rand_fire$slope <- slope_rand$slope
Rand_fire$aspect <- aspect_rand$aspect
Rand_fire$topo_position <- topo_rand$topo_position
Rand_fire$elevation <- elev_rand$elevation
Rand_fire$BVG <- round(BVG_rand$BVG)
head(Rand_fire); dim(Rand_fire)


unique(Rand_fire$QPWS_rand_firefreq) # Check this looks right
Rand_fire$Sentinel_rand_firefreq <- round(Rand_fire$Sentinel_rand_firefreq)
unique(Rand_fire$Sentinel_rand_firefreq)
unique(is.na(Rand_fire)) # We still have some NA values which was can replace using nearestLand.



write.csv(Rand_fire, './00_Data/Fire_data/Outputs/Random_points_data/Fire_frequency_random_environmental_pres_IBRA.csv', row.names = F)






Rand_fire <- read.csv('./00_Data/Fire_data/Outputs/Random_points_data/Fire_frequency_random_environmental_pres_IBRA.csv', header = T)
head(Rand_fire)


# NOTE the following was only run for the dataset saved as resampled data
# 4.2 Replace NA values for the environmental data ----
# The following function is from the SEEG-Oxford/seegSDM GitHub repository https://github.com/SEEG-Oxford/seegSDM/tree/master. While this is a package readily available for installation, this library was unable to be installed on the current version fo R with rgeos also having been deprecated at the end of 2023. The function was downloaded from GitHub on the 1st of May 2024.
nearestLand <- function (points, raster, max_distance) {
  nearest <- function (lis, raster) {
    neighbours <- matrix(lis[[1]], ncol = 2)
    point <- lis[[2]]
    land <- !is.na(neighbours[, 2])
    if (!any(land)) {
      return (c(NA, NA))
    } else{
      coords <- xyFromCell(raster, neighbours[land, 1])   
      if (nrow(coords) == 1) {
        return (coords[1, ])
      }
      dists <- sqrt((coords[, 1] - point[1]) ^ 2 +
                      (coords[, 2] - point[2]) ^ 2)
      return (coords[which.min(dists), ])
    }
  }
  neighbour_list <- raster::extract(raster, points,
                                    buffer = max_distance,
                                    cellnumbers = TRUE)
  neighbour_list <- lapply(1:nrow(points),
                           function(i) {
                             list(neighbours = neighbour_list[[i]],
                                  point = as.numeric(points[i, ]))
                           })
  return (t(sapply(neighbour_list, nearest, raster)))
}


# 4.2.1 Read in the environmental data ----- 
twi <- raster('./00_Data/Environmental_data/Outputs/TWI/SEQ_IBRA_TWI_cropped_focal_masked.tif')
temp <- raster('./00_Data/Environmental_data/Outputs/SILO_Temperature/Average_Temp_seasonality_SEQ_IBRA_reproj_cropped_focal_masked.tif')
precip <- raster('./00_Data/Environmental_data/Outputs/SILO_Rainfall/Average_precipseason_SEQ_IBRA_reproj_cropped_focal_masked.tif')
fpc <- raster('./00_Data/Environmental_data/Outputs/FPC/Average_FPC_SEQ_IBRA_cropped_focal_masked.tif')
ndvi <- raster('./00_Data/Environmental_data/Outputs/BOM_NDVI/Average_NDVI_cropped_focal_masked.tif')
soil <- raster('./00_Data/Environmental_data/Outputs/Soil_clay/SEQ_IBRA_soilclay_cropped_focal_masked.tif')
slp <- raster('./00_Data/Environmental_data/Outputs/DEM/SEQ_IBRA_slope_cropped_focal_masked.tif')
asp <- raster('./00_Data/Environmental_data/Outputs/DEM/SEQ_IBRA_aspect_cropped_focal_masked.tif')
topo <- raster('00_Data/Environmental_data/Outputs/DEM/SEQ_IBRA_TPI_cropped_focal_masked.tif')
eleva <- raster('./00_Data/Environmental_data/Outputs/DEM/SEQ_IBRA_DEM_cropped_focal_masked.tif')
bvg <- raster('./00_Data/Environmental_data/Outputs/BVG/BVG_cropped_masked.tif')


# 4.2.2 Extract the rows for each environmental data raster that have NA values ----
pt_twi <- Rand_fire[is.na(Rand_fire$TWI), c(3:4)]
colnames(pt_twi) <- c('x', 'y')
pt_soil <- Rand_fire[is.na(Rand_fire$soil_clay), c(3:4)]
colnames(pt_soil) <- c('x', 'y')
pt_slp <- Rand_fire[is.na(Rand_fire$slope), c(3:4)]
colnames(pt_slp) <- c('x', 'y')
pt_asp <- Rand_fire[is.na(Rand_fire$aspect), c(3:4)]
colnames(pt_asp) <- c('x', 'y')
pt_topo <- Rand_fire[is.na(Rand_fire$topo_position), c(3:4)]
colnames(pt_topo) <- c('x', 'y')
pt_elev <- Rand_fire[is.na(Rand_fire$elevation), c(3:4)]
colnames(pt_elev) <- c('x', 'y')
pt_bvg <- Rand_fire[is.na(Rand_fire$BVG), c(3:4)]
colnames(pt_bvg) <- c('x', 'y')
pt_fpc <- Rand_fire[is.na(Rand_fire$FPC), c(3:4)]
colnames(pt_fpc) <- c('x', 'y')
pt_ndvi <- Rand_fire[is.na(Rand_fire$NDVI), c(3:4)]
colnames(pt_ndvi) <- c('x', 'y')
pt_temp <- Rand_fire[is.na(Rand_fire$temp_season), c(3:4)]
colnames(pt_temp) <- c('x', 'y')
pt_precip <- Rand_fire[is.na(Rand_fire$precip_season), c(3:4)]
colnames(pt_precip) <- c('x', 'y')





# 4.2.3 Find nearest raster cell with non-NA value coordinates ----
# Specify a sufficiently large distance to search, the function uses the closest non-NA cells
nearest.twi <- nearestLand(pt_twi, twi, 558)
nearest.soil <- nearestLand(pt_soil, soil, 6000)
nearest.slp <- nearestLand(pt_slp, slp, 558)
nearest.asp <- nearestLand(pt_asp, asp, 558)
nearest.topo <- nearestLand(pt_topo, topo, 558)
nearest.elev <- nearestLand(pt_elev, eleva, 558)
nearest.bvg <- nearestLand(pt_bvg, bvg, 200)
nearest.fpc <- nearestLand(pt_fpc, fpc, 200)
nearest.ndvi <- nearestLand(pt_ndvi, ndvi, 5000)
nearest.temp <- nearestLand(pt_temp, temp, 7000)
nearest.precip <- nearestLand(pt_precip, precip, 7000)


# 4.2.4 Extract the values for these data points ----
twi.na <- terra::extract(twi, nearest.twi) %>% 
  cbind(nearest.twi)
colnames(twi.na) <- c('twi', 'x', 'y')
soil.na <- terra::extract(soil_clay, nearest.soil) %>% 
  cbind(nearest.soil)
slp.na <- terra::extract(slope, nearest.slp) %>% 
  cbind(nearest.slp)
asp.na <- terra::extract(aspect, nearest.asp) %>% 
  cbind(nearest.asp)
tpi.na <- terra::extract(topo_position, nearest.topo) %>% 
  cbind(nearest.topo)
elev.na <- terra::extract(elev, nearest.elev) %>% 
  cbind(nearest.elev)
bvg.na <- terra::extract(bvg, nearest.bvg) %>% 
  cbind(nearest.bvg)
bvg_na <- as.data.frame(bvg.na)
colnames(bvg_na) <- c('BVG', 'x', 'y')
fpc.na <- terra::extract(FPC, nearest.fpc) %>% 
  cbind(nearest.fpc)
ndvi.na <- terra::extract(NDVI, nearest.ndvi) %>% 
  cbind(nearest.ndvi)
temp.na <- terra::extract(tempseason, nearest.temp) %>% 
  cbind(nearest.temp)
precip.na <- terra::extract(precipseason, nearest.precip) %>% 
  cbind(nearest.precip)

# 4.2.5 Combine the new values for NAs with the original data frame ----
# Replace any NA values with the nearest raster value for that point
Rand_fire$TWI <- ifelse(is.na(Rand_fire$TWI), twi.na[,1], Rand_fire$TWI)
Rand_fire$soil_clay <- ifelse(is.na(Rand_fire$soil_clay), soil.na$Percent_clay, Rand_fire$soil_clay)
Rand_fire$slope <- ifelse(is.na(Rand_fire$slope), slp.na$Slope, Rand_fire$slope)
Rand_fire$aspect <- ifelse(is.na(Rand_fire$aspect), asp.na$Aspect, Rand_fire$aspect)
Rand_fire$topo_position <- ifelse(is.na(Rand_fire$topo_position), tpi.na$Topo_position_index, Rand_fire$topo_position)
Rand_fire$elevation <- ifelse(is.na(Rand_fire$elevation), elev.na$Elevation, Rand_fire$elevation)
Rand_fire$BVG <- ifelse(is.na(Rand_fire$BVG), bvg_na$BVG, Rand_fire$BVG)
Rand_fire$FPC <- ifelse(is.na(Rand_fire$FPC), fpc.na$Avg_FPC, Rand_fire$FPC)
Rand_fire$NDVI <- ifelse(is.na(Rand_fire$NDVI), ndvi.na$Avg_NDVI, Rand_fire$NDVI)
Rand_fire$temp_season <- ifelse(is.na(Rand_fire$temp_season), temp.na$Avg_Temperature_seasonality, Rand_fire$temp_season)
Rand_fire$precip_season <- ifelse(is.na(Rand_fire$precip_season), precip.na$Avg_Precipitation_seasonality, Rand_fire$precip_season)
unique(is.na(Rand_fire))
head(Rand_fire); dim(Rand_fire)


# Save the output
write.csv(Rand_fire, './00_Data/Fire_data/Outputs/Random_points_data/Fire_frequency_random_environmental_pres_IBRA_resampled.csv', row.names = F)




# 5. Create background points or psuedo absences ----
# Valavi et al 2021 paper - Predictive performance of presence-only species distribution models: a benchmark study with reproducible code - recent statistics papers suggest sampling background points irrespective of the presence points. They also suggest that there should be many more background points than presence points. 

# 5.1 Generate the random points across the SEQ study region
# As our aim is to improve estimates of fire frequency for areas outside of public land, we will restrict model training to areas of public land where we have a better idea whether or not the land has burnt. So we restrict background point creation to protected land in SEQ


# Read in the data needed for cropping and erasing
canal <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Canal_SEQ_IBRA.gpkg')
lake <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Lakes_SEQ_IBRA.gpkg')
pond <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Ponds_SEQ_IBRA.gpkg')
reservoir <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Reservoirs_SEQ_IBRA.gpkg')
watercourse <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Watercourses_SEQ_IBRA.gpkg')

# Remove areas that are not land from the protected land polygon for spatial point creation
QLD <- vect('./00_Data/Australia_shapefile/STE11aAust.shp')
coast <- subset(QLD, QLD$STATE_NAME == "Queensland") %>%
  project('EPSG:3577')

SEQ2 <- mask(SEQ, coast) 
plet(SEQ2)

protected_land <- vect('./00_Data/Protected_areas/Protected_areas.shp') %>% 
  project('EPSG:3577') %>%
  crop(SEQ2) %>% 
  erase(canal) %>% 
  erase(lake) %>% 
  erase(pond) %>% 
  erase(reservoir) %>% 
  erase(watercourse)
plet(protected_land)

# Compare spatial extent to the other datasets
protected_land
tempseason
# The spatial extents slightly differ which could be problematic for sampling these points and replacing the non-NA values. So instead of using the protected_land polygons to create the samples, we will use this to mask one of the spatial rasters for point creation

bg_ext <- mask(tempseason, protected_land)
plet(bg_ext)
bg_ext


rm(list = setdiff(ls(), c("SEQ", "SEQ2", 'protected_land', "aspect", "elev", "TWI", "tempseason", "precipseason", "FPC", "soil_clay", "slope", "topo_position", "QPWS_SEQ_ff", "QPWS_rand", "Sentinel_ff", "Rand_fire", "nearestLand", "BVG", "NDVI", "bg_ext", 'twi', 'asp', 'bvg', 'eleva', 'fpc', 'ndvi', 'precip', 'slp', 'soil', 'temp', 'topo', 'twi'))); gc()


# Create the background points in protected areas of SEQ
# This step is quite computationally intensive, would recommend if encountering issues either splitting the region of interest into smaller blocks and creating a subset of the background points (e.g., splitting region into 8 blocks and creating 10,000 random background points within each block) or without splitting repeating spatSample process to create a small number of random points multiple times across the study region to more than the intended number (e.g., creating 95,000 points when needing only 80,000) then removing points which have duplicated coordinates.

set.seed(183)
bg_rand <- spatSample(bg_ext, size = 162383, method = "random", na.rm = T, as.points = T) 
plot(bg_rand); head(bg_rand); length(bg_rand)
bg_rand

bg_rand2 <- bg_rand
bg_rand2$ID <- 1:nrow(bg_rand2)
bg_rand2$x <- crds(bg_rand2)[,1]
bg_rand2$y <- crds(bg_rand2)[,2]
head(bg_rand2)
bg_rand2 <- bg_rand2[, 3:4]
head(bg_rand2)
bg_rand2

writeVector(bg_rand2, './00_Data/Fire_data/Outputs/Background_points_data/bg_rand_IBRA.gpkg', overwrite = T)


bg_rand <- vect('./00_Data/Fire_data/Outputs/Background_points_data/bg_rand_IBRA.gpkg')
head(bg_rand); dim(bg_rand); plet(bg_rand)



# 6. Add the fire and environmental data to the background points ----
# 6.1 Fire data ----
Sentinel_bg_ff <- terra::extract(Sentinel_ff, bg_rand)
unique(is.na(Sentinel_bg_ff)) # There are some NA values, we need to look at these more closely to determine which values have NA
Sentinel_bg_ff <- cbind(Sentinel_bg_ff, bg_rand) # Add the point coordinates to the dataframe

# Replace NA values using the nearestLand function
pt_sent <- Sentinel_bg_ff[is.na(Sentinel_bg_ff$Sentinel_fire_freq), c(4,5)]
head(pt_sent);dim(pt_sent)
pts <- vect(pt_sent, geom = c('x', 'y'), crs = 'EPSG:3577') # This point should definitely have a value, we did have some random diagonal lines across with no data that this may have fallen within so we will find the closest value to assign to this point
plet(pts)

sent <- raster('./00_Data/Fire_data/Outputs/Sentinel/Sentinel_ff_hydrographical_mask_SEQ_focal_cropped.tif')


nearest.sent <- nearestLand(pt_sent, sent, 22) # Find the closest non-NA values
nearest.sent # Check this has worked
sent.na <- terra::extract(sent, nearest.sent) 
sent.na
sent.na <- cbind(sent.na, nearest.sent) # Bring the coordinates and new values together.
Sentinel_bg_ff$Sentinel_fire_freq <- ifelse(is.na(Sentinel_bg_ff$Sentinel_fire_freq), sent.na[,1], Sentinel_bg_ff$Sentinel_fire_freq)
unique(is.na(Sentinel_bg_ff)) # There are no more NA values
Sentinel_bg_ff <- Sentinel_bg_ff[c(1,2)]


QPWS_bg_ff <- terra::extract(QPWS_SEQ_ff, bg_rand)
unique(QPWS_bg_ff$QPWS_SEQ_ff) # We expect to get NA values as there are locations with no fire history (e.g., any area outside QPWS estates), we will replace all NA values with 0.
QPWS_bg_ff[is.na(QPWS_bg_ff$QPWS_SEQ_ff),] <- 0
unique(QPWS_bg_ff$QPWS_SEQ_ff) # Check this has worked


# 6.2 Environmental data ----

TWI_ran <- terra::extract(TWI, bg_rand)
colnames(TWI_ran) <- c("ID", "TWI")



tempseason_ran <- terra::extract(tempseason, bg_rand)


precipseason_ran <- terra::extract(precipseason, bg_rand)


FPC_ran <- terra::extract(FPC, bg_rand) 

NDVI_ran <- terra::extract(NDVI, bg_rand)

soil_ran <- terra::extract(soil_clay, bg_rand)


sloperan <- terra::extract(slope, bg_rand)
colnames(sloperan) <- c("ID", "slope")


aspectran <- terra::extract(aspect, bg_rand)
colnames(aspectran) <- c("ID", "aspect")


topo_ran <- terra::extract(topo_position, bg_rand)
colnames(topo_ran) <- c("ID", "topo_position")

elev_ran <- terra::extract(elev, bg_rand)
colnames(elev_ran) <- c('ID', 'elevation')

BVG_ran <- terra::extract(BVG, bg_rand)
colnames(BVG_ran) <- c("ID", "BVG")


# 6.3 Create a dataset with the background points fire and environmental data ----
Background_data <- as.data.frame(1:80000)
Background_data$Sentinel_rand_firefreq <- Sentinel_bg_ff$Sentinel_fire_freq
Background_data$QPWS_rand_firefreq <- QPWS_bg_ff$QPWS_SEQ_ff
Background_data$Lon <- crds(bg_rand)[,1]
Background_data$Lat <- crds(bg_rand)[,2]
Background_data$TWI <- TWI_ran$TWI
Background_data$temp_season <- tempseason_ran$Avg_Temperature_seasonality
Background_data$precip_season <- precipseason_ran$Avg_Precipitation_seasonality
Background_data$FPC <- FPC_ran$Avg_FPC
Background_data$NDVI <- NDVI_ran$Avg_NDVI
Background_data$soil_clay <- soil_ran$Percent_clay
Background_data$slope <- sloperan$slope
Background_data$aspect <- aspectran$aspect
Background_data$topo_position <- topo_ran$topo_position
Background_data$elevation <- elev_ran$elevation
Background_data$BVG <- BVG_ran$BVG

# Check how this looks
head(Background_data); tail(Background_data)
str(Background_data)
unique(is.na(Background_data)) # We do have NA values for all the environmental data so will need to fix these
unique(Background_data$Sentinel_rand_firefreq)
Background_data$Sentinel_rand_firefreq <- round(Background_data$Sentinel_rand_firefreq)
Background_data <- Background_data[, c(2:ncol(Background_data))]
head(Background_data)

write.csv(Background_data, './00_Data/Fire_data/Outputs/Background_points_data/Fire_frequency_background_environmental_data_IBRA.csv', row.names = F)


# 6.4 Replace NA values ----
# 6.4.1 Extract the rows for each environmental data raster that have NA values -----
pt_twi <- Background_data[is.na(Background_data$TWI), c(3,4)]
colnames(pt_twi) <- c('x', 'y')


pt_soil <- Background_data[is.na(Background_data$soil_clay), c(3,4)]
colnames(pt_soil) <- c('x', 'y')

pt_slope <- Background_data[is.na(Background_data$slope), c(3,4)]
colnames(pt_slope) <- c('x', 'y')

pt_aspect <- Background_data[is.na(Background_data$aspect), c(3,4)]
colnames(pt_aspect) <- c('x', 'y')

pt_topo <- Background_data[is.na(Background_data$topo_position), c(3,4)]
colnames(pt_topo) <- c('x', 'y')

pt_elev <- Background_data[is.na(Background_data$elevation), c(3,4)]
colnames(pt_elev) <- c('x', 'y')

pt_bvg <- Background_data[is.na(Background_data$BVG), c(3,4)]
colnames(pt_bvg) <- c('x', 'y')

pt_temp_season <- Background_data[is.na(Background_data$temp_season), c(3,4)]
colnames(pt_temp_season) <- c('x', 'y')

pt_precip_season <- Background_data[is.na(Background_data$precip_season), c(3,4)]
colnames(pt_precip_season) <- c('x', 'y')

pt_fpc <- Background_data[is.na(Background_data$FPC), c(3,4)]
colnames(pt_fpc) <- c('x', 'y')

pt_ndvi <- Background_data[is.na(Background_data$NDVI), c(3,4)]
colnames(pt_ndvi) <- c('x', 'y')





# 6.4.2 Find nearest raster cell with non-NA value coordinates ----
nearest.twi <- nearestLand(pt_twi, twi, 55)
unique(is.na(nearest.twi))
nearest.soil <- nearestLand(pt_soil, soil, 7000) # Most are replaced at smaller distances but few NAs still remain
unique(is.na(nearest.soil))
nearest.slope <- nearestLand(pt_slope, slp, 40)
unique(is.na(nearest.slope))
nearest.aspect <- nearestLand(pt_aspect, asp, 40)
unique(is.na(nearest.aspect))
nearest.topo <- nearestLand(pt_topo, topo, 200)
unique(is.na(nearest.topo))
nearest.elev <- nearestLand(pt_elev, eleva, 40)
unique(is.na(nearest.elev))
nearest.bvg <- nearestLand(pt_bvg, bvg, 120)
unique(is.na(nearest.bvg))
nearest.temp <- nearestLand(pt_temp, temp, 7000)
unique(is.na(nearest.temp))
nearest.precip <- nearestLand(pt_precip, precip, 7000)
unique(is.na(nearest.precip))
nearest.ndvi <- nearestLand(pt_ndvi, ndvi, 5000)
unique(is.na(nearest.ndvi))
nearest.fpc <- nearestLand(pt_fpc, fpc, 200)
unique(is.na(nearest.fpc))

# 6.4.3 Extract the values for these data points ----
twi.na <- terra::extract(TWI, nearest.twi)
soil.na <- terra::extract(soil_clay, nearest.soil)
slope.na <- terra::extract(slope, nearest.slope)
aspect.na <- terra::extract(aspect, nearest.aspect)
topo.na <- terra::extract(topo_position, nearest.topo)
elev.na <- terra::extract(elev, nearest.elev)
bvg.na <- terra::extract(bvg, nearest.bvg)

# 6.4.5 Combine the new values for NAs with the original data frame ----
# Firstly we need to get the point coordinates so this can be used to bind the data back in to the correct values
twi.na <- cbind(twi.na, nearest.twi)
soil.na <- cbind(soil.na, nearest.soil)
slope.na <- cbind(slope.na, nearest.slope)
aspect.na <- cbind(aspect.na, nearest.aspect)
topo.na <- cbind(topo.na, nearest.topo)
elev.na <- cbind(elev.na, nearest.elev)
bvg.na <- cbind(bvg.na, nearest.bvg)
bvg_na <- as.data.frame(bvg.na)
colnames(bvg_na) <- c('BVG', 'x', 'y')

# Replace any NA values with the nearest raster value for that point
Background_data$TWI <- ifelse(is.na(Background_data$TWI), twi.na$Topo_wetness_index, Background_data$TWI)
Background_data$soil_clay <- ifelse(is.na(Background_data$soil_clay), soil.na$Percent_clay, Background_data$soil_clay)
Background_data$slope <- ifelse(is.na(Background_data$slope), slope.na$Slope, Background_data$slope)
Background_data$aspect <- ifelse(is.na(Background_data$aspect), aspect.na$Aspect, Background_data$aspect)
Background_data$topo_position <- ifelse(is.na(Background_data$topo_position), topo.na$Topo_position_index, Background_data$topo_position)
Background_data$elevation <- ifelse(is.na(Background_data$elevation), elev.na$Elevation, Background_data$elevation)
Background_data$BVG <- ifelse(is.na(Background_data$BVG), bvg_na$BVG, Background_data$BVG)
Background_data$temp_season <- ifelse(is.na(Background_data$temp_season), temp.na$Avg_Temperature_seasonality, Background_data$temp_season)
Background_data$precip_season <- ifelse(is.na(Background_data$precip_season), precip.na$Avg_Precipitation_seasonality, Background_data$precip_season)
Background_data$FPC <- ifelse(is.na(Background_data$FPC), fpc.na$Avg_FPC, Background_data$FPC)
Background_data$NDVI <- ifelse(is.na(Background_data$NDVI), ndvi.na$Avg_NDVI, Background_data$NDVI)

unique(is.na(Background_data)) # Make sure all NA values have been replaced.
head(Background_data)

## FINAL OUTPUT
# Save the output
write.csv(Background_data, './00_Data/Fire_data/Outputs/Background_points_data/Fire_frequency_background_environmental_data_resampled_IBRA.csv', row.names = F)



