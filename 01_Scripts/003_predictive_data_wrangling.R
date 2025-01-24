# Written by Felicity Charles
# Date:1/08/2023

##### Fire frequency analysis ----
# This script begins to produce the data needed for producing the predictive model for fire frequency

# R version 4.3.1

# 1. Load required packages ----
library(terra) # terra_1.7-78 
library(dplyr) # dplyr_1.1.4
library(raster) # raster_3.6-23
library(tidyverse) # tidyverse_2.0.0
library(sf) # sf_1.0-14 
library(tmap) # tmap_3.3-3

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
QPWS_SEQ_ff <- rast('./00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_freq_hydrographical_mask.tif')
unique(QPWS_SEQ_ff$QPWS_SEQ_freq_raster)
plet(QPWS_SEQ_ff)
QPWS_rand <- vect('./00_Data/Fire_data/Outputs/QPWS_random.gpkg')

Sentinel_ff <- rast('./00_Data/Fire_data/Outputs/Sentinel/Sentinel_ff_hydrographical_mask_SEQ_focal_cropped.tif')
Sentinel_ff <- round(Sentinel_ff$focal_mean)
plot(Sentinel_ff)
unique(Sentinel_ff$focal_mean)


TWI <- rast('./00_Data/Environmental_data/Outputs/TWI/SEQ_TWI_cropped_focal.tif')
tempseason <- rast('./00_Data/Environmental_data/Outputs/BioClim/Tempseason_SEQ_cropped_focal.tif')
precipseason <- rast('./00_Data/Environmental_data/Outputs/BioClim/precipseason_SEQ_cropped_focal.tif')
diurnal_temp <- rast('./00_Data/Environmental_data/Outputs/BioClim/Diurnal_temp_meanSEQ_cropped_focal.tif')
solar_radiation <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/Solar_radiation_seq_cropped_focal.tif')
FPC <- rast('./00_Data/Environmental_data/Outputs/FPC/FPC_all_cropped_focal.tif')
soil_clay <- rast('./00_Data/Environmental_data/Outputs/Soil_clay/SEQ_soilclay_cropped_focal.tif')
slope <- rast('./00_Data/Environmental_data/Outputs/DEM/SEQ_slope_cropped_focal.tif')
aspect <- rast('./00_Data/Environmental_data/Outputs/DEM/SEQaspect_cropped_focal.tif')
topo_position <- rast('00_Data/Environmental_data/Outputs/DEM/SEQ_TPI_cropped_focal.tif')
elev <- rast('./00_Data/Environmental_data/Outputs/DEM/SEQ_DEM_reproj_cropped_focal.tif')


# 2.2 Read in the 10,000 random points to act as presences ----
Rand_fire <- read.csv('./00_Data/Fire_data/Outputs/Random_points_data/Fire_frequency_random_points_QPWS_pres.csv', header = T)
head(Rand_fire); dim(Rand_fire)
Rand_fire <- Rand_fire[, c(2:6)]

# 3. Extract environmental data for random points ----
TWI_rand <- terra::extract(TWI, QPWS_rand)
colnames(TWI_rand) <- c("ID", "TWI")


tempseason_rand <- terra::extract(tempseason, QPWS_rand) # Getting some NA values
colnames(tempseason_rand) <- c("ID", "tempseason")

precipseason_rand <- terra::extract(precipseason, QPWS_rand)# Getting some NA values
colnames(precipseason_rand) <- c("ID", "precipseason")


diurnal_temprand <- terra::extract(diurnal_temp, QPWS_rand)# Getting some NA values
colnames(diurnal_temprand) <- c("ID", "diurnal_temp")


solar_rad_rand <- terra::extract(solar_radiation, QPWS_rand)
colnames(solar_rad_rand) <- c("ID", "solar_radiation")


FPC_rand <- terra::extract(FPC, QPWS_rand) 
colnames(FPC_rand) <- c("ID", "FPC")


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

# 4.1 Add the environmental data to the fire dataframe ----
Rand_fire$TWI <- TWI_rand$TWI
Rand_fire$tempseason <- tempseason_rand$tempseason
Rand_fire$precipseason <- precipseason_rand$precipseason
Rand_fire$diurnal_temp <- diurnal_temprand$diurnal_temp
Rand_fire$solar_radiation <- solar_rad_rand$solar_radiation
Rand_fire$FPC <- FPC_rand$FPC
Rand_fire$soil_clay <- soil_clay_rand$soil_clay
Rand_fire$slope <- slope_rand$slope
Rand_fire$aspect <- aspect_rand$aspect
Rand_fire$topo_position <- topo_rand$topo_position
Rand_fire$elevation <- elev_rand$elevation
Rand_fire$Lon <- Rand_fire$Lon
Rand_fire$Lat <- Rand_fire$Lat
head(Rand_fire); dim(Rand_fire)
unique(Rand_fire$QPWS_rand_firefreq) # Check this looks right
Rand_fire$Sentinel_rand_firefreq <- round(Rand_fire$Sentinel_rand_firefreq)
unique(Rand_fire$Sentinel_rand_firefreq)
unique(is.na(Rand_fire)) # We still have some NA values which was can replace using nearestLand.



write.csv(Rand_fire, './00_Data/Fire_data/Outputs/Random_points_data/Fire_frequency_random_environmental_pres.csv')

View(Rand_fire) 





# Lets have a look at where the NA values may be falling
Rand_fire <- read.csv('./00_Data/Fire_data/Outputs/Random_points_data/Fire_frequency_random_environmental_pres.csv', header = T)
Rand_fire <- Rand_fire[, c(2:17)]
head(Rand_fire)
pt_precip <- Rand_fire[is.na(Rand_fire$precipseason), c(4:5)]
pt_precip_sf <- st_as_sf(pt_precip, coords = c(1:2), crs = 'EPSG:3577')
precipseason <- rast('./00_Data/Environmental_data/Outputs/BioClim/precipseason_SEQ_cropped_focal.tif')


tmap_mode("view")

# NOTE the mapping does not seem to be working entirely correctly as some 'NA' values appear to be in cells with a value. This gives the general idea of what we are looking at though. 
tm_shape(precipseason)+
  tm_raster(style = "cont")+
  tm_shape(pt_precip_sf)+
  tm_dots()





# NOTE the following was only run for the resampled data
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
twi <- raster('./00_Data/Environmental_data/Outputs/TWI/SEQ_TWI_cropped_focal.tif')
temp <- raster('./00_Data/Environmental_data/Outputs/BioClim/Tempseason_SEQ_cropped_focal.tif')
precip <- raster('./00_Data/Environmental_data/Outputs/BioClim/precipseason_SEQ_cropped_focal.tif')
diurnal <- raster('./00_Data/Environmental_data/Outputs/BioClim/Diurnal_temp_meanSEQ_cropped_focal.tif')
solar <- raster('./00_Data/Environmental_data/Outputs/Solar_radiation/Solar_radiation_seq_cropped_focal.tif')
fpc <- raster('./00_Data/Environmental_data/Outputs/FPC/FPC_all_cropped_focal.tif')
soil <- raster('./00_Data/Environmental_data/Outputs/Soil_clay/SEQ_soilclay_cropped_focal.tif')
slp <- raster('./00_Data/Environmental_data/Outputs/DEM/SEQ_slope_cropped_focal.tif')
asp <- raster('./00_Data/Environmental_data/Outputs/DEM/SEQaspect_cropped_focal.tif')
topo <- raster('00_Data/Environmental_data/Outputs/DEM/SEQ_TPI_cropped_focal.tif')
eleva <- raster('./00_Data/Environmental_data/Outputs/DEM/SEQ_DEM_reproj_cropped_focal.tif')


# 4.2.2 Extract the rows for each environmental data raster that have NA values -----
pt_temp <- Rand_fire[is.na(Rand_fire$tempseason), c(4:5)]
colnames(pt_temp) <- c('x', 'y')
pt_precip <- Rand_fire[is.na(Rand_fire$precipseason), c(4:5)]
colnames(pt_precip) <- c('x', 'y')
pt_diurnal <- Rand_fire[is.na(Rand_fire$diurnal_temp), c(4:5)]
colnames(pt_diurnal) <- c('x', 'y')
pt_solar <- Rand_fire[is.na(Rand_fire$solar_radiation), c(4:5)]
colnames(pt_solar) <- c('x', 'y')
pt_FPC <- Rand_fire[is.na(Rand_fire$FPC), c(4:5)]
colnames(pt_FPC) <- c('x', 'y')
pt_soil <- Rand_fire[is.na(Rand_fire$soil_clay), c(4:5)]
colnames(pt_soil) <- c('x', 'y')



# 4.2.3 Find nearest raster cell with non-NA value coordinates ----
nearest.temp <- nearestLand(pt_temp, temp, 558)
nearest.precip <- nearestLand(pt_precip, precip, 558)
nearest.diurnal <- nearestLand(pt_diurnal, diurnal, 558)
nearest.solar <- nearestLand(pt_solar, solar, 1857) # This is quite far as there are a couple of islands off the coast that have no solar radiation data, the nearest raster value is quite a distance from the furthest point
nearest.fpc <- nearestLand(pt_FPC, fpc, 88)
nearest.soil <- nearestLand(pt_soil, soil, 996)

# 4.2.4 Extract the values for these data points ----
temp.na <- terra::extract(tempseason, nearest.temp)
precip.na <- terra::extract(precipseason, nearest.precip)
diurnal.na <- terra::extract(diurnal_temp, nearest.diurnal)
solar.na <- terra::extract(solar_radiation, nearest.solar)
fpc.na <- terra::extract(FPC, nearest.fpc)
soil.na <- terra::extract(soil_clay, nearest.soil)

# 4.2.5 Combine the new values for NAs with the original data frame ----
# Firstly we need to get the point coordinates so this can be used to bind the data back in to the correct values
temp.na <- cbind(temp.na, nearest.temp)
precip.na <- cbind(precip.na, nearest.precip)
diurnal.na <- cbind(diurnal.na, nearest.diurnal)
solar.na <- cbind(solar.na, nearest.solar)
fpc.na <- cbind(fpc.na, nearest.fpc)
soil.na <- cbind(soil.na, nearest.soil)

# Replace any NA values with the nearest raster value for that point
Rand_fire$tempseason <- ifelse(is.na(Rand_fire$tempseason), temp.na$temperature_seasonality, Rand_fire$tempseason)
Rand_fire$precipseason <- ifelse(is.na(Rand_fire$precipseason), precip.na$precipitation_seasonality, Rand_fire$precipseason)
Rand_fire$diurnal_temp <- ifelse(is.na(Rand_fire$diurnal_temp), diurnal.na$diurnal_temp_seasonality, Rand_fire$diurnal_temp)
Rand_fire$solar_radiation <- ifelse(is.na(Rand_fire$solar_radiation), solar.na$Avg_solar_radiation, Rand_fire$solar_radiation)
Rand_fire$FPC <- ifelse(is.na(Rand_fire$FPC), fpc.na$Foliage_proj_cover, Rand_fire$FPC)
Rand_fire$soil_clay <- ifelse(is.na(Rand_fire$soil_clay), soil.na$percent_soil_clay, Rand_fire$soil_clay)
unique(is.na(Rand_fire))
head(Rand_fire); dim(Rand_fire)

# Save the output
write.csv(Rand_fire, './00_Data/Fire_data/Outputs/Random_points_data/Fire_frequency_random_environmental_pres_resampled.csv')





# 5. Create background points or psuedo absences ----
# Valavi et al 2021 paper - Predictive performance of presence-only species distribution models: a benchmark study with reproducible code - recent statistics papers suggest sampling background points irrespective of the presence points. They also suggest that there should be many more background points than presence points. So we will skip the following code and increase the number of background points we produce

# 5.1 Generate the random points within QPWS estates that are land areas
# Read in the data needed for cropping and erasing
SEQ <- vect('./00_Data/SEQ_bound/SEQ.gpkg')
Aus <- vect('./00_Data/Australia_shapefile/STE11aAust.shp') %>% 
  project("EPSG:3577")
canal <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Canal_SEQ.gpkg')
lake <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Lakes_SEQ.gpkg')
pond <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Ponds_SEQ.gpkg')
reservoir <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Reservoirs_SEQ.gpkg')
watercourse <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Watercourses_SEQ.gpkg')

# Crop and erase protected areas so that it only includes our area of interest and areas of land
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

# Create the background points
set.seed(183)
bg_rand <- spatSample(protected_land, size = 80000, method = "random") 
plet(bg_rand)
head(bg_rand)
length(bg_rand)
bg_rand$ID <- 1:nrow(bg_rand)
coords <- crds(bg_rand)
bg_rand$x <- coords[,1]
bg_rand$y <- coords[,2]
head(bg_rand)
bg_rand <- bg_rand[, c('ID', 'x', 'y')]
writeVector(bg_rand, './00_Data/Fire_data/Outputs/Background_points_data/bg_rand.gpkg', overwrite = T)


bg_rand <- vect('./00_Data/Fire_data/Outputs/Background_points_data/bg_rand.gpkg')
head(bg_rand); dim(bg_rand)

# 6. Add the fire and environmental data to the background points ----
# 6.1 Fire data ----
Sentinel_bg_ff <- terra::extract(Sentinel_ff, bg_rand)
unique(is.na(Sentinel_bg_ff)) # There are some NA values, we need to look at these more closely to determine which values have NA
View(Sentinel_bg_ff) # There are only a few points that have NA values
Sentinel_bg_ff <- cbind(Sentinel_bg_ff, bg_rand) # Add the point coordinates to the dataframe

# NOTE, the following was only run for the resampled data
# Replace NA values using the nearestLand function
pt_sent <- Sentinel_bg_ff[is.na(Sentinel_bg_ff$focal_mean), c(4,5)]
head(pt_sent);dim(pt_sent)
pts <- vect(pt_sent, geom = c('x', 'y'), crs = 'EPSG:3577') # This point should definitely have a value, we did have some random diagonal lines across with no data that this may have fallen within so we will find the closest value to assign to this point
plet(pts)

sent <- raster('./00_Data/Fire_data/Outputs/Sentinel/Sentinel_ff_hydrographical_mask_SEQ_focal_cropped.tif')


nearest.sent <- nearestLand(pt_sent, sent, 22) # Find the closest non-NA values
nearest.sent # Check this has worked
sent.na <- terra::extract(sent, nearest.sent) 
sent.na
sent.na <- cbind(sent.na, nearest.sent) # Bring the coordinates and new values together.
Sentinel_bg_ff$focal_mean <- ifelse(is.na(Sentinel_bg_ff$focal_mean), sent.na[,1], Sentinel_bg_ff$focal_mean)
unique(is.na(Sentinel_bg_ff)) # There are no more NA values
Sentinel_bg_ff <- Sentinel_bg_ff[c(1,2)]



QPWS_bg_ff <- terra::extract(QPWS_SEQ_ff, bg_rand)
unique(QPWS_bg_ff$QPWS_SEQ_freq_raster) # We expect to get NA values as there are locations with no fire history (e.g., any area outside QPWS estates), we will replace all NA values with 0.
QPWS_bg_ff[is.na(QPWS_bg_ff$QPWS_SEQ_freq_raster),] <- 0
unique(QPWS_bg_ff$QPWS_SEQ_freq_raster) # Check this has worked


# 6.2 Environmental data ----

TWI_rand <- terra::extract(TWI, bg_rand)
colnames(TWI_rand) <- c("ID", "TWI")



tempseason_rand <- terra::extract(tempseason, bg_rand)
colnames(tempseason_rand) <- c("ID", "tempseason")


precipseason_rand <- terra::extract(precipseason, bg_rand)
colnames(precipseason_rand) <- c("ID", "precipseason")


diurnal_tempran <- terra::extract(diurnal_temp, bg_rand)
colnames(diurnal_tempran) <- c("ID", "diurnal_temp")


solar_radran <- terra::extract(solar_radiation, bg_rand)
colnames(solar_radran) <- c("ID", "solar_radiation")


FPC_rand <- terra::extract(FPC, bg_rand) 
colnames(FPC_rand) <- c("ID", "FPC")


soil_rand <- terra::extract(soil_clay, bg_rand)
colnames(soil_rand) <- c("ID", "soil_clay")


sloperan <- terra::extract(slope, bg_rand)
colnames(sloperan) <- c("ID", "slope")


aspectran <- terra::extract(aspect, bg_rand)
colnames(aspectran) <- c("ID", "aspect")


topo_ran <- terra::extract(topo_position, bg_rand)
colnames(topo_ran) <- c("ID", "topo_position")

elev_ran <- terra::extract(elev, bg_rand)
colnames(elev_ran) <- c('ID', 'elevation')



# 6.3 Create a dataset with the background points fire and environmental data ----
Background_data <- QPWS_bg_ff
colnames(Background_data) <- c('ID', 'QPWS_rand_firefreq')
Background_data$Sentinel_rand_firefreq <- Sentinel_bg_ff$focal_mean
Background_data$Lon <- coords[,1]
Background_data$Lat <- coords[,2]
Background_data$TWI <- TWI_rand$TWI
Background_data$tempseason <- tempseason_rand$tempseason
Background_data$precipseason <- precipseason_rand$precipseason
Background_data$diurnal_temp <- diurnal_tempran$diurnal_temp
Background_data$solar_radiation <- solar_radran$solar_radiation
Background_data$FPC <- FPC_rand$FPC
Background_data$soil_clay <- soil_rand$soil_clay
Background_data$slope <- sloperan$slope
Background_data$aspect <- aspectran$aspect
Background_data$topo_position <- topo_ran$topo_position
Background_data$elevation <- elev_ran$elevation

# Check how this looks
head(Background_data)
tail(Background_data)
str(Background_data)
unique(is.na(Background_data)) # We do have NA values for all the environmental data so will need to fix these
unique(Background_data$Sentinel_rand_firefreq)
write.csv(Background_data, './00_Data/Fire_data/Outputs/Background_points_data/Fire_frequency_background_environmental_data.csv')


# 6.4 Replace NA values ----
# 6.4.1 Extract the rows for each environmental data raster that have NA values -----
pt_twi <- Background_data[is.na(Background_data$TWI), c(4:5)]
colnames(pt_twi) <- c('x', 'y')
pt_temp <- Background_data[is.na(Background_data$tempseason), c(4:5)]
colnames(pt_temp) <- c('x', 'y')
pt_precip <- Background_data[is.na(Background_data$precipseason), c(4:5)]
colnames(pt_precip) <- c('x', 'y')
pt_diurnal <- Background_data[is.na(Background_data$diurnal_temp), c(4:5)]
colnames(pt_diurnal) <- c('x', 'y')
pt_solar <- Background_data[is.na(Background_data$solar_radiation), c(4:5)]
colnames(pt_solar) <- c('x', 'y')
pt_FPC <- Background_data[is.na(Background_data$FPC), c(4:5)]
colnames(pt_FPC) <- c('x', 'y')
pt_soil <- Background_data[is.na(Background_data$soil_clay), c(4:5)]
colnames(pt_soil) <- c('x', 'y')
pt_slope <- Background_data[is.na(Background_data$slope), c(4:5)]
colnames(pt_slope) <- c('x', 'y')
pt_aspect <- Background_data[is.na(Background_data$aspect), c(4:5)]
colnames(pt_aspect) <- c('x', 'y')
pt_topo <- Background_data[is.na(Background_data$topo_position), c(4:5)]
colnames(pt_topo) <- c('x', 'y')
pt_elev <- Background_data[is.na(Background_data$elevation), c(4:5)]

# 6.4.2 Find nearest raster cell with non-NA value coordinates ----
nearest.twi <- nearestLand(pt_twi, twi, 22)
unique(is.na(nearest.twi))
nearest.temp <- nearestLand(pt_temp, temp, 830)
unique(is.na(nearest.temp))
nearest.precip <- nearestLand(pt_precip, precip, 830)
unique(is.na(nearest.precip))
nearest.diurnal <- nearestLand(pt_diurnal, diurnal, 830)
unique(is.na(nearest.diurnal))
nearest.solar <- nearestLand(pt_solar, solar, 9000) # This is quite far as there are a couple of islands off the coast that have no solar radiation data, the nearest raster value is quite a distance from the furthest point
unique(is.na(nearest.solar))
nearest.fpc <- nearestLand(pt_FPC, fpc, 400)
unique(is.na(nearest.fpc))
nearest.soil <- nearestLand(pt_soil, soil, 5000)
unique(is.na(nearest.soil))
nearest.slope <- nearestLand(pt_slope, slp, 25)
unique(is.na(nearest.slope))
nearest.aspect <- nearestLand(pt_aspect, asp, 25)
unique(is.na(nearest.aspect))
nearest.topo <- nearestLand(pt_topo, topo, 45)
unique(is.na(nearest.topo))
nearest.elev <- nearestLand(pt_elev, eleva, 25)

# 6.4.3 Extract the values for these data points ----
twi.na <- terra::extract(TWI, nearest.twi)
temp.na <- terra::extract(tempseason, nearest.temp)
precip.na <- terra::extract(precipseason, nearest.precip)
diurnal.na <- terra::extract(diurnal_temp, nearest.diurnal)
solar.na <- terra::extract(solar_radiation, nearest.solar)
fpc.na <- terra::extract(FPC, nearest.fpc)
soil.na <- terra::extract(soil_clay, nearest.soil)
slope.na <- terra::extract(slope, nearest.slope)
aspect.na <- terra::extract(aspect, nearest.aspect)
topo.na <- terra::extract(topo_position, nearest.topo)
elev.na <- terra::extract(elev, nearest.elev)

# 6.4.5 Combine the new values for NAs with the original data frame ----
# Firstly we need to get the point coordinates so this can be used to bind the data back in to the correct values
twi.na <- cbind(twi.na, nearest.twi)
colnames(twi.na) <- c('TWI', 'Lon','Lat')
temp.na <- cbind(temp.na, nearest.temp)
precip.na <- cbind(precip.na, nearest.precip)
diurnal.na <- cbind(diurnal.na, nearest.diurnal)
solar.na <- cbind(solar.na, nearest.solar)
fpc.na <- cbind(fpc.na, nearest.fpc)
soil.na <- cbind(soil.na, nearest.soil)
slope.na <- cbind(slope.na, nearest.slope)
aspect.na <- cbind(aspect.na, nearest.aspect)
topo.na <- cbind(topo.na, nearest.topo)
elev.na <- cbind(elev.na, nearest.elev)

# Replace any NA values with the nearest raster value for that point
Background_data$TWI <- ifelse(is.na(Background_data$TWI), twi.na$TWI, Background_data$TWI)
Background_data$tempseason <- ifelse(is.na(Background_data$tempseason), temp.na$temperature_seasonality, Background_data$tempseason)
Background_data$precipseason <- ifelse(is.na(Background_data$precipseason), precip.na$precipitation_seasonality, Background_data$precipseason)
Background_data$diurnal_temp <- ifelse(is.na(Background_data$diurnal_temp), diurnal.na$diurnal_temp_seasonality, Background_data$diurnal_temp)
Background_data$solar_radiation <- ifelse(is.na(Background_data$solar_radiation), solar.na$Avg_solar_radiation, Background_data$solar_radiation)
Background_data$FPC <- ifelse(is.na(Background_data$FPC), fpc.na$Foliage_proj_cover, Background_data$FPC)
Background_data$soil_clay <- ifelse(is.na(Background_data$soil_clay), soil.na$percent_soil_clay, Background_data$soil_clay)
Background_data$slope <- ifelse(is.na(Background_data$slope), slope.na$slope, Background_data$slope)
Background_data$aspect <- ifelse(is.na(Background_data$aspect), aspect.na$aspect, Background_data$aspect)
Background_data$topo_position <- ifelse(is.na(Background_data$topo_position), topo.na$Topo_position_index, Background_data$topo_position)
Background_data$elevation <- ifelse(is.na(Background_data$elevation), elev.na$elevation, Background_data$elevation)
unique(is.na(Background_data)) # Make sure all NA values have been replaced.
head(Background_data)

## FINAL OUTPUT
# Save the output
write.csv(Background_data, './00_Data/Fire_data/Outputs/Background_points_data/Fire_frequency_background_environmental_data_resampled.csv')



