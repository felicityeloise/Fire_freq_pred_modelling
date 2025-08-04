# Written by Felicity Charles
# Date:1/08/2023
# Update: 30/07/2025

##### Fire frequency analysis ----
# This script gathers together environmental data needed for predicting  fire frequency

# R version 4.3.1

# 1. Load required packages ----
library(terra) # terra_1.7-78 
library(dplyr) # dplyr_1.1.4 
library(landform) # landform_0.2


# Install gdalUtilities for working with large datasets
#library(devtools)
#install_github("JoshOBrien/gdalUtilities")
library(gdalUtilities) # gdalUtilities_1.2.5


# 1. Read in, crop and aggregate environmental data ----
# Create an SEQ vector 
SEQ <- vect('./00_Data/SEQ_bound/SEQ_IBRA.gpkg')

# Use SEQ to create a template raster for resampling 
template <- rast(SEQ)
res(template) <- c(30,30)
ext(template) <- ext(SEQ)


# Investigated alternative FPC data, TERN has same data I already am using. QSpatial also has SLATS but prior to 2017/2018 only mapping land cover change so not able to get yearly data for woody vegetation extent.
# Also investigated ANUClimate data from TERN for climate data but limited temporal coverage depending on the dataset not covering most recent years of the study period. BOM has gridded climatology data but this is also temporally limited, cost associated with requesting data. May be even coarser than 5km resolution



# 1.1 Rainfall ----
#### THIS IS WORSE SPATIAL RESOULTION THAN WORLDCLIM DATA AT 5km BUT BETTER TEMPORAL COVERAGE. NOT SURE IF THIS IS THE BEST WAY TO APPROACH THIS, GOING TO BE SLOW
# 1.1.1 Load daily rainfall data
load_SILO_rain <- function(year, dest_dir = './00_Data/Environmental_data/SILO_Rainfall'){
  if(!dir.exists(dest_dir)) dir.create(dest_dir)
  
  base_url <- "https://s3-ap-southeast-2.amazonaws.com/silo-open-data/Official/annual/daily_rain"
  file_name <- paste0(year, ".daily_rain.nc")
  url <- file.path(base_url, file_name)
  local_file <- file.path(dest_dir, file_name)
  
  # Download file
  if(!file.exists(local_file)){
    message("Downloading: ", year)
    download.file(url, local_file, mode = 'wb')
  }
  else{
    message("Using cached file: ", year)
  }
  
  # Load raster
  r <- rast(local_file) %>% 
    project('EPSG:3577') %>% 
    crop(SEQ)
  assign(paste0("rain_", year), r, envir = .GlobalEnv)
}

# Loop over 1987 to 2023
years <- 1987:2023
invisible(lapply(years, load_SILO_rain))
rain_1987; plot(rain_1987)



# IF using this data we can continue to do the following - will need to rescale but maybe leave the rescaling until after these values are calculated
# 1.1.2 Calculate rainfall seasonality
# Calculate mean monthly rainfall
rain_list <- list(rain_1987, rain_1988, rain_1989, rain_1990, rain_1991, rain_1992, rain_1993, rain_1994, rain_1995, rain_1996, rain_1997, rain_1998, rain_1999, rain_2000, rain_2001, rain_2002, rain_2003, rain_2004, rain_2005, rain_2006, rain_2007, rain_2008, rain_2009, rain_2010, rain_2011, rain_2012, rain_2013, rain_2014, rain_2015, rain_2016, rain_2017, rain_2018, rain_2019, rain_2020, rain_2021, rain_2022, rain_2023)
years <- 1987:2023


for(i in seq_along(rain_list)){
  rain <- rain_list[[i]]
  year <- years[i]
  dates <- time(rain) # Get time vector 
  month_groups <- format(dates, '%Y-%m') # Create grouping by month
  
  monthly_avg <- tapp(rain, index = month_groups, fun = mean, na.rm = T) # Calculate monthly average rainfall per cell
  assign(paste0('avg_monthly_rain', year), monthly_avg)
}
avg_monthly_rain1987 # Check

# Caculate the seasonality index
# The following equation is adapted from Walsh, R. and Lawler, D. (1981). Rainfall seasonality: Description, spatial patterns and change through time (British Isles, Africa). Weather, 36(7), 201-208. doi:10.1002/j.1477-8696.1981.tb05400.x. 
# This index can theoretically vary from 0 (when all months have the same      #
# rainfall) to 1.83 (when all the rainfall ocurrs in a single month)           #
# A qualitative classification of degrees of seasonality is the following:     #
# -----------------------------------------------------------------------------#
#  si values   |                    Rainfall regime                            #
# -----------------------------------------------------------------------------#
#     <= 0.19  | Very equable                                                  #
# 0.20 - 0.39  | Equable but with a definite wetter season                     #
# 0.40 - 0.59  | Rather seasonal with a short drier season                     #
# 0.60 - 0.79  | Seasonal                                                      #
# 0.80 - 0.99  | Markedly seasonal with a long drier season                    #
# 1.00 - 1.19  | Most rain in 3 months or less                                 #
#     >= 1.20  | Extreme, almost all rain in 1-2 months                        #
################################################################################
for(i in years){
  month_rast_name <- paste0("avg_monthly_rain", i)
  month_rast <- get(month_rast_name)
  
  R <- sum(month_rast) # Calculate mean annual rainfall from mean monthly rainfall
  
  si_rain <- (1/R) * sum(abs(month_rast - R/12))
  assign(paste0("rain_seasonality_", i), si_rain)
}
# Would also need to save the output
rain_seasonality_1987; rain_seasonality_2023 # Check these seem reasonable


# 1.2 Temperature ----
# 1.2.1 Min temperature
load_SILO_mintemp <- function(year, dest_dir = './00_Data/Environmental_data/SILO_Min_temp'){
  if(!dir.exists(dest_dir)) dir.create(dest_dir)
  
  base_url <- "https://s3-ap-southeast-2.amazonaws.com/silo-open-data/Official/annual/min_temp"
  file_name <- paste0(year, ".min_temp.nc")
  url <- file.path(base_url, file_name)
  local_file <- file.path(dest_dir, file_name)
  
  # Download file
  if(!file.exists(local_file)){
    message("Downloading: ", year)
    download.file(url, local_file, mode = 'wb')
  }
  else{
    message("Using cached file: ", year)
  }
  
  # Load raster
  r <- rast(local_file) %>% 
    project('EPSG:3577') %>% 
    crop(SEQ)
  assign(paste0("mintemp_", year), r, envir = .GlobalEnv)
}

# Loop over 1987 to 2023
years <- 1987:2023
invisible(lapply(years, load_SILO_mintemp))
mintemp_1987; plot(mintemp_1987) # Values are per day




# 1.2.1 Max temperature
load_SILO_maxtemp <- function(year, dest_dir = './00_Data/Environmental_data/SILO_Max_temp'){
  if(!dir.exists(dest_dir)) dir.create(dest_dir)
  
  base_url <- "https://s3-ap-southeast-2.amazonaws.com/silo-open-data/Official/annual/max_temp"
  file_name <- paste0(year, ".max_temp.nc")
  url <- file.path(base_url, file_name)
  local_file <- file.path(dest_dir, file_name)
  
  # Download file
  if(!file.exists(local_file)){
    message("Downloading: ", year)
    download.file(url, local_file, mode = 'wb')
  }
  else{
    message("Using cached file: ", year)
  }
  
  # Load raster
  r <- rast(local_file) %>% 
    project('EPSG:3577') %>% 
    crop(SEQ)
  assign(paste0("maxtemp_", year), r, envir = .GlobalEnv)
}

# Loop over 1987 to 2023
years <- 1987:2023
invisible(lapply(years, load_SILO_maxtemp))
maxtemp_1987; plot(maxtemp_1987) # Values are per day

# 1.2.2 Calculate mean daily temperature
for (i in years ){
  cat("Processing year: ", i, "\n")
  
  tmin <- get(paste0("mintemp_", i)) # Get minimum temperature raster for year i
  tmax <- get(paste0("maxtemp_", i)) # Get maximum temperature raster for year i
  
  tavg <- (tmin + tmax) / 2 # Calculate daily mean temperature
  
  dates <- time(tavg) 
  month_groups <- format(dates, "%Y-%m") # Create month groups
  
  monthly_avg <- tapp(tavg, index = month_groups, fun = mean, na.rm = T) # Calculate monthly average temperature
  assign(paste0('avg_monthly_temp', i), monthly_avg)
}

avg_monthly_temp1987; plot(avg_monthly_temp1987) # Check

# 1.2.3 Calculate temperature seasonality
for(i in years){
  month_rast_name <- paste0("avg_monthly_temp", i)
  month_rast <- get(month_rast_name)
  
  R <- sum(month_rast) # Calculate mean annual temperature from mean monthly temperature
  
  si_temp <- (1/R) * sum(abs(month_rast - R/12))
  assign(paste0("temp_seasonality_", i), si_temp)
}
# Would also need to save the output
temp_seasonality_1987; temp_seasonality_2023 # Check these seem reasonable
plot(temp_seasonality_1987)

#### For rescaling, should actually use bilinear interpolation (method = bilinear) as this is continuous data.

rm(list = setdiff(ls(), "SEQ", "template")); gc()


# 1.1 Topographic Wetness Index ----

TWI <- rast("/vsicurl/https://s3.data.csiro.au/dapprd/000005588v002/data/TopographicWetnessIndex_1_arcsecond_resolution/mosaic/twi_1s.tif?response-content-disposition=attachment%3B%20filename%3D%22twi_1s.tif%22&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20250731T110706Z&X-Amz-SignedHeaders=host&X-Amz-Expires=172800&X-Amz-Credential=GSLZCSPGJ4EIL6KMFTJR%2F20250731%2FCDC%2Fs3%2Faws4_request&X-Amz-Signature=723defc55e80019ed8c3e0ca15045c0bcf13a5b779fa2c6ad2ae336cf9fc249c") # To get this link need to click download on the file from https://data.csiro.au/collection/csiro:5588 and then as it is downloading open up the downloads drop down > right click on the file > copy download link. Need to get a new link each time this line is run otherwise it will not work.

writeRaster(TWI, "./00_Data/Environmental_data/Outputs/TWI/TWI.tif")

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/TWI/TWI.tif', # Source file
         dstfile = './00_Data/Environmental_data/Outputs/TWI/TWI_reproj.tif', # Destination file
         t_srs = 'EPSG:3577', # CRS to be transformed to  
         tr = c(30,30), # Resolution to be transformed to
         r = 'bilinear') # Use bilinear interpolation for continuous data 

TWI <- rast('./00_Data/Environmental_data/Outputs/TWI/TWI_reproj.tif')
TWI2 <- crop(TWI, SEQ)
writeRaster(TWI2, "./00_Data/Environmental_data/Outputs/TWI/SEQ_TWI.tif")


# 1.2 Climate data ----
# BIO2 - Mean Diurnal range BIO4 - temperature seasonality and BIO15  - precipitation seasonality

gdalwarp(srcfile = "./00_Data/Environmental_data/BioClim/wc2.1_30s_bio_2.tif",
         dstfile = "./00_Data/Environmental_data/Outputs/BioClim/Diurnal_temp_mean.tif",
         t_srs = 'EPSG:3577')


diurnal_temp_mean <- rast("./00_Data/Environmental_data/Outputs/BioClim/Diurnal_temp_mean.tif")
dtm <- crop(diurnal_temp_mean, SEQ)
dtm1 <- resample(dtm, template, method = 'bilinear') # More appropriate to use resample for precise resolution changes with continuous data. Disaggregation will just duplicate values which is better for categorical data. Aggregation is still appropriate for other times as we can control how values are upscaled (e.g., mean, max, or sum)
#dtm1 <- disagg(dtm, fact = 30)
plot(dtm1) # Check how this looks
dtm1 # Check crs and resolution

writeRaster(dtm1, './00_Data/Environmental_data/Outputs/BioClim/Diurnal_temp_meanSEQ.tif')


gdalwarp(srcfile = './00_Data/Environmental_data/BioClim/wc2.1_30s_bio_4.tif',
         dstfile = '00_Data/Environmental_data/Outputs/BioClim/Tempseason.tif',
         t_srs = 'EPSG:3577')


tempseason <- rast('./00_Data/Environmental_data/Outputs/BioClim/Tempseason.tif')
range(unique(tempseason$Tempseason)) 

temp <- crop(tempseason, SEQ)
plot(temp)
tempr <- resample(temp, template, method = 'bilinear')
#tempr <- disagg(temp, fact = 30)
plot(tempr) # Check to see if this looks right
tempr # We have the right coordinate reference system and resolution is nearly correct. 

writeRaster(tempr, '00_Data/Environmental_data/Outputs/BioClim/Tempseason_SEQ_IBRA.tif')


gdalwarp(srcfile = './00_Data/Environmental_data/BioClim/wc2.1_30s_bio_15.tif',
         dstfile = '00_Data/Environmental_data/Outputs/BioClim/precipseason.tif',
         t_srs = 'EPSG:3577')

# Need to crop and change resolution
precipseason <- rast('00_Data/Environmental_data/Outputs/BioClim/precipseason.tif')
precip <- crop(precipseason, SEQ)
precipr <- resample(precip, template, method = 'bilinear')
#precipr <- disagg(precip, fact = 30)
plot(precipr)
precipr

writeRaster(precipr, '00_Data/Environmental_data/Outputs/BioClim/precipseason_SEQ_IBRA.tif')




# 1.3 Foliage projective cover ----
FPC14 <- rast('./00_Data/Environmental_data/FPC/DP_QLD_FPC2014.tif')
unique(FPC14$DP_QLD_FPC2014)
# Need to do some more adjustments to this data - metadata states that data ranges between 100-200 which is equivalent to 0-100% FPC. values erroneously predicted above 100% or below 0% have been classed as above 200 and below 100 respectively. Zero values indicate NULL data. The data actually seems to be ranging between 88-213. Post 2014, values range between 0-100 which would denote the % cover without any further changes being required. Let's take a look at the data in ArcGIS as well to make sure this is true for the 2014 dataset.


# Create matrices for reclassification

A = matrix(
  c(88:99, 201:213),
  nrow = 25,
  ncol = 2)
A[,2] <- 0

B = matrix(
  c(100:200),
  nrow = 101,
  ncol = 1
)
B <- cbind(B, 0:100)


reclas <- rbind(A, B)

# Now reclassify FPC14
FPC14r <- classify(FPC14, rcl = reclas)
FPC14r # Check how this looks
plot(FPC14r)

FPC14_seq <- crop(FPC14r, SEQ)
FPC14seq <- project(FPC14_seq, 'EPSG:3577')
FPC14seq
plot(FPC14seq)

writeRaster(FPC14seq, './00_Data/Environmental_data/Outputs/FPC/FPC2014_SEQ.tif')


# 1.3.1 Add in the data from more recent years post 2014 ----

# Firstly, look at the new data to see what needs to be changed
FPC18 <- rast('./00_Data/Environmental_data/FPC/DP_QLD_S2_WOODY_FPC_2018.tif')

# Need to aggregate the data to a coarser resolution, from 10m to 30m, and then we also need to crop the data to SEQ

FPC18_seq <- crop(FPC18, SEQ)

FPC18seq <- terra::aggregate(FPC18_seq, fact = 3, fun = 'mean') # Calculate the average FPC for the cell when aggregating. 
FPC18seq # Check how this looks
plot(FPC18seq)

writeRaster(FPC18seq, './00_Data/Environmental_data/Outputs/FPC/FPC18_SEQ.tif')


FPC19 <- rast('./00_Data/Environmental_data/FPC/DP_QLD_S2_WOODY_FPC_2019.tif')
FPC19_seq <- crop(FPC19, SEQ)
FPC19seq <- terra::aggregate(FPC19_seq, fact = 3, fun = 'mean')
FPC19seq
plot(FPC19seq)
writeRaster(FPC19seq, './00_Data/Environmental_data/Outputs/FPC/FPC19_SEQ.tif')


FPC20 <- rast('./00_Data/Environmental_data/FPC/DP_QLD_S2_WOODY_FPC_2020.tif')
FPC20_seq <- crop(FPC20, SEQ)
FPC20seq <- terra::aggregate(FPC20_seq, fact = 3, fun = 'mean')
FPC20seq
plot(FPC20seq)
writeRaster(FPC20seq, './00_Data/Environmental_data/Outputs/FPC/FPC20_SEQ.tif')


FPC21 <- rast('./00_Data/Environmental_data/FPC/DP_QLD_S2_FPC_2021.tif') # The coordinate reference system has not been read in the same manner as the others so we will need to fix this

FPC21_seq <- crop(FPC21, SEQ)
FPC21seq <- project(FPC21_seq, 'EPSG:3577')
FPC21seq <- terra::aggregate(FPC21_seq, fact = 3, fun = 'mean')
FPC21seq
plot(FPC21seq)
writeRaster(FPC21seq, './00_Data/Environmental_data/Outputs/FPC/FPC21_SEQ.tif')


# 1.3.2 Combine the FPC data into one raster ----
FPC14seq <- resample(FPC14seq, FPC18seq, method = 'bilinear') # Need the extents to match
writeRaster(FPC14seq, './00_Data/Environmental_data/Outputs/FPC/FPC14_SEQ_resampled.tif')

#stack <- c(FPC14,FPC18,FPC19,FPC20,FPC21)
# This produces a raster with each year as a separate raster layer. But what we want, because FPC2014 is from 1988-2014 and the others are separate years, we need an average FPC value across the years.

#FPC <- terra::mean(FPC14seq, FPC18seq, FPC19seq, FPC20seq, FPC21seq)
#FPC
#plot(FPC)


writeRaster(FPC, './00_Data/Environmental_data/Outputs/FPC/FPC_all.tif')




# 1.4 Elevation data ----
# Download the DEM ad access straight from download folder

DEM <- rast('./00_Data/Environmental_data/69816/srtm-1sec-dem-v1-COG.tif') # Work from the original data


# Need to cut down to smaller area
Aus <- vect('./00_Data/Australia_shapefile/STE11aAust.shp')
QLD <- subset(Aus, Aus$STATE_NAME == "Queensland")
QLD <- project(QLD, 'EPSG:4326')
DEMqld <- crop(DEM, QLD)

SEQ <- project(SEQ, 'EPSG:4326')
DEMseq <- crop(DEM, SEQ)

writeRaster(DEMseq, './00_Data/Environmental_data/Outputs/DEM/SEQ_DEM.tif')
gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/DEM/SEQ_DEM.tif',
         dstfile = './00_Data/Environmental_data/Outputs/DEM/SEQ_DEM_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30),
         r = 'bilinear')
DEM <- rast('./00_Data/Environmental_data/Outputs/DEM/SEQ_DEM_reproj.tif')



# 1.4.1 Calculate slope, aspect and topographic position index (TPI) ----

slope <- terrain(DEMseq, v = "slope", unit = "degrees")
writeRaster(slope, './00_Data/Environmental_data/Outputs/DEM/slope.tif') # Save output


# Fix projection and change resolution so it can be used with the other data
gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/DEM/slope.tif',
         dstfile = './00_Data/Environmental_data/Outputs/DEM/slope_reproj.tif',
         tr = c(30,30),
         r = 'bilinear')

slope <- rast('./00_Data/Environmental_data/Outputs/DEM/slope_reproj.tif') %>% 
  project('EPSG:3577')
slope #Check
slope <- crop(slope, SEQ)
plot(slope)
slope

writeRaster(slope, './00_Data/Environmental_data/Outputs/DEM/SEQ_slope.tif')


aspect <- terrain(DEMseq, v = "aspect")
writeRaster(aspect, './00_Data/Environmental_data/Outputs/DEM/aspect.tif')

# Fix projection and change resolution so it can be used with the other data
gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/DEM/aspect.tif',
         dstfile = './00_Data/Environmental_data/Outputs/DEM/aspect_reproj.tif',
         tr = c(30,30),
         r = 'bilinear')

aspect <- rast('./00_Data/Environmental_data/Outputs/DEM/aspect_reproj.tif') %>% 
  project('EPSG:3577')
aspect <- crop(aspect, SEQ)
plot(aspect)
aspect

writeRaster(aspect, './00_Data/Environmental_data/Outputs/DEM/SEQaspect.tif')


TP1 <- landform(DEMqld, class.type = "slope.position")
TPI <- TP1$all

# Here TPI values are as follows:
# 1 = Valley
# 2 = Lower slope
# 3 = Flat slope
# 4 = Middle slope
# 5 = Upper slope
# 6 = Ridge

writeRaster(TPI, "./00_Data/Environmental_data/Outputs/DEM/TPI.tif")

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/DEM/TPI.tif',
         dstfile = './00_Data/Environmental_data/Outputs/DEM/TPI_reproj.tif',
         t_srs ='EPSG:3577',
         tr = c(30,30),
         r = 'near') #TPI is categorical numeric

TPI <- rast('./00_Data/Environmental_data/Outputs/DEM/TPI_reproj.tif')
TPI <- crop(TPI, SEQ)
plot(TPI)
TPI

writeRaster(TPI, './00_Data/Environmental_data/Outputs/DEM/SEQ_TPI.tif')





# 1.5 Soil nutrients using soil % clay ----
# Download data from https://esoil.io/TERNLandscapes/Public/Pages/SLGA/GetData-COGSDataStore_SLGA.html
# This variable is included as nutrients influence plant growth which would then influence the occurrence of fire

# 0 to 0.05m
clay0to0.05 <- rast("./00_Data/Environmental_data/Soil_clay/CLY_000_005_EV_N_P_AU_TRN_N_20210902.tif")
clay0to0.05 <- crop(clay0to0.05, QLD)
plot(clay0to0.05) # Check this worked correctly

writeRaster(clay0to0.05, "./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_0-05m.tif")

gdalwarp(srcfile = "./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_0-05m.tif",
         dstfile = "./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_0-05m_reproj.tif",
         t_srs = 'EPSG:3577',
         tr = c(30,30),
         r = 'bilinear')

clay0to0.05 <- rast("./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_0-05m_reproj.tif")
clay0to0.05 <- crop(clay0to0.05, SEQ)
plot(clay0to0.05) # Check this looks right

writeRaster(clay0to0.05, "./00_Data/Environmental_data/Outputs/Soil_clay/SEQ_Clay_0-05m.tif")



# 0.05 to 0.15m
clay2 <- rast("./00_Data/Environmental_data/Soil_clay/CLY_005_015_EV_N_P_AU_TRN_N_20210902.tif")
clay2 <- crop(clay2, QLD)
plot(clay2) # Check this worked correctly

writeRaster(clay2, "./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_0-015m.tif")

gdalwarp(srcfile = "./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_0-015m.tif",
         dstfile = "./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_0-015m_reproj.tif",
         t_srs = 'EPSG:3577',
         tr = c(30,30),
         r = 'bilinear')

clay2 <- rast("./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_0-015m_reproj.tif")
clay2 <- crop(clay2, SEQ)
plot(clay2) # Check this looks right

writeRaster(clay2, "./00_Data/Environmental_data/Outputs/Soil_clay/SEQ_Clay_0-015m.tif")



# 0.15 to 0.3m

clay3 <- rast("./00_Data/Environmental_data/Soil_clay/CLY_015_030_EV_N_P_AU_TRN_N_20210902.tif")
clay3 <- crop(clay3, QLD)
plot(clay3) # Check this worked correctly

writeRaster(clay3, "./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_015-03m.tif")

gdalwarp(srcfile = "./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_015-03.tif",
         dstfile = "./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_015-03_reproj.tif",
         t_srs = 'EPSG:3577',
         tr = c(30,30),
         r = 'bilinear')

clay3 <- rast("./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_015-03m_reproj.tif")
clay3 <- crop(clay3, SEQ)
plot(clay3) # Check this looks right

writeRaster(clay3, "./00_Data/Environmental_data/Outputs/Soil_clay/SEQ_Clay_015-03m.tif")




# 0.3 to 0.6m 

clay4 <- rast("./00_Data/Environmental_data/Soil_clay/CLY_030_060_EV_N_P_AU_TRN_N_20210902.tif")
clay4 <- crop(clay4, QLD)
plot(clay4) # Check this worked correctly

writeRaster(clay4, "./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_03-06m.tif")

gdalwarp(srcfile = "./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_03-06m.tif",
         dstfile = "./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_03-06m_reproj.tif",
         t_srs = 'EPSG:3577',
         tr = c(30,30),
         r = 'bilinear')

clay4 <- rast("./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_03-06m_reproj.tif")
clay4 <- crop(clay4, SEQ)
plot(clay4) # Check this looks right

writeRaster(clay4, "./00_Data/Environmental_data/Outputs/Soil_clay/SEQ_Clay_03-06m.tif")



# 0.6 to 1m

clay5 <- rast("./00_Data/Environmental_data/Soil_clay/CLY_060_100_EV_N_P_AU_TRN_N_20210902.tif")
clay5 <- crop(clay5, QLD)
plot(clay5) # Check this worked correctly

writeRaster(clay5, "./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_06-1m.tif")

gdalwarp(srcfile = "./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_06-1m.tif",
         dstfile = "./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_06-1m_reproj.tif",
         t_srs = 'EPSG:3577',
         tr = c(30,30),
         r = 'bilinear')

clay5 <- rast("./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_06-1m_reproj.tif")
clay5 <- crop(clay5, SEQ)
plot(clay5) # Check this looks right

writeRaster(clay5, "./00_Data/Environmental_data/Outputs/Soil_clay/SEQ_Clay_06-1m.tif")



# 1 to 2m

clay6 <- rast("./00_Data/Environmental_data/Soil_clay/CLY_100_200_EV_N_P_AU_TRN_N_20210902.tif")
clay6 <- crop(clay6, QLD)
plot(clay6) # Check this worked correctly

writeRaster(clay6, "./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_1-2m.tif")

gdalwarp(srcfile = "./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_1-2m.tif",
         dstfile = "./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_1-2m_reproj.tif",
         t_srs = 'EPSG:3577',
         tr = c(30,30),
         r = 'bilinear')

clay6 <- rast("./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_1-2m_reproj.tif")
clay6 <- crop(clay6, SEQ)
plot(clay6) # Check this looks right

writeRaster(clay6, "./00_Data/Environmental_data/Outputs/Soil_clay/SEQ_Clay_1-2m.tif")



# For our purposes, an average of percent soil clay is sufficient as this would give us an idea of the nutrients contained within the soil but we do not need super detailed information of the nutrients at the particular depth 
soil_clay <- terra::mean(clay0to0.5, clay2, clay3, clay4, clay5, clay6)
soil_clay # Check
plot(soil_clay) # Check this looks right
names(soil_clay) <- "%_clay"

writeRaster(soil_clay, "./00_Data/Environmental_data/Outputs/Soil_clay/SEQ_soilclay.tif")




# 2. Crop data to the same extent and resolve resolution issues for BioClim data ----
# To be able to work with this data all together we need to fix the extents as they are not the same despite all having been cropped by the same spatial extent of SEQ. 
# When viewing these layers overlayed in ArcGIS we can see that  slope top extent is smallest, FPC bottom extent is smallest, and the BioClim variables have the smallest x extents


ext(SEQ) #1881028, 2121562, -3247627, -2660998


# Using terra::crop() is not cropping these layers correctly so use gdalwarp instead

gdalwarp(srcfile = './00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_IBRA_freq_hydrographical_mask.tif',
         './00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_IBRA_freq_hydrographical_mask_cropped.tif',
         te = c(1881028, -3247627, 2121562, -2660998),
         tr = c(30,30),
         r = 'bilinear')



gdalwarp(srcfile = './00_Data/Fire_data/Outputs/Sentinel/Sentinel_ff_hydrographical_mask_SEQ_IBRA_focal.tif',
         './00_Data/Fire_data/Outputs/Sentinel/Sentinel_ff_hydrographical_mask_SEQ_IBRA_focal_cropped.tif',
         te = c(1881028, -3247627, 2121562, -2660998))


gdalwarp(src = './00_Data/Environmental_data/Outputs/FPC/FPC_all.tif',
         './00_Data/Environmental_data/Outputs/FPC/FPC_all_cropped.tif',
         te = c(1881028, -3247627, 2121562, -2660998))


gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/BioClim/Diurnal_temp_meanSEQ.tif',
         './00_Data/Environmental_data/Outputs/BioClim/Diurnal_temp_meanSEQ_cropped.tif',
         te = c(1881028, -3247627, 2121562, -2660998),
         tr = c(30,30),
         r = 'bilinear') # Also fix the resolution issue we know occurs for this data


gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/BioClim/precipseason_SEQ_IBRA.tif',
         './00_Data/Environmental_data/Outputs/BioClim/precipseason_SEQ_IBRA_cropped.tif',
         te = c(1881028, -3247627, 2121562, -2660998),
         tr = c(30,30),
         r = 'bilinear') # Also fix the resolution issue we know occurs for this data


gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/BioClim/tempseason_SEQ_IBRA.tif',
         './00_Data/Environmental_data/Outputs/BioClim/tempseason_SEQ_IBRA_cropped.tif',
         te = c(1881028, -3247627, 2121562, -2660998),
         tr = c(30,30),
         r = 'bilinear') # Also fix the resolution issue we know occurs for this data


gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/TWI/SEQ_TWI.tif',
         './00_Data/Environmental_data/Outputs/TWI/SEQ_TWI_cropped.tif',
         te = c(1881028, -3247627, 2121562, -2660998))



gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Soil_clay/SEQ_soilclay.tif',
         './00_Data/Environmental_data/Outputs/Soil_clay/SEQ_soilclay_cropped.tif',
         te = c(1881028, -3247627, 2121562, -2660998))

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/DEM/SEQ_DEM_reproj.tif',
         './00_Data/Environmental_data/Outputs/DEM/SEQ_DEM_reproj_cropped.tif',
         te = c(1881028, -3247627, 2121562, -2660998),
         tr = c(30,30),
         r = 'bilinear')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/DEM/SEQ_slope.tif',
         './00_Data/Environmental_data/Outputs/DEM/SEQ_slope_cropped.tif',
         te = c(1881028, -3247627, 2121562, -2660998))


gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/DEM/SEQaspect.tif',
         './00_Data/Environmental_data/Outputs/DEM/SEQaspect_cropped.tif',
         te = c(1881028, -3247627, 2121562, -2660998))


gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/DEM/SEQ_TPI.tif',
         './00_Data/Environmental_data/Outputs/DEM/SEQ_TPI_cropped.tif',
         te = c(1881028, -3247627, 2121562, -2660998))




# 3. Interpolate missing data
# Each raster has missing data, we would like to fill these gaps for any spatial layers where these gaps are small. Basically this means every layer but FPC as there are very large gaps. We could do FPC but would nee to have a good look at what is produced


# Read in the newly cropped environmental predictors data and give each raster a meaningful name

TWI <- rast('./00_Data/Environmental_data/Outputs/TWI/SEQ_TWI_cropped.tif')
names(TWI) <- "Topo_wetness_index"

tempseason <- rast('./00_Data/Environmental_data/Outputs/BioClim/Tempseason_SEQ_IBRA_cropped.tif')
names(tempseason) <- "temp_seasonality"

precipseason <- rast('./00_Data/Environmental_data/Outputs/BioClim/precipseason_SEQ_IBRA_cropped.tif')
names(precipseason) <- "precip_seasonality"

diurnal_temp <- rast('./00_Data/Environmental_data/Outputs/BioClim/Diurnal_temp_meanSEQ_cropped.tif')
names(diurnal_temp) <- "diurnal_temp_seasonality"


FPC <- rast('./00_Data/Environmental_data/Outputs/FPC/FPC_all_cropped.tif')
names(FPC) <- "average_foliage_proj_cover"

soil_clay <- rast('./00_Data/Environmental_data/Outputs/Soil_clay/SEQ_soilclay_cropped.tif')
names(soil_clay) <- "percent_clay"

slope <- rast('./00_Data/Environmental_data/Outputs/DEM/SEQ_slope_cropped.tif')

aspect <- rast('./00_Data/Environmental_data/Outputs/DEM/SEQaspect_cropped.tif')

topo_position <- rast('00_Data/Environmental_data/Outputs/DEM/SEQ_TPI_cropped.tif')
names(topo_position) <- "topo_position"

elev <- rast('./00_Data/Environmental_data/Outputs/DEM/SEQ_DEM_reproj_cropped.tif')
names(elev) <- "elevation"

# Use focal to fill in only the NA values

TWI_foc <- focal(TWI, fun = "mean", na.policy = "only", na.rm = T)
range(unique(TWI$Topo_wetness_index))
TWI_foc # The same range of values
names(TWI_foc) <- "Topo_wetness_index"
writeRaster(TWI_foc, './00_Data/Environmental_data/Outputs/TWI/SEQ_TWI_cropped_focal.tif')


soil_cl_foc <- focal(soil_clay, fun = "mean", na.poly = "only", na.rm = T)
range(unique(soil_clay$percent_clay))
unique(is.na(soil_clay$percent_clay))
soil_cl_foc # The minimum has increased. 
plet(soil_cl_foc)
names(soil_cl_foc) <- "percent_soil_clay"
writeRaster(soil_cl_foc, './00_Data/Environmental_data/Outputs/Soil_clay/SEQ_soilclay_cropped_focal.tif')

slop_foc <- focal(slope, fun = "mean", na.policy = "only", na.rm = T)
range(unique(slope$slope))
slop_foc
names(slop_foc) <- "slope"
writeRaster(slop_foc, './00_Data/Environmental_data/Outputs/DEM/SEQ_slope_cropped_focal.tif')

aspect_foc <- focal(aspect, fun = "mean", na.policy = "only", na.rm = T)
range(unique(aspect$aspect))
aspect_foc
names(aspect_foc) <- "aspect"
writeRaster(aspect_foc, './00_Data/Environmental_data/Outputs/DEM/SEQaspect_cropped_focal.tif')

topo_foc <- focal(topo_position, fun = "modal", na.policy = "only", na.rm = T) # This is the only case where we will not use the average of nearby cells, we will instead use the modal 
range(unique(topo_position$topo_position))
topo_foc
names(topo_foc) <- "Topo_position_index"
writeRaster(topo_foc, './00_Data/Environmental_data/Outputs/DEM/SEQ_TPI_cropped_focal.tif')


elev_foc <- focal(elev, fun = "mean", na.policy = "only", na.rm = T)
range(unique(elev_foc$focal_mean))
elev_foc
names(elev_foc) <- "elevation"
writeRaster(elev_foc, './00_Data/Environmental_data/Outputs/DEM/SEQ_DEM_reproj_cropped_focal.tif')

temp_foc <- focal(tempseason, fun = "mean", na.policy = "only", na.rm = T)
range(unique(tempseason$temp_seasonality))
temp_foc
names(temp_foc) <- "temperature_seasonality"
writeRaster(temp_foc, './00_Data/Environmental_data/Outputs/BioClim/tempseason_SEQ_IBRA_cropped_focal.tif')

precip_foc <- focal(precipseason, fun = "mean", na.policy = "only", na.rm = T)
range(unique(precipseason$precip_seasonality))
precip_foc
names(precip_foc) <- "precipitation_seasonality"
writeRaster(precip_foc, './00_Data/Environmental_data/Outputs/BioClim/precipseason_SEQ_IBRA_cropped_focal.tif')

diurnal_foc <- focal(diurnal_temp, fun = "mean", na.policy = "only", na.rm = T)
range(unique(diurnal_temp$diurnal_temp_seasonality))
diurnal_foc
names(diurnal_foc) <- "diurnal_temp_seasonality"
writeRaster(diurnal_foc, './00_Data/Environmental_data/Outputs/BioClim/Diurnal_temp_meanSEQ_cropped_focal.tif')


FPC_foc <- focal(FPC, fun = 'mean', na.policy = 'only', na.rm = T)
range(unique(FPC$average_foliage_proj_cover))
FPC_foc
plot(FPC_foc)
plot(FPC)
names(FPC_foc) <- "Foliage_proj_cover"
writeRaster(FPC_foc, './00_Data/Environmental_data/Outputs/FPC/FPC_all_cropped_focal.tif')

# While FPC does appear to have large areas of no data, using focal is not replacing all these areas but just expanding to fill in some of these regions with data. White patches could be related to water bodies or some other land cover 




# 4. Mask raster environmental data by hydrological features, we do not want predictions or to train the model using data where we know a fire would not burn
canal <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Canal_SEQ_IBRA.gpkg')
lake <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Lakes_SEQ_IBRA.gpkg')
pond <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Ponds_SEQ_IBRA.gpkg')
reservoir <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Reservoirs_SEQ_IBRA.gpkg')
watercourse <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Watercourses_SEQ_IBRA.gpkg')

# Also in some cases will need to mask by the coastline
Aus <- vect('./00_Data/Australia_shapefile/STE11aAust.shp')
coast <- subset(Aus, Aus$STATE_NAME == "Queensland")
coast <- project(coast, 'EPSG:3577')
coast <- crop(coast, e)

TWI <- rast('./00_Data/Environmental_data/Outputs/TWI/SEQ_TWI_cropped_focal.tif')
tempseason <- rast('./00_Data/Environmental_data/Outputs/BioClim/tempseason_SEQ_IBRA_cropped_focal.tif')    
precipseason <- rast('./00_Data/Environmental_data/Outputs/BioClim/precipseason_SEQ_IBRA_cropped_focal.tif')  
diurnal_temp <- rast('./00_Data/Environmental_data/Outputs/BioClim/Diurnal_temp_meanSEQ_cropped_focal.tif')  
FPC <- rast('./00_Data/Environmental_data/Outputs/FPC/FPC_all_cropped_focal.tif')  
soil_clay <- rast('./00_Data/Environmental_data/Outputs/Soil_clay/SEQ_soilclay_cropped_focal.tif')  
slope <- rast('./00_Data/Environmental_data/Outputs/DEM/SEQ_slope_cropped_focal.tif')  
aspect <- rast('./00_Data/Environmental_data/Outputs/DEM/SEQaspect_cropped_focal.tif')  
topo_position <- rast('./00_Data/Environmental_data/Outputs/DEM/SEQ_TPI_cropped_focal.tif') 
elev <- rast('./00_Data/Environmental_data/Outputs/DEM/SEQ_DEM_reproj_cropped_focal.tif')

TWI <- mask(TWI, canal, inverse = T)
TWI <- mask(TWI, lake, inverse = T)
TWI <- mask(TWI, pond, inverse = T)
TWI <- mask(TWI, reservoir, inverse = T)
TWI <- mask(TWI, watercourse, inverse = T)
plot(TWI)


tempseason <- mask(tempseason, canal, inverse = T)
tempseason <- mask(tempseason, lake, inverse = T)
tempseason <- mask(tempseason, pond, inverse = T)
tempseason <- mask(tempseason, reservoir, inverse = T)
tempseason <- mask(tempseason, watercourse, inverse = T)


precipseason <- mask(precipseason, canal, inverse = T)
precipseason <- mask(precipseason, lake, inverse = T)
precipseason <- mask(precipseason, pond, inverse = T)
precipseason <- mask(precipseason, reservoir, inverse = T)
precipseason <- mask(precipseason, watercourse, inverse = T)


diurnal_temp <- mask(diurnal_temp, canal, inverse = T)
diurnal_temp <- mask(diurnal_temp, lake, inverse = T)
diurnal_temp <- mask(diurnal_temp, pond, inverse = T)
diurnal_temp <- mask(diurnal_temp, reservoir, inverse = T)
diurnal_temp <- mask(diurnal_temp, watercourse, inverse = T)



FPC <- mask(FPC, canal, inverse = T)
FPC <- mask(FPC, lake, inverse = T)
FPC <- mask(FPC, pond, inverse = T)
FPC <- mask(FPC, reservoir, inverse = T)
FPC <- mask(FPC, watercourse, inverse = T)


soil_clay <- mask(soil_clay, canal, inverse = T)
soil_clay <- mask(soil_clay, lake, inverse = T)
soil_clay <- mask(soil_clay, pond, inverse = T)
soil_clay <- mask(soil_clay, reservoir, inverse = T)
soil_clay <- mask(soil_clay, watercourse, inverse = T)


slope <- mask(slope, canal, inverse = T)
slope <- mask(slope, lake, inverse = T)
slope <- mask(slope, pond, inverse = T)
slope <- mask(slope, reservoir, inverse = T)
slope <- mask(slope, watercourse, inverse = T)
slope <- mask(slope, coast)


aspect <- mask(aspect, canal, inverse = T)
aspect <- mask(aspect, lake, inverse = T)
aspect <- mask(aspect, pond, inverse = T)
aspect <- mask(aspect, reservoir, inverse = T)
aspect <- mask(aspect, watercourse, inverse = T)
aspect <- mask(aspect, coast)


topo_position <- mask(topo_position, canal, inverse = T)
topo_position <- mask(topo_position, lake, inverse = T)
topo_position <- mask(topo_position, pond, inverse = T)
topo_position <- mask(topo_position, reservoir, inverse = T)
topo_position <- mask(topo_position, watercourse, inverse = T)
topo_position <- mask(topo_position, coast)


elev <- mask(elev, canal, inverse = T)
elev <- mask(elev, lake, inverse = T)
elev <- mask(elev, pond, inverse = T)
elev <- mask(elev, reservoir, inverse = T)
elev <- mask(elev, watercourse, inverse = T)
elev <- mask(elev, coast)


# 5. Combine the data into a single raster for later use ----
# Make sure all the names are meaningful


names(TWI) <- "Topo_wetness_index"

names(tempseason) <- "Temp_seasonality"

names(precipseason) <- "Precip_seasonality"

names(diurnal_temp) <- "Diurnal_temp_seasonality"

names(FPC) <- "Average_foliage_proj_cover"

names(soil_clay) <- 'Percent_clay'

names(slope) <- "Slope"

names(aspect) <- "Aspect"

names(topo_position) <- "Topo_position_index"

names(elev) <- "Elevation"


# Combine the data and save the output
# Need to make sure QPWS fire data has the same extent as this other data
library(gdalUtilities)
gdalwarp('./00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_IBRA_freq_hydrographical_mask_cropped.tif',
         './00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_IBRA_freq_hydrographical_mask_cropped_reproj.tif',
         te = c(1881028, -3247627, 2121562, -2660998))

QPWS_SEQ_ff <- rast('./00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_IBRA_freq_hydrographical_mask_cropped_reproj.tif')
names(QPWS_SEQ_ff) <- "QPWS_firefreq"


# GAMs and GLMs will not produce predictions for areas where any layer has NA values. We need to replace NA values in the QPWS dataset to be 0s despite the fact that these may not be true absences. We will be replacing these NA values with zeros when we extract values from each spatial layer to create the dataset for feeding into our modelling algorithms.

QPWS_SEQ_ff[is.na(QPWS_SEQ_ff)] <- 0
unique(QPWS_SEQ_ff)

predictors <- c(QPWS_SEQ_ff, TWI, tempseason, precipseason, diurnal_temp, FPC, soil_clay, slope, aspect, topo_position, elev)
writeRaster(predictors, './00_Data/SDM_data/predictors.tif', overwrite = T)



