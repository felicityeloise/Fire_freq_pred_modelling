# Written by Felicity Charles
# Date:1/08/2023
# Update: 30/07/2025

##### Fire frequency analysis ----
# This script gathers together environmental data needed for predicting  fire frequency

# R version 4.5.1

# 1. Load required packages ----
library(terra) # terra_1.7-78 # Updated to 1.8-60
library(dplyr) # dplyr_1.1.4 
library(landform) # landform_0.2
library(httr) # httr_1.4.7
library(raster) # raster 3.6-23 # Updated to 3.6-32
library(sf) # s3 1.0-14 # Updated to 1.0-21
library(httr) # httr 1.4.7

# Instalhttr# Install gdalUtilities for working with large datasets
#library(devtools)
#install_github("JoshOBrien/gdalUtilities")
library(gdalUtilities) # gdalUtilities_1.2.5


# 1. Read in, crop and aggregate environmental data ----
# Create an SEQ vector 
SEQ <- vect('./00_Data/SEQ_bound/SEQ_IBRA.gpkg')
SEQ2 <- project(SEQ, 'EPSG:4326')


# Use SEQ to create a template raster for resampling 
template <- rast(SEQ)
res(template) <- c(30,30)
ext(template) <- ext(SEQ)


# 1.1 Rainfall ----
# 1.1.1 Load daily rainfall data
load_SILO_rain <- function(year, dest_dir = './00_Data/Environmental_data/SILO_Rainfall'){
  if(!dir.exists(dest_dir)) dir.create(dest_dir)
  
  base_url <- "https://s3-ap-southeast-2.amazonaws.com/silo-open-data/Official/annual/daily_rain"
  file_name <- paste0(year, ".daily_rain.nc")
  url <- file.path(base_url, file_name)
  local_file <- file.path(dest_dir, file_name)
  
  # Download file
  download.file(url, local_file, mode = 'wb') # Method = 'wb' preserves the byte content of the original file
  
  # Load raster
  r <- rast(local_file) %>% 
    crop(aoi)
  assign(paste0("rain_", year), r, envir = .GlobalEnv) # Add rasters to the environment as they are processed
}
# Loop over 1987 to 2023
options(timeout = 400)
aoi <- SEQ2 
years <- 1987:2023 # When an error is thrown and the connection fails, just update this code to start from the year which failed and re-run L57-L58
invisible(lapply(years, load_SILO_rain))
rain_1987; plot(rain_1987); gc()

# 1.1.2 Calculate rainfall seasonality
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
avg_monthly_rain1987; gc() # Check

# Calculate seasonality 
for(i in years){
  month_rast_name <- paste0("avg_monthly_rain", i)
  month_rast <- get(month_rast_name)
  
  # Calculate seasonality in same manner as BioClim BIO4: standard deviation of monthly temps × 100
  si_rain <- app(month_rast, fun = function(x) sd(x, na.rm = TRUE) * 100)
  names(si_rain) <- paste0("rain_seasonality_", i)
  assign(paste0("rain_seasonality_", i), si_rain)
}
rain_seasonality_1987; rain_seasonality_2023; gc() # Check these seem reasonable

# Combine into a raster stack
rain_seasonality <- c(rain_seasonality_1987, rain_seasonality_1988, rain_seasonality_1989, rain_seasonality_1990, rain_seasonality_1991, rain_seasonality_1992, rain_seasonality_1993, rain_seasonality_1994, rain_seasonality_1995, rain_seasonality_1996, rain_seasonality_1997, rain_seasonality_1998, rain_seasonality_1999, rain_seasonality_2000, rain_seasonality_2001, rain_seasonality_2002, rain_seasonality_2003, rain_seasonality_2004, rain_seasonality_2005, rain_seasonality_2006, rain_seasonality_2007, rain_seasonality_2008, rain_seasonality_2009, rain_seasonality_2010, rain_seasonality_2011, rain_seasonality_2012, rain_seasonality_2013, rain_seasonality_2014, rain_seasonality_2015, rain_seasonality_2016, rain_seasonality_2017, rain_seasonality_2018, rain_seasonality_2019, rain_seasonality_2020, rain_seasonality_2021, rain_seasonality_2022, rain_seasonality_2023)
rain_seasonality

# Reproject 
rain_seasonality[rain_seasonality < 0] <- NA # SILO uses values like -32768 as values for missing data which means these values are not recognised as NAs during reprojection.
writeRaster(rain_seasonality, './00_Data/Environmental_data/Outputs/SILO_Rainfall/Rainfall_seasonality_SEQ_IBRA.tif', overwrite = T)

# Rescale
gdalwarp('./00_Data/Environmental_data/Outputs/SILO_Rainfall/Rainfall_seasonality_SEQ_IBRA.tif',
         './00_Data/Environmental_data/Outputs/SILO_Rainfall/Rainfall_seasonality_SEQ_IBRA_reproj.tif', 
         t_srs = 'EPSG:3577',
         tr = c(30,30), 
         r = 'bilinear')
# Check
rain_seasonalityr <- rast('./00_Data/Environmental_data/Outputs/SILO_Rainfall/Rainfall_seasonality_SEQ_IBRA_reproj.tif')


# Calculate long-term average
rain_seasonality_avg <- mean(rain_seasonalityr, na.rm = T)
rain_seasonality_avg
plet(rain_seasonality_avg)

writeRaster(rain_seasonality_avg, './00_Data/Environmental_data/Outputs/SILO_Rainfall/Average_Rain_seasonality_SEQ_IBRA.tif', overwrite = T)
rm(list = setdiff(ls(), c("SEQ", "SEQ2", "template"))); gc()



# 1.2 Temperature ----
# 1.2.1 Min temperature
load_SILO_mintemp <- function(year, dest_dir = './00_Data/Environmental_data/SILO_Min_temp'){
  if(!dir.exists(dest_dir)) dir.create(dest_dir)
  
  base_url <- "https://s3-ap-southeast-2.amazonaws.com/silo-open-data/Official/annual/min_temp"
  file_name <- paste0(year, ".min_temp.nc")
  url <- file.path(base_url, file_name)
  local_file <- file.path(dest_dir, file_name)
  
  # Download file
  download.file(url, local_file, mode = 'wb') # Method = 'wb' preserves the byte content of the original file
  
  # Load raster
  r <- rast(local_file) %>% 
    crop(aoi)
  assign(paste0("mintemp_", year), r, envir = .GlobalEnv) # Add rasters to the environment as they are processed
}
# Loop over 1987 to 2023
options(timeout = 400)
aoi <- SEQ2
years <- 1987:2023 # When an error is thrown and the connection fails, just update this code to start from the year which failed and re-run L124-L125
invisible(lapply(years, load_SILO_mintemp))
mintemp_1987; plot(mintemp_1987); gc() # Values are per day

# 1.2.1 Max temperature
load_SILO_maxtemp <- function(year, dest_dir = './00_Data/Environmental_data/SILO_Max_temp'){
  if(!dir.exists(dest_dir)) dir.create(dest_dir)
  
  base_url <- "https://s3-ap-southeast-2.amazonaws.com/silo-open-data/Official/annual/max_temp"
  file_name <- paste0(year, ".max_temp.nc")
  url <- file.path(base_url, file_name)
  local_file <- file.path(dest_dir, file_name)
  
  # Download file
  download.file(url, local_file, mode = 'wb') # Method = 'wb' preserves the byte content of the original file
  
  # Load raster
  r <- rast(local_file) %>% 
    crop(aoi)
  assign(paste0("maxtemp_", year), r, envir = .GlobalEnv) # Add rasters to the environment as they are processed
}
# Loop over 1987 to 2023
aoi <- SEQ2
years <- 1987:2023 # When an error is thrown and the connection fails, just update this code to start from the year which failed and re-run L147-L148
invisible(lapply(years, load_SILO_maxtemp))
maxtemp_1987; plot(maxtemp_1987) # Values are per day

rm(list = setdiff(ls(), c("SEQ", "template", 'years'))); gc()

# 1.2.2 Calculate mean daily temperature
years <- 1987:2023
for (i in years ){
  
  tmin_file <- file.path(paste0('./00_Data/Environmental_data/SILO_Min_temp/',i,".min_temp.nc")) # Get the raster file for minimum temperature raster for year i 
  tmax_file <- file.path(paste0('./00_Data/Environmental_data/SILO_Max_temp/',i,".max_temp.nc")) # Get the raster file for maximum temperature raster for year i
  
  # Load rasters
  tmin <- rast(tmin_file)
  tmax <- rast(tmax_file)
  
  # Calculate daily mean temperature
  tavg <- (tmin + tmax) / 2 
  
  # Create month groups
  dates <- time(tavg) # Create a time series of times when sampling occurred
  month_groups <- format(dates, "%Y-%m") # Produce the month groups based on 4-digit year and month
  
  # Calculate monthly average temperature
  monthly_avg <- tapp(tavg, index = month_groups, fun = mean, na.rm = T) 
  assign(paste0('avg_monthly_temp', i), monthly_avg, envir = .GlobalEnv) # Add rasters to the environment
}

avg_monthly_temp1987; plot(avg_monthly_temp1987) # Check

# 1.2.3 Calculate temperature seasonality
# NOTE: to interpret values in terms of °C, divide the value by 100 (e.g., 350/100 = 3.5°C variation in annual temperature range)

for(i in years){
  month_rast_name <- paste0("avg_monthly_temp", i)
  month_rast <- get(month_rast_name)
  
  # Calculate seasonality in same manner as BioClim BIO4: standard deviation of monthly temps × 100
  si_temp <- app(month_rast, fun = function(x) sd(x, na.rm = TRUE) * 100)
  names(si_temp) <- paste0("temp_seasonality_", i)
  assign(paste0("temp_seasonality_", i), si_temp, envir = .GlobalEnv)
}
temp_seasonality_1987; temp_seasonality_2023 # Check these seem reasonable
plot(temp_seasonality_1987)

# Combine into a raster stack - note this will mean that we have a seperate value for each year.
temp_seasonality <- c(temp_seasonality_1987, temp_seasonality_1988, temp_seasonality_1989, temp_seasonality_1990, temp_seasonality_1991, temp_seasonality_1992, temp_seasonality_1993, temp_seasonality_1994, temp_seasonality_1995, temp_seasonality_1996, temp_seasonality_1997, temp_seasonality_1998, temp_seasonality_1999, temp_seasonality_2000, temp_seasonality_2001, temp_seasonality_2002, temp_seasonality_2003, temp_seasonality_2004, temp_seasonality_2005, temp_seasonality_2006, temp_seasonality_2007, temp_seasonality_2008, temp_seasonality_2009, temp_seasonality_2010, temp_seasonality_2011, temp_seasonality_2012, temp_seasonality_2013, temp_seasonality_2014, temp_seasonality_2015, temp_seasonality_2016, temp_seasonality_2017, temp_seasonality_2018, temp_seasonality_2019, temp_seasonality_2020, temp_seasonality_2021, temp_seasonality_2022, temp_seasonality_2023)
temp_seasonality
range(unique(temp_seasonality))
writeRaster(temp_seasonality, './00_Data/Environmental_data/Outputs/SILO_Temperature/Temp_seasonality_SEQ_IBRA.tif', overwrite = T)
# Temperature seasonality data is now structured such that we have 36 years of data which allows us to capture year-to-year variations across not just long-term averages (e.g., BioClim BIO4).

# Calculate the long-term average
temp_seasonality[temp_seasonality < 0] <- NA # SILO uses values like -32768 as values for missing data which means these values are not recognised as NAs during reprojection.
temp_seasonality_avg <- mean(temp_seasonality, na.rm = T)
temp_seasonality_avg

# Too memory intensive to follow same process as for the rainfall data
temp_seasonality_avgr <- project(temp_seasonality_avg, 'EPSG:3577', res = 30, method = 'bilinear')
temp_seasonality_avgr
plet(temp_seasonality_avgr)

writeRaster(temp_seasonality_avgr, './00_Data/Environmental_data/Outputs/SILO_Temperature/Average_Temp_seasonality_SEQ_IBRA.tif', overwrite = T)
gc()

# 1.3 Topographic Wetness Index ----
# Use CSIRO web coverage service for Topographic Wetness Index, no other download methods (e.g., downlaod.file, vsicurl, gdal_translate or wms_url will work on Mac OS)

# Queensland bounding box (xmin, ymin, xmax, ymax)
qld_bbox <- c(138.0, -29.2, 153.6, -9.0)

# Download TWI data via WCS
wcs_url <- paste0(
  "https://www.asris.csiro.au/arcgis/services/TERN/SRTM_attributes_3s_ACLEP_AU/MapServer/WCSServer",
  "?service=WCS&version=1.0.0&request=GetCoverage&coverage=14&format=GeoTIFF",
  "&bbox=", paste(qld_bbox, collapse=","),
  "&crs=EPSG:4326&width=2000&height=1500"
)


# Load raster
TWI <- rast(wcs_url)
NAflag(TWI) # Check - NaN

writeRaster(TWI, './00_Data/Environmental_data/TWI/TWI.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/TWI/TWI.tif', # Source file
         dstfile = './00_Data/Environmental_data/Outputs/TWI/TWI_reproj_v2.tif', # Destination file
         t_srs = 'EPSG:3577', # CRS to be transformed to  
         tr = c(30,30), # Resolution to be transformed to
         r = 'bilinear') # Use bilinear interpolation for continuous data 

TWI <- rast('./00_Data/Environmental_data/Outputs/TWI/TWI_reproj_v2.tif')

TWI <- crop(TWI, SEQ)
plot(TWI)
writeRaster(TWI, './00_Data/Environmental_data/Outputs/TWI/TWI_SEQ_IBRA_reproj.tif')


# 1.4 Foliage projective cover ----
# Download FPC data from Queensland Spatial Catalogue and load into R
# For the earlier datasets, 2012, 2013, and 2014, we need to do some more adjustments to this data - metadata states that data ranges between 100-200 which is equivalent to 0-100% FPC. Values erroneously predicted above 100% or below 0% have been classed as above 200 and below 100 respectively. Zero values indicate NULL data. The 2014 data actually seems to be ranging between 88-213. Post 2014, values range between 0-100 which would denote the % cover without any further changes being required. 

FPC12 <- rast('./00_Data/Environmental_data/FPC/DP_QLD_FPC2012_WOODED_VEG.tif')
FPC13 <- rast('./00_Data/Environmental_data/FPC/QLD_FPC2013_WOODED_VEG.tif')
FPC14 <- rast('./00_Data/Environmental_data/FPC/DP_QLD_FPC2014.tif')

unique(FPC12$DP_QLD_FPC2012_WOODED_VEG)
NAflag(FPC12) # NaN
unique(FPC13$QLD_FPC2013_WOODED_VEG)
NAflag(FPC13) #NaN
unique(FPC14$DP_QLD_FPC2014)
NAflag(FPC14) #NaN

# Create matrices for reclassification
# For 2012 and 2013, image values are FPC + 100 such that FPC of 5% has a value of 105. So in these cases a value of 100 = 0% FPC
# 2012 data
reclas12 = matrix(
  c(100:200),
  nrow = 101,
  ncol = 1
) %>% 
  cbind(0:100)

# 2013 data
reclas13 = matrix(
  c(100:200),
  nrow = 101,
  ncol = 1
) %>%
  cbind(0:100)


# 2014 data
A14 = matrix(
  c(88:99, 201:213),
  nrow = 25,
  ncol = 2)
A14[,2] <- 0

B14 = matrix(
  c(100:200),
  nrow = 101,
  ncol = 1
)
B14 <- cbind(B14, 0:100)

reclas14 <- rbind(A14, B14)


# Now reclassify these datasets
FPC12r <- classify(FPC12, rcl = reclas12)
FPC12r; plot(FPC12r)
FPC13r <- classify(FPC13, rcl = reclas13)
FPC13r; plot(FPC13r)
FPC14r <- classify(FPC14, rcl = reclas14)
FPC14r; plot(FPC14r) # Check how this looks


# Crop to SEQ
FPC12seq <- crop(FPC12r, SEQ) %>% 
  project('EPSG:3577') %>% 
  resample(template, method = 'bilinear')
FPC12seq;plot(FPC12seq)
writeRaster(FPC12seq, './00_Data/Environmental_data/Outputs/FPC/FPC12_SEQ_IBRA.tif')

FPC13seq <- crop(FPC13r, SEQ) %>% 
  project('EPSG:3577')
writeRaster(FPC13seq, './00_Data/Environmental_data/Outputs/FPC/FPC13_SEQ_IBRA.tif')

FPC14seq <- crop(FPC14r, SEQ) %>% 
  project('EPSG:3577')
FPC14seq; plot(FPC14seq)

writeRaster(FPC14seq, './00_Data/Environmental_data/Outputs/FPC/FPC2014_SEQ_IBRA.tif')


# Check the data
FPC18 <- rast('./00_Data/Environmental_data/FPC/DP_QLD_S2_WOODY_FPC_2018.tif')
NAflag(FPC18) # NaN

FPC18 <- crop(FPC18, SEQ) %>% 
  aggregate(fact = 3, fun = 'mean') # Calculate average FPC for the cell when aggregating to a coarser resolution
FPC18; plot(FPC18) # Check how this looks

writeRaster(FPC18, './00_Data/Environmental_data/Outputs/FPC/FPC18_SEQ_IBRA.tif')


FPC19 <- rast('./00_Data/Environmental_data/FPC/DP_QLD_S2_WOODY_FPC_2019.tif')
NAflag(FPC19) # NaN
FPC19 <- crop(FPC19, SEQ) %>% 
  aggregate(fact = 3, fun = 'mean')
FPC19; plot(FPC19)
writeRaster(FPC19, './00_Data/Environmental_data/Outputs/FPC/FPC19_SEQ_IBRA.tif')


FPC20 <- rast('./00_Data/Environmental_data/FPC/DP_QLD_S2_WOODY_FPC_2020.tif')
NAflag(FPC20) # NaN
FPC20 <- crop(FPC20, SEQ) %>% 
  aggregate(fact = 3, fun = 'mean')
FPC20; plot(FPC20)
writeRaster(FPC20, './00_Data/Environmental_data/Outputs/FPC/FPC20_SEQ_IBRA.tif')


FPC21 <- rast('./00_Data/Environmental_data/FPC/DP_QLD_S2_FPC_2021.tif')
NAflag(FPC21) # NaN
FPC21 <- crop(FPC21, SEQ) %>% 
  aggregate(fact = 3, fun = 'mean') %>% 
  project('EPSG:3577') # The coordinate reference system has not been read in the same manner as the others so we will need to fix this
FPC21; plot(FPC21)
writeRaster(FPC21, './00_Data/Environmental_data/Outputs/FPC/FPC21_SEQ_IBRA.tif')


FPC22 <- rast('./00_Data/Environmental_data/FPC/DP_QLD_S2_FPC_2022.tif')
NAflag(FPC22) # NaN
FPC22 <- crop(FPC22, SEQ) %>% 
  aggregate(fact = 3, fun = 'mean')
FPC22; plot(FPC22)
writeRaster(FPC22, './00_Data/Environmental_data/Outputs/FPC/FPC22_SEQ_IBRA.tif')


FPC23 <- rast('./00_Data/Environmental_data/FPC/DP_QLD_S2_FPC_2023.tif')
NAflag(FPC23) # NaN
FPC23 <- crop(FPC23, SEQ) %>% 
  aggregate(fact = 3, fun = 'mean')
FPC23; plot(FPC23)
writeRaster(FPC23, './00_Data/Environmental_data/Outputs/FPC/FPC23_SEQ_IBRA.tif')


# 1.4.2 Combine the FPC data into one raster ----
FPC12seq <- resample(FPC12seq, FPC18, method = 'bilinear') # Need the extents to match
FPC13seq <- resample(FPC13seq, FPC18, method = 'bilinear') # Need the extents to match
FPC14seq <- resample(FPC14seq, FPC18, method = 'bilinear') # Need the extents to match
# Check
FPC12seq; FPC13seq; FPC14seq

writeRaster(FPC12seq, './00_Data/Environmental_data/Outputs/FPC/FPC12_SEQ_IBRA_resampled.tif')
writeRaster(FPC13seq, './00_Data/Environmental_data/Outputs/FPC/FPC13_SEQ_IBRA_resampled.tif')
writeRaster(FPC14seq, './00_Data/Environmental_data/Outputs/FPC/FPC14_SEQ_IBRA_resampled.tif')

# To conduct a sensitivity analysis, lets create few different FPC datasets. Rather than averaging to get a psuedo-long-term average, lets just create raster stacks. For now as we will compare FPC and NDVI, we will just use the long-term average
FPC_all <- c(FPC12seq, FPC13seq, FPC14seq, FPC18, FPC19, FPC20, FPC21, FPC22, FPC23)
FPC_all; plot(FPC_all)
writeRaster(FPC_all, './00_Data/Environmental_data/Outputs/FPC/FPC_ALL_SEQ_IBRA.tif')

FPC <- mean(FPC_all, na.rm = T)
writeRaster(FPC, './00_Data/Environmental_data/Outputs/FPC/Average_FPC_SEQ_IBRA.tif')



# We will test using NDVI data instead of FPC to see which is better. NDVI data is at a much coarser scale, 5km than the intial scale of FPC data which may disadvantage this data but NDVI has a better temporal scale relative to the fire data than FPC.
# Download BOM NDVI data for May 1992 to December 2018 from https://www.bom.gov.au/jsp/awap/ndvi/index.jsp?colour=colour&map=ndviave&year=2018&month=8&period=month&area=nat

aoi <- SEQ
years <- 1992:2018
months <- c("January", "February", "March", "April", "May", "June", 
            "July", "August", "September", "October", "November", "December")


terraOptions(memfrac = 0.8, progress = 2)
for (year in years) {
  monthly_rasters <- list()
  
  for (month in months) {
    # Skip months before May 1992
    if (year == 1992 & month %in% c("January", "February", "March", "April")) next 
    
    # Create file path
    file_path <- file.path("./00_Data/Environmental_data/BOM_NDVI", paste0(year, "_", month))
    
    # Check if file exists and process
    if (file.exists(file_path)) {
      # Load raster from first extracted file and crop
      r <- rast(file_path)
      r[r < -100] <- NA # NA values coded as -999 or -9999
      monthly_rasters[[month]] <- crop(project(r, template, method = 'bilinear'), aoi)
    } 
  }
  # Calculate yearly average from monthly rasters and save to disk
  if (length(monthly_rasters) > 0) {
    yearly_mean <- mean(rast(monthly_rasters), na.rm = TRUE)
    writeRaster(yearly_mean, paste0("./00_Data/Environmental_data/Outputs/BOM_NDVI/NDVI_", year, "_average.tif"), overwrite = TRUE)
  }
  
  tmpFiles(remove = TRUE)
  gc()
}

# Create a raster stack of yearly average NDVI rasters
years <- 1992:2018
annual_rasters <- rast(paste0("./00_Data/Environmental_data/Outputs/BOM_NDVI/NDVI_", years, "_average.tif"))
names(annual_rasters) <- paste0("NDVI_", years)
annual_rasters

# Calcualte long-term average
NDVI <- mean(annual_rasters, na.rm = T)
plot(NDVI)
NDVI
writeRaster(NDVI, './00_Data/Environmental_data/Outputs/BOM_NDVI/Average_NDVI.tif', overwrite = T)

# The rsi package could be used to access satellite imagery for 2019 onwards, however the below was tested for a singular month and did not finish running in ~24 hours.
# For 2019 to 2023 data, use the rsi package
#SEQ_sf <- st_as_sf(SEQ) # Convert our area of interest to a sf object


#sent19 <- get_stac_data(
#aoi = SEQ_sf,
#start_date = "2019-01-01",
#end_date = "2019-01-31", 
#pixel_x_size = 30,
#pixel_y_size = 30,
#stac_source = "https://planetarycomputer.microsoft.com/api/stac/v1",
#collection = "sentinel-2-l2a",
#bands = c("B04", "B08"),
#cloud_mask = FALSE,
#composite_function = "median")  # Combine scenes into one image with median removing outliers/clouds



# 1.5 Elevation data ----
# Download the DEM data file from https://data.gov.au/data/dataset/9a9284b6-eb45-4a13-97d0-91bf25f1187b/resource/7653b920-5334-4267-8df3-5d55b21f05ec and unzip the file using an alternative file management software other than that installed as part of computing system such as WinZip.
# Create the QLD polygon for cropping
QLD <- vect('./00_Data/Australia_shapefile/STE11aAust.shp')
QLD <- subset(QLD, QLD$STATE_NAME == "Queensland") %>%
  project('EPSG:4326')
# Read in while cropping the DEM data to QLD 
DEM <- rast('./00_Data/Environmental_data/69816/srtm-1sec-dem-v1-COG.tif') %>% 
  crop(QLD)
writeRaster(DEM, './00_Data/Environmental_data/Outputs/DEM/QLD_DEM.tif')


gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/DEM/QLD_DEM.tif',
         dstfile = './00_Data/Environmental_data/Outputs/DEM/QLD_DEM_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30),
         r = 'bilinear')
DEM <- rast('./00_Data/Environmental_data/Outputs/DEM/QLD_DEM_reproj.tif') %>% 
  crop(SEQ)
plot(DEM)
writeRaster(DEM, './00_Data/Environmental_data/Outputs/DEM/SEQ_IBRA_DEM_reproj.tif')

# 1.5.1 Calculate slope, aspect and topographic position index (TPI) ----
slope <- terrain(DEM, v = "slope", unit = "degrees")
slope; plot(slope) #Check
writeRaster(slope, './00_Data/Environmental_data/Outputs/DEM/slope_SEQ_IBRA.tif') # Save output


aspect <- terrain(DEM, v = "aspect")
plot(aspect); aspect
writeRaster(aspect, './00_Data/Environmental_data/Outputs/DEM/aspect_SEQ_IBRA.tif')


TP1 <- landform(DEM, class.type = "slope.position")
TPI <- TP1$all

# Here TPI values are as follows:
# 1 = Valley
# 2 = Lower slope
# 3 = Flat slope
# 4 = Middle slope
# 5 = Upper slope
# 6 = Ridge
plot(TPI);TPI
writeRaster(TPI, "./00_Data/Environmental_data/Outputs/DEM/TPI_SEQ_IBRA.tif")




# 1.6 Soil nutrients using soil % clay ----
# This variable is included as nutrients influence plant growth which would then influence the occurrence of fire
# Go to https://portal.tern.org.au/ and sign-in then generate API key by clicking on your Profile > TERN ACCOUNT > Create API key on side menu > API key name = 'Soils_access' > Request API Key
# Also, have to use the httr package to read the file. 


my_api_key <- "paste API key here"
get_slga_clay <- function(depth_code, aoi, api_key) {
  url <- paste0("https://data.tern.org.au/model-derived/slga/NationalMaps/SoilAndLandscapeGrid/CLY/v2/CLY_", 
                depth_code, "_EV_N_P_AU_TRN_N_20210902.tif")
  
  temp_file <- tempfile(fileext = ".tif")
  
  # Simple download with long timeout
  GET(url,
      authenticate("apikey", api_key, type = "basic"),
      write_disk(temp_file, overwrite = TRUE),
      timeout(7200),  # 2 hours - plenty of time
      progress())
  
  # Load, crop, and save
  r <- rast(temp_file)
  r[r < -100] <- NA # NA values stored as -9999
  r <- project(r, 'EPSG:3577') %>% 
    crop(aoi)
  assign(paste0("soil_clay_", depth_code), r, envir = .GlobalEnv)
}
depth_codes <- c("000_005", "005_015", "015_030", "030_060", "060_100", "100_200")
for(i in depth_codes) {
  get_slga_clay(i, aoi = SEQ, my_api_key)
}

writeRaster(soil_clay_000_005, './00_Data/Environmental_data/Soil_clay/SEQ_IBRA_soil_clay1.tif')
writeRaster(soil_clay_005_015, './00_Data/Environmental_data/Soil_clay/SEQ_IBRA_soil_clay2.tif')
writeRaster(soil_clay_015_030, './00_Data/Environmental_data/Soil_clay/SEQ_IBRA_soil_clay3.tif')
writeRaster(soil_clay_030_060, './00_Data/Environmental_data/Soil_clay/SEQ_IBRA_soil_clay4.tif')
writeRaster(soil_clay_060_100, './00_Data/Environmental_data/Soil_clay/SEQ_IBRA_soil_clay5.tif')
writeRaster(soil_clay_100_200, './00_Data/Environmental_data/Soil_clay/SEQ_IBRA_soil_clay6.tif')

# For our purposes, an average of percent soil clay is sufficient as this would give us an idea of the nutrients contained within the soil but we do not need super detailed information of the nutrients at each particular depth 
soil_clay_000_005 <- rast('./00_Data/Environmental_data/Soil_clay/SEQ_IBRA_soil_clay1.tif')
NAflag(soil_clay_000_005) # NaN
soil_clay_005_015 <- rast('./00_Data/Environmental_data/Soil_clay/SEQ_IBRA_soil_clay2.tif')
soil_clay_015_030 <- rast('./00_Data/Environmental_data/Soil_clay/SEQ_IBRA_soil_clay3.tif')
soil_clay_030_060 <- rast('./00_Data/Environmental_data/Soil_clay/SEQ_IBRA_soil_clay4.tif')
soil_clay_060_100 <- rast('./00_Data/Environmental_data/Soil_clay/SEQ_IBRA_soil_clay5.tif')
soil_clay_100_200 <- rast('./00_Data/Environmental_data/Soil_clay/SEQ_IBRA_soil_clay6.tif')

soil_clay <- terra::mean(soil_clay_000_005, soil_clay_005_015, soil_clay_015_030, soil_clay_030_060, soil_clay_060_100, soil_clay_100_200)
soil_clay; plot(soil_clay) # Check this looks right

soil_clay <- resample(soil_clay, template, method = 'bilinear')

names(soil_clay) <- "%_clay_avg"

writeRaster(soil_clay, "./00_Data/Environmental_data/Outputs/Soil_clay/SEQ_IBRA_soilclay.tif")



# 1.7 Broad Vegetation groups ----

BVG <- vect('./00_Data/Environmental_data/Remnant_2021_broad_veg_groups/Remnant_broad_vegetation_groups.shp') %>%
  project('EPSG:3577') %>% 
  crop(SEQ)
BVG # Check if there is anything else we need to remove?

BVG <- BVG[, 5] # Only keep the relevant columns
unique(BVG$dbvg5m) # We can remove water
BVG <- BVG[!BVG$dbvg5m == 'water']


# We need to rasterise this data for use with the other datasets
unique(BVG$dbvg5m) # rasterisation will require all dbvg5m codes to be numeric

# Assign numeric codes
BVG$dbvg5m_coded <- ifelse(BVG$dbvg5m == "non-remnant", 17,
                           ifelse(BVG$dbvg5m == "plantation", 18,
                                  as.numeric(BVG$dbvg5m)))
unique(BVG$dbvg5m_coded)

rtemp <- rast(xmin = 1881028, xmax = 2121594, ymin = -3442593, ymax = -2656210, res = 30, crs = 'EPSG:3577')

BVG_rast <- rasterize(BVG, rtemp, field = 'dbvg5m_coded', fun = 'min') # Rasterize the data
BVG_rast # Check how this looks
writeRaster(BVG_rast, './00_Data/Environmental_data/Outputs/BVG/BVG.tif')

BVG_rast <- rast('./00_Data/Environmental_data/Outputs/BVG/BVG.tif')
unique(BVG_rast$BVG)



# 2. Reinforce spatial extent and resolution ----
# To be able to work with this data all together we need to fix issues with the spatial extents as they are not the same despite all having been cropped by the same spatial extent of SEQ and also resolution.

# Re-load all data
SEQ

Sentinel_ff <- rast('./00_Data/Fire_data/Outputs/Sentinel/SEQ_IBRA/Sentinel_ff_hydrographical_mask_SEQ_IBRA_focal.tif') %>% 
  print()

QPWS_ff <- rast('./00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_IBRA_freq_hydrographical_mask.tif') %>% 
  print()

TWI <- rast('./00_Data/Environmental_data/Outputs/TWI/TWI_SEQ_IBRA_reproj.tif') %>% 
  print()

temp_season <- rast('./00_Data/Environmental_data/Outputs/SILO_Temperature/Average_Temp_seasonality_SEQ_IBRA.tif') %>% 
  print()

precip_season <- rast('./00_Data/Environmental_data/Outputs/SILO_Rainfall/Average_Rain_seasonality_SEQ_IBRA.tif') %>% 
  print()

FPC <- rast('./00_Data/Environmental_data/Outputs/FPC/Average_FPC_SEQ_IBRA.tif') %>% 
  print()

NDVI <- rast('./00_Data/Environmental_data/Outputs/BOM_NDVI/Average_NDVI.tif') %>% 
  print()

elevation <- rast('./00_Data/Environmental_data/Outputs/DEM/SEQ_IBRA_DEM_reproj.tif') %>% 
  print()

slope <- rast('./00_Data/Environmental_data/Outputs/DEM/slope_SEQ_IBRA.tif') %>% 
  print()

aspect <- rast('./00_Data/Environmental_data/Outputs/DEM/aspect_SEQ_IBRA.tif') %>% 
  print()

TPI <- rast('./00_Data/Environmental_data/Outputs/DEM/TPI_SEQ_IBRA.tif') %>% 
  print()

soil_clay <- rast('./00_Data/Environmental_data/Outputs/Soil_clay/SEQ_IBRA_soilclay.tif') %>% 
  print()

BVG_rast <- rast('./00_Data/Environmental_data/Outputs/BVG/BVG.tif') %>% 
  print()

# While all the datasets were cropped using terra::crop(), extents do not match which is required for being able to run models in subsequent steps

gdalwarp(srcfile = './00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_IBRA_freq_hydrographical_mask.tif',
         './00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_IBRA_freq_hydrographical_mask_cropped.tif',
         te = c(1881028, -3247627, 2121562, -2660998),
         tr = c(30,30),
         r = 'bilinear')


gdalwarp(srcfile = './00_Data/Fire_data/Outputs/Sentinel/SEQ_IBRA/Sentinel_ff_hydrographical_mask_SEQ_IBRA_focal.tif',
         './00_Data/Fire_data/Outputs/Sentinel/Sentinel_ff_hydrographical_mask_SEQ_IBRA_focal_cropped.tif',
         te = c(1881028, -3247627, 2121562, -2660998))


gdalwarp(src = './00_Data/Environmental_data/Outputs/FPC/Average_FPC_SEQ_IBRA.tif',
         './00_Data/Environmental_data/Outputs/FPC/Average_FPC_SEQ_IBRA_cropped.tif',
         te = c(1881028, -3247627, 2121562, -2660998))



gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/SILO_Rainfall/Average_Rain_seasonality_SEQ_IBRA.tif',
         './00_Data/Environmental_data/Outputs/SILO_Rainfall/Average_Rain_seasonality_SEQ_IBRA_reproj_cropped.tif',
         te = c(1881028, -3247627, 2121562, -2660998)) 

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/SILO_Temperature/Average_Temp_seasonality_SEQ_IBRA.tif',
         './00_Data/Environmental_data/Outputs/SILO_Temperature/Average_Temp_seasonality_SEQ_IBRA_reproj_cropped.tif',
         te = c(1881028, -3247627, 2121562, -2660998)) 



gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/TWI/TWI_SEQ_IBRA_reproj.tif',
         './00_Data/Environmental_data/Outputs/TWI/TWI_SEQ_IBRA_reproj_cropped.tif',
         te = c(1881028, -3247627, 2121562, -2660998))



gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Soil_clay/SEQ_IBRA_soilclay.tif',
         './00_Data/Environmental_data/Outputs/Soil_clay/SEQ_IBRA_soilclay_cropped.tif',
         te = c(1881028, -3247627, 2121562, -2660998))


gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/BOM_NDVI/Average_NDVI.tif',
         './00_Data/Environmental_data/Outputs/BOM_NDVI/Average_NDVI_reproj_cropped.tif',
         te = c(1881028, -3247627, 2121562, -2660998))


gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/DEM/SEQ_IBRA_DEM_reproj.tif',
         './00_Data/Environmental_data/Outputs/DEM/SEQ_IBRA_DEM_reproj_cropped.tif',
         te = c(1881028, -3247627, 2121562, -2660998),
         tr = c(30,30),
         r = 'bilinear')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/DEM/slope_SEQ_IBRA.tif',
         './00_Data/Environmental_data/Outputs/DEM/slope_SEQ_IBRA_cropped.tif',
         te = c(1881028, -3247627, 2121562, -2660998))


gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/DEM/aspect_SEQ_IBRA.tif',
         './00_Data/Environmental_data/Outputs/DEM/aspect_SEQ_IBRA_cropped.tif',
         te = c(1881028, -3247627, 2121562, -2660998))


gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/DEM/TPI_SEQ_IBRA.tif',
         './00_Data/Environmental_data/Outputs/DEM/TPI_SEQ_IBRA_cropped.tif',
         te = c(1881028, -3247627, 2121562, -2660998))


gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/BVG/BVG.tif',
         './00_Data/Environmental_data/Outputs/BVG/BVG_cropped.tif',
         te = c(1881028, -3247627, 2121562, -2660998))


# 3. Interpolate missing data
# Each raster has missing data, we would like to fill these gaps for any spatial layers where these gaps are small. Basically this means every layer but FPC as there are very large gaps. We could do FPC but would need to have a good look at what is produced


# Read in the newly cropped environmental predictors data and give each raster a meaningful name
QPWS_ff <- rast('./00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_IBRA_freq_hydrographical_mask_cropped.tif') %>% 
  print()

TWI <- rast('./00_Data/Environmental_data/Outputs/TWI/TWI_SEQ_IBRA_reproj_cropped.tif')
names(TWI) <- "Topo_wetness_index"

tempseason <- rast('./00_Data/Environmental_data/Outputs/SILO_Temperature/Average_Temp_seasonality_SEQ_IBRA_reproj_cropped.tif')
names(tempseason) <- "Avg_Temperature_seasonality"


precipseason <- rast('./00_Data/Environmental_data/Outputs/SILO_Rainfall/Average_Rain_seasonality_SEQ_IBRA_reproj_cropped.tif')
names(precipseason) <- "Avg_Precipitation_seasonality"

FPC <- rast('./00_Data/Environmental_data/Outputs/FPC/Average_FPC_SEQ_IBRA_cropped.tif')
names(FPC) <- "Avg_FPC"

NDVI <- rast('./00_Data/Environmental_data/Outputs/BOM_NDVI/Average_NDVI_reproj_cropped.tif')
names(NDVI) <- "Avg_NDVI"

soil_clay <- rast('./00_Data/Environmental_data/Outputs/Soil_clay/SEQ_IBRA_soilclay_cropped.tif')
names(soil_clay) <- "percent_clay"

slope <- rast('./00_Data/Environmental_data/Outputs/DEM/slope_SEQ_IBRA_cropped.tif')

aspect <- rast('./00_Data/Environmental_data/Outputs/DEM/aspect_SEQ_IBRA_cropped.tif')

topo_position <- rast('00_Data/Environmental_data/Outputs/DEM/TPI_SEQ_IBRA_cropped.tif')
names(topo_position) <- "topo_position"

elev <- rast('./00_Data/Environmental_data/Outputs/DEM/SEQ_IBRA_DEM_reproj_cropped.tif')
names(elev) <- "elevation"

# Use focal to fill in only the NA values

TWI_foc <- focal(TWI, fun = "mean", na.policy = "only", na.rm = T)
range(unique(TWI$Topo_wetness_index))
TWI_foc # The same range of values
names(TWI_foc) <- "Topo_wetness_index"
writeRaster(TWI_foc, './00_Data/Environmental_data/Outputs/TWI/SEQ_IBRA_TWI_cropped_focal.tif')


soil_cl_foc <- focal(soil_clay, fun = "mean", na.poly = "only", na.rm = T)
range(unique(soil_clay$percent_clay)); unique(is.na(soil_clay$percent_clay)); soil_cl_foc 
names(soil_cl_foc) <- "percent_soil_clay"
writeRaster(soil_cl_foc, './00_Data/Environmental_data/Outputs/Soil_clay/SEQ_IBRA_soilclay_cropped_focal.tif')

slop_foc <- focal(slope, fun = "mean", na.policy = "only", na.rm = T)
range(unique(slope$slope)); slop_foc
names(slop_foc) <- "slope"
writeRaster(slop_foc, './00_Data/Environmental_data/Outputs/DEM/SEQ_IBRA_slope_cropped_focal.tif')

aspect_foc <- focal(aspect, fun = "mean", na.policy = "only", na.rm = T)
range(unique(aspect$aspect)); aspect_foc
names(aspect_foc) <- "aspect"
writeRaster(aspect_foc, './00_Data/Environmental_data/Outputs/DEM/SEQ_IBRA_aspect_cropped_focal.tif')

topo_foc <- focal(topo_position, fun = "modal", na.policy = "only", na.rm = T) # This is the only case where we will not use the average of nearby cells, we will instead use the modal 
range(unique(topo_position$topo_position)); topo_foc
names(topo_foc) <- "Topo_position_index"
writeRaster(topo_foc, './00_Data/Environmental_data/Outputs/DEM/SEQ_IBRA_TPI_cropped_focal.tif')


elev_foc <- focal(elev, fun = "mean", na.policy = "only", na.rm = T)
range(unique(elev_foc$focal_mean)); elev_foc
names(elev_foc) <- "elevation"
writeRaster(elev_foc, './00_Data/Environmental_data/Outputs/DEM/SEQ_IBRA_DEM_reproj_cropped_focal.tif')

temp_foc <- focal(tempseason, fun = "mean", na.policy = "only", na.rm = T)
temp_foc
writeRaster(temp_foc, './00_Data/Environmental_data/Outputs/SILO_Temperature/Average_tempseason_SEQ_IBRA_reproj_cropped_focal.tif')

precip_foc <- focal(precipseason, fun = "mean", na.policy = "only", na.rm = T)
writeRaster(precip_foc, './00_Data/Environmental_data/Outputs/SILO_Rainfall/Average_precipseason_SEQ_IBRA_reproj_cropped_focal.tif')


NDVI_foc <- focal(NDVI, fun = "mean", na.policy = "only", na.rm = T)
writeRaster(NDVI_foc, './00_Data/Environmental_data/Outputs/BOM_NDVI/Average_NDVI_SEQ_IBRA_cropped_focal.tif', overwrite = T)

# While FPC does appear to have large areas of no data, using focal is not replacing all these areas but just expanding to fill in some of these regions with data. White patches could be related to water bodies or some other land cover 
FPC_foc <- focal(FPC, fun = 'mean', na.policy = 'only', na.rm = T)
FPC; FPC_foc
names(FPC_foc) <- "Avg_FPC"
writeRaster(FPC_foc, './00_Data/Environmental_data/Outputs/FPC/Average_FPC_IBRA_cropped_focal.tif')


# 4. Mask environmental raster data by hydrological features, we do not want predictions or to train the model using data where we know a fire would not burn
canal <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Canal_SEQ_IBRA.gpkg')
lake <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Lakes_SEQ_IBRA.gpkg')
pond <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Ponds_SEQ_IBRA.gpkg')
reservoir <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Reservoirs_SEQ_IBRA.gpkg')
watercourse <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Watercourses_SEQ_IBRA.gpkg')

# Also in some cases will need to mask by the coastline
Aus <- vect('./00_Data/Australia_shapefile/STE11aAust.shp')
coast <- subset(Aus, Aus$STATE_NAME == "Queensland") %>% 
  project("EPSG:3577") %>% 
  crop(SEQ)

TWI <- rast('./00_Data/Environmental_data/Outputs/TWI/SEQ_IBRA_TWI_cropped_focal.tif')
tempseason <- rast('./00_Data/Environmental_data/Outputs/SILO_Temperature/Average_tempseason_SEQ_IBRA_reproj_cropped_focal.tif')
precipseason <- rast('./00_Data/Environmental_data/Outputs/SILO_Rainfall/Average_precipseason_SEQ_IBRA_reproj_cropped_focal.tif')
FPC <- rast('./00_Data/Environmental_data/Outputs/FPC/FPC_ALL_IBRA_cropped_focal.tif')  
NDVI <- rast('./00_Data/Environmental_data/Outputs/BOM_NDVI/Average_NDVI_SEQ_IBRA_cropped_focal.tif')
soil_clay <- rast('./00_Data/Environmental_data/Outputs/Soil_clay/SEQ_IBRA_soilclay_cropped_focal.tif')  
slope <- rast('./00_Data/Environmental_data/Outputs/DEM/SEQ_IBRA_slope_cropped_focal.tif')  
aspect <- rast('./00_Data/Environmental_data/Outputs/DEM/SEQ_IBRA_aspect_cropped_focal.tif')  
topo_position <- rast('./00_Data/Environmental_data/Outputs/DEM/SEQ_IBRA_TPI_cropped_focal.tif') 
elev <- rast('./00_Data/Environmental_data/Outputs/DEM/SEQ_IBRA_DEM_reproj_cropped_focal.tif')
BVG <- rast('./00_Data/Environmental_data/Outputs/BVG/BVG_cropped.tif')

TWI <- mask(TWI, canal, inverse = T)
TWI <- mask(TWI, lake, inverse = T)
TWI <- mask(TWI, pond, inverse = T)
TWI <- mask(TWI, reservoir, inverse = T)
TWI <- mask(TWI, watercourse, inverse = T)
TWI <- mask(TWI, coast)
plot(TWI)
names(TWI) <- "Topo_wetness_index"
writeRaster(TWI, './00_Data/Environmental_data/Outputs/TWI/SEQ_IBRA_TWI_cropped_focal_masked.tif', overwrite = T)

tempseason <- mask(tempseason, canal, inverse = T)
tempseason <- mask(tempseason, lake, inverse = T)
tempseason <- mask(tempseason, pond, inverse = T)
tempseason <- mask(tempseason, reservoir, inverse = T)
tempseason <- mask(tempseason, watercourse, inverse = T)
plot(tempseason)
tempseason <- resample(tempseason, template, method = 'bilinear') %>% 
  print()
tempseason <- mask(tempseason, coast)
names(tempseason) <- "Avg_Temperature_seasonality"
writeRaster(tempseason, './00_Data/Environmental_data/Outputs/SILO_Temperature/Average_Temp_seasonality_SEQ_IBRA_reproj_cropped_focal_masked.tif', overwrite = T)



precipseason <- mask(precipseason, canal, inverse = T)
precipseason <- mask(precipseason, lake, inverse = T)
precipseason <- mask(precipseason, pond, inverse = T)
precipseason <- mask(precipseason, reservoir, inverse = T)
precipseason <- mask(precipseason, watercourse, inverse = T)
precipseason <- resample(precipseason, template, method = 'bilinear') %>% 
  print()
precipseason <- mask(precipseason, coast)
names(precipseason) <-  "Avg_Precipitation_seasonality"
writeRaster(precipseason, './00_Data/Environmental_data/Outputs/SILO_Rainfall/Average_precipseason_SEQ_IBRA_reproj_cropped_focal_masked.tif', overwrite = T)


FPC <- mask(FPC, canal, inverse = T)
FPC <- mask(FPC, lake, inverse = T)
FPC <- mask(FPC, pond, inverse = T)
FPC <- mask(FPC, reservoir, inverse = T)
FPC <- mask(FPC, watercourse, inverse = T)
FPC <- mask(FPC, coast)
plot(FPC)
names(FPC) <- "Avg_FPC"
writeRaster(FPC, './00_Data/Environmental_data/Outputs/FPC/Average_FPC_SEQ_IBRA_cropped_focal_masked.tif', overwrite = T)


NDVI <-  mask(NDVI, canal, inverse = T)
NDVI <- mask(NDVI, lake, inverse = T)
NDVI <- mask(NDVI, pond, inverse = T)
NDVI <- mask(NDVI, reservoir, inverse = T)
NDVI <- mask(NDVI, watercourse, inverse = T)
NDVI <- mask(NDVI, coast)
plot(NDVI)
names(NDVI) <- "Avg_NDVI"
writeRaster(NDVI, './00_Data/Environmental_data/Outputs/BOM_NDVI/Average_NDVI_cropped_focal_masked.tif', overwrite = T)


soil_clay <- mask(soil_clay, canal, inverse = T)
soil_clay <- mask(soil_clay, lake, inverse = T)
soil_clay <- mask(soil_clay, pond, inverse = T)
soil_clay <- mask(soil_clay, reservoir, inverse = T)
soil_clay <- mask(soil_clay, watercourse, inverse = T)
soil_clay <- mask(soil_clay, coast)
plot(soil_clay)
names(soil_clay) <- 'Percent_clay'
writeRaster(soil_clay, './00_Data/Environmental_data/Outputs/Soil_clay/SEQ_IBRA_soilclay_cropped_focal_masked.tif', overwrite = T)


slope <- mask(slope, canal, inverse = T)
slope <- mask(slope, lake, inverse = T)
slope <- mask(slope, pond, inverse = T)
slope <- mask(slope, reservoir, inverse = T)
slope <- mask(slope, watercourse, inverse = T)
slope <- mask(slope, coast)
plot(slope)
names(slope) <- "Slope"
writeRaster(slope, './00_Data/Environmental_data/Outputs/DEM/SEQ_IBRA_slope_cropped_focal_masked.tif', overwrite = T)

aspect <- mask(aspect, canal, inverse = T)
aspect <- mask(aspect, lake, inverse = T)
aspect <- mask(aspect, pond, inverse = T)
aspect <- mask(aspect, reservoir, inverse = T)
aspect <- mask(aspect, watercourse, inverse = T)
aspect <- mask(aspect, coast)
plot(aspect)
names(aspect) <- "Aspect"
writeRaster(aspect, './00_Data/Environmental_data/Outputs/DEM/SEQ_IBRA_aspect_cropped_focal_masked.tif', overwrite = T)


topo_position <- mask(topo_position, canal, inverse = T)
topo_position <- mask(topo_position, lake, inverse = T)
topo_position <- mask(topo_position, pond, inverse = T)
topo_position <- mask(topo_position, reservoir, inverse = T)
topo_position <- mask(topo_position, watercourse, inverse = T)
topo_position <- mask(topo_position, coast)
plot(topo_position)
names(topo_position) <- "Topo_position_index"
writeRaster(topo_position, './00_Data/Environmental_data/Outputs/DEM/SEQ_IBRA_TPI_cropped_focal_masked.tif', overwrite = T)


elev <- mask(elev, canal, inverse = T)
elev <- mask(elev, lake, inverse = T)
elev <- mask(elev, pond, inverse = T)
elev <- mask(elev, reservoir, inverse = T)
elev <- mask(elev, watercourse, inverse = T)
elev <- mask(elev, coast)
elev <- resample(elev, template, method = 'bilinear') %>% 
  print()
plot(elev)
names(elev) <- "Elevation"
writeRaster(elev, './00_Data/Environmental_data/Outputs/DEM/SEQ_IBRA_DEM_cropped_focal_masked.tif', overwrite = T)


BVG <- mask(BVG, canal, inverse = T)
BVG <- mask(BVG, lake, inverse = T)
BVG <- mask(BVG, pond, inverse = T)
BVG <- mask(BVG, reservoir, inverse = T)
BVG <- mask(BVG, watercourse, inverse = T)
BVG <- mask(BVG, coast)
BVG <- resample(BVG, template, method = 'bilinear') %>% 
  print()
plot(BVG)
writeRaster(BVG, './00_Data/Environmental_data/Outputs/BVG/BVG_cropped_masked.tif', overwrite = T)




# 5. Combine the data into a single raster for later use ----
QPWS_SEQ_ff <- rast('./00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_IBRA_freq_hydrographical_mask_cropped.tif') %>% 
  resample(template, method = 'bilinear')
names(QPWS_SEQ_ff) <- "QPWS_firefreq"

QPWS_SEQ_ff[is.na(QPWS_SEQ_ff)] <- 0 # GAMs and GLMs will not produce predictions for areas where any layer has NA values so we need to reaplce NA values in the QPWS data with 0s despite the fact that these may not be true absences.
range(unique(QPWS_SEQ_ff))
plot(QPWS_SEQ_ff)

QPWS_SEQ_ff <- mask(QPWS_SEQ_ff, coast)
plot(QPWS_SEQ_ff)
writeRaster(QPWS_SEQ_ff, './00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_IBRA_freq_hydrographical_mask_cropped.tif', overwrite = T)

predictors <- c(QPWS_SEQ_ff, TWI, tempseason, precipseason, FPC, soil_clay, slope, aspect, topo_position, elev, NDVI, BVG)
predictors
writeRaster(predictors, './00_Data/SDM_data/predictors_SEQ_IBRA.tif', overwrite = T)



