# Written by Felicity Charles
# Date:1/08/2023


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


# To improve our predictive accuracy for QPWS based fire frequency outside QPWS estates lets include environmental data
# 1. Read in, crop and aggregate environmental data ----
# Create an SEQ vector 
e <- ext(1902033, 2111776, -3257627, -2954985)
SEQ <- as.polygons(e, 'EPSG:3577')
writeVector(SEQ, './00_Data/SEQ_bound/SEQ.gpkg')

SEQ <- vect('./00_Data/SEQ_bound/SEQ.gpkg')

# 1.1 Topographic Wetness Index ----

TWI <- rast("/vsicurl/https://s3.data.csiro.au/dapprd/000005588v002/data/TopographicWetnessIndex_1_arcsecond_resolution/mosaic/twi_1s.tif?response-content-disposition=attachment%3B%20filename%3D%22twi_1s.tif%22&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20240311T015651Z&X-Amz-SignedHeaders=host&X-Amz-Expires=172800&X-Amz-Credential=WJU445XD13UTAMAHK6LU%2F20240311%2FCDC%2Fs3%2Faws4_request&X-Amz-Signature=83bcbd4b9498bff94952d508f45207fae833ccb1df4656744d60cab81e602b5f") # To get this link need to click download on the file from https://data.csiro.au/collection/csiro:5588 and then as it is downloading open up the downloads drop down > right click on the file > copy download link. Need to get a new link each time this line is run otherwise it will not work.

writeRaster(TWI, overwrite = T, "./00_Data/Environmental_data/Outputs/TWI/TWI.tif")

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/TWI/TWI.tif', # Source file
         dstfile = './00_Data/Environmental_data/Outputs/TWI/TWI_reproj.tif', # Destination file
         t_srs = 'EPSG:3577', # CRS to be transformed to  
         tr = c(30,30)) # Resolution to be transformed to

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
dtm1 <- disagg(dtm, fact = 30)
dtm2 <- crop(dtm1, SEQ)
plot(dtm2) # Check how this looks
dtm2 # Check crs and resolution

writeRaster(dtm2, './00_Data/Environmental_data/Outputs/BioClim/Diurnal_temp_meanSEQ.tif')


gdalwarp(srcfile = './00_Data/Environmental_data/BioClim/wc2.1_30s_bio_4.tif',
         dstfile = '00_Data/Environmental_data/Outputs/BioClim/Tempseason.tif',
         t_srs = 'EPSG:3577')


tempseason <- rast('./00_Data/Environmental_data/Outputs/BioClim/Tempseason.tif')
range(unique(tempseason$Tempseason)) 

temp <- crop(tempseason, SEQ)
plot(temp)
tempr <- disagg(temp, fact = 30)
tempr1 <- crop(tempr, SEQ)
plot(tempr1) # Check to see if this looks right
tempr1 # We have the right coordinate reference system and resolution is nearly correct. 

writeRaster(tempr1, '00_Data/Environmental_data/Outputs/BioClim/Tempseason_SEQ.tif')


gdalwarp(srcfile = './00_Data/Environmental_data/BioClim/wc2.1_30s_bio_15.tif',
         dstfile = '00_Data/Environmental_data/Outputs/BioClim/precipseason.tif',
         t_srs = 'EPSG:3577')

# Need to crop and change resolution
precipseason <- rast('00_Data/Environmental_data/Outputs/BioClim/precipseason.tif')
precip <- crop(precipseason, SEQ)
precipr <- disagg(precip, fact = 30)
precipr1 <- crop(precipr, SEQ)
plot(precipr1)
precipr1

writeRaster(precipr1, '00_Data/Environmental_data/Outputs/BioClim/precipseason_SEQ.tif')




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

# Need to aggregate the data to a coarser resoltuion, from 10m to 30m, and then we also need to crop the data to SEQ

FPC18_seq <- crop(FPC18, SEQ)

FPC18seq <- terra::aggregate(FPC18_seq, fact = 3)
FPC18seq # Check how this looks
plot(FPC18seq)

writeRaster(FPC18seq, './00_Data/Environmental_data/Outputs/FPC/FPC18_SEQ.tif')


FPC19 <- rast('./00_Data/Environmental_data/FPC/DP_QLD_S2_WOODY_FPC_2019.tif')
FPC19_seq <- crop(FPC19, SEQ)
FPC19seq <- terra::aggregate(FPC19_seq, fact = 3)
FPC19seq
plot(FPC19seq)
writeRaster(FPC19seq, './00_Data/Environmental_data/Outputs/FPC/FPC19_SEQ.tif')


FPC20 <- rast('./00_Data/Environmental_data/FPC/DP_QLD_S2_WOODY_FPC_2020.tif')
FPC20_seq <- crop(FPC20, SEQ)
FPC20seq <- terra::aggregate(FPC20_seq, fact = 3)
FPC20seq
plot(FPC20seq)
writeRaster(FPC20seq, './00_Data/Environmental_data/Outputs/FPC/FPC20_SEQ.tif')


FPC21 <- rast('./00_Data/Environmental_data/FPC/DP_QLD_S2_FPC_2021.tif') # The coordinate reference system has not been read in the same manner as the others so we will need to fix this

FPC21_seq <- crop(FPC21, SEQ)
FPC21seq <- project(FPC21_seq, 'EPSG:3577')
FPC21seq <- terra::aggregate(FPC21_seq, fact = 3)
FPC21seq
plot(FPC21seq)
writeRaster(FPC21seq, './00_Data/Environmental_data/Outputs/FPC/FPC21_SEQ.tif')


# 1.3.2 Combine the FPC data into one raster ----
FPC14seq <- resample(FPC14seq, FPC18seq) # Need the extents to match

#stack <- c(FPC14,FPC18,FPC19,FPC20,FPC21)
# This produces a raster with each year as a separate raster layer. But what we want, because FPC2014 is from 1988-2014 and the others are separate years, we need an average FPC value across the years.

FPC <- terra::mean(FPC14seq, FPC18seq, FPC19seq, FPC20seq, FPC21seq)
FPC
plot(FPC)


writeRaster(FPC, './00_Data/Environmental_data/Outputs/FPC/FPC_all.tif')




# 1.4 Elevation data ----
# Download the DEM ad access straight from download folder

DEM <- rast('./00_Data/Environmental_data/69816/srtm-1sec-dem-v1-COG.tif') # Work from the original data


# Need to cut down to smaller area
Aus <- vect('./00_Data/Australia_shapefile/STE11aAust.shp')
QLD <- subset(Aus, Aus$STATE_NAME == "Queensland")
QLD <- project(QLD, 'EPSG:4326')
DEMqld <- crop(DEM, QLD)

SEQ <- vect('./00_Data/SEQ_bound/SEQ.gpkg') # Make sure we have the original SEQ file here
SEQ <- project(SEQ, 'EPSG:4326')
DEMseq <- crop(DEM, SEQ)

writeRaster(DEMseq, './00_Data/Environmental_data/Outputs/DEM/SEQ_DEM.tif')
gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/DEM/SEQ_DEM.tif',
         dstfile = './00_Data/Environmental_data/Outputs/DEM/SEQ_DEM_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))
DEM <- rast('./00_Data/Environmental_data/Outputs/DEM/SEQ_DEM_reproj.tif')



# 1.4.1 Calculate slope, aspect and topographic position index (TPI) ----

slope <- terrain(DEMqld, v = "slope", unit = "degrees")
writeRaster(slope, './00_Data/Environmental_data/Outputs/DEM/slope.tif') # Save output


# Fix projection and change resolution so it can be used with the other data
gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/DEM/slope.tif',
         dstfile = './00_Data/Environmental_data/Outputs/DEM/slope_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

slope <- rast('./00_Data/Environmental_data/Outputs/DEM/slope_reproj.tif')
slope <- crop(slope, SEQ)
plot(slope)
slope

writeRaster(slope, './00_Data/Environmental_data/Outputs/DEM/SEQ_slope.tif')


aspect <- terrain(DEMqld, v = "aspect")
writeRaster(aspect, './00_Data/Environmental_data/Outputs/DEM/aspect.tif')

# Fix projection and change resolution so it can be used with the other data
gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/DEM/aspect.tif',
         dstfile = './00_Data/Environmental_data/Outputs/DEM/aspect_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

aspect <- rast('./00_Data/Environmental_data/Outputs/DEM/aspect_reproj.tif')
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
         tr = c(30,30)) 

TPI <- rast('./00_Data/Environmental_data/Outputs/DEM/TPI_reproj.tif')
TPI <- crop(TPI, SEQ)
plot(TPI)
TPI

writeRaster(TPI, './00_Data/Environmental_data/Outputs/DEM/SEQ_TPI.tif')





# 1.5 Soil nutrients using soil % clay ----
# Download data from https://esoil.io/TERNLandscapes/Public/Pages/SLGA/GetData-COGSDataStore_SLGA.html
# This variable is included as nutrients influence plant growth which would then influence the occurrence of fire
# As with solar radiation we will first crop to QLD and then run gdalwarp and then we can crop to SEQ

# 0 to 0.05m
clay0to0.05 <- rast("./00_Data/Environmental_data/Soil_clay/CLY_000_005_EV_N_P_AU_TRN_N_20210902.tif")
clay0to0.05 <- crop(clay0to0.05, QLD)
plot(clay0to0.05) # Check this worked correctly
 
writeRaster(clay0to0.05, "./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_0-05m.tif")

gdalwarp(srcfile = "./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_0-05m.tif",
         dstfile = "./00_Data/Environmental_data/Outputs/Soil_clay/QLD_Clay_0-05m_reproj.tif",
         t_srs = 'EPSG:3577',
         tr = c(30,30))

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
         tr = c(30,30))

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
         tr = c(30,30))

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
         tr = c(30,30))

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
         tr = c(30,30))

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
         tr = c(30,30))

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

e <- c(1902030, 2111790, -3257630, -2954990)


# Using terra::crop() is not cropping these layers correctly so use gdalwarp instead

gdalwarp(srcfile = './00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_freq_hydrographical_mask.tif',
         './00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_freq_hydrographical_mask_cropped.tif',
         te = c(1902030, -3257630, 2111790, -2954990),
         tr = c(30,30))



gdalwarp(srcfile = './00_Data/Fire_data/Outputs/Sentinel/Sentinel_ff_hydrographical_mask_SEQ_focal.tif',
         './00_Data/Fire_data/Outputs/Sentinel/Sentinel_ff_hydrographical_mask_SEQ_focal_cropped.tif',
         te = c(1902030, -3257630, 2111790, -2954990))


gdalwarp(src = './00_Data/Environmental_data/Outputs/FPC/FPC_all.tif',
         './00_Data/Environmental_data/Outputs/FPC/FPC_all_cropped.tif',
         te = c(1902030, -3257630, 2111790, -2954990))


gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/BioClim/Diurnal_temp_meanSEQ.tif',
         './00_Data/Environmental_data/Outputs/BioClim/Diurnal_temp_meanSEQ_cropped.tif',
         te = c(1902030, -3257630, 2111790, -2954990),
         tr = c(30,30)) # Also fix the resolution issue we know occurs for this data


gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/BioClim/precipseason_SEQ.tif',
         './00_Data/Environmental_data/Outputs/BioClim/precipseason_SEQ_cropped.tif',
         te = c(1902030, -3257630, 2111790, -2954990),
         tr = c(30,30)) # Also fix the resolution issue we know occurs for this data


gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/BioClim/tempseason_SEQ.tif',
         './00_Data/Environmental_data/Outputs/BioClim/tempseason_SEQ_cropped.tif',
         te = c(1902030, -3257630, 2111790, -2954990),
         tr = c(30,30)) # Also fix the resolution issue we know occurs for this data


gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/TWI/SEQ_TWI.tif',
         './00_Data/Environmental_data/Outputs/TWI/SEQ_TWI_cropped.tif',
         te = c(1902030, -3257630, 2111790, -2954990))



gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Soil_clay/SEQ_soilclay.tif',
         './00_Data/Environmental_data/Outputs/Soil_clay/SEQ_soilclay_cropped.tif',
         te = c(1902030, -3257630, 2111790, -2954990))

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/DEM/SEQ_DEM_reproj.tif',
         './00_Data/Environmental_data/Outputs/DEM/SEQ_DEM_reproj_cropped.tif',
         te = c(1902030, -3257630, 2111790, -2954990),
         tr = c(30,30))

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/DEM/SEQ_slope.tif',
         './00_Data/Environmental_data/Outputs/DEM/SEQ_slope_cropped.tif',
         te = c(1902030, -3257630, 2111790, -2954990))


gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/DEM/SEQaspect.tif',
         './00_Data/Environmental_data/Outputs/DEM/SEQaspect_cropped.tif',
         te = c(1902030, -3257630, 2111790, -2954990))


gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/DEM/SEQ_TPI.tif',
         './00_Data/Environmental_data/Outputs/DEM/SEQ_TPI_cropped.tif',
         te = c(1902030, -3257630, 2111790, -2954990))




# 3. Interpolate missing data
# Each raster has missing data, we would like to fill these gaps for any spatial layers where these gaps are small. Basically this means every layer but FPC as there are very large gaps. We could do FPC but would nee to have a good look at what is produced


# Read in the newly cropped environmental predictors data and give each raster a meaningful name

TWI <- rast('./00_Data/Environmental_data/Outputs/TWI/SEQ_TWI_cropped.tif')
names(TWI) <- "Topo_wetness_index"

tempseason <- rast('./00_Data/Environmental_data/Outputs/BioClim/Tempseason_SEQ_cropped.tif')
names(tempseason) <- "temp_seasonality"

precipseason <- rast('./00_Data/Environmental_data/Outputs/BioClim/precipseason_SEQ_cropped.tif')
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
writeRaster(temp_foc, './00_Data/Environmental_data/Outputs/BioClim/tempseason_SEQ_cropped_focal.tif')

precip_foc <- focal(precipseason, fun = "mean", na.policy = "only", na.rm = T)
range(unique(precipseason$precip_seasonality))
precip_foc
names(precip_foc) <- "precipitation_seasonality"
writeRaster(precip_foc, './00_Data/Environmental_data/Outputs/BioClim/precipseason_SEQ_cropped_focal.tif')

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
canal <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Canal_SEQ.gpkg')
lake <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Lakes_SEQ.gpkg')
pond <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Ponds_SEQ.gpkg')
reservoir <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Reservoirs_SEQ.gpkg')
watercourse <- vect('./00_Data/Environmental_data/Outputs/Hydrographic_features/Watercourses_SEQ.gpkg')

# Also in some cases will need to mask by the coastline
Aus <- vect('./00_Data/Australia_shapefile/STE11aAust.shp')
coast <- subset(Aus, Aus$STATE_NAME == "Queensland")
coast <- project(coast, 'EPSG:3577')
coast <- crop(coast, e)

TWI <- rast('./00_Data/Environmental_data/Outputs/TWI/SEQ_TWI_cropped_focal.tif')
tempseason <- rast('./00_Data/Environmental_data/Outputs/BioClim/tempseason_SEQ_cropped_focal.tif')    
precipseason <- rast('./00_Data/Environmental_data/Outputs/BioClim/precipseason_SEQ_cropped_focal.tif')  
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
gdalwarp('./00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_freq_hydrographical_mask_cropped.tif',
         './00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_freq_hydrographical_mask_cropped_reproj.tif',
         te = c(1902030, -3257630, 2111790, -2954990))

QPWS_SEQ_ff <- rast('./00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_freq_hydrographical_mask_cropped_reproj.tif')
names(QPWS_SEQ_ff) <- "QPWS_firefreq"


# GAMs and GLMs will not produce predictions for areas where any layer has NA values. We need to replace NA values in the QPWS dataset to be 0s despite the fact that these may not be true absences. We will be replacing these NA values with zeros when we extract values from each spatial layer to create the dataset for feeding into our modelling algorithms.

QPWS_SEQ_ff[is.na(QPWS_SEQ_ff)] <- 0
unique(QPWS_SEQ_ff)

predictors <- c(QPWS_SEQ_ff, TWI, tempseason, precipseason, diurnal_temp, FPC, soil_clay, slope, aspect, topo_position, elev)
writeRaster(predictors, './00_Data/SDM_data/predictors.tif', overwrite = T)



