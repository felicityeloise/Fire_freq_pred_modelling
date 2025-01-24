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
# Need to do some more adjustments to this data - metadata states that data ranges between 100-200 which is equivalent to 0-100% FPC. values erroneously predicted above 100% or below 0% have been classed as 200 and 100 respectively. Zero values indicate NULL data. The data actually seems to be ranging between 88-213. Post 2014, values range between 0-100 which would denote the % cover without any further changes being required. Let's take a look at the data in ArcGIS as well to make sure this is true for the 2014 dataset. Looking in ArcGIS we can see that the actual values do range between 100-200 with additional 0 for NULL data.


# Create matrices for reclassification

A = matrix(
  c(88:99),
  nrow = 12,
  ncol = 2)
A[,2] <- 0

B = matrix(
  c(100:200),
  nrow = 101,
  ncol = 1
)
B <- cbind(B, 0:100)

C = matrix(
  c(201:213),
  nrow = 13,
  ncol = 2
)
C[,2] <- 100


reclas <- rbind(A, B, C)

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
FPC14 <- resample(FPC14, FPC18) # Need the extents to match

#stack <- c(FPC14,FPC18,FPC19,FPC20,FPC21)
# This produces a raster with each year as a separate raster layer. But what we want, because FPC2014 is from 1988-2014 and the others are separate years, we need an average FPC value across the years.

FPC <- terra::mean(FPC14, FPC18, FPC19, FPC20, FPC21)
FPC
plot(FPC)


writeRaster(mean_FPC, './00_Data/Environmental_data/Outputs/FPC/FPC_all.tif')




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





# 1.5 Solar radiation ----
# Download .nc files from https://s3-ap-southeast-2.amazonaws.com/silo-open-data/Official/annual/index.html


# First crop the data down to QLD then run gdalwarp and then crop to SEQ to reduce run times
Aus <- vect('./00_Data/Australia_shapefile/STE11aAust.shp')
QLD <- subset(Aus, Aus$STATE_NAME == "Queensland")
QLD <- project(QLD, "EPSG:4326")

SEQ <- vect('./00_Data/SEQ_bound/SEQ.gpkg') # Make sure we have the original SEQ file here


# Note that the Averaged files and reproj files were deleted due to their large space requirements
# These steps take ~ 1 hour
r87 <- rast('./00_Data/Environmental_data/Solar_radiation/1987.radiation.nc')
r87 <- crop(r87, QLD)
rad_avg87 <- mean(r87)
plot(rad_avg87)
writeRaster(rad_avg87, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad87_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad87_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r87avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r87 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r87avg_reproj.tif')
plot(r87)
r87seq <- crop(r87, SEQ)
plot(r87seq) # Slight rotation issue but I am unsure of why

writeRaster(r87seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad87_seq.tif")



r88 <- rast('./00_Data/Environmental_data/Solar_radiation/1988.radiation.nc')
r88 <- crop(r88, QLD)
rad_avg88 <- mean(r88)
plot(rad_avg88)
writeRaster(rad_avg88, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad88_avg.tif')


gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad88_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r88avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r88 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r88avg_reproj.tif')

r88seq <- crop(r88, SEQ)
plot(r88seq) 

writeRaster(r88seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad88_seq.tif")


r89 <- rast('./00_Data/Environmental_data/Solar_radiation/1989.radiation.nc')
r89 <- crop(r89, QLD)
rad_avg89 <- mean(r89)
writeRaster(rad_avg89, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad89_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad89_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r89avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r89 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r89avg_reproj.tif')

r89seq <- crop(r89, SEQ)
plot(r89seq) # Slight rotation issue but I am unsure of why

writeRaster(r89seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad89_seq.tif")







r90 <- rast('./00_Data/Environmental_data/Solar_radiation/1990.radiation.nc')
r90 <- crop(r90, QLD)
rad_avg90 <- mean(r90)
writeRaster(rad_avg90, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad90_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad90_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r90avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r90 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r90avg_reproj.tif')

r90seq <- crop(r90, SEQ)
plot(r90seq) # Slight rotation issue but I am unsure of why

writeRaster(r90seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad90_seq.tif")




r91 <- rast('./00_Data/Environmental_data/Solar_radiation/1991.radiation.nc')
r91 <- crop(r91, QLD)
rad_avg91 <- mean(r91)
writeRaster(rad_avg91, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad91_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad91_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r91avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r91 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r91avg_reproj.tif')

r91seq <- crop(r91, SEQ)
plot(r91seq) # Slight rotation issue but I am unsure of why

writeRaster(r91seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad91_seq.tif")




r92 <- rast('./00_Data/Environmental_data/Solar_radiation/1992.radiation.nc')
r92 <- crop(r92, QLD)
rad_avg92 <- mean(r92)
writeRaster(rad_avg92, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad92_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad92_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r92avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r92 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r92avg_reproj.tif')

r92seq <- crop(r92, SEQ)
plot(r92seq) # Slight rotation issue but I am unsure of why

writeRaster(r92seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad92_seq.tif")




r93 <- rast('./00_Data/Environmental_data/Solar_radiation/1993.radiation.nc')
r93 <- crop(r93, QLD)
rad_avg93 <- mean(r93)
writeRaster(rad_avg93, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad93_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad93_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r93avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r93 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r93avg_reproj.tif')

r93seq <- crop(r93, SEQ)
plot(r93seq) # Slight rotation issue but I am unsure of why

writeRaster(r93seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad93_seq.tif")




r94 <- rast('./00_Data/Environmental_data/Solar_radiation/1994.radiation.nc')
r94 <- crop(r94, QLD)
rad_avg94 <- mean(r94)
writeRaster(rad_avg94, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad94_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad94_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r94avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))


r94 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r94avg_reproj.tif')

r94seq <- crop(r94, SEQ)
plot(r94seq) # Slight rotation issue but I am unsure of why

writeRaster(r94seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad94_seq.tif")





r95 <- rast('./00_Data/Environmental_data/Solar_radiation/1995.radiation.nc')
r95 <- crop(r95, QLD)
rad_avg95 <- mean(r95)
writeRaster(rad_avg95, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad95_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad95_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r95avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r95 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r95avg_reproj.tif')

r95seq <- crop(r95, SEQ)
plot(r95seq) # Slight rotation issue but I am unsure of why

writeRaster(r95seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad95_seq.tif")




r96 <- rast('./00_Data/Environmental_data/Solar_radiation/1996.radiation.nc')
r96 <- crop(r96, QLD)
rad_avg96 <- mean(r96)
writeRaster(rad_avg96, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad96_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad96_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r96avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r96 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r96avg_reproj.tif')

r96seq <- crop(r96, SEQ)
plot(r96seq) # Slight rotation issue but I am unsure of why

writeRaster(r96seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad96_seq.tif")




r97 <- rast('./00_Data/Environmental_data/Solar_radiation/1997.radiation.nc')
r97 <- crop(r97, QLD)
rad_avg97 <- mean(r97)
writeRaster(rad_avg97, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad97_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad97_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r97avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r97 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r97avg_reproj.tif')

r97seq <- crop(r97, SEQ)
plot(r97seq) # Slight rotation issue but I am unsure of why

writeRaster(r97seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad97_seq.tif")




r98 <- rast('./00_Data/Environmental_data/Solar_radiation/1998.radiation.nc')
r98 <- crop(r98, QLD)
rad_avg98 <- mean(r98)
writeRaster(rad_avg98, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad98_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad98_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r98avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r98 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r98avg_reproj.tif')

r98seq <- crop(r98, SEQ)
plot(r98seq) # Slight rotation issue but I am unsure of why

writeRaster(r98seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad98_seq.tif")





r99 <- rast('./00_Data/Environmental_data/Solar_radiation/1999.radiation.nc')
r99 <- crop(r99, QLD)
rad_avg99 <- mean(r99)
writeRaster(rad_avg99, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad99_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad99_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r99avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r99 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r99avg_reproj.tif')

r99seq <- crop(r99, SEQ)
plot(r99seq) # Slight rotation issue but I am unsure of why

writeRaster(r99seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad99_seq.tif")




r00 <- rast('./00_Data/Environmental_data/Solar_radiation/2000.radiation.nc')
r00 <- crop(r00, QLD)
rad_avg00 <- mean(r00)
writeRaster(rad_avg00, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad00_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad00_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r00avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r00 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r00avg_reproj.tif')

r00seq <- crop(r00, SEQ)
plot(r00seq) # Slight rotation issue but I am unsure of why

writeRaster(r00seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad00_seq.tif")




r01 <- rast('./00_Data/Environmental_data/Solar_radiation/2001.radiation.nc')
r01 <- crop(r01, QLD)
rad_avg01 <- mean(r01)
writeRaster(rad_avg01, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad01_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad01_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r01avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r01 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r01avg_reproj.tif')

r01seq <- crop(r01, SEQ)
plot(r01seq) # Slight rotation issue but I am unsure of why

writeRaster(r01seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad01_seq.tif")



r02 <- rast('./00_Data/Environmental_data/Solar_radiation/2002.radiation.nc')
r02 <- crop(r02, QLD)
rad_avg02 <- mean(r02)
writeRaster(rad_avg02, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad02_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad02_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r02avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r02 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r02avg_reproj.tif')

r02seq <- crop(r02, SEQ)
plot(r02seq) # Slight rotation issue but I am unsure of why

writeRaster(r02seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad02_seq.tif")





r03 <- rast('./00_Data/Environmental_data/Solar_radiation/2003.radiation.nc')
r03 <- crop(r03, QLD)
rad_avg03 <- mean(r03)
writeRaster(rad_avg03, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad03_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad03_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r03avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))


r03 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r03avg_reproj.tif')

r03seq <- crop(r03, SEQ)
plot(r03seq) # Slight rotation issue but I am unsure of why

writeRaster(r03seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad03_seq.tif")




r04 <- rast('./00_Data/Environmental_data/Solar_radiation/2004.radiation.nc')
r04 <- crop(r04, QLD)
rad_avg04 <- mean(r04)
writeRaster(rad_avg04, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad04_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad04_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r04avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))


r04 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r04avg_reproj.tif')

r04seq <- crop(r04, SEQ)
plot(r04seq) # Slight rotation issue but I am unsure of why

writeRaster(r04seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad04_seq.tif")




r05 <- rast('./00_Data/Environmental_data/Solar_radiation/2005.radiation.nc')
r05 <- crop(r05, QLD)
rad_avg05 <- mean(r05)
writeRaster(rad_avg05, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad05_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad05_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r05avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r05 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r05avg_reproj.tif')

r05seq <- crop(r05, SEQ)
plot(r05seq) # Slight rotation issue but I am unsure of why

writeRaster(r05seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad05_seq.tif")




r06 <- rast('./00_Data/Environmental_data/Solar_radiation/2006.radiation.nc')
r06 <- crop(r06, QLD)
rad_avg06 <- mean(r06)
writeRaster(rad_avg06, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad06_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad06_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r06avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r06 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r06avg_reproj.tif')

r06seq <- crop(r06, SEQ)
plot(r06seq) # Slight rotation issue but I am unsure of why

writeRaster(r06seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad06_seq.tif")




r07 <- rast('./00_Data/Environmental_data/Solar_radiation/2007.radiation.nc')
r07 <- crop(r07, QLD)
rad_avg07 <- mean(r07)
writeRaster(rad_avg07, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad07_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad07_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r07avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r07 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r07avg_reproj.tif')

r07seq <- crop(r07, SEQ)
plot(r07seq) # Slight rotation issue but I am unsure of why

writeRaster(r07seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad07_seq.tif")




r08 <- rast('./00_Data/Environmental_data/Solar_radiation/2008.radiation.nc')
r08 <- crop(r08, QLD)
rad_avg08 <- mean(r08)
writeRaster(rad_avg08, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad08_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad08_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r08avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r08 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r08avg_reproj.tif')

r08seq <- crop(r08, SEQ)
plot(r08seq) # Slight rotation issue but I am unsure of why

writeRaster(r08seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad08_seq.tif")




r09 <- rast('./00_Data/Environmental_data/Solar_radiation/2009.radiation.nc')
r09 <- crop(r09, QLD)
rad_avg09 <- mean(r09)
writeRaster(rad_avg09, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad09_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad09_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r09avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r09 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r09avg_reproj.tif')

r09seq <- crop(r09, SEQ)
plot(r09seq) # Slight rotation issue but I am unsure of why

writeRaster(r09seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad09_seq.tif")




r10 <- rast('./00_Data/Environmental_data/Solar_radiation/2010.radiation.nc')
r10 <- crop(r10, QLD)
rad_avg10 <- mean(r10)
writeRaster(rad_avg10, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad10_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad10_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r10avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))



r10 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r10avg_reproj.tif')

r10seq <- crop(r10, SEQ)
plot(r10seq) # Slight rotation issue but I am unsure of why

writeRaster(r10seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad10_seq.tif")




r11 <- rast('./00_Data/Environmental_data/Solar_radiation/2011.radiation.nc')
r11 <- crop(r11, QLD)
rad_avg11 <- mean(r11)
writeRaster(rad_avg11, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad11_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad11_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r11avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r11 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r11avg_reproj.tif')

r11seq <- crop(r11, SEQ)
plot(r11seq) # Slight rotation issue but I am unsure of why

writeRaster(r11seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad11_seq.tif")



r12 <- rast('./00_Data/Environmental_data/Solar_radiation/2012.radiation.nc')
r12 <- crop(r12, QLD)
rad_avg12 <- mean(r12)
writeRaster(rad_avg12, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad12_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad12_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r12avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))


r12 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r12avg_reproj.tif')

r12seq <- crop(r12, SEQ)
plot(r12seq) # Slight rotation issue but I am unsure of why

writeRaster(r12seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad12_seq.tif")




r13 <- rast('./00_Data/Environmental_data/Solar_radiation/2013.radiation.nc')
r13 <- crop(r13, QLD)
rad_avg13 <- mean(r13)
writeRaster(rad_avg13, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad13_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad13_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r13avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))


r13 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r13avg_reproj.tif')

r13seq <- crop(r13, SEQ)
plot(r13seq) # Slight rotation issue but I am unsure of why

writeRaster(r13seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad13_seq.tif")



r14 <- rast('./00_Data/Environmental_data/Solar_radiation/2014.radiation.nc')
r14 <- crop(r14, QLD)
rad_avg14 <- mean(r14)
writeRaster(rad_avg14, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad14_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad14_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r14avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r14 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r14avg_reproj.tif')

r14seq <- crop(r14, SEQ)
plot(r14seq) # Slight rotation issue but I am unsure of why

writeRaster(r14seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad14_seq.tif")




r15 <- rast('./00_Data/Environmental_data/Solar_radiation/2015.radiation.nc')
r15 <- crop(r15, QLD)
rad_avg15 <- mean(r15)
writeRaster(rad_avg15, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad15_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad15_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r15avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r15 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r15avg_reproj.tif')

r15seq <- crop(r15, SEQ)
plot(r15seq) # Slight rotation issue but I am unsure of why

writeRaster(r15seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad15_seq.tif")




r16 <- rast('./00_Data/Environmental_data/Solar_radiation/2016.radiation.nc')
r16 <- crop(r16, QLD)
rad_avg16 <- mean(r16)
writeRaster(rad_avg16, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad16_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad16_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r16avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r16 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r16avg_reproj.tif')

r16seq <- crop(r16, SEQ)
plot(r16seq) # Slight rotation issue but I am unsure of why

writeRaster(r16seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad16_seq.tif")




r17 <- rast('./00_Data/Environmental_data/Solar_radiation/2017.radiation.nc')
r17 <- crop(r17, QLD)
rad_avg17 <- mean(r17)
writeRaster(rad_avg17, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad17_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad17_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r17avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r17 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r17avg_reproj.tif')

r17seq <- crop(r17, SEQ)
plot(r17seq) # Slight rotation issue but I am unsure of why

writeRaster(r17seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad17_seq.tif")




r18 <- rast('./00_Data/Environmental_data/Solar_radiation/2018.radiation.nc')
r18 <- crop(r18, QLD)
rad_avg18 <- mean(r18)
writeRaster(rad_avg18, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad18_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad18_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r18avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r18 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r18avg_reproj.tif')

r18seq <- crop(r18, SEQ)
plot(r18seq) # Slight rotation issue but I am unsure of why

writeRaster(r18seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad18_seq.tif")





r19 <- rast('./00_Data/Environmental_data/Solar_radiation/2019.radiation.nc')
r19 <- crop(r19, QLD)
rad_avg19 <- mean(r19)
writeRaster(rad_avg19, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad19_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad19_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r19avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r19 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r19avg_reproj.tif')

r19seq <- crop(r19, SEQ)
plot(r19seq) # Slight rotation issue but I am unsure of why

writeRaster(r19seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad19_seq.tif")




r20 <- rast('./00_Data/Environmental_data/Solar_radiation/2020.radiation.nc')
r20 <- crop(r20, QLD)
rad_avg20 <- mean(r20)
writeRaster(rad_avg20, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad20_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad20_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r20avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r20 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r20avg_reproj.tif')

r20seq <- crop(r20, SEQ)
plot(r20seq) # Slight rotation issue but I am unsure of why

writeRaster(r20seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad20_seq.tif")




r21 <- rast('./00_Data/Environmental_data/Solar_radiation/2021.radiation.nc')
r21 <- crop(r21, QLD)
rad_avg21 <- mean(r21)
writeRaster(rad_avg21, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad21_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad21_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r21avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r21 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r21avg_reproj.tif')

r21seq <- crop(r21, SEQ)
plot(r21seq) # Slight rotation issue but I am unsure of why

writeRaster(r21seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad21_seq.tif")


r22 <- rast('./00_Data/Environmental_data/Solar_radiation/2022.radiation.nc')
r22 <- crop(r22, QLD)
rad_avg22 <- mean(r22)
writeRaster(rad_avg22, './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad22_avg.tif')

gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Averaged/rad22_avg.tif',
         dstfile = './00_Data/Environmental_data/Outputs/Solar_radiation/r22avg_reproj.tif',
         t_srs = 'EPSG:3577',
         tr = c(30,30))

r22 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/r22avg_reproj.tif')

r22seq <- crop(r22, SEQ)
plot(r22seq) # Slight rotation issue but I am unsure of why

writeRaster(r22seq, "./00_Data/Environmental_data/Outputs/Solar_radiation/rad22_seq.tif")



# 1.5.1 Read in the solar radiation data for SEQ and combine into one raster ----

r87 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad87_seq.tif')
r88 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad88_seq.tif')
r89 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad89_seq.tif')
r90 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad90_seq.tif')
r91 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad91_seq.tif')
r92 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad92_seq.tif')
r93 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad93_seq.tif')
r94 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad94_seq.tif')
r95 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad95_seq.tif')
r96 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad96_seq.tif')
r97 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad97_seq.tif')
r98 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad98_seq.tif')
r99 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad99_seq.tif')
r00 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad00_seq.tif')
r01 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad01_seq.tif')
r02 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad02_seq.tif')
r03 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad03_seq.tif')
r04 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad04_seq.tif')
r05 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad05_seq.tif')
r06 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad06_seq.tif')
r07 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad07_seq.tif')
r08 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad08_seq.tif')
r09 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad09_seq.tif')
r10 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad10_seq.tif')
r11 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad11_seq.tif')
r12 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad12_seq.tif')
r13 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad13_seq.tif')
r14 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad14_seq.tif')
r15 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad15_seq.tif')
r16 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad16_seq.tif')
r17 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad17_seq.tif')
r18 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad18_seq.tif')
r19 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad19_seq.tif')
r20 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad20_seq.tif')
r21 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad21_seq.tif')
r22 <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/rad22_seq.tif')

# As with the other environmental variables, to work with this data we need it to be averaged
solar_rad <- terra::mean(r87, r88, r89, r90, r91, r92, r93, r94, r95, r96, r97, r98, r99, r00, r01, r02, r03, r04, r05, r06, r07, r08, r09, r10, r11, r12, r13, r14, r15, r16, r17, r18, r19, r20, r21, r22)

plot(solar_rad)
solar_rad # Check this has worked



writeRaster(solar_rad, './00_Data/Environmental_data/Outputs/Solar_radiation/Solar_radiation_seq.tif', overwrite = T)





# 1.6 Soil nutrients using soil % clay ----
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


gdalwarp(srcfile = './00_Data/Environmental_data/Outputs/Solar_radiation/Solar_radiation_seq.tif',
         './00_Data/Environmental_data/Outputs/Solar_radiation/Solar_radiation_seq_cropped.tif',
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

solar_radiation <- rast('./00_Data/Environmental_data/Outputs/Solar_radiation/Solar_radiation_seq_cropped.tif')
names(solar_radiation) <- "average_solar_rad."

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

solar_rad_foc <- focal(solar_radiation, fun = "mean", na.policy = "only", na.rm = T)
range(unique(solar_radiation$average_solar_rad.))
solar_rad_foc
names()
names(solar_rad_foc) <- "Avg_solar_radiation"
writeRaster(solar_rad_foc, './00_Data/Environmental_data/Outputs/Solar_radiation/Solar_radiation_seq_cropped_focal.tif')

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
tempseason <- rast('E:/PhD/R_analysis/Fire_freq/00_Data/Environmental_data/Outputs/BioClim/tempseason_SEQ_cropped_focal.tif')    
precipseason <- rast('E:/PhD/R_analysis/Fire_freq/00_Data/Environmental_data/Outputs/BioClim/precipseason_SEQ_cropped_focal.tif')  
diurnal_temp <- rast('E:/PhD/R_analysis/Fire_freq/00_Data/Environmental_data/Outputs/BioClim/Diurnal_temp_meanSEQ_cropped_focal.tif')  
solar_radiation <- rast('E:/PhD/R_analysis/Fire_freq/00_Data/Environmental_data/Outputs/Solar_radiation/Solar_radiation_seq_cropped_focal.tif')  
FPC <- rast('E:/PhD/R_analysis/Fire_freq/00_Data/Environmental_data/Outputs/FPC/FPC_all_cropped_focal.tif')  
soil_clay <- rast('E:/PhD/R_analysis/Fire_freq/00_Data/Environmental_data/Outputs/Soil_clay/SEQ_soilclay_cropped_focal.tif')  
slope <- rast('E:/PhD/R_analysis/Fire_freq/00_Data/Environmental_data/Outputs/DEM/SEQ_slope_cropped_focal.tif')  
aspect <- rast('E:/PhD/R_analysis/Fire_freq/00_Data/Environmental_data/Outputs/DEM/SEQaspect_cropped_focal.tif')  
topo_position <- rast('E:/PhD/R_analysis/Fire_freq/00_Data/Environmental_data/Outputs/DEM/SEQ_TPI_cropped_focal.tif') 
elev <- rast('E:/PhD/R_analysis/Fire_freq/00_Data/Environmental_data/Outputs/DEM/SEQ_DEM_reproj_cropped_focal.tif')

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


solar_radiation <- mask(solar_radiation, canal, inverse = T)
solar_radiation <- mask(solar_radiation, lake, inverse = T)
solar_radiation <- mask(solar_radiation, pond, inverse = T)
solar_radiation <- mask(solar_radiation, reservoir, inverse = T)
solar_radiation <- mask(solar_radiation, watercourse, inverse = T)
solar_radiation <- mask(solar_radiation, coast)


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

names(solar_radiation) <- "Average_solar_rad"

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

predictors <- c(QPWS_SEQ_ff, TWI, tempseason, precipseason, diurnal_temp, solar_radiation, FPC, soil_clay, slope, aspect, topo_position, elev)
writeRaster(predictors, './00_Data/SDM_data/predictors.tif', overwrite = T)



