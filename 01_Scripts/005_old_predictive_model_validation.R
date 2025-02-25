# Written by Felicity Charles
# Date: 06/10/2024

##### Fire frequency analysis ----
# This script validates model predictions, checking correlations between the original dataset and predicted values 

# R version 4.3.1 

# 1. Load required packages ----

library(terra) # terra_1.7-78
library(tidyterra) # tidyterra_0.6.1
library(ggspatial) # ggspatial_1.1.9
library(dplyr) # dplyr_1.1.4
library(sf) # sf_1.0-14 
library(ggplot2) # ggplot2_3.5.1
library(cowplot) # cowplot_1.1.1
library(ggforce) # ggforce_0.4.2 

# 2. Load original data, predictive model data, and environmental data
unweighted_pred <- rast('./04_Results/Prediction_rasters/Unweighted_pred.tif')
down_wt_pred <- rast('./04_Results/Prediction_rasters/Downweighted_pred.tif')
IWLR_pred <- rast('./04_Results/Prediction_rasters/IWLR_pred.tif')
gam_pred <- rast('./04_Results/Prediction_rasters/GAM_pred.tif')
glm_pred <- rast('./04_Results/Prediction_rasters/GLM_pred.tif')
Sentinel_ff <- rast('./00_Data/Fire_data/Outputs/Sentinel/Sentinel_ff_hydrographical_mask_SEQ_focal_cropped.tif')


environmental_preds <- rast('./00_Data/SDM_data/predictors.tif')
BVG <- vect('./00_Data/Environmental_data/Remnant_2021_broad_veg_groups/Remnant_broad_vegetation_groups.shp') %>%
  project('EPSG:3577') %>% 
  crop(gam_pred)
RE <- vect('./00_Data/Environmental_data/Regional_ecosystem/Biodiversity_status_of_remnant_regional_ecosystems.shp') %>% 
  project('EPSG:3577') %>% 
  crop(gam_pred)

QPWS_rand <- vect('./00_Data/Fire_data/Outputs/QPWS_random.gpkg')
QPWS_ff <- rast('./00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_freq_hydrographical_mask_cropped_reproj.tif')
SEQ <- vect('./00_Data/Australia_shapefile/STE11aAust.shp') %>% 
  project('EPSG:3577') %>% 
  crop(gam_pred)
protected_land <- vect('./00_Data/Protected_areas/Protected_areas.shp') %>% 
  project('EPSG:3577') %>% 
  crop(Sentinel_ff)
HV <- vect('./00_Data/Spatial data/HV.kml') %>% 
  project('EPSG:3577')
Gillies <- vect('./00_Data/Spatial data/nature-refuge-gillies-ridge-nature-refuge.kml') %>% 
  project('EPSG:3577')
Bulimbah <- vect('./00_Data/Spatial data/nature-refuge-bulimbah-nature-refuge.kml') %>% 
  project('EPSG:3577')
Bartopia <- vect('./00_Data/Spatial data/nature-refuge-bartopia-nature-refuge.kml') %>% 
  project('EPSG:3577')

# 3. Validate model predictions ----
# We know when we look at the each raster the min and max value are incorrect as this is not what gets plotted on a map.

range(Sentinel_ff$focal_mean) # Maximum is 26
environmental_preds$QPWS_firefreq # Maximum is 12
gam_pred # Maximum is 15
down_wt_pred # Maximum is 72 so there is definitely over-estimation. We could undertake post-processing to remove this over-estimation
plot(down_wt_pred)
unweighted_pred # Maximum is 64, overestimating
range(unweighted_pred$lyr1)
summary(unweighted_pred$lyr1)
IWLR_pred # Maximum is 8, underestimating
glm_pred # Maximum is 16, underestimating

# So just looking at the maximum values, most model underpredict fire frequency. The GLM, GAM and down weighted BRT overpredict QPWS fire frequency but only the down-weighted BRT overpredicts Sentinel fire frequency as well. SO let's take a look at the distribution of the predictions to really see what is going on as overprediction may be limited to a few locations and relatively few fire frequencies.

# 3.1 Check correlation of predictive outputs with QPWS data -----
QPWS_ff_rand <- extract(environmental_preds$QPWS_firefreq, QPWS_rand)

# Original Sentinel data
Sentinel_rand <- extract(Sentinel_ff, QPWS_rand)
sent_cor <- cor.test(QPWS_ff_rand$QPWS_firefreq, Sentinel_rand$focal_mean) # Correlation = 0.2517431  


# Unweighted model
unweighted_rand <- extract(unweighted_pred, QPWS_rand)
unwt_cor <- cor.test(QPWS_ff_rand$QPWS_firefreq, unweighted_rand$lyr1) # Correlation = 0.3389155  
# Slight improvement of correlation between from Sentinel data


# Downweighted model
down_rand <- extract(down_wt_pred, QPWS_rand)
down_cor <- cor.test(QPWS_ff_rand$QPWS_firefreq, down_rand$lyr1) # Correlation = 0.3375106   


# IWLR weighted model
IWLR_rand <- extract(IWLR_pred, QPWS_rand)
IWLR_cor <- cor.test(QPWS_ff_rand$QPWS_firefreq, IWLR_rand$lyr1) # Correlation = 0.03899995    


# GAM 
gam_rand <- extract(gam_pred, QPWS_rand)
gam_cor <- cor.test(QPWS_ff_rand$QPWS_firefreq, gam_rand$lyr1) # correlation = 0.4738458    

# GLM
glm_rand <- extract(glm_pred, QPWS_rand)
glm_cor <- cor.test(QPWS_ff_rand$QPWS_firefreq, glm_rand$lyr1) # correlation = 0.7452386  


# The GLM has the highest correlation with QPWS fire frequency data, followed by the GAM and then the unweighted BRT. The down-weighted BRT which has the best model fit is only slightly worse than the unweighted BRT. We need to make some decisions about which is the best model to use. Let's take a look at the correlation with the Sentinel data as this will also guide our choices. 


# 3.2 Check correlation of predictive outputs with Sentinel data
set.seed(480)
SEQ_pts <- spatSample(Sentinel_ff, 10000, na.rm = T, xy = T)
SEQ_pts <- SEQ_pts[,1:2]

Sent_SEQ <- extract(Sentinel_ff, SEQ_pts)

# Unweighted model
unweighted_SEQ <- extract(unweighted_pred, SEQ_pts)
unweighted_SEQ <- round(unweighted_SEQ)
unwt_SEQ_cor <- cor.test(Sent_SEQ$focal_mean, unweighted_SEQ$lyr1) # Correlation = 0.326055  

# Down-weighted model
down_SEQ <- round(extract(down_wt_pred, SEQ_pts))
dwt_SEQ_cor <- cor.test(Sent_SEQ$focal_mean, down_SEQ$lyr1) # Correlation = 0.3205573  


# IWLR weighted model
IWLR_SEQ <- round(extract(IWLR_pred, SEQ_pts))
IWLR_SEQ_cor <- cor.test(Sent_SEQ$focal_mean, IWLR_SEQ$lyr1) # Correlation = 0.3082806   

# GAM model
GAM_SEQ <- round(extract(gam_pred, SEQ_pts))
GAM_SEQ_cor <- cor.test(Sent_SEQ$focal_mean, GAM_SEQ$lyr1) # Correlation = 0.316026   
GAM_SEQ[is.na(GAM_SEQ)] <- 0

# GLM model
GLM_SEQ <- round(extract(glm_pred, SEQ_pts))
GLM_SEQ_cor <- cor.test(Sent_SEQ$focal_mean, GLM_SEQ$lyr1) # Correlation = 0.2051413  
GLM_SEQ[is.na(GAM_SEQ)] <- 0

# While the GLM performed best for the QPWS data, the correlation with Sentinel data is worse than correlation between QPWS and Sentinel. This is not the best model to be used. The model with the highest correlation to Sentinel data is the GAM, this also had the highest correlation to the QPWS data. The down-weighted and unweighted BRT models correlation to Sentinel data is not much worse than the GAM so could feasibly be considered good models. In reality, when we consider the AUC ROC and precision-recall statistics, among the other model evaulation statistics, the GAM does not perform much worse than either of these BRT models. Our best model should be the model that has the highest correlation with the data that is already available to us, considering all of these aspects, we would say that the GAM model is actually our best model.



# 3.3 Compare predictions to QPWS and Sentinel data for locations which I have sampled -----
# 3.3.1 Transect based comparison
transects <- read.csv('C:/Users/s4590925/OneDrive - The University of Queensland/Desktop/GitHub/Fire_recruit/00_Data/Transect_location_data.csv', header = T, stringsAsFactors = T)
head(transects)
View(transects)

# Convert to spatial
transects_sf <- st_as_sf(transects, coords = c("Latitude", "Longitude"), crs = 'EPSG:4326')


transects_v <- vect(transects_sf) %>% 
  project('EPSG:3577')
plet(transects_v)


transects_rand <- extract(QPWS_ff, transects_v) # QPWS
transects_sent <- extract(Sentinel_ff, transects_v) 
transects_unwt <- extract(unweighted_pred, transects_v)
transects_dwt <- extract(down_wt_pred, transects_v)
transects_IWLR <- extract(IWLR_pred, transects_v)
transect_gam <- extract(gam_pred, transects_v)
transect_glm <- extract(glm_pred, transects_v)

transects_fire <- as.data.frame(transects$Transect)
transects_fire$Sentinel <- transects_sent$focal_mean
transects_fire$QPWS <- transects_rand$QPWS_SEQ_freq_raster
transects_fire$GAM <- transect_gam$lyr1
transects_fire$Downweighted <- transects_dwt$lyr1
transects_fire$Unweighted <- transects_unwt$lyr1
transects_fire$IWLR <- transects_IWLR$lyr1
transects_fire$GLM <- transect_glm$lyr1


transects_fire$GAM <- round(transect_gam$lyr1)
transects_fire$Downweighted <- round(transects_dwt$lyr1)
transects_fire$Unweighted <- round(transects_unwt$lyr1)
transects_fire$IWLR <- round(transects_IWLR$lyr1)
transects_fire$GLM <- round(transect_glm$lyr1)


View(transects_fire)

# All models tend to underpredict for QPWS but the GAM produces predictions that are more close to those for Sentinel where there is no QPWS data compared to the down-weighted model.


# 3.3.2 Genetics based comparison ----
genetics <- read.table('C:/Users/s4590925/OneDrive - The University of Queensland/Desktop/GitHub/Allocasuarina__genfire/00_Data/Individual_site_data.txt', header = T, stringsAsFactors = T)
head(genetics)

# Convert to spatial
genetics_sf <- st_as_sf(genetics, coords = c("Latitude", "Longitude"), crs = 'EPSG:4326')

genetics_v <- vect(genetics_sf) %>% 
  project('EPSG:3577')
plet(genetics_v)


genetics_qpws <- extract(QPWS_ff, genetics_v)
genetics_sent <- extract(Sentinel_ff, genetics_v)
genetics_unwt <- extract(unweighted_pred, genetics_v)
genetics_dwt <- extract(down_wt_pred, genetics_v)
genetics_IWLR <- extract(IWLR_pred, genetics_v)
genetics_gam <- extract(gam_pred, genetics_v)
genetics_glm <- extract(glm_pred, genetics_v)

genetics_fire <- as.data.frame(genetics$ind)
genetics_fire$Sentinel <- genetics_sent$focal_mean
genetics_fire$QPWS <- genetics_qpws$QPWS_SEQ_freq_raster
genetics_fire$GAM <- round(genetics_gam$lyr1)
genetics_fire$Downweighted <- round(genetics_dwt$lyr1)
genetics_fire$Unweighted <- round(genetics_unwt$lyr1)
genetics_fire$IWLR <- round(genetics_IWLR$lyr1)
genetics_fire$glm <- round(genetics_glm$lyr1)

View(genetics_fire)

# As above, the GAM is more likely to predict fire frequency that is closer to the value from Sentinel for areas of no QPWS data than our down-weighted BRT. 






# 4. Create histograms for the distribution of fire frequency on QPWS land which we will plot as inset map -----
# 4.1 Create maps to use with these histograms ----

unweighted <-ggplot() + 
  geom_spatvector(data = SEQ, fill = 'transparent', col = 'black')+
  geom_spatraster(data = unweighted_pred) +
  theme_cowplot(font_size = 12)+
  scale_fill_viridis_c(na.value = 'transparent', limits = c(1,30), breaks = c(1, 5,10,15,20,25,30), direction = 1) +
  labs(fill = 'Fire frequency') +
  annotation_scale(location = "bl", style = 'ticks', pad_y = unit(0.5, 'cm'), pad_x = unit(15, 'cm'), text_cex = 2)+
  annotation_north_arrow(location = "bl", which_north = T, height = unit(2, "cm"), width = unit(1.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(25, 'cm'), style = north_arrow_fancy_orienteering) +
  theme(legend.key.height = unit(2.5, 'cm'),
        legend.key.width = unit(1.75, 'cm'),
        legend.title = element_text(face = 'bold', size = 25),
        legend.text = element_text(size = 20))



downweighted <- ggplot() + 
  geom_spatvector(data = SEQ, fill = 'transparent', col = 'black')+
  geom_spatraster(data = down_wt_pred) +
  theme_cowplot(font_size = 12)+  
  scale_fill_viridis_c(na.value = 'transparent', limits = c(1,30), breaks = c(1, 5,10,15,20,25,30), direction = 1) +
  labs(fill = 'Fire frequency') +
  annotation_scale(location = "bl", style = 'ticks', pad_y = unit(0.5, 'cm'), pad_x = unit(15, 'cm'), text_cex = 2)+
  annotation_north_arrow(location = "bl", which_north = T, height = unit(2, "cm"), width = unit(1.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(25, 'cm'), style = north_arrow_fancy_orienteering) +
  theme(legend.key.height = unit(2.5, 'cm'),
        legend.key.width = unit(1.75, 'cm'),
        legend.title = element_text(face = 'bold', size = 25),
        legend.text = element_text(size = 20))




IWLR <- ggplot() + 
  geom_spatvector(data = SEQ, fill = 'transparent', col = 'black')+
  geom_spatraster(data = IWLR_pred) +
  theme_cowplot(font_size = 12)+  
  scale_fill_viridis_c(na.value = 'transparent', limits = c(1,30), breaks = c(1,5,10,15,20,25,30), direction = 1) +
  labs(fill = 'Fire frequency') +
  annotation_scale(location = "bl", style = 'ticks', pad_y = unit(0.5, 'cm'), pad_x = unit(15, 'cm'), text_cex = 2)+
  annotation_north_arrow(location = "bl", which_north = T, height = unit(2, "cm"), width = unit(1.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(25, 'cm'), style = north_arrow_fancy_orienteering) +
  theme(legend.key.height = unit(2.5, 'cm'),
        legend.key.width = unit(1.75, 'cm'),
        legend.title = element_text(face = 'bold', size = 25),
        legend.text = element_text(size = 20))


QPWS <- ggplot() +
  geom_spatvector(data = SEQ, fill = 'transparent', col = 'black')+
  geom_spatraster(data = environmental_preds$QPWS_firefreq) +
  theme_cowplot(font_size = 12)+  
  scale_fill_viridis_c(na.value = 'transparent', limits = c(1,30), breaks = c(1, 5,10,15,20,25,30), direction = 1) +
  labs(fill = 'Fire frequency')+
  annotation_scale(location = "bl", style = 'ticks', pad_y = unit(0.5, 'cm'), pad_x = unit(15, 'cm'), text_cex = 2)+
  annotation_north_arrow(location = "bl", which_north = T, height = unit(2, "cm"), width = unit(1.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(25, 'cm'), style = north_arrow_fancy_orienteering) +
  theme(legend.key.height = unit(2.5, 'cm'),
        legend.key.width = unit(1.75, 'cm'),
        legend.title = element_text(face = 'bold', size = 25),
        legend.text = element_text(size = 20))





GAM_m <- ggplot() +
  geom_spatvector(data = SEQ, fill = 'transparent', col = 'black')+
  geom_spatraster(data = gam_pred) +
  theme_cowplot(font_size = 12)+  
  scale_fill_viridis_c(na.value = 'transparent', limits = c(1,30), breaks = c(1,5,10,15,20,25,30), direction = 1) +
  labs(fill = 'Fire frequency')+
  annotation_scale(location = "bl", style = 'ticks', pad_y = unit(0.5, 'cm'), pad_x = unit(15, 'cm'), text_cex = 2)+
  annotation_north_arrow(location = "bl", which_north = T, height = unit(2, "cm"), width = unit(1.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(25, 'cm'), style = north_arrow_fancy_orienteering) +
  theme(legend.key.height = unit(1.75, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(face = 'bold', size = 25),
        legend.text = element_text(size = 20))






GLM_m <- ggplot() +
  geom_spatvector(data = SEQ, fill = 'transparent', col = 'black')+
  geom_spatraster(data = glm_pred) +
  theme_cowplot(font_size = 14)+  
  scale_fill_viridis_c(na.value = 'transparent', limits = c(1,30), breaks = c(1,5,10,15,20,25,30), direction = 1) +
  labs(fill = 'Fire frequency')+
  annotation_scale(location = "bl", style = 'ticks', pad_y = unit(0.5, 'cm'), pad_x = unit(15, 'cm'), text_cex = 2)+
  annotation_north_arrow(location = "bl", which_north = T, height = unit(2, "cm"), width = unit(1.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(25, 'cm'), style = north_arrow_fancy_orienteering) +
  theme(legend.key.height = unit(2.5, 'cm'),
        legend.key.width = unit(1.75, 'cm'),
        legend.title = element_text(face = 'bold', size = 25),
        legend.text = element_text(size = 20))





Sent <- ggplot() +
  geom_spatraster(data = environmental_preds$QPWS_firefreq) +
  geom_spatraster(data = Sentinel_ff) +
  theme_cowplot(font_size = 14)+  
  scale_fill_viridis_c(na.value = 'white', limits = c(1,30), breaks = c(1, 5,10,15,20,25,30), direction = 1) +
  geom_spatvector(data = SEQ, fill = 'transparent', col = 'black')+
  labs(fill = 'Fire frequency')+
  annotation_scale(location = "bl", style = 'ticks', pad_y = unit(0.5, 'cm'), pad_x = unit(15, 'cm'), text_cex = 2)+
  annotation_north_arrow(location = "bl", which_north = T, height = unit(2, "cm"), width = unit(1.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(25, 'cm'), style = north_arrow_fancy_orienteering) +
  theme(legend.key.height = unit(1.75, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(face = 'bold', size = 25),
        legend.text = element_text(size = 20))



Sent + theme(plot.background = element_rect(fill = 'gray93'), legend.key.width = unit(0, "cm"), legend.text = element_blank(), legend.title = element_blank())

# Produce plot with all maps

r_map_p <- plot_grid(GAM_m + theme(legend.position = 'none'), downweighted + theme(legend.position = "none"), unweighted+ theme(legend.position = "none"), IWLR+ theme(legend.position = "none"), nrow = 1, rel_widths = c(0.2,0.2,0.2,0.2))




# 4.2 Prepare data for histogram plots ----
QPWS_pres <- extract(QPWS_ff, QPWS_rand)
QPWS_pres[is.na(QPWS_pres)] <- 0
summary(QPWS_pres)
colnames(QPWS_pres) <- c('ID', 'lyr1')
QPWS_pres$Dataset <- 'QPWS'
QPWS_pres[is.na(QPWS_pres)] <- 0

Sent_pres <- extract(Sentinel_ff, QPWS_rand)
colnames(Sent_pres) <- c('ID', 'lyr1')
Sent_pres$lyr1 <- round(Sent_pres$lyr1)
Sent_pres$Dataset <- 'Sentinel'
head(Sent_pres)


uwt_pres <- extract(unweighted_pred, QPWS_rand)
uwt_pres <- as.data.frame(uwt_pres)
uwt_pres$ID <- 1:10000
colnames(uwt_pres) <- c('ID', 'lyr1')
uwt_pres$lyr1 <- round(uwt_pres$lyr1)
uwt_pres$Dataset <- 'Unweighted BRT'

dwt_pres <- extract(down_wt_pred, QPWS_rand)
dwt_pres <- as.data.frame(dwt_pres)
dwt_pres$ID <- 1:10000
colnames(dwt_pres) <- c('ID', 'lyr1')
dwt_pres$lyr1 <- round(dwt_pres$lyr1)
dwt_pres$Dataset <- 'Down-weighted BRT'
head(dwt_pres)

IWLR_pres <- extract(IWLR_pred, QPWS_rand)
IWLR_pres <- as.data.frame(IWLR_pres)
IWLR_pres$ID <- 1:10000
colnames(IWLR_pres) <- c('ID', 'lyr1')
IWLR_pres$lyr1 <- round(IWLR_pres$lyr1)
IWLR_pres$Dataset <- 'IWLR BRT'
head(IWLR_pres)

gam_pres <- extract(gam_pred, QPWS_rand)
gam_pres <- as.data.frame(gam_pres)
gam_pres$ID <- 1:10000
colnames(gam_pres) <- c('ID', 'lyr1')
gam_pres$lyr1 <- round(gam_pres$lyr1)
gam_pres$Dataset <- 'GAM'
head(gam_pres)

glm_pres <- extract(glm_pred, QPWS_rand)
glm_pres <- as.data.frame(glm_pres)
glm_pres$ID <- 1:10000
colnames(glm_pres) <- c('ID', 'lyr1')
glm_pres$lyr1 <- round(glm_pres$lyr1)
glm_pres$Dataset <- 'GLM'
head(glm_pres)


# Create the dataframes for histograms
Sent_pres <- rbind(QPWS_pres, Sent_pres)
Sent_pres$Dataset <- as.factor(Sent_pres$Dataset)
Sent_pres$Dataset <- relevel(Sent_pres$Dataset, "Sentinel")

uwt_pres <- rbind(uwt_pres, QPWS_pres)
uwt_pres$Dataset <- as.factor(uwt_pres$Dataset)
uwt_pres$Dataset <- relevel(uwt_pres$Dataset, "Unweighted BRT")

dwt_pres <- rbind(dwt_pres, QPWS_pres)
dwt_pres$Dataset <- as.factor(dwt_pres$Dataset)
dwt_pres$Dataset <- relevel(dwt_pres$Dataset, "Down-weighted BRT")
head(dwt_pres); tail(dwt_pres)

IWLR_pres <- rbind(IWLR_pres, QPWS_pres)
IWLR_pres$Dataset <- as.factor(IWLR_pres$Dataset)

gam_pres <- rbind(gam_pres,QPWS_pres)
gam_pres[is.na(gam_pres)] <- 0

glm_pres <- rbind(glm_pres, QPWS_pres)
glm_pres[is.na(glm_pres)] <- 0


# 4.3 Produce histograms ----




Sent_pres1 <- ggplot(Sent_pres[Sent_pres$lyr1 <=7,], aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1) +
  scale_fill_manual(values = c('steelblue', 'gray50'))+
  scale_x_continuous('Fire frequency', breaks = c(0,1,2,3,4,5,6,7))+
  scale_y_continuous('Count') +
  theme(legend.title = element_text(size = 15, face = 'bold'),
        axis.title = element_text(size = 15, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())


Sent_pres2 <- ggplot(Sent_pres[Sent_pres$lyr1 >7,], aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', bins = 9)+
  scale_fill_manual(values = c('steelblue', 'gray50'))+
  scale_x_continuous('Fire frequency', breaks = c(8,9,10,11,12,13,14,15,16))+
  scale_y_continuous(' Count', limits = c(0,50))+
  theme(legend.title = element_text(size = 15, face = 'bold'),
        axis.title = element_text(size = 15, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())




dwt_pres1 <- ggplot(dwt_pres[dwt_pres$lyr1<= 6, ], aes(x = lyr1, fill = Dataset))+
         geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('#579C97', 'gray50'))+
  scale_x_continuous('Fire frequency', breaks = c(0,1,2,3,4,5,6))+
  scale_y_continuous('Count', limits = c(0, 4500)) +
  theme(legend.title = element_text(size = 15, face = 'bold'),
        axis.title = element_text(size = 15, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())


dwt_pres2 <- ggplot(dwt_pres[dwt_pres$lyr1 > 6, ], aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('#579C97', 'gray50'))+
  scale_x_continuous('Fire frequency', breaks = c(7,8,9,10,11,12,13,14,15))+
  scale_y_continuous('Count', limits = c(0, 45)) +
  theme(legend.title = element_text(size = 15, face = 'bold'),
        axis.title = element_text(size = 15, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())





uwt_pres1 <- ggplot(uwt_pres[uwt_pres$lyr1 <= 6, ], aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('#2A6D7A', 'gray50'))+
  scale_x_continuous('Fire frequency', breaks = c(0,1,2,3,4,5,6))+
  scale_y_continuous('Count', limits = c(0, 4500)) +
  theme(legend.title = element_text(size = 15, face = 'bold'),
        axis.title = element_text(size = 15, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

uwt_pres2 <- ggplot(uwt_pres[uwt_pres$lyr1 > 6, ], aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('#2A6D7A', 'gray50'))+
  scale_x_continuous('Fire frequency', breaks = c(7,8,9,10,11,12,13,14,15,16))+
  scale_y_continuous('Count', limits = c(0, 55), breaks = seq(from = 0, to = 50, by = 10)) +
  theme(legend.title = element_text(size = 15, face = 'bold'),
        axis.title = element_text(size = 15, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())




IWLR_pres1 <- ggplot(IWLR_pres[IWLR_pres$lyr1 <= 6, ], aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('#8FCCB4', 'gray50'))+
  scale_x_continuous('Fire frequency', breaks = c(0,1,2,3,4,5,6))+
  scale_y_continuous('Count', limits = c(0, 7000)) +
  theme(legend.title = element_text(size = 15, face = 'bold'),
        axis.title = element_text(size = 15, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())




gam_pres1 <- ggplot(gam_pres[gam_pres$lyr1 <=6, ], aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c("#AAA970", 'gray50'))+
  scale_x_continuous('Fire frequency', breaks = c(0,1,2,3,4,5,6))+
  scale_y_continuous('Count') +
  theme(legend.title = element_text(size = 15, face = 'bold'),
        axis.title = element_text(size = 15, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())


gam_pres2 <- ggplot(gam_pres[gam_pres$lyr1 >6, ], aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c("#AAA970", 'gray50'))+
  scale_x_continuous('Fire frequency', breaks = seq(7,11,1))+
  scale_y_continuous('Count', limits = c(0,40)) +
  theme(legend.title = element_text(size = 15, face = 'bold'),
        axis.title = element_text(size = 15, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())





glm_pres1 <- ggplot(glm_pres[glm_pres$lyr1 <=6, ], aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c("#492050", 'gray50'))+
  scale_x_continuous('Fire frequency', breaks = c(0,1,2,3,4,5,6))+
  scale_y_continuous('Count', limits = c(0, 7000)) +
  theme(legend.title = element_text(size = 15, face = 'bold'),
        axis.title = element_text(size = 15, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())


glm_pres2 <- ggplot(glm_pres[glm_pres$lyr1 >6, ], aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c("#492050", 'gray50'))+
  scale_x_continuous('Fire frequency', breaks = seq(7,16,1))+
  scale_y_continuous('Count', limits = c(0,40)) +
  theme(legend.title = element_text(size = 15, face = 'bold'),
        axis.title = element_text(size = 15, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())





# New versions
glm_pres_plot <- 
  ggplot(glm_pres, aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c("#492050", 'gray30'),limits = c('GLM'))+
                      scale_x_continuous('Fire frequency', breaks = seq(0,16,1), limits = c(0,16))+
                      scale_y_continuous('Count', limits = c(0, 7000)) +
                      theme(legend.title = element_text(size = 20, face = 'bold'),
                            axis.title = element_text(size = 20, face = 'bold'),
                            axis.text = element_text(size = 15),
                            legend.text = element_text(size = 20),
                            legend.key.size = unit(1, 'cm'),
                            panel.grid.major.y = element_line(colour = 'gray80'),
                            panel.grid.minor.y = element_line(colour = 'gray80'),
                            axis.ticks.y = element_line(colour = 'gray80'),
                            panel.grid.major.x = element_blank(),
                            panel.grid.minor.x = element_blank(),
                            title = element_text(size = 20),
                            panel.background = element_blank(),
                            axis.line = element_line(colour = 'black'))+
                      labs(title = expression(bold("(b)"))) +
                      facet_zoom(ylim = c(0,80), zoom.size = 1)


Sent_pres_plot <- ggplot(Sent_pres, aes(x = lyr1, fill = Dataset))+
                      geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1) +
                      scale_fill_manual(values = c('steelblue', 'gray30'), labels = c('Sentinel', "Public land"))+
                      scale_x_continuous('Fire frequency', breaks = seq(0,16,1), limits = c(0,16))+
                      scale_y_continuous('Count', limits = c(0, 7000)) +
                      theme(legend.title = element_text(size = 20, face = 'bold'),
                            axis.title = element_text(size = 20, face = 'bold'),
                            axis.text = element_text(size = 15),
                            legend.text = element_text(size = 20),
                            legend.key.size = unit(1, 'cm'),
                            panel.grid.major.y = element_line(colour = 'gray80'),
                            panel.grid.minor.y = element_line(colour = 'gray80'),
                            axis.ticks.y = element_line(colour = 'gray80'),
                            panel.grid.major.x = element_blank(),
                            panel.grid.minor.x = element_blank(),
                            title = element_text(size = 20),
                            panel.background = element_blank(),
                            axis.line = element_line(colour = 'black'))+
                      labs(title = expression(bold("(a)"))) +
                      facet_zoom(ylim = c(0,80), zoom.size = 1)

gam_pres_plot <- ggplot(gam_pres, aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c("#AAA970", 'gray30'), limits = c('GAM'))+
  scale_x_continuous('Fire frequency', breaks = seq(0,16,1), limits = c(0,16))+
  scale_y_continuous('Count', limits = c(0, 7000)) +
  theme(legend.title = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size = 20, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.y = element_line(colour = 'gray80'),
        panel.grid.minor.y = element_line(colour = 'gray80'),
        axis.ticks.y = element_line(colour = 'gray80'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        title = element_text(size = 20),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black'))+
  labs(title = expression(bold("(c)"))) +
  facet_zoom(ylim = c(0,80), zoom.size = 1)




dwt_pres_plot <- ggplot(dwt_pres, aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('#579C97', 'gray30'), limits = ("Down-weighted BRT"))+
  scale_x_continuous('Fire frequency', breaks = seq(0,16,1), limits = c(0,16))+
  scale_y_continuous('Count', limits = c(0, 7000)) +
  theme(legend.title = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size = 20, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.y = element_line(colour = 'gray80'),
        panel.grid.minor.y = element_line(colour = 'gray80'),
        axis.ticks.y = element_line(colour = 'gray80'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        title = element_text(size = 20),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black'))+
  labs(title = expression(bold("(d)"))) +
  facet_zoom(ylim = c(0,80), zoom.size = 1)


uwt_pres_plot <- ggplot(uwt_pres, aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('#2A6D7A', 'gray30'),limits = c('Unweighted BRT'))+
  scale_x_continuous('Fire frequency', breaks = seq(0,16,1), limits = c(0,16))+
  scale_y_continuous('Count', limits = c(0, 7000)) +
  theme(legend.title = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size = 20, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.y = element_line(colour = 'gray80'),
        panel.grid.minor.y = element_line(colour = 'gray80'),
        axis.ticks.y = element_line(colour = 'gray80'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        title = element_text(size = 20),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black'))+
  labs(title = expression(bold("(e)"))) +
  facet_zoom(ylim = c(0,80), zoom.size = 1)

IWLR_pres_plot <-   ggplot(IWLR_pres[IWLR_pres$lyr1 <= 6, ], aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('#8FCCB4', 'gray30'), limits = ("IWLR BRT"))+
  scale_x_continuous('Fire frequency', breaks = seq(0,16,1), limits = c(-1,16))+
  scale_y_continuous('Count', limits = c(0, 7000)) +
  theme(legend.title = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size = 20, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.y = element_line(colour = 'gray80'),
        panel.grid.minor.y = element_line(colour = 'gray80'),
        axis.ticks.y = element_line(colour = 'gray80'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        title = element_text(size = 20),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black'))+
  labs(title = expression(bold("(f)"))) +
  facet_zoom(ylim = c(0, 80), zoom.size = 1)



# 4.4 Create plots for documents ----

p1 <- plot_grid(Sent_pres_plot + theme(axis.title.x = element_blank(), legend.position = "bottom"), glm_pres_plot + theme(legend.position = "bottom"), nrow = 2)
p1
#1500x1000



p2 <- plot_grid(gam_pres_plot + theme(legend.position = "bottom", axis.title.x = element_blank()), 
                dwt_pres_plot + theme(legend.position = "bottom", legend.title = element_blank(), axis.title = element_blank()), 
                uwt_pres_plot + theme(legend.position = "bottom"), 
                IWLR_pres_plot+ theme(legend.position = "bottom", legend.title = element_blank(), axis.title.y = element_blank()), nrow = 2)

p2
#2200x1000

unweighted_presence_plot <- plot_grid(unweighted, uwt_pres_plot,
                                      ncol = 1,
                                      nrow = 2,
                                      rel_heights = c(4,1))
# 1500x1100
downweighted_presence_plot <- plot_grid(downweighted, dwt_pres_plot,
                                        ncol = 1,
                                        nrow = 2,
                                        rel_heights = c(4,1))
# 1500x1100

IWLR_presence_plot <- plot_grid(IWLR, IWLR_pres_plot,
                                ncol = 1,
                                nrow = 2,
                                rel_heights = c(4,1))
# 1000x1100

GAM_presence_plot <- plot_grid(GAM_m, gam_pres_plot,
                               ncol = 1, nrow = 2, rel_heights = c(4,1))
#1200x1100

GLM_presence_plot <- plot_grid(GLM_m, glm_pres_plot,
                               ncol = 1, nrow = 2, rel_heights = c(4,1))
#1400x1100


Sent_qpws_plot <- plot_grid(Sent, Sent_pres_plot,
                            ncol = 1, nrow = 2, rel_heights = c(4,1))
#1400x1100


# 5. Examine model performance compared to Sentinel -----
# 5.1 Produce points outside of QPWS land -----
plet(Sentinel_ff)
Sentinel_ff_m <- mask(Sentinel_ff, QPWS_ff, inverse = T)
plet(Sentinel_ff_m)

#set.seed(480)
#SEQ_pts <- spatSample(Sentinel_ff, 10000, na.rm = T, xy = T)
#SEQ_pts <- SEQ_pts[,1:2]

Sentinel_bg <- extract(Sentinel_ff, SEQ_pts)
Sentinel_bg$Dataset <- "Sentinel"
colnames(Sentinel_bg) <- c("ID", "lyr1", "Dataset")
Sentinel_bg$lyr1 <- round(Sentinel_bg$lyr1)


Sent_bg <- extract(Sentinel_ff, SEQ_pts)
Sent_bg$Dataset <- "Sentinel"
colnames(Sent_bg) <- c("ID", "lyr1", "Dataset")
Sent_bg$lyr1 <- round(Sent_bg$lyr1)


QPWS_bg <- extract(QPWS_ff, SEQ_pts)
colnames(QPWS_bg) <- c("ID", "lyr1")
QPWS_bg$Dataset <- "QPWS"
QPWS_bg[is.na(QPWS_bg)] <- 0 

dwt_bg <- extract(down_wt_pred, SEQ_pts)
colnames(dwt_bg) <- c("ID", "lyr1")
dwt_bg$Dataset <- "Down-weighted BRT"
dwt_bg$lyr1 <- round(dwt_bg$lyr1)

uwt_bg <- extract(unweighted_pred, SEQ_pts)
colnames(uwt_bg) <- c("ID", "lyr1")
uwt_bg$Dataset <- 'Unweighted BRT'
uwt_bg$lyr1 <- round(uwt_bg$lyr1)

IWLR_bg <- extract(IWLR_pred, SEQ_pts)
colnames(IWLR_bg) <- c("ID", "lyr1")
IWLR_bg$Dataset <- "IWLR BRT"
IWLR_bg$lyr1<- round(IWLR_bg$lyr1)

gam_bg <- extract(gam_pred, SEQ_pts)
colnames(gam_bg) <- c("ID", "lyr1")
gam_bg$Dataset <- "GAM"
gam_bg$lyr1 <- round(gam_bg$lyr1)
gam_bg[is.na(gam_bg)] <- 0


glm_bg <- extract(glm_pred, SEQ_pts)
colnames(glm_bg) <- c('ID', 'lyr1')
glm_bg$Dataset <- "GLM"
glm_bg$lyr1 <- round(glm_bg$lyr1)
glm_bg[is.na(glm_bg)] <- 0

# 5.2 Calculate the correlations ----
dwt_bg_cor <- cor.test(Sentinel_bg$lyr1, dwt_bg$lyr1)
dwt_bg_cor # 0.3206938  

uwt_bg_cor <- cor.test(Sentinel_bg$lyr1, uwt_bg$lyr1)
uwt_bg_cor # 0.3253162 

IWLR_bg_cor <- cor.test(Sentinel_bg$lyr1, IWLR_bg$lyr1)
IWLR_bg_cor # 0.3073885 

gam_bg_cor <- cor.test(Sentinel_bg$lyr1, gam_bg$lyr1)
gam_bg_cor # 0.3549376  

glm_bg_cor <- cor.test(Sentinel_bg$lyr1, glm_bg$lyr1)
glm_bg_cor # 0.2946344  


qpws_sent_cor <- cor.test(Sent_bg$lyr1, QPWS_bg$lyr1)
qpws_sent_cor # 0.1163371 



# Combine the data
Sent_bg <- rbind(Sent_bg, QPWS_bg)
head(Sent_bg)
Sent_bg$Dataset <- as.factor(Sent_bg$Dataset)
Sent_bg$Dataset <- relevel(Sent_bg$Dataset, "Sentinel")


dwt_bg <- rbind(Sentinel_bg, dwt_bg)
dwt_bg$Dataset <- as.factor(dwt_bg$Dataset)
dwt_bg$Dataset <- relevel(dwt_bg$Dataset, "Sentinel")


uwt_bg <- rbind(Sentinel_bg, uwt_bg)
uwt_bg$Dataset <- as.factor(uwt_bg$Dataset)
uwt_bg$Dataset <- relevel(uwt_bg$Dataset, "Sentinel")


IWLR_bg <- rbind(Sentinel_bg, IWLR_bg)
IWLR_bg$Dataset <- as.factor(IWLR_bg$Dataset)
IWLR_bg$Dataset <- relevel(IWLR_bg$Dataset, "Sentinel")


gam_bg <- rbind(Sentinel_bg, gam_bg)
gam_bg$Dataset <- as.factor(gam_bg$Dataset)
gam_bg$Dataset <- relevel(gam_bg$Dataset, "Sentinel")


glm_bg <- rbind(Sentinel_bg, glm_bg)
glm_bg$Dataset <- as.factor(glm_bg$Dataset)
glm_bg$Dataset <- relevel(glm_bg$Dataset, "Sentinel")




# The correlation with Sentinel for any model is low but any of these predictions are better than the correlation to QPWS data.


# 5.3 Produce histograms ----


Sent_bg1 <- ggplot(Sent_bg[Sent_bg$lyr1 <=6,], aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1) +
  scale_fill_manual(values = c('steelblue', 'gray50'))+
  scale_x_continuous('Fire frequency', breaks = c(0,1,2,3,4,5,6))+
  scale_y_continuous('Count', breaks = seq(from = 0, to = 10000, by = 2000)) +
  theme(legend.title = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size = 20, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())


Sent_bg2 <- ggplot(Sent_bg[Sent_bg$lyr1 >6,], aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('steelblue', 'gray50'))+
  scale_x_continuous('Fire frequency', breaks = seq(from = 7, to = 21, by = 1))+
  scale_y_continuous(' Count', limits = c(0,35))+
  theme(legend.title = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size = 20, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())




dwt_bg_p <- ggplot(dwt_bg, aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('steelblue','#579C97'))+
  scale_x_continuous('Fire frequency', breaks = seq(0, 21,1))+
  scale_y_continuous('Count') +
  theme(legend.title = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size = 20, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        title = element_text(size = 20))+
  labs(title = expression(bold("(b)"))) +
  facet_zoom(ylim = c(0, 35), zoom.size = 1)
  
dwt_bg1 <- ggplot(dwt_bg[dwt_bg$lyr1<= 4, ], aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('steelblue','#579C97'))+
  scale_x_continuous('Fire frequency', breaks = c(0,1,2,3,4))+
  scale_y_continuous('Count') +
  theme(legend.title = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size = 20, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())


dwt_bg2 <- ggplot(dwt_bg[dwt_bg$lyr1 > 4, ], aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('steelblue', '#579C97'))+
  scale_x_continuous('Fire frequency', breaks = seq(from = 5, to = 13, by = 1), limits = c(4.5, 13.5))+
  scale_y_continuous('Count') +
  theme(legend.title = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size = 20, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())






uwt_bg1 <- ggplot(uwt_bg[uwt_bg$lyr1 <= 4, ], aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('steelblue', '#2A6D7A'))+
  scale_x_continuous('Fire frequency', breaks = c(0,1,2,3,4))+
  scale_y_continuous('Count') +
  theme(legend.title = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size = 20, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

uwt_bg2 <- ggplot(uwt_bg[uwt_bg$lyr1 > 4, ], aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('steelblue', '#579C97'))+
  scale_x_continuous('Fire frequency', breaks = seq(from = 5, to = 13, by = 1), limits = c(4.5, 13.5))+
  scale_y_continuous('Count', limits = c(0,450)) +
  theme(legend.title = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size = 20, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())




IWLR_bg1 <- ggplot(IWLR_bg[IWLR_bg$lyr1 <= 3, ], aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('steelblue', '#8FCCB4'))+
  scale_x_continuous('Fire frequency', breaks = c(0,1,2,3))+
  scale_y_continuous('Count', breaks = seq(0, 20000, 1000)) +
  theme(legend.title = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size = 20, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())


IWLR_bg2 <- ggplot(IWLR_bg[IWLR_bg$lyr1 > 3, ], aes(x = lyr1, fill = Dataset))+
   geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
   scale_fill_manual(values = c('steelblue', '#8FCCB4'))+
   scale_x_continuous('Fire frequency', breaks = c(4,5,6,7,8), limits = c(3.5,8.5))+
   scale_y_continuous('Count', limits = c(0, 800)) +
  theme(legend.title = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size = 20, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
# Predictions up to 8 but only sampled up to 5




gam_bg1 <- ggplot(gam_bg[gam_bg$lyr1 <=3, ], aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('steelblue', "#AAA970"))+
  scale_x_continuous('Fire frequency', breaks = c(0,1,2,3))+
  scale_y_continuous('Count', breaks = seq(0, 15000, 1000)) +
  theme(legend.title = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size = 20, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

gam_bg2 <- ggplot(gam_bg[gam_bg$lyr1 >3, ], aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('steelblue', "#AAA970"))+
  scale_x_continuous('Fire frequency', limits = c(3.5, 15.5), breaks = seq(4,15,1))+
  scale_y_continuous('Count', limits = c(0,800)) +
  theme(legend.title = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size = 20, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())# Predictions up to 15 but only sampled up to 10





glm_bg1 <- ggplot(glm_bg[glm_bg$lyr1 <=2, ], aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('steelblue', "#492050"))+
  scale_x_continuous('Fire frequency')+
  scale_y_continuous('Count', limits = c(0, 7000)) +
  theme(legend.title = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size = , face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

glm_bg2 <- ggplot(glm_bg[glm_bg$lyr1 >2, ], aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('steelblue', "#492050"))+
  scale_x_continuous('Fire frequency', limits = c(2.5, 17.5), breaks = seq(3,17,1))+
  scale_y_continuous('Count') +
  theme(legend.title = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size = , face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

# Predicts up to 17 but only sampled up to 9



# New versions
glm_bg_p <- ggplot(glm_bg, aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c("#492050", 'steelblue'), limits = 'GLM')+
  scale_x_continuous('Fire frequency', breaks = seq(0,21,1))+
  scale_y_continuous('Count', limits = c(0,10000)) +
  theme(legend.title = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size = , face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.y = element_line(colour = 'gray80'),
        panel.grid.minor.y = element_line(colour = 'gray80'),
        axis.ticks.y = element_line(colour = 'gray80'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        title = element_text(size = 20),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black'))+
  labs(title = expression(bold("(b)"))) +
  facet_zoom(ylim = c(0, 200), zoom.size = 1)

Sent_bg_p <- ggplot(Sent_bg, aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1) +
  scale_fill_manual(values = c('steelblue', 'gray50'), labels = c('Sentinel', "Public land"))+
  scale_x_continuous('Fire frequency', breaks = seq(0,21,1))+
  scale_y_continuous('Count', limits = c(0,10000)) +
  theme(legend.title = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size = 20, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.y = element_line(colour = 'gray80'),
        panel.grid.minor.y = element_line(colour = 'gray80'),
        axis.ticks.y = element_line(colour = 'gray80'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        title = element_text(size = 20),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black'))+
  labs(title = expression(bold("(a)"))) +
  facet_zoom(ylim = c(0,200), zoom.size = 1)

gam_bg_p <- ggplot(gam_bg, aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c("#AAA970", 'steelblue'), limits = "GAM")+
  scale_x_continuous('Fire frequency', breaks = seq(0,21,1))+
  scale_y_continuous('Count', limits = c(0,10000)) +
  theme(legend.title = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size = 20, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.y = element_line(colour = 'gray80'),
        panel.grid.minor.y = element_line(colour = 'gray80'),
        axis.ticks.y = element_line(colour = 'gray80'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        title = element_text(size = 20),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black'))+
  labs(title = expression(bold("(c)"))) +
  facet_zoom(ylim = c(0, 200), zoom.size = 1)


dwt_bg_p <- ggplot(dwt_bg, aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('#579C97', 'steelblue'), limits = 'Down-weighted BRT')+
  scale_x_continuous('Fire frequency', breaks = seq(0, 21,1))+
  scale_y_continuous('Count', limits = c(0,10000)) +
  theme(legend.title = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size = 20, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.y = element_line(colour = 'gray80'),
        panel.grid.minor.y = element_line(colour = 'gray80'),
        axis.ticks.y = element_line(colour = 'gray80'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        title = element_text(size = 20),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black'))+
  labs(title = expression(bold("(d)"))) +
  facet_zoom(ylim = c(0, 200), zoom.size = 1)

uwt_bg_p <- ggplot(uwt_bg, aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('#2A6D7A', 'steelblue'), limits = "Unweighted BRT")+
  scale_x_continuous('Fire frequency', breaks = seq(0,21,1))+
  scale_y_continuous('Count', limits = c(0,10000)) +
  theme(legend.title = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size = 20, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.y = element_line(colour = 'gray80'),
        panel.grid.minor.y = element_line(colour = 'gray80'),
        axis.ticks.y = element_line(colour = 'gray80'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        title = element_text(size = 20),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black'))+
  labs(title = expression(bold("(e)"))) +
  facet_zoom(ylim = c(0,200), zoom.size = 1)

IWLR_bg_p <- ggplot(IWLR_bg, aes(x = lyr1, fill = Dataset))+
  geom_histogram(colour = '#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('#8FCCB4', 'steelblue'), limits = "IWLR BRT")+
  scale_x_continuous('Fire frequency', breaks = seq(0,21,1))+
  scale_y_continuous('Count', limits = c(0,10000)) +
  theme(legend.title = element_text(size = 20, face = 'bold'),
        axis.title = element_text(size = 20, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.y = element_line(colour = 'gray80'),
        panel.grid.minor.y = element_line(colour = 'gray80'),
        axis.ticks.y = element_line(colour = 'gray80'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        title = element_text(size = 20),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black'))+
  labs(title = expression(bold("(f)"))) +
  facet_zoom(ylim = c(0,200), zoom.size = 1)




p1_bg <- plot_grid(Sent_bg_p + theme(axis.title.x = element_blank(), legend.position = "bottom"), glm_bg_p + theme(legend.position = "bottom"), nrow = 2)
p1_bg
#1500x1200



p2_bg <- plot_grid(gam_bg_p + theme(legend.position = "bottom", axis.title.x = element_blank()), 
                dwt_bg_p + theme(legend.position = "bottom", legend.title = element_blank(), axis.title = element_blank()), 
                uwt_bg_p + theme(legend.position = "bottom"), 
                IWLR_bg_p+ theme(legend.position = "bottom", legend.title = element_blank(), axis.title.y = element_blank()), nrow = 2)

p2_bg
#3200X1200




# 5.4 Create plots for documents ----

unweighted_bg_plot <- plot_grid(unweighted, uwt_bg_p,
                                      ncol = 1,
                                      nrow = 2,
                                      rel_heights = c(4,1))
# 1900x1100

downweighted_bg_plot <- plot_grid(downweighted, dwt_bg_p,
                                        ncol = 1,
                                        nrow = 2,
                                        rel_heights = c(4,1))
# 1900x1100

IWLR_bg_plot <- plot_grid(IWLR, IWLR_bg_p,
                                ncol = 1,
                                nrow = 2,
                                rel_heights = c(4,1))
# 1900x1100

GAM_bg_plot <- plot_grid(GAM_m, gam_bg_p,
                               ncol = 1, nrow = 2, rel_heights = c(4,1))
# 1900x1100

GLM_bg_plot <- plot_grid(GLM_m, glm_bg_p,
                               ncol = 1, nrow = 2, rel_heights = c(4,1))
# 1900x1100


Sent_bg_plot <- plot_grid(Sent, Sent_bg_p,
                            ncol = 1, nrow = 2, rel_heights = c(4,1))
# 1900x1100



# Final maps -----

r_unwt_p <- plot_grid(uwt_pres_plot + theme(axis.title.x = element_blank()), uwt_bg_p,
                      nrow = 2, ncol = 1)
unweighted_plot <- plot_grid(unweighted, r_unwt_p, nrow = 2, ncol = 1, rel_heights = c(5,2))
  
# 1700x2200


r_dwt_p <- plot_grid(dwt_pres_plot + theme(axis.title.x = element_blank()), dwt_bg_p, nrow = 2, ncol = 1)
downweighted_plot <- plot_grid(downweighted, r_dwt_p, nrow = 2, ncol = 1, rel_heights = c(5,2))


r_iwlr_p <- plot_grid(IWLR_pres_plot + theme(axis.title.x = element_blank()), IWLR_bg_p, nrow = 2, ncol = 1)
IWLR_plot <- plot_grid(IWLR, r_iwlr_p, nrow = 2, ncol = 1, rel_heights = c(5,2))


r_gam_p <- plot_grid(gam_pres_plot + theme(axis.title.x = element_blank()), gam_bg_p, nrow = 2, ncol = 1)
gam_plot <- plot_grid(GAM_m, r_gam_p, nrow = 2, ncol = 1, rel_heights = c(5,2))


r_glm_p <- plot_grid(glm_pres_plot + theme(axis.title.x = element_blank()), glm_bg_p, nrow = 2, ncol = 1)
glm_plot <- plot_grid(GLM_m, r_glm_p, nrow = 2, ncol = 1, rel_heights = c(5,2))



r_Sent_p <- plot_grid(Sent_pres_plot + theme(axis.title.x = element_blank()), Sent_bg_p, nrow = 2, ncol = 1)
Sent_plot <- plot_grid(Sent, r_Sent_p, nrow = 2, ncol = 1, rel_heights = c(5,2))










