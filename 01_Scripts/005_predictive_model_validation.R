# Written by Felicity Charles
# Date: 06/10/2024
# Updated: 12/11/2025

##### Fire frequency analysis ----
# This script validates model predictions, checking correlations between the original dataset and predicted values 

# R version 4.3.1 

# 1. Load required packages ----

library(terra) # terra_1.7-78 updated to 1.8-70
library(tidyterra) # tidyterra_0.6.1 updated to 0.7.2
library(ggspatial) # ggspatial_1.1.9 updated to 1.1.10
library(dplyr) # dplyr_1.1.4
library(sf) # sf_1.0-14 updated to 1.0-21
library(ggplot2) # ggplot2_3.5.1 updated to 4.0.0
library(cowplot) # cowplot_1.1.1 updated to 1.2.0
library(ggforce) # ggforce_0.4.2 updated to 0.5.0

# 2. Load original data, predictive model data, and environmental data
unweighted_pred_fpc <- rast('./04_Results/Prediction_rasters/SEQ_IBRA_Unweighted_pred_FPC.tif')
down_wt_pred_fpc <- rast('./04_Results/Prediction_rasters/SEQ_IBRA_Downweighted_FPC_pred.tif')
IWLR_pred_fpc <- rast('./04_Results/Prediction_rasters/SEQ_IBRA_IWLR_FPC_pred.tif')
gam_pred_fpc <- rast('./04_Results/Prediction_rasters/SEQ_IBRA_GAM_FPC_pred.tif')
glm_pred_fpc <- rast('./04_Results/Prediction_rasters/SEQ_IBRA_GLM_FPC_pred.tif')

unweighted_pred_ndvi <- rast('./04_Results/Prediction_rasters/SEQ_IBRA_Unweighted_pred_NDVI.tif')
down_wt_pred_ndvi <- rast('./04_Results/Prediction_rasters/SEQ_IBRA_Downweighted_NDVI_pred.tif')
IWLR_pred_ndvi <- rast('./04_Results/Prediction_rasters/SEQ_IBRA_IWLR_NDVI_pred.tif')
gam_pred_ndvi <- rast('./04_Results/Prediction_rasters/SEQ_IBRA_GAM_NDVI_pred.tif')
glm_pred_ndvi <- rast('./04_Results/Prediction_rasters/SEQ_IBRA_GLM_NDVI_pred.tif')

Sentinel_ff <- rast('./00_Data/Fire_data/Outputs/Sentinel/Sentinel_ff_hydrographical_mask_SEQ_IBRA_focal_cropped.tif')
QPWS_rand <- vect('./00_Data/Fire_data/Outputs/QPWS_random_SEQ_IBRA.gpkg')
QPWS_ff <- rast('./00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_IBRA_freq_hydrographical_mask_cropped.tif')
SEQ_land <- vect('./00_Data/SEQ_bound/SEQ_IBRA.gpkg')
protected_land <- vect('./00_Data/Protected_areas/Protected_areas.shp') %>% 
  project('EPSG:3577') %>% 
  crop(SEQ_land)
plot(protected_land)

Sentinel_ff_m <- mask(Sentinel_ff, protected_land, inverse = T) %>% 
  mask(SEQ_land)
plet(Sentinel_ff_m)
QPWS <- ifel(QPWS_ff == 0, NA, QPWS_ff) # Extract the cells where QPWS_ff is not zero
SEQ_ff <- merge(QPWS, Sentinel_ff_m)
plot(SEQ_ff)


Sentinel_ff <- mask(Sentinel_ff, SEQ_land)
plet(Sentinel_ff)

# 3. Sentinel_ff# 3. Validate model predictions ----
# We know when we look at the each raster the min and max value are incorrect as this is not what gets plotted on a map.

range(Sentinel_ff$Sentinel_fire_freq) # Maximum is 29
range(QPWS_ff$QPWS_firefreq) # Maximum is 12

# FPC models
gam_pred_fpc # Maximum is 40
down_wt_pred_fpc; range(down_wt_pred_fpc); plot(down_wt_pred_fpc) # Maximum is 71 probably just some outliers as the majority is ~30 fires, so very limited overestimation if any
unweighted_pred_fpc; plot(unweighted_pred_fpc) # Maximum is 115, again outliers are present as plot shows ~55 fires, so also overestimating
IWLR_pred_fpc; plot(IWLR_pred_fpc) # Maximum is 9, underestimating
glm_pred_fpc # Maximum is 29

# NDVI models
gam_pred_ndvi; plot(gam_pred_ndvi) # Maximum is 18
down_wt_pred_ndvi; range(down_wt_pred_ndvi); plot(down_wt_pred_ndvi) # Maximum is 75 so some outliers as the majority is ~40 fires, so some overestimation
unweighted_pred_ndvi; plot(unweighted_pred_ndvi) # Maximum is 135, again outliers are present as plot shows ~55 fires, so also overestimating
IWLR_pred_ndvi; plot(IWLR_pred_ndvi) # Maximum is 9, underestimating
glm_pred_ndvi # Maximum is 22


# So just looking at the maximum values, most model result in under predictions of fire frequency. The FPC GLM has the best prediction, matching Sentinel, the FPC GAM is not bad but may be overestimating. The NDVI GAM and GLM are ok, some underestimation. SO let's take a look at the distribution of the predictions to really see what is going on as overprediction may be limited to a few locations and relatively few fire frequencies.

# 3.1 Check correlation of predictive outputs with QPWS data -----
QPWS_ff_rand <- extract(QPWS_ff, QPWS_rand)

# Original Sentinel data
Sentinel_rand <- extract(Sentinel_ff, QPWS_rand)
sent_cor <- cor.test(QPWS_ff_rand$QPWS_firefreq, Sentinel_rand$Sentinel_fire_freq) # Correlation = 0.3309308 

# 3.1.1 FPC models ----
# Unweighted model
unweighted_rand_fpc <- extract(unweighted_pred_fpc, QPWS_rand)
unwt_cor_fpc <- cor.test(QPWS_ff_rand$QPWS_firefreq, unweighted_rand_fpc$lyr1) # Correlation = 0.3767418
# Slight improvement of correlation between from Sentinel data


# Downweighted model
down_rand_fpc <- extract(down_wt_pred_fpc, QPWS_rand)
down_cor_fpc <- cor.test(QPWS_ff_rand$QPWS_firefreq, down_rand_fpc$lyr1) # Correlation = 0.4104306


# IWLR weighted model
IWLR_rand_fpc <- extract(IWLR_pred_fpc, QPWS_rand)
IWLR_cor_fpc <- cor.test(QPWS_ff_rand$QPWS_firefreq, IWLR_rand_fpc$lyr1) # Correlation = 0.03987862


# GAM 
gam_rand_fpc <- extract(gam_pred_fpc, QPWS_rand)
gam_cor_fpc <- cor.test(QPWS_ff_rand$QPWS_firefreq, gam_rand_fpc$lyr1) # correlation = 0.5260704

# GLM
glm_rand_fpc <- extract(glm_pred_fpc, QPWS_rand)
glm_cor_fpc <- cor.test(QPWS_ff_rand$QPWS_firefreq, glm_rand_fpc$lyr1) # correlation = 0.5763089


# The GLM has the highest correlation with QPWS fire frequency data, followed by the GAM for the models with FPC


# 3.1.2 NDVI models ----
# Unweighted model
unweighted_rand_ndvi <- extract(unweighted_pred_ndvi, QPWS_rand)
unwt_cor_ndvi <- cor.test(QPWS_ff_rand$QPWS_firefreq, unweighted_rand_ndvi$lyr1) # Correlation = 0.3799667
# Slight improvement of correlation between from Sentinel data


# Downweighted model
down_rand_ndvi <- extract(down_wt_pred_ndvi, QPWS_rand)
down_cor_ndvi <- cor.test(QPWS_ff_rand$QPWS_firefreq, down_rand_ndvi$lyr1) # Correlation = 0.4174479


# IWLR weighted model
IWLR_rand_ndvi <- extract(IWLR_pred_ndvi, QPWS_rand)
IWLR_cor_ndvi <- cor.test(QPWS_ff_rand$QPWS_firefreq, IWLR_rand_ndvi$lyr1) # Correlation = -0.009195731 


# GAM 
gam_rand_ndvi <- extract(gam_pred_ndvi, QPWS_rand)
gam_cor_ndvi <- cor.test(QPWS_ff_rand$QPWS_firefreq, gam_rand_ndvi$lyr1) # correlation = 0.58607289

# GLM
glm_rand_ndvi <- extract(glm_pred_ndvi, QPWS_rand)
glm_cor_ndvi <- cor.test(QPWS_ff_rand$QPWS_firefreq, glm_rand_ndvi$lyr1) # correlation = 0.7394897 

# The GLM has the highest correlation with QPWS fire frequency data for the NDVI models. All NDVI models, excluding the IWLR BRT have better correlations with QPWS data than the FPC models.


# 4. Create maps to use with these histograms ----
Sent <- ggplot() +
  geom_spatraster(data = Sentinel_ff) +
  theme_cowplot(font_size = 17)+  
  scale_fill_viridis_c(na.value = 'transparent', limits = c(1,30), breaks = c(1, 5, 10, 15, 20, 25, 30), direction = 1) +
  geom_spatvector(data = SEQ_land, fill = 'transparent', col = 'black')+
  labs(fill = 'Fire frequency')+
  annotation_scale(location = "bl", style = 'ticks', pad_y = unit(0.5, 'cm'), pad_x = unit(15, 'cm'), text_cex = 2)+
  annotation_north_arrow(location = "bl", which_north = T, height = unit(2, "cm"), width = unit(1.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(25, 'cm'), style = north_arrow_fancy_orienteering) +
  theme(legend.key.height = unit(1.75, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(face = 'bold', size = 25),
        legend.text = element_text(size = 20),
        legend.position = "none")


QPWS <- ggplot()+
  geom_spatraster(data = QPWS_ff) +
  theme_cowplot(font_size = 20)+  
  scale_fill_viridis_c(na.value = 'transparent', limits = c(1,30), breaks = c(1, 5, 10, 15, 20, 25, 30), direction = 1) +
  geom_spatvector(data = SEQ_land, fill = 'transparent', col = 'black')+
  labs(fill = 'Fire frequency')+
  annotation_north_arrow(location = "br", which_north = T, height = unit(2, "cm"), width = unit(1.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(-0.3, 'cm'), style = north_arrow_fancy_orienteering) +
  theme(legend.key.height = unit(2.5, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(face = 'bold', size = 25),
        legend.text = element_text(size = 20))
#1837 x 1200


Observed <- ggplot()+
  geom_spatraster(data = SEQ_ff) +
  theme_cowplot(font_size = 20)+  
  scale_fill_viridis_c(na.value = 'transparent', limits = c(1,30), breaks = c(1, 5, 10, 15, 20, 25, 30), direction = 1) +
  geom_spatvector(data = SEQ_land, fill = 'transparent', col = 'black')+
  labs(fill = 'Fire frequency')+
  annotation_north_arrow(location = "br", which_north = T, height = unit(2, "cm"), width = unit(1.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(-0.3, 'cm'), style = north_arrow_fancy_orienteering) +
  theme(legend.key.height = unit(2.5, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(face = 'bold', size = 25),
        legend.text = element_text(size = 20))

# 4.1 FPC maps ---- 
unweighted_pred_fpc <- mask(unweighted_pred_fpc, SEQ_land) %>% 
  round()
unweighted_pred_fpc <- ifel(unweighted_pred_fpc$lyr1 >30, 30, unweighted_pred_fpc$lyr1)
unweighted_pred_fpc <- ifel(unweighted_pred_fpc$lyr1 == 0, NA, unweighted_pred_fpc$lyr1)
unweighted_pred_fpc
unweighted_fpc <- ggplot() + 
  geom_spatvector(data = SEQ_land, fill = 'transparent', col = 'black')+
  geom_spatraster(data = unweighted_pred_fpc) +
  theme_cowplot(font_size = 17)+
  scale_fill_viridis_c(na.value = 'transparent', limits = c(1,30), breaks = c(1, 5, 10, 15, 20, 25, 30), direction = 1) +
  labs(fill = 'Fire frequency') +
  annotation_scale(location = "bl", style = 'ticks', pad_y = unit(0.5, 'cm'), pad_x = unit(15, 'cm'), text_cex = 2)+
  annotation_north_arrow(location = "bl", which_north = T, height = unit(2, "cm"), width = unit(1.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(25, 'cm'), style = north_arrow_fancy_orienteering) +
  theme(legend.key.height = unit(2.5, 'cm'),
        legend.key.width = unit(1.75, 'cm'),
        legend.title = element_text(face = 'bold', size = 25),
        legend.text = element_text(size = 20))


down_wt_pred_fpc <- mask(down_wt_pred_fpc, SEQ_land) %>% 
  round()
down_wt_pred_fpc <- ifel(down_wt_pred_fpc$lyr1 >30, 30, down_wt_pred_fpc$lyr1)
down_wt_pred_fpc <- ifel(down_wt_pred_fpc$lyr1 == 0, NA, down_wt_pred_fpc$lyr1)
downweighted_fpc <- ggplot() + 
  geom_spatvector(data = SEQ_land, fill = 'transparent', col = 'black')+
  geom_spatraster(data = down_wt_pred_fpc) +
  theme_cowplot(font_size = 17)+  
  scale_fill_viridis_c(na.value = 'transparent', limits = c(1,30), breaks = c(1, 5, 10, 15, 20, 25, 30), direction = 1) +
  labs(fill = 'Fire frequency') +
  annotation_scale(location = "bl", style = 'ticks', pad_y = unit(0.5, 'cm'), pad_x = unit(15, 'cm'), text_cex = 2)+
  annotation_north_arrow(location = "bl", which_north = T, height = unit(2, "cm"), width = unit(1.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(25, 'cm'), style = north_arrow_fancy_orienteering) +
  theme(legend.key.height = unit(2.5, 'cm'),
        legend.key.width = unit(1.75, 'cm'),
        legend.title = element_text(face = 'bold', size = 25),
        legend.text = element_text(size = 20))


IWLR_pred_fpc <- mask(IWLR_pred_fpc, SEQ_land) %>% 
  round()
IWLR_pred_fpc <- ifel(IWLR_pred_fpc$lyr1 >30, 30, IWLR_pred_fpc$lyr1)
IWLR_pred_fpc <- ifel(IWLR_pred_fpc$lyr1 == 0, NA, IWLR_pred_fpc$lyr1)
IWLR_fpc <- ggplot() + 
  geom_spatvector(data = SEQ_land, fill = 'transparent', col = 'black')+
  geom_spatraster(data = IWLR_pred_fpc) +
  theme_cowplot(font_size = 17)+  
  scale_fill_viridis_c(na.value = 'transparent', limits = c(0.91,30), breaks = c(1, 5, 10, 15, 20, 25, 30), direction = 1) +
  labs(fill = 'Fire frequency') + 
  annotation_scale(location = "bl", style = 'ticks', pad_y = unit(0.5, 'cm'), pad_x = unit(15, 'cm'), text_cex = 2)+
  annotation_north_arrow(location = "bl", which_north = T, height = unit(2, "cm"), width = unit(1.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(25, 'cm'), style = north_arrow_fancy_orienteering) +
  theme(legend.key.height = unit(2.5, 'cm'),
        legend.key.width = unit(1.75, 'cm'),
        legend.title = element_text(face = 'bold', size = 25),
        legend.text = element_text(size = 20))

gam_pred_fpc <- round(gam_pred_fpc)
gam_pred_fpc <- ifel(gam_pred_fpc$lyr1 >30, 30, gam_pred_fpc$lyr1)
gam_pred_fpc <- ifel(gam_pred_fpc$lyr1 == 0, NA, gam_pred_fpc$lyr1)
GAM_fpc <- ggplot() +
  geom_spatvector(data = SEQ_land, fill = 'transparent', col = 'black')+
  geom_spatraster(data = gam_pred_fpc) +
  theme_cowplot(font_size = 17)+  
  scale_fill_viridis_c(na.value = 'transparent', limits = c(1,30), breaks = c(1, 5, 10, 15, 20, 25, 30), direction = 1) +
  labs(fill = 'Fire frequency')+
  #annotation_north_arrow(location = "br", which_north = T, height = unit(2, "cm"), width = unit(1.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(-0.3, 'cm'), style = north_arrow_fancy_orienteering) +
  theme(legend.key.height = unit(2.5, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(face = 'bold', size = 25),
        legend.text = element_text(size = 20),
        legend.position = "none")
#1680 x 1100

glm_pred_fpc <- round(glm_pred_fpc)
glm_pred_fpc <- ifel(glm_pred_fpc$lyr1 >30, 30, glm_pred_fpc$lyr1)
glm_pred_fpc <- ifel(glm_pred_fpc$lyr1 == 0, NA, glm_pred_fpc$lyr1)
GLM_fpc <- ggplot() +
  geom_spatvector(data = SEQ_land, fill = 'transparent', col = 'black')+
  geom_spatraster(data = glm_pred_fpc) +
  theme_cowplot(font_size = 17)+  
  scale_fill_viridis_c(na.value = 'transparent', limits = c(1,30), breaks = c(1, 5, 10, 15, 20, 25, 30), direction = 1) +
  labs(fill = 'Fire frequency')+
  annotation_scale(location = "bl", style = 'ticks', pad_y = unit(0.5, 'cm'), pad_x = unit(15, 'cm'), text_cex = 2)+
  annotation_north_arrow(location = "bl", which_north = T, height = unit(2, "cm"), width = unit(1.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(25, 'cm'), style = north_arrow_fancy_orienteering) +
  theme(legend.key.height = unit(2.5, 'cm'),
        legend.key.width = unit(1.75, 'cm'),
        legend.title = element_text(face = 'bold', size = 25),
        legend.text = element_text(size = 20),
        legend.position = "none")


# Produce plot with all maps
r_map_p <- plot_grid(downweighted_fpc + theme(legend.position = "none"), unweighted_fpc + theme(legend.position = "none"), IWLR_fpc + theme(legend.position = "none"), nrow = 1, rel_widths = c(0.2,0.2,0.2))
# 1650 x 1075



# 4.1 NDVI maps ---- 
unweighted_pred_ndvi <- mask(unweighted_pred_ndvi, SEQ_land) %>% 
  round()
unweighted_pred_ndvi <- ifel(unweighted_pred_ndvi$lyr1 >30, 30, unweighted_pred_ndvi$lyr1)
unweighted_pred_ndvi <- ifel(unweighted_pred_ndvi$lyr1 == 0, NA, unweighted_pred_ndvi$lyr1)
unweighted_ndvi <- ggplot() + 
  geom_spatvector(data = SEQ_land, fill = 'transparent', col = 'black')+
  geom_spatraster(data = unweighted_pred_ndvi) +
  theme_cowplot(font_size = 17)+
  scale_fill_viridis_c(na.value = 'transparent', limits = c(1,30), breaks = c(1, 5, 10, 15, 20, 25, 30), direction = 1) +
  labs(fill = 'Fire frequency') +
  annotation_scale(location = "bl", style = 'ticks', pad_y = unit(0.5, 'cm'), pad_x = unit(15, 'cm'), text_cex = 2)+
  annotation_north_arrow(location = "bl", which_north = T, height = unit(2, "cm"), width = unit(1.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(25, 'cm'), style = north_arrow_fancy_orienteering) +
  theme(legend.key.height = unit(2.5, 'cm'),
        legend.key.width = unit(1.75, 'cm'),
        legend.title = element_text(face = 'bold', size = 25),
        legend.text = element_text(size = 20))


down_wt_pred_ndvi <- mask(down_wt_pred_ndvi, SEQ_land) %>% 
  round()
down_wt_pred_ndvi <- ifel(down_wt_pred_ndvi$lyr1 >30, 30, down_wt_pred_ndvi$lyr1)
down_wt_pred_ndvi <- ifel(down_wt_pred_ndvi$lyr1 == 0, NA, down_wt_pred_ndvi$lyr1)
downweighted_ndvi <- ggplot() + 
  geom_spatvector(data = SEQ_land, fill = 'transparent', col = 'black')+
  geom_spatraster(data = down_wt_pred_ndvi) +
  theme_cowplot(font_size = 17)+  
  scale_fill_viridis_c(na.value = 'transparent', limits = c(1,30), breaks = c(1, 5, 10, 15, 20, 25, 30), direction = 1) +
  labs(fill = 'Fire frequency') +
  annotation_scale(location = "bl", style = 'ticks', pad_y = unit(0.5, 'cm'), pad_x = unit(15, 'cm'), text_cex = 2)+
  annotation_north_arrow(location = "bl", which_north = T, height = unit(2, "cm"), width = unit(1.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(25, 'cm'), style = north_arrow_fancy_orienteering) +
  theme(legend.key.height = unit(2.5, 'cm'),
        legend.key.width = unit(1.75, 'cm'),
        legend.title = element_text(face = 'bold', size = 25),
        legend.text = element_text(size = 20))


IWLR_pred_ndvi <- mask(IWLR_pred_ndvi, SEQ_land) %>% 
  round()
IWLR_pred_ndvi <- ifel(IWLR_pred_ndvi$lyr1 >30, 30, IWLR_pred_ndvi$lyr1)
IWLR_pred_ndvi <- ifel(IWLR_pred_ndvi$lyr1 == 0, NA, IWLR_pred_ndvi$lyr1)
IWLR_ndvi <- ggplot() + 
  geom_spatvector(data = SEQ_land, fill = 'transparent', col = 'black')+
  geom_spatraster(data = IWLR_pred_ndvi) +
  theme_cowplot(font_size = 17)+  
  scale_fill_viridis_c(na.value = 'transparent', limits = c(1,30), breaks = c(1, 5, 10, 15, 20, 25, 30), direction = 1) +
  labs(fill = 'Fire frequency') + 
  annotation_scale(location = "bl", style = 'ticks', pad_y = unit(0.5, 'cm'), pad_x = unit(15, 'cm'), text_cex = 2)+
  annotation_north_arrow(location = "bl", which_north = T, height = unit(2, "cm"), width = unit(1.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(25, 'cm'), style = north_arrow_fancy_orienteering) +
  theme(legend.key.height = unit(2.5, 'cm'),
        legend.key.width = unit(1.75, 'cm'),
        legend.title = element_text(face = 'bold', size = 25),
        legend.text = element_text(size = 20))


gam_pred_ndvi <- round(gam_pred_ndvi)
gam_pred_ndvi <- ifel(gam_pred_ndvi$lyr1 >30, 30, gam_pred_ndvi$lyr1)
gam_pred_ndvi <- ifel(gam_pred_ndvi$lyr1 == 0, NA, gam_pred_ndvi$lyr1)
GAM_ndvi <- ggplot() +
  geom_spatvector(data = SEQ_land, fill = 'transparent', col = 'black')+
  geom_spatraster(data = gam_pred_ndvi) +
  theme_cowplot(font_size = 17)+  
  scale_fill_viridis_c(na.value = 'transparent', limits = c(1,30), breaks = c(1, 5, 10, 15, 20, 25, 30), direction = 1) +
  labs(fill = 'Fire frequency')+
  #annotation_north_arrow(location = "br", which_north = T, height = unit(2, "cm"), width = unit(1.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(-0.3, 'cm'), style = north_arrow_fancy_orienteering) +
  theme(legend.key.height = unit(2.5, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(face = 'bold', size = 25),
        legend.text = element_text(size = 20),
        legend.position = "none")
#1686 x 1100


glm_pred_ndvi <- round(glm_pred_ndvi)
glm_pred_ndvi <- ifel(glm_pred_ndvi$lyr1 >30, 30, glm_pred_ndvi$lyr1)
glm_pred_ndvi <- ifel(glm_pred_ndvi$lyr1 == 0, NA, glm_pred_ndvi$lyr1)
GLM_ndvi <- ggplot() +
  geom_spatvector(data = SEQ_land, fill = 'transparent', col = 'black')+
  geom_spatraster(data = glm_pred_ndvi) +
  theme_cowplot(font_size = 17)+  
  scale_fill_viridis_c(na.value = 'transparent', limits = c(1,30), breaks = c(1, 5, 10, 15, 20, 25, 30), direction = 1) +
  labs(fill = 'Fire frequency')+
  annotation_scale(location = "bl", style = 'ticks', pad_y = unit(0.5, 'cm'), pad_x = unit(15, 'cm'), text_cex = 2)+
  annotation_north_arrow(location = "bl", which_north = T, height = unit(2, "cm"), width = unit(1.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(25, 'cm'), style = north_arrow_fancy_orienteering) +
  theme(legend.key.height = unit(2.5, 'cm'),
        legend.key.width = unit(1.75, 'cm'),
        legend.title = element_text(face = 'bold', size = 25),
        legend.text = element_text(size = 20),
        legend.position = "none")


# Produce plot with all maps
r_map_p_ndvi <- plot_grid(downweighted_ndvi + theme(legend.position = "none"), unweighted_ndvi + theme(legend.position = "none"), IWLR_ndvi + theme(legend.position = "none"), nrow = 1, rel_widths = c(0.2,0.2,0.2))
# 1650 x 1075



# 4.2 Prepare data for histogram plots ----
QPWS_pres <- extract(QPWS_ff, QPWS_rand)
QPWS_pres[is.na(QPWS_pres)] <- 0
summary(QPWS_pres)
colnames(QPWS_pres) <- c('ID', 'lyr1')
QPWS_pres[is.na(QPWS_pres)] <- 0

Sent_pres <- extract(Sentinel_ff, QPWS_rand)
colnames(Sent_pres) <- c('ID', 'lyr1')
Sent_pres$Dataset <- 'Sentinel'
head(Sent_pres)

# 4.2.1 FPC data ----
uwt_fpc_pres <- extract(unweighted_pred_fpc, QPWS_rand)
uwt_fpc_pres <- as.data.frame(uwt_fpc_pres)
uwt_fpc_pres$ID <- 1:10000
colnames(uwt_fpc_pres) <- c('ID', 'lyr1')
head(uwt_fpc_pres)

dwt_fpc_pres <- extract(down_wt_pred_fpc, QPWS_rand)
dwt_fpc_pres <- as.data.frame(dwt_fpc_pres)
dwt_fpc_pres$ID <- 1:10000
colnames(dwt_fpc_pres) <- c('ID', 'lyr1')
head(dwt_fpc_pres)

IWLR_fpc_pres <- extract(IWLR_pred_fpc, QPWS_rand)
IWLR_fpc_pres <- as.data.frame(IWLR_fpc_pres)
IWLR_fpc_pres$ID <- 1:10000
colnames(IWLR_fpc_pres) <- c('ID', 'lyr1')
head(IWLR_fpc_pres)

gam_fpc_pres <- extract(gam_pred_fpc, QPWS_rand)
gam_fpc_pres <- as.data.frame(gam_fpc_pres)
gam_fpc_pres$ID <- 1:10000
colnames(gam_fpc_pres) <- c('ID', 'lyr1')
head(gam_fpc_pres)
gam_fpc_pres[is.na(gam_fpc_pres)] <- 0


glm_fpc_pres <- extract(glm_pred_fpc, QPWS_rand)
glm_fpc_pres <- as.data.frame(glm_fpc_pres)
glm_fpc_pres$ID <- 1:10000
colnames(glm_fpc_pres) <- c('ID', 'lyr1')
head(glm_fpc_pres)
glm_fpc_pres[is.na(glm_fpc_pres)] <- 0



# 4.2.2 NDVI data ----
uwt_ndvi_pres <- extract(unweighted_pred_ndvi, QPWS_rand)
uwt_ndvi_pres <- as.data.frame(uwt_ndvi_pres)
uwt_ndvi_pres$ID <- 1:10000
colnames(uwt_ndvi_pres) <- c('ID', 'lyr1')
head(uwt_ndvi_pres)

dwt_ndvi_pres <- extract(down_wt_pred_ndvi, QPWS_rand)
dwt_ndvi_pres <- as.data.frame(dwt_ndvi_pres)
dwt_ndvi_pres$ID <- 1:10000
colnames(dwt_ndvi_pres) <- c('ID', 'lyr1')
head(dwt_ndvi_pres)

IWLR_ndvi_pres <- extract(IWLR_pred_ndvi, QPWS_rand)
IWLR_ndvi_pres <- as.data.frame(IWLR_ndvi_pres)
IWLR_ndvi_pres$ID <- 1:10000
colnames(IWLR_ndvi_pres) <- c('ID', 'lyr1')
head(IWLR_ndvi_pres)

gam_ndvi_pres <- extract(gam_pred_ndvi, QPWS_rand)
gam_ndvi_pres <- as.data.frame(gam_ndvi_pres)
gam_ndvi_pres$ID <- 1:10000
colnames(gam_ndvi_pres) <- c('ID', 'lyr1')
head(gam_ndvi_pres)
gam_ndvi_pres[is.na(gam_ndvi_pres)] <- 0


glm_ndvi_pres <- extract(glm_pred_ndvi, QPWS_rand)
glm_ndvi_pres <- as.data.frame(glm_ndvi_pres)
glm_ndvi_pres$ID <- 1:10000
colnames(glm_ndvi_pres) <- c('ID', 'lyr1')
head(glm_ndvi_pres)
glm_ndvi_pres[is.na(glm_ndvi_pres)] <- 0


# 4.3 Produce histograms ----
# 4.3.1 Get RGB code for hexadecimal colours we want to use, note the col2rgb value needs to be divded by 255 to give the [0,1] required by rgb.
gb_c <- 'gray80'
col2rgb(gb_c, alpha = F)

gl_c <- "#492050"
col2rgb(gl_c, alpha = F)

s_c <- 'steelblue'
col2rgb(s_c, alpha = F)

ga_c <- '#AAA970'
col2rgb(ga_c, alpha = F)

uwt_c <- '#579C97'
col2rgb(uwt_c, alpha = F)

dwt_c <- '#2A6D7A'
col2rgb(dwt_c, alpha = F)

iwlr_c <- '#8FCCB4'
col2rgb(iwlr_c, alpha = F)



# 4.3.2 Produce histogram for each plot separately
gb_p <- hist(round(QPWS_pres$lyr1), breaks = seq(-1,16,1))
s_p <- hist(round(Sent_pres$lyr1), breaks = seq(-1,16,1))

# 4.3.2.1 FPC models
# Main plots 
gl_p_fpc <- hist(round(glm_fpc_pres$lyr1), breaks = seq(-1,16,1))
ga_p_fpc <- hist(round(gam_fpc_pres$lyr1),  breaks = seq(-1,16,1))

# Subplots
dwt_p_fpc <- hist(round(dwt_fpc_pres$lyr1[dwt_fpc_pres$lyr1 >=0 & dwt_fpc_pres$lyr1 <17]),  breaks = seq(-1,16,1))
uwt_p_fpc <- hist(round(uwt_fpc_pres$lyr1[uwt_fpc_pres$lyr1 >=0 & uwt_fpc_pres$lyr1 < 17]),  breaks = seq(-1,16,1))
iwlr_p_fpc <- hist(round(IWLR_fpc_pres$lyr1),  breaks = seq(-1,16,1))

# 4.3.2.1 NDVI models
# Main plots 
gl_p_ndvi <- hist(round(glm_ndvi_pres$lyr1), breaks = seq(-1,16,1))
ga_p_ndvi <- hist(round(gam_ndvi_pres$lyr1),  breaks = seq(-1,16,1))

# Subplots
dwt_p_ndvi <- hist(round(dwt_ndvi_pres$lyr1[dwt_ndvi_pres$lyr1 >=0 & dwt_ndvi_pres$lyr1 <17]),  breaks = seq(-1,16,1))
uwt_p_ndvi <- hist(round(uwt_ndvi_pres$lyr1[uwt_ndvi_pres$lyr1 >=0 & uwt_ndvi_pres$lyr1 < 17]),  breaks = seq(-1,16,1))
iwlr_p_ndvi <- hist(round(IWLR_ndvi_pres$lyr1),  breaks = seq(-1,16,1))




# 4.3.3 Create the subplots with the reduced y axis
gb_ps <- hist(round(QPWS_pres$lyr1), ylim = c(0,100), breaks = seq(-1,16,1))
s_ps <- hist(round(Sent_pres$lyr1), ylim = c(0,100), breaks = seq(-1,16,1))

# 4.3.3.1 FPC models 
# Main plots
gl_ps_fpc <- hist(round(glm_fpc_pres$lyr1), ylim = c(0,100), breaks = seq(-1,16,1))
ga_ps_fpc <- hist(round(gam_fpc_pres$lyr1), ylim = c(0,100), breaks = seq(-1,16,1))

# Subplots
dwt_ps_fpc <- hist(round(dwt_fpc_pres$lyr1[dwt_fpc_pres$lyr1 >=0 & dwt_fpc_pres$lyr1 <17]), ylim = c(0,100), breaks = seq(-1,16,1))
uwt_ps_fpc <- hist(round(uwt_fpc_pres$lyr1[uwt_fpc_pres$lyr1 >=0 & uwt_fpc_pres$lyr1 < 17]), ylim = c(0,100), breaks = seq(-1,16,1))
iwlr_ps_fpc <- hist(round(IWLR_fpc_pres$lyr1), ylim = c(0,100), breaks = seq(-1,16,1))



# 4.3.3.1 NDVI models 
# Main plots
gl_ps_ndvi <- hist(round(glm_ndvi_pres$lyr1), ylim = c(0,100), breaks = seq(-1,16,1))
ga_ps_ndvi <- hist(round(gam_ndvi_pres$lyr1), ylim = c(0,100), breaks = seq(-1,16,1))

# Subplots
dwt_ps_ndvi <- hist(round(dwt_ndvi_pres$lyr1[dwt_ndvi_pres$lyr1 >=0 & dwt_ndvi_pres$lyr1 <17]), ylim = c(0,100), breaks = seq(-1,16,1))
uwt_ps_ndvi <- hist(round(uwt_ndvi_pres$lyr1[uwt_ndvi_pres$lyr1 >=0 & uwt_ndvi_pres$lyr1 < 17]), ylim = c(0,100), breaks = seq(-1,16,1))
iwlr_ps_ndvi <- hist(round(IWLR_ndvi_pres$lyr1), ylim = c(0,100), breaks = seq(-1,16,1))



# 4.3.4 Combine histograms for each fire data grouping ----
# 4.3.4.1 FPC models
dev.new(width = 15, height = 10, res = 300, dpi = 80, noRStudioGD = T)
par(mfrow = c(3, 2), mar = c(4,7,4,2), oma = c(1,3,1,22), mgp = c(1, 1.5, 0)) # Check if this helps with the axis labels compared to the tick marks, may need to fiddle with the line for labels further as well
# Observed data
plot(gb_p, col = rgb(204/255, 204/255,204/255, 1), 
     ylab = "",  las = 1, ylim = c(0,7000), yaxt = "n",
     xlab = "", xlim = c(0,16), xaxt = "n", breaks = seq(0, 16, 1), 
     border = 'white', main = "", cex.axis = 2.4, cex.lab = 2.8)
plot(s_p, col = rgb(70/255, 130/255, 180/255, 0.4),  ylim = c(0,7000),
     ylab = expression(bold("")),  las = 1,
     xlab = "", xlim = c(0,16), xaxt = "n",  
     border = 'white',bty= 'l', main = "", 
     add = T)
axis(side = 1, at = seq(-1,16,1), cex.axis = 2.1, line = -0.4)
axis(side = 1, at = c(10,12,14,16), cex.axis = 2.1, line = -0.4)
axis(side = 2, at = seq(0, 7000, 1000), cex.axis = 2.1, line = 0.3, las = 1)
mtext(expression(bold("Count of cells")), side = 2, cex = 1.7, line = 6.5)
mtext(expression(bold("Fire frequency")), side = 1, line = 4, cex = 1.7)
mtext(expression(bold("(a) Observed")), side = 3, line = 0.7, at = 19, cex = 2.5)
mtext(expression(paste("Pearson's ", italic("r"), " = 0.331")), line = -2.5, at = 11, cex = 2)


plot(gb_ps, col = rgb(204/255, 204/255, 204/255, 1), 
     ylab = "", las = 1,  ylim = c(0,100), yaxt = "n",
     xlab = "", xlim = c(0,16), xaxt = "n", breaks = seq(0, 16, 1),
     border = 'white', main = "", cex.axis = 2.4, cex.lab = 2.8)
plot(s_ps, col = rgb(70/255, 130/255, 180/255, 0.5),
     ylab = "", las = 1,  ylim = c(0,100), yaxt = "n",
     xlab = "", xlim = c(0,16), xaxt = "n",
     border = 'white',  main = "", 
     add = T)
axis(side = 1, at = seq(-1,16,1), cex.axis = 2.1, line = -0.3)
axis(side = 1, at = c(10,12,14,16), cex.axis = 2.1, line = -0.3)
axis(side = 2, at = seq(0, 100, 20), cex.axis = 2.1, line = 0.4, las = 1 )
axis(side = 2, at = seq(0, 100, 10), labels = F, line = 0.4)
mtext(expression(bold("Fire frequency")), side = 1, line = 4, cex = 1.7)

#mtext(expression(bold("Density")), side = 2, cex = 2.2, line = 3)

# GLM data
plot(gb_p, col = rgb(204/255, 204/255, 204/255, 1), 
     ylab = "",  las = 1, ylim = c(0,7000), yaxt = "n",
     xlab = "", xlim = c(0,16),  xaxt = "n",
     border = 'white',bty= 'l', main = "", cex.axis = 2.4, cex.lab = 2.8)
plot(gl_p_fpc, col = rgb(73/255, 32/255, 80/255, 0.5),
     ylab = "",  las = 1, ylim = c(0,7000), yaxt = "n",
     xlab = "", xlim = c(0,16),  xaxt = "n",
     border = 'white', main = "", cex.axis = 2.4, cex.lab = 2.8,
     add = T)
axis(side = 1, at = seq(-1,16,1), cex.axis = 2.1, line = -0.4)
axis(side = 1, at = c(10,12,14,16), cex.axis = 2.1, line = -0.4)
axis(side = 2, at = seq(0, 7000, 1000), cex.axis = 2.1, line = 0.3, las = 1)
mtext(expression(bold("Count of cells")), side = 2, cex = 1.7, line = 6.5)
mtext(expression(bold("Fire frequency")), side = 1, line = 4, cex = 1.7)
mtext(expression(bold("(b) GLM")), side = 3, line = 0.1, at = 19, cex = 2.5)
mtext(expression(paste("Pearson's ", italic("r"), " = 0.577")), line = -2.5, at = 11, cex = 2)

plot(gb_ps, col = rgb(204/255, 204/255, 204/255, 1),
     ylab = "", las = 1,  ylim = c(0,100), yaxt = "n",
     xlab = "", xlim = c(0,16), xaxt = "n",
     border = 'white', main = "", cex.axis = 2.4, cex.lab = 2.8)
plot(gl_ps_fpc, col = rgb(73/255, 32/255, 80/255, 0.5),
     ylab = "", las = 1,  ylim = c(0,100), yaxt = "n",
     xlab = "", xlim = c(0,16), xaxt = "n",
     border = 'white', main = "", cex.lab = 2.8, 
     add = T)
axis(side = 1, at = seq(-1,16,1), cex.axis = 2.1, line = -0.3)
axis(side = 1, at = c(10,12,14,16), cex.axis = 2.1, line = -0.3)
axis(side = 2, at = seq(0, 100, 20), cex.axis = 2.1, line = 0.4, las = 1)
axis(side = 2, at = seq(0, 100, 10), labels = F, line = 0.4)
mtext(expression(bold("Fire frequency")), side = 1, line = 4, cex = 1.7)


# GAM
plot(gb_p, col = rgb(204/255, 204/255, 204/255, 1), 
     ylab = "",  las = 1, ylim = c(0,7000), yaxt = "n",
     xlab = "", xlim = c(0,16),  xaxt = "n",
     border = 'white',bty= 'l', main = "", cex.axis = 2.4, cex.lab = 2.8)
plot(ga_p_fpc, col = rgb(170/255, 179/255, 112/255, 0.5),
     ylab = "",  las = 1, ylim = c(0,7000), yaxt = "n",
     xlab = "", xlim = c(0,16),  xaxt = "n",
     border = 'white', main = "", cex.axis = 2.4, cex.lab = 2.8,
     add = T)
axis(side = 1, at = seq(-1,16,1), cex.axis = 2.1, line = -0.4)
axis(side = 1, at = c(10,12,14,16), cex.axis = 2.1, line = -0.4)
axis(side = 2, at = seq(0, 7000, 1000), cex.axis = 2.1, line = 0.3, las = 1)
mtext(expression(bold("Count of cells")), side = 2, cex = 1.7, line = 6.5)
mtext(expression(bold("Fire frequency")), side = 1, line = 4, cex = 1.7)
mtext(expression(bold("(c) GAM")), side = 3, line = 0.1, at = 19, cex = 2.5)
mtext(expression(paste("Pearson's ", italic("r"), " = 0.526")), line = -2.5, at = 11, cex = 2)

plot(gb_ps, col = rgb(204/255, 204/255, 204/255, 1),
     ylab = "", las = 1,  ylim = c(0,100), yaxt = "n",
     xlab = "", xlim = c(0,16), xaxt = "n",
     border = 'white', main = "", cex.axis = 2.4, cex.lab = 2.8)
plot(ga_ps_fpc, col = rgb(170/255, 179/255, 112/255, 0.5),
     ylab = "", las = 1,  ylim = c(0,100), yaxt = "n",
     xlab = "", xlim = c(0,16), xaxt = "n",
     border = 'white', main = "", cex.lab = 2.8, 
     add = T)
axis(side = 1, at = seq(-1,16,1), cex.axis = 2.1, line = -0.3)
axis(side = 1, at = c(10,12,14,16), cex.axis = 2.1, line = -0.3)
axis(side = 2, at = seq(0, 100, 20), cex.axis = 2.1, line = 0.4, las = 1)
axis(side = 2, at = seq(0, 100, 10), labels = F, line = 0.4)
mtext(expression(bold("Fire frequency")), side = 1, line = 4, cex = 1.7)

par(xpd = NA)  # Allow drawing outside the plot region
legend(x = 15.3, y = 455, 
       legend = c("Public land", "Satellite", "GLM", "GAM", "Down-weighted BRT", "Unweighted BRT", "Infinite BRT"),
       fill = c(rgb(204/255, 204/255, 204/255, 1),
                rgb(70/255, 130/255, 180/255, 0.5),
                rgb(73/255, 32/255, 80/255, 0.5),
                rgb(170/255, 179/255, 112/255, 0.5),
                rgb(42/255, 109/255, 122/255, 0.5),
                rgb(87/255, 156/255, 151/255, 0.5),
                rgb(143/255, 204/255, 180/255, 0.5)),
       border = "white",
       bty = "n",
       cex = 2.5)





### Subplots
dev.new(width = 20, height = 8, res = 300, dpi = 80, noRStudioGD = T)
par(mfrow = c(2, 4), mar = c(6,4,4,2), oma = c(1,3,1,0), mgp = c(1,1.5,0))


# DWT
plot(gb_p, col = rgb(204/255, 204/255, 204/255, 1), 
     ylab = "",  las = 1, ylim = c(0,7000), yaxt = "n",
     xlab = "", xlim = c(-1,16),  xaxt = "n",
     border = 'white',bty= 'l', main = "", cex.axis = 2.4, cex.lab = 2.8)
plot(dwt_p_fpc, col = rgb(42/255, 109/255, 122/255, 0.5), 
     ylab = "",  las = 1, ylim = c(0,7000), yaxt = "n",
     xlab = "", xlim = c(-1,16),  xaxt = "n",
     border = 'white', main = "", cex.axis = 2.4, cex.lab = 2.8,
     add = T)
axis(side = 1, at = c(0,5,10,16), cex.axis = 2.1, line = -0.4, tick = F)
axis(side = 1, at = seq(0,16, 1), labels = F, line = -0.4)
axis(side = 2, at = seq(0, 7000, 1000), cex.axis = 2.1, line = -0.8, las = 1)
axis(side = 1, at = c(-1, 16), labels = F, line = -0.4, lwd.ticks = 0)
mtext(expression(bold("Count of cells")), side = 2, cex = 1.7, line = 5)
mtext(expression(bold("Fire frequency")), side = 1, line = 3.5, cex = 1.7)
mtext(expression(bold("(d) Down-weighted BRT")), side = 3, line = 1, at = 19, cex = 2.5)
mtext(expression(paste("Pearson's ", italic("r"), " = 0.373")), line = -2.5, at = 9, cex = 2)

plot(gb_ps, col = rgb(204/255, 204/255, 204/255, 1),
     ylab = "", las = 1,  ylim = c(0,100), yaxt = "n",
     xlab = "", xlim = c(-1,16), xaxt = "n",
     border = 'white', main = "", cex.axis = 2.4, cex.lab = 2.8)
plot(dwt_ps_fpc, col = rgb(42/255, 109/255, 122/255, 0.5),
     ylab = "", las = 1,  ylim = c(0,100), yaxt = "n",
     xlab = "", xlim = c(-1,16), xaxt = "n",
     border = 'white', main = "", cex.lab = 2.8, 
     add = T)
axis(side = 1, at = c(0,5,10,16), cex.axis = 2.1, line = -0.4, tick = F)
axis(side = 1, at = seq(0,16, 1), labels = F, line = -0.4)
axis(side = 1, at = c(-1, 16), labels = F, line = -0.4, lwd.ticks = 0)
axis(side = 2, at = seq(0, 100, 20), cex.axis = 2.1, line = -0.8, las = 1)
mtext(expression(bold("Fire frequency")), side = 1, line = 3.5, cex = 1.7)



# UWT
plot(gb_p, col = rgb(204/255, 204/255, 204/255, 1), 
     ylab = "",  las = 1, ylim = c(0,7000), yaxt = "n",
     xlab = "", xlim = c(-1,16),  xaxt = "n",
     border = 'white',bty= 'l', main = "", cex.axis = 2.4, cex.lab = 2.8)
plot(uwt_p_fpc, col = rgb(87/255, 156/255, 151/255, 0.5),
     ylab = "",  las = 1, ylim = c(0,7000), yaxt = "n",
     xlab = "", xlim = c(-1,16),  xaxt = "n",
     border = 'white', main = "", cex.axis = 2.4, cex.lab = 2.8,
     add = T)
axis(side = 1, at = c(0,5,10,16), cex.axis = 2.1, line = -0.4, tick = F)
axis(side = 1, at = seq(0,16, 1), labels = F, line = -0.4)
axis(side = 2, at = seq(0, 7000, 1000), cex.axis = 2.1, line = -0.8, las = 1)
axis(side = 1, at = c(-1, 16), labels = F, line = -0.4, lwd.ticks = 0)
mtext(expression(bold("Count of cells")), side = 2, cex = 1.7, line = 5)
mtext(expression(bold("Fire frequency")), side = 1, line = 3.5, cex = 1.7)
mtext(expression(bold("(e) Unweighted BRT")), side = 3, line = 1, at = 19, cex = 2.4)
mtext(expression(paste("Pearson's ", italic("r"), " = 0.377")), line = -2.5, at = 9, cex = 2)

plot(gb_ps, col = rgb(204/255, 204/255, 204/255, 1),
     ylab = "", las = 1,  ylim = c(0,100), yaxt = "n",
     xlab = "", xlim = c(-1,16), xaxt = "n",
     border = 'white', main = "", cex.axis = 2.4, cex.lab = 2.8)
plot(uwt_ps_fpc, col = rgb(87/255, 156/255, 151/255, 0.5),
     ylab = "", las = 1,  ylim = c(0,100), yaxt = "n",
     xlab = "", xlim = c(-1,16), xaxt = "n",
     border = 'white', main = "", cex.lab = 2.8, 
     add = T)
axis(side = 1, at = c(0,5,10,16), cex.axis = 2.1, line = -0.4, tick = F)
axis(side = 1, at = seq(0,16, 1), labels = F, line = -0.4)
axis(side = 1, at = c(-1, 16), labels = F, line = -0.4, lwd.ticks = 0)
axis(side = 2, at = seq(0, 100, 20), cex.axis = 2.1, line = -0.8, las = 1)
mtext(expression(bold("Fire frequency")), side = 1, line = 3.5, cex = 1.7)




# IWLR
plot(gb_p, col = rgb(204/255, 204/255, 204/255, 1), 
     ylab = "",  las = 1, ylim = c(0,7000), yaxt = "n",
     xlab = "", xlim = c(-1,16),  xaxt = "n",
     border = 'white',bty= 'l', main = "", cex.axis = 2.4, cex.lab = 2.8)
plot(iwlr_p_fpc, col = rgb(143/255, 204/255, 180/255, 0.5),
     ylab = "",  las = 1, ylim = c(0,7000), yaxt = "n",
     xlab = "", xlim = c(-1,16),  xaxt = "n",
     border = 'white', main = "", cex.axis = 2.4, cex.lab = 2.8,
     add = T)
axis(side = 1, at = c(0,5,16), cex.axis = 2.1, line = -0.4, tick = F)
axis(side = 1, at = seq(0,16, 1), labels = F, line = -0.4)
axis(side = 2, at = seq(0, 7000, 1000), cex.axis = 2.1, line = -0.8, las = 1)
axis(side = 1, at = c(-1, 16), labels = F, line = -0.4, lwd.ticks = 0)
mtext(expression(bold("Count of cells")), side = 2, cex = 1.7, line = 5)
mtext(expression(bold("Fire frequency")), side = 1, line = 3.5, cex = 1.7)
mtext(expression(bold("(f) Infinite BRT")), side = 3, line = 1, at = 19, cex = 2.5)
mtext(expression(paste("Pearson's ", italic("r"), " = 0.040")), line = -2.5, at = 9, cex = 2)

plot(gb_ps, col = rgb(204/255, 204/255, 204/255, 1),
     ylab = "", las = 1,  ylim = c(0,100), yaxt = "n",
     xlab = "", xlim = c(-1,16), xaxt = "n",
     border = 'white', main = "", cex.axis = 2.4, cex.lab = 2.8)
plot(iwlr_ps_fpc, col = rgb(143/255, 204/255, 180/255, 0.5),
     ylab = "", las = 1,  ylim = c(0,100), yaxt = "n",
     xlab = "", xlim = c(-1,16), xaxt = "n",
     border = 'white', main = "", cex.lab = 2.8, 
     add = T)
axis(side = 1, at = c(0,5,10,16), cex.axis = 2.1, line = -0.4, tick = F)
axis(side = 1, at = seq(0,16, 1), labels = F, line = -0.4)
axis(side = 1, at = c(-1, 16), labels = F, line = -0.4, lwd.ticks = 0)
axis(side = 2, at = seq(0, 100, 20), cex.axis = 2.1, line = -0.8, las = 1)
mtext(expression(bold("Fire frequency")), side = 1, line = 3.5, cex = 1.7)

save.image('./02_Workspaces/005_predictive_model_validations.RData')
#load('./02_Workspaces/005_predictive_model_validations.RData')




# NOW NDVI models


## Lets also look for 18
# We also want to determine what % of cells we have taken out of maps for Satllite, Unweighted BRT, and Downweighted BRT that were above 18 fires.
plot(Sentinel_ff$Fire_freq)
unique(terra::cells(Sentinel_ff, c(19:26))) # Number of cells with ff above 18 = 12371
cells(Sentinel_ff) # Total number of cells = 69976534
# Number of cells excluded 
(12371/69976534)*100 # % of cells with Fire freq >18

cells(unweighted_pred) # Total cells =  70535296. Number unchanged by adding round
cells(round(unweighted_pred), c(19:31)) # 238 cells. Needed to add round due to the decimals in predictions
(238/70535296)*100 # % of cells with fire freq >18

cells(down_wt_pred) # Total cells = 70535296
cells(round(down_wt_pred), c(19:27)) #185  cells
(185/70535296)*100 # % of cells with fire freq >18
