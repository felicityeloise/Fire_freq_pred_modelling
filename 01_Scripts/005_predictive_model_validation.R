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


QPWS_rand <- vect('./00_Data/Fire_data/Outputs/QPWS_random.gpkg')
QPWS_ff <- rast('./00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_freq_hydrographical_mask_cropped_reproj.tif')
SEQ <- vect('./00_Data/Australia_shapefile/STE11aAust.shp') %>% 
  project('EPSG:3577') %>% 
  crop(gam_pred)

Sentinel_ff_m <- mask(Sentinel_ff, QPWS_ff, inverse = T)
plet(Sentinel_ff_m)


# Join polygons for estates with the same name
# This was done in ArcGIS by using the aggregate polygons function with an aggregation distance of 150m and aggregate field of NAMEABBREV.
protected_land <- vect('./00_Data/Protected_areas/Protected_areas_dissolved.shp') %>% 
  project('EPSG:3577') %>% 
  crop(SEQ)

plet(protected_land)


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
QPWS_ff_rand <- extract(QPWS_ff, QPWS_rand)

# Original Sentinel data
Sentinel_rand <- extract(Sentinel_ff, QPWS_rand)
sent_cor <- cor.test(QPWS_ff_rand$QPWS_SEQ_freq_raster, Sentinel_rand$focal_mean) # Correlation = 0.2517431  


# Unweighted model
unweighted_rand <- extract(unweighted_pred, QPWS_rand)
unwt_cor <- cor.test(QPWS_ff_rand$QPWS_SEQ_freq_raster, unweighted_rand$lyr1) # Correlation = 0.3379137  
# Slight improvement of correlation between from Sentinel data


# Downweighted model
down_rand <- extract(down_wt_pred, QPWS_rand)
down_cor <- cor.test(QPWS_ff_rand$QPWS_SEQ_freq_raster, down_rand$lyr1) # Correlation = 0.336561   


# IWLR weighted model
IWLR_rand <- extract(IWLR_pred, QPWS_rand)
IWLR_cor <- cor.test(QPWS_ff_rand$QPWS_SEQ_freq_raster, IWLR_rand$lyr1) # Correlation = 0.03820975 


# GAM 
gam_rand <- extract(gam_pred, QPWS_rand)
gam_cor <- cor.test(QPWS_ff_rand$QPWS_SEQ_freq_raster, gam_rand$lyr1) # correlation = 0.4726568     

# GLM
glm_rand <- extract(glm_pred, QPWS_rand)
glm_cor <- cor.test(QPWS_ff_rand$QPWS_SEQ_freq_raster, glm_rand$lyr1) # correlation = 0.7445872   


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


pred_ff <- mosaic(glm_pred, gam_pred, fun = "max")
pred_ff

transects_pred <- extract(pred_ff, transects_v)
transects_rand <- extract(QPWS_ff, transects_v) # QPWS
transects_fire$QPWS <- transects_rand$QPWS_SEQ_freq_raster
transects_fire$pred <- round(transects_pred$lyr1)



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






# 4. Create maps to use with these histograms ----
 

unweighted <- ggplot() + 
  geom_spatvector(data = SEQ, fill = 'transparent', col = 'black')+
  geom_spatraster(data = unweighted_pred) +
  theme_cowplot(font_size = 17)+
  scale_fill_viridis_c(na.value = 'transparent', limits = c(0.51,17), breaks = seq(1,18,1), direction = 1) +
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
  theme_cowplot(font_size = 17)+  
  scale_fill_viridis_c(na.value = 'transparent', limits = c(0.51,17), breaks = seq(1,18,1), direction = 1) +
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
  theme_cowplot(font_size = 17)+  
  scale_fill_viridis_c(na.value = 'transparent', limits = c(0.51,17), breaks = seq(1,18,1), direction = 1) +
  labs(fill = 'Fire frequency') + 
  annotation_scale(location = "bl", style = 'ticks', pad_y = unit(0.5, 'cm'), pad_x = unit(15, 'cm'), text_cex = 2)+
  annotation_north_arrow(location = "bl", which_north = T, height = unit(2, "cm"), width = unit(1.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(25, 'cm'), style = north_arrow_fancy_orienteering) +
  theme(legend.key.height = unit(2.5, 'cm'),
        legend.key.width = unit(1.75, 'cm'),
        legend.title = element_text(face = 'bold', size = 25),
        legend.text = element_text(size = 20))





GAM_m <- ggplot() +
  geom_spatvector(data = SEQ, fill = 'transparent', col = 'black')+
  geom_spatraster(data = gam_pred) +
  theme_cowplot(font_size = 17)+  
  scale_fill_viridis_c(na.value = 'transparent', limits = c(0.51,17), breaks = seq(1,18,1), direction = 1) +
  labs(fill = 'Fire frequency')+
  annotation_north_arrow(location = "br", which_north = T, height = unit(2, "cm"), width = unit(1.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(-0.3, 'cm'), style = north_arrow_fancy_orienteering) +
  theme(legend.key.height = unit(2.5, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(face = 'bold', size = 25),
        legend.text = element_text(size = 20),
        legend.position = "none")






GLM_m <- ggplot() +
  geom_spatvector(data = SEQ, fill = 'transparent', col = 'black')+
  geom_spatraster(data = glm_pred) +
  theme_cowplot(font_size = 17)+  
  scale_fill_viridis_c(na.value = 'transparent', limits = c(0.51,17), breaks = seq(1,18,1), direction = 1) +
  labs(fill = 'Fire frequency')+
  annotation_scale(location = "bl", style = 'ticks', pad_y = unit(0.5, 'cm'), pad_x = unit(15, 'cm'), text_cex = 2)+
  annotation_north_arrow(location = "bl", which_north = T, height = unit(2, "cm"), width = unit(1.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(25, 'cm'), style = north_arrow_fancy_orienteering) +
  theme(legend.key.height = unit(2.5, 'cm'),
        legend.key.width = unit(1.75, 'cm'),
        legend.title = element_text(face = 'bold', size = 25),
        legend.text = element_text(size = 20),
        legend.position = "none")





Sent <- ggplot() +
  geom_spatraster(data = Sentinel_ff) +
  theme_cowplot(font_size = 17)+  
  scale_fill_viridis_c(na.value = 'transparent', limits = c(1,17), breaks = seq(1,18,1), direction = 1) +
  geom_spatvector(data = SEQ, fill = 'transparent', col = 'black')+
  labs(fill = 'Fire frequency')+
  annotation_scale(location = "bl", style = 'ticks', pad_y = unit(0.5, 'cm'), pad_x = unit(15, 'cm'), text_cex = 2)+
  annotation_north_arrow(location = "bl", which_north = T, height = unit(2, "cm"), width = unit(1.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(25, 'cm'), style = north_arrow_fancy_orienteering) +
  theme(legend.key.height = unit(1.75, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(face = 'bold', size = 25),
        legend.text = element_text(size = 20),
        legend.position = "none")


QPWS <- ggplot()+
  geom_spatraster(data = environmental_preds$QPWS_firefreq) +
  theme_cowplot(font_size = 17)+  
  scale_fill_viridis_c(na.value = 'transparent', limits = c(1,17), breaks = seq(1,18,1), direction = 1) +
  geom_spatvector(data = SEQ, fill = 'transparent', col = 'black')+
  labs(fill = 'Fire frequency')+
  annotation_north_arrow(location = "br", which_north = T, height = unit(2, "cm"), width = unit(1.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(-0.3, 'cm'), style = north_arrow_fancy_orienteering) +
  theme(legend.key.height = unit(2.5, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(face = 'bold', size = 25),
        legend.text = element_text(size = 20))

# Produce plot with all maps

r_map_p <- plot_grid(unweighted+ theme(legend.position = "none"), downweighted + theme(legend.position = "none"), IWLR+ theme(legend.position = "none"), nrow = 1, rel_widths = c(0.2,0.2,0.2))




# 4.1 Prepare data for histogram plots ----

QPWS_pres <- extract(QPWS_ff, QPWS_rand)
QPWS_pres[is.na(QPWS_pres)] <- 0
summary(QPWS_pres)
colnames(QPWS_pres) <- c('ID', 'lyr1')
QPWS_pres[is.na(QPWS_pres)] <- 0

Sent_pres <- extract(Sentinel_ff, QPWS_rand)
colnames(Sent_pres) <- c('ID', 'lyr1')
Sent_pres$Dataset <- 'Sentinel'
head(Sent_pres)


uwt_pres <- extract(unweighted_pred, QPWS_rand)
uwt_pres <- as.data.frame(uwt_pres)
uwt_pres$ID <- 1:10000
colnames(uwt_pres) <- c('ID', 'lyr1')

dwt_pres <- extract(down_wt_pred, QPWS_rand)
dwt_pres <- as.data.frame(dwt_pres)
dwt_pres$ID <- 1:10000
colnames(dwt_pres) <- c('ID', 'lyr1')
head(dwt_pres)

IWLR_pres <- extract(IWLR_pred, QPWS_rand)
IWLR_pres <- as.data.frame(IWLR_pres)
IWLR_pres$ID <- 1:10000
colnames(IWLR_pres) <- c('ID', 'lyr1')
head(IWLR_pres)

gam_pres <- extract(gam_pred, QPWS_rand)
gam_pres <- as.data.frame(gam_pres)
gam_pres$ID <- 1:10000
colnames(gam_pres) <- c('ID', 'lyr1')
head(gam_pres)
gam_pres[is.na(gam_pres)] <- 0


glm_pres <- extract(glm_pred, QPWS_rand)
glm_pres <- as.data.frame(glm_pres)
glm_pres$ID <- 1:10000
colnames(glm_pres) <- c('ID', 'lyr1')
head(glm_pres)
glm_pres[is.na(glm_pres)] <- 0





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

dwt_c <- '#579C97'
col2rgb(dwt_c, alpha = F)

uwt_c <- '#2A6D7A'
col2rgb(uwt_c, alpha = F)

iwlr_c <- '#8FCCB4'
col2rgb(iwlr_c, alpha = F)



# 4.3.2 Produce histogram for each plot separately
# Main plots 
gb_p <- hist(round(QPWS_pres$lyr1), breaks = seq(-1,16,1))
gl_p <- hist(round(glm_pres$lyr1), breaks = seq(-1,16,1))
s_p <- hist(round(Sent_pres$lyr1), breaks = seq(-1,16,1))



# Subplots
ga_p <- hist(round(gam_pres$lyr1),  breaks = seq(-1,16,1))
dwt_p <- hist(round(dwt_pres$lyr1),  breaks = seq(-1,16,1))
uwt_p <- hist(round(uwt_pres$lyr1),  breaks = seq(-1,16,1))
iwlr_p <- hist(round(IWLR_pres$lyr1),  breaks = seq(-1,16,1))





# 4.3.3 Create the subplots with the recuded y axis
# Main plots
gb_ps <- hist(round(QPWS_pres$lyr1), ylim = c(0,100), breaks = seq(-1,16,1))
gl_ps <- hist(round(glm_pres$lyr1), ylim = c(0,100), breaks = seq(-1,16,1))
s_ps <- hist(round(Sent_pres$lyr1), ylim = c(0,100), breaks = seq(-1,16,1))

# Subplots
ga_ps <- hist(round(gam_pres$lyr1), ylim = c(0,100), breaks = seq(-1,16,1))
dwt_ps <- hist(round(dwt_pres$lyr1), ylim = c(0,100), breaks = seq(-1,16,1))
uwt_ps <- hist(round(uwt_pres$lyr1), ylim = c(0,100), breaks = seq(-1,16,1))
iwlr_ps <- hist(round(IWLR_pres$lyr1), ylim = c(0,100), breaks = seq(-1,16,1))





# 4.3.4 Combine histograms for each fire data grouping
dev.new(width = 60, height = 40, res = 300, dpi = 80, noRStudioGD = T)
par(mfrow = c(3, 2), mar = c(4,4,4,2), oma = c(1,3,1,17))
# Observed data
plot(gb_p, col = rgb(204/255, 204/255,204/255, 1), 
     ylab = "",  las = 1, ylim = c(0,0.7), freq = F, yaxt = "n",
     xlab = expression(bold('Fire frequency')), xlim = c(0,16), xaxt = "n", 
     border = 'white', main = "", cex.axis = 2, cex.lab = 2.6)
mtext(expression(bold("Density")), side = 2, cex = 2.2, line = 3.5)
plot(s_p, col = rgb(70/255, 130/255, 180/255, 0.4),  ylim = c(0,0.7), freq = F,
     ylab = expression(bold("")),  las = 1, , yaxt = "n",
     xlab = expression(bold("Fire frequency")), xlim = c(0,16), xaxt = "n",  
     border = 'white',bty= 'l', main = "", 
     add = T)
axis(side = 1, at = seq(-1,16,1), cex.axis = 1.8, line = -0.4)
axis(side = 1, at = c(10,12,14,16), cex.axis = 1.8, line = -0.4)
axis(side = 2, at = seq(0, 0.7, 0.1), cex.axis = 1.8, line = 0.3, las = 1)
mtext(expression(bold("(a) Observed")), side = 3, line = 0.8, at = 16, cex = 3)
mtext(expression(paste("Pearson's ", italic("r"), " = 0.252")), line = -2, at = 12, cex = 1.5)


plot(gb_ps, col = rgb(204/255, 204/255, 204/255, 1), 
     ylab = "", las = 1,  ylim = c(0,0.02), freq = F, yaxt = "n",
     xlab = expression(bold("Fire frequency")), xlim = c(0,16), xaxt = "n",
     border = 'white', main = "", cex.axis = 2, cex.lab = 2.6)
plot(s_ps, col = rgb(70/255, 130/255, 180/255, 0.5),
     ylab = expression(bold("")), las = 1, ylim = c(0,0.02), freq = F, yaxt = "n",
     xlab = expression(bold("Fire frequency")), xlim = c(0,16), xaxt = "n",
     border = 'white',  main = "", 
     add = T)
axis(side = 1, at = seq(-1,16,1), cex.axis = 1.8, line = -0.4)
axis(side = 1, at = c(10,12,14,16), cex.axis = 1.8, line = -0.4)
axis(side = 2, at = seq(0, 0.02, 0.01), cex.axis = 1.8, line = 0.3, las = 1)

#mtext(expression(bold("Density")), side = 2, cex = 2.2, line = 3)

par(xpd = NA)
legend('topright', inset = c(-0.52, -0.35), fill = c('gray80', 'steelblue', "#492050", '#AAA970', '#2A6D7A','#579C97','#8FCCB4'), legend = c("Public", "Satellite", "GLM", "GAM", "Unweigthed BRT", "Down-weighted BRT", "IWLR BRT"), cex = 2.5, bty = 'n', border = 'transparent')
par(xpd = F)



# GLM data
plot(gb_p, col = rgb(204/255, 204/255, 204/255, 1), 
     ylab = "", las = 1,  ylim = c(0,0.7), freq = F, yaxt = "n",
     xlab = expression(bold("Fire frequency")), xlim = c(0,16),  xaxt = "n",
     border = 'white',bty= 'l', main = "", cex.axis = 2, cex.lab = 2.6)
plot(gl_p, col = rgb(73/255, 32/255, 80/255, 0.5),
     ylab = expression(bold("Density")),  las = 1, ylim = c(0,0.7), freq = F, yaxt = "n",
     xlab = expression(bold("Fire frequency")), xlim = c(0,16),  xaxt = "n",
     border = 'white', main = "", cex.axis = 2, cex.lab = 2.6,
     add = T)
axis(side = 1, at = seq(-1,16,1), cex.axis = 1.8, line = -0.4)
axis(side = 1, at = c(10,12,14,16), cex.axis = 1.8, line = -0.4)
axis(side = 2, at = seq(0, 0.7, 0.1), cex.axis = 1.8, line = 0.3, las = 1)
mtext(expression(bold("Density")), side = 2, cex = 2.2, line = 3.5)
mtext(expression(bold("(b) GLM")), side = 3, line = 0.1, at = 16, cex = 3)
mtext(expression(paste("Pearson's ", italic("r"), " = 0.745")), line = -2, at = 12, cex = 1.5)

plot(gb_ps, col = rgb(204/255, 204/255, 204/255, 1),
     ylim = c(0,0.02), freq = F, ylab = "", las = 1, yaxt = "n",
     xlab = expression(bold("Fire frequency")), xlim = c(0,16), xaxt = "n",
     border = 'white', main = "", cex.axis = 2, cex.lab = 2.6)
plot(gl_ps, col = rgb(73/255, 32/255, 80/255, 0.5),
     ylim = c(0,0.02), freq = F, ylab = expression(bold("")), las = 1, yaxt = "n", 
     xlab = expression(bold("Fire frequency")), xlim = c(0,16), xaxt = "n",
     border = 'white', main = "", cex.lab = 2.6, 
     add = T)
axis(side = 1, at = seq(-1,16,1), cex.axis = 1.8, line = -0.4)
axis(side = 1, at = c(10,12,14,16), cex.axis = 1.8, line = -0.4)
axis(side = 2, at = seq(0, 0.02, 0.01), cex.axis = 1.8, line = 0.3, las = 1)
#mtext(expression(bold("Density")), side = 2, cex = 2.2, line = 3)


# GAM
plot(gb_p, col = rgb(204/255, 204/255, 204/255, 1), 
     ylab = "", las = 1,  ylim = c(0,0.7), freq = F,yaxt = "n",
     xlab = expression(bold("Fire frequency")), xlim = c(0,16),  xaxt = "n",
     border = 'white',bty= 'l', main = "", cex.axis = 2, cex.lab = 2.6)
plot(ga_p, col = rgb(170/255, 169/255, 112/255, 0.5),
     ylab = expression(bold("Density")),  las = 1, ylim = c(0,0.7), freq = F,
     xlab = expression(bold("Fire frequency")), xlim = c(0,16),  xaxt = "n",
     border = 'white', main = "", cex.axis = 2, cex.lab = 2.6,
     add = T)
axis(side = 1, at = seq(-1,16,1), cex.axis = 1.8, line = -0.4)
axis(side = 1, at = c(10,12,14,16), cex.axis = 1.8, line = -0.4)
axis(side = 2, at = seq(0, 0.7, 0.1), cex.axis = 1.8, line = 0.3, las = 1)
mtext(expression(bold("Density")), side = 2, cex = 2.2, line = 3.5)
mtext(expression(bold("(c) GAM")), side = 3, line = 0.1, at = 16, cex = 3)
mtext(expression(paste("Pearson's ", italic("r"), " = 0.473")), line = -2, at = 12, cex = 1.5)

plot(gb_ps, col = rgb(204/255, 204/255, 204/255, 1),
     ylim = c(0,0.02), freq = F, ylab = "", las = 1,yaxt = "n",
     xlab = expression(bold("Fire frequency")), xlim = c(0,16), xaxt = "n",
     border = 'white', main = "", cex.axis = 2, cex.lab = 2.6)
plot(ga_ps, col = rgb(170/255, 169/255, 112/255, 0.5),
     ylim = c(0,0.02), freq = F, ylab = expression(bold("")), las = 1, yaxt = "n",
     xlab = expression(bold("Fire frequency")), xlim = c(0,16), xaxt = "n",
     border = 'white', main = "", cex.lab = 2.6, 
     add = T)
axis(side = 1, at = seq(-1,16,1), cex.axis = 1.8, line = -0.4)
axis(side = 1, at = c(10,12,14,16), cex.axis = 1.8, line = -0.4)
axis(side = 2, at = seq(0, 0.02, 0.01), cex.axis = 1.8, line = 0.3, las = 1)
#mtext(expression(bold("Density")), side = 2, cex = 2.2, line = 3)






### Subplots
dev.new(width = 20, height = 8, res = 300, dpi = 80, noRStudioGD = T)
par(mfrow = c(2, 4), mar = c(6,3,4,2), oma = c(1,3,0,0))



# UWT
plot(gb_p, col = rgb(204/255, 204/255, 204/255, 1), 
     ylab = "", las = 1,  ylim = c(0,0.7), freq = F,yaxt = "n",
     xlab = expression(bold("Fire frequency")), xlim = c(-1,16),  xaxt = "n",
     border = 'white',bty= 'l', main = "", cex.axis = 2, cex.lab = 2.6)
plot(uwt_p, col = rgb(42/255, 109/255, 122/255, 0.5),
     ylab = expression(bold("Density")),  las = 1, ylim = c(0,0.7), freq = F,yaxt = "n",
     xlab = expression(bold("Fire frequency")), xlim = c(-1,16),  xaxt = "n",
     border = 'white', main = "", cex.axis = 2, cex.lab = 2.6,
     add = T)
axis(side = 1, at = c(0,5,10,16), cex.axis = 1.8, line = -0.4)
axis(side = 1, at = seq(0,16, 1), labels = F, line = -0.4)
axis(side = 2, at = seq(0, 0.7, 0.1), cex.axis = 1.8, line = -0.8, las = 1)
axis(side = 1, at = c(-1, 16), labels = F, line = -0.4)
mtext(expression(bold("Density")), side = 2, cex = 2.2, line = 2.5)
mtext(expression(bold("(d) Unweighted BRT")), side = 3, line = 1, at = 16, cex = 2)
mtext(expression(paste("Pearson's ", italic("r"), " = 0.338")), line = -2, at = 9, cex = 1.75)

plot(gb_ps, col = rgb(204/255, 204/255, 204/255, 1),
     ylim = c(0,0.02), freq = F, ylab = "", las = 1, yaxt = "n",
     xlab = expression(bold("Fire frequency")), xlim = c(-1,16), xaxt = "n",
     border = 'white', main = "", cex.axis = 2, cex.lab = 2.6)
plot(uwt_ps, col = rgb(42/255, 109/255, 122/255, 0.5),
     ylim = c(0,0.02), freq = F, ylab = expression(bold("")), las = 1, yaxt = "n",
     xlab = expression(bold("Fire frequency")), xlim = c(-1,16), xaxt = "n",
     border = 'white', main = "", cex.lab = 2.6, 
     add = T)
axis(side = 1, at = c(0,5,10,16), cex.axis = 1.8, line = -0.4)
axis(side = 1, at = seq(0,16, 1), labels = F, line = -0.4)
axis(side = 1, at = c(-1, 16), labels = F, line = -0.4)
axis(side = 2, at = seq(0, 0.02, 0.01), cex.axis = 1.8, line = -0.8, las = 1)
#mtext(expression(bold("Density")), side = 2, cex = 2.2, line = 3)



# DWT
plot(gb_p, col = rgb(204/255, 204/255, 204/255, 1), 
     ylab = "", las = 1,  ylim = c(0,0.7), freq = F,yaxt = "n",
     xlab = expression(bold("Fire frequency")), xlim = c(-1,16),  xaxt = "n",
     border = 'white',bty= 'l', main = "", cex.axis = 2, cex.lab = 2.6)
plot(dwt_p, col = rgb(87/255, 156/255, 151/255, 0.5), 
     ylab = expression(bold("Density")),  las = 1, ylim = c(0,0.7), freq = F,yaxt = "n",
     xlab = expression(bold("Fire frequency")), xlim = c(-1,16),  xaxt = "n",
     border = 'white', main = "", cex.axis = 2, cex.lab = 2.6,
     add = T)
axis(side = 1, at = c(0,5,10,16), cex.axis = 1.8, line = -0.4)
axis(side = 1, at = seq(0,16, 1), labels = F, line = -0.4)
axis(side = 2, at = seq(0, 0.7, 0.1), cex.axis = 1.8, line = -0.8, las = 1)
axis(side = 1, at = c(-1, 16), labels = F, line = -0.4)
mtext(expression(bold("Density")), side = 2, cex = 2.2, line = 2.5)
mtext(expression(bold("(e) Down-weighted BRT")), side = 3, line = 1, at = 16, cex = 2)
mtext(expression(paste("Pearson's ", italic("r"), " = 0.337")), line = -2, at = 9, cex = 1.75)

plot(gb_ps, col = rgb(204/255, 204/255, 204/255, 1),
     ylim = c(0,0.02), freq = F, ylab = "", las = 1, yaxt = "n",
     xlab = expression(bold("Fire frequency")), xlim = c(-1,16), xaxt = "n",
     border = 'white', main = "", cex.axis = 2, cex.lab = 2.6)
plot(dwt_ps, col = rgb(87/255, 156/255, 151/255, 0.5),
     ylim = c(0,0.02), freq = F, ylab = expression(bold("")), las = 1, yaxt = "n",
     xlab = expression(bold("Fire frequency")), xlim = c(-1,16), xaxt = "n",
     border = 'white', main = "", cex.lab = 2.6, 
     add = T)
axis(side = 1, at = c(0,5,10,16), cex.axis = 1.8, line = -0.4)
axis(side = 1, at = seq(0,16, 1), labels = F, line = -0.4)
axis(side = 1, at = c(-1, 16), labels = F, line = -0.4)
axis(side = 2, at = seq(0, 0.02, 0.01), cex.axis = 1.8, line = -0.8, las = 1)
#mtext(expression(bold("Density")), side = 2, cex = 2.2, line = 3)



# IWLR
plot(gb_p, col = rgb(204/255, 204/255, 204/255, 1), 
     ylab = "", las = 1,  ylim = c(0,0.7), freq = F,yaxt = "n",
     xlab = expression(bold("Fire frequency")), xlim = c(-1,16),  xaxt = "n",
     border = 'white',bty= 'l', main = "", cex.axis = 2, cex.lab = 2.6)
plot(iwlr_p, col = rgb(143/255, 204/255, 180/255, 0.5),
     ylab = expression(bold("Density")),  las = 1, ylim = c(0,0.7), freq = F,yaxt = "n",
     xlab = expression(bold("Fire frequency")), xlim = c(-1,16),  xaxt = "n",
     border = 'white', main = "", cex.axis = 2, cex.lab = 2.6,
     add = T)
axis(side = 1, at = c(0,5,10,16), cex.axis = 1.8, line = -0.4)
axis(side = 1, at = seq(0,16, 1), labels = F, line = -0.4)
axis(side = 2, at = seq(0, 0.7, 0.1), cex.axis = 1.8, line = -0.8, las = 1)
axis(side = 1, at = c(-1, 16), labels = F, line = -0.4)
mtext(expression(bold("Density")), side = 2, cex = 2.2, line = 2.5)
mtext(expression(bold("(f) IWLR BRT")), side = 3, line = 1, at = 16, cex = 2)
mtext(expression(paste("Pearson's ", italic("r"), " = 0.038")), line = -2, at = 9, cex = 1.75)

plot(gb_ps, col = rgb(204/255, 204/255, 204/255, 1),
     ylim = c(0,0.02), freq = F, ylab = "", las = 1, yaxt = "n",
     xlab = expression(bold("Fire frequency")), xlim = c(-1,16), xaxt = "n",
     border = 'white', main = "", cex.axis = 2, cex.lab = 2.6)
plot(iwlr_ps, col = rgb(143/255, 204/255, 180/255, 0.5),
     ylim = c(0,0.02), freq = F, ylab = expression(bold("")), las = 1, yaxt = "n",
     xlab = expression(bold("Fire frequency")), xlim = c(-1,16), xaxt = "n",
     border = 'white', main = "", cex.lab = 2.6, 
     add = T)
axis(side = 1, at = c(0,5,10,16), cex.axis = 1.8, line = -0.4)
axis(side = 1, at = seq(0,16, 1), labels = F, line = -0.4)
axis(side = 1, at = c(-1, 16), labels = F, line = -0.4)
axis(side = 2, at = seq(0, 0.02, 0.01), cex.axis = 1.8, line = -0.8, las = 1)
#mtext(expression(bold("Density")), side = 2, cex = 2.2, line = 3)

save.image('./02_Workspaces/005_predictive_model_validations.RData')
#load('./02_Workspaces/005_predictive_model_validations.RData')




# We also want to determine what % of cells we have taken out of maps for Satllite, Unweighted BRT, and Downweighted BRT that were above 16 fires.
plot(Sentinel_ff$Fire_freq)
unique(terra::cells(Sentinel_ff, c(17:26))) # Number of cells with ff above 16 = 32945
cells(Sentinel_ff) # Total number of cells = 69976534
       
# Number of cells excluded 
(32945/69976534)*100 # 0.04708007% of cells with Fire freq >16

cells(unweighted_pred) # Total cells =  70535296. Number unchanged by adding round
cells(round(unweighted_pred), c(17:64)) # 2012 cells. Needed to add round due to the decimals in predictions
(2021/70535296)*100 # 0.002865232% of cells with fire freq >16

cells(down_wt_pred) # Total cells = 70535296
cells(round(down_wt_pred), c(17:72)) # 1413 cells
(1413/70535296)*100 # 0.002003252% of cells with fire freq >16
