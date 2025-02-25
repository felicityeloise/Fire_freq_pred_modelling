# Written by Felicity Charles
# Date: 06/10/2024

##### Fire frequency analysis ----
# This script validates model predictions, ensuring predictions are reasonable given the vegetation type and recommended fire regime. 

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



# Further ground truthing/model validation can be performed to ensure the predictions seem realistic considering the broad vegetation group, regional ecosystem fire regime management suggestions, and land use.
# 6. Examine model predictions based on vegetation type -----
# Group BVGs 2-14, keeping 1. Rainforest and scrubs, 15. Wetlands, and 16. Mangroves and saltmarshes separate as these areas are less likely to burn than forests or woodlands
# Want to take DBVG5M, the column of interest from BVG 


# 6.1 Extract from BVGs any areas that correspond to DBVG5M - 1, 15, or 16 ----
lowburn_BVGs <- BVG %>% 
  filter(BVG$dbvg5m == "1" | BVG$dbvg5m == "15" | BVG$dbvg5m == "16")


# 6.2 Read in presence background dataframes and point locations ----
Rand_fire <- read.csv('./00_Data/Fire_data/Outputs/Random_points_data/Fire_frequency_random_environmental_pres_resampled.csv', header = T)
head(Rand_fire); dim(Rand_fire)
Rand_fire <- Rand_fire[, c(4, 5:6)]
pres_pts <-  vect('./00_Data/Fire_data/Outputs/QPWS_random.gpkg')


Background_data <- read.csv('./00_Data/Fire_data/Outputs/Background_points_data/Fire_frequency_background_environmental_data_resampled.csv', header = T)
head(Background_data); dim(Background_data)
Background_data <- Background_data[, c(3, 5:6)]
bg_pts <- vect('./00_Data/Fire_data/Outputs/Background_points_data/bg_rand.gpkg')


# 6.3 Extract information for presence and background points ---- 
lowb_BVG <- extract(lowburn_BVGs, pres_pts)

lowb_BVG_bg <- extract(lowburn_BVGs, bg_pts)

fire_pred <- extract(gam_pred, pres_pts)

fire_predbg <- extract(gam_pred, bg_pts)

bg_QPWS <- extract(QPWS_ff, bg_pts)

bg_Sentinel <- extract(Sentinel_ff, bg_pts)
bg_Sentinel[is.na(bg_Sentinel$focal_mean),] <- 0

# 6.3.1 Add extracted information for presence and background points to the dataframes ----
# Want to code if it is a high burn BVG as 0s first, all else will be NA. Then all others we want to take the BVG code to add to this column. Similar to what we did for nearest land. 
Rand_fire$BVG_firecat <- lowb_BVG$dbvg5m
Rand_fire$fire_pred <- round(fire_pred$lyr1)
head(Rand_fire)
tail(Rand_fire)
unique(Rand_fire$BVG_firecat)



Background_data$BVG_firecat <- lowb_BVG_bg$dbvg5m
Background_data$fire_pred <- round(fire_predbg$lyr1)
Background_data$QPWS_NAtoSent <- bg_QPWS$QPWS # Add new column for the QPWS data that still contains NA values so that we can replace any rows with NA values with the Sentinel value.
Background_data$QPWS_NAtoSent <- ifelse(is.na(Background_data$QPWS_NAtoSent), round(bg_Sentinel$focal_mean), Background_data$QPWS_NAtoSent)
head(Background_data)
tail(Background_data)
unique(Background_data$BVG_firecat)
unique(Background_data$QPWS_NAtoSent)

# 6.3.2 Remove rows with NA values for BVG fire category ----
# As we are only interested in knowing how the predictions are distributed for the low burn BVGs, we will remove any rows of data that does not correspond to these BVGs.
Rand_fire <- Rand_fire[!is.na(Rand_fire$BVG_firecat),]
head(Rand_fire)

Background_data <- Background_data[!is.na(Background_data$BVG_firecat),]
head(Rand_fire)


# 6.4 What are the fire frequencies that we see for BVGs that we do not expect to burn? ----
# Format the data in such a way that we can produce histograms showing the distribution of fire frequencies for presence points and background points. We can see whether predictions match expectations for QPWS land, which we know should be accurate. We can also determine if our expectations are met outside of QPWS land, with comparisons made to the Sentinel fire frequency data where QPWS data was not avaialble for the background points to confirm that prediction accuracy is not degraded outside of QPWS estates for these vegetation groups. 

# 6.4.1 Reformat the data ----
lowburn_BVG_pres_gam <- Rand_fire[, c(5,4)]
colnames(lowburn_BVG_pres_gam) <- c("fire_freq", "BVG")
lowburn_BVG_pres_gam$Dataset <- "GAM"

lowburn_BVG_pres_QPWS <- Rand_fire[, c(1,4)]
colnames(lowburn_BVG_pres_QPWS) <- c("fire_freq", "BVG")
lowburn_BVG_pres_QPWS$Dataset <- "QPWS"

lowburn_BVG_pres <- rbind(lowburn_BVG_pres_QPWS, lowburn_BVG_pres_gam)
head(lowburn_BVG_pres)
str(lowburn_BVG_pres)
lowburn_BVG_pres$Dataset <- as.factor(lowburn_BVG_pres$Dataset)
levels(lowburn_BVG_pres$Dataset) <- c("QPWS-Sentinel", "GAM")
lowburn_BVG_pres$Dataset <- relevel(lowburn_BVG_pres$Dataset, "GAM")
head(lowburn_BVG_pres)
tail(lowburn_BVG_pres)
unique(lowburn_BVG_pres)
lowburn_BVG_pres[is.na(lowburn_BVG_pres$fire_freq), 1] <- 0
unique(lowburn_BVG_pres)

# Change BVG level names
str(lowburn_BVG_pres)
lowburn_BVG_pres$BVG <- factor(lowburn_BVG_pres$BVG)
levels(lowburn_BVG_pres$BVG) <- c("1. Rainforest and scrubs", "15. Wetlands", "16. Mangroves and saltmarshes")



head(Background_data)
lowburn_BVG_bg_gam <- Background_data[, c(5,4)]
colnames(lowburn_BVG_bg_gam) <- c("fire_freq", "BVG")
lowburn_BVG_bg_gam$Dataset <- "GAM"

lowburn_BVG_bg_QPWS <- Background_data[, c(6,4)]
colnames(lowburn_BVG_bg_QPWS) <- c("fire_freq", "BVG")
lowburn_BVG_bg_QPWS$Dataset <- "QPWS"

lowburn_BVG_bg <- rbind(lowburn_BVG_bg_QPWS, lowburn_BVG_bg_gam)
colnames(lowburn_BVG_bg) <- c('fire_freq', 'BVG', 'Dataset')
str(lowburn_BVG_bg)
lowburn_BVG_bg$Dataset <- as.factor(lowburn_BVG_bg$Dataset)
levels(lowburn_BVG_bg$Dataset) <- c('QPWS-Sentinel', 'GAM')
lowburn_BVG_bg$Dataset <- relevel(lowburn_BVG_bg$Dataset, "GAM")
head(lowburn_BVG_bg)
tail(lowburn_BVG_bg)
unique(lowburn_BVG_bg) 
lowburn_BVG_bg[is.na(lowburn_BVG_bg$fire_freq), 1] <- 0
unique(lowburn_BVG_bg) 

# Change BVG level names
str(lowburn_BVG_bg)
lowburn_BVG_bg$BVG <- factor(lowburn_BVG_bg$BVG)
levels(lowburn_BVG_bg$BVG) <- c("1. Rainforest and scrubs", "15. Wetlands", "16. Mangroves and saltmarshes")


# 6.4.2 Produce histogram plots -----
rainforest_pres <- ggplot(lowburn_BVG_pres[lowburn_BVG_pres$BVG == '1. Rainforest and scrubs', ], aes(x = fire_freq, fill = Dataset))+
  geom_histogram(colour ='#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('#AAA970', 'slategray3')) +
  scale_x_continuous('', breaks = seq(0,6,1)) +
  scale_y_continuous('Count', limits = c(0,700)) +
  theme(legend.title = element_text(size = 15, face = 'bold'),
        axis.title = element_text(size = 15, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(size = 15))+
  labs(title = expression(bold("(a)")))+
  facet_zoom(ylim = c(0,50), zoom.size = 1)



rainforest_bg <- ggplot(lowburn_BVG_bg[lowburn_BVG_bg$BVG == '1. Rainforest and scrubs',  ], aes(x = fire_freq, fill = Dataset))+
  geom_histogram(colour ='#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('#AAA970', 'slategray3')) +
  scale_x_continuous('Fire frequency', breaks = seq(0,13,1)) +
  scale_y_continuous('', limits = c(0,8000)) +
  theme(legend.title = element_text(size = 15, face = 'bold'),
        axis.title = element_text(size = 15, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(size = 15))+
  labs(title = expression(bold("(b)")))+
  facet_zoom(ylim = c(0,200), zoom.size = 1)



non_rainforest_pres <-  ggplot(lowburn_BVG_pres[lowburn_BVG_pres$BVG == '15. Wetlands' | lowburn_BVG_pres$BVG ==  "16. Mangroves and saltmarshes",], aes(x = fire_freq, fill = Dataset))+
  geom_histogram(colour ='#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('#AAA970', 'slategray3')) +
  scale_x_continuous('Fire frequency', breaks = seq(0,6,1)) +
  scale_y_continuous('Count') +
  theme(legend.title = element_text(size = 15, face = 'bold'),
        axis.title = element_text(size = 15, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(size = 15))+
  labs(title = expression(bold("(a)")))


nonrainforest_bg <- ggplot(lowburn_BVG_bg[lowburn_BVG_bg$BVG != '1. Rainforest and scrubs',  ], aes(x = fire_freq, fill = Dataset))+
  geom_histogram(colour ='#e9ecef', alpha = 0.6, position = 'identity', binwidth = 1)+
  scale_fill_manual(values = c('#AAA970', 'slategray3')) +
  scale_x_continuous('Fire frequency', breaks = seq(0,13,1)) +
  scale_y_continuous('', limits = c(0,450)) +
  theme(legend.title = element_text(size = 15, face = 'bold'),
        axis.title = element_text(size = 15, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.key.size = unit(1, 'cm'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        plot.title = element_text(size = 15))+
  labs(title = expression(bold("(b)")))+
  facet_zoom(ylim = c(0,70), zoom.size = 1)




legend <- get_legend(non_rainforest_pres + theme(legend.box.margin = margin(0,0,0,15)))


rainforest_plot <- plot_grid(rainforest_pres+ theme(legend.position = "none"), rainforest_bg+ theme(legend.position = "none"), nrow = 2)
rainforest_plot <- plot_grid(rainforest_plot, legend, rel_heights = c(4,0.6), ncol = 1)


non_rainforest_plot <- plot_grid(non_rainforest_pres + theme(legend.position = "none"), nonrainforest_bg + theme(legend.position = "none"), nrow = 2)
non_rainforest_plot <- plot_grid(non_rainforest_plot, legend, rel_heights = c(5, 0.6), ncol = 1)

#save.image('./02_Workspaces/005_predictive_model_validations.RData')

# 7. RE. Compare estimates of fire frequency based on fire management guidelines to the predictions ----
# Here we will compare fire frequency estimates to the actual values recorded for both QPWS and Sentinel, and also determine how our predicted values compare to fire management guidelines. 

# Maximum fire frequency estimate column based on the interval_min, minimum fire frequency estimate based on the interval_max

# 7.1 Read in the data and remove unnecessary information ----
RE_fire <- read.csv('./00_Data/Environmental_data/Regional_ecosystem_fire_management_guidelines/fire-management-guidelines-v13.1.csv', header = T, stringsAsFactors = T)
RE_fire <- RE_fire[, c(2, 9, 13, 14)]
head(RE_fire)
str(RE_fire)

# 7.1.1 Remove any rows with NA values for fire return interval information
RE_fire <- RE_fire[!is.na(RE_fire$INTERVAL_MIN),]
str(RE_fire)
unique(RE_fire$INTERVAL_MIN)
unique(RE_fire$INTERVAL_MAX)  



# 7.2 Calculate fire frequency estimates for maximum and minimum fire frequency ---
RE_fire$max_firefreq <- round(35/RE_fire$INTERVAL_MIN)
RE_fire$min_firefreq <- round(35/RE_fire$INTERVAL_MAX)

# 7.3 Extract regional ecosystem for a subset of presence background points, then we can use this to pattern match to this dataset for extraction of minimum and maximum estimated fire frequency. 
# Look at the data
RE # We are interested in the column re1
unique(RE$re1) # We need to remove any letters and anything appearing after this from re1
length(unique(RE$re1)) #303 REs


# 7.3.1 Perform pattern replacement ----
# Some REs contain an x followed by a number meaning we cannot just remove a-z from the character string as the RE may then be misidentified due to the following number. Firstly, identify which REs have an x for pattern replacement of these REs. Then use pattern replacement to remove any alphabetical characters from the RE value.
RE2 <- as.data.frame(RE)
RE2$x <- grepl('x', RE2$re1, fixed = T) # Find which rows have an x in the character string for re1
RE2 <- RE2[RE2$x == T, 2]
unique(RE2)

RE$re1 <- gsub('12.11.9x1', '12.11.9', RE$re1)
RE$re1 <- gsub('12.12.19x5', '12.12.19', RE$re1)
RE$re1 <- gsub('12.12.19x2', '12.12.19', RE$re1)
RE$re1 <- gsub('12.12.19x3', '12.12.19', RE$re1)
RE$re1 <- gsub('12.9-10.1x1', '12.9-10', RE$re1)
unique(RE$re1)
length(unique(RE$re1)) # 299 REs


RE$re1 <- gsub("[a-z]", "", RE$re1)
unique(RE$re1)


# 7.3.3 Add columns to RE spatial data frame for minimum estimated fire frequency and maximum estimated fire frequency ----
# Check the dataframes
head(RE_fire);head(RE)
# Change the column names of RE_fire to match those in RE as well
colnames(RE_fire) <- c("re1", "dbvg5m", "interval_min", "interval_max", "max_ff_est", "min_ff_est")


# 7.3.3.1 Add fire interval data to the spatial regional ecosystem data ----
# Create a new dataset that contains the rows from RE
RE_df <- left_join(as.data.frame(RE),
                   RE_fire,
                   by = 're1')
head(RE_df) # Check how this looks

# Take the columns we are interested in from this dataframe and add this to the spatial RE dataset
RE
RE$interval_min <- RE_df$interval_min
head(RE)
RE$interval_max <- RE_df$interval_max
RE$max_ff_est <- RE_df$max_ff_est
RE$min_ff_est <- RE_df$min_ff_est
head(RE)

# 7.3.4 Create a random sample of 150 points across SEQ for validation purposes ----
set.seed(436)
rand_pts <- spatSample(gam_pred, 5000, "random", as.points = T)
rand_pts
plet(rand_pts)

set.seed(436)
sm_rand_pts <- spatSample(gam_pred, 150, "random", as.points = T)

# 7.3.5 Extract fire frequency and regional ecosystem information for points ----
RE_pts <- extract(RE, rand_pts) # This contains the regional ecosystem, BVG, and estimated maximum and minimum fire frequencies
QPWS_pts <- extract(QPWS_ff, rand_pts)
pred_pts <- extract(gam_pred, rand_pts)
unique(round(pred_pts$lyr1))
Sent_pts <- extract(Sentinel_ff, rand_pts)


RE_pts <- extract(RE, sm_rand_pts) # This contains the regional ecosystem, BVG, and estimated maximum and minimum fire frequencies
QPWS_pts <- extract(QPWS_ff, sm_rand_pts)
pred_pts <- extract(gam_pred, sm_rand_pts)
unique(round(pred_pts$lyr1))
Sent_pts <- extract(Sentinel_ff, sm_rand_pts)
# 7.3.6 Investigate fire frequency observed and predictions against minimum and maximum fire frequency estimates from fire management guidelines ----

RE_randpt_fire <- QPWS_pts

# Replace any NA values for QPWS with Sentinel data
RE_randpt_fire$QPWS_NAtoSent <- ifelse(is.na(RE_randpt_fire$QPWS_SEQ_freq_raster), Sent_pts$focal_mean, RE_randpt_fire$QPWS_SEQ_freq_raster)

RE_randpt_fire$GAM_ff <- pred_pts$lyr1
RE_randpt_fire$RE <- RE_pts$re1
RE_randpt_fire$BVG <- RE_pts$dbvg5m.x
RE_randpt_fire$max_ff_est <- RE_pts$max_ff_est
RE_randpt_fire$min_ff_est <- RE_pts$min_ff_est

head(RE_randpt_fire) # Check how this looks
RE_randpt_fire <- RE_randpt_fire[!is.na(RE_randpt_fire$max_ff_est),]
RE_randpt_fire # Look at the data



sm_RE_randpt_fire <- QPWS_pts

# Replace any NA values for QPWS with Sentinel data
sm_RE_randpt_fire$QPWS_NAtoSent <- ifelse(is.na(sm_RE_randpt_fire$QPWS_SEQ_freq_raster), Sent_pts$focal_mean, sm_RE_randpt_fire$QPWS_SEQ_freq_raster)

sm_RE_randpt_fire$GAM_ff <- pred_pts$lyr1
sm_RE_randpt_fire$RE <- RE_pts$re1
sm_RE_randpt_fire$BVG <- RE_pts$dbvg5m.x
sm_RE_randpt_fire$max_ff_est <- RE_pts$max_ff_est
sm_RE_randpt_fire$min_ff_est <- RE_pts$min_ff_est

head(sm_RE_randpt_fire) # Check how this looks
sm_RE_randpt_fire <- sm_RE_randpt_fire[!is.na(sm_RE_randpt_fire$max_ff_est),]
sm_RE_randpt_fire # Look at the data


# Produce histograms of the predictions against the maximum and minimum intervals
plot(round(RE_randpt_fire$GAM_ff), RE_randpt_fire$max_ff_est)

unique(RE_randpt_fire$RE)

# Not sure if this will be something to plot or not because each RE would have different requirements.