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
library(ggforce) # ggforce_0.4.2 
library(geomtextpath) # geomtextpath_0.1.5

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


# Further ground truthing/model validation can be performed to ensure the predictions seem realistic considering the broad vegetation group, regional ecosystem fire regime management suggestions, and land use.

# 3. Compare estimates of fire frequency based on fire management guidelines to the predictions ----
# Here we will compare fire frequency estimates to the actual values recorded for both QPWS and Sentinel, and also determine how our predicted values compare to fire management guidelines. 

# Maximum fire frequency estimate column based on the interval_min, minimum fire frequency estimate based on the interval_max

# 3.1 Read regional ecosystem fire management data and remove unnecessary information ----
RE_fire <- read.csv('./00_Data/Environmental_data/Regional_ecosystem_fire_management_guidelines/fire-management-guidelines-v13.1.csv', header = T, stringsAsFactors = T)
RE_fire <- RE_fire[, c(2, 9, 13, 14)]
head(RE_fire)

RE_fire <- RE_fire[!is.na(RE_fire$INTERVAL_MIN),]
str(RE_fire)
unique(RE_fire$INTERVAL_MIN)
unique(RE_fire$INTERVAL_MAX)  



# 3.2 Calculate fire frequency estimates for maximum and minimum fire frequency ---
RE_fire$max_firefreq <- round(35/RE_fire$INTERVAL_MIN)
RE_fire$min_firefreq <- round(35/RE_fire$INTERVAL_MAX)

# 3.3 Clean up regional ecosystems dataset
# Look at the data
RE # We are interested in the column re1
unique(RE$re1) # We need to remove any letters and anything appearing after this from re1
length(unique(RE$re1)) #303 REs


# 3.3.1 Perform pattern replacement ----
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


# Remove RE categories for:  plantation, water, canal, estuary, ocean, sand, and non-remnant
RE <- RE[RE$re1 != "plantation", ]
RE <- RE[RE$re1 != "water",]
RE <- RE[RE$re1 != "non-remnant",]
RE <- RE[RE$re1 != "canal",]
RE <- RE[RE$re1 != "estuary",]
RE <- RE[RE$re1 != "ocean",]
RE <- RE[RE$re1 != "sand",]

unique(RE$re1)
head(RE)

RE$re1 <- gsub("[a-z]", "", RE$re1) # Now remove the a-z letters that appear at the end of some numbered REs
unique(RE$re1) 

# 3.3.4 Add fire interval information to RE dataset  ----
# Check the dataframes
head(RE_fire);head(RE)
# Change the column names of RE_fire to match those in RE
colnames(RE_fire) <- c("re1", "dbvg5m", "interval_min", "interval_max", "max_ff_est", "min_ff_est")

RE_df <- left_join(as.data.frame(RE),
                   RE_fire,
                   by = 're1')
head(RE_df) # Check how this looks

RE
RE$interval_min <- RE_df$interval_min
head(RE)
RE$interval_max <- RE_df$interval_max
RE$max_ff_est <- RE_df$max_ff_est
RE$min_ff_est <- RE_df$min_ff_est
head(RE)

str(RE)

# 4. Extract fire frequency and regional ecosystem information for points 1000 random points for each BVG aggregation ----
Rainforest <- RE[RE$dbvg5m == 1, ]
Rainforest$Agg <- "Rainforest"

Sclerophyll <- RE[RE$dbvg5m == 2 | RE$dbvg5m == 3 | RE$dbvg5m == 4 | RE$dbvg5m == 5 | RE$dbvg5m == 6 | RE$dbvg5m == 7 | RE$dbvg5m == 8 | RE$dbvg5m == 9 | RE$dbvg5m == 10 | RE$dbvg5m == 11, ]
Sclerophyll$Agg <- "Sclerophyll"
unique(Sclerophyll$dbvg5m)

Shrubland <- RE[RE$dbvg5m == 12 | RE$dbvg5m == 13 | RE$dbvg5m == 14, ]
Shrubland$Agg <- "Grassland & Shrubland"

Other <- RE[RE$dbvg5m == 15 | RE$dbvg5m == 16, ]
Other$Agg <- "Other"


set.seed(480)
Rain_pts <- spatSample(Rainforest, 1000, "random") 
plet(Rain_pts)
dim(Rain_pts)

set.seed(480)
Scle_pts <- spatSample(Sclerophyll, 1000, "random")
plet(Scle_pts)
dim(Scle_pts)

set.seed(480)
Shrub_pts <- spatSample(Shrubland, 1000, "random")
plet(Shrub_pts)
dim(Shrub_pts)

set.seed(480)
Other_pts <- spatSample(Other, 1663, "random") # Need to increse the number of points for this group to produce 1000 points
plet(Other_pts)
dim(Other_pts)


# Combine these points into one
RE_pts <- rbind(Rain_pts, Scle_pts, Shrub_pts, Other_pts)
head(RE_pts); tail(RE_pts); dim(RE_pts) # Make sure this looks right
unique(RE_pts$Agg)


# Create a raster stack to extract information
names(QPWS_ff) <- "Ground based"
names(Sentinel_ff) <- "Satellite"
names(glm_pred) <- "GLM"

glm_pred <- round(glm_pred)
glm_pred[is.na(glm_pred)] <- 0 

enviro_info <- c(QPWS_ff, Sentinel_ff, glm_pred)

# Extract information
enviro_info_pts <- extract(enviro_info, RE_pts)
head(enviro_info_pts)
unique(is.na(enviro_info_pts))
dim(enviro_info_pts)



RE_ext <- extract(RE, RE_pts) # This contains the regional ecosystem, BVG, and estimated maximum and minimum fire frequencies
head(RE_ext)
RE_ext_subset <- RE_ext[, c(1,3,29,32:ncol(RE_ext))]
head(RE_ext_subset)
unique(RE_ext_subset$re1)

# 5. Create dataframe from extracted information ----
colnames(enviro_info_pts) <- c("id.y", "Ground_based", "Satellite", "GLM")

RE_randpt_fire <- RE_ext_subset
RE_randpt_fire <- left_join(RE_ext_subset, enviro_info_pts, by = 'id.y')
head(RE_randpt_fire)
unique(RE_randpt_fire$re1)

# Add BVG aggregation information
RE_pts_df <- as.data.frame(RE_pts[, 35])
RE_pts_df$id.y <- 1:4000
head(RE_pts_df)
RE_randpt_fire <- left_join(RE_randpt_fire, RE_pts_df, by = 'id.y')
head(RE_randpt_fire)
unique(RE_randpt_fire$Agg)

# Overwrite NA values

RE_randpt_fire$Satellite[is.na(RE_randpt_fire$Satellite)] <- 0
unique(is.na(RE_randpt_fire$Satellite))

RE_randpt_fire$Ground_based[is.na(RE_randpt_fire$Ground_based)] <- 0
unique(is.na(RE_randpt_fire))
unique(RE_randpt_fire$Agg)
head(RE_randpt_fire)

# Create new dataset with average and standard deviation information for each dataset
Rainforest_fire <- RE_randpt_fire[RE_randpt_fire$Agg == "Rainforest", ]
Sclero_fire <- RE_randpt_fire[RE_randpt_fire$Agg == "Sclerophyll", ]
Shrub_fire <- RE_randpt_fire[RE_randpt_fire$Agg == "Grassland & Shrubland", ]
Other_fire <- RE_randpt_fire[RE_randpt_fire$Agg == "Other", ]


Rainforest_fire_df <- as.data.frame(rbind(mean(Rainforest_fire$Ground_based), mean(Rainforest_fire$Satellite), mean(Rainforest_fire$GLM)))
colnames(Rainforest_fire_df) <- "fire_freq_mean"
Rainforest_fire_df$Dataset <- as.factor(c("Ground_based", "Satellite", "GLM"))
Rainforest_fire_df$sdv <- rbind(sd(Rainforest_fire$Ground_based), sd(Rainforest_fire$Satellite), sd(Rainforest_fire$GLM))
Rainforest_fire_df$Min_ff <- round(mean(Rainforest_fire$min_ff_est))
Rainforest_fire_df$Max_ff <- round(mean(Rainforest_fire$max_ff_est))
Rainforest_fire_df$Agg <- as.factor("Rainforest")
head(Rainforest_fire_df)




Sclero_fire_df <- as.data.frame(rbind(mean(Sclero_fire$Ground_based), mean(Sclero_fire$Satellite), mean(Sclero_fire$GLM)))
colnames(Sclero_fire_df) <- "fire_freq_mean"
Sclero_fire_df$Dataset <- as.factor(c("Ground_based", "Satellite", "GLM"))
Sclero_fire_df$sdv <- rbind(sd(Sclero_fire$Ground_based), sd(Sclero_fire$Satellite), sd(Sclero_fire$GLM))
Sclero_fire_df$Min_ff <- round(mean(Sclero_fire$min_ff_est))
Sclero_fire_df$Max_ff <- round(mean(Sclero_fire$max_ff_est))
Sclero_fire_df$Agg <- as.factor("Sclero")
head(Sclero_fire_df)



Shrub_fire_df <- as.data.frame(rbind(mean(Shrub_fire$Ground_based), mean(Shrub_fire$Satellite), mean(Shrub_fire$GLM)))
colnames(Shrub_fire_df) <- "fire_freq_mean"
Shrub_fire_df$Dataset <- as.factor(c("Ground_based", "Satellite", "GLM"))
Shrub_fire_df$sdv <- rbind(sd(Shrub_fire$Ground_based), sd(Shrub_fire$Satellite), sd(Shrub_fire$GLM))
Shrub_fire_df$Min_ff <- round(mean(Shrub_fire$min_ff_est))
Shrub_fire_df$Max_ff <- round(mean(Shrub_fire$max_ff_est))
Shrub_fire_df$Agg <- as.factor("Grassland & Shrubland")
head(Shrub_fire_df)




Other_fire_df <- as.data.frame(rbind(mean(Other_fire$Ground_based), mean(Other_fire$Satellite), mean(Other_fire$GLM)))
colnames(Other_fire_df) <- "fire_freq_mean"
Other_fire_df$Dataset <- as.factor(c("Ground_based", "Satellite", "GLM"))
Other_fire_df$sdv <- rbind(sd(Other_fire$Ground_based), sd(Other_fire$Satellite), sd(Other_fire$GLM))
Other_fire_df$Min_ff <- round(mean(Other_fire$min_ff_est))
Other_fire_df$Max_ff <- mean(Other_fire$max_ff_est)
Other_fire_df$Agg <- as.factor("Other")
head(Other_fire_df)
# The average maximum fire frequency estimate returns a very small number here which seems unusual
Other_fire$max_ff_est
# The data is zero-inflated and we can see that we have only three other values of "4". Let's manually calculate the average instead
Other_fire_df$Max_ff <- (0+4)/2
head(Other_fire_df)



RE_avg_fire <- rbind(Rainforest_fire_df, Sclero_fire_df, Shrub_fire_df, Other_fire_df)
RE_avg_fire
str(RE_avg_fire)
RE_avg_fire$Dataset <- factor(RE_avg_fire$Dataset, levels = c("Ground_based", "Satellite", "GLM"))

# For all BVG aggregations the minimum fire frequency estimate based on fire return interval guidelines for the past 35 years is 0. We won't draw this on the plot, just make note of this in the caption, we will only plot the maximum bar for each aggregation. Need to note the average maximum fire frequency estimate for the other BVG aggregate was manually calculated due to the zero-inflated data. 


# 7. Produce plot ----
# We have to manually plot the horizontal lines
ggplot(RE_avg_fire, aes(x = Agg, y = fire_freq_mean, fill = Dataset)) +
  geom_col(position = position_dodge()) +
  geom_textsegment(aes(x = 0.55, xend = 1.45, y = 0, yend = 0, label = "Recommended frequency"), colour = 'gray60', vjust = 1.2)+
  geom_segment(aes(x = 0.55, xend = 1.45, 0), colour = 'gray60', size = 1.5) +

  geom_textsegment(aes(x = 1.55, xend = 2.45, y = 9, yend = 9, label = "Max. recommended"), colour = 'gray60', vjust = 1.2)+
  geom_segment(aes(x = 1.55, xend = 2.45, 9), colour = 'gray60', size = 1.5)+
  geom_textsegment(aes(x = 1.55, xend = 2.45, y = 3, yend = 3, label = "Min. recommended"), colour = 'gray60', vjust = -0.2)+
  geom_segment(aes(x = 1.55, xend = 2.45, 3), colour = 'gray60', size = 1.5)+
  
  geom_textsegment(aes(x = 2.55, xend = 3.45, y = 9, yend = 9, label = "Max. recommended"), colour = 'gray60', vjust = 1.2)+
  geom_segment(aes(x = 2.55, xend = 3.45, 9), colour = 'gray60', size = 1.5)+
  geom_textsegment(aes(x = 2.55, xend = 3.45, y = 2, yend = 2, label = "Min. recommended"), colour = 'gray60', vjust = -0.2)+
  geom_segment(aes(x = 2.55, xend = 3.45, 2), colour = 'gray60', size = 1.5)+
  
  geom_textsegment(aes(x = 3.55, xend = 4.45, y= 2, yend = 2, label = "Max. recommended"), colour = 'gray60', vjust = -0.2)+
  geom_segment(aes(x = 3.55, xend = 4.45, y = 2), colour = 'gray60', size = 1.5) +
  geom_textsegment(aes(x = 3.55, xend = 4.45, y = 0, yend = 0, label = "Min. recommended"), colour = 'gray60', vjust = 1.2) +
  geom_segment(aes(x = 3.55, xend = 4.45, y = 0), colour = 'gray60', size = 1.5) +
    
  geom_errorbar(aes(ymin = ifelse(fire_freq_mean - sdv <0, 0, fire_freq_mean - sdv), ymax = round(fire_freq_mean + sdv)), position = position_dodge(0.9), width = 0.1)+
  scale_y_continuous(limits = c(0, 9.1), breaks = seq(0,9,1))+
  scale_fill_manual(values = c('gray80', 'steelblue', '#492050'), labels = c("Ground based", "Sentinel", "GLM"))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.title = element_text(size = 20, face = 'bold'),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 20),
        legend.position = "bottom")+
  labs(x = expression(bold("Broad vegetation aggregation")), y = expression(bold("Fire frequency")))+
  scale_x_discrete(labels = c("Rainforest \n BVG 1-7", "Sclerophyll \n BVG 8-27", "Grassland & Shrubland \n BVG 28-33", "Other \n BVG 34-35"))
  

  



save.image('./02_Workspaces/005_predictive_model_validations_vegetation.RData')
#load('./02_Workspaces/005_predictive_model_validations_vegetation.RData')

