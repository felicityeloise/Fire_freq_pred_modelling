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
Other$Agg <- "Wetland, Mangrove, Saltmarsh"


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
names(QPWS_ff) <- "Fire_freq"
names(Sentinel_ff) <- "Fire_freq"
names(glm_pred) <- "Fire_freq"

glm_pred <- round(glm_pred)
glm_pred[is.na(glm_pred)] <- 0 


# Extract information
QPWS_pts <- extract(QPWS_ff, RE_pts)
QPWS_pts$Dataset <- "Public land"
Sent_pts <- extract(Sentinel_ff, RE_pts)
Sent_pts$Dataset <- "Satellite"
glm_pts <- extract(glm_pred, RE_pts)
glm_pts$Dataset <- "GLM"

RE_ext <- extract(RE, RE_pts) # This contains the regional ecosystem, BVG, and estimated maximum and minimum fire frequencies
head(RE_ext)
RE_ext_subset <- RE_ext[, c(1,3,29,32:ncol(RE_ext))]
head(RE_ext_subset)
colnames(RE_ext_subset) <- c("ID", "re1", "dbvg5m", 'interval_min', "interval_max", "max_ff_est", "min_ff_est")
unique(RE_ext_subset$re1)



# Now we are going to combine each fire dataset with the RE information, then this will be stacked so that we have one column for fire frequency
RE_pts_df <- as.data.frame(RE_pts[, 35])
RE_pts_df$ID <- 1:4000
head(RE_pts_df)
RE_ext_subset <- left_join(RE_ext_subset, RE_pts_df, by = 'ID')
head(RE_ext_subset)

QPWS_RE <- left_join(RE_ext_subset, QPWS_pts, by = 'ID')
head(QPWS_RE)

Sent_RE <- left_join(RE_ext_subset, Sent_pts, by = "ID")
head(Sent_RE)

GLM_RE <- left_join(RE_ext_subset, glm_pts, by = "ID")
head(GLM_RE)


RE_randpt_fire <- rbind(QPWS_RE, Sent_RE, GLM_RE)
head(RE_randpt_fire); tail(RE_randpt_fire); dim(RE_randpt_fire)
unique(is.na(RE_randpt_fire))

# Overwrite NA values for fire frequency
RE_randpt_fire$Fire_freq[is.na(RE_randpt_fire$Fire_freq)] <- 0
unique(is.na(RE_randpt_fire))
str(RE_randpt_fire)

RE_randpt_fire$Dataset <- factor(RE_randpt_fire$Dataset, levels = c("GLM", "Public land", "Satellite"))
RE_randpt_fire$Dataset <- factor(RE_randpt_fire$Dataset, levels = c("Public land", "Satellite", "GLM"))

str(RE_randpt_fire)
unique(RE_randpt_fire$Fire_freq)
RE_randpt_fire$Fire_freq <- round(RE_randpt_fire$Fire_freq)


RE_randpt_fire$Agg <- factor(RE_randpt_fire$Agg)
RE_randpt_fire$Agg <- factor(RE_randpt_fire$Agg, levels = c("Rainforest", "Sclerophyll", "Grassland & Shrubland", "Wetland, Mangrove, Saltmarsh"))


# Labels for facetting, to add some extra information 
Agg_labs <- c("Rainforest" = "(a) Rainforest \n BVG 1-7",
              "Sclerophyll" = "(b) Sclerophyll \n BVG 8-27",
              "Grassland & Shrubland" = '(c) Grassland & Shrubland \n BVG 28-33',
            "Wetland, Mangrove, Saltmarsh" = "(d) Wetland, Mangrove, Saltmarsh \n BVG 34-35")

ggplot(data = RE_randpt_fire, aes(x = Fire_freq, fill = Dataset))+
  geom_histogram(bins = 11, position = position_dodge(0.5))+
  scale_fill_manual(values = c("gray80", "steelblue", "#492050"), labels = c("Public", "Satellite", "GLM"))+
  theme_bw()+
  scale_y_continuous(expression(bold("Count")), limits = c(0,900), breaks = seq(0,900,100)) +
  scale_x_continuous(expression(bold("Fire frequency")), breaks = seq(0,10,1)) +
  facet_wrap(vars(RE_randpt_fire$Agg), labeller = as_labeller(Agg_labs), scales = 'fixed', axes = "all")+
  theme(axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = 'bold', size = 20),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.position = "bottom")


save.image('./02_Workspaces/005_predictive_model_validations_vegetation.RData')
#load('./02_Workspaces/005_predictive_model_validations_vegetation.RData')


# Lets do this again but add predictions from the GAM as well
names(gam_pred) <- "Fire_freq"

gam_pred <- round(gam_pred)
gam_pred[is.na(gam_pred)] <- 0 


# Extract information
GAM_pts <- extract(gam_pred, RE_pts)
GAM_pts$Dataset <- "GAM"



GAM_RE <- left_join(RE_ext_subset, GAM_pts, by = "ID")
head(GAM_RE)


RE_randpt_fire <- rbind(QPWS_RE, Sent_RE, GLM_RE, GAM_RE)
head(RE_randpt_fire); tail(RE_randpt_fire); dim(RE_randpt_fire)
unique(is.na(RE_randpt_fire))

# Overwrite NA values for fire frequency
RE_randpt_fire$Fire_freq[is.na(RE_randpt_fire$Fire_freq)] <- 0
unique(is.na(RE_randpt_fire))
str(RE_randpt_fire)

RE_randpt_fire$Dataset <- factor(RE_randpt_fire$Dataset, levels = c("GAM", "GLM", "Public land", "Satellite"))
RE_randpt_fire$Dataset <- factor(RE_randpt_fire$Dataset, levels = c("Public land", "Satellite", "GLM", "GAM"))

str(RE_randpt_fire)
unique(RE_randpt_fire$Fire_freq)
RE_randpt_fire$Fire_freq <- round(RE_randpt_fire$Fire_freq)


RE_randpt_fire$Agg <- factor(RE_randpt_fire$Agg)
RE_randpt_fire$Agg <- factor(RE_randpt_fire$Agg, levels = c("Rainforest", "Sclerophyll", "Grassland & Shrubland", "Wetland, Mangrove, Saltmarsh"))



ggplot(data = RE_randpt_fire, aes(x = Fire_freq, fill = Dataset))+
  geom_histogram(bins = 11, position = position_dodge(0.7))+
  scale_fill_manual(values = c("gray80", "steelblue", "#492050", "#AAA970"), labels = c("Public land", "Satellite", "GLM", "GAM"))+
  theme_bw()+
  scale_y_continuous(expression(bold("Count of cells")), limits = c(0,900), breaks = seq(0,900,100)) +
  scale_x_continuous(expression(bold("Fire frequency")), breaks = seq(0,10,1)) +
  facet_wrap(vars(RE_randpt_fire$Agg), labeller = as_labeller(Agg_labs), scales = 'fixed', axes = "all")+
  theme(axis.line = element_line(colour = 'black'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = 'bold', size = 20),
        axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.position = "bottom")

# 1100 x 1000
save.image('./02_Workspaces/005_predictive_model_validations_vegetation.RData')

table(RE_randpt_fire$Dataset, RE_randpt_fire$Fire_freq)
