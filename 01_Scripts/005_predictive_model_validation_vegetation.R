# Written by Felicity Charles
# Date: 06/10/2024
# Updated: 25/05/2026

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
unweighted_pred_fpc <- rast('./04_Results/Prediction_rasters/SEQ_IBRA_Unweighted_pred_FPC.tif')
down_wt_pred_fpc <- rast('./04_Results/Prediction_rasters/SEQ_IBRA_Downweighted_FPC_pred.tif')
IWLR_pred_fpc <- rast('./04_Results/Prediction_rasters/SEQ_IBRA_IWLR_FPC_pred.tif')
gam_pred_fpc <- rast('./04_Results/Prediction_rasters/SEQ_IBRA_GAM_FPC_pred.tif')
glm_pred_fpc <- rast('./04_Results/Prediction_rasters/SEQ_IBRA_GLM_FPC_pred.tif')


Sentinel_ff <- rast('./00_Data/Fire_data/Outputs/Sentinel/Sentinel_ff_hydrographical_mask_SEQ_IBRA_focal_cropped.tif')


environmental_preds <- rast('./00_Data/SDM_data/predictors.tif')
#BVG <- vect('./00_Data/Environmental_data/Remnant_2021_broad_veg_groups/Remnant_broad_vegetation_groups.shp') %>%
#project('EPSG:3577') %>% 
#crop(gam_pred_fpc)
#RE <- vect('./00_Data/Environmental_data/Regional_ecosystem/Biodiversity_status_of_remnant_regional_ecosystems.shp') %>% 
#project('EPSG:3577') %>% 
#crop(gam_pred_fpc)

QPWS_rand <- vect('./00_Data/Fire_data/Outputs/QPWS_random.gpkg')
QPWS_ff <- rast('./00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_freq_hydrographical_mask_cropped_reproj.tif')
SEQ <- vect('./00_Data/SEQ_bound/SEQ_IBRA.gpkg')



# Further ground truthing/model validation can be performed to ensure the predictions seem realistic considering the broad vegetation group, regional ecosystem fire regime management suggestions, and land use.

# 3. Compare estimates of fire frequency based on fire management guidelines to the predictions ----
# Here we will compare fire frequency estimates to the actual values recorded for both QPWS and Sentinel, and also determine how our predicted values compare to fire management guidelines. 
fire_reg <- vect('./00_Data/Environmental_data/DP_SEQ_VEG_FIRE_RGM_100K_A/SEQ_VEG_FIRE_RGM_100K_A.shp')

fireg <- fire_reg[!is.na(fire_reg$FREQUENCY),]
table(fireg$FREQUENCY, fireg$FVG)

head(fireg [5,])

table(fireg$FREQUENCY, fireg$FVG)
# Ecosystems lacking frequency information are:

# Mangroves and saltmarsh - do not burn intentionally or no fire
# Rainforests - no fire
# Riparian - no fire
# Wet tall open - no fire
# Plantations - fire frequency as required

table(fireg$FREQUENCY == "Do not intentionally burn", fireg$FVG)
table(fireg$FREQUENCY == "As required", fireg$FVG)


fireg <- subset(fireg, fireg$FVG != "Water" & fireg$FVG != "Sand" & fireg$FVG != "Non-remnant" & fireg$FVG != "Plantations")
unique(fireg$FVG)

# Assign negative values for do not intentionally burn and as required and 0 for no fire
fireg$FREQUENCY <- ifelse(fireg$FREQUENCY == "Do not intentionally burn", -1, fireg$FREQUENCY)
fireg$FREQUENCY <- ifelse(fireg$FREQUENCY == "As required", -2, fireg$FREQUENCY)
fireg$FREQUENCY <- ifelse(fireg$FREQUENCY == "No fire", 0 , fireg$FREQUENCY)
unique(fireg$FREQUENCY)





# Calculate minimum and maximum fire frequency for 36 years
# Extract minimum interval (number before dash, 0 for non-numeric)
fireg$min_int <- ifelse(grepl("\\d+-\\d+", fireg$FREQUENCY),
                        as.numeric(sub("(\\d+)-.*", "\\1", fireg$FREQUENCY)),
                        0)

# Extract maximum interval (number after dash, 0 for non-numeric)  
fireg$max_int <- ifelse(grepl("\\d+-\\d+", fireg$FREQUENCY),
                        as.numeric(sub(".*-(\\d+).*", "\\1", fireg$FREQUENCY)),
                        0)

fireg$max_ff <- ifelse(round(36/fireg$min_int) == "Inf", 0 , round(36/fireg$min_int))
fireg$min_ff <- ifelse(round(36/fireg$max_int) == "Inf", 0 , round(36/fireg$max_int))



# Now we want to get the fire frequency information for each dataset for the polygons of fireg
fireg <- project(fireg, 'EPSG:3577')


# Create 1000 points for each FVG
# Will rename FVGs as follows
# Rainforest, vine forest, and brigalow
# Open forests/woodlands
# Grasslands
# Coastal forests and headlands
# Melaleuca communities
# Heath communities
# Mangroves and saltmarsh
# Riparian, foredune, beach ridge communities
# Wet tall open forests

set.seed(42); Rain <- spatSample(fireg[fireg$FVG == "Rainforests, dry vine forests and brigalow communities"],  1000, 'random'); head(Rain); dim(Rain)
set.seed(42); open <- spatSample(fireg[fireg$FVG == "Open forests/woodlands"],  1000, 'random'); dim(open)
set.seed(42); grass <- spatSample(fireg[fireg$FVG == "Grasslands"],  15580, 'random'); dim(grass)
set.seed(42); coast <- spatSample(fireg[fireg$FVG == "Coastal fringing forests and headlands"],  30690, 'random'); dim(coast);
set.seed(42); mel <- spatSample(fireg[fireg$FVG == "Melaleuca communities"],  1540, 'random'); dim(mel)
set.seed(42); heath <- spatSample(fireg[fireg$FVG == "Heath communities"],  1000, 'random'); dim(heath)
set.seed(42); mang <- spatSample(fireg[fireg$FVG == "Mangroves and saltmarsh"],  1370, 'random'); dim(mang)
mang <- mang[1:1000,]; dim(mang) # Unable to only produce 1000 so just remove the extra point
set.seed(42); rip <- spatSample(fireg[fireg$FVG == "Riparian, foredune, coral cay island and beach ridge communities"],  1291, 'random'); dim(rip)
set.seed(42); wet <- spatSample(fireg[fireg$FVG == "Wet tall open forests"],  1279, 'random'); dim(wet)

Rain$FVG_name <- "Rainforests, vine forests & brigalow"
open$FVG_name <- open$FVG
grass$FVG_name <- grass$FVG
coast$FVG_name <- "Coastal forests & headlands"
mel$FVG_name <- mel$FVG
heath$FVG_name <- heath$FVG
mang$FVG_name <- "Mangroves & saltmarsh"
rip$FVG_name <- "Riparian, foredune & beach ridges"
wet$FVG_name <- wet$FVG
head(wet)

FVGrec_pts <- rbind(open, grass, coast, mel, heath)


# Extract predictions for Satellite data -----
FVG_recs_sent <- extract(Sentinel_ff, FVGrec_pts)
colnames(FVG_recs_sent) <-  c('ID', 'Fire_frequency')
FVG_recs_sent$Dataset <- "Satellite"
FVG_recs_sent$frequency_status <- ifelse(!is.na(FVG_recs_sent$Fire_frequency),
                                         ifelse(FVG_recs_sent$Fire_frequency > FVGrec_pts$max_ff, "Higher",
                                                ifelse(FVG_recs_sent$Fire_frequency < FVGrec_pts$min_ff, "Lower",
                                                       "Within")),
                                         NA)



rain_sent <- extract(Sentinel_ff, Rain)
colnames(rain_sent) <- c('ID', 'Fire_frequency')
rain_sent$Dataset <- "Satellite"
rain_sent$frequency_status <- "No fire"
rain_sent <- as.data.frame(cbind(Rain[ , c(26:28, ncol(Rain))], rain_sent))

mang_sent <- extract(Sentinel_ff, mang)
colnames(mang_sent) <- c('ID', 'Fire_frequency')
mang_sent$Dataset <- 'Satellite'
mang_sent$frequency_status <- ifelse(!is.na(mang_sent$Fire_frequency), 
                                     ifelse(mang$FREQUENCY == -1, 
                                            'Do not intentionally burn', 
                                            'No fire'), 
                                     'No fire')
mang_sent <- as.data.frame(cbind(mang[ , c(26:28, ncol(mang))], mang_sent))

rip_sent <- extract(Sentinel_ff, rip)
colnames(rip_sent) <- c('ID', 'Fire_frequency')
rip_sent$Dataset <- 'Satellite'
rip_sent$frequency_status <- "No fire"
rip_sent <- as.data.frame(cbind(rip[ , c(26:28, ncol(rip))], rip_sent))


wet_sent <- extract(Sentinel_ff, wet)
colnames(wet_sent) <- c('ID', 'Fire_frequency')
wet_sent$Dataset <- 'Satellite'
wet_sent$frequency_status <- 'As required'
wet_sent <- as.data.frame(cbind(wet[ , c(26:28, ncol(wet))], wet_sent))

# Extract predictions for public land data ----
FVG_recs_pub <- extract(QPWS_ff, FVGrec_pts)
colnames(FVG_recs_pub) <-  c('ID', 'Fire_frequency')
FVG_recs_pub$Dataset <- "Public land"
FVG_recs_pub$frequency_status <- ifelse(!is.na(FVG_recs_pub$Fire_frequency),
                                        ifelse(FVG_recs_pub$Fire_frequency > FVGrec_pts$max_ff, "Higher",
                                               ifelse(FVG_recs_pub$Fire_frequency < FVGrec_pts$min_ff, "Lower",
                                                      "Within")),
                                        NA)


rain_pub <- extract(QPWS_ff, Rain)
colnames(rain_pub) <- c('ID', 'Fire_frequency')
rain_pub$Dataset <- "Public land"
rain_pub$frequency_status <- "No fire"
rain_pub <- as.data.frame(cbind(Rain[ , c(26:28, ncol(Rain))], rain_pub))

mang_pub <- extract(QPWS_ff, mang)
colnames(mang_pub) <- c('ID', 'Fire_frequency')
mang_pub$Dataset <- 'Public land'
mang_pub$frequency_status <- ifelse(!is.na(mang_pub$Fire_frequency), 
                                    ifelse(mang$FREQUENCY == -1, 
                                           'Do not intentionally burn', 
                                           'No fire'), 
                                    'No fire')
mang_pub <- as.data.frame(cbind(mang[ , c(26:28, ncol(mang))], mang_pub))

rip_pub <- extract(QPWS_ff, rip)
colnames(rip_pub) <- c('ID', 'Fire_frequency')
rip_pub$Dataset <- 'Public land'
rip_pub$frequency_status <- "No fire"
rip_pub <- as.data.frame(cbind(rip[ , c(26:28, ncol(rip))], rip_pub))


wet_pub <- extract(QPWS_ff, wet)
colnames(wet_pub) <- c('ID', 'Fire_frequency')
wet_pub$Dataset <- 'Public land'
wet_pub$frequency_status <- 'As required'
wet_pub <- as.data.frame(cbind(wet[ , c(26:28, ncol(wet))], wet_pub))


# Extract predictions for GAM FPC data -----
FVG_recs_gam_fpc <- extract(gam_pred_fpc, FVGrec_pts)
colnames(FVG_recs_gam_fpc) <-  c('ID', 'Fire_frequency')
FVG_recs_gam_fpc$Dataset <- "GAM"
FVG_recs_gam_fpc$frequency_status <- ifelse(!is.na(FVG_recs_gam_fpc$Fire_frequency),
                                            ifelse(FVG_recs_gam_fpc$Fire_frequency > FVGrec_pts$max_ff, "Higher",
                                                   ifelse(FVG_recs_gam_fpc$Fire_frequency < FVGrec_pts$min_ff, "Lower",
                                                          "Within")),
                                            NA)



rain_gam_fpc <- extract(gam_pred_fpc, Rain)
colnames(rain_gam_fpc) <- c('ID', 'Fire_frequency')
rain_gam_fpc$Dataset <- "GAM"
rain_gam_fpc$frequency_status <- "No fire"
rain_gam_fpc <- as.data.frame(cbind(Rain[ , c(26:28, ncol(Rain))], rain_gam_fpc))

mang_gam_fpc <- extract(gam_pred_fpc, mang)
colnames(mang_gam_fpc) <- c('ID', 'Fire_frequency')
mang_gam_fpc$Dataset <- 'GAM'
mang_gam_fpc$frequency_status <- ifelse(!is.na(mang_gam_fpc$Fire_frequency), 
                                        ifelse(mang$FREQUENCY == -1, 
                                               'Do not intentionally burn', 
                                               'No fire'), 
                                        'No fire')
mang_gam_fpc <- as.data.frame(cbind(mang[ , c(26:28, ncol(mang))], mang_gam_fpc))

rip_gam_fpc <- extract(gam_pred_fpc, rip)
colnames(rip_gam_fpc) <- c('ID', 'Fire_frequency')
rip_gam_fpc$Dataset <- 'GAM'
rip_gam_fpc$frequency_status <- "No fire"
rip_gam_fpc <- as.data.frame(cbind(rip[ , c(26:28, ncol(rip))], rip_gam_fpc))


wet_gam_fpc <- extract(gam_pred_fpc, wet)
colnames(wet_gam_fpc) <- c('ID', 'Fire_frequency')
wet_gam_fpc$Dataset <- 'GAM'
wet_gam_fpc$frequency_status <- 'As required'
wet_gam_fpc <- as.data.frame(cbind(wet[ , c(26:28, ncol(wet))], wet_gam_fpc))




# Extract predictions for GLM FPC data -----
FVG_recs_glm_fpc <- extract(glm_pred_fpc, FVGrec_pts)
colnames(FVG_recs_glm_fpc) <-  c('ID', 'Fire_frequency')
FVG_recs_glm_fpc$Dataset <- "GLM"
FVG_recs_glm_fpc$frequency_status <- ifelse(!is.na(FVG_recs_glm_fpc$Fire_frequency),
                                            ifelse(FVG_recs_glm_fpc$Fire_frequency > FVGrec_pts$max_ff, "Higher",
                                                   ifelse(FVG_recs_glm_fpc$Fire_frequency < FVGrec_pts$min_ff, "Lower",
                                                          "Within")),
                                            NA)



rain_glm_fpc <- extract(glm_pred_fpc, Rain)
colnames(rain_glm_fpc) <- c('ID', 'Fire_frequency')
rain_glm_fpc$Dataset <- "GLM"
rain_glm_fpc$frequency_status <- "No fire"
rain_glm_fpc <- as.data.frame(cbind(Rain[ , c(26:28, ncol(Rain))], rain_glm_fpc))

mang_glm_fpc <- extract(glm_pred_fpc, mang)
colnames(mang_glm_fpc) <- c('ID', 'Fire_frequency')
mang_glm_fpc$Dataset <- 'GLM'
mang_glm_fpc$frequency_status <- ifelse(!is.na(mang_glm_fpc$Fire_frequency), 
                                        ifelse(mang$FREQUENCY == -1, 
                                               'Do not intentionally burn', 
                                               'No fire'), 
                                        'No fire')
mang_glm_fpc <- as.data.frame(cbind(mang[ , c(26:28, ncol(mang))], mang_glm_fpc))

rip_glm_fpc <- extract(glm_pred_fpc, rip)
colnames(rip_glm_fpc) <- c('ID', 'Fire_frequency')
rip_glm_fpc$Dataset <- 'GLM'
rip_glm_fpc$frequency_status <- "No fire"
rip_glm_fpc <- as.data.frame(cbind(rip[ , c(26:28, ncol(rip))], rip_glm_fpc))


wet_glm_fpc <- extract(glm_pred_fpc, wet)
colnames(wet_glm_fpc) <- c('ID', 'Fire_frequency')
wet_glm_fpc$Dataset <- 'GLM'
wet_glm_fpc$frequency_status <- 'As required'
wet_glm_fpc <- as.data.frame(cbind(wet[ , c(26:28, ncol(wet))], wet_glm_fpc))







# Combine datasets for plotting
FVGs_withrecs_glm_fpc <- cbind(FVGrec_pts[,c(26:28, ncol(FVGrec_pts))], FVG_recs_glm_fpc)
FVGs_withrecs_gam_fpc <- cbind(FVGrec_pts[,c(26:28, ncol(FVGrec_pts))], FVG_recs_gam_fpc)
FVGs_withrecs_sent <- cbind(FVGrec_pts[,c(26:28, ncol(FVGrec_pts))], FVG_recs_sent)
FVGs_withrecs_pub <- cbind(FVGrec_pts[,c(26:28, ncol(FVGrec_pts))], FVG_recs_pub)
FVGs_withrecs <- rbind(FVGs_withrecs_glm_fpc, FVGs_withrecs_gam_fpc, FVGs_withrecs_sent, FVGs_withrecs_pub)
head(FVGs_withrecs)

# Reorder the frequency status 
FVGs_withrecs$frequency_status <- factor(FVGs_withrecs$frequency_status,
                                         levels = c("Lower", "Within", "Higher"))


# Reorder FVG_name factor
FVGs_withrecs$FVG_name <- factor(FVGs_withrecs$FVG_name, 
                                 levels = c("Open forests/woodlands",
                                            "Melaleuca communities", 
                                            "Heath communities",
                                            "Grasslands",
                                            "Coastal forests & headlands"))

FVGs_withrecs_df <- as.data.frame(FVGs_withrecs)
FVGs_withrecs_df$Dataset <- as.factor(FVGs_withrecs_df$Dataset)
str(FVGs_withrecs_df)
FVGs_withrecs_df$Dataset <- factor(FVGs_withrecs_df$Dataset, levels = c("Public land", "Satellite", "GAM", "GLM"))



Mangrove <- rbind(mang_pub, mang_sent, mang_glm_fpc, mang_gam_fpc)
Riparian <- rbind(rip_pub, rip_sent, rip_glm_fpc, rip_gam_fpc)
Wet <- rbind(wet_pub, wet_sent, wet_glm_fpc, wet_gam_fpc)
Rainforest <- rbind(rain_pub, rain_sent, rain_glm_fpc, rain_gam_fpc)

FVGs_nrecs <- rbind(Mangrove, Riparian, Wet, Rainforest)
head(FVGs_nrecs)
str(FVGs_nrecs)
unique(FVGs_nrecs$frequency_status)
FVGs_nrecs_df <- as.data.frame(FVGs_nrecs)
FVGs_nrecs_df$Dataset <- as.factor(FVGs_nrecs_df$Dataset)
str(FVGs_nrecs_df)
FVGs_nrecs_df$Dataset <- factor(FVGs_nrecs_df$Dataset, levels = c("Public land", "Satellite", "GAM", "GLM"))




library(ggridges); library(ggplot2)

# Labels for facetting
Agg_labs <- c("Open forests/woodlands" = "(a) Open forests/\nwoodlands",
              "Melaleuca communities" = "(b) Melaleuca\ncommunities",
              "Heath communities" = "(c) Heath\ncommunities",
              "Grasslands" = "(d) Grasslands",
              "Coastal forests & headlands" = "(e) Coastal forests\n& headlands",
              "Mangroves & saltmarsh" = "(f) Mangroves & saltmarsh",
              "Riparian, foredune & beach ridges" = "(g) Riparian, foredune & beach ridges",
              "Wet tall open forests" = "(h) Wet tall open forests",
              "Rainforests, vine forests & brigalow" = "(i) Rainforests, vine forests & brigalow")



p1 <- ggplot(data = FVGs_withrecs_df[!is.na(FVGs_withrecs_df$Fire_frequency),], 
             aes(x = Fire_frequency, 
                 y = factor(frequency_status, levels = c("Lower", 
                                                         "Within", 
                                                         "Higher")), 
                 fill = Dataset,
                 colour = Dataset)) +
  geom_density_ridges(alpha = 0.4, scale = 0.95, stat = 'binline', breaks = seq(0,25,1),draw_baseline = F, linewidth = 0.5) + # Plot as a histogram with 31 bins (1 bin per fire frequency), set scale as <1 to ensure histograms do not overlap, do not draw baseline to remove trailing lines 
  facet_wrap(vars(FVG_name), scales = "free", labeller = as_labeller(Agg_labs)) + # Allow x axis to vary between communities based on the highest fire frequency for that community, allowing easier interrogation
  scale_fill_manual(values = c("gray70", "steelblue", "#492050", "#AAA970"))+
  scale_colour_manual(values = c("gray70", "steelblue", "#492050", "#AAA970"))+
  scale_y_discrete(drop = FALSE, expand = c(0.01,0)) +  # This ensures all factor levels are shown
  scale_x_continuous(
    breaks = function(x) {
      max_val <- min(ceiling(max(x)/5) * 5, 25)  # Cap at 25
      seq(0, max_val, by = 5)
    },
    limits = function(x) {
      max_val <- min(ceiling(max(x)/5) * 5, 25)  # Cap at 25
      c(0, max_val)
    },
    expand = c(0.01,0)) +  
  labs(x = "Fire Frequency", y = "Fire regime status count")+
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(colour = 'black'),
        strip.text = element_text(face = "bold", size = 20),
        legend.title = element_blank(),
        axis.title = element_text(face = 'bold', size = 20),
        axis.text.x = element_text(vjust = 0),
        axis.text.y = element_text(vjust = 0),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 18),
        legend.position = c(0.8, 0.3))
p1
#1300x851


filt_FVGs_nrecs_df <- FVGs_nrecs_df[FVGs_nrecs_df$Fire_frequency != 0, ]
table(filt_FVGs_nrecs_df$FVG_name, round(filt_FVGs_nrecs_df$Fire_frequency))

p2 <- ggplot(data = FVGs_nrecs_df[!is.na(FVGs_nrecs_df$Fire_frequency),], 
             aes(x = Fire_frequency, 
                 y = factor(frequency_status, levels = c("No fire", 
                                                         "Do not intentionally burn", 
                                                         "As required")), 
                 fill = Dataset, colour = Dataset)) +
  geom_density_ridges(alpha = 0.6, scale = 0.95, stat = 'binline', breaks = seq(1,25,1), draw_baseline = F) + # Plat as a histogram with 21 bins (1 bin per fire frequency), set scael as <1 to ensure histograms do not overlap, do not draw baseline to remove trailing lines 
  facet_wrap(vars(FVG_name), scales = "free", labeller = as_labeller(Agg_labs), drop = T, nrow = 2) + # Allow x axis to vary between communities based on the highest fire frequency for that commubnity, allowing easier interrogation
  scale_fill_manual(values = c("gray80", "steelblue", "#492050", "#AAA970"), labels = c("Public land", "Satellite", "GLM", "GAM"))+
  scale_colour_manual(values = c("gray80", "steelblue", "#492050", "#AAA970"), labels = c("Public land", "Satellite", "GLM", "GAM"))+
  scale_y_discrete(expand = c(0.01,0), labels = c("No fire", "Do not\nintentionally\nburn", "As required")) +  
  scale_x_continuous(
    breaks = c(1, 5, 10, 15, 20),
    limits = function(x) {
      max_val <- min(ceiling(max(x)/5) * 5, 30)
      c(1, max_val)
    },
    expand = c(0.01,0)) +  
  labs(x = "Fire Frequency", y = "Fire regime status count")+
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(colour = 'black'),
        strip.text = element_text(face = "bold", size = 20),
        legend.position = "none",
        axis.title = element_text(face = 'bold', size = 20),
        axis.text.x = element_text(vjust = 0),
        axis.text.y = element_text(vjust = 0),
        axis.text = element_text(size = 14)) 
p2
#1400 x 1100


# Create an inset plot to show the distribution of zeros for these vegetation types 
veg_types <- unique(FVGs_nrecs_df$FVG) # Get vegetation types
FVGs_nrecs_0s <- FVGs_nrecs_df[FVGs_nrecs_df$Fire_frequency == 0 | is.na(FVGs_nrecs_df$Fire_frequency), ] # Get the zeros
FVGs_nrecs_0s$Fire_frequency[is.na(FVGs_nrecs_0s$Fire_frequency)] <- 0
# Check the data
table(FVGs_nrecs_0s$FVG, FVGs_nrecs_0s$Dataset)

riparian <- ggplot(data = FVGs_nrecs_0s[FVGs_nrecs_0s$FVG == "Riparian, foredune, coral cay island and beach ridge communities", ],
                   aes(x = 1, fill = Dataset, colour = Dataset)) +  # Changed x to constant 1
  geom_bar(stat = "count", position = position_dodge(width = 0.1), width = 0.2)+
  scale_fill_manual(values = c("Public land" = "gray80", "Satellite" = "steelblue", 
                               "GLM" = "#492050", "GAM" = "#AAA970"))+
  scale_colour_manual(values = c("Public land" = "gray80", "Satellite" = "steelblue", 
                                 "GLM" = "#492050", "GAM" = "#AAA970"))+
  scale_y_continuous(limits = c(0,1000))+
  labs(y = "Count") +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(colour = 'black'),
        strip.text = element_text(face = "bold", size = 20),
        legend.position = "none",
        axis.title.y = element_text(face = 'bold', size = 20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        axis.text.y = element_text(vjust = 0),
        axis.text = element_text(size = 14))
riparian

wet <- ggplot(data = FVGs_nrecs_0s[FVGs_nrecs_0s$FVG == "Wet tall open forests", ],
              aes(x = 1, fill = Dataset, colour = Dataset)) +  # Changed x to constant 1
  geom_bar(stat = "count", position = position_dodge(width = 0.1), width = 0.2)+
  scale_fill_manual(values = c("Public land" = "gray80", "Satellite" = "steelblue", 
                               "GLM" = "#492050", "GAM" = "#AAA970"))+
  scale_colour_manual(values = c("Public land" = "gray80", "Satellite" = "steelblue", 
                                 "GLM" = "#492050", "GAM" = "#AAA970"))+
  scale_y_continuous(limits = c(0,1000))+
  labs(y = "Count") +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(colour = 'black'),
        strip.text = element_text(face = "bold", size = 20),
        legend.position = "none",
        axis.title.y = element_text(face = 'bold', size = 20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        axis.text.y = element_text(vjust = 0),
        axis.text = element_text(size = 14))
wet

rain <- ggplot(data = FVGs_nrecs_0s[FVGs_nrecs_0s$FVG == "Rainforests, dry vine forests and brigalow communities", ],
               aes(x = 1, fill = Dataset, colour = Dataset)) +  # Changed x to constant 1
  geom_bar(stat = "count", position = position_dodge(width = 0.1), width = 0.2)+
  scale_fill_manual(values = c("Public land" = "gray80", "Satellite" = "steelblue", 
                               "GLM" = "#492050", "GAM" = "#AAA970"))+
  scale_colour_manual(values = c("Public land" = "gray80", "Satellite" = "steelblue", 
                                 "GLM" = "#492050", "GAM" = "#AAA970"))+
  scale_y_continuous(limits = c(0,1000))+
  labs(y = "Count") +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(colour = 'black'),
        strip.text = element_text(face = "bold", size = 20),
        legend.position = "none",
        axis.title.y = element_text(face = 'bold', size = 20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        axis.text.y = element_text(vjust = 0),
        axis.text = element_text(size = 14))
rain

mang <- ggplot(data = FVGs_nrecs_0s[FVGs_nrecs_0s$FVG == "Mangroves and saltmarsh", ],
               aes(x = 1, fill = Dataset, colour = Dataset)) +  # Changed x to constant 1
  geom_bar(stat = "count", position = position_dodge(width = 0.1), width = 0.2)+
  scale_fill_manual(values = c("Public land" = "gray80", "Satellite" = "steelblue", 
                               "GLM" = "#492050", "GAM" = "#AAA970"))+
  scale_colour_manual(values = c("Public land" = "gray80", "Satellite" = "steelblue", 
                                 "GLM" = "#492050", "GAM" = "#AAA970"))+
  labs(y = "Count") +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(colour = 'black'),
        strip.text = element_text(face = "bold", size = 20),
        legend.position = "none",
        axis.title.y = element_text(face = 'bold', size = 20),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        axis.text.y = element_text(vjust = 0),
        axis.text = element_text(size = 14))
mang



p2_with_insets <- ggdraw(p2) +
  # Mangroves (top left) - adjust coordinates as needed
  draw_plot(mang, x = 0.39, y = 0.73, width = 0.12, height = 0.2) +
  # Rainforest (top right)
  draw_plot(rain, x = 0.87, y = 0.73, width = 0.12, height = 0.2) +
  # Riparian (bottom left)
  draw_plot(riparian, x = 0.39, y = 0.25, width = 0.12, height = 0.2) +
  # Wet tall (bottom right)
  draw_plot(wet, x = 0.87, y = 0.25, width = 0.12, height = 0.2)
p2_with_insets
ggsave("./04_Results/Plots/Vegetation_validation/Sub_FPC_insets.png", p2_with_insets, width = 14, height = 11, dpi = 300)

######## NOW WE WANT TO PLOT P2 under P1










##### OLD 



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
