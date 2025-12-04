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
library(ggridges) # ggridges_0.5.7
library(cowplot)


# 2. Load original data, predictive model data, and environmental data ----
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





# 4. Calculate minimum and maximum fire frequency for 36 years -----
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





# 5. Create 1000 points for each FVG ----
fireg <- project(fireg, 'EPSG:3577')

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
head(wet);dim(wet)

FVGrec_pts <- rbind(open, grass, coast, mel, heath)

# 6. Extract fire frequencies for each FVG ----
# 6.1 Extract predictions for Satellite data
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

# 6.2 Extract predictions for public land data ----
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


# 6.3 Extract predictions for GAM FPC data -----
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




# 6.4 Extract predictions for GLM FPC data -----
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







# 6.5 Combine fire frequency FVG information for each dataset for plotting -----
# 6.5.1 FVGs with recommendations ----
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


# 6.5.2 FVGs without recommendations ----
Mangrove <- rbind(mang_pub, mang_sent, mang_glm_fpc, mang_gam_fpc); dim (Mangrove)
Riparian <- rbind(rip_pub, rip_sent, rip_glm_fpc, rip_gam_fpc); dim (Riparian)
Wet <- rbind(wet_pub, wet_sent, wet_glm_fpc, wet_gam_fpc); dim (Wet)
Rainforest <- rbind(rain_pub, rain_sent, rain_glm_fpc, rain_gam_fpc); dim(Rainforest)

FVGs_nrecs <- rbind(Mangrove, Riparian, Wet, Rainforest)
head(FVGs_nrecs)
str(FVGs_nrecs)
unique(FVGs_nrecs$frequency_status)
FVGs_nrecs_df <- as.data.frame(FVGs_nrecs)
FVGs_nrecs_df$Dataset <- as.factor(FVGs_nrecs_df$Dataset)
str(FVGs_nrecs_df)
FVGs_nrecs_df$Dataset <- factor(FVGs_nrecs_df$Dataset, levels = c("Public land", "Satellite", "GAM", "GLM"))


# 7. Create plots ----
# 7.1 Create labels for facets
Agg_labs <- c("Open forests/woodlands" = "(a) Open forests/\nwoodlands",
              "Melaleuca communities" = "(b) Melaleuca\ncommunities",
              "Heath communities" = "(c) Heath\ncommunities",
              "Grasslands" = "(d) Grasslands",
              "Coastal forests & headlands" = "(e) Coastal forests\n& headlands",
              "Mangroves & saltmarsh" = "(f) Mangroves\n & saltmarsh",
              "Riparian, foredune & beach ridges" = "(h) Riparian, foredune\n & beach ridges",
              "Wet tall open forests" = "(i) Wet tall\n open forests",
              "Rainforests, vine forests & brigalow" = "(g) Rainforests, vine\n forests & brigalow")


# 7.2 Create main plot for vegetation with recommendations
# Calculate maximum cell count for each FVG based on the frequency status and dataset 
dat <- FVGs_withrecs_df[!is.na(FVGs_withrecs_df$Fire_frequency),]
counts_by_dataset <- aggregate(Fire_frequency ~ FVG_name + frequency_status + Dataset, 
                               data = dat, 
                               FUN = length)
names(counts_by_dataset)[4] <- "count"

# Select the maximum cell count for the FVG and frequency status
max_counts <- aggregate(count ~ FVG_name + frequency_status, 
                        data = counts_by_dataset, 
                        FUN = max)
names(max_counts)[3] <- "max_count"

# Calculate upper limits
max_counts$upper_limit <- ceiling(max_counts$max_count / 10) * 10

max_counts$frequency_status <- factor(max_counts$frequency_status, 
                                      levels = c("Lower", "Within", "Higher"))

# Add numeric y position for plotting
max_counts$y_pos <- as.numeric(max_counts$frequency_status) + 0.8

# Plot
p1 <- ggplot(data = dat, 
             aes(x = Fire_frequency, 
                 y = factor(frequency_status, levels = c("Lower", 
                                                         "Within", 
                                                         "Higher")), 
                 fill = Dataset,
                 colour = Dataset)) +
  geom_density_ridges(alpha = 0.4, scale = 0.95, stat = 'binline', binwidth = 1, boundary = 0, draw_baseline = F, linewidth = 0.5, panel_scaling = T) + # Plat as a histogram with 21 bins (1 bin per fire frequency), set scale as <1 to ensure histograms do not overlap, do not draw baseline to remove trailing lines 
  geom_text(data = max_counts, aes(x = -1, y = y_pos, label = upper_limit), inherit.aes = FALSE, hjust = 1.4, size = 5,  color = "black") +
  facet_wrap(vars(FVG_name), scales = "free", labeller = as_labeller(Agg_labs), ncol = 2) + # Allow x axis to vary between communities based on the highest fire frequency for that community, allowing easier interrogation
  scale_fill_manual(values = c("gray70", "steelblue", "#492050", "#AAA970"))+
  scale_colour_manual(values = c("gray70", "steelblue", "#492050", "#AAA970"))+
  scale_y_discrete(drop = FALSE, expand = expansion(add = c(0.2, 0.8))) +  # This ensures all factor levels are shown
  scale_x_continuous(breaks = c(0,5,10,15,20), limits = c(-1,21), expand = c(0.01,0)) +  
  labs(x = "Fire Frequency", y = "Fire regime status (number of cells)")+
  coord_cartesian(clip = "off") +  # Allow drawing outside plot area
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
        axis.text = element_text(size = 16),
        legend.text = element_text(size = 20),
        legend.position = c(0.6, 0.2))
p1
#1500x851

# Subplot for fire sensitive vegetation
filt_FVGs_nrecs_df <- FVGs_nrecs_df[!is.na(FVGs_nrecs_df$Fire_frequency) & FVGs_nrecs_df$Fire_frequency != 0, ]
table(filt_FVGs_nrecs_df$FVG_name, round(filt_FVGs_nrecs_df$Fire_frequency))
filt_FVGs_nrecs_df$frequency_status <- factor(filt_FVGs_nrecs_df$frequency_status,
                                              levels = c("No fire", 
                                                         "Do not intentionally burn", 
                                                         "As required"))

# Calculate maximum cell count for each FVG based on the frequency status and dataset 
counts_by_dataset2 <- aggregate(Fire_frequency ~ FVG_name + frequency_status + Dataset, 
                                data = filt_FVGs_nrecs_df, 
                                FUN = length)
names(counts_by_dataset2)[4] <- "count"

# Select the maximum cell count for the FVG and frequency status
max_counts2 <- aggregate(count ~ FVG_name + frequency_status, 
                         data = counts_by_dataset2, 
                         FUN = max)
names(max_counts2)[3] <- "max_count"

# Calculate upper limits
max_counts2$upper_limit <- ceiling(max_counts2$max_count / 10) * 10

max_counts2$frequency_status <- factor(max_counts2$frequency_status, 
                                       levels = c("No fire", 
                                                  "Do not intentionally burn", 
                                                  "As required"))

# Add numeric y position for plotting
max_counts2$y_pos <- ave(rep(1, nrow(max_counts2)), 
                         max_counts2$FVG_name, 
                         FUN = seq_along) + 0.92



p2 <- ggplot(data = filt_FVGs_nrecs_df, 
             aes(x = Fire_frequency, 
                 y = frequency_status,
                 fill = Dataset, colour = Dataset)) +
  geom_density_ridges(alpha = 0.4, scale = 0.95, stat = 'binline', binwidth = 1, boundary = 1, draw_baseline = F, linewidth = 0.5, panel_scaling = T) + # Plat as a histogram with 21 bins (1 bin per fire frequency), set scale as <1 to ensure histograms do not overlap, do not draw baseline to remove trailing lines 
  geom_text(data = max_counts2, aes(x = 0, y = y_pos, label = upper_limit), inherit.aes = FALSE, hjust = 1.7, size = 4,  color = "black") +
  facet_wrap(vars(FVG_name), scales = "free_y", labeller = as_labeller(Agg_labs), drop = T, nrow = 1) + # Allow x axis to vary between communities based on the highest fire frequency for that community, allowing easier interrogation
  scale_fill_manual(values = c("gray80", "steelblue", "#492050", "#AAA970"), labels = c("Public land", "Satellite", "GLM", "GAM"))+
  scale_colour_manual(values = c("gray80", "steelblue", "#492050", "#AAA970"), labels = c("Public land", "Satellite", "GLM", "GAM"))+
  scale_y_discrete(expand = c(0.010), drop = T, labels = c("No fire" = "No fire", "Do not intentionally burn" = "Do not\nintentionally\nburn", "As required" = "As required")) +  
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20),limits = c(0,21), expand = c(0.01,0)) +  
  labs(x = "Fire Frequency", y = "Fire regime status\n (number of cells)")+
  coord_cartesian(clip = "off") +  # Allow drawing outside plot area
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
  geom_bar(stat = "count", position = position_dodge(width = 0.26), width = 0.2)+
  scale_fill_manual(values = c("Public land" = "gray80", "Satellite" = "steelblue", 
                               "GLM" = "#492050", "GAM" = "#AAA970"))+
  scale_colour_manual(values = c("Public land" = "gray80", "Satellite" = "steelblue", 
                                 "GLM" = "#492050", "GAM" = "#AAA970"))+
  scale_x_continuous(breaks = 1, labels = "0") +
  scale_y_continuous(limits = c(0,1000))+
  labs(y = "Number\n of cells") +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(colour = 'black'),
        strip.text = element_text(face = "bold", size = 16),
        legend.position = "none",
        axis.title.y = element_text(face = 'bold', size = 16),
        axis.title.x = element_blank(),
        axis.text.y = element_text(vjust = 0),
        axis.text = element_text(size = 12),
        plot.background = element_blank(),
        axis.ticks.x = element_blank())
riparian

wet <- ggplot(data = FVGs_nrecs_0s[FVGs_nrecs_0s$FVG == "Wet tall open forests", ],
              aes(x = 1, fill = Dataset, colour = Dataset)) +  # Changed x to constant 1
  geom_bar(stat = "count", position = position_dodge(width = 0.26), width = 0.2)+
  scale_fill_manual(values = c("Public land" = "gray80", "Satellite" = "steelblue", 
                               "GLM" = "#492050", "GAM" = "#AAA970"))+
  scale_colour_manual(values = c("Public land" = "gray80", "Satellite" = "steelblue", 
                                 "GLM" = "#492050", "GAM" = "#AAA970"))+
  scale_y_continuous(limits = c(0,1000))+
  scale_x_continuous(breaks = 1, labels = "0") +
  labs(y = "Number\n of cells") +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(colour = 'black'),
        strip.text = element_text(face = "bold", size = 16),
        legend.position = "none",
        axis.title.y = element_text(face = 'bold', size = 16),
        axis.title.x = element_blank(),
        axis.text.y = element_text(vjust = 0),
        axis.text = element_text(size = 12),
        plot.background = element_blank(),
        axis.ticks.x = element_blank())
wet

rain <- ggplot(data = FVGs_nrecs_0s[FVGs_nrecs_0s$FVG == "Rainforests, dry vine forests and brigalow communities", ],
               aes(x = 1, fill = Dataset, colour = Dataset)) +  # Changed x to constant 1
  geom_bar(stat = "count", position = position_dodge(width = 0.26), width = 0.2)+
  scale_fill_manual(values = c("Public land" = "gray80", "Satellite" = "steelblue", 
                               "GLM" = "#492050", "GAM" = "#AAA970"))+
  scale_colour_manual(values = c("Public land" = "gray80", "Satellite" = "steelblue", 
                                 "GLM" = "#492050", "GAM" = "#AAA970"))+
  scale_y_continuous(limits = c(0,1000))+
  scale_x_continuous(breaks = 1, labels = "0") +
  labs(y = "Number\n of cells") +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(colour = 'black'),
        strip.text = element_text(face = "bold", size = 16),
        legend.position = "none",
        axis.title.y = element_text(face = 'bold', size = 16),
        axis.title.x = element_blank(),
        axis.text.y = element_text(vjust = 0),
        axis.text = element_text(size = 12),
        plot.background = element_blank(),
        axis.ticks.x = element_blank())
rain

mang <- ggplot(data = FVGs_nrecs_0s[FVGs_nrecs_0s$FVG == "Mangroves and saltmarsh", ],
               aes(x = 1, fill = Dataset, colour = Dataset)) +  # Changed x to constant 1
  geom_bar(stat = "count", position = position_dodge(width = 0.26), width = 0.2)+
  scale_fill_manual(values = c("Public land" = "gray80", "Satellite" = "steelblue", 
                               "GLM" = "#492050", "GAM" = "#AAA970"))+
  scale_colour_manual(values = c("Public land" = "gray80", "Satellite" = "steelblue", 
                                 "GLM" = "#492050", "GAM" = "#AAA970"))+
  labs(y = "Number\n of cells") +
  scale_x_continuous(breaks = 1, labels = "0") +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(colour = 'black'),
        strip.text = element_text(face = "bold", size = 16),
        legend.position = "none",
        axis.title.y = element_text(face = 'bold', size = 16),
        axis.title.x = element_blank(),
        axis.text.y = element_text(vjust = 0),
        axis.text = element_text(size = 12),
        plot.background = element_blank(),
        axis.ticks.x = element_blank())
mang



p2_with_insets <- ggdraw(p2) +
  # Mangroves
  draw_plot(mang, x = 0.18, y = 0.48, width = 0.12, height = 0.31) +
  # Rainforest
  draw_plot(rain, x = 0.40, y = 0.48, width = 0.12, height = 0.31) +
  # Riparian
  draw_plot(riparian, x = 0.63, y = 0.48, width = 0.12, height = 0.31) +
  # Wet tall open forest
  draw_plot(wet, x = 0.88, y = 0.48, width = 0.12, height = 0.31)
p2_with_insets


combined_plot <- plot_grid(p1, p2_with_insets, ncol = 1, rel_heights = c(4, 1))
#combined_plot

ggsave("./04_Results/Plots/Vegetation_validation/Veg_fire_regime.tiff", combined_plot, width = 16, height = 18, dpi = 300)


table(round(FVGs_withrecs_df$Fire_frequency), FVGs_withrecs_df$Dataset)





