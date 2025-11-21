# Written by Felicity Charles
# Date: 26/02/2025



#####UPDATE WITH IBRA CROPPED IF REQUIRED
##### Fire frequency analysis ----
# General study area map showing the extent of remnant vegetation cover, protected areas, and landmarks in SEQ with inset map of Queensland


# 1. Load packages -----
library(terra)
library(dplyr)
library(ggplot2)
library(tidyterra)
library(ggspatial)
library(cowplot)
library(sf)


# 2. Read and format spatial data ----
SEQ_bound <- vect('./00_Data/SEQ_bound/SEQ_IBRA.gpkg')
Aus <- download.file("https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs-edition-3/jul2021-jun2026/access-and-downloads/digital-boundary-files/STE_2021_AUST_SHP_GDA2020.zip", destfile = './00_Data/Spatial data/Australia.zip', mode = "wb", cacheOK = F)
unzip(zipfile = './00_Data/Spatial data/Australia.zip', exdir = './00_Data/Spatial data/Australia')
SEQ <- vect('./00_Data/Spatial data/Australia/STE_2021_AUST_GDA2020.shp') %>%
  project("EPSG:3577") %>% 
  crop(SEQ_bound)

Sentinel_ff <- rast('./00_Data/Fire_data/Outputs/Sentinel/Sentinel_ff_hydrographical_mask_SEQ_IBRA_focal_cropped.tif') %>% 
  mask(SEQ)
QPWS_ff <- rast('./00_Data/Fire_data/Outputs/SEQ/QPWS_SEQ_IBRA_freq_hydrographical_mask_cropped.tif')
# Combine fire frequency datasets into one, keeping the value for QPWS first where it is not 0

QPWS_ff <- ifel(QPWS_ff == 0, NA, QPWS_ff) # Extract the cells where QPWS_ff is not zero
Sent_ff <- mask(Sentinel_ff, QPWS_ff, inverse = T) # Mask areas of Sentinel data where there is QPWS data
SEQ_ff <- merge(QPWS_ff, Sent_ff)
par(mfrow = c(1,3)); plot(QPWS_ff); plot(SEQ_ff); plot(Sent_ff)


protected_land <- vect('./00_Data/Protected_areas/Protected_areas.shp') %>% 
  project('EPSG:3577') %>% 
  crop(SEQ)
protected_land$Public_estate <- "Public estate"


# 3. Produce main study area map
p1 <-
  ggplot() + 
  geom_spatraster(data = SEQ_ff) +
  scale_fill_viridis_c(option = 'viridis', na.value = 'transparent', breaks = c(1,5,10,15,20,25, 30), limits = c(1,30)) +
  labs(fill = "Fire frequency")+
  theme_bw() +
  theme_cowplot(font_size = 15) +
  geom_spatvector(data = protected_land, aes(colour = 'Public estate'), fill = 'gray', alpha = 0.5, show.legend = T) +
  scale_colour_manual(values = 'gray', name = "") +
  annotation_north_arrow(location = "bl", which_north = T, height = unit(1, "cm"), width = unit(0.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(9, "cm"), style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl", text_cex = 1)+
  labs(x = "", y = "") +
  theme(legend.key.height = unit(1, 'cm'),
        legend.text = element_text(size = 15),
        legend.title = element_text(face = 2, size = 18),
        plot.background = element_rect(fill = "white"))


p1


Aus <- vect('./00_Data/Spatial data/Australia/STE_2021_AUST_GDA2020.shp') %>%
  project("EPSG:3577")


inset <-
  ggplot()+
  geom_spatvector(data = Aus, fill = NA, col = 'midnightblue')+
  geom_spatvector(data = SEQ_bound, col = 'black', fill = NA, lwd = 1) +
  geom_spatraster(data = SEQ_ff) +
  scale_fill_viridis_c(option = 'viridis', na.value = 'transparent', breaks = c(0,5,10,15,20,25, 30), limits = c(0,30))+
  theme_void() +
  theme(legend.position = "none") +
  annotation_scale(location = "bl", text_cex = 1, pad_x = unit(2, "cm"))



combined <- ggdraw()+
  draw_plot(p1, x = 0, y = 0, width = 1, height = 0.85) +  # p1 takes bottom 85%
  draw_plot(inset, x = 0.22, y = 0.65, width = 1.2, height = 0.16)  
ggsave("./04_Results/Study_area_fire_map.png", combined, width = 650, height = 1176, units = "px", dpi = 96)









# OLD Code to be edited
BVG <- vect('./00_Data/Environmental_data/Remnant_2021_broad_veg_groups/Remnant_broad_vegetation_groups.shp') %>%
  project('EPSG:3577') %>% 
  crop(SEQ) 

protected_land <- vect('./00_Data/Protected_areas/Protected_areas.shp') %>% 
  project('EPSG:3577') %>% 
  crop(SEQ)


places <- vect('./00_Data/Environmental_data/Place_names/QSC_Extracted_Data_20250226_155128183750-19836/Place_names_gazetteer.shp') %>% 
  project('EPSG:3577') %>% 
  crop(SEQ)

places <- places[places$place_name == "Warwick" | places$place_name =="Brisbane"| places$place_name =="Toowoomba"| places$place_name =="Gladstone"| places$place_name =="Bundaberg" | places$place_name =="Gold Coast"| places$place_name =="Sunshine Coast"]

places <- st_as_sf(places)

length(places$place_name)
unique(places$place_name)
places$place_name

places <- places[!duplicated(places$place_name),]
length(places)

# Remove areas that are water and assign plantation and non-remnant as NA values
rem_nat_veg <- BVG[BVG$dbvg5m != "water" & BVG$dbvg5m != "plantation" & BVG$dbvg5m != "non-remnant"]
unique(rem_nat_veg$dbvg5m)
rem_nat_veg$Remnant_veg <- "Remnant native vegetation" # Add column for colour

protected_land$Public_estate <- "Public estate"




# We want to calculate the area of SEQ land, Rainforest vegetation, and Sclerophyll vegetation. Then we want to calculate area of rainforest burnt in 2019-2020
Rainforest <- BVG[BVG$dbvg5m == 1, ]
Sclerophyll <- BVG[BVG$dbvg5m == 2 | BVG$dbvg5m == 3 | BVG$dbvg5m == 4 | BVG$dbvg5m == 5 | BVG$dbvg5m == 6 | BVG$dbvg5m == 7 | BVG$dbvg5m == 8 | BVG$dbvg5m == 9 | BVG$dbvg5m == 10 | BVG$dbvg5m == 11, ]


Rainforest$ha <- expanse(Rainforest)/10000
Sclerophyll$ha <- expanse(Sclerophyll)/10000
SEQ$ha <- expanse(SEQ)/10000
rem_nat_veg$ha <- expanse(rem_nat_veg)/10000


QPWS1987 <- vect('./00_Data/Fire_data/Outputs/QPWS_fire_hist_1987.gpkg') %>% 
  crop(SEQ)


Rainforest_fire <- crop(QPWS1987, Rainforest)
Rainforest_wf <- Rainforest_fire[Rainforest_fire$TYPE == "WF"]

Rainforest_wf_ha <- expanse(Rainforest_wf)/10000
Rainforest_fire_ha <- expanse(Rainforest_fire)/10000


# Lets do how much cover is rainforest of remnant veg cover and sclero of remnant veg cover

((sum(Rainforest$ha))/(sum(rem_nat_veg$ha)))*100

(sum(Sclerophyll$ha)/sum(rem_nat_veg$ha))*100

(sum(Rainforest_wf_ha)/sum(Rainforest_fire_ha))*100



# 3. Produce the main study area map ----

library(ggrepel)

places_coords <- st_coordinates(places)
places_df <- data.frame(
  place_name = places$place_name,
  x = places_coords[,1],
  y = places_coords[,2])

p1 <-
  ggplot() + 
  geom_spatvector(data = rem_nat_veg, aes(colour = Remnant_veg), fill = 'gray') +
  theme_bw() +
  theme_cowplot(font_size = 15)+  
  geom_spatvector(data = protected_land, aes(colour = Public_estate), fill = NA, show.legend = T) +
  scale_colour_manual(values = c('black', 'gray')) +
  geom_spatvector(data = SEQ, colour = "midnightblue", fill = NA) +
  annotation_north_arrow(location = "bl", which_north = T, height = unit(1, "cm"), width = unit(0.75, "cm"), pad_y = unit(0.1, "cm"), pad_x = unit(9, "cm"), style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl", text_cex = 1)+
  geom_text_repel(data = places_df, aes(x = x, y = y, label = place_name), 
                  show.legend = F, fontface = "bold",
                  box.padding = 1.5, 
                  point.padding = 1,
                  min.segment.length = 0,
                  segment.color = "gray50",
                  segment.size = 0.3,
                  force = 10,
                  max.overlaps = Inf)+
  labs(x = "", y = "") +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 10))


p1



# Create the inset map
# Read in the data required for the inset ----
Aus <- vect('./00_Data/Spatial data/Australia/STE_2021_AUST_GDA2020.shp') %>%
  project("EPSG:3577")


inset <-
  ggplot()+
  geom_spatvector(data = Aus, fill = NA, col = 'midnightblue')+
  geom_spatvector(data = SEQ_bound, col = 'black', fill = NA, lwd = 1.5) +
  theme_void() +
  theme(legend.position = "none") +
  annotation_scale(location = "bl", text_cex = 1, pad_x = unit(3, "cm"))







ggdraw()+
  draw_plot(p1)+
  draw_plot(inset, height = 0.4, 
            x = 0.15, y = 0.59)



