# Written by Felicity Charles
# Date: 26/02/2025

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
SEQ_bound <- vect('./00_Data/SEQ_bound/SEQ.gpkg')
Aus <- download.file("https://www.abs.gov.au/statistics/standards/australian-statistical-geography-standard-asgs-edition-3/jul2021-jun2026/access-and-downloads/digital-boundary-files/STE_2021_AUST_SHP_GDA2020.zip", destfile = './00_Data/Spatial data/Australia.zip', mode = "wb", cacheOK = F)
unzip(zipfile = './00_Data/Spatial data/Australia.zip', exdir = './00_Data/Spatial data/Australia')
SEQ <- vect('./00_Data/Spatial data/Australia/STE_2021_AUST_GDA2020.shp') %>%
  project("EPSG:3577") %>% 
  crop(SEQ_bound)
BVG <- vect('./00_Data/Environmental_data/Remnant_2021_broad_veg_groups/Remnant_broad_vegetation_groups.shp') %>%
  project('EPSG:3577') %>% 
  crop(SEQ)

protected_land <- vect('./00_Data/Protected_areas/Protected_areas.shp') %>% 
  project('EPSG:3577') %>% 
  crop(SEQ)


places <- vect('./00_Data/Environmental_data/Place_names/QSC_Extracted_Data_20250226_155128183750-19836/Place_names_gazetteer.shp') %>% 
  project('EPSG:3577') %>% 
  crop(SEQ)

places <- places[places$place_name == "Warwick" | places$place_name =="Brisbane"| places$place_name =="Toowoomba"| places$place_name =="Kingaroy"| places$place_name =="Gympie"|  places$place_name =="Caboolture"| places$place_name =="Redcliffe"| places$place_name =="Gold Coast"| places$place_name =="Sunshine Coast"]

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
  geom_sf_text(data = places, aes(label = place_name, geometry = geometry), show.legend = F, fontface = "bold")+
  labs(x = "", y = "") +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 10))
  
  




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
            


