library(tidyverse)
library(here)
library(sf)
library(terra)
library(cowplot)
library(ggsflabel)
# library(ggspatial)

#Shoreline shapefile
wa_outline <- here("data", "spatial", "WA_state_outline", "WA_state_outline.shp") %>% 
  read_sf()

shoreline <- here("data","spatial", "shorezone_shoreline_only", "shorezone_shoreline_only.shp") %>% 
  read_sf(crs = 2927) %>%  #Washington State Plane South (ft) / NAD83
  st_transform(crs = 4326)

# shoreline_crop <- st_crop(shoreline, ext(c(-124, -121, 47, 49))) #crop to the region of interest
# # plot(shoreline_crop)
# 
# shoreline_crop_wide <- st_crop(shoreline, ext(c(-124.2, -120, 47, 49))) #crop to the region of interest
# plot(shoreline_crop_wide)

#GPS locations for our survey stations with each ipa
SOS_core_sites <- factor(c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"), 
                         levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"))

SOS_sites <- here("data", "spatial","reupdatingthearmoringgeodatabase", "Shoreline_armoring_shore_origin_sites_UTM.shp") %>% 
  read_sf() %>% #UTM zone 10 32610
  st_transform(crs = 4326) %>% #transform site coordinates into the same crs as the shoreline layer
  filter(!Site_name %in% c("Lost Lake", "Titlow Park", "Penrose Point", "Waterman Shoreline Preserve", "Howarth Park", "Maylor Point")) %>% 
  mutate(site = case_when(Site_name == "Dockton Park" ~ "DOK",
                          Site_name == "Cornet Bay" ~ "COR",
                          Site_name == "Edgewater Beach" ~ "EDG",
                          Site_name == "Family Tides" ~ "FAM",
                          Site_name == "Seahurst Park" ~ "SHR", 
                          Site_name == "Turn Island" ~ "TUR"), .after = "Site_name") %>% 
  mutate(IPA = ifelse(site == "TUR" & IPA == "Restored", "Natural", IPA)) %>% 
  mutate(IPA = factor(IPA), site = factor(site, levels = SOS_core_sites), coords = st_coordinates(.)) %>% 
  arrange(site, IPA) %>% 
  mutate(site_ID = row_number(), .after = "site")


#load 30m resolution land cover data
ccap_2016lc <- rast(here("data", "spatial", "2016_CCAP_Job987470", "2016_CCAP_J987470.tif"))

#classify land cover data based on NOAA codes
cover <- c("background", "unclassified", "high_intensity_developed",
           "medium_intensity_developed", "low_intensity_developed",
           "developed_open_space", "cultivated_land", "pasture.hay", 
           "grassland", "deciduous_forest", "evergreen_forest",
           "mixed_forest", "shrub.scrub", "palustrine_forested_wetland",
           "palustrine_scrub.shrub_wetland", "palustrine_emergent_wetland",
           "estuarine_forested_wetland", "estuarine_scrub.shrub_wetland",
           "estuarine_emergent_wetland", "unconsolidated_shore", "bare_land",
           "open_water", "palustrine_aquatic_bed", "estuarine_aquatic_bed","tundra", "snow.ice")

levels(ccap_2016lc) <- data.frame(id=0:25, cover=cover)
# plot(ccap_2016lc)

simple_cover <- matrix(
  c(0,0, 1,1, 2,2, 3,2, 4,2, 5,2, 6,3, 7,3, 8,4, 9,5, 10,5, 11,5, 12,6, 
    13,7, 14,7, 15,7, 16,7, 17,7, 18,7, 19,8, 20,10, 21,9, 22,8, 23,8, 24,10, 25,10),
  ncol = 2, 
  byrow = TRUE
)

ccap_2016lc_simple <- terra::classify(ccap_2016lc, rcl = simple_cover)

#classify land cover data based on NOAA groupings of the codes
cover_simple <- c("background", "unclassified", "developed",
           "agricultural", "grassland",
           "forest", "scrub", "wetland",
           "shoreline","water", "barren")

levels(ccap_2016lc_simple) <- data.frame(id=0:10, cover=cover_simple)


# plot(ccap_2016lc_simple)
simple_ccap <- project(ccap_2016lc_simple, "EPSG:4326")
# plot(simple_ccap)

ccap_clip <- crop(simple_ccap, ext(c(-123.42, -122, 47, 49))) #crop to the region of interest
# plot(ccap_clip)

ccap_lower_res <- terra::aggregate(ccap_clip, fact = 10, fun = mean)
# plot(ccap_lower_res)

ccap_df <- as.data.frame(ccap_lower_res, xy = TRUE) %>%  # Keeps x, y coordinates
  mutate(cover_group = case_when(cover < 1 ~ "background",
                             cover >= 1 & cover < 2 ~ "unclassified",
                             cover >= 2 & cover < 3 ~ "developed",
                             cover >= 3 & cover < 4 ~ "agricultural",
                             cover >= 4 & cover < 5 ~ "grassland",
                             cover >= 5 & cover < 6 ~ "forest",
                             cover >= 6 & cover < 7 ~ "scrub",
                             cover >= 7 & cover < 8 ~ "wetland",
                             cover >= 8 & cover < 9 ~ "shoreline",
                             cover >= 9 & cover < 10 ~ "water",
                             cover >= 10 ~ "barren")) 
  
ccap_full_df <- as.data.frame(ccap_clip, xy = TRUE)
ggplot() +
  geom_raster(data = ccap_full_df, aes(x = x, y = y, fill = as.factor(cover))) +
  geom_sf(data = wa_outline, fill = "transparent") +
  scale_fill_manual(values = c("lightgrey", "yellow", "darkslategrey", "darkgreen",  "lightgreen", "lightyellow","tan","white", "lightblue")) +
  coord_sf(xlim = c(-123.42, -122), ylim = c(47, 49)) 
  
land_cover_cols <- c("lightgrey","brown","lightgreen","darkgreen","lightyellow","lightblue","tan","white","darkslategrey")

region_lc <- ggplot() +
  geom_sf(data = SOS_sites, aes(shape = IPA), color = "magenta", size = 3) + #size = 7
  geom_raster(data = ccap_full_df, aes(x = x, y = y, fill = as.factor(cover))) + 
  # geom_sf(data = wa_outline, fill = "transparent", linewidth = 0.5) +
  scale_fill_manual(values = land_cover_cols) +
  # geom_sf(data = wa_outline, fill = "lightgrey", linewidth = 0.5) +
  coord_sf(xlim = c(-123.42, -122), ylim = c(47, 49), clip = 'off') +
  # coord_sf(xlim = c(-124, -121), ylim = c(47, 49)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.box = "horizontal") +
  guides(shape = guide_legend(order = 1),
         fill = guide_legend(order = 2)) +
  labs(shape = "Category", fill = "Land cover")

lc_legend <- get_legend(region_lc)

#whole region
region <- ggplot() +
  # geom_raster(data = ccap_full_df, aes(x = x, y = y, fill = as.factor(cover))) + 
  # geom_sf(data = wa_outline, fill = "transparent", linewidth = 0.5) +
  # scale_fill_manual(values = land_cover_cols) +
  geom_sf(data = wa_outline, fill = "lightgrey", linewidth = 0.5) +
  geom_sf(data = SOS_sites, aes(shape = IPA), color = "magenta", size = 3) + #size = 7
  coord_sf(xlim = c(-123.42, -122), ylim = c(47, 49), clip = 'off') +
  # coord_sf(xlim = c(-124, -121), ylim = c(47, 49)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none") +
  labs(fill = "Land cover")

#san juans
# sji <- ggplot() +
#   geom_raster(data = ccap_full_df, aes(x = x, y = y, fill = as.factor(cover))) + 
#   geom_sf(data = wa_outline, fill = "transparent", linewidth = 0.5) +
#   scale_fill_manual(values = land_cover_cols) +
#   geom_sf_label_repel(data = filter(SOS_sites, site %in% c("FAM", "TUR")), 
#                       aes(label = site_ID), color = "black", size = 5,  
#                nudge_y = ifelse(SOS_sites$site == "FAM", 0, .005),
#                nudge_x = ifelse(SOS_sites$site == "FAM", 0.005, 0)) +
#   geom_sf(data = SOS_sites, aes(shape = IPA), color = "magenta", size = 4) +
#   coord_sf(xlim = c(-123.03, -122.91), ylim = c(48.52, 48.62)) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text.x  = element_blank(),
#         axis.text.y  = element_blank(),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(), 
#         legend.position = "none") 

fam <- ggplot() +
  geom_raster(data = ccap_full_df, aes(x = x, y = y, fill = as.factor(cover))) + 
  geom_sf(data = wa_outline, fill = "transparent", linewidth = 0.5) +
  scale_fill_manual(values = land_cover_cols) +
  geom_sf_label_repel(data = filter(SOS_sites, site %in% c("FAM")), 
                      aes(label = site_ID), color = "black", size = 5,  
                      nudge_x =  0.005) +
  geom_sf(data = SOS_sites, aes(shape = IPA), color = "magenta", size = 4) +
  coord_sf(xlim = c(-123.02, -122.94), ylim = c(48.575, 48.64)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.ticks = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        legend.position = "none",
        plot.margin=grid::unit(c(0,0,0,0), "mm")) 

tur <- ggplot() +
  geom_raster(data = ccap_full_df, aes(x = x, y = y, fill = as.factor(cover))) + 
  geom_sf(data = wa_outline, fill = "transparent", linewidth = 0.5) +
  scale_fill_manual(values = land_cover_cols) +
  geom_sf_label_repel(data = filter(SOS_sites, site %in% c("TUR")), 
                      aes(label = site_ID), color = "black", size = 5,  
                      nudge_x =  0.005) +
  geom_sf(data = SOS_sites, aes(shape = IPA), color = "magenta", size = 4) +
  coord_sf(xlim = c(-123.03, -122.91), ylim = c(48.50, 48.56)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.ticks = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(), 
        legend.position = "none",
        plot.margin=grid::unit(c(0,0,0,0), "mm")) 


#cornet
cor <- ggplot() +
  geom_raster(data = ccap_full_df, aes(x = x, y = y, fill = as.factor(cover))) + 
  geom_sf(data = wa_outline, fill = "transparent", linewidth = 0.5) +
  scale_fill_manual(values = land_cover_cols) +
  geom_sf_label_repel(data = filter(SOS_sites, site %in% c("COR")), 
                      aes(label = site_ID), color = "black", size = 5,  
                      nudge_x =  0.005) +
  geom_sf(data = SOS_sites, aes(shape = IPA), color = "magenta", size = 4) +
  coord_sf(xlim = c(-122.675, -122.585), ylim = c(48.38, 48.425)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.ticks = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.margin=grid::unit(c(0,0,0,0), "mm"))

#seahurst and dockton
# ggplot() +
#   geom_sf(data = wa_outline) +
#   geom_sf(data = SOS_sites, pch = 18, size = 3) + #size = 7
#   coord_sf(xlim = c(-122.5, -122.3), ylim = c(47.35, 47.5)) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.border = element_blank(),
#         axis.ticks = element_blank(),
#         axis.text.x  = element_blank(),
#         axis.text.y  = element_blank(),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank())

#seahurst
shr <- ggplot() +
  geom_raster(data = ccap_full_df, aes(x = x, y = y, fill = as.factor(cover))) + 
  geom_sf(data = wa_outline, fill = "transparent", linewidth = 0.5) +
  scale_fill_manual(values = land_cover_cols) +
  geom_sf_label_repel(data = filter(SOS_sites, site %in% c("SHR")), 
                      aes(label = site_ID), color = "black", size = 5,  
                      nudge_x =  -0.005) +
  geom_sf(data = SOS_sites, aes(shape = IPA), color = "magenta", size = 5) +
  coord_sf(xlim = c(-122.38, -122.35), ylim = c(47.47, 47.49)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        axis.ticks = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        plot.margin=grid::unit(c(0,0,0,0), "mm"))

#dockton
dok <- ggplot() +
  geom_raster(data = ccap_full_df, aes(x = x, y = y, fill = as.factor(cover))) + 
  geom_sf(data = wa_outline, fill = "transparent", linewidth = 0.5) +
  scale_fill_manual(values = land_cover_cols) +
  geom_sf(data = SOS_sites, aes(shape = IPA), color = "magenta", size = 3) +
  coord_sf(xlim = c(-122.5, -122.4), ylim = c(47.35, 47.4)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

#edgewater
edg <- ggplot() +
  geom_raster(data = ccap_full_df, aes(x = x, y = y, fill = as.factor(cover))) + 
  geom_sf(data = wa_outline, fill = "transparent", linewidth = 0.5) +
  scale_fill_manual(values = land_cover_cols) +
  geom_sf(data = SOS_sites, aes(shape = IPA), color = "magenta", size = 3) +
  coord_sf(xlim = c(-123, -122.85), ylim = c(47.12, 47.2)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none")

# ggdraw() + 
#   draw_plot(region) + 
#   draw_plot(sji, x = 0.05, y = 0.49,
#             width = 0.4, height = 0.5) +
#   draw_plot(cor, x = 0.64, y = 0.69,
#             width = 0.3, height = 0.3) +
#   draw_plot(shr, x = 0.64, y = 0.35, 
#             width = 0.3, height = 0.3) +
#   draw_plot(dok, x = 0.64, y = 0.01,
#             width = 0.3, height = 0.3) + 
#   draw_plot(edg, x = 0.05, y = 0.01,
#             width = 0.3, height = 0.3)

# library(patchwork)
# 
# sji + edg + region + cor + shr + dok +
#   plot_layout(ncol = 3) +
#   plot_annotation(title = "Land cover and survey sites in Washington State",
#                   theme = theme(plot.title = element_text(hjust = 0.5, size = 16))) +
#   plot_layout(guides = "collect") & theme(legend.position = "bottom")


# region_bbox <- st_bbox(shoreline_crop_wide)

# sji_shoreline_crop <- st_crop(shoreline_crop, ext(c(-123.03, -122.91), c(48.52, 48.62))) #crop to the san juans
# sji_bbox <- st_bbox(sji_shoreline_crop)
# sji_w <- sji_bbox$xmax - sji_bbox$xmin
# sji_h <- sji_bbox$ymax - sji_bbox$ymin

fam_shoreline_crop <- st_crop(shoreline_crop, ext(c(-123.02, -122.94), c(48.575, 48.64))) #crop to fam
fam_bbox <- st_bbox(fam_shoreline_crop)
fam_w <- fam_bbox$xmax - fam_bbox$xmin
fam_h <- fam_bbox$ymax - fam_bbox$ymin

tur_shoreline_crop <- st_crop(shoreline_crop, ext(c(-123.03, -122.91), c(48.50, 48.56))) #crop to tur
tur_bbox <- st_bbox(tur_shoreline_crop)
tur_w <- tur_bbox$xmax - tur_bbox$xmin
tur_h <- tur_bbox$ymax - tur_bbox$ymin

cor_shoreline_crop <- st_crop(shoreline_crop, ext(c(-122.675, -122.585), c(48.38, 48.425))) #crop to cornet bay
cor_bbox <- st_bbox(cor_shoreline_crop)
cor_w <- cor_bbox$xmax - cor_bbox$xmin
cor_h <- cor_bbox$ymax - cor_bbox$ymin

shr_shoreline_crop <- st_crop(shoreline_crop, ext(c(-122.38, -122.35), c(47.47, 47.49))) #crop to seahurst
shr_bbox <- st_bbox(shr_shoreline_crop)
shr_w <- shr_bbox$xmax - shr_bbox$xmin
shr_h <- shr_bbox$ymax - shr_bbox$ymin
  
regionymax <- 50
regionymin <- 47
regionxmax <- -119.5
regionxmin <- -124.2

region + 
  annotation_custom(grob = ggplotGrob(fam),
                    xmax = regionxmin + fam_w*11,
                    xmin = regionxmin,
                    ymin = regionymax - fam_h*11,
                    ymax = regionymax) +
  annotation_custom(grob = ggplotGrob(tur),
                    xmax = regionxmin + tur_w*11,
                    xmin = regionxmin,
                    ymin = 46.7) +
  annotation_custom(grob = ggplotGrob(cor),
                    xmax = regionxmax,
                    ymin = regionymax - cor_h*11,
                    ymax = regionymax) +
  annotation_custom(grob = ggplotGrob(shr),
                    xmax = regionxmax,
                    xmin = regionxmax - shr_w*125,
                    ymin = 46.7)
  # annotation_custom(grob = lc_legend, 
  #                   xmax = regionxmax,
  #                   ymin = 47.1)

## make a regional layer that includes canada and us/ca border?  
