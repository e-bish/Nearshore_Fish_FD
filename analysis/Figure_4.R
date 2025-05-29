#load libraries
library(tidyverse)
library(here)

load(here("data","mFD_results.Rdata"))
load(here("data", "fish_L_sample.Rdata")) 

#### Plot the assemblages weighted by abundance (do sites have more generalists - points in the middle of the trait space?)
spp_weights <- fish_L_sample %>% 
  rownames_to_column(var = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa", "year")) %>% 
  pivot_longer( !c(site, ipa, year), names_to = "species") %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")),
         region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"), 
         veg = ifelse(site %in% c("TUR", "COR", "SHR"), "present", "absent"), .after = site)

spp_weights_site <- spp_weights %>% 
  group_by(site, species) %>% 
  summarize(avg_abund = mean(value))

fish_L_pa <- decostand(fish_L_sample, method = "pa")

fish_pa_df <- fish_L_pa %>% 
  rownames_to_column(var = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa", "year")) %>% 
  pivot_longer( !c(site, ipa, year), names_to = "species") %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")),
         region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"), 
         veg = ifelse(site %in% c("TUR", "COR", "SHR"), "present", "absent"), .after = site) %>% 
  filter(!value == 0) %>% 
  dplyr::select(!value)

spp_coords <- trait_space %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "species")

spp_coords_df <- fish_pa_df %>% 
  left_join(spp_coords)

site_abund_df <- spp_coords_df %>% 
  dplyr::select(!10:16) %>% 
  left_join(spp_weights_site)

FRic_df <- fish_pa_df %>% 
  left_join(spp_coords)

calculate_hull <- function(df) {
  df[chull(df$PC1, df$PC2), ]
}

FRic.hulls_site <- FRic_df %>% 
  group_by(site) %>% 
  group_split() %>%
  map_dfr(~ calculate_hull(.), .id = "hull_no") 

veg_colors <- c("lightgrey", "lightgreen", "lightgreen", "lightgreen", "lightgrey", "lightgrey")

site_abund_df %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_polygon(data = FRic.hulls_site, aes(x = PC1, y = PC2, fill = hull_no, group = site), 
               show.legend = FALSE, alpha = 0.75, inherit.aes = FALSE) +
  scale_fill_manual(values = veg_colors) +
  geom_point(aes(size = avg_abund)) +
  theme_classic() +
  theme(strip.text = element_text(face = "bold"), 
        strip.background = element_rect(fill = NA, colour = NA),
        legend.title = element_text(hjust = 0.5)) + 
  labs(size = "Average\nAbundance", x = "PCoA axis 1", y = "PCoA axis 2") + 
  facet_wrap(~site)

#### region ####
spp_weights_region <- spp_weights %>% 
  group_by(region, species) %>% 
  summarize(avg_abund = mean(value)) %>% 
  filter(avg_abund > 0)

fish_pa_df_region <- fish_L_pa %>% 
  rownames_to_column(var = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa", "year")) %>% 
  pivot_longer( !c(site, ipa, year), names_to = "species") %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")),
         region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"), 
         veg = ifelse(site %in% c("TUR", "COR", "SHR"), "present", "absent"), .after = site) %>% 
  filter(!value == 0) %>% 
  dplyr::select(region, species) %>% 
  distinct()

spp_coords_df.r <- fish_pa_df_region %>% 
  left_join(spp_coords)

region_abund_df <- spp_coords_df.r %>% 
  dplyr::select(!5:12) %>% 
  left_join(spp_weights_region)

FRic_df.r <- fish_pa_df_region %>% 
  left_join(spp_coords)

FRic.hulls_region <- FRic_df.r %>% 
  group_by(region) %>% 
  group_split() %>%
  map_dfr(~ calculate_hull(.), .id = "hull_no") 

region_abund_df %>% 
  ggplot(aes(x = PC1, y = PC2)) +
  geom_polygon(data = FRic.hulls_region, aes(x = PC1, y = PC2, fill = hull_no, group = region), 
               show.legend = FALSE, alpha = 0.75, inherit.aes = FALSE) +
  # scale_fill_manual(values = veg_colors) +
  geom_point(aes(size = avg_abund)) +
  theme_classic() +
  theme(strip.text = element_text(face = "bold"), 
        strip.background = element_rect(fill = NA, colour = NA),
        legend.title = element_text(hjust = 0.5)) + 
  labs(size = "Average\nAbundance", x = "PCoA axis 1", y = "PCoA axis 2") + 
  facet_wrap(~region)


mFD::alpha.multidim.plot(output_alpha_fd_multidim = alpha_fd_indices,
                         )
