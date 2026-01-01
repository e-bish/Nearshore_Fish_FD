library(tidyverse)
library(here)
library(vegan)
library(janitor)
library(pairwiseAdonis)
library(ggrepel)
library(patchwork)


set.seed(2025)
load(here("data", "net_core.Rdata")) #created in "02_tidy_data.R"
load(here("data", "fish_L_full.Rdata")) #created in "03_create_matrices.R"
fish_Q <- read_csv(here("data", "fish_traits.csv")) #created in "03_create_matrices.R"
load(here("data", "wq_tb.Rdata")) #created in "Figure_S1_water_quality.R"

site_colors <- rev(c("#8c510a","#d8b365", 
                     "lightgoldenrod",
                     "lightblue",
                     "#5ab4ac", "#01665e"))

#### summarize catch ####
total_catch <- sum(net_core$species_count)

net_core %>% 
  group_by(tax_group) %>% 
  summarize(sum = sum(species_count), prop = sum/total_catch) %>% 
  arrange(-prop)

#eelgrass community
net_core %>% 
  filter(tax_group %in% c("Stichaeid", "Sygnathid", "Aulorhynchus", "Pholidae")) %>% 
  summarize(sum = sum(species_count), prop = sum/total_catch)

#catch by station
catch_by_station <- net_core %>% 
  group_by(station) %>% 
  summarize(catch_sum = sum(species_count)) %>% 
  mutate(prop = catch_sum/sum(total_catch))

# all_catch_plot <- catch_by_station %>% 
#   ggplot(aes(x = factor(station, labels = c("shallow", "middle", "deep")), 
#                         y = catch_sum)) +
#   geom_bar(stat = 'identity') +
#   theme_classic() + 
#   labs(x = "Station", y = "Sum of Catch", title = "All species catch")
  
#demersal catch by station
demersal_spp <- fish_Q %>% 
  filter(DemersPelag %in% "demersal")

demersal_catch_by_station <- net_core %>% 
  filter(ComName %in% demersal_spp$ComName) %>% 
  group_by(station) %>% 
  summarize(catch_sum = sum(species_count)) %>% 
  mutate(prop = catch_sum/sum(catch_sum))

net_core %>% 
  mutate(Demers = ifelse(ComName %in% demersal_spp$ComName, "Demersal", "Not Demersal")) %>% 
  group_by(Demers, station) %>% 
  summarize(catch_sum = sum(species_count)) %>% 
  ggplot(aes(x = factor(station, labels = c("shallow", "middle", "deep")), 
             y = catch_sum,
             fill = factor(Demers, levels = c("Not Demersal", "Demersal")))) +
  geom_bar(stat = 'identity', color = "black") +
  scale_fill_manual(values = c("white", "lightgrey")) +
  theme_classic() +
  labs(x = "Station", y = "Sum of catch", 
       fill = "Position in the\nwater column")

ggsave(here("figures", "Fig_S1.png"), 
       width = 6, height = 5.5, dpi = 300)

#### nmds 1 ####
fish_L_all <- fish_L_full %>% 
  mutate(ipa = ifelse(ipa == "Natural2", "Natural", ipa)) %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")),
         ipa = factor(ipa, levels = c("Natural",
                                      "Armored",
                                      "Restored")))

fish_dist <- vegdist(sqrt(fish_L_all[4:45]), method = "bray")

nmds <- metaMDS(fish_dist, k= 3, trymax=1000, plot = FALSE)

points <- data.frame(nmds$points) %>% 
  cbind(fish_L_all %>% select(site:year)) %>% 
  mutate(yr_abb = str_sub(year, start = 3))

hulls <- points %>%
  group_by(site) %>% 
  slice(chull(MDS1,MDS2))

nmds_all <- ggplot() +
  geom_point(data = points, aes(x = MDS1, y = MDS2, 
                                color = site, shape = ipa), 
             size = 2) + 
  # geom_text_repel(data = points, 
  #                 aes(x = MDS1, y = MDS2, color = site, label = yr_abb), 
  #                 max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) +
  theme_classic(base_size = 10) +
  annotate("text", x = 0.7, y = 1.4, size = 3,
           label = paste("Stress = ", round(nmds$stress, digits = 3))) +
  annotate("text", label = "A", 
           x = -2, y = 1.4,
           size = 6, fontface = "bold") +
  labs(color = "Site", shape = "Condition\ncategory")

#### nmds 2 ####

wq_tb_exp <- wq_tb %>%
  pivot_longer(c(secchi_depth_m, do_mg_l, salinity_ppm, temperature), names_to = "metric") %>%
  group_by(shoreline, year, metric) %>%
  summarize(mean = round(mean(value, na.rm = TRUE), 2), min = min(value, na.rm = TRUE), max = max(value, na.rm = TRUE)) %>%
  ungroup() %>% 
  mutate(site = str_sub(shoreline, 1, 3),
         ipa = str_sub(shoreline, 4), .before = shoreline) %>% 
  select(!c(shoreline,min,max)) %>% 
  pivot_wider(names_from = metric, values_from = mean) %>% 
  mutate(year = factor(year))

fish_L_abb <- fish_L_full %>% 
  filter(!year %in% c(2018, 2019)) %>% 
  full_join(wq_tb_exp) %>% 
  mutate(ipa = ifelse(ipa == "Natural2", "Natural", ipa)) %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")),
         ipa = factor(ipa, levels = c("Natural",
                                      "Armored",
                                      "Restored")))
  
fish_dist_abb <- vegdist(sqrt(fish_L_abb[4:45]), method = "bray")

nmds_abb <- metaMDS(fish_dist_abb, k= 3, trymax=1000, plot = FALSE)

points_abb <- data.frame(nmds_abb$points) %>% 
  cbind(fish_L_abb %>% select(site:year)) %>% 
  mutate(yr_abb = str_sub(year, start = 3))

hulls_abb <- points_abb %>%
  group_by(site) %>% 
  slice(chull(MDS1,MDS2))

env.test <- envfit(nmds_abb, fish_L_abb[46:49], permutations = 9999)
env.test 

wq_vecs <- as.data.frame(scores(env.test, "vectors")) * ordiArrowMul(env.test)

wq_vec_df <- wq_vecs %>% 
  rownames_to_column(var = "wq_metric") %>% 
  filter(!wq_metric == c("do_mg_l", "salinity_ppm")) %>% 
  mutate(wq_labels = c("Secchi", "Temp"))

#### Fig 2 ####
nmds_abb_plot <- ggplot() +
  geom_point(data = points_abb, 
             aes(x = MDS1, y = MDS2, 
                 color = site, shape = ipa), 
             size = 2) + 
  geom_text_repel(data = points_abb, 
                  aes(x = MDS1, y = MDS2, color = site, label = yr_abb),
                  max.overlaps = 20, show.legend = FALSE) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               data = wq_vec_df, linewidth =0.75, 
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text(data = wq_vec_df, aes(x = NMDS1, y = NMDS2), colour = "grey30",
            label = wq_vec_df$wq_labels, size = 4, nudge_y = -0.05) +
  scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) +
  theme_classic(base_size = 10) +
  annotate("text", x = 0.2, y = 1.4, size = 3,
           label = paste("Stress = ", round(nmds_abb$stress, digits = 3))) +
  annotate("text", label = "B", 
           x = -2, y = 1.4,
           size = 6, fontface = "bold") +
  labs(color = "Site", shape = "Condition\ncategory", y = "")

nmds_all + nmds_abb_plot + plot_layout(guides = "collect")
# + plot_annotation(tag_levels = 'A')

ggsave(here("figures", "Fig_2.png"), 
       width = 174, height = 100, units = "mm", dpi = 300)
