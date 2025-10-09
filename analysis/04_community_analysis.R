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
load(here("data", "wq_tb.Rdata")) #created in "Figure_S1_water_quality.R"

site_colors <- rev(c("#8c510a","#d8b365", 
                     "lightgoldenrod",
                     "lightblue",
                     "#5ab4ac", "#01665e"))

#### nmds ####
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
  geom_point(data = points, aes(x = MDS1, y = MDS2, color = site, shape = ipa), size = 3) + 
  geom_text_repel(data = points, 
                  aes(x = MDS1, y = MDS2, color = site, label = yr_abb), 
                  max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) +
  theme_classic() +
  annotate("text", x = -1, y = 1.4, 
           label = paste("Stress = ", round(nmds$stress, 3))) +
  labs(color = "Site", shape = "Condition Category")

#### nmds2 ####

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
nmds_abb <- ggplot() +
  geom_point(data = points_abb, aes(x = MDS1, y = MDS2, color = site, shape = ipa), size = 3) + 
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
  theme_classic() +
  annotate("text", x = -1.4, y = 1.4, 
           label = paste("Stress = ", round(nmds_abb$stress, digits = 3))) +
  labs(color = "Site", shape = "Condition Category", y = "")

nmds_all + nmds_abb + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A')

ggsave(here("figures", "Fig_2.png"), 
       width = 6, height = 5.5, dpi = 300)
