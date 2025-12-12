library(here)
library(tidyverse)
library(FD)
library(RColorBrewer)
library(patchwork)

#load abundance and trait data
load(here("data", "fish.list.Rdata")) #created in "03_create_matrices.R"

SOS_core_sites <- factor(c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"), 
                         levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"))

#### look at traits by shoreline ####
fish_cwm <- functcomp(x = fish.list$trait,
                      a = as.matrix(fish.list$abund),
                      CWM.type = "dom")

fish_cwm_df <- fish_cwm %>% 
  rownames_to_column( "sample") %>% 
  separate_wider_delim(cols = sample, delim = "_", names = c("site", "ipa", "year")) %>% 
  mutate(site = factor(site, levels = SOS_core_sites),
         region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"), 
         veg = ifelse(site %in% c("TUR", "COR", "SHR"), "present", "absent"), .after = site,
         ipa = ifelse(ipa == "Natural2", "Natural", ipa), #combine the two natural sites at TUR 
         shoreline = paste0(site, case_when(ipa == "Armored" ~ "A", 
                                            ipa == "Restored" ~ "R", 
                                            ipa == "Natural2" ~ "N2",
                                            TRUE ~ "N")))

fish_cwm_long <- fish_cwm_df %>% 
  select(!mean_length_mm) %>% 
  pivot_longer(!c(site, ipa, year, region, shoreline, veg), names_to = "trait_group", values_to = "traits") %>% 
  mutate(across(where(is.character), as.factor)) %>% 
  arrange(trait_group, traits)

bs_col <- brewer.pal(n = 4, "Blues")
dp_col <- brewer.pal(n = 4, "YlGn")
mi_col <- brewer.pal(n = 3, "Reds")
fg_col <- brewer.pal(n = 3, "RdPu")
trait_cols <- c(bs_col, dp_col,
                mi_col,
                fg_col)

#by site
trait_wrap_plot <- fish_cwm_long %>% 
  ggplot(aes(x = site, fill = factor(traits, levels = unique(traits)))) +
  geom_bar(stat = "count", position = "fill") + 
  theme_classic() + 
  scale_fill_manual(values = trait_cols) +
  theme(strip.background = element_rect(fill = NA, colour = NA),
        strip.text.x = element_text(size = 10)) +
  facet_wrap(~factor(trait_group, 
                     labels = c("Transverse shape",
                                "Vertical distribition",
                                "Feeding guild", 
                                "Migration behavior"))) +
  labs(fill = "Trait value", y = "Relative abundance of\nCWM trait values", x = "Site")

length_cwm <- fish_cwm_df %>% 
  select(!body_shape_i:feeding_guild)

length_plot <- length_cwm %>% 
  ggplot(aes(x = site, y = mean_length_mm)) +
  geom_boxplot() + 
  geom_point(alpha = 0.4) +
  theme_classic() + 
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10)) +
  labs(y = "CWM of\nmean lengths (mm)", title = "Size")

length_plot + trait_wrap_plot + plot_layout(ncol = 1, heights = c(0.5,1)) + 
  plot_layout(guides = "collect") 

ggsave(here("figures", "Fig_3.png"), 
       width = 174, height = 160, units = "mm", dpi = 300)

#by ipa
ipa_trait_wrap <- fish_cwm_long %>% 
  ggplot(aes(x = ipa, fill = factor(traits, levels = unique(traits)))) +
  geom_bar(stat = "count", position = "fill") + 
  theme_classic() + 
  scale_fill_manual(values = trait_cols) +
  theme(strip.background = element_rect(fill = NA, colour = NA),
        strip.text.x = element_text(size = 10)) +
  facet_wrap(~factor(trait_group, 
                     labels = c("Transverse shape",
                                "Vertical distribition",
                                "Feeding guild", 
                                "Migration behavior"))) +
  labs(fill = "Trait value", y = "Relative abundance of\nCWM trait values", x = "Condition Category")

length_cwm <- fish_cwm_df %>% 
  select(!body_shape_i:feeding_guild)

ipa_length_plot <- length_cwm %>% 
  ggplot(aes(x = ipa, y = mean_length_mm)) +
  geom_boxplot() + 
  geom_point(alpha = 0.4) +
  theme_classic() + 
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10)) +
  labs(y = "CWM of\nmean lengths (mm)", title = "Size")

ipa_length_plot + ipa_trait_wrap + plot_layout(ncol = 1, heights = c(0.5,1)) + plot_layout(guides = "collect") 

ggsave(here("figures", "Fig_S3.png"), 
       width = 174, height = 160, units = "mm", dpi = 300)

#by year

year_trait_wrap <- fish_cwm_long %>% 
  ggplot(aes(x = year, fill = factor(traits, levels = unique(traits)))) +
  geom_bar(stat = "count", position = "fill") + 
  theme_classic() + 
  scale_fill_manual(values = trait_cols) +
  theme(strip.background = element_rect(fill = NA, colour = NA),
        strip.text.x = element_text(size = 10)) +
  facet_wrap(~factor(trait_group, 
                     labels = c("Transverse shape",
                                "Vertical distribition",
                                "Feeding guild", 
                                "Migration behavior"))) +
  labs(fill = "Trait value", y = "Relative abundance of\nCWM trait values", x = "Year")


length_cwm <- fish_cwm_df %>% 
  select(!body_shape_i:feeding_guild)

year_length_plot <- length_cwm %>% 
  ggplot(aes(x = year, y = mean_length_mm)) +
  geom_boxplot() + 
  geom_point(alpha = 0.4) +
  theme_classic() + 
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10)) +
  labs(y = "CWM of\nmean lengths (mm)", title = "Size")

year_length_plot + year_trait_wrap + plot_layout(ncol = 1, heights = c(0.5,1)) + plot_layout(guides = "collect") 

ggsave(here("figures", "Fig_S4.png"), 
       width = 174, height = 160, units = "mm", dpi = 300)
