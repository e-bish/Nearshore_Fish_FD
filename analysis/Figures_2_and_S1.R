library(here)
library(tidyverse)
library(FD)
library(RColorBrewer)

#load abundance and trait data
load(here("data", "fish.list.Rdata")) #last saved 3/21/25

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
trait_cols <- c(bs_col, dp_col, mi_col, fg_col)

#by ipa
fish_cwm_long %>% 
  ggplot(aes(x = ipa, fill = factor(traits, levels = unique(traits)))) +
  geom_bar(stat = "count", position = "fill") + 
  theme_classic() + 
  scale_fill_manual(values = trait_cols) +
  facet_wrap(~trait_group) +
  labs(fill = "traits", y = "proportion of sampling points where dominant")

ipa_trait_wrap <- fish_cwm_long %>% 
  ggplot(aes(x = ipa, fill = factor(traits, levels = unique(traits)))) +
  geom_bar(stat = "count", position = "fill") + 
  theme_classic() + 
  scale_fill_manual(values = trait_cols) +
  theme(strip.background = element_rect(fill = NA, colour = NA),
        strip.text.x = element_text(size = 10)) +
  facet_wrap(~factor(trait_group, 
                     labels = c("Body Shape",
                                "Water Column Position",
                                "Feeding Guild",
                                "Migrations"))) +
  labs(fill = "Traits", y = "Proportion of sampling points where dominant", x = "Condition category")


length_cwm <- fish_cwm_df %>% 
  select(!body_shape_i:feeding_guild)

ipa_length_plot <- length_cwm %>% 
  ggplot(aes(x = ipa, y = mean_length_mm)) +
  geom_boxplot() + 
  geom_point(alpha = 0.4) +
  theme_classic() + 
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 10)) +
  labs(y = "CWM of mean lengths (mm)", title = "Body Size")

ipa_length_plot + ipa_trait_wrap + plot_layout(ncol = 1, heights = c(0.5,1))

ggsave(here("figures", "Figure_2.png"), 
       width = 6, height = 8, dpi = 300)
