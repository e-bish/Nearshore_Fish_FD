library(vegan)
library(tidyverse)
library(janitor)
library(here)
library(ggrepel)


load(here("data", "net_core.Rdata")) #created in "02_tidy_data.R"

set.seed(2025)

spp_sum_df <- net_core %>% 
  group_by(site, ComName) %>% 
  summarize(spp_sum = sum(species_count)) %>% 
  ungroup() %>% 
  filter(!is.na(spp_sum)) %>% 
  pivot_wider(names_from = ComName, values_from = spp_sum, values_fill = 0) %>% 
  # clean_names() %>% 
  select(!site)

#run all of the models
rad <- radfit(spp_sum_df)

#view results
plot(rad)

#check support for using a single model for all sites
rad[["1"]][["models"]]
#preemption (129.165) mandelbrot (133.165), delta AIC = 4

rad[["6"]][["models"]]
#preemption (257.250) mandelbrot (261.25), delta AIC = 4

#using best fit model for each site
FAM_mod <- rad.preempt(spp_sum_df[1,])
TUR_mod <- rad.zipfbrot(spp_sum_df[2,])
COR_mod <- rad.zipfbrot(spp_sum_df[3,])
SHR_mod <- rad.zipfbrot(spp_sum_df[4,])
DOK_mod <- rad.zipfbrot(spp_sum_df[5,])
EDG_mod <- rad.preempt(spp_sum_df[6,])

plot(DOK_mod)


FAM_points <- points(FAM_mod)[[1]] %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "species") %>% 
  mutate(site = "FAM")
TUR_points <- points(TUR_mod)[[1]] %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "species") %>% 
  mutate(site = "TUR")
COR_points <- points(COR_mod)[[1]] %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "species") %>% 
  mutate(site = "COR")
SHR_points <- points(SHR_mod)[[1]] %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "species") %>% 
  mutate(site = "SHR")
DOK_points <- points(DOK_mod)[[1]] %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "species") %>% 
  mutate(site = "DOK")
EDG_points <- points(EDG_mod)[[1]] %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "species") %>% 
  mutate(site = "EDG")

combined_points <- FAM_points %>% 
  bind_rows(TUR_points, COR_points, SHR_points, DOK_points, EDG_points) %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")))

site_colors <- rev(c("#8c510a","#d8b365", 
                     "lightgoldenrod",
                     "lightblue",
                     "#5ab4ac", "#01665e"))

combined_points %>% 
  ggplot(aes(x = rnk, y = poi, shape = site, color = site)) +
  geom_point(size = 3) +
  theme_classic() +
  geom_text_repel(aes(label=ifelse(poi>250,as.character(species),'')),
                  show.legend = FALSE, max.overlaps = 15) + 
  labs(x = "Rank", y = "Abundance", color = "Site", shape = "Site") +
  xlim(0,15) +
  scale_color_manual(values = site_colors)

ggsave(here("figures", "Fig_S3.png"), 
       width = 6, height = 5.5, dpi = 300)
