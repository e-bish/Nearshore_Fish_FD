library(tidyverse)
library(here)
library(vegan)
library(pairwiseAdonis)
library(ggrepel)

set.seed(2025)
load(here("data", "net_core.Rdata")) 
load(here("data", "fish_L_full.Rdata"))

site_colors <- rev(c("#8c510a","#d8b365", 
                     # "#f6e8c1",
                     "lightgoldenrod",
                     # "#c7eae8",
                     "lightblue",
                     "#5ab4ac", "#01665e"))

### RAD ###

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

combined_points %>% 
  ggplot(aes(x = rnk, y = poi, shape = site, color = site)) +
  geom_point(size = 3) +
  theme_classic() +
  geom_text_repel(aes(label=ifelse(poi>250,as.character(species),'')),
                  show.legend = FALSE, max.overlaps = 15) + 
  labs(x = "Rank", y = "Abundance", color = "Site", shape = "Site") +
  xlim(0,15) +
  scale_color_manual(values = site_colors)

#### rda for site ####

fish_L_full <- fish_L_full %>% 
  mutate(ipa = ifelse(ipa == "Natural2", "Natural", ipa)) %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")),
         ipa = factor(ipa, levels = c("Natural",
                                      "Armored",
                                      "Restored")))

fish_dist <- vegdist(sqrt(fish_L_full[4:45]), method = "bray")

mod1 <- dbrda(fish_dist ~ ipa + site + year, 
              data = fish_L_full)

anova(mod1, by = "margin", model = "reduced", permutations = 9999)

pairwise.adonis2(fish_dist ~ site, 
                 data = fish_L_full, by = "margin")

site_beta <- betadisper(vegdist(fish_dist), 
                     group = fish_L_full$site, type = "median")

permutest(site_beta, pairwise = TRUE)

#visualize the ordination
rda_scores1 <- scores(mod1)
sites_scores1 <- as.data.frame(rda_scores1$sites)
biplot_scores1 <- as.data.frame(rda_scores1$biplot) %>% 
  rownames_to_column(var = "fac") %>% 
  filter(grepl("ipa", fac))

points1 <- sites_scores1 %>% 
  cbind(fish_L_full %>% select(site:year))

hulls <- points1 %>%
  group_by(site) %>% 
  slice(chull(dbRDA1,dbRDA2))

ggplot(data = points1, aes(x = dbRDA1, y = dbRDA2)) +
  geom_point(aes(color = site, shape = ipa), size = 3) + 
  geom_polygon(data = hulls, 
               aes(x = dbRDA1, 
                   y = dbRDA2, 
                   fill = site), 
               alpha = 0.2) +
  # stat_ellipse(aes(group = site, color = site),
  #              linetype = "dashed", show.legend = FALSE) +
  # geom_segment(data=biplot_scores1,
  #              aes(x = 0, y = 0, xend=dbRDA1, yend=dbRDA2),
  #              arrow=arrow(length = unit(0.01, "npc")),
  #              lwd=0.75) +
  # geom_text(data=biplot_scores1,
  #           aes(x=dbRDA1*0.9,
  #               y=dbRDA2*0.9,
  #               label=fac),
  #           nudge_x = c(-0.2, -0.3, -0.2, -0.2),
  #           nudge_y = c(0.1, -0.05, 0, 0),
  #           size=4) +
  scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) +
  labs(shape = "Shoreline Condition", fill = "Site", color = "Site") +
  theme_classic() 
