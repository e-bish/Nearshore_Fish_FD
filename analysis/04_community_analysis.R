library(tidyverse)
library(here)
library(vegan)

load(here("data", "fish_L_full.Rdata"))

fish_L_full <- fish_L_full %>% 
  mutate(ipa = ifelse(ipa == "Natural2", "Natural", ipa))

#### rda for site ####
fish_L_full <- fish_L_full %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")),
         ipa = factor(ipa, levels = c("Natural",
                                      "Armored",
                                      "Restored")))

mod1 <- dbrda(fish_L_full[4:45] ~ ipa + site + year, 
              data = fish_L_full, distance = "bray")

anova(mod1, by = "margin")

pairwise.adonis2(fish_L_full[4:45] ~ site, 
                 data = mFD_results, distance = "bray")

plot(mod1)

rda_scores1 <- scores(mod1)
sites_scores1 <- as.data.frame(rda_scores1$sites)
biplot_scores1 <- as.data.frame(rda_scores1$biplot) %>% 
  rownames_to_column(var = "fac") %>% 
  filter(grepl("ipa", fac))

# gg_ordiplot(ord = biplot_scores1, #for some reason the scale gets weird if you don't specify this
#             groups = mFD_results$site,
#             ellipse = TRUE,
#             hull = FALSE,
#             spiders = FALSE)
# 
# gg_ordiplot(ord = biplot_scores2, #for some reason the scale gets weird if you don't specify this
#             groups = mFD_sub$site,
#             ellipse = TRUE,
#             hull = FALSE,
#             spiders = FALSE)

points1 <- sites_scores1 %>% 
  cbind(fish_L_full %>% select(site:year))

hulls <- points1 %>%
  group_by(site) %>% 
  slice(chull(dbRDA1,dbRDA2))

site_colors <- rev(c("#8c510a","#d8b365", 
                     # "#f6e8c1",
                     "lightgoldenrod",
                     # "#c7eae8",
                     "lightblue",
                     "#5ab4ac", "#01665e"))

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
