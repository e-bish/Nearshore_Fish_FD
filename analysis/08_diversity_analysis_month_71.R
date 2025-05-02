#FD analysis

#load libraries
library(tidyverse)
library(here)
library(ggrepel)
library(mFD)
library(ggordiplots)
library(patchwork)
library(FD)

#### Start here for FD analysis ####

load(here("data", "fish.list.month.Rdata"))

fish.list <- fish.list.month

#### run with mFD ####
traits_cat <- data.frame(trait_name = colnames(fish.list$trait),
                         trait_type = c("Q", "N", "N", "N", "N"))

#Species trait summary
traits_summary <- sp.tr.summary(tr_cat = traits_cat, 
                                sp_tr = fish.list$trait, 
                                stop_if_NA = T)
traits_summary

#create the trait space
dist_mat <- funct.dist(sp_tr = fish.list$trait, 
                       tr_cat = traits_cat,
                       metric = "gower",
                       weight_type = "equal",
                       stop_if_NA = TRUE)

#examine the quality of the potential functional spaces 
space_quality <- quality.fspaces(sp_dist = dist_mat,
                                 maxdim_pcoa = 10,
                                 deviation_weighting = "absolute",
                                 fdist_scaling = FALSE,
                                 fdendro = "ward.D2")

round(space_quality$"quality_fspaces",3) #lowest value is the best (<.1 is good), meaning species pairs are accurately represented
#aka the distances in euclidean space are accurately reflecting the gowers distances

trait_space <- space_quality$"details_fspaces"$"sp_pc_coord"

#test correlations between traits and axes
traits_correlations <- traits.faxes.cor(
  sp_tr          = fish.list$trait, 
  sp_faxes_coord = trait_space[ , c("PC1", "PC2", "PC3")], 
  plot           = TRUE)

traits_correlations$"tr_faxes_stat"[which(traits_correlations$"tr_faxes_stat"$"p.value" < 0.05), ]

traits_correlations$"tr_faxes_plot" + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#calculate FD Indices
alpha_fd_indices <- alpha.fd.multidim(sp_faxes_coord = trait_space[ , c("PC1", "PC2", "PC3")],
                                      asb_sp_w = data.matrix(fish.list$abund),
                                      ind_vect = c("fdis", "feve", "fric", "fdiv"),
                                      scaling = TRUE,
                                      check_input = TRUE,
                                      details_returned = TRUE)

#the mFD package uses ape::pcoa() which automatically removes negative eigenvalues rather than applying a correction

mFD_values <- alpha_fd_indices$"functional_diversity_indices"

metrics <- c("Species_Richness", "FDis", "FEve", "FRic", "FDiv")
metrics_clean <- c("Species Richness", "FDis", "FEve", "FRic", "FDiv")

colnames(mFD_values)[1:5] <- metrics

# load(here("data", "rows_w_few_spp_month.Rdata"))
# 
# add_small_samples_back <- rows_w_few_spp %>% 
#   rownames_to_column(var = "sample") %>% 
#   separate_wider_delim(sample, delim = "_", names = c("site", "ipa", "month")) %>% 
#   mutate(Species_Richness = rowSums(across(where(is.numeric)))) %>% 
#   select(site, ipa, month, Species_Richness)

mFD_results <- mFD_values %>%
  select(!6:8) %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa", "month"), cols_remove = TRUE) %>% 
  # full_join(add_small_samples_back) %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")),
         region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"), 
         site_pair = case_when(site %in% c("EDG", "SHR") ~ "beach", 
                               site %in% c("DOK", "COR") ~ "marina",
                               site %in% c("TUR", "FAM") ~ "rocky_coast"),
         month = factor(month, levels = c("Apr", "May", "Jun", "Jul", "Aug", "Sept")),
         veg = ifelse(site %in% c("TUR", "COR", "SHR"), "present", "absent"), 
         season = ifelse(month %in% c("May", "Jun", "Jul"), "peak", "shoulder"), .after = site,
         season1 = case_when(month %in% c("Apr", "Sept") ~ "tail",
                             month %in% c("Jun", "Jul") ~ "peak",
                                          TRUE ~ "shoulder"),
         shoreline = paste0(site, case_when(ipa == "Armored" ~ "A", 
                                            ipa == "Restored" ~ "R", 
                                            ipa == "Natural2" ~ "N2",
                                            TRUE ~ "N")),
         across(where(is.numeric), ~replace_na(., 0))) %>% 
  mutate(ipa = ifelse(ipa == "Natural2", "Natural", ipa)) 


# mFD_results.month <- mFD_results

# save(mFD_results.month, file = here("data","mFD_results.month.Rdata")) #last saved 4/29/25
# load(here("data","mFD_results.Rdata"))


effort <- mFD_results %>% 
  group_by(site, ipa) %>% 
  summarize(n())

### Figure 3 ###
site_colors <- rev(c("#8c510a","#d8b365", 
                     # "#f6e8c1",
                     "lightgoldenrod",
                     # "#c7eae8",
                     "lightblue",
                     "#5ab4ac", "#01665e"))

plot_prep <- mFD_results %>%
  mutate(ipa2 = ifelse(shoreline == "TURN2", "alt", "no_alt")) %>% 
  group_by(site, ipa, ipa2) %>% 
  mutate(id = factor(cur_group_id())) %>%
  ungroup() %>% 
  pivot_longer(!c(shoreline, site, region, site_pair, ipa, ipa2, month, season, season1, veg, id), 
               names_to = "metric", values_to = "value") %>%
  group_by(id, metric) %>% 
  mutate(min = min(value), max = max(value)) %>% 
  ungroup()  


p1 <- plot_prep %>% 
  filter(metric == "Species_Richness") %>% 
  # filter(!value == 0) %>% 
  ggplot(aes(x = id, y = value)) +
  geom_boxplot(aes(fill = site)) +
  geom_point(aes(shape = ipa), size = 2, alpha = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        legend.box = "horizontal", 
        legend.title = element_text(hjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 10)) +
  scale_fill_manual(values = site_colors) +
  scale_alpha_manual(values = c(1, 0.6, 0.3)) +
  labs(title = "Species Richness", 
       y = "Value", fill = "Site", shape = "Shoreline\ncondition") 

p2 <- plot_prep %>% 
  filter(!metric == "Species_Richness") %>% 
  # filter(!value == 0) %>% 
  ggplot(aes(x = id, y = value)) +
  geom_boxplot(aes(fill = site)) +
  geom_point(aes(shape = ipa), size = 2, alpha = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.75, hjust=1),
        legend.box = "horizontal", 
        legend.title = element_text(hjust = 0.5),
        strip.placement = "outside",
        strip.text.x = element_text(size = 10),
        strip.background = element_rect(fill = NA,
                                        colour = NA)) +
  scale_fill_manual(values = site_colors) +
  labs(x = "Shoreline ID", y = "Value", fill = "Site", shape = "Shoreline\ncondition") +
  facet_wrap(~factor(metric, levels = c("FDis", "FRic", "FEve", "FDiv")), 
             scales = "free_y", axes = "all_x", axis.labels = "margins")

layout <- '
AB
CC
'

free(p1) + guide_area() + p2 + 
  plot_layout(guides = 'collect', 
              design = layout,
              heights = c(1.1,2)) &
  theme(legend.direction = 'vertical',
        legend.box = 'horizontal')

#### summary stats ####


######################
#### figure S2 ####
mFD_results %>% 
  pivot_longer(!c(site, ipa, month, season, season1, shoreline, region, veg), names_to = "metric", values_to = "value") %>% 
  ggplot(aes(x = ipa, y = value)) +
  geom_boxplot() +
  geom_point(show.legend = FALSE, alpha = 0.5) +
  theme_classic() +
  facet_wrap(~factor(metric, levels = c("Species_Richness", "FDis", "FRic", "FEve", "FDiv"),
                     labels = c("Species Richness", "FDis", "FRic", "FEve", "FDiv")),
             scales = "free_y") + 
  labs(x = "Condition category", y = "Value") + 
  theme(strip.background = element_rect(fill = NA, colour = NA))

#month
mFD_results %>% 
  pivot_longer(!c(site, ipa, month, season, season1, shoreline, region, veg), names_to = "metric", values_to = "value") %>% 
  ggplot(aes(x = month, y = value)) +
  geom_boxplot() +
  geom_point(show.legend = FALSE) +
  theme_classic() +
  facet_wrap(~factor(metric, levels = c("Species_Richness", "FDis", "FRic", "FEve", "FDiv"),
                     labels = c("Species Richness", "FDis", "FRic", "FEve", "FDiv")),
             scales = "free_y") + 
  labs(x = "Month", y = "Value") + 
  theme(strip.background = element_rect(fill = NA, colour = NA))

#### permanova and permdisp ####
#specify the permutations
#in permutation terminology sites are PLOTS and shorelines are SAMPLES
#so to compare sites, we want only whole plots randomized. To compare shorelines, we want only split plots randomized
#specify the permutations
split_plot_shuffle <- how(within = Within(type = "free"),
                          plots = Plots(strata = mFD_results$site, type = "none"),
                          nperm = 9999, 
                          observed = TRUE)

#check permutation structure
#here, within a site we are randomizing ipa and month
# head(mFD_results[, c("site", "ipa", "month")], 15)
# check(mFD_results, control =split_plot_shuffle)
# head(mFD_results[shuffle(nrow(mFD_results), control = split_plot_shuffle), c("site", "ipa", "month")], 15)

FD_dist <- vegdist(mFD_results[,c("FDis", "FEve", "FRic", "FDiv")], method = "euc")

adonis2(FD_dist ~ ipa*site*month, data = mFD_results, permutations = split_plot_shuffle, by = "margin")
#not significant

adonis2(FD_dist ~ ipa*site*month - ipa:site:month, data = mFD_results, permutations = split_plot_shuffle, by = "margin")
#significant interaction of ipa*site, no difference in site*month or ipa*month

adonis2(FD_dist ~ ipa+site*month - site:month, data = mFD_results, permutations = split_plot_shuffle, by = "margin")
#confirmed no ipa differences

adonis2(FD_dist ~ ipa+site+month, data = mFD_results, permutations = split_plot_shuffle, by = "margin")


ipa.disp <- betadisper(FD_dist, mFD_results$ipa, type = c("median"), bias.adjust = FALSE,
                       sqrt.dist = FALSE, add = FALSE)

boxplot(ipa.disp)
permutest(ipa.disp, pairwise = TRUE)
#restored is significantly lower dispersion

#check to see if sites are different if we allow permutations of ipas

#specify the permutations
whole_plot_shuffle <- how(within = Within(type = "none"),
                          plots = Plots(strata = mFD_results$site, type = "free"),
                          nperm = 720, #max
                          observed = TRUE)

#check permutation structure
# head(mFD_results[, c("site", "ipa", "month")], 15)
# check(mFD_results, control =whole_plot_shuffle)
# head(mFD_results[shuffle(nrow(mFD_results), control = whole_plot_shuffle), c("site", "ipa", "month")], 15)

adonis2(FD_dist ~ site * month, #including site and month so we can parse out that variance
        data = mFD_results, permutations = whole_plot_shuffle, by = "margin")
#another, less conservative, approach would be to do unrestricted permutations at this stage because there wasn't a significant difference in 
#ipas. This approach is suggested in Anderson and Braak, 2003. It is, however, debated whether this is an appropriate approach
#and the documentation for PRIMER-E suggests that discarding an aspect of your study design is something to think twice about because it 
#may not be wise to assume there aren't differences just because you didn't find them. 

#again we see a significant interaction between site and month

#### NOT TOTALLY SURE ABOUT THIS CODE YET ####
library(pairwiseAdonis)
pairwise.adonis2(FD_dist ~ month, nperm = 9999, data = mFD_results)
#june is different from every other month
pairwise.adonis2(FD_dist ~ site, nperm = 9999, data = mFD_results)
#cornet is different from any other site

?pairwise.adonis2

####

site.disp <- betadisper(FD_dist, mFD_results$site, type = c("median"), bias.adjust = FALSE,
                        sqrt.dist = FALSE, add = FALSE)

boxplot(site.disp)
permutest(site.disp, pairwise = TRUE)

month.disp <- betadisper(FD_dist, mFD_results$month, type = c("median"), bias.adjust = FALSE,
                        sqrt.dist = FALSE, add = FALSE)

boxplot(month.disp)
permutest(month.disp, pairwise = TRUE)

#no difference in site dispersion, dispersion in May is higher than August

#playing around with ordinations
nmds <- metaMDS(mFD_results[8:11], distance = "euc", 
                                  autotransform = FALSE,
                                  engine = "monoMDS",
                                  k = 3,
                                  weakties = TRUE,
                                  model = "global",
                                  maxit = 400,
                                  try = 40,
                                  trymax = 100, 
                                  trace = FALSE)
plot(nmds)

nmds_points <- data.frame(nmds$points)
nmds_points <- bind_cols(mFD_results[1:6], nmds_points) 

ggplot(data=nmds_points,
       aes(x=MDS1, y=MDS2,
           color= month)) + 
  geom_point(size = 3) +
  # stat_ellipse(aes(group = depth, color = depth), 
  #              linetype = "dashed", show.legend = FALSE) +
  theme(axis.line = element_blank(), 
        axis.ticks = element_blank(),
        axis.text =  element_blank()) +
  theme_classic() +
  # labs(fill = "Site", shape = "Year") + 
  # scale_shape_manual(values = c(21,24,22)) +
  # scale_fill_manual(values = depth_colors) + 
  # scale_color_manual(values = depth_colors) +
  # guides(fill=guide_legend(override.aes=list(color=c(depth_colors)))) +
  # annotate("text", x = -1.2, y = 1.4, 
  #          label = paste("Stress = ", round(nonmotile.nmds$stress, 3))) +
  theme(text = element_text(size = 14))

#just june and september, to compare
plot_prep %>% 
  filter(month %in% c("Jun", "Sept")) %>% 
  ggplot(aes(x = site, y = value, fill = month)) + 
  geom_boxplot() +
  theme_classic() + 
  facet_wrap(~metric, scales = "free_y")

#by season
plot_prep %>% 
  ggplot(aes(x = site, y = value, fill = season1)) + 
  geom_boxplot() +
  theme_classic() + 
  facet_wrap(~metric, scales = "free_y")


#by season and region
plot_prep %>% 
  # filter(!value == 0) %>% 
  ggplot(aes(x = region, y = value, fill = season1)) + 
  geom_boxplot() +
  theme_classic() + 
  facet_wrap(~metric, scales = "free_y")


#### Figure 3 by region ####
p1 <- plot_prep %>% 
  filter(metric == "Species_Richness") %>% 
  ggplot(aes(x = month, y = value, fill = region)) +
  geom_boxplot() +
  # geom_point(aes(shape = ipa), size = 2, alpha = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        legend.box = "horizontal", 
        legend.title = element_text(hjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 10)) +
  scale_fill_manual(values = site_colors) +
  scale_alpha_manual(values = c(1, 0.6, 0.3)) +
  labs(title = "Species Richness", 
       y = "Value", fill = "Region") 

p2 <- plot_prep %>% 
  filter(!metric == "Species_Richness") %>% 
  ggplot(aes(x = month, y = value, fill = region)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.75, hjust=1),
        legend.box = "horizontal", 
        legend.title = element_text(hjust = 0.5),
        strip.placement = "outside",
        strip.text.x = element_text(size = 10),
        strip.background = element_rect(fill = NA,
                                        colour = NA)) +
  scale_fill_manual(values = site_colors) +
  labs(x = "Shoreline ID", y = "Value", fill = "Region") +
  facet_wrap(~factor(metric, levels = c("FDis", "FRic", "FEve", "FDiv")), 
             scales = "free_y", axes = "all_x", axis.labels = "margins")

free(p1) + guide_area() + p2 + 
  plot_layout(guides = 'collect', 
              design = layout,
              heights = c(1.1,2)) &
  theme(legend.direction = 'vertical',
        legend.box = 'horizontal')

plot_prep %>% 
  ggplot(aes(x = month, y = value, color = site, shape = ipa)) + 
  geom_point() + 
  theme_classic() + 
  facet_wrap(~metric, scales = "free_y")

#### Figure 3 by site ipa ####
p1 <- plot_prep %>% 
  filter(metric == "Species_Richness") %>% 
  ggplot(aes(x = site, y = value, fill = ipa)) +
  geom_boxplot() +
  # geom_point(aes(shape = ipa), size = 2, alpha = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        legend.box = "horizontal", 
        legend.title = element_text(hjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 10)) +
  scale_fill_manual(values = site_colors) +
  scale_alpha_manual(values = c(1, 0.6, 0.3)) +
  labs(title = "Species Richness", 
       y = "Value", fill = "Shoreline Condition") 

p2 <- plot_prep %>% 
  filter(!metric == "Species_Richness") %>% 
  ggplot(aes(x = site, y = value, fill = ipa)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.75, hjust=1),
        legend.box = "horizontal", 
        legend.title = element_text(hjust = 0.5),
        strip.placement = "outside",
        strip.text.x = element_text(size = 10),
        strip.background = element_rect(fill = NA,
                                        colour = NA)) +
  scale_fill_manual(values = site_colors) +
  labs(x = "Site", y = "Value", fill = "Shoreline Condition") +
  facet_wrap(~factor(metric, levels = c("FDis", "FRic", "FEve", "FDiv")), 
             scales = "free_y", axes = "all_x", axis.labels = "margins")

free(p1) + guide_area() + p2 + 
  plot_layout(guides = 'collect', 
              design = layout,
              heights = c(1.1,2)) &
  theme(legend.direction = 'vertical',
        legend.box = 'horizontal')


#### Figure 3 by site_pair ipa ####
p1 <- plot_prep %>% 
  filter(metric == "Species_Richness") %>% 
  ggplot(aes(x = site_pair, y = value, fill = ipa)) +
  geom_boxplot() +
  # geom_point(aes(shape = ipa), size = 2, alpha = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        legend.box = "horizontal", 
        legend.title = element_text(hjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 10)) +
  scale_fill_manual(values = site_colors) +
  scale_alpha_manual(values = c(1, 0.6, 0.3)) +
  labs(title = "Species Richness", 
       y = "Value", fill = "Shoreline Condition") 

p2 <- plot_prep %>% 
  filter(!metric == "Species_Richness") %>% 
  ggplot(aes(x = site_pair, y = value, fill = ipa)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.75, hjust=1),
        legend.box = "horizontal", 
        legend.title = element_text(hjust = 0.5),
        strip.placement = "outside",
        strip.text.x = element_text(size = 10),
        strip.background = element_rect(fill = NA,
                                        colour = NA)) +
  scale_fill_manual(values = site_colors) +
  labs(x = "Site Pair", y = "Value", fill = "Shoreline Condition") +
  facet_wrap(~factor(metric, levels = c("FDis", "FRic", "FEve", "FDiv")), 
             scales = "free_y", axes = "all_x", axis.labels = "margins")

free(p1) + guide_area() + p2 + 
  plot_layout(guides = 'collect', 
              design = layout,
              heights = c(1.1,2)) &
  theme(legend.direction = 'vertical',
        legend.box = 'horizontal')

plot_prep %>% 
  ggplot(aes(x = month, y = value, color = site, shape = ipa)) + 
  geom_point() + 
  theme_classic() + 
  facet_wrap(~metric, scales = "free_y")


#SR v FR

mFD_values %>% 
  ggplot(aes(x = Species_Richness, y = FRic)) +
  geom_point() + 
  theme_classic()
