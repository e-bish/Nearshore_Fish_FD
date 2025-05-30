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

load(here("data", "fish.list.season.Rdata"))

fish.list <- fish.list.season

fish.list$trait <- fish.list$trait %>% dplyr::select(!migrations)

#### run with mFD ####
traits_cat <- data.frame(trait_name = colnames(fish.list$trait),
                         trait_type = c("Q", "N", "N", "N"))

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

traits_correlations$"tr_faxes_plot"

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

load(here("data", "rows_w_few_spp_season1.Rdata"))

add_small_samples_back <- rows_w_few_spp %>% 
  rownames_to_column(var = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa", "season")) %>% 
  mutate(Species_Richness = rowSums(across(where(is.numeric)))) %>% 
  select(site, ipa, season, Species_Richness)

mFD_results <- mFD_values %>%
  select(!6:8) %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa", "season"), cols_remove = TRUE) %>% 
  full_join(add_small_samples_back) %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")),
         region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"), 
         veg = ifelse(site %in% c("TUR", "COR", "SHR"), "present", "absent"), .after = site,
         shoreline = paste0(site, case_when(ipa == "Armored" ~ "A", 
                                            ipa == "Restored" ~ "R", 
                                            ipa == "Natural2" ~ "N2",
                                            TRUE ~ "N")),
         season = factor(season, levels = c("Apr-May", "Jun-Jul", "Aug-Sept")),
         across(where(is.numeric), ~replace_na(., 0))) %>% 
  mutate(ipa = ifelse(ipa == "Natural2", "Natural", ipa)) 

# save(mFD_results, file = here("data","mFD_results_season.Rdata")) #last saved 5/27/25
load(here("data","mFD_results_season.Rdata"))

#plot with full range
plot_prep <- mFD_results %>%
  mutate(ipa2 = ifelse(shoreline == "TURN2", "alt", "no_alt")) %>% 
  group_by(site, ipa, ipa2) %>% 
  mutate(id = factor(cur_group_id())) %>%
  ungroup() %>% 
  pivot_longer(!c(shoreline, site, ipa, ipa2, season, region, veg, id), 
               names_to = "metric", values_to = "value") %>%
  group_by(id, metric) %>% 
  mutate(min = min(value), max = max(value)) %>% 
  ungroup() 


p1 <- plot_prep %>% 
  filter(metric == "Species_Richness") %>% 
  ggplot(aes(x = id, y = value, color = site, shape = ipa), show.legend = TRUE) +
  geom_point(size = 3, alpha = 0.5) +
  geom_linerange(aes(ymin=min,ymax=max),linetype=1, show.legend = FALSE)+
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        # axis.text.x = element_text(angle = 60, vjust = 0.75, hjust=1),
        legend.box = "horizontal", 
        legend.title = element_text(hjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 10)) +
  labs(title = "Species Richness", 
       # x = "Shoreline ID", 
       y = "Value", color = "Site", shape = "Shoreline\ntype") 

p2 <- plot_prep %>% 
  filter(!metric == "Species_Richness") %>% 
  ggplot(aes(x = id, y = value, color = site, shape = ipa), show.legend = TRUE) +
  geom_point(size = 3, alpha = 0.5) +
  geom_linerange(aes(ymin=min,ymax=max),linetype=1, show.legend = FALSE)+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.75, hjust=1),
        legend.box = "horizontal", 
        legend.title = element_text(hjust = 0.5),
        strip.placement = "outside",
        strip.text.x = element_text(size = 10),
        strip.background = element_rect(fill = NA,
                                        colour = NA)) +
  labs(x = "Shoreline ID", y = "Value", color = "Site", shape = "Shoreline\ntype") +
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

#### Figure 3 ###
# plot with means
mFD_averages <- mFD_results %>% 
  pivot_longer(!c(shoreline, site, ipa, season, region, veg), names_to = "metric", values_to = "value") %>% 
  group_by(site, ipa, shoreline, metric) %>% 
  summarize(mean = mean(value), sd = sd(value), max = max(value), min = min(value))

plot_prep <- mFD_results %>%
  mutate(ipa2 = ifelse(shoreline == "TURN2", "alt", "no_alt")) %>% 
  group_by(site, ipa, ipa2) %>% 
  mutate(id = factor(cur_group_id())) %>%
  ungroup() %>% 
  pivot_longer(!c(shoreline, site, ipa, ipa2, season, region, veg, id), 
               names_to = "metric", values_to = "value") %>%
  group_by(site, ipa, id, metric) %>% 
  summarize(avg = mean(value), min = min(value), max = max(value)) %>% 
  ungroup() 

site_colors <- rev(c("#8c510a","#d8b365", 
                     # "#f6e8c1",
                     "lightgoldenrod",
                     # "#c7eae8",
                     "lightblue",
                     "#5ab4ac", "#01665e"))

p1 <- plot_prep %>% 
  filter(metric == "Species_Richness") %>% 
  ggplot(aes(x = id, y = avg, color = site, shape = ipa), show.legend = TRUE) +
  geom_point(size = 3) +
  geom_linerange(aes(ymin=min,ymax=max),linetype=1, show.legend = FALSE)+
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        # axis.text.x = element_text(angle = 60, vjust = 0.75, hjust=1),
        legend.box = "horizontal", 
        legend.title = element_text(hjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 10)) +
  scale_color_manual(values = site_colors) +
  labs(title = "Species Richness", 
       # x = "Shoreline ID", 
       y = "Value", color = "Site", shape = "Shoreline\ntype") 

p2 <- plot_prep %>% 
  filter(!metric == "Species_Richness") %>% 
  ggplot(aes(x = id, y = avg, color = site, shape = ipa), show.legend = TRUE) +
  geom_point(size = 3) +
  geom_linerange(aes(ymin=min,ymax=max),linetype=1, show.legend = FALSE)+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, vjust = 0.75, hjust=1),
        legend.box = "horizontal", 
        legend.title = element_text(hjust = 0.5),
        strip.placement = "outside",
        strip.text.x = element_text(size = 10),
        strip.background = element_rect(fill = NA,
                                        colour = NA)) +
  scale_color_manual(values = site_colors) +
  labs(x = "Shoreline ID", y = "Value", color = "Site", shape = "Shoreline\ntype") +
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
  pivot_longer(!c(site, ipa, season, shoreline, region, veg), names_to = "metric", values_to = "value") %>% 
  ggplot(aes(x = ipa, y = value)) +
  geom_boxplot() +
  geom_point(show.legend = FALSE) +
  theme_classic() +
  facet_wrap(~factor(metric, levels = c("Species_Richness", "FDis", "FRic", "FEve", "FDiv"),
                     labels = c("Species Richness", "FDis", "FRic", "FEve", "FDiv")),
             scales = "free_y") + 
  labs(x = "Condition category", y = "Value") + 
  theme(strip.background = element_rect(fill = NA, colour = NA))

#region
mFD_results %>% 
  pivot_longer(!c(site, ipa, season, shoreline, region, veg), names_to = "metric", values_to = "value") %>% 
  ggplot(aes(x = region, y = value)) +
  geom_boxplot() +
  geom_point(show.legend = FALSE) +
  theme_classic() +
  facet_wrap(~factor(metric, levels = c("Species_Richness", "FDis", "FRic", "FEve", "FDiv"),
                     labels = c("Species Richness", "FDis", "FRic", "FEve", "FDiv")),
             scales = "free_y") + 
  labs(x = "Season", y = "Value") + 
  theme(strip.background = element_rect(fill = NA, colour = NA))

#season
mFD_results %>% 
  pivot_longer(!c(site, ipa, season, shoreline, region, veg), names_to = "metric", values_to = "value") %>% 
  ggplot(aes(x = season, y = value)) +
  geom_boxplot() +
  geom_point(show.legend = FALSE) +
  theme_classic() +
  facet_wrap(~factor(metric, levels = c("Species_Richness", "FDis", "FRic", "FEve", "FDiv"),
                     labels = c("Species Richness", "FDis", "FRic", "FEve", "FDiv")),
             scales = "free_y") + 
  labs(x = "Season", y = "Value") + 
  theme(strip.background = element_rect(fill = NA, colour = NA))

#### permanova ####
#test ipa and site
adonis2(mFD_results[8:11] ~ ipa*site*season, data = mFD_results, permutations = 9999, by = "terms", method = "euclidean")
adonis2(mFD_results[8:11] ~ ipa*site*season - ipa:site:season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[8:11] ~ ipa+site+season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#site + season

adonis2(mFD_results['Species_Richness'] ~ ipa*site, data = mFD_results, permutations = 9999, by = "terms", method = "euclidean")
adonis2(mFD_results['Species_Richness'] ~ site*ipa, data = mFD_results, permutations = 9999, by = "terms", method = "euclidean")

#Species_Richness
adonis2(mFD_results[7] ~ season*site*ipa, data = mFD_results, permutations = 9999, by = "terms", method = "euclidean")
adonis2(mFD_results[7] ~ season*ipa*site, data = mFD_results, permutations = 9999, by = "terms", method = "euclidean")
#site; order does not impact our conclusions


# adonis2(mFD_results[7] ~ ipa*site*season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
# adonis2(mFD_results[7] ~ ipa*site*season - ipa:site:season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
# #ipa:site and site:season
# adonis2(mFD_results[7] ~ ipa+site+season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#site + season

#FRic
adonis2(mFD_results['FRic'] ~ season*site*ipa, data = mFD_results, permutations = 9999, by = "terms", method = "euclidean")
adonis2(mFD_results['FRic'] ~ season*ipa*site, data = mFD_results, permutations = 9999, by = "terms", method = "euclidean")
#site; order slightly affects ipa conclusions


# adonis2(mFD_results[10] ~ ipa*site*season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
# adonis2(mFD_results[10] ~ ipa*site*season - ipa:site:season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
# #site:season
# adonis2(mFD_results[10] ~ ipa+site+season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#site + season 

#FEve
adonis2(mFD_results['FEve'] ~ season*site*ipa, data = mFD_results, permutations = 9999, by = "terms", method = "euclidean")
adonis2(mFD_results['FEve'] ~ season*ipa*site, data = mFD_results, permutations = 9999, by = "terms", method = "euclidean")
#none, order doesn't affect conclusions

# adonis2(mFD_results[9] ~ ipa*site*season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
# adonis2(mFD_results[9] ~ ipa*site*season - ipa:site:season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
# adonis2(mFD_results[9] ~ ipa+site+season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")

#FDiv
adonis2(mFD_results['FDiv'] ~ season*site*ipa, data = mFD_results, permutations = 9999, by = "terms", method = "euclidean")
adonis2(mFD_results['FDiv'] ~ season*ipa*site, data = mFD_results, permutations = 9999, by = "terms", method = "euclidean")
#order doesn't affect conclusions

# adonis2(mFD_results[11] ~ ipa*site*season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
# adonis2(mFD_results[11] ~ ipa*site*season - ipa:site:season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
# adonis2(mFD_results[11] ~ ipa+site+season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#site

#FDis
adonis2(mFD_results['FDis'] ~ season*site*ipa, data = mFD_results, permutations = 9999, by = "terms", method = "euclidean")
adonis2(mFD_results['FDis'] ~ season*ipa*site, data = mFD_results, permutations = 9999, by = "terms", method = "euclidean")
#order doesn't affect conclusions

# adonis2(mFD_results[8] ~ ipa*site*season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
# adonis2(mFD_results[8] ~ ipa*site*season - ipa:site:season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
# adonis2(mFD_results[8] ~ ipa+site+season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")

#### region permanova ###
whole_plot_shuffle <- how(within = Within(type = "none"),
                          plots = Plots(strata = mFD_results$site, type = "free"),
                          nperm = 999)

mod <- rda(mFD_results[8:11] ~ region+veg, data = mFD_results) #condition by site doesn't work here
anova(mod, permutations = whole_plot_shuffle, by = "margin")


adonis2(mFD_results[8:11] ~ region:veg, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")
adonis2(mFD_results[8:11] ~ region+veg, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")
#region

#Species_Richness
adonis2(mFD_results['Species_Richness'] ~ region*veg, data = mFD_results, permutations = whole_plot_shuffle, by = "terms", method = "euclidean")
adonis2(mFD_results['Species_Richness'] ~ veg*region, data = mFD_results, permutations = whole_plot_shuffle, by = "terms", method = "euclidean")


# adonis2(mFD_results['Species_Richness'] ~ region*veg, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")
# adonis2(mFD_results['Species_Richness'] ~ region+veg, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")


#FRic
adonis2(mFD_results['FRic'] ~ region*veg, data = mFD_results, permutations = whole_plot_shuffle, by = "terms", method = "euclidean")
adonis2(mFD_results['FRic'] ~ veg*region, data = mFD_results, permutations = whole_plot_shuffle, by = "terms", method = "euclidean")

# adonis2(mFD_results['FRic'] ~ region*veg, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")
# adonis2(mFD_results['FRic'] ~ region+veg, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")


#FEve
adonis2(mFD_results['FEve'] ~ region*veg, data = mFD_results, permutations = whole_plot_shuffle, by = "terms", method = "euclidean")
adonis2(mFD_results['FEve'] ~ veg*region, data = mFD_results, permutations = whole_plot_shuffle, by = "terms", method = "euclidean")


# adonis2(mFD_results[9] ~ region:veg, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")
# adonis2(mFD_results['FEve'] ~ region+veg, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")


#FDiv
adonis2(mFD_results['FDiv'] ~ region*veg, data = mFD_results, permutations = whole_plot_shuffle, by = "terms", method = "euclidean")
adonis2(mFD_results['FDiv'] ~ veg*region, data = mFD_results, permutations = whole_plot_shuffle, by = "terms", method = "euclidean")

# adonis2(mFD_results[11] ~ region:veg, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")
# adonis2(mFD_results['FDiv'] ~ region+veg, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")

#FDis
adonis2(mFD_results['FDis'] ~ region*veg, data = mFD_results, permutations = whole_plot_shuffle, by = "terms", method = "euclidean")
adonis2(mFD_results['FDis'] ~ veg*region, data = mFD_results, permutations = whole_plot_shuffle, by = "terms", method = "euclidean")


# adonis2(mFD_results[8] ~ region:veg, data = mFD_results, permutations =whole_plot_shuffle, by = "margin", method = "euclidean")
# adonis2(mFD_results['FDis'] ~ region+veg, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")


########## compare sites
ipa_shuffle <- how(within = Within(type = "free"),
                      plots = Plots(strata = mFD_results$ipa, type = "none"),
                      nperm = 9999)

mod <- rda(mFD_results['Species_Richness'] ~ site + Condition(ipa), data = mFD_results)
plot(mod)
summary(mod, display = NULL)
res <- anova(mod, permutations = ipa_shuffle, by = "margin")

#### feeling good about this ####
season_shuffle <- how(within = Within(type = "none"),
                      plots = Plots(strata = mFD_results$season, 
                                    type = "free"),
                      block = mFD_results$shoreline,
                      nperm = 9999)

#condition removes the between block variation 
mod <- rda(mFD_results[8:11] ~ season + Condition(shoreline), data = mFD_results) 
summary(mod)
anova(mod, permutations = season_shuffle, by = "margin")

mod <- rda(mFD_results['Species_Richness'] ~ season + Condition(shoreline), data = mFD_results) 
anova(mod, permutations = season_shuffle, by = "margin")

mod <- rda(mFD_results['FRic'] ~ season + Condition(shoreline), data = mFD_results) 
anova(mod, permutations = season_shuffle, by = "margin")

mod <- rda(mFD_results['FEve'] ~ season + Condition(shoreline), data = mFD_results) 
anova(mod, permutations = season_shuffle, by = "margin")

mod <- rda(mFD_results['FDiv'] ~ season + Condition(shoreline), data = mFD_results) 
anova(mod, permutations = season_shuffle, by = "margin")

mod <- rda(mFD_results['FDis'] ~ season + Condition(shoreline), data = mFD_results) 
anova(mod, permutations = season_shuffle, by = "margin")

####

## RDA
mod <- rda(mFD_results[8:11] ~ ipa + site + season, data = mFD_results)
plot(mod)

rda_scores <- scores(mod)
sites_scores <- as.data.frame(rda_scores[[1]])
biplot_scores <- as.data.frame(rda_scores[[2]])

RDA_ordiplot <- gg_ordiplot(ord = biplot_scores, #for some reason the scale gets weird if you don't specify this
                            groups = mFD_results$site,
                            ellipse = TRUE,
                            hull = FALSE,
                            spiders = FALSE)

points <- biplot_scores %>% 
  cbind(mFD_results %>% dplyr::select(site:season))

ggplot(data = points, aes(x = RDA1, y = RDA2)) +
  geom_point(aes(color = site, shape = season), size = 3) + 
  stat_ellipse(aes(group = site, color = site),
               linetype = "dashed", show.legend = FALSE) +
  # geom_text_repel(data = points, aes(x = MDS1, y = MDS2, color = site, label = month)) +
  # geom_polygon(data = hulls, aes(x = MDS1, y = MDS2, fill = site), alpha = 0.2) +
  # annotate("text", x = Inf, y = Inf, label = paste("stress = ", S), vjust = 2, hjust = 2) +
  # annotate("text", x = Inf, y = Inf, label = paste("k = ", K), vjust = 2, hjust = 2) +
  theme_minimal() 

## play with nmds
FD.nmds <- metaMDS(mFD_results[8:11], distance = "euc", 
                          autotransform = FALSE,
                          engine = "monoMDS",
                          k = 3,
                          weakties = TRUE,
                          model = "global",
                          maxit = 400,
                          try = 40,
                          trymax = 100, 
                          trace = FALSE)

nmds_points <- data.frame(FD.nmds$points)
nmds_points <- bind_cols(mFD_results[1:6], nmds_points) 

ggplot(data=nmds_points,
       aes(x=MDS1, y=MDS2)) + 
  geom_point(size = 3) +
  stat_ellipse(aes(group = season, color = season),
               linetype = "dashed", show.legend = TRUE) +
  theme(axis.line = element_blank(), 
        axis.ticks = element_blank(),
        axis.text =  element_blank()) +
  theme_classic() +
  labs(fill = "Site", shape = "Year") + 
  # guides(fill=guide_legend(override.aes=list(color=c(depth_colors)))) +
  # annotate("text", x = -1.2, y = 1.4, 
  #          label = paste("Stress = ", round(nonmotile.nmds$stress, 3))) +
  theme(text = element_text(size = 14))

### temporal beta diversity ###
library(adespatial)

apr_may_mat <- fish.list$abund %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa", "season"), cols_remove = TRUE) %>% 
  filter(season == "Apr-May") %>% 
  add_row(site = "EDG", ipa = "Armored") %>% 
  add_row(site = "TUR", ipa = "Natural") %>% 
  arrange(site, ipa) %>% 
  mutate(across(where(is.numeric), ~replace_na(., 0))) %>% 
  dplyr::select(!site:season) %>% 
  decostand(method = "pa")

jun_jul_mat <- fish.list$abund %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa", "season"), cols_remove = TRUE) %>% 
  filter(season == "Jun-Jul") %>% 
  arrange(site, ipa) %>% 
  dplyr::select(!site:season) %>%
  decostand(method = "pa")
  

aug_sept_mat <- fish.list$abund %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa", "season"), cols_remove = TRUE) %>% 
  filter(season == "Aug-Sept") %>% 
  add_row(site = "EDG", ipa = "Armored") %>% 
  add_row(site = "EDG", ipa = "Natural") %>% 
  add_row(site = "DOK", ipa = "Natural") %>% 
  mutate(across(where(is.numeric), ~replace_na(., 0))) %>% 
  arrange(site, ipa) %>% 
  dplyr::select(!site:season) %>% 
  decostand(method = "pa")

step_1 <- TBI(apr_may_mat, jun_jul_mat, method = "%difference", nperm = 999, test.t.perm = TRUE)
summary(step_1)
step_1$BCD.mat

step_2 <- TBI(jun_jul_mat, aug_sept_mat, method = "%difference", nperm = 999, test.t.perm = TRUE)
summary(step_2)
step_2$BCD.mat
step_2$t.test_B.C

apr_may_mat2 <- fish.list$abund %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa", "season"), cols_remove = TRUE) %>% 
  filter(season == "Apr-May") %>% 
  add_row(site = "TUR", ipa = "Natural") %>% 
  arrange(site, ipa) %>% 
  mutate(across(where(is.numeric), ~replace_na(., 0))) %>% 
  dplyr::select(!site:season) %>% 
  decostand(method = "pa")

aug_sept_mat2 <- fish.list$abund %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa", "season"), cols_remove = TRUE) %>% 
  filter(season == "Aug-Sept") %>% 
  add_row(site = "EDG", ipa = "Natural") %>% 
  add_row(site = "DOK", ipa = "Natural") %>% 
  mutate(across(where(is.numeric), ~replace_na(., 0))) %>% 
  arrange(site, ipa) %>% 
  dplyr::select(!site:season) %>% 
  decostand(method = "pa")

step_3 <- TBI(aug_sept_mat2, apr_may_mat2, method = "%difference", nperm = 999, test.t.perm = TRUE)
summary(step_3)
step_3$BCD.mat
step_3$t.test_B.C
