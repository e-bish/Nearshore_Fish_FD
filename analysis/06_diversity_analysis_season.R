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

fish.list$trait <- fish.list$trait %>% select(!migrations)

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
adonis2(mFD_results[8:11] ~ ipa*site*season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[8:11] ~ ipa*site*season - ipa:site:season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[8:11] ~ ipa+site+season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#site + season

#Species_Richness
adonis2(mFD_results[7] ~ ipa*site*season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[7] ~ ipa*site*season - ipa:site:season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#ipa:site and site:season
adonis2(mFD_results[7] ~ ipa+site+season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#site + season

#FRic
adonis2(mFD_results[10] ~ ipa*site*season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[10] ~ ipa*site*season - ipa:site:season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#site:season
adonis2(mFD_results[10] ~ ipa+site+season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#site + season 

#FEve
adonis2(mFD_results[9] ~ ipa*site*season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[9] ~ ipa*site*season - ipa:site:season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[9] ~ ipa+site+season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")

#FDiv
adonis2(mFD_results[11] ~ ipa*site*season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[11] ~ ipa*site*season - ipa:site:season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[11] ~ ipa+site+season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#site

#FDis
adonis2(mFD_results[8] ~ ipa*site*season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[8] ~ ipa*site*season - ipa:site:season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[8] ~ ipa+site+season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")

#### region permanova ####
#approach #1
adonis2(mFD_results[8:11] ~ ipa*region*season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[8:11] ~ ipa*region*season - ipa:region:season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[8:11] ~ region+ipa+season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#region + season

#Species_Richness
adonis2(mFD_results['Species_Richness'] ~ ipa*region*season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results['Species_Richness'] ~ ipa*region*season - ipa:region:season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results['Species_Richness'] ~ region+ipa+season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#region + season

#FRic
adonis2(mFD_results['FRic'] ~ ipa*region*season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results['FRic'] ~ ipa*region*season - ipa:region:season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results['FRic'] ~ region+ipa+season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#region + season

#FEve
adonis2(mFD_results[9] ~ ipa*region*season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[9] ~ ipa*region*season - ipa:region:season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results['FEve'] ~ region+ipa+season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")

#FDiv
adonis2(mFD_results[11] ~ ipa*region*season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[11] ~ ipa*region*season - ipa:region:season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results['FDiv'] ~ region+ipa+season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")

#FDis
adonis2(mFD_results[8] ~ ipa*region*season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[8] ~ ipa*region*season - ipa:region:season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results['FDis'] ~ region+ipa+season, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
