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

load(here("data", "fish.list.Rdata"))

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

load(here("data", "rows_w_few_spp.Rdata"))

add_small_samples_back <- rows_w_few_spp %>% 
  rownames_to_column(var = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa", "year")) %>% 
  mutate(Species_Richness = rowSums(across(where(is.numeric)))) %>% 
  select(site, ipa, year, Species_Richness)

mFD_results <- mFD_values %>%
  select(!6:8) %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa", "year"), cols_remove = TRUE) %>% 
  full_join(add_small_samples_back) %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")),
         region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"), 
         veg = ifelse(site %in% c("TUR", "COR", "SHR"), "present", "absent"), .after = site,
         shoreline = paste0(site, case_when(ipa == "Armored" ~ "A", 
                                            ipa == "Restored" ~ "R", 
                                            ipa == "Natural2" ~ "N2",
                                            TRUE ~ "N")),
         across(where(is.numeric), ~replace_na(., 0))) %>% 
  mutate(ipa = ifelse(ipa == "Natural2", "Natural", ipa)) 

save(mFD_results, file = here("data","mFD_results.Rdata")) #last saved 3/21/25
load(here("data","mFD_results.Rdata"))

#plot with full range
plot_prep <- mFD_results %>%
  mutate(ipa2 = ifelse(shoreline == "TURN2", "alt", "no_alt")) %>% 
  group_by(site, ipa, ipa2) %>% 
  mutate(id = factor(cur_group_id())) %>%
  ungroup() %>% 
  pivot_longer(!c(shoreline, site, ipa, ipa2, year, region, veg, id), 
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
  pivot_longer(!c(shoreline, site, ipa, year, region, veg), names_to = "metric", values_to = "value") %>% 
  group_by(site, ipa, shoreline, metric) %>% 
  summarize(mean = mean(value), sd = sd(value), max = max(value), min = min(value))

plot_prep <- mFD_results %>%
  mutate(ipa2 = ifelse(shoreline == "TURN2", "alt", "no_alt")) %>% 
  group_by(site, ipa, ipa2) %>% 
  mutate(id = factor(cur_group_id())) %>%
  ungroup() %>% 
  pivot_longer(!c(shoreline, site, ipa, ipa2, year, region, veg, id), 
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

ggsave(here("figures", "figure_3.png"), 
       width = 8, height = 6, dpi = 300) 

#### summary stats ####


######################
#### figure S2 ####
mFD_results %>% 
  pivot_longer(!c(site, ipa, year, shoreline, region, veg), names_to = "metric", values_to = "value") %>% 
  ggplot(aes(x = ipa, y = value)) +
  geom_boxplot() +
  geom_point(show.legend = FALSE) +
  theme_classic() +
  facet_wrap(~factor(metric, levels = c("Species_Richness", "FDis", "FRic", "FEve", "FDiv"),
                     labels = c("Species Richness", "FDis", "FRic", "FEve", "FDiv")),
             scales = "free_y") + 
  labs(x = "Condition category", y = "Value") + 
  theme(strip.background = element_rect(fill = NA, colour = NA))

ggsave(here("figures", "figure_S2.png"), 
       width = 8, height = 6, dpi = 300) 


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
# head(mFD_results[, c("site", "ipa", "year")], 15)
# check(mFD_results, control =split_plot_shuffle)
# head(mFD_results[shuffle(nrow(mFD_results), control = split_plot_shuffle), c("site", "ipa", "year")], 15)

FD_dist <- vegdist(mFD_results[,c("FDis", "FEve", "FRic", "FDiv")], method = "euc")

adonis2(FD_dist ~ ipa + site + year, data = mFD_results, permutations = split_plot_shuffle, by = "margin")

ipa.disp <- betadisper(FD_dist, mFD_results$ipa, type = c("median"), bias.adjust = FALSE,
                       sqrt.dist = FALSE, add = FALSE)

boxplot(ipa.disp)
permutest(ipa.disp, pairwise = TRUE)

#check to see if sites are different 

#specify the permutations
whole_plot_shuffle <- how(within = Within(type = "none"),
                          plots = Plots(strata = mFD_results$site, type = "free"),
                          nperm = 999, #max is 720 for our sample size
                          observed = TRUE)

#check permutation structure
# head(mFD_results[, c("site", "ipa", "year")], 15)
# check(mFD_results, control =whole_plot_shuffle)
# head(mFD_results[shuffle(nrow(mFD_results), control = whole_plot_shuffle), c("site", "ipa", "year")], 15)

adonis2(FD_dist ~ site + year, #including site and year so we can parse out that variance
        data = mFD_results, permutations = whole_plot_shuffle, by = "margin")
#another, less conservative, approach would be to do unrestricted permutations at this stage because there wasn't a significant difference in 
#ipas. This approach is suggested in Anderson and Braak, 2003. It is, however, debated whether this is an appropriate approach
#and the documentation for PRIMER-E suggests that discarding an aspect of your study design is something to think twice about because it 
#may not be wise to assume there aren't differences just because you didn't find them. 


site.disp <- betadisper(FD_dist, mFD_results$site, type = c("median"), bias.adjust = FALSE,
                        sqrt.dist = FALSE, add = FALSE)

boxplot(site.disp)
permutest(site.disp, pairwise = TRUE)