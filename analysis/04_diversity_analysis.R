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

set.seed(2025)

load(here("data", "fish.list.Rdata"))

fish.list$trait <- fish.list$trait %>%
  select(!migrations)

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

# save(mFD_results, file = here("data","mFD_results.Rdata")) #last saved 5/27/25 without migrations
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

#### FD by site ####
mFD_results %>% 
  pivot_longer(!c(site, ipa, year, shoreline, region, veg), names_to = "metric", values_to = "value") %>% 
  ggplot(aes(x = site, y = value)) +
  geom_boxplot() +
  geom_point(show.legend = FALSE) +
  theme_classic() +
  facet_wrap(~factor(metric, levels = c("Species_Richness", "FDis", "FRic", "FEve", "FDiv"),
                     labels = c("Species Richness", "FDis", "FRic", "FEve", "FDiv")),
             scales = "free_y") + 
  labs(x = "Condition category", y = "Value") + 
  theme(strip.background = element_rect(fill = NA, colour = NA))


#### figure S2, FD by ipa ####
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

#### FD by region ####
mFD_results %>% 
  pivot_longer(!c(site, ipa, year, shoreline, region, veg), names_to = "metric", values_to = "value") %>% 
  ggplot(aes(x = region, y = value)) +
  geom_boxplot() +
  geom_point(show.legend = FALSE) +
  theme_classic() +
  facet_wrap(~factor(metric, levels = c("Species_Richness", "FDis", "FRic", "FEve", "FDiv"),
                     labels = c("Species Richness", "FDis", "FRic", "FEve", "FDiv")),
             scales = "free_y") + 
  labs(x = "Region", y = "Value") + 
  theme(strip.background = element_rect(fill = NA, colour = NA))

#### FD by year ####
mFD_results %>% 
  pivot_longer(!c(site, ipa, year, shoreline, region, veg), names_to = "metric", values_to = "value") %>% 
  ggplot(aes(x = year, y = value)) +
  geom_boxplot() +
  geom_point(show.legend = FALSE) +
  theme_classic() +
  facet_wrap(~factor(metric, levels = c("Species_Richness", "FDis", "FRic", "FEve", "FDiv"),
                     labels = c("Species Richness", "FDis", "FRic", "FEve", "FDiv")),
             scales = "free_y") + 
  labs(x = "Condition category", y = "Value") + 
  theme(strip.background = element_rect(fill = NA, colour = NA))

#### permanova ####

adonis2(mFD_results[8:11] ~ ipa*site*year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[8:11] ~ ipa*site*year - ipa:site:year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[8:11] ~ ipa+site+year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#site

#Species_Richness
adonis2(mFD_results[7] ~ ipa*site*year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[7] ~ ipa*site*year - ipa:site:year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#ipa:site
adonis2(mFD_results[7] ~ ipa+site+year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#site + year

#FRic
adonis2(mFD_results[10] ~ ipa*site*year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#significant with migrations
adonis2(mFD_results[10] ~ ipa*site*year - ipa:site:year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[10] ~ ipa+site+year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#site + year 

#FEve
adonis2(mFD_results[9] ~ ipa*site*year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[9] ~ ipa*site*year - ipa:site:year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[9] ~ ipa+site+year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")

#FDiv
adonis2(mFD_results[11] ~ ipa*site*year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[11] ~ ipa*site*year - ipa:site:year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[11] ~ ipa+site+year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")

#FDis
adonis2(mFD_results[8] ~ ipa*site*year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[8] ~ ipa*site*year - ipa:site:year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[8] ~ ipa+site+year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#site

#pairwise tests
library(pairwiseAdonis)

pairwise.adonis2(mFD_results['Species_Richness'] ~ site, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean" )
pairwise.adonis2(mFD_results['Species_Richness'] ~ year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean" )
pairwise.adonis2(mFD_results['FDis'] ~ site, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean" )


#### permdisp ####
FD_dist <- vegdist(mFD_results[8:11], method = "euclidean")

ipa.disp <- betadisper(FD_dist, mFD_results$ipa, type = c("median"), bias.adjust = FALSE,
                       sqrt.dist = FALSE, add = FALSE)

boxplot(ipa.disp)
permutest(ipa.disp, pairwise = TRUE, permutations = 9999)


site.disp <- betadisper(FD_dist, mFD_results$site, type = c("median"), bias.adjust = FALSE,
                        sqrt.dist = FALSE, add = FALSE)

boxplot(site.disp)
permutest(site.disp, pairwise = TRUE, permutations = 9999)

#### rda for site ####
mod <- rda(mFD_results['FRic'] ~ ipa + site + year, data = mFD_results)
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
  cbind(mFD_results %>% select(site:year))

ggplot(data = points, aes(x = RDA1, y = RDA2)) +
  geom_point(aes(color = ipa, shape = site), size = 3) + 
  # stat_ellipse(aes(group = site, color = site), 
  #              linetype = "dashed", show.legend = FALSE) +
  # geom_text_repel(data = points, aes(x = MDS1, y = MDS2, color = site, label = month)) +
  # geom_polygon(data = hulls, aes(x = MDS1, y = MDS2, fill = site), alpha = 0.2) +
  # annotate("text", x = Inf, y = Inf, label = paste("stress = ", S), vjust = 2, hjust = 2) +
  # annotate("text", x = Inf, y = Inf, label = paste("k = ", K), vjust = 2, hjust = 2) +
  theme_minimal() 

#### permanova region ####
#approach #1
adonis2(mFD_results[8:11] ~ ipa*region*year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[8:11] ~ ipa*region*year - ipa:region:year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[8:11] ~ ipa+region+year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#region

#Species_Richness
adonis2(mFD_results[7] ~ ipa*region*year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[7] ~ ipa*region*year - ipa:region:year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[7] ~ ipa+region+year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#region

#FRic
adonis2(mFD_results[10] ~ ipa*region*year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[10] ~ ipa*region*year - ipa:region:year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[10] ~ ipa+region+year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#region

#FEve
adonis2(mFD_results[9] ~ ipa*region*year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[9] ~ ipa*region*year - ipa:region:year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[9] ~ ipa+region+year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#region

#FDiv
adonis2(mFD_results[11] ~ ipa*region*year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[11] ~ ipa*region*year - ipa:region:year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[11] ~ ipa+region+year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#region

#FDis
adonis2(mFD_results[8] ~ ipa*region*year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[8] ~ ipa*region*year - ipa:region:year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results[8] ~ ipa+region+year, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")


#approach #2 (similar results to approach 1)
mFD_results.region <- mFD_results %>% 
  pivot_longer(!c(shoreline, site, ipa, year, region, veg), 
               names_to = "metric", values_to = "value") %>% 
  group_by(site, region, year, metric) %>% 
  summarize(site_avg = mean(value)) %>% 
  pivot_wider(names_from = metric, values_from = site_avg)

adonis2(mFD_results.region[4:7] ~ region*year, data = mFD_results.region, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results.region[4:7] ~ region+year, data = mFD_results.region, permutations = 9999, by = "margin", method = "euclidean")
#region

#Species_Richness
adonis2(mFD_results.region['Species_Richness'] ~ region*year, data = mFD_results.region, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results.region['Species_Richness'] ~ region+year, data = mFD_results.region, permutations = 9999, by = "margin", method = "euclidean")
#region

#FRic
adonis2(mFD_results.region['FRic'] ~ region*year, data = mFD_results.region, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results.region["FRic"] ~ region+year, data = mFD_results.region, permutations = 9999, by = "margin", method = "euclidean")
#region

#FEve
adonis2(mFD_results.region['FEve'] ~ region*year, data = mFD_results.region, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results.region['FEve'] ~ region+year, data = mFD_results.region, permutations = 9999, by = "margin", method = "euclidean")
#region

#FDiv
adonis2(mFD_results.region['FDiv'] ~ region*year, data = mFD_results.region, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results.region['FDiv'] ~ region+year, data = mFD_results.region, permutations = 9999, by = "margin", method = "euclidean")


#FDis
adonis2(mFD_results.region['FDis'] ~ region*year, data = mFD_results.region, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results.region['FDis'] ~ region+year, data = mFD_results.region, permutations = 9999, by = "margin", method = "euclidean")

## approach #3 (year is more often significant than region)
whole_plot_shuffle <- how(within = Within(type = "none"),
                          plots = Plots(strata = mFD_results$site, type = "free"),
                          nperm = 999) #max is 720 for 6 sites
                          
adonis2(mFD_results[8:11] ~ region*year, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")
adonis2(mFD_results[8:11] ~ region+year+ipa, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")
#year

adonis2(mFD_results['Species_Richness'] ~ region*year, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")
adonis2(mFD_results['Species_Richness'] ~ region+year, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")
#year

adonis2(mFD_results['FRic'] ~ region*year, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")
adonis2(mFD_results['FRic'] ~ region+year, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")


adonis2(mFD_results['FEve'] ~ region*year, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")
adonis2(mFD_results['FEve'] ~ region+year, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")
#region

adonis2(mFD_results['FDiv'] ~ region*year, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")
adonis2(mFD_results['FDiv'] ~ region+year, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")
#year

adonis2(mFD_results['FDis'] ~ region*year, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")
adonis2(mFD_results['FDis'] ~ region+year, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")


#### rda for region ####
mod2 <- rda(mFD_results.region[4:7] ~ region + year, data = mFD_results.region)
plot(mod2)

rda_scores <- scores(mod2)
sites_scores <- as.data.frame(rda_scores[[1]])
biplot_scores <- as.data.frame(rda_scores[[2]])

# RDA_ordiplot <- gg_ordiplot(ord = biplot_scores, #for some reason the scale gets weird if you don't specify this
#                             groups = mFD_results$site,
#                             ellipse = TRUE,
#                             hull = FALSE,
#                             spiders = FALSE)

points <- biplot_scores %>% 
  cbind(mFD_results.region %>% select(site:year))

ggplot(data = points, aes(x = RDA1, y = RDA2, color = region)) +
  geom_text(size = 3, aes(label = site)) + 
  stat_ellipse(aes(group = region, color = region), 
               linetype = "dashed", show.legend = FALSE) +
  # geom_text_repel(data = points, aes(x = MDS1, y = MDS2, color = site, label = month)) +
  # geom_polygon(data = hulls, aes(x = MDS1, y = MDS2, fill = site), alpha = 0.2) +
  # annotate("text", x = Inf, y = Inf, label = paste("stress = ", S), vjust = 2, hjust = 2) +
  # annotate("text", x = Inf, y = Inf, label = paste("k = ", K), vjust = 2, hjust = 2) +
  theme_minimal() 

region_dist <- vegdist(mFD_results.region[4:8], method = "euclidean")

region.disp <- betadisper(region_dist, mFD_results.region$region, type = c("median"), bias.adjust = FALSE,
                          sqrt.dist = FALSE, add = FALSE)

boxplot(region.disp)
permutest(region.disp, pairwise = TRUE, permutations = 9999)
