#FD analysis

#load libraries
library(tidyverse)
library(here)
library(ggrepel)
library(mFD)
library(ggordiplots)
library(patchwork)
library(FD)
library(pairwiseAdonis)

#### Start here for FD analysis ####

set.seed(2025)

load(here("data", "fish.list.Rdata"))

fish.list$trait <- fish.list$trait %>%
  mutate(migrations = ifelse(migrations == "non-migratory", "resident", "migratory")) %>% 
  mutate(migrations = factor(migrations))

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
  dplyr::select(!6:8) %>% 
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
  mutate(ipa = ifelse(ipa == "Natural2", "Natural", ipa)) %>% 
  arrange(shoreline, year)

# save(mFD_results, file = here("data","mFD_results.Rdata")) #last saved 5/30/25 with simplified migrations
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
# 
# ggsave(here("figures", "figure_3.png"), 
#        width = 8, height = 6, dpi = 300) 

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
  labs(x = "Site", y = "Value") + 
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
# 
# ggsave(here("figures", "figure_S2.png"), 
#        width = 8, height = 6, dpi = 300) 

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

#### FD by veg ####
mFD_results %>% 
  pivot_longer(!c(site, ipa, year, shoreline, region, veg), names_to = "metric", values_to = "value") %>% 
  ggplot(aes(x = veg, y = value)) +
  geom_boxplot() +
  geom_point(show.legend = FALSE) +
  theme_classic() +
  facet_wrap(~factor(metric, levels = c("Species_Richness", "FDis", "FRic", "FEve", "FDiv"),
                     labels = c("Species Richness", "FDis", "FRic", "FEve", "FDiv")),
             scales = "free_y") + 
  labs(x = "Eelgrass", y = "Value") + 
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
plot_shuffle <- how(within = Within(type = "series"),
                          plots = Plots(strata = mFD_results$shoreline, type = "free"),
                          nperm = 9999)

#check permutation structure
head(mFD_results[, c("site", "ipa", "year")], 15)
check(mFD_results, control =plot_shuffle)
head(mFD_results[shuffle(nrow(mFD_results), control = plot_shuffle), c("site", "ipa", "year")], 15)

adonis2(mFD_results[8:11] ~ ipa*site, data = mFD_results, permutations = plot_shuffle, by = "margin", method = "euclidean")
adonis2(mFD_results[8:11] ~ ipa+site, data = mFD_results, permutations = plot_shuffle, by = "margin", method = "euclidean")

FD_dist <- vegdist(mFD_results[8:11], method = "euclidean")

FD.ipa.disp <- betadisper(FD_dist, mFD_results$ipa, type = c("median"), bias.adjust = FALSE,
                          sqrt.dist = FALSE, add = FALSE)

boxplot(FD.ipa.disp)
permutest(FD.ipa.disp, pairwise = TRUE, permutations = plot_shuffle)

## Species_Richness
adonis2(mFD_results['Species_Richness'] ~ ipa*site, data = mFD_results, permutations = plot_shuffle, by = "margin", method = "euclidean")
adonis2(mFD_results['Species_Richness'] ~ ipa+site, data = mFD_results, permutations = plot_shuffle, by = "margin", method = "euclidean")

pairwise.adonis2(mFD_results['Species_Richness'] ~ site, data = mFD_results, permutations = 999, by = "margin", method = "euclidean")
#A: COR
#B: DOK, FAM, SHR, TUR
#C: DOK, EDG, FAM, SHR

# dispersion
SR_dist <- vegdist(mFD_results['Species_Richness'], method = "euclidean")

SR.ipa.disp <- betadisper(SR_dist, mFD_results$ipa, type = c("median"), bias.adjust = FALSE,
                          sqrt.dist = FALSE, add = FALSE)

boxplot(SR.ipa.disp)
permutest(SR.ipa.disp, pairwise = TRUE, permutations = plot_shuffle)

SR.site.disp <- betadisper(SR_dist, mFD_results$site, type = c("median"), bias.adjust = FALSE,
                          sqrt.dist = FALSE, add = FALSE)

boxplot(SR.site.disp)
permutest(SR.site.disp, pairwise = TRUE, permutations = plot_shuffle)
#A: TUR, FAM, SHR, EDG
#B: COR, TUR
#C: FAM, DOK, SHR, EDG

## FRic
adonis2(mFD_results['FRic'] ~ ipa*site, data = mFD_results, permutations = plot_shuffle, by = "margin", method = "euclidean")
adonis2(mFD_results['FRic'] ~ ipa+site, data = mFD_results, permutations = plot_shuffle, by = "margin", method = "euclidean")
#site

pairwise.adonis2(mFD_results['FRic'] ~ site, data = mFD_results, permutations = 999, by = "margin", method = "euclidean")
#A: COR
#B: DOK, EDG, TUR
#C: EDG, FAM, TUR, SHR

#dispersion
FRic_dist <- vegdist(mFD_results['FRic'], method = "euclidean")

FRic.ipa.disp <- betadisper(FRic_dist, mFD_results$ipa, type = c("median"), bias.adjust = FALSE,
                          sqrt.dist = FALSE, add = FALSE)

boxplot(FRic.ipa.disp)
permutest(FRic.ipa.disp, pairwise = TRUE, permutations = plot_shuffle)

FRic.site.disp <- betadisper(FRic_dist, mFD_results$site, type = c("median"), bias.adjust = FALSE,
                           sqrt.dist = FALSE, add = FALSE)

boxplot(FRic.site.disp)
permutest(FRic.site.disp, pairwise = TRUE, permutations = plot_shuffle)
#A: TUR, FAM, EDG
#B: COR
#C: SHR, FAM, TUR
#D: TUR, DOK, EDG

## FEve
adonis2(mFD_results['FEve'] ~ ipa*site, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results['FEve'] ~ ipa+site, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")

FEve_dist <- vegdist(mFD_results['FEve'], method = "euclidean")

FEve.ipa.disp <- betadisper(FEve_dist, mFD_results$ipa, type = c("median"), bias.adjust = FALSE,
                            sqrt.dist = FALSE, add = FALSE)

boxplot(FEve.ipa.disp)
permutest(FEve.ipa.disp, pairwise = TRUE, permutations = plot_shuffle)

FEve.site.disp <- betadisper(FEve_dist, mFD_results$site, type = c("median"), bias.adjust = FALSE,
                             sqrt.dist = FALSE, add = FALSE)

boxplot(FEve.site.disp)
permutest(FEve.site.disp, pairwise = TRUE, permutations = plot_shuffle)
#COR different from EDG

## FDiv
adonis2(mFD_results['FDiv'] ~ ipa*site, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
adonis2(mFD_results['FDiv'] ~ ipa+site, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")


FDiv_dist <- vegdist(mFD_results['FDiv'], method = "euclidean")

FDiv.ipa.disp <- betadisper(FDiv_dist, mFD_results$ipa, type = c("median"), bias.adjust = FALSE,
                            sqrt.dist = FALSE, add = FALSE)

boxplot(FDiv.ipa.disp)
permutest(FDiv.ipa.disp, pairwise = TRUE, permutations = plot_shuffle)

FDiv.site.disp <- betadisper(FDiv_dist, mFD_results$site, type = c("median"), bias.adjust = FALSE,
                             sqrt.dist = FALSE, add = FALSE)

boxplot(FDiv.site.disp)
permutest(FDiv.site.disp, pairwise = TRUE, permutations = plot_shuffle)
#A: FAM, TUR, SHR, DOK, EDG
#B: COR, TUR, SHR, EDG

## FDis
adonis2(mFD_results['FDis'] ~ ipa*site, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#ipa:site (barely)
adonis2(mFD_results['FDis'] ~ ipa+site, data = mFD_results, permutations = 9999, by = "margin", method = "euclidean")
#site

pairwise.adonis2(mFD_results['FDis'] ~ site, data = mFD_results, permutations = 999, by = "margin", method = "euclidean")


FDis_dist <- vegdist(mFD_results['FDis'], method = "euclidean")

FDis.ipa.disp <- betadisper(FDis_dist, mFD_results$ipa, type = c("median"), bias.adjust = FALSE,
                            sqrt.dist = FALSE, add = FALSE)

boxplot(FDis.ipa.disp)
permutest(FDis.ipa.disp, pairwise = TRUE, permutations = plot_shuffle)

FDis.site.disp <- betadisper(FDis_dist, mFD_results$site, type = c("median"), bias.adjust = FALSE,
                             sqrt.dist = FALSE, add = FALSE)

boxplot(FDis.site.disp)
permutest(FDis.site.disp, pairwise = TRUE, permutations = plot_shuffle)
#no differences

#### rda for site ####
mod <- rda(mFD_results[8:11] ~ ipa + site + year, data = mFD_results)
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
  geom_point(aes(color = site, shape = ipa), size = 3) + 
  # stat_ellipse(aes(group = site, color = site), 
  #              linetype = "dashed", show.legend = FALSE) +
  # geom_text_repel(data = points, aes(x = MDS1, y = MDS2, color = site, label = month)) +
  # geom_polygon(data = hulls, aes(x = MDS1, y = MDS2, fill = site), alpha = 0.2) +
  # annotate("text", x = Inf, y = Inf, label = paste("stress = ", S), vjust = 2, hjust = 2) +
  # annotate("text", x = Inf, y = Inf, label = paste("k = ", K), vjust = 2, hjust = 2) +
  theme_minimal() 

#### permanova region ####
whole_plot_shuffle <- how(within = Within(type = "none"),
                          plots = Plots(strata = mFD_results$site, type = "free"),
                          nperm = 999)

#check permutation structure
head(mFD_results[, c("site", "ipa", "year")], 15)
check(mFD_results, control =whole_plot_shuffle)
head(mFD_results[shuffle(nrow(mFD_results), control = whole_plot_shuffle), c("site", "ipa", "year")], 15)

#Species_Richness
adonis2(mFD_results['Species_Richness'] ~ region*veg, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")
adonis2(mFD_results['Species_Richness'] ~ region+veg, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")

SR.region.disp <- betadisper(SR_dist, mFD_results$region, type = c("median"), bias.adjust = FALSE,
                            sqrt.dist = FALSE, add = FALSE)

boxplot(SR.region.disp)
permutest(SR.region.disp, pairwise = TRUE, permutations = whole_plot_shuffle)
#south has lower dispersion

SR.veg.disp <- betadisper(SR_dist, mFD_results$veg, type = c("median"), bias.adjust = FALSE,
                             sqrt.dist = FALSE, add = FALSE)

boxplot(SR.veg.disp)
permutest(SR.veg.disp, pairwise = TRUE, permutations = plot_shuffle)
#higher dispersion when present

#FRic
adonis2(mFD_results['FRic'] ~ region*veg, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")
adonis2(mFD_results['FRic'] ~ region+veg, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")

FRic.region.disp <- betadisper(FRic_dist, mFD_results$region, type = c("median"), bias.adjust = FALSE,
                             sqrt.dist = FALSE, add = FALSE)

boxplot(FRic.region.disp)
permutest(FRic.region.disp, pairwise = TRUE, permutations = whole_plot_shuffle)
#south has lower dispersion

FRic.veg.disp <- betadisper(FRic_dist, mFD_results$veg, type = c("median"), bias.adjust = FALSE,
                          sqrt.dist = FALSE, add = FALSE)

boxplot(FRic.veg.disp)
permutest(FRic.veg.disp, pairwise = TRUE, permutations = plot_shuffle)
#higher dispersion when present

#FEve
adonis2(mFD_results['FEve'] ~ region*veg, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")
adonis2(mFD_results['FEve'] ~ region+veg, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")

FEve.region.disp <- betadisper(FEve_dist, mFD_results$region, type = c("median"), bias.adjust = FALSE,
                               sqrt.dist = FALSE, add = FALSE)

boxplot(FEve.region.disp)
permutest(FEve.region.disp, pairwise = TRUE, permutations = whole_plot_shuffle)

FEve.veg.disp <- betadisper(FEve_dist, mFD_results$veg, type = c("median"), bias.adjust = FALSE,
                            sqrt.dist = FALSE, add = FALSE)

boxplot(FEve.veg.disp)
permutest(FEve.veg.disp, pairwise = TRUE, permutations = plot_shuffle)
#lower dispersion when present

#FDiv
adonis2(mFD_results['FDiv'] ~ region*veg, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")
adonis2(mFD_results['FDiv'] ~ region+veg, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")

FDiv.region.disp <- betadisper(FDiv_dist, mFD_results$region, type = c("median"), bias.adjust = FALSE,
                               sqrt.dist = FALSE, add = FALSE)

boxplot(FDiv.region.disp)
permutest(FDiv.region.disp, pairwise = TRUE, permutations = whole_plot_shuffle)

FDiv.veg.disp <- betadisper(FDiv_dist, mFD_results$veg, type = c("median"), bias.adjust = FALSE,
                            sqrt.dist = FALSE, add = FALSE)

boxplot(FDiv.veg.disp)
permutest(FDiv.veg.disp, pairwise = TRUE, permutations = plot_shuffle)
#lower dispersion when present

#FDis
adonis2(mFD_results['FDis'] ~ region*veg, data = mFD_results, permutations =whole_plot_shuffle, by = "margin", method = "euclidean")
adonis2(mFD_results['FDis'] ~ region+veg, data = mFD_results, permutations = whole_plot_shuffle, by = "margin", method = "euclidean")

FDis.region.disp <- betadisper(FDis_dist, mFD_results$region, type = c("median"), bias.adjust = FALSE,
                               sqrt.dist = FALSE, add = FALSE)

boxplot(FDis.region.disp)
permutest(FDis.region.disp, pairwise = TRUE, permutations = whole_plot_shuffle)

FDis.veg.disp <- betadisper(FDis_dist, mFD_results$veg, type = c("median"), bias.adjust = FALSE,
                            sqrt.dist = FALSE, add = FALSE)

boxplot(FDis.veg.disp)
permutest(FDis.veg.disp, pairwise = TRUE, permutations = plot_shuffle)

## Annual test
annual_shuffle <- how(within = Within(type = "free"),
                      plots = Plots(strata = mFD_results$shoreline, type = "none"),
                      nperm = 9999)

#use rda because it allows you to condition on shoreline, which removes the effects
#of variations between shorelines

mod <- rda(mFD_results['Species_Richness'] ~ year + Condition(shoreline), data = mFD_results)
anova(mod, permutations = annual_shuffle, by = "margin")

mod <- rda(mFD_results['FRic'] ~ year + Condition(shoreline), data = mFD_results)
anova(mod, permutations = annual_shuffle, by = "margin")

mod <- rda(mFD_results['FEve'] ~ year + Condition(shoreline), data = mFD_results)
anova(mod, permutations = annual_shuffle, by = "margin")

mod <- rda(mFD_results['FDiv'] ~ year + Condition(shoreline), data = mFD_results)
anova(mod, permutations = annual_shuffle, by = "margin")

mod <- rda(mFD_results['FDis'] ~ year + Condition(shoreline), data = mFD_results)
anova(mod, permutations = annual_shuffle, by = "margin")

#### rda for region ####

mod1 <- rda(mFD_results[8:11] ~ ipa + site, data = mFD_results)
mod2 <- rda(mFD_results[8:11] ~ region + veg, data = mFD_results)

rda_scores <- scores(mod1)
sites_scores <- as.data.frame(rda_scores[[1]])
biplot_scores <- as.data.frame(rda_scores[[2]])

# RDA_ordiplot <- gg_ordiplot(ord = biplot_scores, #for some reason the scale gets weird if you don't specify this
#                             groups = mFD_results$site,
#                             ellipse = TRUE,
#                             hull = FALSE,
#                             spiders = FALSE)

points <- biplot_scores %>% 
  cbind(mFD_results %>% select(site:year))

ggplot(data = points, aes(x = RDA1, y = RDA2)) +
  geom_text(size = 3, aes(label = paste0(ipa, year), color = site)) + 
  # stat_ellipse(aes(group = region, color = site),
  #              linetype = "dashed", show.legend = TRUE, level = .90) +
  theme_minimal() 

alt_mFD_results <- mFD_results %>% 
  filter(!Species_Richness < 4)

mod3 <- rda(alt_mFD_results[8:11] ~ ipa + site, data = alt_mFD_results)

rda_scores <- scores(mod3)
sites_scores <- as.data.frame(rda_scores[[1]])
biplot_scores <- as.data.frame(rda_scores[[2]])

# RDA_ordiplot <- gg_ordiplot(ord = biplot_scores, 
#                             groups = alt_mFD_results$ipa,
#                             ellipse = FALSE,
#                             hull = TRUE,
#                             spiders = FALSE)



points <- biplot_scores %>% 
  cbind(mFD_results %>% select(site:year))

ggplot(data = points, aes(x = RDA1, y = RDA2)) +
  geom_point(size = 3, aes(shape = ipa, fill = site), color = "black") + 
  scale_shape_manual(values = c(21,24,22)) +
  stat_ellipse(aes(group = site, color = site),
               linetype = "dashed", show.legend = FALSE, level = .90) +
  labs(fill = "Site", shape = "Shoreline Condition") +
  scale_fill_manual(values = site_colors) +
  scale_color_manual(values = site_colors) +
  guides(fill=guide_legend(override.aes=list(color=c(site_colors)))) +
  # geom_text_repel(data = points, aes(x = RDA1, y = RDA2, color = site, label = year)) +
  theme_minimal()

# veg + region results are likely not driven by the samples with FD = 0, they're
# driven by the fact that dock and edg are similar and they're both in south sound with no eelgrass.



