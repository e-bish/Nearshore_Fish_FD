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
                                 deviation_weighting = c("absolute", "squared"),
                                 fdist_scaling = TRUE,
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

trait_plot <- traits_correlations$"tr_faxes_plot"

trait_plot + scale_x_discrete(guide = guide_axis(angle = 90)) 

#calculate FD Indices
alpha_fd_indices <- alpha.fd.multidim(sp_faxes_coord = trait_space[ , c("PC1", "PC2", "PC3")],
                                      asb_sp_w = data.matrix(fish.list$abund),
                                      ind_vect = c("fric", "feve", "fdiv", "fdis"),
                                      scaling = TRUE,
                                      check_input = TRUE,
                                      details_returned = TRUE)

#the mFD package uses ape::pcoa() which automatically removes negative eigenvalues rather than applying a correction

mFD_values <- alpha_fd_indices$"functional_diversity_indices"

metrics <- c("Species_Richness", "FRic", "FEve", "FDiv", "FDis")
metrics_clean <- c("Species Richness", "FRic", "FEve", "FDiv", "FDis")

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
  ggplot(aes(x = site, y = value, fill = site)) +
  geom_boxplot() +
  geom_point(show.legend = FALSE) +
  theme_classic() +
  facet_wrap(~factor(metric, levels = c("Species_Richness", "FDis", "FRic", "FEve", "FDiv"),
                     labels = c("Species Richness", "FDis", "FRic", "FEve", "FDiv")),
             scales = "free_y") + 
  labs(x = "Site", y = "Value") + 
  scale_fill_manual(values = site_colors) +
  theme(strip.background = element_rect(fill = NA, colour = NA))


#### figure S2, FD by ipa ####
mFD_results %>% 
  pivot_longer(!c(site, ipa, year, shoreline, region, veg), names_to = "metric", values_to = "value") %>% 
  ggplot(aes(x = ipa, y = value)) +
  geom_boxplot() +
  geom_point(aes(color = site)) +
  theme_classic() +
  facet_wrap(~factor(metric, levels = c("Species_Richness", "FDis", "FRic", "FEve", "FDiv"),
                     labels = c("Species Richness", "FDis", "FRic", "FEve", "FDiv")),
             scales = "free_y") + 
  labs(x = "Condition category", y = "Value") + 
  scale_color_manual(values = site_colors) +
  theme(strip.background = element_rect(fill = NA, colour = NA))
# 
# ggsave(here("figures", "figure_S2.png"), 
#        width = 8, height = 6, dpi = 300) 

#### FD by region ####
mFD_results %>% 
  pivot_longer(!c(site, ipa, year, shoreline, region, veg), names_to = "metric", values_to = "value") %>% 
  ggplot(aes(x = region, y = value)) +
  geom_boxplot() +
  geom_point(aes(color = site)) +
  theme_classic() +
  facet_wrap(~factor(metric, levels = c("Species_Richness", "FDis", "FRic", "FEve", "FDiv"),
                     labels = c("Species Richness", "FDis", "FRic", "FEve", "FDiv")),
             scales = "free_y") + 
  scale_color_manual(values = site_colors) +
  labs(x = "Region", y = "Value") + 
  theme(strip.background = element_rect(fill = NA, colour = NA))

#### FD by veg ####
mFD_results %>% 
  pivot_longer(!c(site, ipa, year, shoreline, region, veg), names_to = "metric", values_to = "value") %>% 
  ggplot(aes(x = veg, y = value)) +
  geom_boxplot() +
  geom_point(aes(color = site)) +
  theme_classic() +
  facet_wrap(~factor(metric, levels = c("Species_Richness", "FDis", "FRic", "FEve", "FDiv"),
                     labels = c("Species Richness", "FDis", "FRic", "FEve", "FDiv")),
             scales = "free_y") + 
  scale_color_manual(values = site_colors) +
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

## shuffle shorelines to assess site and shoreline condition variables ##
plot_shuffle <- how(within = Within(type = "series"),
                          plots = Plots(strata = mFD_results$shoreline, type = "free"),
                          nperm = 9999)

#check permutation structure
head(mFD_results[, c("site", "ipa", "year")], 15)
check(mFD_results, control =plot_shuffle)
head(mFD_results[shuffle(nrow(mFD_results), control = plot_shuffle), c("site", "ipa", "year")], 15)

#test interactions first 
calc_int.plot_permanova <- function(metric) {
  result <- adonis2(mFD_results[metric] ~ ipa*site*year - ipa:site:year, 
                    data = mFD_results, 
                    permutations = plot_shuffle, 
                    by = "margin", 
                    method = "euclidean")
  
  colnames(result) <- paste(colnames(result), metric, sep = '_')
  
  output <- as.data.frame(result)
  
  return(output)
}


int.plot_permanova_list <- map(metrics, calc_int.plot_permanova)
int.plot_permanova_df <- list_cbind(int.plot_permanova_list, name_repair = "minimal")
int.plot_permanova_df <- int.plot_permanova_df %>% 
  rownames_to_column(var = "X")

write_csv(int.plot_permanova_df, here("data", "int.plot_permanova_df.csv"))

#test single factors 
calc_plot_permanova <- function(metric) {
  result <- adonis2(mFD_results[metric] ~ ipa+site+year, 
                    data = mFD_results, 
                    permutations = plot_shuffle, 
                    by = "margin", 
                    method = "euclidean")
  
  colnames(result) <- paste(colnames(result), metric, sep = '_')
  
  output <- as.data.frame(result) 
  
  return(output)
}

plot_permanova_list <- map(metrics, calc_plot_permanova)
plot_permanova_df <- list_cbind(plot_permanova_list, name_repair = "minimal") 
plot_permanova_df <- plot_permanova_df %>% 
  rownames_to_column(var = "X")

write_csv(plot_permanova_df, here("data", "plot_permanova_df.csv"))

## follow up pairwise tests 

pairwise.adonis2(mFD_results['Species_Richness'] ~ site, data = mFD_results, permutations = 999, by = "margin", method = "euclidean")
#A: COR
#B: DOK, FAM, SHR, TUR
#C: DOK, EDG, FAM, SHR

pairwise.adonis2(mFD_results['FRic'] ~ site, data = mFD_results, permutations = 999, by = "margin", method = "euclidean")
#A: COR
#B: DOK, EDG, TUR
#C: EDG, FAM, TUR, SHR

# SR dispersion
library(patchwork)

compare_dispersion <- function(metric, scale) {
  
  dist <- vegdist(mFD_results[[metric]], method = "euclidean")
  
  disp <- betadisper(dist, mFD_results[[scale]], type = c("median"), bias.adjust = FALSE,
                     sqrt.dist = FALSE, add = FALSE)
  
  output <- as.data.frame(disp$distances)
  
  names(output) <- metric
  
  return(output)
  
}

ipa_dispersions <- mFD_results %>% 
  select(site:year) %>% 
  cbind(map2_dfc(metrics, "ipa", compare_dispersion)) %>% 
  pivot_longer(cols = !site:year, names_to = "metric")

ggplot(data = ipa_dispersions) +
  geom_boxplot(aes(x = ipa, y = value)) +
  theme_classic() + 
  facet_wrap(~metric)

site_dispersions <- mFD_results %>% 
  select(site:year) %>% 
  cbind(map2_dfc(metrics, "site", compare_dispersion)) %>% 
  pivot_longer(cols = !site:year, names_to = "metric") %>% 
  mutate(metric = factor(metric, levels = metrics))

site_full <- ggplot(data = site_dispersions) +
  geom_boxplot(aes(x = site, y = value)) +
  theme_classic() + 
  labs(y = "Dispersion", title = "Full") +
  facet_wrap(~metric, ncol = 1, scales = "free")


compare_sub_dispersion <- function(metric, scale) {
  
  dist <- vegdist(mFD_sub[[metric]], method = "euclidean")
  
  disp <- betadisper(dist, mFD_sub[[scale]], type = c("median"), bias.adjust = FALSE,
                     sqrt.dist = FALSE, add = FALSE)
  
  output <- as.data.frame(disp$distances)
  
  names(output) <- metric
  
  return(output)
  
}

site_sub_dispersions <- mFD_sub %>% 
  select(site:year) %>% 
  cbind(map2_dfc(metrics, "site", compare_sub_dispersion)) %>% 
  pivot_longer(cols = !site:year, names_to = "metric") %>% 
  mutate(metric = factor(metric, levels = metrics))

site_sub <- ggplot(data = site_sub_dispersions) +
  geom_boxplot(aes(x = site, y = value)) +
  theme_classic() + 
  labs(y = "Dispersion", title = "Sub") +
  facet_wrap(~metric, ncol = 1, scales = "free")

site_full + site_sub + plot_layout(axis_titles = "collect")

compare_disp_pval <- function(metric, scale, shuffle) {
  
  dist <- vegdist(mFD_results[[metric]], method = "euclidean")
  
  disp <- betadisper(dist, mFD_results[[scale]], type = c("median"), bias.adjust = FALSE,
                     sqrt.dist = FALSE, add = FALSE)
  
  permu <- permutest(disp, pairwise = TRUE, permutations = shuffle)
  
  pvals <- as.data.frame(permu[["pairwise"]]["permuted"]) %>% 
    rownames_to_column(var = "pair")
  
  names(pvals)[2] <- metric
  
  return(pvals)
  
}

ipa_arg_list <- list(metrics, rep("ipa", times = 5), rep(list(plot_shuffle), times = 5))

ipa_disp_pvals <- pmap_dfc(ipa_arg_list, compare_disp_pval) %>% 
  rename(pairs = "pair...1") %>% 
  select(!contains("..."))


site_arg_list <- list(metrics, rep("site", times = 5), rep(list(plot_shuffle), times = 5))

site_disp_pvals <- pmap_dfc(site_arg_list, compare_disp_pval) %>% 
  rename(pairs = "pair...1") %>% 
  select(!contains("..."))

#### old code ####
SR_dist <- vegdist(mFD_results['Species_Richness'], method = "euclidean")

SR.ipa.disp <- betadisper(SR_dist, mFD_results$ipa, type = c("median"), bias.adjust = FALSE,
                          sqrt.dist = FALSE, add = FALSE)

boxplot(SR.ipa.disp)
test <- permutest(SR.ipa.disp, pairwise = TRUE, permutations = plot_shuffle)

SR.site.disp <- betadisper(SR_dist, mFD_results$site, type = c("median"), bias.adjust = FALSE,
                          sqrt.dist = FALSE, add = FALSE)

boxplot(SR.site.disp)
permutest(SR.site.disp, pairwise = TRUE, permutations = plot_shuffle)
#A: TUR, FAM, SHR, EDG
#B: COR, TUR
#C: FAM, DOK, SHR, EDG

## FRic dispersion
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

### test without sites that have zeros ###
mFD_sub <- mFD_results %>% 
  filter(!Species_Richness < 4)

sub.plot_shuffle <- how(within = Within(type = "series"),
                    plots = Plots(strata = mFD_sub$shoreline, type = "free"),
                    nperm = 9999)

FRic.sub_dist <- vegdist(mFD_sub['FRic'], method = "euclidean")

FRic.sub.disp <- betadisper(FRic.sub_dist, mFD_sub$site, type = c("median"), bias.adjust = FALSE,
                             sqrt.dist = FALSE, add = FALSE)

boxplot(FRic.sub.disp)
#can't test with permutations because the design is no longer balanced, but it looks like patterns are largely the same

## FEve dispersion
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

## FDiv dispersion
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

## FDis dispersion
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

mod1 <- rda(mFD_results['FDiv'] ~ ipa + site + year, data = mFD_results)
mod2 <- rda(mFD_sub['FDiv'] ~ ipa + site + year, data = mFD_sub)
plot(mod1)
plot(mod2)

rda_scores1 <- scores(mod1)
sites_scores1 <- as.data.frame(rda_scores1[[1]])
biplot_scores1 <- as.data.frame(rda_scores1[[2]])

rda_scores2 <- scores(mod2)
sites_scores2 <- as.data.frame(rda_scores2[[1]])
biplot_scores2 <- as.data.frame(rda_scores2[[2]])

gg_ordiplot(ord = biplot_scores1, #for some reason the scale gets weird if you don't specify this
                              groups = mFD_results$site,
                              ellipse = TRUE,
                              hull = FALSE,
                              spiders = FALSE)

gg_ordiplot(ord = biplot_scores2, #for some reason the scale gets weird if you don't specify this
            groups = mFD_sub$site,
            ellipse = TRUE,
            hull = FALSE,
            spiders = FALSE)


points1 <- biplot_scores1 %>% 
  cbind(mFD_results %>% select(site:year))

points2 <- biplot_scores2 %>% 
  cbind(mFD_sub %>% select(site:year))

ggplot(data = points1, aes(x = RDA1, y = PC1)) +
  geom_point(aes(color = site, shape = ipa), size = 3) + 
  stat_ellipse(aes(group = site, color = site),
               linetype = "dashed", show.legend = FALSE) +
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


#test interactions first 
calc_int.whole_permanova <- function(metric) {
  result <- adonis2(mFD_results[metric] ~ region*veg, 
                    data = mFD_results, 
                    permutations = whole_plot_shuffle, 
                    by = "margin", 
                    method = "euclidean")
  
  colnames(result) <- paste(colnames(result), metric, sep = '_')
  
  output <- as.data.frame(result)
  
  return(output)
}


int.whole_permanova_list <- map(metrics, calc_int.whole_permanova)
int.whole_permanova_df <- list_cbind(int.whole_permanova_list, name_repair = "minimal")
int.whole_permanova_df <- int.whole_permanova_df %>% 
  rownames_to_column(var = "X")

write_csv(int.whole_permanova_df, here("data", "int.whole_permanova_df.csv"))

#test single factors 
calc_whole_permanova <- function(metric) {
  result <- adonis2(mFD_results[metric] ~ region + veg, 
                    data = mFD_results, 
                    permutations = whole_plot_shuffle, 
                    by = "margin", 
                    method = "euclidean")
  
  colnames(result) <- paste(colnames(result), metric, sep = '_')
  
  output <- as.data.frame(result) 
  
  return(output)
}

whole_permanova_list <- map(metrics, calc_whole_permanova)
whole_permanova_df <- list_cbind(whole_permanova_list, name_repair = "minimal")
whole_permanova_df <- whole_permanova_df %>% 
  rownames_to_column(var = "X")

write_csv(whole_permanova_df, here("data", "whole_permanova_df.csv"))

#test dispersion
region_dispersions <- mFD_results %>% 
  select(site:year) %>% 
  cbind(map2_dfc(metrics, "region", compare_dispersion)) %>% 
  pivot_longer(cols = !site:year, names_to = "metric")

ggplot(data = region_dispersions) +
  geom_boxplot(aes(x = region, y = value)) +
  theme_classic() + 
  facet_wrap(~metric, scales = "free")

region_arg_list <- list(metrics, rep("region", times = 5), rep(list(whole_plot_shuffle), times = 5))

region_disp_pvals <- pmap_dfc(region_arg_list, compare_disp_pval) %>% 
  rename(pairs = "pair...1") %>% 
  select(!contains("..."))

veg_dispersions <- mFD_results %>% 
  select(site:year) %>% 
  cbind(map2_dfc(metrics, "veg", compare_dispersion)) %>% 
  pivot_longer(cols = !site:year, names_to = "metric")

ggplot(data = veg_dispersions, aes(x = veg, y = value)) +
  geom_boxplot() +
  geom_point(aes(color = site)) +
  theme_classic() + 
  scale_color_manual(values = site_colors) +
  facet_wrap(~metric, scales = "free")

veg_arg_list <- list(metrics, rep("veg", times = 5), rep(list(whole_plot_shuffle), times = 5))

veg_disp_pvals <- pmap_dfc(veg_arg_list, compare_disp_pval) %>% 
  rename(pairs = "pair...1") %>% 
  select(!contains("..."))

## Annual test
annual_shuffle <- how(within = Within(type = "series"),
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



