#Group comparisons after FD analysis

#load libraries
library(tidyverse)
library(here)
library(vegan)
library(patchwork)
library(pairwiseAdonis)


load(here("data","mFD_results.Rdata")) #created in 04_diversity_analysis
metrics <- c("Species_Richness", "FRic", "FEve", "FDiv", "FDis")

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

# write_csv(int.plot_permanova_df, here("data", "int.plot_permanova_df.csv"))
# int.plot_permanova_df <- read_csv(here("data", "int.plot_permanova_df.csv"))

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

# write_csv(plot_permanova_df, here("data", "plot_permanova_df.csv"))
# plot_permanova_df <- read_csv(here("data", "plot_permanova_df.csv"))

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
  dplyr::select(!contains("..."))

SR_disp <- site_disp_pvals %>% 
  pull(Species_Richness, pairs)
FRic_disp <- site_disp_pvals %>% 
  pull(FRic, pairs)
FEve_disp <- site_disp_pvals %>% 
  pull(FEve, pairs)
FDiv_disp <- site_disp_pvals %>% 
  pull(FDiv, pairs)
FDis_disp <- site_disp_pvals %>% 
  pull(FDis, pairs)

SR_disp_letters <- multcompLetters(SR_disp)
FR_disp_letters <- multcompLetters(FRic_disp)
FE_disp_letters <- multcompLetters(FEve_disp)
FV_disp_letters <- multcompLetters(FDiv_disp) 
# FS_disp_letters <- multcompLetters(FDis_disp) no difference

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



