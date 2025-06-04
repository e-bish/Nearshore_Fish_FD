#Group comparisons after FD analysis

#load libraries
library(tidyverse)
library(here)
library(vegan)
library(patchwork)
library(pairwiseAdonis)
library(multcompView)

#load custom functions
source(here("analysis", "group_comparison_functions.R"))

#load data 
load(here("data","mFD_results.Rdata")) #created in 04_diversity_analysis
metrics <- c("Species_Richness", "FRic", "FEve", "FDiv", "FDis")

#### Permute factors at the shoreline level ####

## shuffle shorelines to assess site and shoreline condition variables ##
plot_shuffle <- how(within = Within(type = "series"),
                    plots = Plots(strata = mFD_results$shoreline, type = "free"),
                    nperm = 9999)

#check permutation structure
# head(mFD_results[, c("site", "ipa", "year")], 15)
# check(mFD_results, control =plot_shuffle)
# head(mFD_results[shuffle(nrow(mFD_results), control = plot_shuffle), c("site", "ipa", "year")], 15)

#test interactions first 
# int.plot_permanova_list <- map(metrics, calc_int.plot_permanova)
# int.plot_permanova_df <- list_cbind(int.plot_permanova_list, name_repair = "minimal")
# int.plot_permanova_df <- int.plot_permanova_df %>% 
#   rownames_to_column(var = "X")

# write_csv(int.plot_permanova_df, here("data", "int.plot_permanova_df.csv"))
int.plot_permanova_df <- read_csv(here("data", "int.plot_permanova_df.csv"))

#test single factors 
# plot_permanova_list <- map(metrics, calc_plot_permanova)
# plot_permanova_df <- list_cbind(plot_permanova_list, name_repair = "minimal") 
# plot_permanova_df <- plot_permanova_df %>% 
#   rownames_to_column(var = "X")

# write_csv(plot_permanova_df, here("data", "plot_permanova_df.csv"))
plot_permanova_df <- read_csv(here("data", "plot_permanova_df.csv"))

## follow up pairwise tests 

pairwise.adonis2(mFD_results['Species_Richness'] ~ site, data = mFD_results, permutations = 999, by = "margin", method = "euclidean")
#A: COR
#B: FAM, SHR, TUR
#C: DOK, EDG, FAM, SHR

pairwise.adonis2(mFD_results['FRic'] ~ site, data = mFD_results, permutations = 999, by = "margin", method = "euclidean")
#A: COR
#B: DOK, EDG
#C: EDG, TUR
#D: SHR, FAM, TUR

mFD_df <- as.data.frame(mFD_results)

pairwise.adonis2(mFD_df['Species_Richness'] ~ year, 
                 data = mFD_df,
                 strata = 'shoreline', 
                 by = "margin", method = "euclidean")

pairwise.adonis2(mFD_df['FRic'] ~ year, 
                 data = mFD_df,
                 strata = 'shoreline', 
                 by = "margin", method = "euclidean")

# dispersion between shoreline conditions
ipa_dispersions <- mFD_results %>% 
  select(site:year) %>% 
  cbind(map2_dfc(metrics, "ipa", compare_dispersion)) %>% 
  pivot_longer(cols = !site:year, names_to = "metric") %>% 
  mutate(metric = factor(metric, levels = metrics))

ggplot(data = ipa_dispersions) +
  geom_boxplot(aes(x = ipa, y = value)) +
  theme_classic() + 
  facet_wrap(~metric, scales = "free")

ipa_arg_list <- list(metrics, rep("ipa", times = 5), rep(list(plot_shuffle), times = 5))

ipa_disp_pvals <- pmap_dfc(ipa_arg_list, compare_disp_pval) %>% 
  rename(pairs = "pair...1") %>% 
  select(!contains("..."))
#no significant differences between shoreline conditions

## site dispersions
site_dispersions <- mFD_results %>% 
  select(site:year) %>% 
  cbind(map2_dfc(metrics, "site", compare_dispersion)) %>% 
  pivot_longer(cols = !site:year, names_to = "metric") %>% 
  mutate(metric = factor(metric, levels = metrics))

site_arg_list <- list(metrics, rep("site", times = 5), rep(list(plot_shuffle), times = 5))

site_disp_pvals <- pmap_dfc(site_arg_list, compare_disp_pval) %>% 
  rename(pairs = "pair...1") %>% 
  dplyr::select(!contains("..."))

site_letters <- map_dfc(metrics, create_letter_df)  %>% 
  rownames_to_column(var = "site") %>% 
  pivot_longer(!site, names_to = "metric", values_to = "letters")

site_signif <- site_dispersions %>% 
  group_by(site, metric) %>% 
  summarize(max = max(value)) %>% 
  full_join(site_letters)

range_df <- site_dispersions %>% 
  group_by(site, metric) %>% 
  summarize(max = max(value)) %>% 
  mutate(ymax = max + 0.35*max) %>% 
  select(!max)

site_dispersions_signif <- full_join(site_dispersions, site_signif) %>% 
  full_join(range_df) %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")),
         metric = factor(metric, levels = metrics)) %>% 
  mutate(letters = ifelse(metric == "FDis", NA, letters))

site_full <- ggplot(data = site_dispersions_signif) +
  geom_boxplot(aes(x = site, y = value)) +
  geom_point(aes(y = ymax, x = site), color = "white", size = 0) +
  geom_text(aes(label = letters, y = max, x = site), vjust = -0.5) +
  theme_classic() + 
  labs(y = "Dispersion", title = "Full") +
  facet_wrap(~metric, ncol = 1, scales = "free")
  
## test for differences in site without zeros
mFD_sub <- mFD_results %>% 
  filter(!Species_Richness < 4)

site_sub_disp <- mFD_sub %>% 
  select(site:year) %>% 
  cbind(map2_dfc(metrics, "site", compare_sub_dispersion)) 
#cannot test with the same permutation structure because the design is no longer balanced

# with(site_sub_disp, (kruskal.test(FRic ~ site)))
# with(site_sub_disp, (kruskal.test(FEve ~ site)))
# with(site_sub_disp, (kruskal.test(FDiv ~ site)))
# with(site_sub_disp, (kruskal.test(FDis ~ site)))
# 
# with(site_sub_disp, (pairwise.wilcox.test(FRic, site)))
# with(site_sub_disp, (pairwise.wilcox.test(FEve, site)))
# with(site_sub_disp, (pairwise.wilcox.test(FDiv, site)))
# with(site_sub_disp, (pairwise.wilcox.test(FDis, site)))

site_sub_dispersions <- site_sub_disp %>% 
  pivot_longer(cols = !site:year, names_to = "metric") %>% 
  mutate(metric = factor(metric, levels = metrics))

site_sub <- ggplot(data = site_sub_dispersions) +
  geom_boxplot(aes(x = site, y = value)) +
  theme_classic() + 
  labs(y = "Dispersion", title = "Sub") +
  facet_wrap(~metric, ncol = 1, scales = "free")

site_full + site_sub + plot_layout(axis_titles = "collect")

#### rda for site ####

mod1 <- rda(mFD_results[8:11] ~ ipa + site + year, data = mFD_results)
mod2 <- rda(mFD_sub[8:11] ~ ipa + site + year, data = mFD_sub)
plot(mod1)
plot(mod2)

rda_scores1 <- scores(mod1)
sites_scores1 <- as.data.frame(rda_scores1$sites)
biplot_scores1 <- as.data.frame(rda_scores1$species)

rda_scores2 <- scores(mod2)
sites_scores2 <- as.data.frame(rda_scores2$sites)
biplot_scores2 <- as.data.frame(rda_scores2$species)

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
  cbind(mFD_results %>% select(site:year))

points2 <- sites_scores2 %>% 
  cbind(mFD_sub %>% select(site:year))

site_colors <- rev(c("#8c510a","#d8b365", 
                     # "#f6e8c1",
                     "lightgoldenrod",
                     # "#c7eae8",
                     "lightblue",
                     "#5ab4ac", "#01665e"))

RDA_full <- ggplot(data = points1, aes(x = RDA1, y = RDA2)) +
  geom_point(aes(color = site, shape = ipa), size = 3) + 
  stat_ellipse(aes(group = site, color = site),
               linetype = "dashed", show.legend = FALSE) +
  geom_segment(data=biplot_scores1,
               aes(x = 0, y = 0, xend=RDA1, yend=RDA2),
               arrow=arrow(length = unit(0.01, "npc")),
               lwd=0.75) +
  geom_text(data=biplot_scores1,
            aes(x=RDA1*0.9,
                y=RDA2*0.9,
                label=metrics[-1]),
            nudge_x = c(-0.2, -0.3, -0.2, -0.2), 
            nudge_y = c(0.1, -0.05, 0, 0),
            size=4) +
  scale_color_manual(values = site_colors) +
  labs(shape = "Shoreline Condition", 
       color = "Site",
       title = "Full") +
  theme_minimal() 


RDA_sub <- ggplot(data = points2, aes(x = RDA1, y = RDA2)) +
  geom_point(aes(color = site, shape = ipa), size = 3) + 
  stat_ellipse(aes(group = site, color = site),
               linetype = "dashed", show.legend = FALSE) +
  geom_segment(data=biplot_scores2,
               aes(x = 0, y = 0, xend=RDA1, yend=RDA2),
               arrow=arrow(length = unit(0.01, "npc")),
               lwd=0.75) +
  geom_text(data=biplot_scores2,
            aes(x=RDA1*0.9,
                y=RDA2*0.9,
                label=metrics[-1]),
            nudge_x = c(-0.15, -0.15, -0.15, -0.15), 
            nudge_y = c(0, -0.13, 0.1, 0),
            size=4) +
  scale_color_manual(values = site_colors) +
  labs(shape = "Shoreline Condition", 
       color = "Site",
       title = "Sub") +
  theme_minimal() 

RDA_full + RDA_sub + plot_layout(guides = 'collect')

# dispersion between years
year_dispersions <- mFD_results %>% 
  select(site:year) %>% 
  cbind(map2_dfc(metrics, "year", compare_dispersion)) %>% 
  pivot_longer(cols = !site:year, names_to = "metric") %>% 
  mutate(metric = factor(metric, levels = metrics))

ggplot(data = year_dispersions) +
  geom_boxplot(aes(x = year, y = value)) +
  theme_classic() + 
  facet_wrap(~metric, scales = "free")

year_arg_list <- list(metrics, rep("year", times = 5), rep(list(plot_shuffle), times = 5))

year_disp_pvals <- pmap_dfc(year_arg_list, compare_disp_pval) %>% 
  rename(pairs = "pair...1") %>% 
  select(!contains("..."))
#no significant differences between years except FDiv dispersion is lower in 2019 and 2022

#### permanova region ####
whole_plot_shuffle <- how(within = Within(type = "none"),
                          plots = Plots(strata = mFD_results$site, type = "free"),
                          nperm = 999)

#check permutation structure
head(mFD_results[, c("site", "ipa", "year")], 15)
check(mFD_results, control =whole_plot_shuffle)
head(mFD_results[shuffle(nrow(mFD_results), control = whole_plot_shuffle), c("site", "ipa", "year")], 15)

#test interactions first 
int.whole_permanova_list <- map(metrics, calc_int.whole_permanova)
int.whole_permanova_df <- list_cbind(int.whole_permanova_list, name_repair = "minimal")
int.whole_permanova_df <- int.whole_permanova_df %>% 
  rownames_to_column(var = "X")

write_csv(int.whole_permanova_df, here("data", "int.whole_permanova_df.csv"))

#test single factors 
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
aov <- anova(mod, permutations = annual_shuffle, by = "margin", model = "reduced")

mod <- rda(mFD_results['FRic'] ~ year + Condition(shoreline), data = mFD_results)
anova(mod, permutations = annual_shuffle, by = "margin", model = "reduced")

mod <- rda(mFD_results['FEve'] ~ year + Condition(shoreline), data = mFD_results)
anova(mod, permutations = annual_shuffle, by = "margin", model = "reduced")

mod <- rda(mFD_results['FDiv'] ~ year + Condition(shoreline), data = mFD_results)
anova(mod, permutations = annual_shuffle, by = "margin", model = "reduced")

mod <- rda(mFD_results['FDis'] ~ year + Condition(shoreline), data = mFD_results)
anova(mod, permutations = annual_shuffle, by = "margin", model = "reduced")

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



