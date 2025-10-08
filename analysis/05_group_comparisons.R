#Group comparisons after FD analysis

#load libraries
library(tidyverse)
library(here)
library(vegan)
library(patchwork)
library(pairwiseAdonis)
library(multcompView)

set.seed(2025)

site_colors <- rev(c("#8c510a","#d8b365", 
                     # "#f6e8c1",
                     "lightgoldenrod",
                     # "#c7eae8",
                     "lightblue",
                     "#5ab4ac", "#01665e"))


#load custom functions
source(here("analysis", "group_comparison_functions.R"))

#load data 
load(here("data","mFD_results.Rdata")) #created in 04_diversity_analysis
metrics <- c("Species_Richness", "FRic", "FEve", "FDiv", "FDis")

mFD_results_long <- mFD_results %>% 
  pivot_longer(!c(site, ipa, year, shoreline, region, veg), names_to = "metric", values_to = "value") %>% 
  mutate(metric = factor(metric, levels = c("Species_Richness", "FRic", "FEve", "FDiv", "FDis"),
         labels = c("Species Richness", "FRic", "FEve", "FDiv", "FDis")))

#### summarize values ####
with(mFD_results, cor.test(Species_Richness, FEve))
with(mFD_results, cor.test(Species_Richness, FDiv))
with(mFD_results, cor.test(Species_Richness, FRic))
with(mFD_results, cor.test(Species_Richness, FDis))


by_site <- mFD_results %>% 
  group_by(site, year) %>% 
  summarize_at(vars(Species_Richness:FDis), mean)

with(by_site, cor.test(Species_Richness, FEve))
with(by_site, cor.test(Species_Richness, FDiv))
with(by_site, cor.test(Species_Richness, FRic))
with(by_site, cor.test(Species_Richness, FDis))


##### summary stats #####

#Species richness 

mFD_results %>% 
  filter(site == "COR") %>% 
  summarize(min = min(Species_Richness), max= max(Species_Richness), avg = mean(Species_Richness), sd = sd(Species_Richness))

mFD_results %>% 
  filter(!site == "COR") %>% 
  summarize(min = min(Species_Richness), max= max(Species_Richness), avg = mean(Species_Richness), sd = sd(Species_Richness))

mFD_results %>% 
  summarize(min = min(FRic), max= max(FRic), avg = mean(FRic), sd = sd(FRic))

mFD_results %>% 
  summarize(min = min(FEve), max= max(FEve), avg = mean(FEve), sd = sd(FEve))

mFD_results %>% 
  summarize(min = min(FDiv), max= max(FDiv), avg = mean(FDiv), sd = sd(FDiv))

mFD_results %>% 
  summarize(min = min(FDis), max= max(FDis), avg = mean(FDis), sd = sd(FDis))


#### visualize differences ####

# plot with means
# mFD_averages <- mFD_results %>% 
#   pivot_longer(!c(shoreline, site, ipa, year, region, veg), names_to = "metric", values_to = "value") %>% 
#   group_by(site, ipa, shoreline, metric) %>% 
#   summarize(mean = mean(value), sd = sd(value), max = max(value), min = min(value)) %>% 
#   filter(!metric %in% c("Species_Richness", "FRic", "FDis"))
# 
# plot_prep <- mFD_results %>%
#   mutate(ipa2 = ifelse(shoreline == "TURN2", "alt", "no_alt")) %>% 
#   group_by(site, ipa, ipa2) %>% 
#   mutate(id = factor(cur_group_id())) %>%
#   ungroup() %>% 
#   pivot_longer(!c(shoreline, site, ipa, ipa2, year, region, veg, id), 
#                names_to = "metric", values_to = "value") %>%
#   group_by(site, ipa, id, metric) %>% 
#   summarize(avg = mean(value), min = min(value), max = max(value)) %>% 
#   ungroup() %>% 
#   filter(!metric %in% c("Species_Richness", "FRic", "FDis"))
# 
# site_colors <- rev(c("#8c510a","#d8b365", 
#                      # "#f6e8c1",
#                      "lightgoldenrod",
#                      # "#c7eae8",
#                      "lightblue",
#                      "#5ab4ac", "#01665e"))
# 
# plot_prep %>% 
#   ggplot(aes(x = id, y = avg, color = site, shape = ipa), show.legend = TRUE) +
#   geom_point(size = 3) +
#   geom_linerange(aes(ymin=min,ymax=max),linetype=1, show.legend = FALSE)+
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 60, vjust = 0.75, hjust=1),
#         legend.box = "horizontal", 
#         legend.title = element_text(hjust = 0.5),
#         strip.placement = "outside",
#         strip.text.x = element_text(size = 10),
#         strip.background = element_rect(fill = NA,
#                                         colour = NA)) +
#   scale_color_manual(values = site_colors) +
#   labs(x = "Shoreline ID", y = "Value", color = "Site", shape = "Shoreline\ntype") +
#   facet_wrap(~metric, 
#              scales = "free_y", axes = "all_x", axis.labels = "margins")
# 
# layout <- '
# AB
# CC
# '
# 
# free(p1) + guide_area() + p2 + 
#   plot_layout(guides = 'collect', 
#               design = layout,
#               heights = c(1.1,2)) &
#   theme(legend.direction = 'vertical',
#         legend.box = 'horizontal')
# # 
# ggsave(here("figures", "figure_3.png"), 
#        width = 8, height = 6, dpi = 300) 


#### figure S2, FD by ipa ####
ipa_FD <- mFD_results_long %>% 
  ggplot(aes(x = ipa, y = value)) +
  geom_boxplot() +
  geom_point(alpha = 0.4, size = 2, show.legend = FALSE) +
  theme_classic() +
  labs(x = "Condition category", y = "Value", color = "Site") + 
  scale_color_manual(values = site_colors) +
  theme(strip.background = element_rect(fill = NA, colour = NA),
        strip.text.x = element_text(" ")) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1, nrow = 5) 
# 
# ggsave(here("figures", "figure_S2.png"), 
#        width = 8, height = 6, dpi = 300) 

#### FD by site ####
site_FD <- mFD_results_long %>% 
  ggplot(aes(x = site, y = value, fill = site)) +
  geom_boxplot() +
  geom_point(alpha = 0.4, size = 2, show.legend = FALSE) +
  theme_classic() +
  labs(x = "Site", y = " ", fill = "Site") + 
  scale_fill_manual(values = site_colors) +
  theme(strip.background = element_rect(fill = NA, colour = NA), 
        strip.text = element_text(face = "bold", size = 12)) +
  facet_wrap(~factor(metric, labels = c("Species Richness",
                                        "Functional Richness",
                                        "Functional Evenness",
                                        "Functional Divergence",
                                        "Functional Dispersion"
    
  )), scales = "free_y", ncol = 1, nrow = 5)

#### FD by year ####
year_FD <- mFD_results_long %>% 
  ggplot(aes(x = year, y = value)) +
  geom_boxplot() +
  geom_point(alpha = 0.4, size = 2, show.legend = FALSE) +
  theme_classic() +
  scale_color_manual(values = site_colors) +
  labs(x = "Year", y = " ") + 
  theme(strip.background = element_rect(fill = NA, colour = NA),
        strip.text.x = element_text(" ")) +
  facet_wrap(~metric, scales = "free_y", ncol = 1) 

ipa_FD + site_FD + year_FD + plot_layout(guides = "collect")


#### FD by region ####
# mFD_results %>%
#   pivot_longer(!c(site, ipa, year, shoreline, region, veg), names_to = "metric", values_to = "value") %>%
#   ggplot(aes(x = region, y = value)) +
#   geom_boxplot() +
#   geom_point(aes(color = site)) +
#   theme_classic() +
#   facet_wrap(~metric,
#              scales = "free_y") +
#   scale_color_manual(values = site_colors) +
#   labs(x = "Region", y = "Value") +
#   theme(strip.background = element_rect(fill = NA, colour = NA))
#consider the effect of COR if revisiting this

# # # #### FD by veg ####
# mFD_results %>%
#   pivot_longer(!c(site, ipa, year, shoreline, region, veg), names_to = "metric", values_to = "value") %>%
#   ggplot(aes(x = veg, y = value)) +
#   geom_boxplot() +
#   geom_point(aes(color = site)) +
#   theme_classic() +
#   facet_wrap(~metric,
#              scales = "free_y") +
#   scale_color_manual(values = site_colors) +
#   labs(x = "Eelgrass", y = "Value") +
#   theme(strip.background = element_rect(fill = NA, colour = NA))
#consider the effect of COR if revisiting this


### try Kruskal- Wallis test
# with(mFD_results, kruskal.test(Species_Richness,ipa))
# with(mFD_results, kruskal.test(SESFRic,ipa))
# with(mFD_results, kruskal.test(FEve,ipa))
# with(mFD_results, kruskal.test(FDiv,ipa))
# with(mFD_results, kruskal.test(SESFDis,ipa))
# 
# with(mFD_results, dunn.test(FDiv, g=ipa, method= "bonferroni", kw=FALSE))
# 
# hist(log(mFD_results$Species_Richness))
# hist(log(mFD_results$SESFRic))
# hist(log(mFD_results$FEve))
# hist(log(mFD_results$FDiv))
# hist(sqrt(mFD_results$SESFDis))


# with(mFD_results, kruskal.test(Species_Richness,site))
# with(mFD_results, kruskal.test(SESFRic,site))
# with(mFD_results, kruskal.test(FEve,site))
# with(mFD_results, kruskal.test(FDiv,site))
# with(mFD_results, kruskal.test(SESFDis,site))
# 
# library(dunn.test)
# 
# with(mFD_results, dunn.test(SESFRic, g=site, method= "bonferroni", kw=FALSE))
# with(mFD_results, dunn.test(FEve, g=site, method= "bonferroni", kw=FALSE))


#### Permute factors at the shoreline level ####

## shuffle shorelines to assess site and shoreline condition variables ##
plot_shuffle <- how(within = Within(type = "free"),
                    plots = Plots(strata = mFD_results$shoreline, type = "free"),
                    nperm = 9999)

#check permutation structure
# head(mFD_results[, c("site", "ipa", "year")], 15)
# check(mFD_results, control =plot_shuffle)
# head(mFD_results[shuffle(nrow(mFD_results), control = plot_shuffle), c("site", "ipa", "year")], 15)

#test interactions first 
int.plot_permanova_list <- map(metrics[-1], calc_int.plot_permanova)
int.plot_permanova_df <- list_cbind(int.plot_permanova_list, name_repair = "minimal")
int.plot_permanova_df <- int.plot_permanova_df %>%
  rownames_to_column(var = "X") %>% 
  mutate(across(where(is.numeric), round, 3)) %>% 
  select(!contains("R2"))

# write_csv(int.plot_permanova_df, here("data", "int.plot_permanova_df.csv")) #last saved 6/27/25
int.plot_permanova_df <- read_csv(here("data", "int.plot_permanova_df.csv"))

#test single factors 
plot_permanova_list <- map(metrics[-1], calc_plot_permanova)
plot_permanova_df <- list_cbind(plot_permanova_list, name_repair = "minimal")
plot_permanova_df <- plot_permanova_df %>%
  rownames_to_column(var = "X") %>% 
  mutate(across(where(is.numeric), round, 4)) %>% 
  select(!contains("R2"))

# write_csv(plot_permanova_df, here("data", "plot_permanova_df.csv")) #last saved 6/27/25
plot_permanova_df <- read_csv(here("data", "plot_permanova_df.csv"))

## follow up pairwise tests 

# pairwise.adonis2(mFD_results['Species_Richness'] ~ site, data = mFD_results,
#                  permutations = 999, by = "margin", method = "euclidean")
# #COR v DOK, EDG, FAM, SHR, TUR
# #DOK & EDG v TUR


## pairwise tests

pairwise.adonis2(mFD_results['FRic'] ~ site, data = mFD_results,
                 permutations = 999, method = "euclidean")
#COR higher than any other site

pairwise.adonis2(mFD_results['FRic'] ~ year, data = mFD_results,
                 permutations = 999, method = "euclidean")
#2018 v 2019 & 2022

pairwise.adonis2(mFD_results['FEve'] ~ site, data = mFD_results,
                 permutations = 999, method = "euclidean")
# DOK vs.FAM, COR, TUR, SHR

# pairwise.adonis2(mFD_results['Species_Richness'] ~ year, 
#                  data = mFD_results,
#                  by = "margin", method = "euclidean")
# #2018 v 2022



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

disp_summary <- site_dispersions %>% 
  group_by(site, metric) %>% 
  summarize(median = median(value)) %>% 
  group_by(metric) %>% 
  summarize(min = min(median), max = max(median)) %>% 
  mutate(diff = max/min)

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
         metric = factor(metric, levels = metrics, 
                         labels = c("Species Richness", "Functional Richness", "Functional Evenness", "Functional Divergence", "Functional Dispersion"))) %>%
  mutate(letters = ifelse(metric %in% c("Functional Richness", "Functional Dispersion") , NA, letters)) 

site_dispersions_signif %>% 
  filter(!metric == "Species Richness") %>% 
  ggplot() +
  geom_boxplot(aes(x = site, y = value)) +
  geom_point(aes(y = ymax, x = site), color = "white", size = 0) +
  geom_text(aes(label = letters, y = max, x = site), vjust = -0.5) +
  theme_classic() + 
  labs(y = "Distance to median", x = "Site") +
  facet_wrap(~metric, ncol = 1, scales = "free")

site_dispersions_signif %>% 
  filter(site == "FAM" | site == "EDG") %>% 
  group_by(site, metric) %>% 
  summarize(mean = mean(value), median = median(value))

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


