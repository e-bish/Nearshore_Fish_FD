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
                     "lightgoldenrod",
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
# with(mFD_results, cor.test(Species_Richness, FEve))
# with(mFD_results, cor.test(Species_Richness, FDiv))
# with(mFD_results, cor.test(Species_Richness, FRic))
# with(mFD_results, cor.test(Species_Richness, FDis))


by_site <- mFD_results %>% 
  group_by(site, year) %>% 
  summarize_at(vars(Species_Richness:FDis), mean)

# with(by_site, cor.test(Species_Richness, FEve))
# with(by_site, cor.test(Species_Richness, FDiv))
# with(by_site, cor.test(Species_Richness, FRic))
# with(by_site, cor.test(Species_Richness, FDis))


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


## visualize differences 
#### figure S2, FD by ipa ####
ipa_FD <- mFD_results_long %>% 
  ggplot(aes(x = ipa, y = value)) +
  geom_boxplot() +
  geom_point(alpha = 0.4, size = 2, show.legend = FALSE) +
  theme_classic() +
  labs(x = "Condition category", y = "Value", color = "Site") + 
  theme(strip.background = element_rect(fill = NA, colour = NA),
        strip.text.x = element_text(" ")) +
  facet_wrap(~ metric, scales = "free_y", ncol = 1, nrow = 5) 

#### FD by site ####
site_FD <- mFD_results_long %>% 
  ggplot(aes(x = site, 
             # fill = site,
             y = value)) +
  geom_boxplot() +
  geom_point(alpha = 0.4, size = 2, show.legend = FALSE) +
  theme_classic() +
  labs(x = "Site", 
       # fill = "Site"
       y = " ") + 
  # scale_fill_manual(values = site_colors) +
  theme(strip.background = element_rect(fill = NA, colour = NA), 
        strip.text = element_text(face = "bold", size = 12)) +
  facet_wrap(~factor(metric, labels = c("Species Richness",
                                        "Functional Richness",
                                        "Functional Evenness",
                                        "Functional Divergence",
                                        "Functional Dispersion"
    
  )), scales = "free_y", ncol = 1, nrow = 5)

#### Fig 4 ####

plot_index_sc <- function(index_name, metric_name) {
  ggplot(mFD_results, aes(y = .data[[index_name]], x = ipa)) +
    geom_boxplot() +
    geom_point(alpha = 0.5) +
    theme_classic() +
    labs(y = metric_name, x = " ") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
          # plot.margin = unit(c(0, 0.5, -5, 0.5), "mm"))
}

plot_index_yr <- function(index) {
  ggplot(mFD_results, aes(y = .data[[index]], x = year)) +
    geom_boxplot() +
    geom_point(alpha = 0.5) +
    theme_classic() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
          # plot.margin = unit(c(0, 0.5, -5, 0.5), "mm")) + 
    labs(y = " ", x = " ") 
}

metric_n_names <- c("Species\nRichness",
                  "Functional\nRichness",
                  "Functional\nEvenness",
                  "Functional\nDivergence",
                  "Functional\nDispersion")

index_names <- c("Species_Richness", "FRic", "FEve", "FDiv", "FDis")

index_plot_sc<- map2(index_names, metric_n_names, plot_index_sc)
index_plot_yr<- map(index_names, plot_index_yr)

index_plot_sc5 <- mFD_results %>%
  ggplot(aes(x = ipa, y = FDis)) +
  geom_boxplot() +
  geom_point(alpha = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x = "Condition\nCategory", y = "Functional\nDispersion") 

index_plot_yr5 <- mFD_results %>%
  ggplot(aes(x = year, y = FDis)) +
  geom_boxplot() +
  geom_point(alpha = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(x = "Year", y = " ") 

(index_plot_sc[[1]] + index_plot_yr[[1]]) /
  (index_plot_sc[[2]] + index_plot_yr[[2]]) /
  (index_plot_sc[[3]] + index_plot_yr[[3]]) /
  (index_plot_sc[[4]] + index_plot_yr[[4]]) / 
  (index_plot_sc5 + index_plot_yr5)

ggsave(here("figures", "Fig_4.png"), 
       width = 84, height = 175, units = "mm", dpi = 300)

# plot_index <- function(metric_name, index) {
#   mFD_results_long %>% 
#     ggplot(aes(x = .data[[index]], y = value)) +
#     geom_boxplot() +
#     geom_point(alpha = 0.5) +
#     theme_classic() +
#     labs(y = metric_name, x = " ") +
#     facet_wrap(~metric, scales = "free_y", ncol = 1) +
#     theme(
#       # axis.title.x=element_blank(),
#       # axis.text.x=element_blank(),
#       # axis.ticks.x=element_blank(),
#       strip.background = element_blank(),
#       strip.text.x = element_blank())
# }

# yr_plots <- map2(metric_names, "ipa", plot_index_sc)

library(cowplot)

site_colors <- rev(c("#8c510a","#d8b365", 
                     "lightgoldenrod",
                     "lightblue",
                     "#5ab4ac", "#01665e"))

metric_names <- c("Species Richness",
                  "Functional Richness",
                  "Functional Evenness",
                  "Functional Divergence",
                  "Functional Dispersion")

site_plots <- mFD_results_long %>%
  mutate(metric_long = case_when(metric == "FRic" ~ "Functional Richness",
                                metric == "FEve" ~ "Functional Evenness",
                                metric == "FDiv" ~ "Functional Divergence",
                                metric == "FDis" ~ "Functional Dispersion",
         TRUE ~ metric)) %>% 
  mutate(metric_long = factor(metric_long, levels = metric_names)) %>% 
  ggplot(aes(x = site, y = value, fill = site)) +
  geom_boxplot() +
  geom_point(alpha = 0.5, show.legend = FALSE) +
  theme_classic() +
  scale_fill_manual(values = site_colors) +
  labs(x = "Site", fill = "Site") +
  facet_wrap(~metric_long, scales = "free_y", 
             strip.position = 'left') +
  ylab(NULL) +
  theme(strip.background = element_blank(),
        strip.placement='outside',
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

sp1 <- site_plots + theme(legend.position = "none")
sp_legend <- get_legend(site_plots)

ggdraw(sp1) +
  draw_grob(sp_legend, x = 0.75, y = 0.1, width = 0.2, height = 0.3)

ggsave(here("figures", "Fig_4.5.png"), 
       width = 174, height = 150, units = "mm", dpi = 300)

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

pairwise.adonis2(mFD_results['FRic'] ~ site, data = mFD_results,
                 permutations = 999, method = "euclidean")
#COR higher than any other site

pairwise.adonis2(mFD_results['FRic'] ~ year, data = mFD_results,
                 permutations = 999, method = "euclidean")
#2018 v 2019 & 2022

pairwise.adonis2(mFD_results['FEve'] ~ site, data = mFD_results,
                 permutations = 999, method = "euclidean")
# DOK vs.FAM, COR, TUR, SHR

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


