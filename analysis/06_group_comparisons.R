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
metric_names <- c("Species Richness",
                  "Functional Richness",
                  "Functional Evenness",
                  "Functional Divergence",
                  "Functional Dispersion")

index_plot1 <- mFD_results %>% 
  ggplot(aes(x = ipa, y = Species_Richness)) +
  geom_boxplot() +
  geom_point(alpha = 0.5) +
  theme_classic() +
  labs(x = " ", y = "Species\nRichness") + 
  theme(plot.margin = unit(c(0, 0, -5, 0.5), "mm"))

plot_index_2 <- function(index, metric) {
  ggplot(mFD_results, aes(y = .data[[index]], x = site)) +
    geom_boxplot() +
    geom_point(alpha = 0.5) +
    theme_classic() +
    labs(y = " ", x = " ") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust=1),
          plot.margin = unit(c(0, 0.5, -5, 0.5), "mm"))
}

plot_index_3 <- function(index, view) {
  ggplot(mFD_results, aes(y = .data[[index]], x = .data[[view]])) +
    geom_boxplot() +
    geom_point(alpha = 0.5) +
    theme_classic() +
    labs(y = " ", x = " ") +
    theme(plot.margin = unit(c(0, 0.5, -5, 0.5), "mm"))
}

index_plot4 <- mFD_results %>% 
  ggplot(aes(x = ipa, y = FDis)) +
  geom_boxplot() +
  geom_point(alpha = 0.5) +
  theme_classic() +
  labs(x = "Shoreline Condition", y = "FDis") + 
  theme(plot.margin = unit(c(0, 0.5, -5, 0.5), "mm"))

index_plot5 <- mFD_results %>% 
  ggplot(aes(x = site, y = FDis)) +
  geom_boxplot() +
  geom_point(alpha = 0.5) +
  theme_classic() +
  labs(x = "Site", y = " ") +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust=1),
        plot.margin = unit(c(0, 0.5, -5, 0.5), "mm"))

index_plot6 <- mFD_results %>% 
  ggplot(aes(x = year, y = FDis)) +
  geom_boxplot() +
  geom_point(alpha = 0.5) +
  theme_classic() +
  labs(x = "Year", y = " ") + 
  theme(plot.margin = unit(c(0, 0.5, -5, 0.5), "mm"))

plot_index_7 <- function(index, view) {
  ggplot(mFD_results, aes(y = .data[[index]], x = .data[[view]])) +
    geom_boxplot() +
    geom_point(alpha = 0.5) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 25, vjust = 0.9, hjust=1),
          plot.margin = unit(c(0, 0.5, -5, 0.5), "mm")) + 
    labs(y = index, x = " ") 
}

# Create all combinations

index_names <- c("Species_Richness", "FRic", "FEve", "FDiv", "FDis")
views <- c("ipa", "site", "year")

index_grid <- expand_grid(index = index_names, view = views)

# Use map2 to iterate properly
index_plot2 <- map2(index_names, metric_names, plot_index_2)
index_plot3 <- map2(index_grid$index, index_grid$view, plot_index_3)
index_plot7 <- map2(index_grid$index, index_grid$view, plot_index_7)

# make headings 
make_heading <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.55, label = label,
             size = 4, fontface = "bold") +
    annotate("segment", x = 0, xend = 1, y = 0.15, yend = 0.15,
             linewidth = 0.8) +
    expand_limits(x = c(0, 1), y = c(0, 1)) +   
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin = margin(-10, 0, 0, 0))
}

titles <- list(
  make_heading("Species Richness"),
  make_heading("Functional Richness"),
  make_heading("Functional Evenness"),
  make_heading("Functional Divergence"),
  make_heading("Functional Dispersion")
)


# Combine all plots
(titles[[1]] /
(index_plot1 + index_plot2[[1]] + index_plot3[[3]] + plot_layout(ncol = 3)) /
  titles[[2]] /
  (index_plot7[[4]] + index_plot2[[2]] + index_plot3[[6]]+ plot_layout(ncol = 3)) / 
  titles[[3]] / 
  (index_plot7[[7]] + index_plot2[[3]] + index_plot3[[9]]+ plot_layout(ncol = 3)) / 
  titles[[4]] / 
  (index_plot7[[10]] + index_plot2[[4]] + index_plot3[[12]]+ plot_layout(ncol = 3)) / 
  titles[[5]] / 
  (index_plot4 + index_plot5 + index_plot6 + plot_layout(ncol = 3))) + 
  plot_layout(heights = c(0.5, 1, 0.5, 1, 0.5, 1, 0.5, 1, 0.5, 1))

ggsave(here("figures", "Fig_4.png"), 
       width = 174, height = 215, units = "mm", dpi = 300)

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


