#Group comparisons after FD analysis

#load libraries
library(tidyverse)
library(here)
library(vegan)
library(patchwork)
library(pairwiseAdonis)
library(multcompView)

set.seed(2025)

#load custom functions
source(here("analysis", "group_comparison_functions.R"))

#load data 
load(here("data","mFD_results.Rdata")) #created in 04_diversity_analysis
load(here("data", "SES_tab.Rdata")) #created in 05_SES_calc
metrics <- c("Species_Richness", "SESFRic", "FEve", "FDiv", "SESFDis")

mFD_results <- mFD_results %>% 
 inner_join(SES_tab)

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
int.plot_permanova_list <- map(metrics, calc_int.plot_permanova)
int.plot_permanova_df <- list_cbind(int.plot_permanova_list, name_repair = "minimal")
int.plot_permanova_df <- int.plot_permanova_df %>%
  rownames_to_column(var = "X") %>% 
  mutate(across(where(is.numeric), round, 3)) %>% 
  select(!contains("R2"))

# write_csv(int.plot_permanova_df, here("data", "int.plot_permanova_df.csv")) #last saved 6/27/25
int.plot_permanova_df <- read_csv(here("data", "int.plot_permanova_df.csv"))

#test single factors 
plot_permanova_list <- map(metrics, calc_plot_permanova)
plot_permanova_df <- list_cbind(plot_permanova_list, name_repair = "minimal")
plot_permanova_df <- plot_permanova_df %>%
  rownames_to_column(var = "X") %>% 
  mutate(across(where(is.numeric), round, 3)) %>% 
  select(!contains("R2"))

# write_csv(plot_permanova_df, here("data", "plot_permanova_df.csv")) #last saved 6/27/25
plot_permanova_df <- read_csv(here("data", "plot_permanova_df.csv"))

## follow up pairwise tests 

pairwise.adonis2(mFD_results['Species_Richness'] ~ site, data = mFD_results,
                 permutations = 999, by = "margin", method = "euclidean")
#COR v DOK, EDG, FAM, SHR, TUR
#DOK & EDG v TUR


pairwise.adonis2(mFD_results['FEve'] ~ site, data = mFD_results,
                 permutations = 999, by = "margin", method = "euclidean")
# DOK vs.FAM, COR, TUR, SHR

pairwise.adonis2(mFD_results['Species_Richness'] ~ year, 
                 data = mFD_results,
                 by = "margin", method = "euclidean")
#2018 v 2022

pairwise.adonis2(mFD_results['SESFRic'] ~ year, data = mFD_results,
                 permutations = 999, by = "margin", method = "euclidean")
#2018 & 2022 v 2019

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
  mutate(letters = ifelse(metric == "SESFDis" , NA, letters))

site_dispersions_signif %>% 
  filter(!metric == "Species_Richness") %>% 
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

##nmds
nmds <- metaMDS(fish_L_full[4:45], 
                distance="bray", k= 3, trymax=1000, plot = FALSE)

wq_tb_exp <- wq_tb_prep %>% 
  mutate(site = str_sub(shoreline, 1, 3),
         ipa = str_sub(shoreline, 4), .before = shoreline) %>% 
  select(!c(shoreline,min,max)) %>% 
  pivot_wider(names_from = metric, values_from = mean)

load(here("data", "fish_L_full.Rdata"))

fish_L_full <- fish_L_full %>% 
  full_join(wq_tb_exp) %>% 
  mutate(ipa = ifelse(ipa == "Natural2", "Natural", ipa))
  
#   mutate(region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"), 
#           veg = ifelse(site %in% c("TUR", "COR", "SHR"), "present", "absent"), .after = year)

env.test <- envfit(nmds, fish_L_full[c(46,47,48,49)], permutations = 9999)
env.test 

wq_vecs <- as.data.frame(scores(env.test, "vectors")) * ordiArrowMul(env.test)

points <- data.frame(nmds$points) %>% 
  cbind(mFD_results %>% select(site:year)) 
  # anti_join(add_small_samples_back)

hulls <- points %>%
  group_by(site) %>% 
  slice(chull(MDS1,MDS2))

ggplot() +
  geom_point(data = points, aes(x = MDS1, y = MDS2, color = site, shape = ipa), size = 3) + 
  # geom_text_repel(data = points, aes(x = MDS1, y = MDS2, color = site, label = year), max.overlaps = 20) +
  geom_polygon(data = hulls, aes(x = MDS1, y = MDS2, fill = site), alpha = 0.2) +
  # annotate("text", x = Inf, y = Inf, label = paste("stress = ", S), vjust = 2, hjust = 2) +
  # annotate("text", x = Inf, y = Inf, label = paste("k = ", K), vjust = 2, hjust = 2) +
   geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = wq_vecs, linewidth =1) +
  geom_text(data = wq_vecs, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            fontface = "bold", label = row.names(wq_vecs)) +
  scale_color_manual(values = site_colors) +
  scale_fill_manual(values = site_colors) +
  theme_minimal() 

simper_site <- simper(comm = fish_L_full[4:45], 
       group = fish_L_full$site, 
       permutations = 9999)
summary(simper_site)
