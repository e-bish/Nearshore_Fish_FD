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
metrics <- c("Species_Richness", "SESFRic", "FEve", "FDiv", "SESFDis")

mFD_results <- mFD_results %>% 
 inner_join(SES_tab)

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
int.plot_permanova_list <- map(metrics, calc_int.plot_permanova)
int.plot_permanova_df <- list_cbind(int.plot_permanova_list, name_repair = "minimal")
int.plot_permanova_df <- int.plot_permanova_df %>%
  rownames_to_column(var = "X")

# write_csv(int.plot_permanova_df, here("data", "int.plot_permanova_df.csv"))
int.plot_permanova_df <- read_csv(here("data", "int.plot_permanova_df.csv"))

#test single factors 
plot_permanova_list <- map(metrics, calc_plot_permanova)
plot_permanova_df <- list_cbind(plot_permanova_list, name_repair = "minimal")
plot_permanova_df <- plot_permanova_df %>%
  rownames_to_column(var = "X")

# write_csv(plot_permanova_df, here("data", "plot_permanova_df.csv"))
plot_permanova_df <- read_csv(here("data", "plot_permanova_df.csv"))

## follow up pairwise tests 

pairwise.adonis2(mFD_results['Species_Richness'] ~ site, data = mFD_results, permutations = 999, by = "margin", method = "euclidean")
#A: COR
#B: FAM, SHR, TUR, DOK
#C: DOK, EDG, FAM, SHR

pairwise.adonis2(mFD_results['FEve'] ~ site, data = mFD_results, permutations = 999, by = "margin", method = "euclidean")
# #A: TUR, SHR, COR, EDG, FAM
# #B: DOK, EDG, FAM

mFD_df <- as.data.frame(mFD_results)

pairwise.adonis2(mFD_df['Species_Richness'] ~ year, 
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
  mutate(letters = ifelse(metric == "SESFDis" , NA, letters))

ggplot(data = site_dispersions_signif) +
  geom_boxplot(aes(x = site, y = value)) +
  geom_point(aes(y = ymax, x = site), color = "white", size = 0) +
  geom_text(aes(label = letters, y = max, x = site), vjust = -0.5) +
  theme_classic() + 
  labs(y = "Dispersion") +
  facet_wrap(~metric, ncol = 1, scales = "free")
  

#### rda for site ####
mod1 <- rda(mFD_results[8:11] ~ region, data = mFD_results)
plot(mod1)

rda_scores1 <- scores(mod1)
sites_scores1 <- as.data.frame(rda_scores1$sites)
biplot_scores1 <- as.data.frame(rda_scores1$species)

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

site_colors <- rev(c("#8c510a","#d8b365", 
                     # "#f6e8c1",
                     "lightgoldenrod",
                     # "#c7eae8",
                     "lightblue",
                     "#5ab4ac", "#01665e"))

ggplot(data = points1, aes(x = RDA1, y = PC1)) +
  geom_point(aes(color = site, shape = ipa), size = 3) + 
  stat_ellipse(aes(group = site, color = site),
               linetype = "dashed", show.legend = FALSE) +
  # geom_segment(data=biplot_scores1,
  #              aes(x = 0, y = 0, xend=RDA1, yend=RDA2),
  #              arrow=arrow(length = unit(0.01, "npc")),
  #              lwd=0.75) +
  # geom_text(data=biplot_scores1,
  #           aes(x=RDA1*0.9,
  #               y=RDA2*0.9,
  #               label=metrics[-1]),
  #           nudge_x = c(-0.2, -0.3, -0.2, -0.2), 
  #           nudge_y = c(0.1, -0.05, 0, 0),
  #           size=4) +
  scale_color_manual(values = site_colors) +
  labs(shape = "Shoreline Condition", 
       color = "Site") +
  theme_minimal() 

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
nmds <- metaMDS(mFD_results[8:11], 
                distance="euc", k= 3, trymax=1000, plot = FALSE)

points <- data.frame(nmds$points) %>% 
  cbind(mFD_results %>% select(site:year)) 
  # anti_join(add_small_samples_back)

hulls <- points %>%
  group_by(site) %>% 
  slice(chull(MDS1,MDS2))

ggplot() +
  geom_point(data = points, aes(x = MDS1, y = MDS2, color = site, shape = ipa), size = 3) + 
  geom_text_repel(data = points, aes(x = MDS1, y = MDS2, color = site, label = year), max.overlaps = 20) +
  # geom_polygon(data = hulls, aes(x = MDS1, y = MDS2, fill = site), alpha = 0.2) +
  # annotate("text", x = Inf, y = Inf, label = paste("stress = ", S), vjust = 2, hjust = 2) +
  # annotate("text", x = Inf, y = Inf, label = paste("k = ", K), vjust = 2, hjust = 2) +
  theme_minimal() 
