#FD analysis

#load libraries
library(tidyverse)
library(here)
library(ggrepel)
library(mFD)
library(ggordiplots)
library(patchwork)

#### Start here for FD analysis ####

set.seed(2025)

load(here("data", "fish.list.Rdata"))
metrics <- c("Species_Richness", "FRic", "FEve", "FDiv", "FDis")


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
save(trait_space, file = here("data", "trait_space.Rda"))

#test correlations between traits and axes
traits_correlations <- traits.faxes.cor(
  sp_tr          = fish.list$trait, 
  sp_faxes_coord = trait_space[ , c("PC1", "PC2", "PC3")], 
  plot           = TRUE)

traits_correlations$"tr_faxes_stat"[which(traits_correlations$"tr_faxes_stat"$"p.value" < 0.05), ]

trait_plot <- traits_correlations$"tr_faxes_plot"


#calculate FD Indices
alpha_fd_indices <- alpha.fd.multidim(sp_faxes_coord = trait_space[ , c("PC1", "PC2", "PC3")],
                                      asb_sp_w = data.matrix(fish.list$abund),
                                      ind_vect = c("fric", "feve", "fdiv", "fdis"),
                                      scaling = FALSE,
                                      check_input = TRUE,
                                      details_returned = TRUE)

#the mFD package uses ape::pcoa() which automatically removes negative eigenvalues rather than applying a correction

mFD_values <- alpha_fd_indices$"functional_diversity_indices"

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

# save(mFD_results, file = here("data","mFD_results.Rdata")) #last saved 6/25/25 
load(here("data","mFD_results.Rdata"))

create_metric_summary_table <- function(metric_ID) {
  tmp <- mFD_results %>% 
    summarize(min = min(.data[[metric_ID]]), 
              max = max(.data[[metric_ID]]), 
              avg = mean(.data[[metric_ID]]), 
              sd = sd(.data[[metric_ID]])) %>% 
    mutate(metric = metric_ID)
  
  return(tmp)
}

summary_stats <- map_dfr(metrics, create_metric_summary_table)

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

#### estimate the contribution of each trait to metrics of FD ####
### based on Stuart-Smith et al. 2013

calculate_FD <- function(trait_mat, trait_vec) {
  
  #### run with mFD ####
  traits_cat <- data.frame(trait_name = colnames(trait_mat),
                           trait_type = trait_vec)
  
  
  dist_mat <- funct.dist(sp_tr = trait_mat, 
                         tr_cat = traits_cat,
                         metric = "gower",
                         weight_type = "equal",
                         stop_if_NA = TRUE)
  
  space_quality <- quality.fspaces(sp_dist = dist_mat,
                                   maxdim_pcoa = 10,
                                   deviation_weighting = c("absolute", "squared"),
                                   fdist_scaling = TRUE,
                                   fdendro = "ward.D2")
  
  trait_space <- space_quality$"details_fspaces"$"sp_pc_coord"
  
  # calculate FD Indices
  alpha_fd_indices <- alpha.fd.multidim(sp_faxes_coord = trait_space[ , c("PC1", "PC2", "PC3")],
                                        asb_sp_w = data.matrix(fish.list$abund),
                                        ind_vect = c("fdis", "feve", "fric", "fdiv"),
                                        scaling = FALSE,
                                        check_input = TRUE,
                                        details_returned = TRUE)
  
  #the mFD package uses ape::pcoa() which automatically removes negative eigenvalues rather than applying a correction
  
  mFD_values <- alpha_fd_indices$"functional_diversity_indices"
  
  return(mFD_values)
}

full_traits <- calculate_FD(fish.list$trait, c("Q", "N", "N", "N", "N"))
minus_length <- calculate_FD(fish.list$trait[,-1], c("N", "N", "N", "N"))
minus_body <- calculate_FD(fish.list$trait[,-2], c("Q", "N", "N", "N"))
minus_pos <- calculate_FD(fish.list$trait[,-3], c("Q", "N", "N", "N"))
minus_fg <- calculate_FD(fish.list$trait[,-5], c("Q", "N", "N", "N"))
minus_mig <- calculate_FD(fish.list$trait[,-4], c("Q", "N", "N", "N"))

mods <- c("minus_length", "minus_body", "minus_pos", "minus_fg", "minus_mig")

extract_R2_tab <- function(metric){
  
  response <- full_traits[,metric]
  
  rsq_df <- data.frame(mods)
  
  rsq_df[1,metric]<- summary(lm(response ~ minus_length[,metric]))$r.squared
  rsq_df[2,metric]<- summary(lm(response ~ minus_mig[,metric]))$r.squared
  rsq_df[3,metric]<- summary(lm(response ~ minus_pos[,metric]))$r.squared
  rsq_df[4,metric]<- summary(lm(response ~ minus_fg[,metric]))$r.squared
  rsq_df[5,metric]<- summary(lm(response ~ minus_body[,metric]))$r.squared
  
  colnames(rsq_df) <- c("mod", paste0(metric, "_rsq"))
  
  return(rsq_df)
}

metrics_short <- c("fric", "feve", "fdiv", "fdis")

combine_R2_tabs <- function(metric) {
  
  tmp <- extract_R2_tab(metric) %>% 
    as_tibble() %>% 
    mutate(contribution = 1 - .[[2]]) %>% 
    mutate_if(is.numeric, round, 2) %>% 
    rename_with(~paste(metric, "cont", sep = "_"), starts_with("contribution"))
  
  return(tmp)
}

combined_df <- map_dfc(metrics_short, combine_R2_tabs) %>% 
  rename(trait = "mod...1") %>% 
  select(-contains("mod")) %>% 
  mutate(trait = c("Body length",
                   "Body transverse shape",
                   "Vertical distribution",
                   "Feeding guild",
                   "Migrations"))

write_csv(combined_df, here("data", "rsq_contributions.csv")) #last saved 6/25/25
