
library(picante)

load(here("data", "fish_L_full.Rdata"))

fish_L_full <- fish_L_full %>% 
  mutate(sample = paste(site, ipa, year, sep = "_"), .after = year) %>% 
  select(!1:3) %>% 
  column_to_rownames(var = "sample")

null_L <- list(fish_L_full)
n_iter <- 1000 #observed data + 999 permutations
#
# #frequency is C4 null model in Gotzenberger et al. that works well for detection of environmental filtering
# #independentswap seems to be the method that other FD papers have used (e.g., Zhang) and is recommended by Swenson
for (i in 2:n_iter) {
  set.seed(i)
  null_L[[i]] <- randomizeMatrix(fish_L_full, null.model ="independentswap", iterations = 1000)
}

# identify_rwfs <- function (df) {
#   rows_w_few_spp <- df %>%
#     decostand(method = "pa") %>%
#     filter(rowSums(.) < 4)
#   
#   return(rows_w_few_spp)
# }
# 
# null_L_dfs <- map(null_L, as.data.frame)
# null_L_rwfs <- map(null_L_dfs, identify_rwfs)

load(here("data", "rows_w_few_spp.Rdata"))

samples_w_few_spp <- rownames(rows_w_few_spp)

remove_rwfs <- function(df) {
  
  df2 <- as.data.frame(df)
  
  out <- df2 %>% filter(!rownames(.) %in% samples_w_few_spp)
  
  return(out)
}

null_L_sample <- map(null_L, remove_rwfs)

load(here("data", "fish.list.Rdata"))

#### run with mFD ####
traits_cat <- data.frame(trait_name = colnames(fish.list$trait),
                         trait_type = c("Q", "N", "N", "N", "N"))

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

trait_space <- space_quality$"details_fspaces"$"sp_pc_coord"

FD_null_results <- list()
for (i in 1:n_iter){
  
  null_fishFD <-  alpha.fd.multidim(sp_faxes_coord = trait_space[ , c("PC1", "PC2", "PC3")],
                                    asb_sp_w = data.matrix(null_L_sample[[i]]),
                                    ind_vect = c("fric","fdis"),
                                    scaling = FALSE,
                                    check_input = TRUE,
                                    details_returned = FALSE)
  
  null_mFD_values <- null_fishFD$"functional_diversity_indices"
  
  FD_null_df <- null_mFD_values %>%
    as_tibble(rownames = "sample")
  
  FD_null_results <- rbind(FD_null_results, FD_null_df)
}

colnames(FD_null_results)[2:4] <- c("Species_Richness", "FRic", "FDis")

add_small_samples_back <- rows_w_few_spp %>% 
  rownames_to_column(var = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa", "year")) %>% 
  mutate(Species_Richness = rowSums(across(where(is.numeric)))) %>% 
  select(site, ipa, year, Species_Richness) %>% 
  rename(Species_Richness_mean = Species_Richness)

FD_null_means <- FD_null_results %>% 
  group_by(sample) %>% 
  summarize(across(where(is.numeric), list(mean = mean, sd = sd))) %>% 
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

# FD_means <- mFD_results %>% 
#   group_by(site) %>% 
#   summarize(across(where(is.numeric), mean)) %>% 
#   relocate(FRic,.after = Species_Richness) %>% 
#   relocate(FDis, .after = FDiv)

SES_tab <- as.data.frame(mFD_results[1:6])
SES_tab[,7] <- FD_null_means$Species_Richness_mean
SES_tab[,8] <- (mFD_results$FRic - FD_null_means$FRic_mean) / FD_null_means$FRic_sd
SES_tab[,9] <- (mFD_results$FDis - FD_null_means$FDis_mean) / FD_null_means$FDis_sd

SES_tab[is.na(SES_tab)] <- 0

colnames(SES_tab)[7:9] <- c("Species_Richness", "SESFRic", "SESFDis")


all(SES_tab$Species_Richness == mFD_results$Species_Richness) #TRUE

