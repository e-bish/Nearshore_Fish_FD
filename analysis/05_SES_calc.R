library(tidyverse)
library(here)
library(mFD)
library(picante)
library(DescTools)

load(here("data", "fish_L_full.Rdata")) #created in 03_create_matrices
load(here("data","mFD_results.Rdata")) #created in 04_diversity_analysis

fish_L_full <- fish_L_full %>% 
  mutate(sample = paste(site, ipa, year, sep = "_"), .after = year) %>% 
  select(!1:3) %>% 
  column_to_rownames(var = "sample")

null_L <- list(fish_L_full)
n_iter <- 1000 #observed data + 999 permutations

# #frequency is C4 null model in Gotzenberger et al. that works well for detection of environmental filtering
# #independentswap seems to be the method that other FD papers have used (e.g., Zhang) and is recommended by Swenson
for (i in 2:n_iter) {
  set.seed(i)
  null_L[[i]] <- randomizeMatrix(fish_L_full, null.model ="independentswap", iterations = 1000)
}

load(here("data", "rows_w_few_spp.Rdata"))

samples_w_few_spp <- rownames(rows_w_few_spp)

remove_rwfs <- function(df) {
  
  df2 <- as.data.frame(df)
  
  out <- df2 %>% filter(!rownames(.) %in% samples_w_few_spp)
  
  return(out)
}

null_L_sample <- map(null_L, remove_rwfs)

#### calculate FD metrics for null matrices with mFD
load(here("data", "fish.list.Rdata"))
load(here("data", "trait_space.Rda"))

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

#check warning on sample df. all dfs should have the same column and row sums
# test_war <- colSums(null_L[[15]]) #warnings appear erroneous

colnames(FD_null_results)[2:4] <- c("Species_Richness", "FRic", "FDis")

hist(FD_null_results$FRic)
hist(sqrt(FD_null_results$FRic))
Skew(sqrt(FD_null_results$FRic))

hist(FD_null_results$FDis)
Skew(FD_null_results$FDis)

#check skewness of the null distributions

check_skew <- function(metric, trans) {
  
  samples <- FD_null_results %>% 
    select(sample) %>% 
    distinct() %>% 
    pull(sample)
  
  skew_vals <- c()
  
  if (trans == FALSE) {
    for (i in samples) {
      df <- FD_null_results %>% 
        filter(sample == i)  
      
      skew_vals[i] <- Skew(df[[metric]])
    }
  } else {
    for (i in samples) {
      df <- FD_null_results %>% 
        filter(sample == i)  
      
      skew_vals[i] <- Skew(sqrt(df[[metric]]))
    }
  }
  return(boxplot(skew_vals))
}

check_skew("FRic", trans = FALSE) #transformation needed 
check_skew("FRic", trans = TRUE)
check_skew("FDis", trans = FALSE) #no transformation needed (skew values between -1 & 1)

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

SES_tab <- as.data.frame(mFD_results[1:6])
SES_tab[,7] <- FD_null_means$Species_Richness_mean
SES_tab[,8] <- (sqrt(mFD_results$FRic) - sqrt(FD_null_means$FRic_mean)) / sqrt(FD_null_means$FRic_sd)
SES_tab[,9] <- (mFD_results$FDis - FD_null_means$FDis_mean) / FD_null_means$FDis_sd

SES_tab[is.na(SES_tab)] <- 0

colnames(SES_tab)[7:9] <- c("Species_Richness", "SESFRic", "SESFDis")

all(SES_tab$Species_Richness == mFD_results$Species_Richness) #TRUE

save(SES_tab, file = here("data", "SES_tab.Rdata"))
