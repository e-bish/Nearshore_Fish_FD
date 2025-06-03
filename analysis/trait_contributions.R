### script to estimate the contribution of each trait to metrics of FD
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
                                        scaling = TRUE,
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
minus_mig <- calculate_FD(fish.list$trait[,-4], c("Q", "N", "N", "N"))
minus_fg <- calculate_FD(fish.list$trait[,-5], c("Q", "N", "N", "N"))

mods <- c("minus_length", "minus_body", "minus_pos", "minus_mig", "minus_fd")

extract_R2_tab <- function(metric){
  
  response <- full_traits[,metric]
  
  rsq_df <- data.frame(mods)
  
  rsq_df[1,metric]<- summary(lm(response ~ minus_length[,metric]))$r.squared
  rsq_df[2,metric]<- summary(lm(response ~ minus_mig[,metric]))$r.squared
  rsq_df[3,metric]<- summary(lm(response ~ minus_pos[,metric]))$r.squared
  rsq_df[4,metric]<- summary(lm(response ~ minus_fd[,metric]))$r.squared
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
                   "Migrations",
                   "Feeding guild"))

write_csv(combined_df, here("data", "rsq_contributions.csv"))
