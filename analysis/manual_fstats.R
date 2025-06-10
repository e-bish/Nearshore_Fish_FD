set.seed(2025)

permutations <- rbind(1:nrow(mFD_results),
                      shuffleSet(n = nrow(mFD_results), 
                                 control = whole_plot_shuffle, nset = 999))

calc_region_pvals <- function (metric) {
  
  mod_result <- adonis2(mFD_results[metric] ~ region + site,
                        data = mFD_results,
                        method = "euclidean",
                        permutations = 0, 
                        by = "terms")
  
  output <- as.data.frame(mod_result)
  
  results <- matrix(nrow = nrow(permutations), ncol = 4)
  
  colnames(results) <- c("region", "site", "residual", "total")
  
  for (i in 1:nrow(permutations)) {
    temp.data <- mFD_results[permutations[i, ], ]
    temp <- adonis2(mFD_results[metric] ~ region + site,
                    data = temp.data,
                    method = "euclidean",
                    permutations = 0, 
                    by = "terms")
    results[i,] <- t(temp$SumOfSqs)
  }
  
  aov_result <- summary(aov(mFD_results[[metric]] ~ region + site + Error(site), data = mFD_results))
  
  results <- results %>% 
    data.frame() %>% 
    mutate(F.region = (region / aov_result[["Error: site"]][[1]][1,1]) / (site / aov_result[["Error: site"]][[1]][2,1]),
           F.site = (site / aov_result[["Error: site"]][[1]][2,1]) / (residual / aov_result[["Error: Within"]][[1]][1,1]))
  
  region_pval <- with(results, sum(F.region >= F.region[1]) / length(F.region))
  
  output[1,'Pr(>F)'] <- region_pval
  output[2,'Pr(>F)'] <- NA
  
  colnames(output) <- paste(colnames(output), metric, sep = '_')
  
  return(output)
}

region_pvals <- map_dfc(metrics, calc_region_pvals) %>% 
  rename(Degrees_Freedom = "Df_Species_Richness") %>% 
  select(!contains("Df"))

calc_veg_pvals <- function (metric) {
  
  mod_result <- adonis2(mFD_results[metric] ~ veg + site,
                        data = mFD_results,
                        method = "euclidean",
                        permutations = 0, 
                        by = "terms")
  
  output <- as.data.frame(mod_result)
  
  results <- matrix(nrow = nrow(permutations), ncol = 4)
  
  colnames(results) <- c("veg", "site", "residual", "total")
  
  for (i in 1:nrow(permutations)) {
    temp.data <- mFD_results[permutations[i, ], ]
    temp <- adonis2(mFD_results[metric] ~ veg + site,
                    data = temp.data,
                    method = "euclidean",
                    permutations = 0, 
                    by = "terms")
    results[i,] <- t(temp$SumOfSqs)
  }
  
  aov_result <- summary(aov(mFD_results[[metric]] ~ veg + site + Error(site), data = mFD_results))
  
  results <- results %>% 
    data.frame() %>% 
    mutate(F.veg = (veg / aov_result[["Error: site"]][[1]][1,1]) / (site / aov_result[["Error: site"]][[1]][2,1]),
           F.site = (site / aov_result[["Error: site"]][[1]][2,1]) / (residual / aov_result[["Error: Within"]][[1]][1,1]))
  
  veg_pval <- with(results, sum(F.veg >= F.veg[1]) / length(F.veg))
  
  output[1,'Pr(>F)'] <- veg_pval
  output[2,'Pr(>F)'] <- NA
  
  colnames(output) <- paste(colnames(output), metric, sep = '_')
  
  return(output)
}

veg_pvals <- map_dfc(metrics, calc_veg_pvals) %>% 
  rename(Degrees_Freedom = "Df_Species_Richness") %>% 
  select(!contains("Df"))
