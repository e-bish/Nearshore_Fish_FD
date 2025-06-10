set.seed(2025)

permutations <- rbind(1:nrow(mFD_results),
                      shuffleSet(n = nrow(mFD_results), 
                                 control = whole_plot_shuffle, nset = 999))

calc_region_pvals <- function (metric) {
  
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
  
  explan_pval <- with(results, sum(F.region >= F.region[1]) / length(F.region))
  site_pval <- with(results, sum(F.site >= F.site[1]) / length(F.site)) #this is the same as the output of adonis2 because the error is the residuals term
  
  results_df <- data.frame(metric = metric, region = explan_pval, site = site_pval)
  
  return(results_df)
}

map_dfr(metrics, calc_region_pvals)


calc_veg_pvals <- function (metric) {
  
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
  
  explan_pval <- with(results, sum(F.veg >= F.veg[1]) / length(F.veg))
  site_pval <- with(results, sum(F.site >= F.site[1]) / length(F.site)) #this is the same as the output of adonis2 because the error is the residuals term
  
  results_df <- data.frame(metric = metric, veg = explan_pval, site = site_pval)
  
  return(results_df)
}

map_dfr(metrics, calc_veg_pvals)

## full table
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
  rename("Df_Species_Richness" = "Degrees_Freedom") %>% 
  select(!contains("Df"))

