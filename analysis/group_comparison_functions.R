
#### calculate plot scale permanova with interaction terms
calc_int.plot_permanova <- function(metric) {
  result <- adonis2(mFD_results[metric] ~ ipa*site*year - ipa:site:year,
                    data = mFD_results,
                    permutations = plot_shuffle,
                    by = "margin",
                    method = "euclidean")
  
  colnames(result) <- paste(colnames(result), metric, sep = '_')
  
  output <- as.data.frame(result)
  
  return(output)
}

#### calculate plot scale permanova with single factors
calc_plot_permanova <- function(metric) {
  result <- adonis2(mFD_results[metric] ~ ipa+site+year,
                    data = mFD_results,
                    permutations = plot_shuffle,
                    by = "margin",
                    method = "euclidean")

  colnames(result) <- paste(colnames(result), metric, sep = '_')

  output <- as.data.frame(result)

  return(output)
}

### calculate whole plot scale permanova for region with manual p values
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



#### calculate dispersion for each group
compare_dispersion <- function(metric, scale) {
  
  dist <- vegdist(mFD_results[[metric]], method = "euclidean")
  
  disp <- betadisper(dist, mFD_results[[scale]], type = c("median"), bias.adjust = FALSE,
                     sqrt.dist = FALSE, add = FALSE)
  
  output <- as.data.frame(disp$distances)
  
  names(output) <- metric
  
  return(output)
  
}

### calculate whole plot scale permanova for veg with manual p values
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

#### calculate dispersion for the results dataframe without zeros
compare_sub_dispersion <- function(metric, scale) {
  
  dist <- vegdist(mFD_sub[[metric]], method = "euclidean")
  
  disp <- betadisper(dist, mFD_sub[[scale]], type = c("median"), bias.adjust = FALSE,
                     sqrt.dist = FALSE, add = FALSE)
  
  output <- as.data.frame(disp$distances)
  
  names(output) <- metric
  
  return(output)
  
}

#### calculate dispersion p values
compare_disp_pval <- function(metric, scale, shuffle) {
  
  dist <- vegdist(mFD_results[[metric]], method = "euclidean")
  
  disp <- betadisper(dist, mFD_results[[scale]], type = c("median"), bias.adjust = FALSE,
                     sqrt.dist = FALSE, add = FALSE)
  
  permu <- permutest(disp, pairwise = TRUE, permutations = shuffle)
  
  pvals <- as.data.frame(permu[["pairwise"]]["permuted"]) %>% 
    rownames_to_column(var = "pair")
  
  names(pvals)[2] <- metric
  
  return(pvals)
  
}

#### extract letters for significant groupings
create_letter_df <- function(metric) {
  
  disp <- site_disp_pvals %>% 
    pull(metric, pairs)
  
  letters <- as.data.frame(multcompLetters(disp)$Letters)
  
  names(letters) <- metric
  
  return(letters)
}
