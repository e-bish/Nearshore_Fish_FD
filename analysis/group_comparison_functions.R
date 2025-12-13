
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


#### extract letters for significant groupings
#### calculate dispersion for each group
compare_dispersion <- function(metric, scale) {
  
  dist <- vegdist(mFD_results[[metric]], method = "euclidean")
  
  disp <- betadisper(dist, mFD_results[[scale]], type = c("median"), bias.adjust = FALSE,
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
create_disp_letter_df <- function(metric) {
  
  disp <- site_disp_pvals %>% 
    pull(FRic, pairs)
  
  letters <- as.data.frame(multcompLetters(disp)$Letters)
  
  names(letters) <- metric
  
  return(letters)
}
