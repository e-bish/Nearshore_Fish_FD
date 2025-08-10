library(vegan)

data(BCI)
S <- specnumber(BCI) # observed number of species
(raremax <- min(rowSums(BCI)))
Srare <- rarefy(BCI, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(BCI, step = 20, sample = raremax, col = "blue", cex = 0.6)


specnumber(fish_L_full[4:45])

spp_sum_mat <- spp_sum_all %>% 
  select(!site)

fish_grouped <- fish_L_full %>% 
  select(!ipa, year) %>% 
  group_by(site) %>% 
  summarise_at(vars(bay_pipefish:whitespotted_greenling), sum)

fish_spec <- specnumber(spp_sum_mat)  
(raremax <- min(rowSums(spp_sum_mat)))
Srare <- rarefy(spp_sum_mat, raremax,  se = TRUE)

plot(fish_spec, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(BCI, step = 20, sample = raremax, col = "blue", cex = 0.6)
