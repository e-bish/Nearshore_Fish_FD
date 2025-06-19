data(BCI)
S <- specnumber(BCI) # observed number of species
(raremax <- min(rowSums(BCI)))
Srare <- rarefy(BCI, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(BCI, step = 20, sample = raremax, col = "blue", cex = 0.6)




#load libraries
library(tidyverse)
library(here)
library(janitor)
library(vegan)
library(rfishbase)
library(GGally)

#load tidy fish data 
load(here("data", "net_core.Rdata")) 

## create the abundance matrix 
# extract unique sampling events to quantify sample effort (slightly unbalanced between depths and shorelines)
sampling_events <- net_core %>% 
  select(year, month, day, site, ipa, station) %>% 
  distinct() %>% 
  filter(!(site == "COR" & year == "2019" & month == "Apr" & day == "30")) %>% #these COR seem like data entry errors because they were only sampled at one station and we already had complete sampling at COR in april and may. No fish recorded in either entry
  filter(!(site == "COR" & year == "2019" & month == "May" & day == "01")) %>% 
  filter(!(site == "TUR" & year == "2018" & month == "Jul" & day == "11")) %>%  #this is an incomplete sampling event. Complete sampling at TUR occured on 7/12/18
  filter(!(site == "TUR" & year == "2018" & month == "Sept" & day == "11")) %>% #another month with repeat sampling; keeping only the second september TUR sampling event
  group_by(year, site, ipa) %>% 
  summarize(no_net_sets = n()) %>% 
  ungroup()

#expand species by sampling event (we're not including day here so we don't have to worry about removing the events we removed above)
expand_species <- net_core %>%
  expand(nesting(year, site, ipa), ComName) %>%
  filter(!is.na(ComName))

#version with catch per set by year
fish_L_prelim <- net_core %>% #L is referring to the RLQ analysis
  filter(!is.na(ComName)) %>% 
  group_by(site, ipa, year, ComName) %>% 
  summarize(spp_sum = sum(species_count)) %>% #sum across depth stations and months within each shoreline type
  ungroup() %>%
  full_join(sampling_events) %>% 
  # mutate(catch_per_set = spp_sum/no_net_sets) %>% #on a couple of occasions we did not sample at all three depths
  filter(!is.na(ComName)) %>% #these are accounted for in the next step
  full_join(expand_species) %>% #add back in all of the events so we can capture 0s
  select(!c(no_net_sets)) %>% 
  mutate(ComName = replace(ComName, ComName == "Pacific Sandfish", "Pacific sandfish")) %>% #if it's capitolized it gets confused about the proper order
  arrange(site, ipa, year, ComName) %>% 
  # mutate(catch_per_set  = replace_na(catch_per_set , 0)) %>% 
  pivot_wider(names_from = ComName, values_from = spp_sum, values_fill = 0) %>% 
  clean_names() %>%
  ungroup() %>% 
  replace(is.na(.), 0) %>% 
  mutate(sum = rowSums(across(where(is.numeric))))

fish_L_prelim %>% 
  filter(sum < 50) %>% 
  ggplot(aes(x = sum)) +
  geom_histogram()

fish_L_rare <- fish_L_prelim %>% 
  filter(!sum < 20) 

rows_removed <- fish_L_prelim %>% 
  filter(sum < 20) 
#we didnt catch enough fish in these samples to say whether the diversity is any lower or higher than other sites

rare_meta <- fish_L_rare %>% 
  select(c(site:year))

fish_L_rare <- fish_L_rare %>% 
  select(!c(site:year, sum))


S <- specnumber(fish_L_rare) # observed number of species
(raremax <- min(rowSums(fish_L_rare)))
Srare <- rarefy(fish_L_rare, 20)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(BCI, step = 20, sample = raremax, col = "blue", cex = 0.6)

fish_L_rarefied <- rrarefy(fish_L_rare, 20)
SR_rarefied <- fish_L_rarefied %>% 
  decostand(method = "pa") %>% 
  as.data.frame() %>% 
  bind_cols(rare_meta) %>% 
  pivot_longer(!site:year, names_to = "species", values_to = "rare_numb") %>% 
  group_by(site, ipa, year) %>% 
  summarize(SR = sum(rare_numb))

SR_rarefied2 <- specnumber(fish_L_rarefied, MARGIN = 1)

rare_meta %>% 
  mutate(SR = SR_rarefied2) %>% 
  ggplot(aes(x = site, y = SR)) +
  geom_boxplot()

SR_rarefied %>% 
  group_by(site) %>% 
  summarize(mean(SR))


mFD_results_fewer <- mFD_results %>% 
  anti_join(rows_removed, by = join_by(site, ipa, year)) #this has one fewer, need to figure out why

#next up, try to run the permanova and see if it still works with this sample size

adonis2(mFD_results_fewer$FRic ~ site + ipa + year, 
        data = mFD_results_fewer, 
        permutations = plot_shuffle,
        method = "euclidean", 
        by = "margin")
#doesn't work if you dont have equal numbers of years