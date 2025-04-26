## Create Matrices 

#load libraries
library(tidyverse)
library(here)
library(janitor)
library(vegan)
library(rfishbase)
library(GGally)

#load tidy fish data 
load(here("data", "net_core.Rdata")) 

net_core <- net_core %>% 
  mutate(season = ifelse(month %in% c("May", "Jun", "Jul"), "peak", "shoulder"), .after = day)

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
  expand(nesting(year, month, season, site, ipa), ComName) %>%
  filter(!is.na(ComName))

#version with catch per set by year
fish_L_full <- net_core %>% #L is referring to the RLQ analysis
  filter(!is.na(ComName)) %>% 
  group_by(year, month, site, ipa, ComName) %>% 
  summarize(spp_sum = sum(species_count)) %>% #sum across depth stations within each shoreline type
  ungroup() %>%
  full_join(sampling_events) %>% 
  mutate(catch_per_set = spp_sum/no_net_sets) %>% #on a couple of occasions we did not sample at all three depths
  filter(!is.na(ComName)) %>% #these are accounted for in the next step
  full_join(expand_species) %>% #add back in all of the events so we can capture 0s
  mutate(ComName = replace(ComName, ComName == "Pacific Sandfish", "Pacific sandfish")) %>% #if it's capitolized it gets confused about the proper order
  arrange(year, month, site, ipa, ComName) %>% 
  mutate(catch_per_set  = replace_na(catch_per_set , 0)) %>% 
  group_by(site, ipa, month, ComName) %>% 
  summarize(avg_cps = mean(catch_per_set)) %>% #average across months and years within a season for each ipa
  pivot_wider(names_from = ComName, values_from = avg_cps, values_fill = 0) %>% 
  clean_names() %>% 
  ungroup() 

fish_L_sample <- fish_L_full %>% 
  mutate(sample = paste(site, ipa, month, sep = "_"), .after = month) %>% 
  select(!1:3) %>% 
  column_to_rownames(var = "sample")

rows_w_few_spp <- fish_L_sample %>%
  decostand(method = "pa") %>%
  filter(rowSums(.) < 4)

save(rows_w_few_spp, file = here("data", "rows_w_few_spp_month.Rdata"))

samples_w_few_spp <- rownames(rows_w_few_spp)

#remove samples that have less than 4 (this number depends on the number of axes we want to retain for the trait space)
fish_L <- fish_L_sample %>%
  filter(!rownames(.) %in% samples_w_few_spp) #removes 37 observations by month

fish_traits <- read_csv("data/fish_traits.csv")

fish_Q <- fish_traits %>% 
  select(-c(Species, ComName)) %>% 
  mutate_if(is.character, as.factor) %>% 
  clean_names() %>% 
  as.data.frame()

rownames(fish_Q) <- colnames(fish_L)

confirm_proper_names <- fish_Q %>% 
  mutate(Species = fish_traits$ComName, .before = mean_length_mm)

fish_meta <- fish_L_sample %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa", "season"), cols_remove = TRUE) %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")),
         region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"), 
         veg = ifelse(site %in% c("TUR", "COR", "SHR"), "present", "absent"), .after = site,
         shoreline = paste0(site, case_when(ipa == "Armored" ~ "A", 
                                            ipa == "Restored" ~ "R", 
                                            ipa == "Natural2" ~ "N2",
                                            TRUE ~ "N")),
         across(where(is.numeric), ~replace_na(., 0))) %>% 
  mutate(ipa = ifelse(ipa == "Natural2", "Natural", ipa)) %>% 
  arrange(site,ipa, shoreline, season) %>% 
  select(site:season)

fish.list.month <- list("trait" = fish_Q, 
                  "abund" = fish_L,
                  "meta" = fish_meta) 

save(fish.list.month, file = here("data", "fish.list.month.Rdata")) #last saved 04/19/25
