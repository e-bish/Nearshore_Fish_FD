#load libraries
library(tidyverse)
library(here)
library(janitor)

#load tidy fish data 
load(here("data", "net_core.Rdata")) 

# extract unique sampling events to quantify sample effort (slightly unbalanced between depths and shorelines)
sampling_events1 <- net_core %>% 
  select(year, month, day, site, ipa, station) %>% 
  distinct() 

# confirm same # of rows by day and by month 
sampling_events1 %>%  group_by(year, month, day) %>% nrow()
sampling_events1 %>% group_by(year, month) %>% nrow()

# we sampled once per month at each site, and confirmed with the previous test that day is redundant with month so it can be excluded
sampling_events <- sampling_events1 %>% 
  group_by(year, site, month, ipa) %>% 
  summarize(no_net_sets = n()) %>% 
  ungroup()

#### No Data Aggregation ####

#expand species by sampling event 
expand_species_events <- net_core %>%
  expand(nesting(year, month, site, ipa, station), ComName) %>%
  filter(!is.na(ComName))

#version with catch per set by sample event
fish_L_full_events <- net_core %>% 
  filter(!is.na(ComName)) %>% 
  group_by(site, ipa, station, year, month, ComName) %>%
  summarize(spp_sum = sum(species_count)) %>% #sum across depth stations and months within each shoreline type
  ungroup() %>%
  full_join(sampling_events) %>% 
  mutate(catch_per_set = spp_sum/no_net_sets) %>% #on a couple of occasions we did not sample at all three depths
  filter(!is.na(ComName)) %>% #these are accounted for in the next step
  full_join(expand_species_events) %>% #add back in all of the events so we can capture 0s
  select(!c(spp_sum, no_net_sets)) %>%
  mutate(ComName = replace(ComName, ComName == "Pacific Sandfish", "Pacific sandfish")) %>% #if it's capitolized it gets confused about the proper order
  arrange(site, ipa, station, year, month, ComName) %>% 
  mutate(catch_per_set  = replace_na(catch_per_set, 0)) %>% 
  pivot_wider(names_from = ComName, values_from = catch_per_set, values_fill = 0) %>% 
  clean_names() %>% 
  ungroup()

fish_L_sample_events <- fish_L_full_events %>% 
  mutate(sample = paste(site, ipa, station, year, month, sep = "_"), .after = month) %>% 
  select(!1:5) %>% 
  column_to_rownames(var = "sample")

#check the abundance matrix for samples with few species
events_rows_w_few_spp <- fish_L_sample_events %>%
  decostand(method = "pa") %>%
  filter(rowSums(.) < 4)

sample_events_w_few_spp <- rownames(events_rows_w_few_spp)

length(sample_events_w_few_spp) / nrow(fish_L_full_events)

#### Aggregate by Depth ####

#expand species by sampling event 
expand_species_depth <- net_core %>%
  expand(nesting(year, month, site, ipa), ComName) %>%
  filter(!is.na(ComName))

#version with catch per set by sample event
fish_L_full_depth <- net_core %>% 
  filter(!is.na(ComName)) %>% 
  group_by(site, ipa, year, month, ComName) %>% 
  summarize(spp_sum = sum(species_count)) %>% 
  ungroup() %>%
  full_join(sampling_events) %>% 
  mutate(catch_per_set = spp_sum/no_net_sets) %>% #on a couple of occasions we did not sample at all three depths
  filter(!is.na(ComName)) %>% #these are accounted for in the next step
  full_join(expand_species_depth) %>% #add back in all of the events so we can capture 0s
  select(!c(spp_sum, no_net_sets)) %>%
  mutate(ComName = replace(ComName, ComName == "Pacific Sandfish", "Pacific sandfish")) %>% #if it's capitolized it gets confused about the proper order
  arrange(site, ipa, year, month, ComName) %>% 
  mutate(catch_per_set  = replace_na(catch_per_set, 0)) %>% 
  pivot_wider(names_from = ComName, values_from = catch_per_set, values_fill = 0) %>% 
  clean_names() %>% 
  ungroup()

fish_L_sample_depth <- fish_L_full_depth%>% 
  mutate(sample = paste(site, ipa, year, month, sep = "_"), .after = month) %>% 
  select(!1:4) %>% 
  column_to_rownames(var = "sample")

#check the abundance matrix for samples with few species
depth_rows_w_few_spp <- fish_L_sample_depth %>%
  decostand(method = "pa") %>%
  filter(rowSums(.) < 4)

sample_depth_w_few_spp <- rownames(depth_rows_w_few_spp)

length(sample_depth_w_few_spp) / nrow(fish_L_full_depth)

