### Code to import data from SOS shoreline restoration project phases 1 (2018/2019) and 2 (2021/2022)

# Load libraries
library(tidyverse)
library(here)

# all of the raw data were downloaded from the survey team's shared google drive and saved as csv files
net_import <- here::here("data","raw", "raw.18.19.csv") %>% read_csv()
net_import2 <- here::here("data","raw", "raw.21.csv") %>% read_csv()
net_import3 <- here::here("data","raw", "raw.22.csv") %>% read_csv()

net_import <-net_import %>%
  separate(date, into = c("year","month", "day"), sep = "-") %>%
  mutate(month = str_pad(month, 2, side = c("left"), pad = "0")) %>%
  mutate(day = str_pad(day, 2, side = c("left"), pad = "0")) %>%
  mutate(tax_group = replace(tax_group, species == "Tube Snout", "Aulorhynchus")) %>% 
  mutate(species = replace(species, species == "Gunnel Sp", "UnID Gunnel")) %>% 
  mutate(species = replace(species, species == "Jellyfish Sp", "UnID Jellyfish")) %>% 
  mutate(ipa = as.factor(ipa)) %>%
  select(-c(transect_notes))

net <- net_import %>%
  mutate(species = replace_na(species, "none"),length_mm = replace_na(length_mm, 0)) %>%
  group_by(year, month, day, site, ipa, station, org_type, tax_group,species) %>% 
  mutate(count = 1) %>%
  summarize(species_count=sum(count), mean_length_mm = mean(length_mm)) %>%
  ungroup() 

#write to csv
write_csv(net, here("data","net_18.19.csv")) 

###############################################################################
## 2021 data

net_2021 <- net_import2 %>%
  mutate(month = str_pad(month, width = 2, pad = "0")) %>% 
  mutate(day = str_pad(day, width = 2, pad = "0")) %>% 
  mutate(length_mm = ifelse(is.na(length_mm), 0, length_mm)) %>% 
  mutate(species = ifelse(is.na(species), "none", species)) %>% 
  group_by(year, month, day, site, ipa, station, org_type, tax_group, species) %>%
  summarize(species_count = n(), mean_length_mm = mean(length_mm)) %>% 
  ungroup() %>% 
  mutate(month = if_else(site == "MA", "06", month)) %>% # we did a July 1st survey at Maylor that we want to count as a June survey
  mutate(ipa = ifelse(site == "TL" & ipa == "Restored", "Natural", ipa)) %>% #fix Titlow misdesignations
  mutate(ipa = ifelse(site == "TL" & ipa == "Armored", "Restored", ipa)) %>% 
  mutate(ipa = ifelse(site == "TL" & ipa == "Armored_2", "Armored", ipa))

#write to csv
write_csv(net_2021, here("data","net_2021.csv"))

###############################################################################
## 2022 data

net_2022 <- net_import3 %>%
  mutate(month = str_pad(month, width = 2, pad = "0")) %>% 
  mutate(day = str_pad(day, width = 2, pad = "0")) %>% 
  mutate(length_mm = ifelse(is.na(length_mm), 0, length_mm)) %>% 
  mutate(species = ifelse(is.na(species), "none", species)) %>% 
  group_by(year, month, day, site, ipa, station, org_type, tax_group, species) %>%
  summarize(species_count = n(), mean_length_mm = mean(length_mm)) %>% 
  ungroup() 

#write to csv
write_csv(net_2022, here("data","net_2022.csv"))
