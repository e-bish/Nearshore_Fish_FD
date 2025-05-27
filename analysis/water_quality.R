# Load libraries
library(tidyverse)
library(here)
library(stringr)


#read csv
wq_import <- here::here("data", "raw","wq_import.csv") %>% read_csv()
SOS_core_sites <- c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")

wq_tb <- wq_import %>% 
  filter(site %in% SOS_core_sites) %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG"))) %>% 
  mutate(ipa = replace(ipa, site == "TUR" & ipa == "Restored", "Natural2")) %>%  #no restoration at Turn Island
  mutate(secchi_depth_m = replace(secchi_depth_m, secchi_depth_m =="NULL", NA)) %>% 
  mutate(secchi_depth_m = as.numeric(unlist(secchi_depth_m))) %>% 
  suppressWarnings() %>% 
  select(!notes) %>% 
  pivot_longer(c(secchi_depth_m, do_mg_l, salinity_ppm, temperature), names_to = "metric") %>% 
  group_by(site, metric) %>% 
  summarize(mean = round(mean(value, na.rm = TRUE), 2), min = min(value, na.rm = TRUE), max = max(value, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(value = paste0(mean, " (", min, "-", max, ")")) %>% 
  select(site, metric, value) %>% 
  pivot_wider(names_from = metric, values_from = value)


write_csv(wq_tb, file = here("data", "wq_tb.csv"))
