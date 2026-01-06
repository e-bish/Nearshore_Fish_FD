# Load libraries
library(tidyverse)
library(here)
library(stringr)

site_colors <- rev(c("#8c510a","#d8b365", 
                     # "#f6e8c1",
                     "lightgoldenrod",
                     # "#c7eae8",
                     "lightblue",
                     "#5ab4ac", "#01665e"))

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
  mutate(shoreline = paste0(site, ipa), .after = ipa)

# save(wq_tb, file = here("data", "wq_tb.Rdata"))

wq_tb_long <- wq_tb %>% 
  pivot_longer(c(secchi_depth_m, do_mg_l, salinity_ppm, temperature), names_to = "metric")

#### Fig S2 ####
wq_tb_long %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")),
    shore_short = ifelse(ipa == "Natural2", "TURN2",
                         str_sub(shoreline, end = 4))) %>% 
  arrange(site, ipa) %>% 
  mutate(shore_short = factor(shore_short, levels = unique(shore_short))) %>% 
  ggplot(aes(x = shore_short, y = value, fill = site)) +
  geom_boxplot() +
  theme_bw() + 
  labs(y = "Value", x = "Shoreline", fill = "Site") +
  scale_fill_manual(values = site_colors) +
  theme(axis.text.x = element_text(angle = 70, vjust = 0.5, hjust=0.5)) +
  facet_grid(rows = vars(metric = factor(metric, 
                                         levels = c("temperature",
                                                    "secchi_depth_m",
                                                    "do_mg_l",
                                                    "salinity_ppm"),
                                         labels = c("Temperature\n(C)",
                                                    "Secchi depth\n(m)",
                                                    "Dissolved\noxygen (mg/L)",
                                                     "Salinity (ppm)"))), 
             # cols = vars(year),
             scales = "free_y")


ggsave(here("figures", "Fig_S2.png"), 
       width = 6, height = 5.5, dpi = 300)

wq_tb_long %>% 
  group_by(metric) %>% 
  summarize(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE), 
    min = min(value, na.rm = TRUE), max = max(value, na.rm = TRUE))

site_wq <- wq_tb_long %>% 
  group_by(site, metric) %>% 
  summarize(mean = mean(value, na.rm = TRUE),
            sd = sd(value, na.rm = TRUE),
            min = min(value, na.rm = TRUE), 
            max = max(value, na.rm = TRUE))

site_wq %>% 
  filter(metric == "do_mg_l")

wq_tb_long %>% 
  mutate(region = ifelse(site %in% c("FAM", "TUR", "COR"), "north", "south")) %>% 
  group_by(region, metric) %>% 
  summarize(mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE))

wq_tb_sd_export <- wq_tb_long %>%
  group_by(site, metric) %>%
  summarize(mean = round(mean(value, na.rm = TRUE), 2), sd = round(sd(value, na.rm = TRUE), 2)) %>%
  ungroup() %>%
  mutate(value = paste0(mean, " (", sd, ")")) %>%
  select(site, metric, value) %>%
  pivot_wider(names_from = metric, values_from = value)

wq_tb_export <- wq_tb_long %>%
  group_by(site, metric) %>%
  summarize(mean = round(mean(value, na.rm = TRUE), 2), min = min(value, na.rm = TRUE), max = max(value, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(value = paste0(mean, " (", min, "-", max, ")")) %>%
  select(site, metric, value) %>%
  pivot_wider(names_from = metric, values_from = value)

# write_csv(wq_tb_export, file = here("data", "wq_tb.csv"))





