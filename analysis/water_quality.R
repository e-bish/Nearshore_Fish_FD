# Load libraries
library(tidyverse)
library(here)
library(stringr)
# library(dunn.test)

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

save(wq_tb, file = here("data", "wq_tb.Rdata"))

wq_tb_long <- wq_tb %>% 
  pivot_longer(c(secchi_depth_m, do_mg_l, salinity_ppm, temperature), names_to = "metric")

#test for normality
# shapiro.test(wq_tb$temperature)
# shapiro.test(wq_tb$salinity_ppm)
# shapiro.test(wq_tb$do_mg_l)
# shapiro.test(wq_tb$secchi_depth_m)
# 
# compare_wq <- function(metric, df) {
#   
#   if (df == "shoreline") {
#     
#     kruskal.results <- kruskal.test(wq_tb[[metric]] ~ wq_tb$shoreline)
#     
#     test.dunn <- dunn.test(wq_tb[[metric]], g = wq_tb$shoreline, method = "bonferroni", kw = FALSE)
#     
#   } else {
#     
#     kruskal.results <- kruskal.test(site_wq[[metric]] ~ site_wq$site)
#     
#     test.dunn <- dunn.test(site_wq[[metric]], g = site_wq$site, method = "bonferroni", kw = FALSE)
#     
#   }
#   
#   dunn.results <- test.dunn[["P.adjusted"]] %>% 
#     as.data.frame() %>% 
#     magrittr::set_rownames(test.dunn[["comparisons"]]) %>% 
#     rename(p_val = ".") %>% 
#     filter(p_val < 0.05) 
#     
#   
#   kruskal.pval <- data.frame(p_val = kruskal.results[["p.value"]], row.names = "kw")
#   
#   export <- dunn.results %>% 
#     bind_rows(kruskal.pval) %>% 
#     rownames_to_column(var = "comparison") %>% 
#     mutate(metric_id = {{metric}}, .before = comparison)
#   
#   return(export)
# }
# 
# wq_metrics <- colnames(wq_tb[7:10])
# wq_results_df <- map2_dfr(wq_metrics, "shoreline", compare_wq)

# site_wq <- wq_tb %>% 
#   group_by(month, year, site) %>% 
#   summarise_at(vars(secchi_depth_m:temperature), mean, na.rm = TRUE)
# 
# wq_results_site <- map2_dfr(wq_metrics, "site", compare_wq)

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
             cols = vars(year),
             scales = "free_y")


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

wq_tb_range_export <- wq_tb_long %>%
  group_by(site, metric) %>%
  summarize(mean = round(mean(value, na.rm = TRUE), 2), min = min(value, na.rm = TRUE), max = max(value, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(value = paste0(mean, " (", min, "-", max, ")")) %>%
  select(site, metric, value) %>%
  pivot_wider(names_from = metric, values_from = value)

# write_csv(wq_tb_export, file = here("data", "wq_tb.csv"))


# wq_tb %>% 
#   pivot_longer(c(secchi_depth_m, do_mg_l, salinity_ppm, temperature), names_to = "metric") %>% 
#   group_by(year, site, metric) %>% 
#   summarize(mean_value = mean(value, na.rm = TRUE), 
#             sd_value = sd(value, na.rm = TRUE)) %>% 
#   ggplot(aes(x = site, y = mean_value, fill = factor(year))) +
#   geom_bar(stat = 'identity', position = 'dodge') +
#   geom_errorbar(aes(ymin = mean_value - sd_value, 
#                     ymax = mean_value + sd_value),
#                 width=.2, position = position_dodge(width = 0.9)) +
#   facet_wrap(~factor(metric, labels = c("Dissolved oxygen (mg/L)",
#                                         "Salinity (ppm)",
#                                         "Secchi depth (m)",
#                                         "Temperature (C)")), scales = "free_y") +
#   theme_bw() + 
#   labs(y = "Mean value", x = "Site", fill = "Year")

# compare_secchi <- wq_tb %>% 
#   group_by(site, month, year) %>% 
#   summarize(mean_val = mean(secchi_depth_m)) %>% 
#   ungroup() %>% 
#   pivot_wider(names_from = site, values_from = mean_val) %>% 
#   relocate(TUR, .after = year) %>% 
#   mutate(across(TUR:EDG, ~  . - TUR)) %>% 
#   pivot_longer(!c(month, year), names_to = "site") %>% 
#   group_by(site) %>% 
#   summarize(mean_diff = mean(value, na.rm = TRUE))
#   
# compare_TUR_temp <- wq_tb %>% 
#   group_by(site, month, year) %>% 
#   summarize(mean_val = mean(temperature)) %>% 
#   ungroup() %>% 
#   pivot_wider(names_from = site, values_from = mean_val) %>% 
#   relocate(TUR, .after = year) %>% 
#   mutate(across(TUR:EDG, ~  . - TUR)) %>% 
#   pivot_longer(!c(month, year), names_to = "site") %>% 
#   group_by(site) %>% 
#   summarize(mean_diff = mean(value, na.rm = TRUE))
# 
# compare_COR_temp <- wq_tb %>% 
#   group_by(site, month, year) %>% 
#   summarize(mean_val = mean(temperature)) %>% 
#   ungroup() %>% 
#   pivot_wider(names_from = site, values_from = mean_val) %>% 
#   relocate(TUR, .after = year) %>% 
#   mutate(across(TUR:EDG, ~  . - COR)) %>% 
#   pivot_longer(!c(month, year), names_to = "site") %>% 
#   group_by(site) %>% 
#   summarize(mean_diff = mean(value, na.rm = TRUE))
#   
  
# wq_tb %>% 
#   pivot_longer(c(secchi_depth_m, do_mg_l, salinity_ppm, temperature), names_to = "metric") %>% 
#   group_by(site, shoreline, month, metric) %>% 
#   mutate(mean_value = mean(value)) %>% 
#   mutate(mean_value = ifelse(is.na(mean_value), value, mean_value)) %>% 
#   mutate(ipa = ifelse(ipa == "Natural2", "Natural", ipa)) %>% 
#   ggplot(aes(x = month, 
#              y = mean_value, 
#              shape = ipa, 
#              color = site)) +
#   geom_point(size = 2) + 
#   geom_line() + 
#   theme_bw() +
#   facet_wrap(~metric, scales = "free_y") +
#   labs(y = "Mean value", x = "Month", color = "Site", shape = "Shoreline\nCondition")
  

