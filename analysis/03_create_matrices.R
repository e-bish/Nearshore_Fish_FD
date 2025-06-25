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
fish_L_full <- net_core %>% #L is referring to the RLQ analysis
  filter(!is.na(ComName)) %>% 
  group_by(site, ipa, year, ComName) %>% 
  summarize(spp_sum = sum(species_count)) %>% #sum across depth stations and months within each shoreline type
  ungroup() %>%
  full_join(sampling_events) %>% 
  mutate(catch_per_set = spp_sum/no_net_sets) %>% #on a couple of occasions we did not sample at all three depths
  filter(!is.na(ComName)) %>% #these are accounted for in the next step
  full_join(expand_species) %>% #add back in all of the events so we can capture 0s
  select(!c(spp_sum, no_net_sets)) %>% 
  mutate(ComName = replace(ComName, ComName == "Pacific Sandfish", "Pacific sandfish")) %>% #if it's capitolized it gets confused about the proper order
  arrange(site, ipa, year, ComName) %>% 
  mutate(catch_per_set  = replace_na(catch_per_set , 0)) %>% 
  pivot_wider(names_from = ComName, values_from = catch_per_set, values_fill = 0) %>% 
  clean_names() %>%
  ungroup()

save(fish_L_full, file = here("data", "fish_L_full.Rdata")) #last saved 6/3/25

fish_L_sample <- fish_L_full %>% 
  mutate(sample = paste(site, ipa, year, sep = "_"), .after = year) %>% 
  select(!1:3) %>% 
  column_to_rownames(var = "sample")

save(fish_L_sample, file = here("data", "fish_L_sample.Rdata")) #last saved 6/3/25

#check the abundance matrix for samples with few species
rows_w_few_spp <- fish_L_sample %>%
  decostand(method = "pa") %>%
  filter(rowSums(.) < 4)

save(rows_w_few_spp, file = here("data", "rows_w_few_spp.Rdata")) #last saved 6/3/25

samples_w_few_spp <- rownames(rows_w_few_spp)

#remove samples that have less than 4 (this number depends on the number of axes we want to retain for the trait space)
fish_L <- fish_L_sample %>%
  filter(!rownames(.) %in% samples_w_few_spp)

#check that all species are represented in at least one sample
colSums(fish_L)

## create the trait matrix 

#extract the names of the species present in the lampara net dataset
spp_names <- net_core %>% 
  distinct(ComName) %>% 
  mutate(Species = NA) 

#link these common names to their scientific names in fishbase
#note that since updating R this takes a very long time to run

sci_names <- vector(mode = 'list', length = length(spp_names))

for (i in 1:nrow(spp_names)) {
  sci_names[[i]] <- rfishbase::common_to_sci(spp_names[i,])
  spp_names[i,2] <- ifelse(nrow(sci_names[[i]]) == 1, sci_names[[i]][[1]], NA)
}

#specify which entry to use for species that have multiple common names
spp_names <- spp_names %>% 
  mutate(Species = case_when(ComName == "Pacific herring" ~ "Clupea pallasii", 
                             ComName == "Pacific Cod" ~ "Gadus macrocephalus",
                             ComName == "Northern Anchovy" ~ "Engraulis mordax",
                             ComName == "Striped Seaperch" ~ "Embiotoca lateralis",
                             ComName == "Rock Sole" ~ "Lepidopsetta bilineata",
                             ComName == "Steelhead trout" ~ "Oncorhynchus mykiss",
                             ComName == "Cutthroat trout" ~ "Oncorhynchus clarkii",
                             ComName == "Tube-snout" ~ "Aulorhynchus flavidus",
                             TRUE ~ Species)) %>% 
  arrange(ComName) %>% 
  filter(!str_detect(ComName, 'UnID'))

write_csv(spp_names, "data/spp_names.csv")
spp_names <- read_csv(here("data", "spp_names.csv"))

#calculate the mean fork length of each species from the subsamples taken in the field
fork_length <- net_core %>% 
  group_by(ComName) %>% 
  summarize(mean_length_mm = mean(mean_length_mm)) %>% 
  inner_join(spp_names) %>% 
  mutate(mean_length_mm = ifelse(ComName == "Tidepool Sculpin", 89.0, mean_length_mm)) #no length in our df so taking the max length from fishbase

#extract other trait data from rfishbase
milieu <- species(spp_names$Species) %>% 
  mutate(DemersPelag = ifelse(Genus == "Oncorhynchus", "epipelagic", DemersPelag)) %>% #juvenile salmon are generally near the surface, 
  #as reviewed in Quinn 2018 and observed by Munsch et al. 2017, and for steelhead cite Daly et al. 2014
  select(Species, BodyShapeI, DemersPelag, AnaCat) %>% 
  mutate(DemersPelag = ifelse(DemersPelag == "pelagic-neritic", "pelagic", DemersPelag)) %>% #simplify this category because there is only one pelagic and two pelagic-neritic species
  rename(migrations = "AnaCat") %>% 
  mutate(migrations = ifelse(migrations == "non-migratory" | is.na(migrations), "resident", "migratory")) %>% 
  mutate(migrations = factor(migrations))

feeding_guild1 <- fooditems(spp_names$Species) 

save(feeding_guild1, file = here("data", "feeding_guild1.Rda")) #save current version in case these hard coded values change later
load(here("data", "feeding_guild1.Rda"))

manual_check <- feeding_guild1 %>% 
  select(Species, FoodI, CommonessII, PredatorStage, Locality2) %>% 
  slice(11, 56, 68, 97, 
        130, 161, 209, 458, #recruits,
        485, 594, 656, 660, 
        696, 780, 817, 997, 
        1006, 1032, 1034,
        1057, 1077, 1100)

omnivores <- feeding_guild1 %>% 
  filter(!Species %in% manual_check$Species) %>% 
  distinct(Species) %>% 
  mutate(FoodI = "omnivorous")

combine_fg <- bind_rows(manual_check, omnivores) %>% 
  mutate(feeding_guild = case_when(FoodI == "zoobenthos" ~ "zoobenthivorous",
                                   FoodI == "zooplankton" ~ "planktivorous", 
                                   TRUE ~ FoodI)) %>% 
  full_join(spp_names) %>% 
  mutate(feeding_guild = case_when(Species == "Blepsias cirrhosus" ~ "zoobenthivorous", #fishbase "Diet"
                                   Species == "Liparis florae" ~ "zoobenthivorous", #based on Liparis pulchellus
                                   TRUE ~ feeding_guild)) %>% 
  select(ComName, Species, feeding_guild) %>% 
  arrange(ComName)

#chinook are omnivorous - duffy et al. 2010
#chinook and chum in eelgrass mostly eat epifaunal invertebrates - Kennedy et al. 2018
#coho and chinook eat some fish, otherwise they (and other species) eat a combination of zoobenthos and zooplankton - Beamish et al. 2003 

fish_traits <- full_join(fork_length, milieu) %>% 
  select(3,1,2,4,5,6) %>% 
  left_join(combine_fg) %>% 
  mutate(ComName = replace(ComName, ComName == "Pacific Sandfish", "Pacific sandfish")) %>%
  arrange(ComName)

write_csv(fish_traits, "data/fish_traits.csv") #last saved 6/25/25

fish_Q <- fish_traits %>% 
  select(-c(Species, ComName)) %>% 
  mutate_if(is.character, as.factor) %>% 
  clean_names() %>% 
  as.data.frame()

rownames(fish_Q) <- colnames(fish_L)

confirm_proper_names <- fish_Q %>% 
  mutate(Species = fish_traits$ComName, .before = mean_length_mm)

# #check the coefficient of variation
CV <- function(x) { 100 * sd(x) / mean(x) }
CV(fish_Q$mean_length_mm)
# # <50 is small, so we don't necessarily need to transform the length data

#save final matrices
fish_meta <- fish_L_sample %>% 
  as_tibble(rownames = "sample") %>% 
  separate_wider_delim(sample, delim = "_", names = c("site", "ipa", "year"), cols_remove = TRUE) %>% 
  mutate(site = factor(site, levels = c("FAM", "TUR", "COR", "SHR", "DOK", "EDG")),
         region = ifelse(site %in% c("FAM", "TUR", "COR"), "North", "South"), 
         veg = ifelse(site %in% c("TUR", "COR", "SHR"), "present", "absent"), .after = site,
         shoreline = paste0(site, case_when(ipa == "Armored" ~ "A", 
                                            ipa == "Restored" ~ "R", 
                                            ipa == "Natural2" ~ "N2",
                                            TRUE ~ "N")),
         across(where(is.numeric), ~replace_na(., 0))) %>% 
  mutate(ipa = ifelse(ipa == "Natural2", "Natural", ipa)) %>% 
  arrange(site,ipa, shoreline, year) %>% 
  select(site:year)

fish.list <- list("trait" = fish_Q, 
                  "abund" = fish_L,
                  "meta" = fish_meta) 

save(fish.list, file = here("data", "fish.list.Rdata")) #last saved 6/25/25
