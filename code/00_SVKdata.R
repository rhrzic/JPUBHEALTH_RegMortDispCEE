# Subnational mortality convergence in Europe by cause of death
# Rok Hrzic (r.hrzic (at) maastrichtuniversity.nl)
# Version March 2025

rm(list = ls())

require(tidyverse)
require(fuzzyjoin)
require(eurostat)

#all cause mortality

SVK_deaths = read.csv(file = 'data/Slovakia/Slovakia_CoD_Rok.csv', sep = ';') %>%
  filter(Year >= 2000 & Age != 'TOT' & Sex != 'b') %>%
  mutate(CauseCategory = case_when(CauseCode == 1 ~ "Cancer",
                                   CauseCode == 2 ~ "Cardiovascular",
                                   CauseCode == 3 ~ "External",
                                   CauseCode == 4 ~ "Other"),
         Age = as.integer(Age),
         Country = 'SVK') %>%
  select(Country, RegionCode, Year, Sex, AgeGroup = Age, CauseCategory, Deaths = Value) 

## No zero counts, assuming zeros are missing, thus require expanding the grid

grid <- expand_grid(
  Country = unique(SVK_deaths$Country),
  RegionCode = unique(SVK_deaths$RegionCode),
  Year = unique(SVK_deaths$Year),
  Sex = unique(SVK_deaths$Sex),
  AgeGroup = unique(SVK_deaths$AgeGroup),
  CauseCategory = unique(SVK_deaths$CauseCategory)
)

SVK_cod_deaths = grid %>%
  left_join(SVK_deaths, by = c("Country", "RegionCode", "Year", "Sex", "AgeGroup", "CauseCategory")) %>%
  mutate(Deaths = replace_na(Deaths, 0))


SVK_all_deaths = SVK_cod_deaths %>%
  group_by(Country, RegionCode, Year, Sex, AgeGroup) %>%
  summarise(CauseCategory = 'All',
            Deaths = sum(Deaths))

rm(SVK_deaths)

SVK_pop = read.csv(file = 'data/Slovakia/Slovakia_Pop_Rok.csv', sep = ';') %>%
  filter(Year >= 2000 & Sex != 'b') %>%
  mutate(Age = as.integer(Age),
         Value = as.integer(Value)) %>%
  select(Country, RegionCode, Year, Sex, AgeGroup = Age, Pop = Value) 


# ESP <- read.csv(file = "data/european_standard_population.csv") %>%
#   mutate(AgeGroup = as.numeric(str_sub(AgeGroup, start = 1, end = 2)),
#          AgeGroup = case_when(AgeGroup %in% c(0, 1, 5) ~ 0,
#                               AgeGroup %in% c(10, 15) ~ 10,
#                               AgeGroup == 20 ~ 20,
#                               AgeGroup == 25 ~ 25,
#                               AgeGroup == 30 ~ 30,
#                               AgeGroup == 35 ~ 35,
#                               AgeGroup == 40 ~ 40,
#                               AgeGroup == 45 ~ 45,
#                               AgeGroup == 50 ~ 50,
#                               AgeGroup == 55 ~ 55,
#                               AgeGroup == 60 ~ 60,
#                               AgeGroup == 65 ~ 65,
#                               AgeGroup == 70 ~ 70,
#                               AgeGroup == 75 ~ 75,
#                               AgeGroup == 80 ~ 80,
#                               AgeGroup %in% c(85, 90, 95) ~ 85)) %>%
#   group_by(AgeGroup) %>%
#   summarise(ESPw = sum(EuropeanStandardPopulation)/100000)


SVK_all_deaths = left_join(SVK_all_deaths, SVK_pop)
  # left_join(., ESP) %>%
  # ungroup() %>%
  # group_by(Country, RegionCode, Year, Sex, CauseCategory) %>%
  # summarise(cdr = 1000*sum(Deaths)/sum(Pop),
  #           sdr = 1000*sum(ESPw*Deaths/Pop),
  #           Pop = sum(Pop))

SVK_cod_deaths = left_join(SVK_cod_deaths, SVK_pop)
  # left_join(., ESP) %>%
  # ungroup() %>%
  # group_by(Country, RegionCode, Year, Sex, CauseCategory) %>%
  # summarise(cdr = 1000*sum(Deaths)/sum(Pop),
  #           sdr = 1000*sum(ESPw*Deaths/Pop),
  #           Pop = sum(Pop))


SVK_deaths = rbind(SVK_cod_deaths, SVK_all_deaths) %>%
  arrange(RegionCode, Year, Sex, CauseCategory, AgeGroup)

saveRDS(SVK_deaths, 'temp/SVK_LAU_for_smooth.rds')

# ggplot(SVK_deaths, aes(x = Year, y = sdr, color = CauseCategory, group = interaction(RegionCode, CauseCategory, Sex)))+
#   geom_line(alpha = 0.5)+
#   facet_grid(. ~ Sex)

## NUTS3

#all cause mortality

SVK_deaths = read.csv(file = 'data/Slovakia/Slovakia_CoD_Rok.csv', sep = ';') %>%
  filter(Year >= 2000 & Age != 'TOT' & Sex != 'b') %>%
  mutate(CauseCategory = case_when(CauseCode == 1 ~ "Cancer",
                                   CauseCode == 2 ~ "Cardiovascular",
                                   CauseCode == 3 ~ "External",
                                   CauseCode == 4 ~ "Other"),
         Age = as.integer(Age),
         Country = 'SVK',
         RegionCode = str_sub(RegionCode, start = 1L, end = 3L)) %>%
  select(Country, RegionCode, Year, Sex, AgeGroup = Age, CauseCategory, Deaths = Value) %>%
  group_by(Country, RegionCode, Year, Sex, AgeGroup, CauseCategory) %>%
  summarise(Deaths = sum(Deaths))

## No zero counts, assuming zeros are missing, thus require expanding the grid

grid <- expand_grid(
  Country = unique(SVK_deaths$Country),
  RegionCode = unique(SVK_deaths$RegionCode),
  Year = unique(SVK_deaths$Year),
  Sex = unique(SVK_deaths$Sex),
  AgeGroup = unique(SVK_deaths$AgeGroup),
  CauseCategory = unique(SVK_deaths$CauseCategory)
)

SVK_cod_deaths = grid %>%
  left_join(SVK_deaths, by = c("Country", "RegionCode", "Year", "Sex", "AgeGroup", "CauseCategory")) %>%
  mutate(Deaths = replace_na(Deaths, 0))


SVK_all_deaths = SVK_cod_deaths %>%
  group_by(Country, RegionCode, Year, Sex, AgeGroup) %>%
  summarise(CauseCategory = 'All',
            Deaths = sum(Deaths))

rm(SVK_deaths)

SVK_pop = read.csv(file = 'data/Slovakia/Slovakia_Pop_Rok.csv', sep = ';') %>%
  filter(Year >= 2000 & Sex != 'b') %>%
  mutate(Age = as.integer(Age),
         Value = as.integer(Value),
         RegionCode = str_sub(RegionCode, start = 1L, end = 3L)) %>%
  select(Country, RegionCode, Year, Sex, AgeGroup = Age, Pop = Value) %>%
  group_by(Country, RegionCode, Year, Sex, AgeGroup) %>%
  summarise(Pop = sum(Pop))


SVK_all_deaths = left_join(SVK_all_deaths, SVK_pop)
  # left_join(., ESP) %>%
  # ungroup() %>%
  # group_by(Country, RegionCode, Year, Sex, CauseCategory) %>%
  # summarise(cdr = 1000*sum(Deaths)/sum(Pop),
  #           sdr = 1000*sum(ESPw*Deaths/Pop),
  #           Pop = sum(Pop))

SVK_cod_deaths = left_join(SVK_cod_deaths, SVK_pop)
  # left_join(., ESP) %>%
  # ungroup() %>%
  # group_by(Country, RegionCode, Year, Sex, CauseCategory) %>%
  # summarise(cdr = 1000*sum(Deaths)/sum(Pop),
  #           sdr = 1000*sum(ESPw*Deaths/Pop),
  #           Pop = sum(Pop))

SVK_deaths = rbind(SVK_cod_deaths, SVK_all_deaths) %>%
  arrange(RegionCode, Year, Sex, CauseCategory, AgeGroup)

saveRDS(SVK_deaths, 'temp/SVK_for_smooth.rds')

# ggplot(SVK_deaths, aes(x = Year, y = sdr, color = CauseCategory, group = interaction(RegionCode, CauseCategory, Sex)))+
#   geom_line(alpha = 0.5)+
#   facet_grid(. ~ Sex)
