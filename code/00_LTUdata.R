# Subnational mortality convergence in Europe by cause of death
# Rok Hrzic (r.hrzic (at) maastrichtuniversity.nl)
# April 2025

rm(list = ls())

require(tidyverse)
require(fuzzyjoin)

#all cause mortality

LTU_deaths = read.csv(file = 'data/Lithuania/Lithuania_CoD_Rok.csv', sep = ';') %>%
  filter(Year >= 2010 & Age != 'TOT' & Sex != 'b') %>%
  mutate(CauseCategory = case_when(CauseCode == 1 ~ "Cancer",
                                   CauseCode == 2 ~ "Cardiovascular",
                                   CauseCode == 3 ~ "External",
                                   CauseCode == 4 ~ "Other"),
         Age = as.integer(Age),
         Country = 'LTU') %>%
  select(Country, RegionCode, Year, Sex, AgeGroup = Age, CauseCategory, Deaths = Value) %>%
  group_by(Country, RegionCode, Year, Sex, AgeGroup, CauseCategory) %>%
  summarise(Deaths = sum(Deaths))

## No zero counts, assuming zeros are missing, thus require expanding the grid

grid <- expand_grid(
  Country = unique(LTU_deaths$Country),
  RegionCode = unique(LTU_deaths$RegionCode),
  Year = unique(LTU_deaths$Year),
  Sex = unique(LTU_deaths$Sex),
  AgeGroup = unique(LTU_deaths$AgeGroup),
  CauseCategory = unique(LTU_deaths$CauseCategory)
)

LTU_cod_deaths = grid %>%
  left_join(LTU_deaths, by = c("Country", "RegionCode", "Year", "Sex", "AgeGroup", "CauseCategory")) %>%
  mutate(Deaths = replace_na(Deaths, 0))  %>%
  arrange(Country, RegionCode, Year, Sex, CauseCategory, AgeGroup)


LTU_all_deaths = LTU_cod_deaths %>%
  group_by(Country, RegionCode, Year, Sex, AgeGroup) %>%
  summarise(CauseCategory = 'All',
            Deaths = sum(Deaths))

rm(LTU_deaths)

LTU_pop = read.csv(file = 'data/Lithuania/Lithuania_Pop_Rok.csv', sep = ';') %>%
  filter(Year >= 2000 & Sex != 'b') %>%
  mutate(Age = as.integer(Age),
         Value = as.integer(Value)) %>%
  select(Country, RegionCode, Year, Sex, AgeGroup = Age, Pop = Value) 


# ESP <- read.csv(file = "data/european_standard_population.csv") %>%
#   mutate(AgeGroup = as.numeric(str_sub(AgeGroup, start = 1, end = 2)),
#          AgeGroup = case_when(AgeGroup == 0 ~ 0,
#                               AgeGroup %in% c(1, 5, 10, 15, 20) ~ 1,
#                               AgeGroup %in% c(25, 30, 35, 40, 45) ~ 25,
#                               AgeGroup %in% c(50, 55, 60, 65) ~ 50,
#                               AgeGroup %in% c(70, 75, 80) ~ 70,
#                               AgeGroup %in% c(85, 90, 95) ~ 85)) %>%
#   group_by(AgeGroup) %>%
#   summarise(ESPw = sum(EuropeanStandardPopulation)/100000)


LTU_all_deaths = left_join(LTU_all_deaths, LTU_pop)

  # left_join(., ESP) %>%
  # ungroup() %>%
  # group_by(Country, RegionCode, Year, Sex, CauseCategory) %>%
  # arrange(AgeGroup) %>%
  # summarise(cdr = 1000*sum(Deaths)/sum(Pop),
  #           sdr = 1000*sum(ESPw*Deaths/Pop),
  #           Pop = sum(Pop))

LTU_cod_deaths = left_join(LTU_cod_deaths, LTU_pop)
  # left_join(., ESP) %>%
  # ungroup() %>%
  # group_by(Country, RegionCode, Year, Sex, CauseCategory) %>%
  # summarise(cdr = 1000*sum(Deaths)/sum(Pop),
  #           sdr = 1000*sum(ESPw*Deaths/Pop),
  #           Pop = sum(Pop))


LTU_deaths = rbind(LTU_cod_deaths, LTU_all_deaths) %>%
  arrange(RegionCode, Year, Sex, CauseCategory, AgeGroup)

saveRDS(LTU_deaths, 'temp/LTU_LAU_for_smooth.rds')

# ggplot(LTU_deaths, aes(x = Year, y = sdr, color = CauseCategory, group = interaction(RegionCode, CauseCategory, Sex)))+
#   geom_line(alpha = 0.5)+
#   facet_grid(. ~ Sex)

## NUTS3

LTU_LAU_NUTS3 = readxl::read_excel('data/Lithuania/LTU.LAU1.xlsx', sheet = 2)

LTU_deaths = read.csv(file = 'data/Lithuania/Lithuania_CoD_Rok.csv', sep = ';') %>%
  filter(Year >= 2010 & Age != 'TOT' & Sex != 'b') %>%
  mutate(CauseCategory = case_when(CauseCode == 1 ~ "Cancer",
                                   CauseCode == 2 ~ "Cardiovascular",
                                   CauseCode == 3 ~ "External",
                                   CauseCode == 4 ~ "Other"),
         Age = as.integer(Age),
         Country = 'LTU') %>%
  select(Country, RegionCode, Year, Sex, AgeGroup = Age, CauseCategory, Deaths = Value) %>%
  left_join(LTU_LAU_NUTS3, by = c('RegionCode' = 'Code LAU')) %>%
  select(-RegionCode) %>%
  rename(RegionCode=`Code NUTS-3`) %>%
  group_by(Country, RegionCode, Year, Sex, AgeGroup, CauseCategory) %>%
  summarise(Deaths = sum(Deaths))

## No zero counts, assuming zeros are missing, thus require expanding the grid

grid <- expand_grid(
  Country = unique(LTU_deaths$Country),
  RegionCode = unique(LTU_deaths$RegionCode),
  Year = unique(LTU_deaths$Year),
  Sex = unique(LTU_deaths$Sex),
  AgeGroup = unique(LTU_deaths$AgeGroup),
  CauseCategory = unique(LTU_deaths$CauseCategory)
)

LTU_cod_deaths = grid %>%
  left_join(LTU_deaths, by = c("Country", "RegionCode", "Year", "Sex", "AgeGroup", "CauseCategory")) %>%
  mutate(Deaths = replace_na(Deaths, 0))  %>%
  arrange(Country, RegionCode, Year, Sex, CauseCategory, AgeGroup)


LTU_all_deaths = LTU_cod_deaths %>%
  group_by(Country, RegionCode, Year, Sex, AgeGroup) %>%
  summarise(CauseCategory = 'All',
            Deaths = sum(Deaths))

rm(LTU_deaths)

LTU_pop = read.csv(file = 'data/Lithuania/Lithuania_Pop_Rok.csv', sep = ';') %>%
  filter(Year >= 2000 & Sex != 'b') %>%
  mutate(Age = as.integer(Age),
         Value = as.integer(Value)) %>%
  select(Country, RegionCode, Year, Sex, AgeGroup = Age, Pop = Value) %>%
  left_join(LTU_LAU_NUTS3, by = c('RegionCode' = 'Code LAU')) %>%
  select(-RegionCode) %>%
  rename(RegionCode=`Code NUTS-3`) %>%
  group_by(Country, RegionCode, Year, Sex, AgeGroup) %>%
  summarise(Pop = sum(Pop))


LTU_all_deaths = left_join(LTU_all_deaths, LTU_pop)
  # left_join(., ESP) %>%
  # ungroup() %>%
  # group_by(Country, RegionCode, Year, Sex, CauseCategory) %>%
  # arrange(AgeGroup) %>%
  # summarise(cdr = 1000*sum(Deaths)/sum(Pop),
  #           sdr = 1000*sum(ESPw*Deaths/Pop),
  #           Pop = sum(Pop))

LTU_cod_deaths = left_join(LTU_cod_deaths, LTU_pop)
  # left_join(., ESP) %>%
  # ungroup() %>%
  # group_by(Country, RegionCode, Year, Sex, CauseCategory) %>%
  # summarise(cdr = 1000*sum(Deaths)/sum(Pop),
  #           sdr = 1000*sum(ESPw*Deaths/Pop),
  #           Pop = sum(Pop))


LTU_deaths = rbind(LTU_cod_deaths, LTU_all_deaths) %>%
  arrange(RegionCode, Year, Sex, CauseCategory, AgeGroup)

saveRDS(LTU_deaths, 'temp/LTU_for_smooth.rds')


