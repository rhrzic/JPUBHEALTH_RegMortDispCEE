# Subnational mortality convergence in Europe by cause of death
# Rok Hrzic (r.hrzic (at) maastrichtuniversity.nl)
# Version September 2024

rm(list = ls())

require(tidyverse)
require(fuzzyjoin)
require(eurostat)

#all cause mortality

CZE_all_deaths_1990_2021 = read.csv(file = 'data/Czechia/CZE_Deaths_NUTS3_1990-2021.csv', sep = ';') %>%
  filter(Year >= 2000 & Age != 'TOT' & Sex != 'b') %>%
  mutate(CauseCategory = 'All',
         Age = as.integer(Age)) %>%
  select(Country, RegionCode, Year, Sex, AgeGroup = Age, CauseCategory, Deaths = Value) 

CZE_all_deaths_2022 = read.csv(file = 'data/Czechia/CZE_Deaths_NUTS3_2022.csv', sep = ';') %>%
  filter(Age != 'TOT' & Sex != 'b') %>%
  filter(Year >= 2000 & Age != 'TOT') %>%
  mutate(CauseCategory = 'All',
         Age = as.integer(Age)) %>%
  select(Country, RegionCode, Year, Sex, AgeGroup = Age, CauseCategory, Deaths = Value) 


CZE_all_deaths = rbind(CZE_all_deaths_1990_2021, CZE_all_deaths_2022)
rm(CZE_all_deaths_1990_2021, CZE_all_deaths_2022)


CZE_pop_1990_2021 = read.csv(file = 'data/Czechia/Czechia_Pop_NUTS3_1990-2021.csv', sep = ';') %>%
  filter(Year >= 2000 & Sex != 'b') %>%
  mutate(Age = as.integer(Age),
         Value = as.integer(Value)) %>%
  select(Country, RegionCode, Year, Sex, AgeGroup = Age, Pop = Value) 

CZE_pop_2022 = read.csv(file = 'data/Czechia/CZE_Mid-year_Pop_NUTS3_2022.csv', sep = ';') %>%
  filter(Age != 'TOT' & Sex != 'b') %>%
  mutate(Age = as.integer(Age),
         Value = str_replace(Value, ',', '.'),
         Value = as.integer(Value)) %>%
  select(Country, RegionCode, Year, Sex, AgeGroup = Age, Pop = Value) 

CZE_pop = rbind(CZE_pop_1990_2021, CZE_pop_2022)
rm(CZE_pop_1990_2021, CZE_pop_2022)

ESP <- read.csv(file = "data/european_standard_population.csv") %>%
  mutate(Age = as.numeric(str_sub(AgeGroup, start = 1, end = 2)),
         ESPw = EuropeanStandardPopulation/100000)

CZE_all_deaths = left_join(CZE_all_deaths, CZE_pop)
  # left_join(., ESP) %>%
  # ungroup() %>%
  # group_by(Country, RegionCode, Year, Sex, CauseCategory) %>%
  # summarise(cdr = 1000*sum(Deaths)/sum(Pop),
  #           sdr = 1000*sum(ESPw*Deaths/Pop),
  #           Pop = sum(Pop))

rm(ESP)

#cod mortality

CZE_cod_deaths <- read.csv(file = "data/Czechia/Czechia_CoD_Rok.csv", sep=";") %>%
  filter(Year >= 2000 & CauseCode != 'all causes') %>%
  mutate(Country = 'CZE') %>%
  select(Country, RegionCode, Year, Sex, AgeGroup = Age, CauseCode, Deaths = Value)


## Categorizion of cod


SIM_list_10 <- read.csv(file = "data/Short-intermediate list Rok.csv", sep=";") %>%
  filter(ICD == 'ICD10_3')

CZE_cod_deaths = CZE_cod_deaths %>%
  regex_left_join(., SIM_list_10, by = c('CauseCode' = 'Regex'))

## I summarise by cause category (not interested in the finer distinctions)

CZE_cod_deaths = CZE_cod_deaths %>%
  group_by(Country, RegionCode, Year, Sex, AgeGroup, CauseCategory) %>%
  summarise(Deaths = sum(Deaths))

# The lack of zero counts causes some age categories to lack the large cause categories, esp CVD and Cancer.

all.causes = expand_grid(Country = 'CZE',
  RegionCode = unique(CZE_cod_deaths$RegionCode), 
            Year = unique(CZE_cod_deaths$Year),
            Sex = unique(CZE_cod_deaths$Sex),
            AgeGroup = unique(CZE_cod_deaths$AgeGroup),
            CauseCategory = unique(CZE_cod_deaths$CauseCategory))


CZE_cod_deaths = CZE_cod_deaths %>%
  left_join(all.causes, .) %>% 
  mutate(Deaths = ifelse(is.na(Deaths), 0, Deaths))


CZE_cod_deaths = CZE_cod_deaths %>%
  mutate(RegionCode = str_sub(RegionCode, 1, -2L)) %>%
  group_by(Country, RegionCode, Year, Sex, AgeGroup, CauseCategory) %>%
  summarise(Deaths = sum(Deaths))


CZE_pop = CZE_pop %>%
  mutate(AgeGroup = as.numeric(str_sub(AgeGroup, start = 1, end = 2)),
         AgeGroup = case_when(AgeGroup == 0 ~ 0,
                              AgeGroup %in% c(1, 5, 10, 15, 20) ~ 1,
                              AgeGroup %in% c(25, 30, 35, 40, 45) ~ 25,
                              AgeGroup %in% c(50, 55, 60, 65) ~ 50,
                              AgeGroup %in% c(70, 75, 80) ~ 70,
                              AgeGroup %in% c(85, 90, 95) ~ 85)) %>%
  group_by(Country, RegionCode, Year, Sex, AgeGroup) %>%
  summarise(Pop = sum(Pop))

# ESP <- read.csv(file = "data/european_standard_population.csv") %>%
#   mutate(AgeGroup = as.numeric(str_sub(AgeGroup, start = 1, end = 2)),
#          Age = case_when(AgeGroup == 0 ~ 0,
#                               AgeGroup %in% c(1, 5, 10, 15, 20) ~ 1,
#                               AgeGroup %in% c(25, 30, 35, 40, 45) ~ 25,
#                               AgeGroup %in% c(50, 55, 60, 65) ~ 50,
#                               AgeGroup %in% c(70, 75, 80) ~ 70,
#                               AgeGroup %in% c(85, 90, 95) ~ 85)) %>%
#   group_by(Age) %>%
#   summarise(ESPw = sum(EuropeanStandardPopulation)/100000)


CZE_cod_deaths = left_join(CZE_cod_deaths, CZE_pop)
  # left_join(., ESP) %>%
  # ungroup() %>%
  # group_by(Country, RegionCode, Year, Sex, CauseCategory) %>%
  # summarise(cdr = 1000*sum(Deaths)/sum(Pop),
  #           sdr = 1000*sum(ESPw*Deaths/Pop),
  #           Pop = sum(Pop))


# ggplot(CZE_cod_deaths, aes(x = Year, y = as.factor(AgeGroup), fill = log(Deaths/Pop)))+
#   geom_tile()+
#   facet_grid(CauseCategory ~ Sex)+
#   scale_fill_viridis_c()


CZE_deaths = rbind(CZE_cod_deaths, CZE_all_deaths) %>%
  arrange(RegionCode, Year, Sex, CauseCategory, AgeGroup)

saveRDS(CZE_deaths, 'temp/CZE_for_smooth.rds')


# ggplot(CZE_deaths, aes(x = Year, y = sdr, color = CauseCategory, group = interaction(RegionCode, CauseCategory, Sex)))+
#          geom_line()+
#          facet_grid(. ~ Sex)
