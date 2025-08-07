# Subnational mortality convergence in Europe by cause of death
# Rok Hrzic (r.hrzic (at) maastrichtuniversity.nl)
# April 2025

rm(list = ls())

require(tidyverse)

##all causes

POL_all_deaths = read.csv(file = 'data/Poland/POL_Deaths_NUTS3_1999-2023.csv', sep = ';') %>%
  filter(Year >= 2002, Age != 'TOT', Sex != 'b') %>%
  select(Country, RegionCode, Year, Sex, Age, Deaths = Value) %>%
  mutate(CauseCategory = 'All',
         AgeGroup = as.integer(Age),
         AgeGroup = ifelse(AgeGroup >= 85, 85, AgeGroup),
         AgeGroup = ifelse(AgeGroup %in% 0:1, 0, AgeGroup)) %>%
  group_by(Country, RegionCode, Year, Sex, AgeGroup, CauseCategory) %>%
  summarise(Deaths = sum(Deaths))


POL_pop = read.csv(file = 'data/Poland/POL_Pop_NUTS3_2002-2023.csv', sep = ';') %>%
  filter(Year >= 2000, Age != 'TOT', Sex != 'b') %>%
  select(Country, RegionCode, Year, Sex, Age, Pop = Value) %>%
  mutate(CauseCategory = 'All',
         AgeGroup = as.integer(Age), 
         AgeGroup = ifelse(Year %in% 2002:2005 & AgeGroup >= 65, 65, AgeGroup))  %>%
  group_by(Country, RegionCode, Year, Sex, AgeGroup) %>%
  summarise(Pop = sum(Pop))

# 
# ESP_65 <- read.csv(file = "data/european_standard_population.csv") %>%
#   mutate(AgeGroup = as.numeric(str_sub(AgeGroup, start = 1, end = 2)),
#          AgeGroup = ifelse(AgeGroup >= 65, 65, AgeGroup),
#          AgeGroup = ifelse(AgeGroup %in% 0:1, 0, AgeGroup)) %>%
#   group_by(AgeGroup) %>%
#   summarise(ESPw = sum(EuropeanStandardPopulation)/100000)
# 
# ESP_65 = expand_grid(Year = 2002:2005,
#                      AgeGroup = unique(ESP_65$AgeGroup)) %>%
#   left_join(., ESP_65)
# 
# ESP_85 <- read.csv(file = "data/european_standard_population.csv") %>%
#   mutate(AgeGroup = as.numeric(str_sub(AgeGroup, start = 1, end = 2)),
#          AgeGroup = ifelse(AgeGroup >= 85, 85, AgeGroup),
#          AgeGroup = ifelse(AgeGroup %in% 0:1, 0, AgeGroup)) %>%
#   group_by(AgeGroup) %>%
#   summarise(ESPw = sum(EuropeanStandardPopulation)/100000)
# 
# ESP_85 = expand_grid(Year = 2006:2023,
#             AgeGroup = unique(ESP_85$AgeGroup)) %>%
#   left_join(., ESP_85)
# 
# ESP = rbind(ESP_65, ESP_85)

POL_all_deaths = left_join(POL_all_deaths, POL_pop)
  # left_join(., ESP) %>%
  # ungroup() %>%
  # group_by(Country, RegionCode, Year, Sex, CauseCategory) %>%
  # summarise(cdr = 1000*sum(Deaths)/sum(Pop),
  #           sdr = 1000*sum(ESPw*Deaths/Pop),
  #           Pop = sum(Pop)) %>%
  # filter(Year %in% 2006:2022)


## cod

## Region code 0263 (Wałbrzych) was removed in 2003 and reintroduced as 0265 in 2013

POL_cod_deaths <- read.csv(file = "data/Poland/Poland_CoD_Rok.csv", sep=";") %>%
  filter(Sex != 'b' & Age != 'TOT' & Year %in% 2000:2019) %>%
  select(Country, RegionCode, Year, Sex, AgeGroup = Age, CauseCategory = CauseCode, Deaths = Value) %>%
  mutate(AgeGroup = as.numeric(AgeGroup),
         CauseCategory = case_when(CauseCategory == 1 ~ 'Other',
                                   CauseCategory == 2 ~ 'Cancer',
                                   CauseCategory == 3 ~ 'Cardiovascular',
                                   CauseCategory == 4 ~ 'External',
                                   CauseCategory == 999 ~ 'Unknown'),
         RegionCode = as.character(RegionCode),
         RegionCode = str_pad(RegionCode, 4, side = 'left', '0'),
         RegionCode = ifelse(RegionCode == "0263", "0265", RegionCode))


## proportional redistribution of missings (category 999)

POL_redistribution_filter = POL_cod_deaths %>%
  filter(CauseCategory == 'Unknown') %>%
  select(Country, RegionCode, Year, Sex, AgeGroup)

POL_to_redistribute = POL_cod_deaths %>%
  left_join(POL_redistribution_filter, .) %>%
  pivot_wider(names_from = CauseCategory, values_from = Deaths) %>%
  replace(is.na(.), 0) %>%
  mutate(Total = Unknown + Other + Cancer +External+Cardiovascular,
         Other = ifelse(Total != Unknown, Other + Unknown*Other/(Total-Unknown), Unknown/4),
         Cancer = ifelse(Total != Unknown, Cancer + Unknown*Cancer/(Total-Unknown), Unknown/4),
         External = ifelse(Total != Unknown,External + Unknown*External/(Total-Unknown), Unknown/4),
         Cardiovascular = ifelse(Total != Unknown,Cardiovascular + Unknown*Cardiovascular/(Total-Unknown), Unknown/4)) %>%
  select(-c(Total, Unknown)) %>%
  pivot_longer(cols = c(Other, Cancer, External, Cardiovascular), names_to = 'CauseCategory', values_to = 'Deaths')

##needs to return POL_deaths with 262133 - 390 = 261743 rows

POL_cod_deaths = POL_cod_deaths %>%
  anti_join(., POL_redistribution_filter)

POL_cod_deaths = rbind(POL_cod_deaths, POL_to_redistribute)


all.causes = expand_grid(Country = 'POL',
                         RegionCode = unique(POL_cod_deaths$RegionCode), 
                         Year = unique(POL_cod_deaths$Year),
                         Sex = unique(POL_cod_deaths$Sex),
                         AgeGroup = unique(POL_cod_deaths$AgeGroup),
                         CauseCategory = unique(POL_cod_deaths$CauseCategory))

POL_cod_deaths = POL_cod_deaths %>%
  right_join(all.causes) %>% 
  mutate(Deaths = ifelse(is.na(Deaths), 0, Deaths))

POL_cod_deaths = POL_cod_deaths %>%
  mutate(AgeGroup = ifelse(AgeGroup %in% 0:1, 0, AgeGroup)) %>%
  group_by(Country, RegionCode, Year, Sex, AgeGroup, CauseCategory) %>%
  summarise(Deaths = sum(Deaths))


## Aggregate to NUTS3


POL_GUS_NUTS3 = readxl::read_xlsx('data/Poland/POL_regions_2019.xlsx') %>%
  select(NUTS3 = `NUTS 3 CODE`, LAU = `LAU CODE`) %>%
  mutate(GUS = str_c(str_sub(LAU,5,6), str_sub(LAU,-4,-3))) %>%
  distinct(NUTS3, GUS)

## This is missing region 0263, which seems to exist in the other dataset. This GUS no longer exists. I need to find when and where the change happened.

#old NUTS3 classification

# POL_regions = readxl::read_xlsx('data/Poland/Poland_regions_2017.xlsx') %>%
#   rename_all(list(~make.names(.))) %>%
#   filter_all(any_vars(complete.cases(.))) %>%
#   fill(1:10, .direction = 'down') %>%
#   mutate(NTS.5..Name = ifelse(!is.na(NTS.5..Name_supplement), str_c(NTS.5..Name, NTS.5..Name_supplement, sep = " "), NTS.5..Name)) %>%
#   select(-NTS.5..Name_supplement) %>%
#   mutate_at(vars(matches("Code")), ~str_remove_all(.x, '\\.')) %>%
#   mutate_at(vars(matches("Code")), ~str_sub(.x, start = 2, end = -1)) %>%
#   mutate(GUS.code = as.character(str_c(str_sub(NTS.2..Code, -2L, -1L), str_sub(NTS.4...Code, -2L, -1L))),
#          LAU2.code = str_sub(NTS.5..Code, 1L, -2L)) %>%
#   select(NTS.3...Name, NTS.3...Code, NTS.4...Name, NTS.4...Code, NTS.5..Name, NTS.5..Code, GUS.code, LAU2.code)
# 
# POL_NUTS = readxl::read_xlsx('data/Poland/EU-28_LAU_2016.xlsx', sheet='PL') %>%
#   select(NUTS_3, LAU2_NAT_CODE, NAME_1)
# 
# POL_regions = left_join(POL_regions, POL_NUTS, by = c('LAU2.code' = 'LAU2_NAT_CODE')) %>%
#   select(GUS.code, NTS.4...Name, NUTS_3, NTS.3...Name) %>%
#   distinct(GUS.code, NTS.4...Name, .keep_all = TRUE)
# test = left_join(POL_regions, POL_GUS_NUTS3, by = c("GUS.code" = 'GUS'))


POL_cod_deaths = left_join(POL_cod_deaths, POL_GUS_NUTS3, by = c('RegionCode' = 'GUS')) %>%
  group_by(Country, NUTS3,  Year, AgeGroup, Sex, CauseCategory) %>%
  summarise(Deaths = sum(Deaths)) %>%
  rename(RegionCode = NUTS3)

POL_pop = read.csv(file = 'data/Poland/POL_Pop_NUTS3_2002-2023.csv', sep = ';') %>%
  filter(Year >= 2000, Age != 'TOT', Sex != 'b') %>%
  select(Country, RegionCode, Year, Sex, Age, Pop = Value) %>%
  mutate(CauseCategory = 'All',
         AgeGroup = as.integer(Age), 
         AgeGroup = case_when(AgeGroup %in% c(0, 5, 10, 15, 20) ~ 0,
                              AgeGroup %in% c(25, 30, 35, 40, 45) ~ 25,
                              AgeGroup %in% c(50, 55, 60, 65) ~ 50,
                              AgeGroup %in% c(70, 75, 80) ~ 70,
                              AgeGroup %in% c(85, 90, 95) ~ 85)) %>%
  group_by(Country, RegionCode, Year, Sex, AgeGroup) %>%
  summarise(Pop = sum(Pop))

# ESP <- read.csv(file = "data/european_standard_population.csv") %>%
#   mutate(AgeGroup = as.numeric(str_sub(AgeGroup, start = 1, end = 2)),
#          AgeGroup = case_when(AgeGroup %in% c(0, 1, 5, 10, 15, 20) ~ 0,
#                               AgeGroup %in% c(25, 30, 35, 40, 45) ~ 25,
#                               AgeGroup %in% c(50, 55, 60, 65) ~ 50,
#                               AgeGroup %in% c(70, 75, 80) ~ 70,
#                               AgeGroup %in% c(85, 90, 95) ~ 85)) %>%
#   group_by(AgeGroup) %>%
#   summarise(ESPw = sum(EuropeanStandardPopulation)/100000)

POL_cod_deaths = left_join(POL_cod_deaths, POL_pop)
  # left_join(., ESP) %>%
  # ungroup() %>%
  # group_by(Country, RegionCode, Year, Sex, CauseCategory) %>%
  # summarise(cdr = 1000*sum(Deaths)/sum(Pop),
  #           sdr = 1000*sum(ESPw*Deaths/Pop),
  #           Pop = sum(Pop))

POL_deaths = rbind(POL_cod_deaths, POL_all_deaths) %>%
  arrange(RegionCode, Year, Sex, CauseCategory, AgeGroup)

saveRDS(POL_deaths, 'temp/POL_for_smooth.rds')

# ggplot(POL_deaths, aes(x = Year, y = sdr, color = CauseCategory, group = interaction(RegionCode, CauseCategory, Sex)))+
#   geom_line()+
#   facet_grid(. ~ Sex)

