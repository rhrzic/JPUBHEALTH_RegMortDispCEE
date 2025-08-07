# Subnational mortality convergence in Europe by cause of death
# Rok Hrzic (r.hrzic (at) maastrichtuniversity.nl)
# April 2025

rm(list = ls())

require(tidyverse)

# all causes

# OAG is 85

ROU_all_deaths = read.csv(file = 'data/Romania/ROU_Deaths_NUTS3_1990-2022.csv', sep = ";") %>%
  filter(Year >= 2000 & Age != 'TOT' & Sex != 'b') %>%
  mutate(CauseCategory = 'All',
         AgeGroup = as.integer(Age)) %>%
  select(Country, RegionCode, Year, Sex, AgeGroup, CauseCategory, Deaths = Value)

ROU_pop = read.csv(file = 'data/Romania/ROU_Mid-year_Pop_NUTS3_1991-2022.csv', sep = ";") %>%
  filter(Year >= 2000 & Age != 'TOT' & Sex != 'b') %>%
  mutate(Age = as.integer(Age), 
         AgeGroup = DemoTools::calcAgeAbr(Age),
         Value = str_replace(Value, ',', '.'),
         Value = as.numeric(Value)) %>%
  group_by(Country, RegionCode, Year, Sex, AgeGroup) %>%
  summarise(Value = sum(Value)) %>%
  select(Country, RegionCode, Year, Sex, AgeGroup, Pop = Value)


ESP <- read.csv(file = "data/european_standard_population.csv") %>%
  mutate(AgeGroup = as.numeric(str_sub(AgeGroup, start = 1, end = 2)),
         AgeGroup = ifelse(AgeGroup >= 85, 85, AgeGroup)) %>%
  group_by(AgeGroup) %>%
  summarise(ESPw = sum(EuropeanStandardPopulation)/100000)
  
ROU_all_deaths = left_join(ROU_all_deaths, ROU_pop)
  # left_join(., ESP) %>%
  # ungroup() %>%
  # group_by(Country, RegionCode, Year, Sex, CauseCategory) %>%
  # summarise(cdr = 1000*sum(Deaths)/sum(Pop),
  #           sdr = 1000*sum(ESPw*Deaths/Pop),
  #           Pop = sum(Pop))

rm(ESP)

# cod

ROU_cod_deaths <- read.csv(file = "data/Romania/Romania_CoD_Rok.csv", sep=";") %>%
  filter(Year %in% 2000:2019 & RegionCode != 'RO' & Sex != 'b' & CauseCode != 'all causes') %>%
  mutate(Value = gsub(',', '.', Value, fixed = T),
         Value = as.numeric(Value),
         CauseCode = as.numeric(CauseCode),
         RegionCode = ifelse(RegionCode %in% c('RO321', 'RO322'), 'RO323', RegionCode)) %>%
  filter(CauseCode <= 14) %>%
  group_by(Country, RegionCode, Year, Sex, Age, CauseCode) %>%
  summarise(Value = sum(Value)) %>%
  select(Country, RegionCode, Year, Sex, AgeGroup = Age, CauseCode, Deaths = Value)

#There seems to be an issue with external mortality causes in RO321 for the oldest men. I'll use the female numbers for now.

# m_RO321 = ROU_deaths %>% filter(RegionCode == 'RO321', CauseCode == 13, Sex == 'f', Age == '85') %>%
#   mutate(Sex = 'm')
# 
# ROU_deaths = ROU_deaths %>% filter(!(RegionCode == 'RO321' & CauseCode == 13 & Sex == 'm' & Age == '85'))
# 
# ROU_deaths = rbind(ROU_deaths, m_RO321)


##There is some sort of mixing up of the TOT age category and the 85 age category for the causes 13 and 14. 
#The 85 category starts picking up the all.cause for all ages value in 2000
#Untrustworthy categories then include: 85 for cause 13 and 14 in men; TOT age category for causes 13 and 14 in men. 

# ROU_deaths = ROU_deaths %>% 
#   mutate(Value = case_when(CauseCode %in% c('13', '14') & !is.na(Value) & Value > 1000 & Age == '85' & lag(Age) == '85' & lag(Age, 2) == '85' & RegionCode == 'RO321' ~ NA, 
#                          TRUE ~ Value),
#          Age = case_when(CauseCode %in% c('13', '14') & !is.na(Value) & Value > 20 & Age == '85' & lag(Age) == '85' & lag(Age, 2) == '85' ~ 'TOT', 
#                          TRUE ~ Age),
#          RegionCode = case_when(CauseCode %in% c('13', '14') & !is.na(Value) & Value > 20 & Age == 'TOT' & lead(Age) == 'TOT' ~ lead(RegionCode), 
#                                 TRUE ~ RegionCode))
# 

# ROU_deaths %>%
#   mutate(Age = factor(Age, levels=c('0', '5', '25', '50', '70', '85', 'TOT')),
#          CauseCode = factor(CauseCode, levels = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '101', '301', '601', '602', 'all.causes'))) %>%
#   group_by(Year, Sex, RegionCode, Age, CauseCode) %>%
#   summarise(missings = is.na(Value)) %>%
#   ggplot(aes(x = Year, y = RegionCode, fill = missings))+
#   geom_tile()+
#   facet_grid(CauseCode ~ Sex+Age)

## Missing values seems MCAR

# ROU_deaths %>%
#   filter(!is.na(Value) & Value <10) %>%
#   ggplot(aes(x = Value)) +
#   geom_histogram(binwidth = 1)
# 
# table(ROU_deaths$Value)

# Values of 1 and 2 are underrepresented, likely designated as NA for good order

# ROU_deaths %>%
#   mutate(Age = factor(Age, levels=c('0', '5', '25', '50', '70', '85', 'TOT')),
#          CauseCode = factor(CauseCode, levels = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '101', '301', '601', '602', 'all.causes'))) %>%
#   group_by(Year, Sex, RegionCode, Age, CauseCode) %>%
#   summarise(datapoints = n()) %>%
#   ggplot(aes(x = Year, y = RegionCode, fill = datapoints))+
#   geom_tile()+
#   facet_grid(CauseCode ~ Sex+Age)

# ROU_deaths %>%
#   filter(CauseCode == 'all.causes' & Age %in% c('0', '5', 'TOT')) %>%
#   mutate(Age = factor(Age, levels=c('0', '5', '25', '50', '70', '85', 'TOT')),
#          CauseCode = factor(CauseCode, levels = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '101', '301', '601', '602', 'all.causes'))) %>%
#   group_by(Year, Sex, RegionCode, Age, CauseCode) %>%
#   summarise(datapoints = n()) %>%
#   ggplot(aes(x = Year, y = RegionCode, fill = datapoints))+
#   geom_tile()+
#   geom_text(aes(label = datapoints))+
#   facet_grid(CauseCode ~ Sex+Age)

## There are data errors in 2008 for women at ages 0, 5, and total. 
## Missing data for RO321 ages 0, 5, tot and  RO422 and RO124 at age 5
## Duplicates for RO422 and RO124 at ages 0 and tot
## Addition: The errors seem to be more widespread throughout the different regions, but always focused on 2008 and women


# ROU_deaths = ROU_deaths %>%
#   filter(!(Age %in% c('0', '5', 'TOT') & Year == 2008 & Sex == 'f'))


# ROU_deaths = ROU_deaths %>%
#   mutate(Value = ifelse(Age == '85' & Sex == 'm' & CauseCode %in% c(13, 14) & Value > 25, NULL, Value))

# ROU_deaths %>%
#   mutate(Age = factor(Age, levels=c('0', '5', '25', '50', '70', '85', 'TOT')),
#          CauseCode = factor(CauseCode, levels = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '101', '301', '601', '602', 'all.causes'))) %>%
#   group_by(Year, Sex, RegionCode, Age, CauseCode) %>%
#   summarise(datapoints = n()) %>%
#   ggplot(aes(x = Year, y = RegionCode, fill = datapoints))+
#   geom_tile()+
#   geom_text(aes(label = datapoints))+
#   facet_grid(CauseCode ~ Sex+Age)

# ROU_deaths = ROU_deaths %>%
#   mutate(Value = ifelse(is.na(Value), 1.5, Value))

# ROU_deaths %>%
#   mutate(Age = factor(Age, levels=c('0', '5', '25', '50', '70', '85', 'TOT')),
#          CauseCode = factor(CauseCode, levels = c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '101', '301', '601', '602', 'all.causes'))) %>%
#   group_by(Year, Sex, RegionCode, Age, CauseCode) %>%
#   summarise(missings = is.na(Value)) %>%
#   ggplot(aes(x = Year, y = RegionCode, fill = missings))+
#   geom_tile()+
#   facet_grid(CauseCode ~ Sex+Age)


# ROU_deaths = ROU_deaths %>%
#   group_by(Country, RegionCode, Year, Sex, Age, CauseCode) %>%
#   summarise(Deaths = sum(Value, na.rm = T)) %>%
#   ungroup()
# 
# fix = expand_grid(Country = unique(ROU_deaths$Country), 
#                   RegionCode = unique(ROU_deaths$RegionCode), 
#                   Year = 2008, Sex = 'f',
#                   Age = c('0', '5', 'TOT'),
#                   CauseCode = unique(ROU_deaths$CauseCode), 
#                   Deaths = NA)
# 
# ROU_deaths = rbind(ROU_deaths, fix)
# 
# ROU_deaths = ROU_deaths %>%
#   group_by(Country, RegionCode, Sex, Age, CauseCode) %>%
#   arrange(Country, RegionCode, Sex, Age, CauseCode, Year) %>%
#   mutate(Deaths = na.approx(Deaths)) %>%
#   ungroup()

## Testing the totals

# ROU_deaths %>%
#   pivot_wider(names_from = Age, values_from = Deaths) %>%
#   mutate(diff=TOT - (`0` + `25` + `5` + `50` + `70` + `85`)) %>%
#   group_by(Year, Sex, RegionCode, CauseCode) %>%
#   summarise(diff = sum(diff)) %>%
#   ggplot(aes(x = Year, y = RegionCode, fill = diff))+
#   geom_tile()+
#   facet_grid(CauseCode ~ Sex)

# ROU_deaths %>%
#   filter(CauseCode == '13' & Sex == 'm') %>%
#   pivot_wider(names_from = Age, values_from = Deaths) %>%
#   mutate(diff=TOT - (`0` + `25` + `5` + `50` + `70` + `85`)) %>%
#   group_by(Year, Sex, RegionCode, CauseCode) %>%
#   summarise(diff = sum(diff)) %>%
#   ggplot(aes(x = Year, y = RegionCode, fill = diff))+
#   geom_tile()+
#   facet_grid(CauseCode ~ Sex)

# 
# ROU_deaths %>%
#   pivot_wider(names_from = CauseCode, values_from = Deaths) %>%
#   mutate(diff= all.causes - (`1`+ `2` + `3` + `4` + `5` + `6` + `7`+ `8` + `9` + `10` + `11` + `12` + `13` + `14`)) %>%
#   group_by(Year, Sex, RegionCode, Age) %>%
#   summarise(diff = sum(diff)/all.causes) %>%
#   ggplot(aes(x = Year, y = RegionCode, fill = diff))+
#   geom_tile()+
#   facet_grid(Age ~ Sex)
# 
# ROU_deaths %>%
#   filter(!CauseCode %in% c('101', '301', '601', '602')) %>%
#   pivot_wider(names_from = c(Age, CauseCode), values_from = Deaths) %>%
#   ggplot(aes(x = Year, y = RegionCode, fill = `70_13`))+
#   geom_tile() +
#   facet_grid(.~Sex)
# 

## No need for this code now.

# fix_coef = ROU_deaths %>%
#   filter(!CauseCode %in% c('101', '301', '601', '602') & Year < 2000 & Sex == 'm') %>%
#   pivot_wider(names_from = c(Age, CauseCode), values_from = Deaths) %>%
#   group_by(RegionCode) %>%
#   summarise(coef_TOT_13 = mean(`TOT_13`/`TOT_all.causes`, na.rm=T),
#          coef_TOT_14 = mean(`TOT_14`/`TOT_all.causes`, na.rm=T),
#          coef_TOT = 
#          coef_85_13 = mean(`85_13`/`TOT_13`, na.rm=T),
#          coef_85_14 = mean(`85_14`/`TOT_14`, na.rm=T)) %>%
#   select(RegionCode, coef_TOT_13, coef_TOT_14, coef_85_13, coef_85_14)
# 
# 
# ROU_deaths_new = ROU_deaths %>%
#   filter(!CauseCode %in% c('101', '301', '601', '602')) %>%
#   left_join(., fix_coef) %>%
#   pivot_wider(names_from = c(Age, CauseCode), values_from = Deaths) %>%
#   mutate(`TOT_13` = ifelse(Year > 1999 & Sex == 'm', `TOT_all.causes`*coef_TOT_13, `TOT_13`),
#          `TOT_14` = ifelse(Year > 1999 & Sex == 'm', `TOT_all.causes`*coef_TOT_14, `TOT_14`),
#          `85_13` = ifelse(Year > 1999 & Sex == 'm', `TOT_13`*coef_85_13, `85_13`),
#          `85_14` = ifelse(Year > 1999 & Sex == 'm', `TOT_14`*coef_85_14, `85_14`))


# ROU_deaths_new %>%
#   ggplot(aes(x = Year, y = RegionCode, fill = `TOT_13`))+
#   geom_tile() +
#   facet_grid(.~Sex)

# ROU_deaths_new = ROU_deaths_new %>%
#   select(-c(coef_TOT_13, coef_TOT_14, coef_85_13, coef_85_14)) %>%
#   pivot_longer(cols = `0_1`:`TOT_all.causes`, names_to = c('Age', 'CauseCode'), values_to = 'Deaths', names_sep = "_")
#   
# 
# rm(ROU_deaths_new, ROU_deaths_fine, ROU_deaths_fix)


## Categorizion of cod
## The provided data has its own specific categorisation of CoDs

ROU_SIM <- read.csv(file = "data/Romania/Romania_SIM_Rok.csv", sep=";")

ROU_cod_deaths = ROU_cod_deaths %>%
  left_join(., ROU_SIM)


rm(ROU_SIM)

## I summarise by cause category (not interested in the finer distinctions)

ROU_cod_deaths = ROU_cod_deaths %>%
  group_by(Country, RegionCode, Year, Sex, AgeGroup, CauseCategory) %>%
  summarise(Deaths = sum(Deaths))

all.causes = expand_grid(Country = 'ROU',
                         RegionCode = unique(ROU_cod_deaths$RegionCode), 
                         Year = unique(ROU_cod_deaths$Year),
                         Sex = unique(ROU_cod_deaths$Sex),
                         AgeGroup = unique(ROU_cod_deaths$AgeGroup),
                         CauseCategory = unique(ROU_cod_deaths$CauseCategory))


ROU_cod_deaths = ROU_cod_deaths %>%
  right_join(all.causes) %>% 
  mutate(Deaths = ifelse(is.na(Deaths), 0, Deaths))


ROU_pop = ROU_pop %>%
  mutate(AgeGroup = as.numeric(str_sub(AgeGroup, start = 1, end = 2)),
         AgeGroup = case_when(AgeGroup %in% c(0,1) ~ 0,
                              AgeGroup %in% c(5, 10, 15, 20) ~ 5,
                              AgeGroup %in% c(25, 30, 35, 40, 45) ~ 25,
                              AgeGroup %in% c(50, 55, 60, 65) ~ 50,
                              AgeGroup %in% c(70, 75, 80) ~ 70,
                              AgeGroup %in% c(85, 90, 95) ~ 85)) %>%
  group_by(Country, RegionCode, Year, Sex, AgeGroup) %>%
  summarise(Pop = sum(Pop))

# ESP <- read.csv(file = "data/european_standard_population.csv") %>%
#   mutate(AgeGroup = as.numeric(str_sub(AgeGroup, start = 1, end = 2)),
#          AgeGroup = case_when(AgeGroup %in% c(0,1) ~ 0,
#                               AgeGroup %in% c(5, 10, 15, 20) ~ 5,
#                               AgeGroup %in% c(25, 30, 35, 40, 45) ~ 25,
#                               AgeGroup %in% c(50, 55, 60, 65) ~ 50,
#                               AgeGroup %in% c(70, 75, 80) ~ 70,
#                               AgeGroup %in% c(85, 90, 95) ~ 85)) %>%
#   group_by(AgeGroup) %>%
#   summarise(ESPw = sum(EuropeanStandardPopulation)/100000)

ROU_cod_deaths = left_join(ROU_cod_deaths, ROU_pop)
  # left_join(., ESP) %>%
  # ungroup() %>%
  # group_by(Country, RegionCode, Year, Sex, CauseCategory) %>%
  # summarise(cdr = 1000*sum(Deaths)/sum(Pop),
  #           sdr = 1000*sum(ESPw*Deaths/Pop),
  #           Pop = sum(Pop))

ROU_deaths = rbind(ROU_cod_deaths, ROU_all_deaths) %>%
  arrange(RegionCode, Year, Sex, CauseCategory, AgeGroup)

saveRDS(ROU_deaths, 'temp/ROU_for_smooth.rds')

# ggplot(ROU_deaths, aes(x = Year, y = sdr, color = CauseCategory, group = interaction(RegionCode, CauseCategory, Sex)))+
#   geom_line()+
#   facet_grid(. ~ Sex)
