# Subnational mortality convergence in Europe by cause of death
# Rok Hrzic (r.hrzic (at) maastrichtuniversity.nl)
# Version September 2024

rm(list = ls())

require(tidyverse)

#all cause

EST_all_deaths = read.csv(file = 'data/Estonia/EST_Deaths_LAU1_1989-2023.csv', sep = ';') %>%
  filter(Year >= 2000, Age != 'TOT', Age != 'UNK', Sex != 'b', RegionCode != 'UNK') %>%
  select(Country, RegionCode, Year, Sex, Age, Deaths = Value) %>%
  mutate(CauseCategory = 'All',
         Age = as.integer(Age), 
         AgeGroup = DemoTools::calcAgeAbr(Age),
         AgeGroup = ifelse(AgeGroup >= 95, 95, AgeGroup)) %>%
  group_by(Country, RegionCode, Year, Sex, AgeGroup, CauseCategory) %>%
  summarise(Deaths = sum(Deaths))
  


EST_pop = read.csv(file = 'data/Estonia/EST_Mid-year_Pop_LAU1_2000-2023.csv', sep = ';') %>%
  filter(Year >= 2000, Age != 'TOT', Age != 'UNK', Sex != 'b', RegionCode != 'UNK') %>%
  select(Country, RegionCode, Year, Sex, Age, Pop = Value) %>%
  mutate(AgeGroup = as.integer(Age),
         Pop = str_replace(Pop, ',', '.'),
         Pop = as.numeric(Pop),
         AgeGroup = ifelse(AgeGroup >= 95, 95, AgeGroup)) %>%
  group_by(Country, RegionCode, Year, Sex, AgeGroup) %>%
  summarise(Pop = sum(Pop))


# ESP <- read.csv(file = "data/european_standard_population.csv") %>%
#   mutate(AgeGroup = as.numeric(str_sub(AgeGroup, start = 1, end = 2))) %>%
#   group_by(AgeGroup) %>%
#   summarise(ESPw = sum(EuropeanStandardPopulation)/100000)

EST_all_deaths = left_join(EST_all_deaths, EST_pop)
  # left_join(., ESP) %>%
  # ungroup() %>%
  # group_by(Country, RegionCode, Year, Sex, CauseCategory) %>%
  # mutate(Pop = ifelse(Pop == 0, 0.5, Pop)) %>%
  # summarise(cdr = 1000*sum(Deaths)/sum(Pop),
  #           sdr = 1000*sum(ESPw*(Deaths/Pop)),
  #           Pop = sum(Pop))


saveRDS(EST_all_deaths, 'temp/EST_LAU_all_for_smooth.rds')



# cod lau

EST_cod_rates <- read.csv(file = "data/Estonia/Estonia maakond SIM list CDRs and SDRs.csv", sep=",") %>%
  filter(CountyID != 999 & Year %in% 2000:2017) %>%
  select(Year, Sex, RegionCode = CountyID, CauseCode = CauseN, cdr = CDR, sdr = SDR2013) %>%
  mutate(Country = 'EST',
         Sex = ifelse(Sex == 1, 'm', 'f'),
         RegionCode = as.character(RegionCode),
         RegionCode = paste0('EE', RegionCode),
         CauseCategory = case_when(CauseCode %in% 7:23 ~ 'Cancer',
                                   CauseCode %in% 31:36 ~ 'Cardiovascular',
                                   CauseCode %in% 51:55 ~ 'External',
                                   TRUE ~ 'Other')) %>%
  group_by(Country, RegionCode, Year, Sex, CauseCategory) %>%
  summarise(cdr = sum(cdr),
            sdr = sum(sdr)/100)


## The number of causes reported varies over time
## No missing values reported, so likely the missing were simply dropped from the dataset
## Which are missing?

all.causes = expand_grid(Country = unique(EST_cod_rates$Country),
                         RegionCode = unique(EST_cod_rates$RegionCode), 
                         Year = unique(EST_cod_rates$Year),
                         Sex = unique(EST_cod_rates$Sex),
                         CauseCategory = unique(EST_cod_rates$CauseCategory))

EST_cod_rates = left_join(all.causes, EST_cod_rates) %>%
  arrange(RegionCode, Year, Sex, CauseCategory) %>%
  mutate(sdr = ifelse(is.na(sdr), 0, sdr),
         cdr = ifelse(is.na(cdr), 0, cdr))

rm(all.causes)

## LAU1

# EST_cod_LAU = rbind(EST_all_deaths, EST_cod_rates) %>%
  # arrange(RegionCode, Year, Sex, CauseCategory, AgeGroup) %>%
  # fill(Pop)

EST_LAU_pop = EST_pop %>%
  ungroup()%>%
  group_by(Country, RegionCode, Year, Sex) %>%
  summarise(Pop = sum(Pop))

EST_cod_rates = left_join(EST_cod_rates, EST_pop)

  
saveRDS(EST_cod_rates, 'temp/EST_LAU_cod_rates.rds')


##NUTS3

EST_regions = read.csv('data/Estonia/EE_LAU_2016.csv', sep=';') %>%
  select(NUTS_3, RegionCode = LAU1_NAT_CODE) %>%
  distinct(NUTS_3, RegionCode, .keep_all = TRUE) %>%
  mutate(RegionCode = as.character(RegionCode),
         RegionCode = paste0('EE', RegionCode)) %>%
  ungroup() %>%
  mutate(Name = case_when(RegionCode == 'EE37' ~ 'Harju county',
                          RegionCode == 'EE39' ~ 'Hiiu county',
                          RegionCode == 'EE44' ~ 'Ida-Viru county',
                          RegionCode == 'EE51' ~ 'Järva county',
                          RegionCode == 'EE49' ~ 'Jõgeva county',
                          RegionCode == 'EE57' ~ 'Lääne county',
                          RegionCode == 'EE59' ~ 'Lääne-Viru county',
                          RegionCode == 'EE67' ~ 'Pärnu county',
                          RegionCode == 'EE65' ~ 'Põlva county',
                          RegionCode == 'EE70' ~ 'Rapla county',
                          RegionCode == 'EE74' ~ 'Saare county',
                          RegionCode == 'EE78' ~ 'Tartu county',
                          RegionCode == 'EE82' ~ 'Valga county',
                          RegionCode == 'EE84' ~ 'Viljandi county',
                          RegionCode == 'EE86' ~ 'Võru county'))

EST_NUTS3_pop = EST_pop %>%
  left_join(., EST_regions) %>%
  group_by(Country, RegionCode, NUTS_3, Year, Sex, AgeGroup) %>%
  summarise(Pop = sum(Pop)) 

EST_NUTS3_all_deaths = EST_all_deaths %>%
  left_join(., EST_regions) %>%
  left_join(., EST_NUTS3_pop) %>%
  group_by(Country, NUTS_3, Year, Sex, CauseCategory, AgeGroup) %>%
  summarise(Deaths = sum(Deaths),
            Pop = sum(Pop)) %>%
  rename(RegionCode=NUTS_3) %>%
  arrange(RegionCode, Year, Sex, CauseCategory, AgeGroup)

saveRDS(EST_NUTS3_all_deaths, 'temp/EST_NUTS3_all_for_smooth.rds')

EST_NUTS3_pop = EST_pop %>%
  ungroup()%>%
  left_join(., EST_regions) %>%
  group_by(Country, NUTS_3, RegionCode, Year, Sex) %>%
  summarise(Pop = sum(Pop))

EST_NUTS3_cod_rates = EST_cod_rates %>%
  left_join(., EST_regions) %>%
  left_join(., EST_NUTS3_pop) %>%
  group_by(Country, NUTS_3, Year, Sex, CauseCategory) %>%
  summarise(cdr = weighted.mean(cdr, Pop),
            sdr = weighted.mean(sdr, Pop),
            Pop = sum(Pop)) %>%
  rename(RegionCode=NUTS_3) %>%
  arrange(RegionCode, Year, Sex, CauseCategory)

saveRDS(EST_NUTS3_cod_rates, 'temp/EST_NUTS3_cod_rates.rds')
