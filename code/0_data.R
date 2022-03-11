# Subnational mortality convergence in Europe by cause of death
# Rok Hrzic (r.hrzic (at) maastrichtuniversity.nl)
# Version March 2022

require(dplyr)

SWE_deaths <- read.csv(file = "data/Sweden/Sweden_CoD_harmonized_Rok.csv", sep=";")

SWE_deaths %>%
  group_by(RegionCode, Year, Sex, Age) %>%
  summarise(Value = sum(Value), CauseCode = 'Total') %>%
  mutate(CauseCode = 'Total') %>%
  rbind(., SWE_deaths)



SWE_pop <- read.csv(file = "data/Sweden/Sweden_CoD_harmonized_Rok.csv", sep=";")
