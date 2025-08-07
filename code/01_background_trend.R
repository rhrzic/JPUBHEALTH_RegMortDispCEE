# Subnational mortality convergence in Europe by cause of death
# Rok Hrzic (r.hrzic (at) maastrichtuniversity.nl)
# August 2025

rm(list = ls())

require(tidyverse)
require(cowplot)

df = readxl::read_xlsx('data/Fig_1_LE_background.xlsx') %>%
  mutate(Sex = ifelse(Sex == 'Females', 'Women', 'Men'))

p1 = ggplot(df)+
  geom_line(aes(x = Year, y = LE, colour = Country))+
  annotate("rect", xmin = 2020, xmax = Inf, ymin = -Inf,   ymax = Inf,fill = "lightgrey", alpha = 0.5)+  
  annotate("text", x = 2020.2, y = 67, label = "COVID-19\npandemic", hjust = 0, vjust = -0.25, size = 3.5)+
  facet_grid(.~Sex)+
  ylab('Life expectancy at birth (years)')+
  xlab('Year')+
  theme_bw()+
  theme(legend.position = 'top',
        legend.title = element_blank(),
        text=element_text(size=14))+
  guides(colour = guide_legend(nrow = 1))

ggsave('plots/fig1.pdf', p1, height = 210, width = 297, units = 'mm')
ggsave('plots/fig1.png', p1, height = 210, width = 297, units = 'mm')

