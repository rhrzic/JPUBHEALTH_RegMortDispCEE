# Subnational mortality convergence in Europe by cause of death
# Rok Hrzic (r.hrzic (at) maastrichtuniversity.nl)
# Version March 2022

require(dplyr)
require(ggplot2)
require(grid)
require(gridExtra)
require(cowplot)
require(directlabels)

national_data_mean <- national_data %>%
  group_by(sex, year, NMS) %>%
  summarise(mean_ple = weighted.mean(ple, pop))

p1 <- ggplot() +
  geom_rect(aes(xmin = 2004, xmax = 2007, ymin = -Inf, ymax = Inf), color = "gray95", fill = "gray95") +
  geom_line(data = filter(national_data, sex == "Men"), aes(x = year, y = ple, group = country, colour = NMS), size = 0.5) +
  stat_summary_bin(data = filter(national_data, sex == "Men"), aes(x = year, y = ple, group = NMS, colour = NMS, label=round(..y..,1)),  
                   fun.min = min, 
                   fun.max = max,
                   geom = "linerange",
                   breaks = c(1990, 1990.5),
                   position = position_nudge(x=-1.5))+
  stat_summary_bin(data = filter(national_data, sex == "Men"), aes(x = year, y = ple, group = NMS, colour = NMS, label=round(..y..,1)),  
                   fun = min, 
                   geom = "text",
                   breaks = c(1990, 1990.5),
                   position = position_nudge(x=-2.5),
                   size = 3)+
  stat_summary_bin(data = filter(national_data, sex == "Men"), aes(x = year, y = ple, group = NMS, colour = NMS, label=round(..y..,1)),  
                   fun = max, 
                   geom = "text",
                   breaks = c(1990, 1990.5),
                   position = position_nudge(x=-2.5),
                   size = 3) +
  stat_summary_bin(data = filter(national_data, sex == "Men"), aes(x = year, y = ple, group = NMS, colour = NMS, label=round(..y..,1)),  
                   fun.min = min, 
                   fun.max = max,
                   geom = "linerange",
                   breaks = c(2004, 2004.1)) +
  stat_summary_bin(data = filter(national_data, sex == "Men"), aes(x = year, y = ple, group = NMS, colour = NMS, label=round(..y..,1)),  
                   fun = min, 
                   geom = "text",
                   breaks = c(2004, 2004.1),
                   position = position_nudge(y=-0.3),
                   size = 3)+
  stat_summary_bin(data = filter(national_data, sex == "Men"), aes(x = year, y = ple, group = NMS, colour = NMS, label=round(..y..,1)),  
                   fun = max, 
                   geom = "text",
                   breaks = c(2004, 2004.1),
                   position = position_nudge(y=0.3),
                   size = 3)+
  stat_summary_bin(data = filter(national_data, sex == "Men"), aes(x = year, y = ple, group = NMS, colour = NMS, label=round(..y..,1)),  
                   fun.min = min, 
                   fun.max = max,
                   geom = "linerange",
                   breaks = c(2007, 2007.1)) +
  stat_summary_bin(data = filter(national_data, sex == "Men"), aes(x = year, y = ple, group = NMS, colour = NMS, label=round(..y..,1)),  
                   fun = min, 
                   geom = "text",
                   breaks = c(2007, 2007.1),
                   position = position_nudge(y=-0.3),
                   size = 3)+
  stat_summary_bin(data = filter(national_data, sex == "Men"), aes(x = year, y = ple, group = NMS, colour = NMS, label=round(..y..,1)),  
                   fun = max, 
                   geom = "text",
                   breaks = c(2007, 2007.1),
                   position = position_nudge(y=0.3),
                   size = 3)+
  stat_summary_bin(data = filter(national_data, sex == "Men"), aes(x = year, y = ple, group = NMS, colour = NMS, label=round(..y..,1)), 
                   geom = "linerange", 
                   fun.min = min, 
                   fun.max = max,
                   breaks = c(2016.5,2017),
                   position = position_nudge(x=1.5)) +
  stat_summary_bin(data = filter(national_data, sex == "Men"), aes(x = year, y = ple, group = NMS, colour = NMS, label=round(..y..,1)),  
                   fun = min, 
                   geom = "text",
                   breaks = c(2016.5,2017),
                   position = position_nudge(x=2.5, y = 0.1),
                   size = 3)+
  stat_summary_bin(data = filter(national_data, sex == "Men"), aes(x = year, y = ple, group = NMS, colour = NMS, label=round(..y..,1)),  
                   fun = max, 
                   geom = "text",
                   breaks = c(2016.5,2017),
                   position = position_nudge(x=2.5, y = -0.1),
                   size = 3)+
  geom_dl(data = filter(national_data, sex == "Men" & country %in% c("SE", "FR", "ES", "DE", "NL", "AT", "PT", "SI", "CZ", "SK", "PL", "HU", "EE", "LT", "LV")), aes(x = year, y = ple, label = country, colour = NMS), method=list("first.points", cex=0.7))+
  geom_dl(data = filter(national_data, sex == "Men" & country %in% c("SE", "FR", "ES", "BE", "NL", "DE", "PT", "SI", "CZ", "SK", "PL", "HU", "EE", "LT", "LV")), aes(x = year, y = ple, label = country, colour = NMS), method=list("last.points", cex=0.7))+
  geom_line(data = filter(national_data_mean, NMS == "New Member States", sex == "Men"), aes(x = year, y = mean_ple, colour = NMS), size = 2)+
  geom_line(data = filter(national_data_mean, NMS == "Old Member States", sex == "Men"), aes(x = year, y = mean_ple, colour = NMS), size = 2)+
  annotate("text", x = 2004.1, y = 62, label = "Short-term\naccession\neffects", hjust = "left", size = 3) +
  facet_wrap(. ~ sex) +
  xlab("Year") +
  ylab("Life expectancy at birth") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6))+
  theme_bw()+
  scale_colour_grey()+
  theme(legend.position = "none")

p2 <- ggplot() +
  geom_rect(aes(xmin = 2004, xmax = 2007, ymin = -Inf, ymax = Inf), color = "gray95", fill = "gray95") +
  geom_line(data = filter(national_data, sex == "Women"), aes(x = year, y = ple, group = country, colour = NMS), size = 0.5) +
  stat_summary_bin(data = filter(national_data, sex == "Women"), aes(x = year, y = ple, group = NMS, colour = NMS, label=round(..y..,1)),  
                   fun.min = min, 
                   fun.max = max,
                   geom = "linerange",
                   breaks = c(1990, 1990.5),
                   position = position_nudge(x=-1.5))+
  stat_summary_bin(data = filter(national_data, sex == "Women"), aes(x = year, y = ple, group = NMS, colour = NMS, label=round(..y..,1)),  
                   fun = min, 
                   geom = "text",
                   breaks = c(1990, 1990.5),
                   position = position_nudge(x=-2.5),
                   size = 3)+
  stat_summary_bin(data = filter(national_data, sex == "Women"), aes(x = year, y = ple, group = NMS, colour = NMS, label=round(..y..,1)),  
                   fun = max, 
                   geom = "text",
                   breaks = c(1990, 1990.5),
                   position = position_nudge(x=-2.5),
                   size = 3)+
  stat_summary_bin(data = filter(national_data, sex == "Women"), aes(x = year, y = ple, group = NMS, colour = NMS, label=round(..y..,1)),  
                   fun.min = min, 
                   fun.max = max,
                   geom = "linerange",
                   breaks = c(2004, 2004.1)) +
  stat_summary_bin(data = filter(national_data, sex == "Women"), aes(x = year, y = ple, group = NMS, colour = NMS, label=round(..y..,1)),  
                   fun = min, 
                   geom = "text",
                   breaks = c(2004, 2004.1),
                   position = position_nudge(y=-0.3),
                   size = 3)+
  stat_summary_bin(data = filter(national_data, sex == "Women"), aes(x = year, y = ple, group = NMS, colour = NMS, label=round(..y..,1)),  
                   fun = max, 
                   geom = "text",
                   breaks = c(2004, 2004.1),
                   position = position_nudge(y=0.3),
                   size = 3)+
  stat_summary_bin(data = filter(national_data, sex == "Women"), aes(x = year, y = ple, group = NMS, colour = NMS, label=round(..y..,1)),  
                   fun.min = min, 
                   fun.max = max,
                   geom = "linerange",
                   breaks = c(2007, 2007.1)) +
  stat_summary_bin(data = filter(national_data, sex == "Women"), aes(x = year, y = ple, group = NMS, colour = NMS, label=round(..y..,1)),  
                   fun = min, 
                   geom = "text",
                   breaks = c(2007, 2007.1),
                   position = position_nudge(y=-0.3),
                   size = 3)+
  stat_summary_bin(data = filter(national_data, sex == "Women"), aes(x = year, y = ple, group = NMS, colour = NMS, label=round(..y..,1)),  
                   fun = max, 
                   geom = "text",
                   breaks = c(2007, 2007.1),
                   position = position_nudge(y=0.3),
                   size = 3)+
  stat_summary_bin(data = filter(national_data, sex == "Women"), aes(x = year, y = ple, group = NMS, colour = NMS, label=round(..y..,1)), 
                   geom = "linerange", 
                   fun.min = min, 
                   fun.max = max,
                   breaks = c(2016.5,2017),
                   position = position_nudge(x=1.5)) +
  stat_summary_bin(data = filter(national_data, sex == "Women"), aes(x = year, y = ple, group = NMS, colour = NMS, label=round(..y..,1)),  
                   fun = min, 
                   geom = "text",
                   breaks = c(2016.5,2017),
                   position = position_nudge(x=2.5),
                   size = 3)+
  stat_summary_bin(data = filter(national_data, sex == "Women"), aes(x = year, y = ple, group = NMS, colour = NMS, label=round(..y..,1)),  
                   fun = max, 
                   geom = "text",
                   breaks = c(2016.5,2017),
                   position = position_nudge(x=2.5),
                   size = 3)+
  geom_dl(data = filter(national_data, sex == "Women" & country %in% c("SE", "FR", "NL", "EL", "BE", "AT", "DE", "DK", "SI", "CZ", "SK", "PL", "HU", "EE", "LT", "LV")), aes(x = year, y = ple, label = country, colour = NMS), method=list("first.points", cex=0.7))+
  geom_dl(data = filter(national_data, sex == "Women" & country %in% c("SE", "ES", "FR", "IT", "PT", "DE", "DK", "SI", "CZ", "SK", "PL", "HU", "EE", "LT", "LV")), aes(x = year, y = ple, label = country, colour = NMS), method=list("last.points", cex=0.7))+
  geom_line(data = filter(national_data_mean, NMS == "New Member States", sex == "Women"), aes(x = year, y = mean_ple, colour = NMS), size = 2)+
  geom_line(data = filter(national_data_mean, NMS == "Old Member States", sex == "Women"), aes(x = year, y = mean_ple, colour = NMS), size = 2)+
  facet_wrap(. ~ sex) +
  xlab("Year") +
  ylab("Life expectancy at birth") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6))+
  theme_bw()+
  scale_colour_grey()+
  theme(legend.position = "none")

DK_women <- national_data %>%
  filter(sex == "Women", country == "DK")

PT_men <- national_data %>%
  filter(sex == "Men", country == "PT")


p3 <- national_data %>%
  filter(sex == "Men" & (NMS == "New Member States" | country == "PT")) %>%
  mutate(ple = ple - PT_men$ple) %>%
  ggplot(aes(x = year, y = ple, group = country), size = 0.8) +
  geom_rect(aes(xmin = 2004, xmax = 2007, ymin = -Inf, ymax = Inf), color = "gray95", fill = "gray95") +
  geom_line() +
  geom_dl(aes(label = country), method=list("last.points", cex=0.7))+
  geom_dl(aes(label = country), method=list("first.points", cex=0.7))+
  facet_wrap(. ~ sex) +
  xlab("Year") +
  ylab("Difference with Portuguese male life expectancy") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6))+
  theme_bw()+
  scale_colour_grey() +  
  labs(colour = "Country group") +
  theme(legend.position = "top", legend.direction = "horizontal")

p4 <- national_data %>%
  filter(sex == "Women" & (NMS == "New Member States" | country == "DK")) %>%
  mutate(ple = ple - DK_women$ple) %>%
  ggplot(aes(x = year, y = ple, group = country), size = 0.8) +
  geom_rect(aes(xmin = 2004, xmax = 2007, ymin = -Inf, ymax = Inf), color = "gray95", fill = "gray95") +
  geom_line() +
  geom_dl(aes(label = country), method=list("last.points", cex=0.7))+
  geom_dl(aes(label = country), method=list("first.points", cex=0.7))+
  facet_wrap(. ~ sex) +
  xlab("Year") +
  ylab("Difference with Danish female life expectancy") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6))+
  theme_bw()+
  scale_colour_grey() +  
  labs(colour = "Country group") +
  theme(legend.position = "top", legend.direction = "horizontal")

fig1 <- plot_grid(p1, p2, p3, p4, nrow = 2, ncol = 2, align = "hv", axis = "tblr", labels = NA)
ggsave('figures/fig1.eps', fig1, scale = 1, width = 12, height = 12, units = "in", device='eps', dpi=700)
ggsave('figures/fig1.png', fig1, scale = 1, width = 12, height = 12, units = "in", device='png', dpi=200)
