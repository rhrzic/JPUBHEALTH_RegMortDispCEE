# Subnational mortality convergence in Europe by cause of death
# Rok Hrzic (r.hrzic (at) maastrichtuniversity.nl)
# Version May 2025

rm(list = ls())

require(tidyverse)
require(cowplot)
require(flextable)
require(officer)

theil_T <- function(x, na.rm = TRUE) {
  x <- if (na.rm) x[!is.na(x)] else x
  mu <- mean(x)
  p  <- x / mu
  mean(p * log(p))
}


### Mortality trend

# df = readRDS('temp/sdr1.RDS') %>%
#   filter(Year %in% 2015:2022 & Country %in% c("CZE", "EST", "POL", "ROU") & CauseCategory == 'All') %>%
#   mutate(Sex = ifelse(Sex == 'm', 'Men', 'Women')) %>%
#   group_by(Country, RegionCode, Year, Sex) %>%
#   summarise(
#     # Point estimate: median of the bootstrap replicates
#     sdr_est = median(smoothed_sdr, na.rm = TRUE),
#     # 95% percentile-based CI
#     sdr_low = quantile(smoothed_sdr, probs = 0.025, na.rm = TRUE),
#     sdr_up = quantile(smoothed_sdr, probs = 0.975, na.rm = TRUE),
#     raw_sdr = first(raw_sdr)) %>%
#   ungroup()
# 
# 
# ggplot(df, aes(x = Year, color = Country, linetype = Sex, shape = Sex, group = interaction(RegionCode, Sex))) +
#   geom_line(aes(y = sdr_est), alpha = 0.5, linewidth = 0.75) +
#   geom_line(aes(y = sdr_low), alpha = 0.3, linewidth = 0.4) +
#   geom_line(aes(y = sdr_up), alpha = 0.3, linewidth = 0.4) +
#   geom_point(aes(y = raw_sdr)) +
#   facet_grid(. ~ Country, scales="free_y") +
#   ylab('Age-standardised death rate (per 1000)')+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 90),
#         legend.position = 'top',
#         text=element_text(size=14))
# 
# 
# df_diff = df %>%
#   group_by(Country, RegionCode, Year, Sex) %>%
#   summarise(sdr_diff = raw_sdr - sdr_est,
#             sdr_up = sdr_up - sdr_est,
#             sdr_low = sdr_low - sdr_est) %>%
#   mutate(annotation = ifelse(Year == 2015 & RegionCode == 'CZ010' & Sex == 'Men', "COVID-19\npandemic", ''))
# 
# ggplot(df_diff, aes(x = Year, group = interaction(RegionCode, Sex))) +
#   geom_rect(aes(xmin = 2019.5, xmax = Inf, ymin = -Inf, ymax =  Inf), fill  = "lightgrey", alpha = 0.1) +
#   geom_line(aes(y = sdr_low, color = Country, linetype = Sex), alpha = 0.3, linewidth = 0.4) +
#   geom_line(aes(y = sdr_up, color = Country, linetype = Sex), alpha = 0.3, linewidth = 0.4) +
#   geom_point(aes(y = sdr_diff, color = Country, shape = Sex)) +
#   geom_text(aes(x = 2020.2, y = Inf, label = annotation, hjust = 0, vjust = 1.5))+
#   facet_grid(Sex ~ Country) +
#   ylab('Difference between observed and forecast ASDR')+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 90),
#         legend.position = 'top',
#         text=element_text(size=14))


df_diff2 = readRDS('temp/sdr.RDS') %>%
  filter(Year %in% 2015:2022 & Country %in% c("CZE", "EST", "POL", "ROU") & CauseCategory == 'All') %>%
  mutate(Sex = ifelse(Sex == 'm', 'Men', 'Women'),sdr_diff = 100*(raw_sdr-smoothed_sdr)/raw_sdr) %>%
  group_by(Country, RegionCode, Year, Sex) %>%
  summarise(sdr_est = median(sdr_diff, na.rm = TRUE),
    # 95% percentile-based CI
    sdr_low = quantile(sdr_diff, probs = 0.025, na.rm = TRUE),
    sdr_up = quantile(sdr_diff, probs = 0.975, na.rm = TRUE)) %>%
  mutate(annotation = ifelse(Year == 2015 & RegionCode == 'CZ010' & Sex == 'Men', "COVID-19\npandemic", ''))

supp_table9 = readRDS('temp/sdr.RDS') %>%
  filter(Year %in% 2015:2022 & Country %in% c("CZE", "EST", "POL", "ROU") & CauseCategory == 'All') %>%
  mutate(Sex = ifelse(Sex == 'm', 'Men', 'Women'),sdr_diff = 100*(raw_sdr-smoothed_sdr)/raw_sdr) %>%
  group_by(Country, Year, Sex) %>%
  summarise(sdr_est = median(sdr_diff, na.rm = TRUE),
            # 95% percentile-based CI
            sdr_low = quantile(sdr_diff, probs = 0.025, na.rm = TRUE),
            sdr_up = quantile(sdr_diff, probs = 0.975, na.rm = TRUE),
            `Median relative difference in ASDR (95% CI)` = sprintf("%.1f (%.1f, %.1f)", sdr_est, sdr_low, sdr_up)) %>%
  select(-c(sdr_est, sdr_low, sdr_up))
  

p4a = ggplot(df_diff2, aes(x = Year)) +
  annotate("rect", xmin = 2020, xmax = Inf, ymin = -Inf,   ymax = Inf,fill = "lightgrey", alpha = 0.5)+  
  annotate("segment", x = -Inf, xend = Inf, y = 0, yend = 0, color = "black", linewidth = 0.5)+
  geom_text(aes(x = 2020.1, y = -Inf, label = annotation, hjust = 0, vjust = -0.5), size = 2)+
  geom_line(aes(y = sdr_est, color = Country,group = interaction(RegionCode, Sex)), alpha = 0.5, linewidth = 0.75) +
  geom_line(aes(y = sdr_low, color = Country,group = interaction(RegionCode, Sex)), alpha = 0.3, linewidth = 0.4) +
  geom_line(aes(y = sdr_up, color = Country,group = interaction(RegionCode, Sex)), alpha = 0.3, linewidth = 0.4) +
    facet_grid(Sex ~ Country) +
  ylab('Rel. difference in forecast vs observed ASDR (95 CI%)')+
  xlab('')+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        legend.position = 'none',
        text=element_text(size=14),
        theme(legend.title = element_blank()))



### Dispersion

df = readRDS('temp/sdr.RDS') %>%
  filter(Year %in% 2015:2022 & Country %in% c("CZE", "EST", "POL", "ROU") & CauseCategory == 'All') %>%
  mutate(Sex = ifelse(Sex == 'm', 'Men', 'Women')) %>%
  ungroup() %>%
  group_by(Country, Sex, Year, CauseCategory, Iteration) %>%
  summarise(SD = sd(smoothed_sdr),
            CV = SD/mean(smoothed_sdr),
            Theil = theil_T(smoothed_sdr))

set.seed(9999) 

df_sd_boot <- readRDS('temp/sdr.RDS') %>%
  filter(Year %in% 2015:2022, Country %in% c("CZE", "EST", "POL", "ROU"), CauseCategory == 'All') %>%
  mutate(Sex = ifelse(Sex == 'm', 'Men', 'Women')) %>%
  group_by(Country, Sex, Year) %>%
  group_modify(~ {
    n <- nrow(.x)
    # for each bootstrap draw, compute sd and cv
    stats_mat <- replicate(
      100,
      {
        samp   <- sample(.x$raw_sdr, size = n, replace = TRUE)
        sd_val <- sd(samp, na.rm = TRUE)
        cv_val <- sd_val / mean(samp, na.rm = TRUE)
        T_val = theil_T(samp, na.rm = TRUE)
        c(sd = sd_val, cv = cv_val, Theil = T_val)
      },
      simplify = "matrix"
    )
    # turn into a tibble
    tibble(
      Iteration = seq_len(100),
      SD_boot    = stats_mat["sd", ],
      CV_boot    = stats_mat["cv", ],
      Theil_boot    = stats_mat["Theil", ]
    )
  }) %>%
  ungroup()

# df_combined = left_join(df, df_sd_boot) %>%
#   group_by(Country, Year, Sex) %>%
#   summarise(
#     SD_est = median(SD, na.rm = TRUE),
#     SD_est_low = quantile(SD, probs = 0.025, na.rm = TRUE),
#     SD_est_up = quantile(SD, probs = 0.975, na.rm = TRUE),
#     SD_obs = median(SD_boot, na.rm = TRUE),
#     SD_obs_low = quantile(SD_boot, probs = 0.025, na.rm = TRUE),
#     SD_obs_up = quantile(SD_boot, probs = 0.975, na.rm = TRUE),
#     CV_est = median(CV, na.rm = TRUE),
#     CV_est_low = quantile(CV, probs = 0.025, na.rm = TRUE),
#     CV_est_up = quantile(CV, probs = 0.975, na.rm = TRUE),
#     CV_obs = median(CV_boot, na.rm = TRUE),
#     CV_obs_low = quantile(CV_boot, probs = 0.025, na.rm = TRUE),
#     CV_obs_up = quantile(CV_boot, probs = 0.975, na.rm = TRUE)) %>%
#   ungroup()
# 
# 
# ggplot(df_combined, aes(x = Year, color = Country, linetype = Sex, shape = Sex, group = interaction(Sex))) +
#   geom_line(aes(y = SD_est), alpha = 0.5, linewidth = 0.75) +
#   geom_line(aes(y = SD_est_low), alpha = 0.3, linewidth = 0.4) +
#   geom_line(aes(y = SD_est_up), alpha = 0.3, linewidth = 0.4) +
#   geom_pointrange(aes(y = SD_obs, ymin = SD_obs_low, ymax = SD_obs_up)) +
#   facet_grid(Sex ~ Country) +
#   ylab('Standard deviation')+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 90),
#         legend.position = 'top',
#         text=element_text(size=14))
# 
# ggplot(df_combined, aes(x = Year, color = Country, linetype = Sex, shape = Sex, group = interaction(Sex))) +

df_combined2 = left_join(df, df_sd_boot) %>%
  mutate(SD_diff = 100*(SD_boot-SD)/SD_boot,
         CV_diff = 100*(CV_boot-CV)/CV_boot,
         Theil_diff = 100*(Theil_boot-Theil)/Theil_boot) %>%
  group_by(Country, Year, Sex) %>%
  summarise(
    SD_diff_est = median(SD_diff, na.rm = TRUE),
    SD_diff_low = quantile(SD_diff, probs = 0.025, na.rm = TRUE),
    SD_diff_up = quantile(SD_diff, probs = 0.975, na.rm = TRUE),
    CV_diff_est = median(CV_diff, na.rm = TRUE),
    CV_diff_low = quantile(CV_diff, probs = 0.025, na.rm = TRUE),
    CV_diff_up = quantile(CV_diff, probs = 0.975, na.rm = TRUE),
    Theil_diff_est = median(Theil_diff, na.rm = TRUE),
    Theil_diff_low = quantile(Theil_diff, probs = 0.025, na.rm = TRUE),
    Theil_diff_up = quantile(Theil_diff, probs = 0.975, na.rm = TRUE)) %>%
  ungroup() %>% 
  mutate(annotation = ifelse(Year == 2015 & Country == 'CZE' & Sex == 'Men', "COVID-19\npandemic", ''))

supp_table9_1 = df_combined2 %>%
  mutate(`Median relative difference in SD (95% CI)` = 
           sprintf("%.1f (%.1f, %.1f)", SD_diff_est, SD_diff_low, SD_diff_up),
         `Median relative difference in CV (95% CI)` = 
           sprintf("%.1f (%.1f, %.1f)", CV_diff_est, CV_diff_low, CV_diff_up)) %>%
  select(Country, Year, Sex, `Median relative difference in SD (95% CI)`,  `Median relative difference in CV (95% CI)`)


p4b =ggplot(df_combined2, aes(x = Year)) +
  annotate("rect", xmin = 2020, xmax = Inf, ymin = -Inf,   ymax = Inf,fill = "lightgrey", alpha = 0.5)+ 
  annotate("segment", x = -Inf, xend = Inf, y = 0, yend = 0, color = "black", linewidth = 0.5)+
  # geom_text(aes(x = 2020.2, y = -Inf, label = annotation, hjust = 0, vjust = -1.5))+
  geom_line(aes(y = SD_diff_est, color = Country, linetype = "Standard deviation")) +
  geom_line(aes(y = SD_diff_low, color = Country, linetype = "Standard deviation"), linewidth = 0.4, alpha = 0.5) +
  geom_line(aes(y = SD_diff_up, color = Country, linetype = "Standard deviation"), linewidth = 0.4, alpha = 0.5) +
  geom_line(aes(y = CV_diff_est, color = Country, linetype = "Coefficient of Variation")) +
  geom_line(aes(y = CV_diff_low, color = Country, linetype = "Coefficient of Variation"),linewidth = 0.4, alpha = 0.5) +
  geom_line(aes(y = CV_diff_up, color = Country, linetype = "Coefficient of Variation"),linewidth = 0.4, alpha = 0.5) +
  # geom_line(aes(y = Theil_diff_est, color = Country, linetype = "Theil's T")) +
  # geom_line(aes(y = Theil_diff_low, color = Country, linetype = "Theil's T"),linewidth = 0.4, alpha = 0.5) +
  # geom_line(aes(y = Theil_diff_up, color = Country, linetype = "Theil's T"),linewidth = 0.4, alpha = 0.5) +
  facet_grid(Sex ~ Country) +
  ylab('Rel. difference in forecast vs observed SD and CoV (95% CI)')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = 'top',
        text=element_text(size=14),
        legend.title = element_blank())

p4 <- plot_grid(
  p4a, p4b,
  ncol      = 1,            # 1 column = vertical stacking
  align     = "v",          # align vertically
  labels    = c("a", "b"),  # captions for each panel
  label_size= 14,           # adjust size as needed
  label_fontface = "bold"   # optional: make labels bold
)

ggsave('plots/fig4.pdf', p4, height = 297, width = 210, units = 'mm')
ggsave('plots/fig4.png', p4, height = 297, width = 210, units = 'mm')



## LAU

df_diff2_lau = readRDS('temp/sdr_LAU.RDS') %>%
  filter(Year %in% 2015:2022 & Country %in% c("EST") & CauseCategory == 'All') %>%
  mutate(Sex = ifelse(Sex == 'm', 'Men', 'Women'),sdr_diff = 100*(raw_sdr-smoothed_sdr)/raw_sdr) %>%
  group_by(Country, RegionCode, Year, Sex) %>%
  summarise(sdr_est = median(sdr_diff, na.rm = TRUE),
            # 95% percentile-based CI
            sdr_low = quantile(sdr_diff, probs = 0.025, na.rm = TRUE),
            sdr_up = quantile(sdr_diff, probs = 0.975, na.rm = TRUE)) %>%
  mutate(annotation = ifelse(Year == 2015 & RegionCode == 'EE37' & Sex == 'Men', "COVID-19\npandemic", ''))

supp_table9a = readRDS('temp/sdr_LAU.RDS') %>%
  filter(Year %in% 2015:2022 & Country %in% c("EST") & CauseCategory == 'All') %>%
  mutate(Sex = ifelse(Sex == 'm', 'Men', 'Women'),sdr_diff = 100*(raw_sdr-smoothed_sdr)/raw_sdr) %>%
  group_by(Country, Year, Sex) %>%
  summarise(sdr_est = median(sdr_diff, na.rm = TRUE),
            # 95% percentile-based CI
            sdr_low = quantile(sdr_diff, probs = 0.025, na.rm = TRUE),
            sdr_up = quantile(sdr_diff, probs = 0.975, na.rm = TRUE),
            `Median relative difference in ASDR (95% CI)` = sprintf("%.1f (%.1f, %.1f)", sdr_est, sdr_low, sdr_up),
            Country = 'EST - LAU') %>%
  select(-c(sdr_est, sdr_low, sdr_up))


supp_fig3a = ggplot(df_diff2_lau, aes(x = Year)) +
  annotate("rect", xmin = 2020, xmax = Inf, ymin = -Inf,   ymax = Inf,fill = "lightgrey", alpha = 0.5)+  
  annotate("segment", x = -Inf, xend = Inf, y = 0, yend = 0, color = "black", linewidth = 0.5)+
  geom_text(aes(x = 2020.2, y = -Inf, label = annotation, hjust = 0, vjust = -0.25))+
  geom_line(aes(y = sdr_est, color = Country,group = interaction(RegionCode, Sex)), alpha = 0.5, linewidth = 0.75) +
  geom_line(aes(y = sdr_low, color = Country,group = interaction(RegionCode, Sex)), alpha = 0.3, linewidth = 0.4) +
  geom_line(aes(y = sdr_up, color = Country,group = interaction(RegionCode, Sex)), alpha = 0.3, linewidth = 0.4) +
  scale_color_manual(values = c('EST' = '#B79F00',
                                'LTU' = '#00BA38',
                                'SVK' = '#F564E3'))+
  facet_grid(Sex ~ Country) +
  ylab('Rel. difference in forecast vs observed ASDR (95 CI%)')+
  xlab('')+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        legend.position = 'none',
        text=element_text(size=14),
        theme(legend.title = element_blank()))



### Dispersion

df_lau = readRDS('temp/sdr_lau.RDS') %>%
  filter(Year %in% 2015:2022 & Country %in% c("EST") & CauseCategory == 'All') %>%
  mutate(Sex = ifelse(Sex == 'm', 'Men', 'Women')) %>%
  ungroup() %>%
  group_by(Country, Sex, Year, CauseCategory, Iteration) %>%
  summarise(SD = sd(smoothed_sdr),
            CV = SD/mean(smoothed_sdr),
            Theil = theil_T(smoothed_sdr))

set.seed(9999) 

df_sd_boot_lau <- readRDS('temp/sdr.RDS') %>%
  filter(Year %in% 2015:2022 & Country %in% c("EST") & CauseCategory == 'All') %>%
  mutate(Sex = ifelse(Sex == 'm', 'Men', 'Women')) %>%
  group_by(Country, Sex, Year) %>%
  group_modify(~ {
    n <- nrow(.x)
    # for each bootstrap draw, compute sd and cv
    stats_mat <- replicate(
      100,
      {
        samp   <- sample(.x$raw_sdr, size = n, replace = TRUE)
        sd_val <- sd(samp, na.rm = TRUE)
        cv_val <- sd_val / mean(samp, na.rm = TRUE)
        T_val = theil_T(samp, na.rm = TRUE)
        c(sd = sd_val, cv = cv_val, Theil = T_val)
      },
      simplify = "matrix"
    )
    # turn into a tibble
    tibble(
      Iteration = seq_len(100),
      SD_boot    = stats_mat["sd", ],
      CV_boot    = stats_mat["cv", ],
      Theil_boot    = stats_mat["Theil", ]
    )
  }) %>%
  ungroup()

df_combined2_lau = left_join(df_lau, df_sd_boot_lau) %>%
  mutate(SD_diff = 100*(SD_boot-SD)/SD_boot,
         CV_diff = 100*(CV_boot-CV)/CV_boot,
         Theil_diff = 100*(Theil_boot-Theil)/Theil_boot) %>%
  group_by(Country, Year, Sex) %>%
  summarise(
    SD_diff_est = median(SD_diff, na.rm = TRUE),
    SD_diff_low = quantile(SD_diff, probs = 0.025, na.rm = TRUE),
    SD_diff_up = quantile(SD_diff, probs = 0.975, na.rm = TRUE),
    CV_diff_est = median(CV_diff, na.rm = TRUE),
    CV_diff_low = quantile(CV_diff, probs = 0.025, na.rm = TRUE),
    CV_diff_up = quantile(CV_diff, probs = 0.975, na.rm = TRUE),
    Theil_diff_est = median(Theil_diff, na.rm = TRUE),
    Theil_diff_low = quantile(Theil_diff, probs = 0.025, na.rm = TRUE),
    Theil_diff_up = quantile(Theil_diff, probs = 0.975, na.rm = TRUE)) %>%
  ungroup() %>% 
  mutate(annotation = ifelse(Year == 2015 & Country == 'EST' & Sex == 'Men', "COVID-19\npandemic", ''))


supp_table9a_1 = df_combined2_lau %>%
  mutate(`Median relative difference in SD (95% CI)` = 
           sprintf("%.1f (%.1f, %.1f)", SD_diff_est, SD_diff_low, SD_diff_up),
         `Median relative difference in CV (95% CI)` = 
           sprintf("%.1f (%.1f, %.1f)", CV_diff_est, CV_diff_low, CV_diff_up),
         Country = 'EST - LAU') %>%
  select(Country, Year, Sex, `Median relative difference in SD (95% CI)`,  `Median relative difference in CV (95% CI)`)


supp_fig3b =ggplot(df_combined2_lau, aes(x = Year)) +
  annotate("rect", xmin = 2020, xmax = Inf, ymin = -Inf,   ymax = Inf,fill = "lightgrey", alpha = 0.5)+  
  annotate("segment", x = -Inf, xend = Inf, y = 0, yend = 0, color = "black", linewidth = 0.5)+
  geom_line(aes(y = SD_diff_est, color = Country, linetype = "Standard deviation")) +
  geom_line(aes(y = SD_diff_low, color = Country, linetype = "Standard deviation"), linewidth = 0.4, alpha = 0.5) +
  geom_line(aes(y = SD_diff_up, color = Country, linetype = "Standard deviation"), linewidth = 0.4, alpha = 0.5) +
  geom_line(aes(y = CV_diff_est, color = Country, linetype = "Coefficient of Variation")) +
  geom_line(aes(y = CV_diff_low, color = Country, linetype = "Coefficient of Variation"),linewidth = 0.4, alpha = 0.5) +
  geom_line(aes(y = CV_diff_up, color = Country, linetype = "Coefficient of Variation"),linewidth = 0.4, alpha = 0.5) +
  # geom_line(aes(y = Theil_diff_est, color = Country, linetype = "Theil's T")) +
  # geom_line(aes(y = Theil_diff_low, color = Country, linetype = "Theil's T"),linewidth = 0.4, alpha = 0.5) +
  # geom_line(aes(y = Theil_diff_up, color = Country, linetype = "Theil's T"),linewidth = 0.4, alpha = 0.5) +
  scale_color_manual(values = c('EST' = '#B79F00',
                                'LTU' = '#00BA38',
                                'SVK' = '#F564E3'))+
  facet_grid(Sex ~ Country) +
  ylab('Rel. difference in forecast vs observed SD and CoV (95% CI)')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = 'top',
        text=element_text(size=14),
        legend.title = element_blank())

supp_fig3 <- plot_grid(
  supp_fig3a, supp_fig3b,
  ncol      = 1,            # 1 column = vertical stacking
  align     = "v",          # align vertically
  labels    = c("a", "b"),  # captions for each panel
  label_size= 14,           # adjust size as needed
  label_fontface = "bold"   # optional: make labels bold
)

ggsave('plots/supp_fig3.png', supp_fig3, height = 297, width = 210, units = 'mm')

supp_table9 = left_join(supp_table9, supp_table9_1)
supp_table9a = left_join(supp_table9a, supp_table9a_1)
supp_table9 = rbind(supp_table9, supp_table9a)

ft <- supp_table9 %>%
  arrange(Country, Sex, Year) %>%
  mutate(Year = as.character(Year)) %>%
  flextable(.)%>%
  set_table_properties(layout = "autofit") %>%
  set_caption("Supplementary Table 9: Relative difference between observed and forecast ASDR, SD, and CV") %>%
  fontsize(size = 8, part = 'all')

ft

save_as_docx(ft, path = 'tables/Supp_table9.docx')
