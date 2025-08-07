# Subnational mortality convergence in Europe by cause of death
# Rok Hrzic (r.hrzic (at) maastrichtuniversity.nl)
# May 2025

rm(list = ls())

require(tidyverse)
require(zoo)
require(cowplot)
require(flextable)
require(officer)

theil_T <- function(x, na.rm = TRUE) {
  x <- if (na.rm) x[!is.na(x)] else x
  mu <- mean(x)
  p  <- x / mu
  mean(p * log(p))
}

df = readRDS('temp/sdr.RDS') %>%
  filter(Year %in% 2000:2019) %>%
  mutate(Sex = ifelse(Sex == 'm', 'Men', 'Women')) %>%
  ungroup() %>%
  group_by(Country, Sex, Year, CauseCategory, Iteration) %>%
  summarise(SD = sd(smoothed_sdr),
            CV = SD/mean(smoothed_sdr),
            Theil = theil_T(smoothed_sdr)) %>%
  group_by(Country, Year, Sex, CauseCategory) %>%
  summarise(
    SD_est = median(SD, na.rm = TRUE),
    SD_low = quantile(SD, probs = 0.025, na.rm = TRUE),
    SD_up = quantile(SD, probs = 0.975, na.rm = TRUE),
    SD_se = sd(CV, na.rm = TRUE),
    CV_est = median(CV, na.rm = TRUE),
    CV_low = quantile(CV, probs = 0.025, na.rm = TRUE),
    CV_up = quantile(CV, probs = 0.975, na.rm = TRUE),
    CV_se = sd(CV, na.rm = TRUE),
    Theil_est = median(Theil, na.rm = TRUE),
    Theil_low = quantile(Theil, probs = 0.025, na.rm = TRUE),
    Theil_up = quantile(Theil, probs = 0.975, na.rm = TRUE),
    Theil_se = sd(Theil, na.rm = TRUE)) %>%
  ungroup()

source("code/segmented.R")

## standard deviation

df_segmented_sd <- df %>%
  mutate(logSD = log(SD_est),
         logSD_se = SD_se/SD_est) %>%
  group_by(Country, Sex, CauseCategory) %>%
  group_modify(~ fit_segmented_boot(.x, response = "logSD", se = 'logSD_se')) %>%
  ungroup()

saveRDS(df_segmented_sd, 'temp/fig2segmented_sd.RDS')

df_segmented_sd = readRDS('temp/fig2segmented_sd.RDS')

df_bp_sd <- df_segmented_sd %>%
  group_by(Country, Sex, CauseCategory) %>%
  slice(1) %>%
  ungroup()


p3a = ggplot(filter(df_segmented_sd, CauseCategory != 'Other'), aes(x = Year, color = Country, linetype = Sex, shape = Sex, group = interaction(Country, Sex))) +
  geom_line(aes(y = SD_est), linewidth = 1) +
  geom_line(aes(y = SD_low), alpha = 0.6, linewidth = 0.4) +
  geom_line(aes(y = SD_up), alpha = 0.6, linewidth = 0.4) +
  geom_line(aes(y = exp(fitted)), color = 'black', linewidth = 0.7) +
  facet_grid(CauseCategory ~ Country, scales="free_y") +
  ylab('Standard deviation')+
  xlab('')+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        legend.position = 'top',
        text=element_text(size=14))


supp_table5 <- df_bp_sd %>%
  dplyr::select(Country, Sex, CauseCategory,
    bp1_est, bp1_lower, bp1_upper,
    bp2_est, bp2_lower, bp2_upper,
    bp3_est, bp3_lower, bp3_upper,
    AAPC1, AAPC1_lower, AAPC1_upper,
    AAPC2, AAPC2_lower, AAPC2_upper,
    AAPC3, AAPC3_lower, AAPC3_upper,
    AAPC4, AAPC4_lower, AAPC4_upper,
    total_pct_change, total_pct_change_lower, total_pct_change_upper,
    AAPC_total, AAPC_total_lower, AAPC_total_upper) %>%
  dplyr::mutate(
    sigA1 = AAPC1_lower * AAPC1_upper > 0,
    sigA2 = AAPC2_lower * AAPC2_upper > 0,
    sigA3 = AAPC3_lower * AAPC3_upper > 0,
    sigA4 = AAPC4_lower * AAPC4_upper > 0,
    sigAT  = AAPC_total_lower * AAPC_total_upper > 0,
    sigT  = total_pct_change_lower * total_pct_change_upper > 0,
    `Breakpoint 1 (95% CI)` = ifelse(is.na(bp1_est), "-", sprintf("%.0f (%.1f, %.1f)", bp1_est, bp1_lower, bp1_upper)),
    `Breakpoint 2 (95% CI)` = ifelse(is.na(bp2_est), "-", sprintf("%.0f (%.1f, %.1f)", bp2_est, bp2_lower, bp2_upper)),
    `Breakpoint 3 (95% CI)` = ifelse(is.na(bp3_est), "-", sprintf("%.0f (%.1f, %.1f)", bp3_est, bp3_lower, bp3_upper)),
    `AAPC 1 (%) (95% CI)` = if_else(is.na(AAPC1), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC1, if_else(sigA1, "*", ""), AAPC1_lower, AAPC1_upper)),
    `AAPC 2 (%) (95% CI)` = if_else(is.na(AAPC2), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC2, if_else(sigA2, "*", ""), AAPC2_lower, AAPC2_upper)),
    `AAPC 3 (%) (95% CI)` = if_else(is.na(AAPC3), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC3,if_else(sigA3, "*", ""), AAPC3_lower, AAPC3_upper)),
    `AAPC 4 (%) (95% CI)` = if_else(is.na(AAPC4), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC4,if_else(sigA4, "*", ""), AAPC4_lower, AAPC4_upper)),
    `AAPC Total (%) (95% CI)` = if_else(is.na(AAPC_total), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC_total,if_else(sigAT, "*", ""), AAPC_total_lower, AAPC_total_upper)),
    `Cumulative change (%) (95% CI)` = if_else(is.na(total_pct_change), "-", sprintf("%.1f%s (%.1f, %.1f)", total_pct_change, if_else(sigT, "*", ""), total_pct_change_lower, total_pct_change_upper)),
  ) %>%
  dplyr::select(Country, Sex, CauseCategory,
                `Breakpoint 1 (95% CI)`,
                `Breakpoint 2 (95% CI)`,
                `Breakpoint 3 (95% CI)`,
                `AAPC 1 (%) (95% CI)`,
                `AAPC 2 (%) (95% CI)`,
                `AAPC 3 (%) (95% CI)`,
                `AAPC 4 (%) (95% CI)`,
                `AAPC Total (%) (95% CI)`,
                `Cumulative change (%) (95% CI)`)


# Create the flextable
ft <- flextable(supp_table5) %>%
  set_table_properties(layout = "autofit") %>%
  set_caption("Supplementary Table 5: Joinpoint regression of standard deviation trends") %>%
  fontsize(size = 8, part = 'all')

ft

# Export the table to a Word document using officer
save_as_docx(ft, path = 'tables/Supp_table5.docx', pr_section = prop_section(page_size(orient = 'landscape')))



## CV

df_segmented_cv <- df %>%
  mutate(logCV = log(CV_est),
         logCV_se = CV_se/CV_est) %>%
  group_by(Country, Sex, CauseCategory) %>%
  group_modify(~ fit_segmented_boot(.x, response = "logCV", se = 'logCV_se')) %>%
  ungroup()

saveRDS(df_segmented_cv, 'temp/fig2segmented_cv.RDS')

df_bp_cv <- df_segmented_cv %>%
  group_by(Country, Sex, CauseCategory) %>%
  slice(1) %>%
  ungroup()


p3b = ggplot(filter(df_segmented_cv, CauseCategory != 'Other'), aes(x = Year, color = Country, linetype = Sex, shape = Sex, group = interaction(Country, Sex))) +
  geom_line(aes(y = CV_est), linewidth = 1) +
  geom_line(aes(y = CV_low), alpha = 0.6, linewidth = 0.4) +
  geom_line(aes(y = CV_up), alpha = 0.6, linewidth = 0.4) +
  geom_line(aes(y = exp(fitted)), color = 'black', linewidth = 0.7) +
  facet_grid(CauseCategory ~ Country, scales="free_y") +
  ylab('Coefficient of variation')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = 'none',
        text=element_text(size=14))

p3 <- plot_grid(
  p3a, p3b,
  ncol      = 1,            # 1 column = vertical stacking
  align     = "v",          # align vertically
  labels    = c("a", "b"),  # captions for each panel
  label_size= 14,           # adjust size as needed
  label_fontface = "bold"   # optional: make labels bold
)


ggsave('plots/fig3.pdf', p3, height = 297, width = 210, units = 'mm')
ggsave('plots/fig3.png', p3, height = 297, width = 210, units = 'mm')


supp_table6 <- df_bp_cv %>%
  dplyr::select(Country, Sex, CauseCategory,
                bp1_est, bp1_lower, bp1_upper,
                bp2_est, bp2_lower, bp2_upper,
                bp3_est, bp3_lower, bp3_upper,
                AAPC1, AAPC1_lower, AAPC1_upper,
                AAPC2, AAPC2_lower, AAPC2_upper,
                AAPC3, AAPC3_lower, AAPC3_upper,
                AAPC4, AAPC4_lower, AAPC4_upper,
                total_pct_change, total_pct_change_lower, total_pct_change_upper,
                AAPC_total, AAPC_total_lower, AAPC_total_upper) %>%
  dplyr::mutate(
    sigA1 = AAPC1_lower * AAPC1_upper > 0,
    sigA2 = AAPC2_lower * AAPC2_upper > 0,
    sigA3 = AAPC3_lower * AAPC3_upper > 0,
    sigA4 = AAPC4_lower * AAPC4_upper > 0,
    sigAT  = AAPC_total_lower * AAPC_total_upper > 0,
    sigT  = total_pct_change_lower * total_pct_change_upper > 0,
    `Breakpoint 1 (95% CI)` = ifelse(is.na(bp1_est), "-", sprintf("%.0f (%.1f, %.1f)", bp1_est, bp1_lower, bp1_upper)),
    `Breakpoint 2 (95% CI)` = ifelse(is.na(bp2_est), "-", sprintf("%.0f (%.1f, %.1f)", bp2_est, bp2_lower, bp2_upper)),
    `Breakpoint 3 (95% CI)` = ifelse(is.na(bp3_est), "-", sprintf("%.0f (%.1f, %.1f)", bp3_est, bp3_lower, bp3_upper)),
    `AAPC 1 (%) (95% CI)` = if_else(is.na(AAPC1), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC1, if_else(sigA1, "*", ""), AAPC1_lower, AAPC1_upper)),
    `AAPC 2 (%) (95% CI)` = if_else(is.na(AAPC2), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC2, if_else(sigA2, "*", ""), AAPC2_lower, AAPC2_upper)),
    `AAPC 3 (%) (95% CI)` = if_else(is.na(AAPC3), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC3,if_else(sigA3, "*", ""), AAPC3_lower, AAPC3_upper)),
    `AAPC 4 (%) (95% CI)` = if_else(is.na(AAPC4), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC4,if_else(sigA4, "*", ""), AAPC4_lower, AAPC4_upper)),
    `AAPC Total (%) (95% CI)` = if_else(is.na(AAPC_total), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC_total,if_else(sigAT, "*", ""), AAPC_total_lower, AAPC_total_upper)),
    `Cumulative change (%) (95% CI)` = if_else(is.na(total_pct_change), "-", sprintf("%.1f%s (%.1f, %.1f)", total_pct_change, if_else(sigT, "*", ""), total_pct_change_lower, total_pct_change_upper)),
  ) %>%
  dplyr::select(Country, Sex, CauseCategory,
                `Breakpoint 1 (95% CI)`,
                `Breakpoint 2 (95% CI)`,
                `Breakpoint 3 (95% CI)`,
                `AAPC 1 (%) (95% CI)`,
                `AAPC 2 (%) (95% CI)`,
                `AAPC 3 (%) (95% CI)`,
                `AAPC 4 (%) (95% CI)`,
                `AAPC Total (%) (95% CI)`,
                `Cumulative change (%) (95% CI)`)


# Create the flextable
ft <- flextable(supp_table6) %>%
  set_table_properties(layout = "autofit") %>%
  set_caption("Supplementary Table 6: Joinpoint regression of trends coefficient of variation trends") %>%
  fontsize(size = 8, part = 'all')

ft

# Export the table to a Word document using officer
save_as_docx(ft, path = 'tables/Supp_table6.docx', pr_section = prop_section(page_size(orient = 'landscape')))

# 
# interpretation = left_join(dplyr::select(supp_table5, Country, Sex, CauseCategory, Change_SD = `AAPC Total (%) (95% CI)`),
#                             dplyr::select(supp_table6, Country, Sex, CauseCategory, Change_CV = `AAPC Total (%) (95% CI)`)) %>%
#   mutate(Change_SD_num = as.numeric(str_extract(Change_SD, "^[^*\\s]+")),
#          Change_CV_num = as.numeric(str_extract(Change_CV, "^[^*\\s]+")),
#          aligned = ifelse(Change_SD_num*Change_CV_num > 0, 'Aligned', 'Opposite')) %>%
#   arrange(CauseCategory, Sex)
# 
# saveRDS(interpretation, 'temp/interpretation.RDS')
# 
# interpretation = readRDS('temp/interpretation.RDS')
# 
# intepretation2 = left_join(interpretation, dplyr::select(supp_table7, Country, Sex, CauseCategory, Change_T = `AAPC Total (%) (95% CI)`)) %>%
#   mutate(Change_T_num = as.numeric(str_extract(Change_T, "^[^*\\s]+")),
#          aligned2 = ifelse(Change_CV_num*Change_CV_num>0, 'Aligned', 'Opposite')) %>%
#   arrange(CauseCategory, Sex)


## Theil index sens

df_segmented_Theil <- df %>%
  mutate(logT = log(Theil_est),
         logTheil_se = Theil_se/Theil_est) %>%
  group_by(Country, Sex, CauseCategory) %>%
  group_modify(~ fit_segmented_boot(.x, response = "logT", se = 'logTheil_se')) %>%
  ungroup()

saveRDS(df_segmented_Theil, 'temp/fig2segmented_Theil.RDS')

df_bp_Theil <- df_segmented_Theil %>%
  group_by(Country, Sex, CauseCategory) %>%
  slice(1) %>%
  ungroup()

supp_table7 <- df_bp_Theil %>%
  dplyr::select(Country, Sex, CauseCategory,
                bp1_est, bp1_lower, bp1_upper,
                bp2_est, bp2_lower, bp2_upper,
                bp3_est, bp3_lower, bp3_upper,
                AAPC1, AAPC1_lower, AAPC1_upper,
                AAPC2, AAPC2_lower, AAPC2_upper,
                AAPC3, AAPC3_lower, AAPC3_upper,
                AAPC4, AAPC4_lower, AAPC4_upper,
                total_pct_change, total_pct_change_lower, total_pct_change_upper,
                AAPC_total, AAPC_total_lower, AAPC_total_upper) %>%
  dplyr::mutate(
    sigA1 = AAPC1_lower * AAPC1_upper > 0,
    sigA2 = AAPC2_lower * AAPC2_upper > 0,
    sigA3 = AAPC3_lower * AAPC3_upper > 0,
    sigA4 = AAPC4_lower * AAPC4_upper > 0,
    sigAT  = AAPC_total_lower * AAPC_total_upper > 0,
    sigT  = total_pct_change_lower * total_pct_change_upper > 0,
    `Breakpoint 1 (95% CI)` = ifelse(is.na(bp1_est), "-", sprintf("%.0f (%.1f, %.1f)", bp1_est, bp1_lower, bp1_upper)),
    `Breakpoint 2 (95% CI)` = ifelse(is.na(bp2_est), "-", sprintf("%.0f (%.1f, %.1f)", bp2_est, bp2_lower, bp2_upper)),
    `Breakpoint 3 (95% CI)` = ifelse(is.na(bp3_est), "-", sprintf("%.0f (%.1f, %.1f)", bp3_est, bp3_lower, bp3_upper)),
    `AAPC 1 (%) (95% CI)` = if_else(is.na(AAPC1), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC1, if_else(sigA1, "*", ""), AAPC1_lower, AAPC1_upper)),
    `AAPC 2 (%) (95% CI)` = if_else(is.na(AAPC2), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC2, if_else(sigA2, "*", ""), AAPC2_lower, AAPC2_upper)),
    `AAPC 3 (%) (95% CI)` = if_else(is.na(AAPC3), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC3,if_else(sigA3, "*", ""), AAPC3_lower, AAPC3_upper)),
    `AAPC 4 (%) (95% CI)` = if_else(is.na(AAPC4), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC4,if_else(sigA4, "*", ""), AAPC4_lower, AAPC4_upper)),
    `AAPC Total (%) (95% CI)` = if_else(is.na(AAPC_total), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC_total,if_else(sigAT, "*", ""), AAPC_total_lower, AAPC_total_upper)),
    `Cumulative change (%) (95% CI)` = if_else(is.na(total_pct_change), "-", sprintf("%.1f%s (%.1f, %.1f)", total_pct_change, if_else(sigT, "*", ""), total_pct_change_lower, total_pct_change_upper)),
  ) %>%
  dplyr::select(Country, Sex, CauseCategory,
                `Breakpoint 1 (95% CI)`,
                `Breakpoint 2 (95% CI)`,
                `Breakpoint 3 (95% CI)`,
                `AAPC 1 (%) (95% CI)`,
                `AAPC 2 (%) (95% CI)`,
                `AAPC 3 (%) (95% CI)`,
                `AAPC 4 (%) (95% CI)`,
                `AAPC Total (%) (95% CI)`,
                `Cumulative change (%) (95% CI)`)


# Create the flextable
ft <- flextable(supp_table7) %>%
  set_table_properties(layout = "autofit") %>%
  set_caption("Supplementary Table 7: Joinpoint regression of Theil index trends") %>%
  fontsize(size = 8, part = 'all')

ft

# Export the table to a Word document using officer
# save_as_docx(ft, path = 'tables/Supp_table7.docx', pr_section = prop_section(page_size(orient = 'landscape')))

## LAU

df_lau = readRDS('temp/sdr_LAU.RDS') %>%
  filter(Year %in% 2000:2019) %>%
  mutate(Sex = ifelse(Sex == 'm', 'Men', 'Women')) %>%
  ungroup() %>%
  group_by(Country, Sex, Year, CauseCategory, Iteration) %>%
  summarise(SD = sd(smoothed_sdr),
            CV = SD/mean(smoothed_sdr),
            Theil = theil_T(smoothed_sdr)) %>%
  group_by(Country, Year, Sex, CauseCategory) %>%
  summarise(
    SD_est = median(SD, na.rm = TRUE),
    SD_low = quantile(SD, probs = 0.025, na.rm = TRUE),
    SD_up = quantile(SD, probs = 0.975, na.rm = TRUE),
    SD_se = sd(CV, na.rm = TRUE),
    CV_est = median(CV, na.rm = TRUE),
    CV_low = quantile(CV, probs = 0.025, na.rm = TRUE),
    CV_up = quantile(CV, probs = 0.975, na.rm = TRUE),
    CV_se = sd(CV, na.rm = TRUE),
    Theil_est = median(Theil, na.rm = TRUE),
    Theil_low = quantile(Theil, probs = 0.025, na.rm = TRUE),
    Theil_up = quantile(Theil, probs = 0.975, na.rm = TRUE),
    Theil_se = sd(Theil, na.rm = TRUE)) %>%
  ungroup()

source("code/segmented.R")

## standard deviation

df_segmented_sd_lau <- df_lau %>%
  mutate(logSD = log(SD_est),
         logSD_se = SD_se/SD_est) %>%
  group_by(Country, Sex, CauseCategory) %>%
  group_modify(~ fit_segmented_boot(.x, response = "logSD", se = 'logSD_se')) %>%
  ungroup()

saveRDS(df_segmented_sd_lau, 'temp/fig2segmented_sd_lau.RDS')

df_bp_sd_lau <- df_segmented_sd_lau %>%
  group_by(Country, Sex, CauseCategory) %>%
  slice(1) %>%
  ungroup()

supp_fig2a = ggplot(filter(df_segmented_sd_lau, CauseCategory != 'Other'), aes(x = Year, color = Country, linetype = Sex, shape = Sex, group = interaction(Country, Sex))) +
  geom_line(aes(y = SD_est), linewidth = 1) +
  geom_line(aes(y = SD_low), alpha = 0.6, linewidth = 0.4) +
  geom_line(aes(y = SD_up), alpha = 0.6, linewidth = 0.4) +
  geom_line(aes(y = exp(fitted)), color = 'black', linewidth = 0.7) +
  scale_color_manual(values = c('EST' = '#B79F00',
                                'LTU' = '#00BA38',
                                'SVK' = '#F564E3'))+
  facet_grid(CauseCategory ~ Country, scales="free_y") +
  ylab('Standard deviation')+
  xlab('')+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        legend.position = 'top',
        text=element_text(size=14))


supp_table8 <- df_bp_sd_lau %>%
  dplyr::select(Country, Sex, CauseCategory,
                bp1_est, bp1_lower, bp1_upper,
                bp2_est, bp2_lower, bp2_upper,
                bp3_est, bp3_lower, bp3_upper,
                AAPC1, AAPC1_lower, AAPC1_upper,
                AAPC2, AAPC2_lower, AAPC2_upper,
                AAPC3, AAPC3_lower, AAPC3_upper,
                AAPC4, AAPC4_lower, AAPC4_upper,
                total_pct_change, total_pct_change_lower, total_pct_change_upper,
                AAPC_total, AAPC_total_lower, AAPC_total_upper) %>%
  dplyr::mutate(
    sigA1 = AAPC1_lower * AAPC1_upper > 0,
    sigA2 = AAPC2_lower * AAPC2_upper > 0,
    sigA3 = AAPC3_lower * AAPC3_upper > 0,
    sigA4 = AAPC4_lower * AAPC4_upper > 0,
    sigAT  = AAPC_total_lower * AAPC_total_upper > 0,
    sigT  = total_pct_change_lower * total_pct_change_upper > 0,
    `Breakpoint 1 (95% CI)` = ifelse(is.na(bp1_est), "-", sprintf("%.0f (%.1f, %.1f)", bp1_est, bp1_lower, bp1_upper)),
    `Breakpoint 2 (95% CI)` = ifelse(is.na(bp2_est), "-", sprintf("%.0f (%.1f, %.1f)", bp2_est, bp2_lower, bp2_upper)),
    `Breakpoint 3 (95% CI)` = ifelse(is.na(bp3_est), "-", sprintf("%.0f (%.1f, %.1f)", bp3_est, bp3_lower, bp3_upper)),
    `AAPC 1 (%) (95% CI)` = if_else(is.na(AAPC1), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC1, if_else(sigA1, "*", ""), AAPC1_lower, AAPC1_upper)),
    `AAPC 2 (%) (95% CI)` = if_else(is.na(AAPC2), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC2, if_else(sigA2, "*", ""), AAPC2_lower, AAPC2_upper)),
    `AAPC 3 (%) (95% CI)` = if_else(is.na(AAPC3), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC3,if_else(sigA3, "*", ""), AAPC3_lower, AAPC3_upper)),
    `AAPC 4 (%) (95% CI)` = if_else(is.na(AAPC4), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC4,if_else(sigA4, "*", ""), AAPC4_lower, AAPC4_upper)),
    `AAPC Total (%) (95% CI)` = if_else(is.na(AAPC_total), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC_total,if_else(sigAT, "*", ""), AAPC_total_lower, AAPC_total_upper)),
    `Cumulative change (%) (95% CI)` = if_else(is.na(total_pct_change), "-", sprintf("%.1f%s (%.1f, %.1f)", total_pct_change, if_else(sigT, "*", ""), total_pct_change_lower, total_pct_change_upper)),
  ) %>%
  dplyr::select(Country, Sex, CauseCategory,
                `Breakpoint 1 (95% CI)`,
                `Breakpoint 2 (95% CI)`,
                `Breakpoint 3 (95% CI)`,
                `AAPC 1 (%) (95% CI)`,
                `AAPC 2 (%) (95% CI)`,
                `AAPC 3 (%) (95% CI)`,
                `AAPC 4 (%) (95% CI)`,
                `AAPC Total (%) (95% CI)`,
                `Cumulative change (%) (95% CI)`)


# Create the flextable
ft <- flextable(supp_table8) %>%
  set_table_properties(layout = "autofit") %>%
  set_caption("Supplementary Table 7: Joinpoint regression of standard deviation trends - LAU level") %>%
  fontsize(size = 8, part = 'all')

ft

# Export the table to a Word document using officer
save_as_docx(ft, path = 'tables/Supp_table7.docx', pr_section = prop_section(page_size(orient = 'landscape')))


## CV

df_segmented_cv_lau <- df_lau %>%
  mutate(logCV = log(CV_est),
         logCV_se = CV_se/CV_est) %>%
  group_by(Country, Sex, CauseCategory) %>%
  group_modify(~ fit_segmented_boot(.x, response = "logCV", se = 'logCV_se')) %>%
  ungroup()

saveRDS(df_segmented_cv_lau, 'temp/fig2segmented_cv_lau.RDS')

df_bp_cv_lau <- df_segmented_cv_lau %>%
  group_by(Country, Sex, CauseCategory) %>%
  slice(1) %>%
  ungroup()


supp_fig2b = ggplot(filter(df_segmented_cv_lau, CauseCategory != 'Other'), aes(x = Year, color = Country, linetype = Sex, shape = Sex, group = interaction(Country, Sex))) +
  geom_line(aes(y = CV_est), linewidth = 1) +
  geom_line(aes(y = CV_low), alpha = 0.6, linewidth = 0.4) +
  geom_line(aes(y = CV_up), alpha = 0.6, linewidth = 0.4) +
  geom_line(aes(y = exp(fitted)), color = 'black', linewidth = 0.7) +
  scale_color_manual(values = c('EST' = '#B79F00',
                                'LTU' = '#00BA38',
                                'SVK' = '#F564E3'))+
  facet_grid(CauseCategory ~ Country, scales="free_y") +
  ylab('Coefficient of variation')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = 'none',
        text=element_text(size=14))

supp_fig2 <- plot_grid(
  supp_fig2a, supp_fig2b,
  ncol      = 1,            # 1 column = vertical stacking
  align     = "v",          # align vertically
  labels    = c("a", "b"),  # captions for each panel
  label_size= 14,           # adjust size as needed
  label_fontface = "bold"   # optional: make labels bold
)

ggsave('plots/supp_fig2.png', supp_fig2, height = 297, width = 210, units = 'mm')


supp_table9 <- df_bp_cv_lau %>%
  dplyr::select(Country, Sex, CauseCategory,
                bp1_est, bp1_lower, bp1_upper,
                bp2_est, bp2_lower, bp2_upper,
                bp3_est, bp3_lower, bp3_upper,
                AAPC1, AAPC1_lower, AAPC1_upper,
                AAPC2, AAPC2_lower, AAPC2_upper,
                AAPC3, AAPC3_lower, AAPC3_upper,
                AAPC4, AAPC4_lower, AAPC4_upper,
                total_pct_change, total_pct_change_lower, total_pct_change_upper,
                AAPC_total, AAPC_total_lower, AAPC_total_upper) %>%
  dplyr::mutate(
    sigA1 = AAPC1_lower * AAPC1_upper > 0,
    sigA2 = AAPC2_lower * AAPC2_upper > 0,
    sigA3 = AAPC3_lower * AAPC3_upper > 0,
    sigA4 = AAPC4_lower * AAPC4_upper > 0,
    sigAT  = AAPC_total_lower * AAPC_total_upper > 0,
    sigT  = total_pct_change_lower * total_pct_change_upper > 0,
    `Breakpoint 1 (95% CI)` = ifelse(is.na(bp1_est), "-", sprintf("%.0f (%.1f, %.1f)", bp1_est, bp1_lower, bp1_upper)),
    `Breakpoint 2 (95% CI)` = ifelse(is.na(bp2_est), "-", sprintf("%.0f (%.1f, %.1f)", bp2_est, bp2_lower, bp2_upper)),
    `Breakpoint 3 (95% CI)` = ifelse(is.na(bp3_est), "-", sprintf("%.0f (%.1f, %.1f)", bp3_est, bp3_lower, bp3_upper)),
    `AAPC 1 (%) (95% CI)` = if_else(is.na(AAPC1), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC1, if_else(sigA1, "*", ""), AAPC1_lower, AAPC1_upper)),
    `AAPC 2 (%) (95% CI)` = if_else(is.na(AAPC2), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC2, if_else(sigA2, "*", ""), AAPC2_lower, AAPC2_upper)),
    `AAPC 3 (%) (95% CI)` = if_else(is.na(AAPC3), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC3,if_else(sigA3, "*", ""), AAPC3_lower, AAPC3_upper)),
    `AAPC 4 (%) (95% CI)` = if_else(is.na(AAPC4), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC4,if_else(sigA4, "*", ""), AAPC4_lower, AAPC4_upper)),
    `AAPC Total (%) (95% CI)` = if_else(is.na(AAPC_total), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC_total,if_else(sigAT, "*", ""), AAPC_total_lower, AAPC_total_upper)),
    `Cumulative change (%) (95% CI)` = if_else(is.na(total_pct_change), "-", sprintf("%.1f%s (%.1f, %.1f)", total_pct_change, if_else(sigT, "*", ""), total_pct_change_lower, total_pct_change_upper)),
  ) %>%
  dplyr::select(Country, Sex, CauseCategory,
                `Breakpoint 1 (95% CI)`,
                `Breakpoint 2 (95% CI)`,
                `Breakpoint 3 (95% CI)`,
                `AAPC 1 (%) (95% CI)`,
                `AAPC 2 (%) (95% CI)`,
                `AAPC 3 (%) (95% CI)`,
                `AAPC 4 (%) (95% CI)`,
                `AAPC Total (%) (95% CI)`,
                `Cumulative change (%) (95% CI)`)



ft <- flextable(supp_table9) %>%
  set_table_properties(layout = "autofit") %>%
  set_caption("Supplementary Table 8: Joinpoint regression of trends coefficient of variation trends - LAU level") %>%
  fontsize(size = 8, part = 'all')

ft

save_as_docx(ft, path = 'tables/Supp_table8.docx', pr_section = prop_section(page_size(orient = 'landscape')))
