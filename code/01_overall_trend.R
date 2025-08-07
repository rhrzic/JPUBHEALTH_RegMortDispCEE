# Subnational mortality convergence in Europe by cause of death
# Rok Hrzic (r.hrzic (at) maastrichtuniversity.nl)
# Version April 2025

rm(list = ls())


require(tidyverse)
require(zoo)
library(flextable)
library(officer)


df = readRDS('temp/sdr.RDS')%>%
  filter(Year %in% 2000:2019) %>%
  mutate(Sex = ifelse(Sex == 'm', 'Men', 'Women')) %>%
  group_by(Country, RegionCode, Year, Sex, CauseCategory) %>%
  summarise(
    # Point estimate: median of the bootstrap replicates
    sdr_est = median(smoothed_sdr, na.rm = TRUE),
    # 95% percentile-based CI
    sdr_low = quantile(smoothed_sdr, probs = 0.025, na.rm = TRUE),
    sdr_up = quantile(smoothed_sdr, probs = 0.975, na.rm = TRUE),
    se = sd(smoothed_sdr, na.rm = TRUE)) %>%
  ungroup()

source("code/segmented.R")


df_segmented <- df %>%
  mutate(logsdr = log(sdr_est),
         logsdr_est = se/sdr_est) %>%
  group_by(Country, Sex, CauseCategory) %>%
  group_modify(~ fit_segmented_boot(.x, response = "logsdr", se = 'logsdr_est')) %>%
  ungroup()

saveRDS(df_segmented, 'temp/fig1segmented1.RDS')

# Create a summary data frame for breakpoint confidence intervals per subgroup.
# We extract one (representative) row per group and also determine the maximum y value
# for placing the CI illustration.
df_bp <- df_segmented %>%
  group_by(Country, Sex, CauseCategory) %>%
  slice(1) %>%
  ungroup()

# Plot the results using ggplot2
p2 <- ggplot(df_segmented, aes(x = Year, color = Country, group = interaction(RegionCode, Sex))) +
  geom_line(aes(y = sdr_est, linetype = Sex), alpha = 0.5, linewidth = 0.75) +
  geom_line(aes(y = sdr_low), alpha = 0.5, linewidth = 0.6, linetype = 'dotted') +
  geom_line(aes(y = sdr_up), alpha = 0.5, linewidth = 0.6, linetype = 'dotted') +
  geom_line(aes(y = exp(fitted), linetype = Sex), color = 'black', linewidth = 0.7) +
  facet_grid(CauseCategory ~ Country, scales="free_y") +
  ylab('Age-standardised death rate (per 1000)')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = 'top',
        text=element_text(size=14))


ggsave('plots/fig2.png', p2, height = 297, width = 210, units = 'mm')
ggsave('plots/fig2.pdf', p2, height = 297, width = 210, units = 'mm')


supp_table3 <- df_bp %>%
  dplyr::select(Country, Sex, CauseCategory,
    bp1_est, bp1_lower, bp1_upper,
    bp2_est, bp2_lower, bp2_upper,
    # bp3_est, bp3_lower, bp3_upper,
    AAPC1, AAPC1_lower, AAPC1_upper,
    AAPC2, AAPC2_lower, AAPC2_upper,
    AAPC3, AAPC3_lower, AAPC3_upper,
    # APC4, APC4_lower, APC4_upper,
    total_pct_change, total_pct_change_lower, total_pct_change_upper,
    AAPC_total, AAPC_total_lower, AAPC_total_upper
  ) %>%
  mutate(
    sigA1 = AAPC1_lower * AAPC1_upper > 0,
    sigA2 = AAPC2_lower * AAPC2_upper > 0,
    sigA3 = AAPC3_lower * AAPC3_upper > 0,
    sigT  = total_pct_change_lower * total_pct_change_upper > 0,
    sigAT  = AAPC_total_lower * AAPC_total_upper > 0,
    `Breakpoint 1 (95% CI)` = if_else(is.na(bp1_est), "-", sprintf("%.0f (%.1f, %.1f)", bp1_est, bp1_lower, bp1_upper)),
    `Breakpoint 2 (95% CI)` = if_else(is.na(bp2_est), "-", sprintf("%.0f (%.1f, %.1f)", bp2_est, bp2_lower, bp2_upper)),
    `AAPC 1 (%) (95% CI)` = if_else(is.na(AAPC1), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC1, if_else(sigA1, "*", ""), AAPC1_lower, AAPC1_upper)),
    `AAPC 2 (%) (95% CI)` = if_else(is.na(AAPC2), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC2, if_else(sigA2, "*", ""), AAPC2_lower, AAPC2_upper)),
    `AAPC 3 (%) (95% CI)` = if_else(is.na(AAPC3), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC3,if_else(sigA3, "*", ""), AAPC3_lower, AAPC3_upper)),
    `AAPC Total (%) (95% CI)` = if_else(is.na(AAPC_total), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC_total,if_else(sigAT, "*", ""), AAPC_total_lower, AAPC_total_upper)),
    `Cumulative change (%) (95% CI)` = if_else(is.na(total_pct_change), "-", sprintf("%.1f%s (%.1f, %.1f)", total_pct_change, if_else(sigT, "*", ""), total_pct_change_lower, total_pct_change_upper)),
    ) %>%
  dplyr::select(Country, Sex, CauseCategory,
         `Breakpoint 1 (95% CI)`,
         `Breakpoint 2 (95% CI)`,
         `AAPC 1 (%) (95% CI)`,
         `AAPC 2 (%) (95% CI)`,
         `AAPC 3 (%) (95% CI)`,
         `AAPC Total (%) (95% CI)`,
         `Cumulative change (%) (95% CI)`)

# Create the flextable
ft <- flextable(supp_table3) %>%
  set_table_properties(layout = "autofit") %>%
  set_caption("Supplementary Table 3: Joinpoint regression of regional mortality trends") %>%
  fontsize(size = 8, part = 'all')
  
ft

# Export the table to a Word document using officer
save_as_docx(ft, path = 'tables/Supp_table3.docx', pr_section = prop_section(page_size(orient = 'landscape')))



## LAU


df_LAU = readRDS('temp/sdr_LAU.RDS') %>%
  filter(Year %in% 2000:2019) %>%
  mutate(Sex = ifelse(Sex == 'm', 'Men', 'Women')) %>%
  group_by(Country, RegionCode, Year, Sex, CauseCategory) %>%
  summarise(
    # Point estimate: median of the bootstrap replicates
    sdr_est = median(smoothed_sdr, na.rm = TRUE),
    # 95% percentile-based CI
    sdr_low = quantile(smoothed_sdr, probs = 0.025, na.rm = TRUE),
    sdr_up = quantile(smoothed_sdr, probs = 0.975, na.rm = TRUE),
    se = sd(smoothed_sdr, na.rm = TRUE)) %>%
  ungroup()

source('code/segmented.R')

  
df_segmented_LAU <- df_LAU %>%
  mutate(logsdr = log(sdr_est),
         logsdr_est = se/sdr_est) %>%
  group_by(Country, Sex, CauseCategory) %>%
  group_modify(~ fit_segmented_boot(.x, response = "logsdr", se = 'logsdr_est')) %>%
  ungroup()

saveRDS(df_segmented_LAU, 'temp/fig1segmented_LAU.RDS')

df_segmented_LAU = readRDS('temp/fig1segmented_LAU.RDS')

# Create a summary data frame for breakpoint confidence intervals per subgroup.
# We extract one (representative) row per group and also determine the maximum y value
# for placing the CI illustration.
df_bp_LAU <- df_segmented_LAU %>%
  group_by(Country, Sex, CauseCategory) %>%
  slice(1) %>%
  ungroup()

supp_p1 = ggplot(df_segmented_LAU, aes(x = Year, y = sdr, shape = Sex, group = interaction(RegionCode, Sex))) +
  geom_line(aes(y = sdr_est, colour = Country, linetype = Sex), alpha = 0.5, linewidth = 0.75) +
  geom_line(aes(y = sdr_low, colour = Country), alpha = 0.5, linewidth = 0.6, linetype = 'dotted') +
  geom_line(aes(y = sdr_up,colour = Country), alpha = 0.5, linewidth = 0.6, linetype = 'dotted') +
  geom_line(aes(y = exp(fitted), linetype = Sex), color = 'black', linewidth = 0.7) +
  facet_wrap(CauseCategory ~ Country, ncol = 3, scales="free_y") +
  ylab('Age-standardised death rate (per 1000)')+
  scale_color_manual(values = c('EST' = '#B79F00',
                                'LTU' = '#00BA38',
                                'SVK' = '#F564E3'))+
  scale_linetype(guide = guide_legend(override.aes = list(color = 'black')))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = 'top',
        text=element_text(size=14))

ggsave('plots/supp_fig1.png', supp_p1, height = 297, width = 200, units = 'mm')

supp_table4 <- df_bp_LAU %>%
  dplyr::select(Country, Sex, CauseCategory,
                bp1_est, bp1_lower, bp1_upper,
                bp2_est, bp2_lower, bp2_upper,
                bp3_est, bp3_lower, bp3_upper,
                AAPC1, AAPC1_lower, AAPC1_upper,
                AAPC2, AAPC2_lower, AAPC2_upper,
                AAPC3, AAPC3_lower, AAPC3_upper,
                AAPC4, AAPC4_lower, AAPC4_upper,
                total_pct_change, total_pct_change_lower, total_pct_change_upper,
                AAPC_total, AAPC_total_lower, AAPC_total_upper
  ) %>%
  mutate(
    sigA1 = AAPC1_lower * AAPC1_upper > 0,
    sigA2 = AAPC2_lower * AAPC2_upper > 0,
    sigA3 = AAPC3_lower * AAPC3_upper > 0,
    sigT  = total_pct_change_lower * total_pct_change_upper > 0,
    sigAT  = AAPC_total_lower * AAPC_total_upper > 0,
    `Breakpoint 1 (95% CI)` = if_else(is.na(bp1_est), "-", sprintf("%.0f (%.1f, %.1f)", bp1_est, bp1_lower, bp1_upper)),
    `Breakpoint 2 (95% CI)` = if_else(is.na(bp2_est), "-", sprintf("%.0f (%.1f, %.1f)", bp2_est, bp2_lower, bp2_upper)),
    `Breakpoint 3 (95% CI)` = if_else(is.na(bp3_est), "-", sprintf("%.0f (%.1f, %.1f)", bp3_est, bp3_lower, bp3_upper)),
    `AAPC 1 (%) (95% CI)` = if_else(is.na(AAPC1), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC1, if_else(sigA1, "*", ""), AAPC1_lower, AAPC1_upper)),
    `AAPC 2 (%) (95% CI)` = if_else(is.na(AAPC2), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC2, if_else(sigA2, "*", ""), AAPC2_lower, AAPC2_upper)),
    `AAPC 3 (%) (95% CI)` = if_else(is.na(AAPC3), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC3,if_else(sigA3, "*", ""), AAPC3_lower, AAPC3_upper)),
    `AAPC 4 (%) (95% CI)` = if_else(is.na(AAPC4), "-", sprintf("%.1f%s (%.1f, %.1f)", AAPC4,if_else(sigA3, "*", ""), AAPC4_lower, AAPC4_upper)),
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
ft <- flextable(supp_table4) %>%
  set_table_properties(layout = "autofit") %>%
  set_caption("Supplementary Table 4: Joinpoint regression of regional mortality trends (LAU level)") %>%
  fontsize(size = 8, part = 'all')

ft

# Export the table to a Word document using officer
save_as_docx(ft, path = 'tables/Supp_table4.docx', pr_section = prop_section(page_size(orient = 'landscape')))


