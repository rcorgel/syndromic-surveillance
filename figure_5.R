################################################################################
# File Name: figure_5                                                          #
#                                                                              #
# Purpose:   Create figure 5 for the manuscript.                               #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Create figure 5                                                #
#                                                                              #
# Project:   Syndromic Surveillance                                            #
# Author:    Ronan Corgel                                                      #
################################################################################

####################
# 1. SET-UP SCRIPT #
####################

# Start with a clear environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(lubridate)
library(assertr)
library(ISOweek)
library(splus2R)
library(scales)
library(sf)
library(geosphere)
library(cowplot)
library(zoo)
library(gsignal)
library(mgcv)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

######################
# 2. CREATE FIGURE 5 #
######################

# Load CDC data
cdc_est <- read_csv('raw/flu-disease-burden-past-season-estimates (1).csv')

cdc_overall <- cdc_est |>
  dplyr::filter(`Age Group` == 'All Ages') |>
  dplyr::select(c(`Flu Season`, Population, `Symptomatic Illnesses`, `Symptomatic Illnesses 95% LL`, `Symptomatic Illnesses 95% UL`)) |>
  mutate(perc_50 = `Symptomatic Illnesses` / Population,
         perc_05 = `Symptomatic Illnesses 95% LL` / Population,
         perc_95 = `Symptomatic Illnesses 95% UL` / Population) |>
  dplyr::select(c(`Flu Season`, perc_50, perc_05, perc_95)) |>
  mutate(Source = 'CDC')


# Set theme
forrest_theme <- theme(legend.position =  "none",
                       legend.position.inside = c(0.873, 0.237),
                       axis.text = element_text(size=14, color = 'black'),
                       axis.title = element_text(size=16),
                       axis.title.y = element_blank(),
                       strip.background = element_blank(),
                       legend.title = element_text(size=16),
                       legend.text = element_text(size=14),
                       plot.title = element_text(size=18, hjust = 0),
                       legend.box="vertical",
                       legend.box.background = element_rect(colour = "black", 
                                                            fill = 'white'))


# Plot
cdc_season_fig <- ggplot(data=cdc_overall, aes(x=`Flu Season`, y=perc_50, ymin=perc_05, ymax=perc_95)) +
  geom_pointrange(aes(group = `Flu Season`, color = `Flu Season`, fill = `Flu Season`), size = 1) +
                  #position = position_dodge2(width = 0.5), size = 1, alpha = 0.7) + 
  geom_hline(yintercept=0.025, lty=2) +  # add a dotted line at x=1 after flip
  geom_hline(yintercept=0.10, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Season") + ylab("Percent of Population") + scale_y_continuous(breaks = seq(0, 0.25, by = 0.05), limits = c(0, 0.25)) +
  ggtitle("Influenza Burden by Season") +
  theme_minimal() +  # use a white background
  forrest_theme

cdc_season_fig 

# Overall predictions
flu_pred_dat_2016 <- readRDS('./tmp/flu_2016_predictions')

flu_county_2016 <- flu_pred_dat_2016 |>
  group_by(county_fips, year_week_dt) |>
  mutate(prediction = sum(pred_case_count)) |>
  distinct(county_fips, year_week_dt, prediction)

# All cause
all_cause <- read.csv('./raw/county_week_ac_v3_imputed.csv', header = T)
all_cause$year_week <- sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", all_cause$year_week)
all_cause$year_week_dt <- ISOweek2date(all_cause$year_week)
all_cause$county_fips <- ifelse(nchar(all_cause$county_fips) == 4, str_pad(all_cause$county_fips, 5, pad = "0"), all_cause$county_fips)


flu_county_2016 <- left_join(flu_county_2016, all_cause[, c(1, 5, 6)],
                             by = c('county_fips', 'year_week_dt'))

flu_county_2016_perc <- flu_county_2016 |>
  group_by(county_fips) |>
  dplyr::filter(year_week_dt > as.Date('2016-10-01') & year_week_dt < as.Date('2017-06-01')) |>
  mutate(prediction = sum(prediction),
         all_cause_sum = sum(all_cause_wtd)) |>
  distinct(county_fips, prediction, all_cause_sum) |>
  mutate(prediction_50 = median(prediction * runif(10000, min = 1/0.58, max = 1/0.42)),
         prediction_05 = quantile(prediction * runif(10000, min = 1/0.58, max = 1/0.42), probs = 0.05),
         prediction_95 = quantile(prediction * runif(10000, min = 1/0.58, max = 1/0.42), probs = 0.95)) |>
  mutate(perc_50 = prediction_50 / all_cause_sum,
         perc_05 = prediction_05 / all_cause_sum,
         perc_95 = prediction_95 / all_cause_sum)

flu_2016_perc <- flu_county_2016_perc |> 
  mutate(Season = '2016-2017') |>
  group_by(Season) |>
  mutate(perc_50 = sum(prediction_50) / sum(all_cause_sum),
         perc_05 = sum(prediction_05) / sum(all_cause_sum),
         perc_95 = sum(prediction_95) / sum(all_cause_sum)) |>
  distinct(Season, perc_50, perc_05, perc_95)





flu_2017_perc <- flu_county_2017 |>
  group_by(year_week_dt, county_fips) |>
  mutate(prediction = sum(prediction),
         all_cause_sum = sum(all_cause_wtd)) |>
  distinct(year_week_dt, prediction, all_cause_sum) |>
  mutate(prediction_50 = median(prediction * runif(10000, min = 1/0.58, max = 1/0.42)),
         prediction_05 = quantile(prediction * runif(10000, min = 1/0.58, max = 1/0.42), probs = 0.05),
         prediction_95 = quantile(prediction * runif(10000, min = 1/0.58, max = 1/0.42), probs = 0.95)) |>
  ungroup() 

flu_2017_perc$year <- '2017'
flu_2016_perc$year <- '2016'

cumsum_prev <- rbind(flu_2017_perc, flu_2016_perc)
ggplot(cumsum_prev) + geom_line(aes(x = year_week_dt, y = perc_50_cumsum, color = year))

flu_2016_perc <- flu_county_2016_perc |> 
  mutate(Season = '2016-2017') |>
  group_by(Season) |>
  mutate(perc_50 = sum(prediction_50) / sum(all_cause_sum),
         perc_05 = sum(prediction_05) / sum(all_cause_sum),
         perc_95 = sum(prediction_95) / sum(all_cause_sum)) |>
  distinct(Season, perc_50, perc_05, perc_95)


flu_pred_dat_2017_all <- flu_pred_dat_2017 |>
  group_by(county_fips, year_week_dt) |>
  mutate(flu_sum = sum(flu_count),
         pred_sum = sum(pred_case_count),
         pred_2_sum = sum(pred_interaction_count),
         pred_3_sum = sum(pred_fever_count)) |>
  distinct(county_fips, year_week_dt, flu_sum, pred_sum, pred_2_sum, pred_3_sum)

ggplot(flu_pred_dat_2017_all[flu_pred_dat_2017_all$county_fips == '36061',]) + 
  geom_line(aes(x = year_week_dt, y = pred_sum), color = 'green') +
  geom_line(aes(x = year_week_dt, y = flu_sum), color = 'blue') + 
  geom_line(aes(x = year_week_dt, y = pred_2_sum), color = 'red') +
  geom_line(aes(x = year_week_dt, y = pred_3_sum), color = 'orange') 




# 2017-18
flu_pred_dat_2017 <- readRDS('./tmp/flu_2017_predictions')

flu_county_2017 <- flu_pred_dat_2017 |>
  group_by(county_fips, year_week_dt) |>
  mutate(prediction = sum(pred_case_count)) |>
  distinct(county_fips, year_week_dt, prediction)

flu_county_2017 <- left_join(flu_county_2017, all_cause[, c(1, 5, 6)],
                             by = c('county_fips', 'year_week_dt'))

flu_county_2017_perc <- flu_county_2017 |>
  group_by(county_fips) |>
  dplyr::filter(year_week_dt > as.Date('2017-10-01') & year_week_dt < as.Date('2018-06-01')) |>
  mutate(prediction = sum(prediction),
         all_cause_sum = sum(all_cause_wtd)) |>
  distinct(county_fips, prediction, all_cause_sum) |>
  mutate(prediction_50 = median(prediction * runif(10000, min = 1/0.58, max = 1/0.42)),
         prediction_05 = quantile(prediction * runif(10000, min = 1/0.58, max = 1/0.42), probs = 0.05),
         prediction_95 = quantile(prediction * runif(10000, min = 1/0.58, max = 1/0.42), probs = 0.95)) |>
  mutate(perc_50 = prediction_50 / all_cause_sum,
         perc_05 = prediction_05 / all_cause_sum,
         perc_95 = prediction_95 / all_cause_sum)

flu_2017_perc <- flu_county_2017_perc |> 
  mutate(Season = '2017-2018') |>
  group_by(Season) |>
  mutate(perc_50 = sum(prediction_50) / sum(all_cause_sum),
         perc_05 = sum(prediction_05) / sum(all_cause_sum),
         perc_95 = sum(prediction_95) / sum(all_cause_sum)) |>
  distinct(Season, perc_50, perc_05, perc_95)


# 2018
flu_pred_dat_2018 <- readRDS('./tmp/flu_3_proc_pred')

flu_pred_dat_2018$pred_count <- ifelse(flu_pred_dat_2018$pred >  0.6593 &
                                         flu_pred_dat_2018$age_grp == '0-4', 1, NA) 
flu_pred_dat_2018$pred_count <- ifelse(flu_pred_dat_2018$pred > 0.7715 &
                                         flu_pred_dat_2018$age_grp == '5-12', 1, flu_pred_dat_2018$pred_count) 
flu_pred_dat_2018$pred_count <- ifelse(flu_pred_dat_2018$pred > 0.5257 &
                                         flu_pred_dat_2018$age_grp == '13-17', 1, flu_pred_dat_2018$pred_count) 
flu_pred_dat_2018$pred_count <- ifelse(flu_pred_dat_2018$pred > 0.4273 &
                                         flu_pred_dat_2018$age_grp == '18-49', 1, flu_pred_dat_2018$pred_count) 
flu_pred_dat_2018$pred_count <- ifelse(flu_pred_dat_2018$pred > 0.332 &
                                         flu_pred_dat_2018$age_grp == '50-64', 1, flu_pred_dat_2018$pred_count) 
flu_pred_dat_2018$pred_count <- ifelse(flu_pred_dat_2018$pred > 0.278 &
                                         flu_pred_dat_2018$age_grp == '>=65', 1, flu_pred_dat_2018$pred_count) 
flu_pred_dat_2018$pred_count <- ifelse(is.na(flu_pred_dat_2018$pred_count), 0, flu_pred_dat_2018$pred_count)

flu_pred_dat_2018$pred_count <- flu_pred_dat_2018$pred_count * flu_pred_dat_2018$patient_count

flu_county_2018 <- flu_pred_dat_2018 |>
  group_by(county_fips, year_week_dt) |>
  mutate(prediction = sum(pred_count)) |>
  distinct(county_fips, year_week_dt, prediction)

flu_county_2018 <- left_join(flu_county_2018, all_cause[, c(1, 5, 6)],
                             by = c('county_fips', 'year_week_dt'))

flu_county_2018_perc <- flu_county_2018 |>
  group_by(county_fips) |>
  dplyr::filter(year_week_dt > as.Date('2018-10-01') & year_week_dt < as.Date('2019-06-01')) |>
  mutate(prediction = sum(prediction),
         all_cause_sum = sum(all_cause_wtd)) |>
  distinct(county_fips, prediction, all_cause_sum) |>
  mutate(prediction_50 = median(prediction * runif(10000, min = 1/0.58, max = 1/0.42)),
         prediction_05 = quantile(prediction * runif(10000, min = 1/0.58, max = 1/0.42), probs = 0.05),
         prediction_95 = quantile(prediction * runif(10000, min = 1/0.58, max = 1/0.42), probs = 0.95)) |>
  mutate(perc_50 = prediction_50 / all_cause_sum,
         perc_05 = prediction_05 / all_cause_sum,
         perc_95 = prediction_95 / all_cause_sum)

flu_2018_perc <- flu_county_2018_perc |> 
  mutate(Season = '2018-2019') |>
  group_by(Season) |>
  mutate(perc_50 = sum(prediction_50) / sum(all_cause_sum),
         perc_05 = sum(prediction_05) / sum(all_cause_sum),
         perc_95 = sum(prediction_95) / sum(all_cause_sum)) |>
  distinct(Season, perc_50, perc_05, perc_95)


estimates <- rbind(flu_2016_perc, flu_2017_perc, flu_2018_perc)
estimates <- estimates |> rename('Flu Season' = 'Season') |>
  mutate(Source = 'Syndromic')

all_estimates <- rbind(cdc_overall, estimates)
all_estimates <- all_estimates |>
  dplyr::filter(`Flu Season` == '2016-2017' | `Flu Season` == '2017-2018' | `Flu Season` == '2018-2019')

all_estimates$Name <- c('CDC Estimate', 'CDC Estimate', 'CDC Estimate', ' ', '  ', 'Syndromic Estimate')

all_estimates$Name <- factor(all_estimates$Name, levels = c('CDC Estimate', ' ', '  ', 'Syndromic Estimate'))

colors <- c('CDC Estimate' = 'darkgray', ' ' = '#3cbb75ff', '  ' = '#9B59B6', 'Syndromic Estimate' = '#d7642c')


forrest_theme <- theme(legend.position =  "bottom",
                       legend.position.inside = c(0.873, 0.237),
                       axis.text = element_text(size=14, color = 'black'),
                       axis.title = element_text(size=16),
                       strip.background = element_blank(),
                       legend.title = element_text(size=16),
                       legend.text = element_text(size=14),
                       plot.title = element_text(size=18, hjust = 0),
                       legend.box="vertical",
                       legend.box.background = element_rect(colour = "white", 
                                                            fill = 'white'))


season_comp_plot <- ggplot(data=all_estimates, aes(x=`Flu Season`, y=perc_50, ymin=perc_05, ymax=perc_95)) +
  
  #geom_hline(yintercept=0.025, lty=2, linewidth = 1.5) +  # add a dotted line at x=1 after flip
  #geom_hline(yintercept=0.10, lty=2, linewidth = 1.5) +  # add a dotted line at x=1 after flip
  geom_pointrange(aes(group = Name, color = Name, fill = Name),
                  position = position_dodge2(width = 0.5), size = 1.5, alpha = 0.75, linewidth = 1.2) + 
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Influenza Season") + ylab("Proportion of Population") + #scale_y_continuous(breaks = seq(0, 0.15, by = 0.05), limits = c(0, 0.17)) +
  ggtitle("Influenza Prevalence by Season") +
  theme_minimal() +  # use a white background
  forrest_theme + ylim(0.065, 0.165) +
  scale_fill_manual('', values = colors, labels = c(" " = "", "  " = "", 'Syndromic Estimat' = "Syndromic Estimat", "CDC Estimate" = 'CDC Estimate')) +
  scale_color_manual('', values = colors, labels = c(" " = "", "  " = "", 'Syndromic Estimat' = "Syndromic Estimat", "CDC Estimate" = 'CDC Estimate')) +
  guides(fill=guide_legend(nrow=1,byrow=FALSE))

season_comp_plot

ggsave('./figs/season_comp_plot.jpg', plot = season_comp_plot, height = 4, width = 7)

# Load map data
usa_albers_state <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_state.rds'))   # convert to sf
usa_albers_county <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_county.rds'))   # convert to sf

# Filter 
flu_county_2018_perc_filt <- flu_county_2018_perc |> dplyr::filter(all_cause_sum > 20000)

# Merge on prev data
usa_albers_county_2018 <- left_join(usa_albers_county, flu_county_2018_perc_filt, by = c('GEOID' = 'county_fips'))

# Set the map theme
map_theme <- theme(legend.position = "right",
                   axis.title = element_text(size=16),
                   strip.background = element_blank(),
                   plot.title = element_blank(),
                   panel.border = element_rect(fill=NA, linewidth = 0.8, color = 'white'),
                   legend.title = element_text(size = 16),
                   legend.text = element_text(size = 12),
                   legend.key.height = unit(0.7, 'cm'),
                   legend.key.width = unit(0.55, "cm"))

map_2018 <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_2018), aes(fill = perc_50, group = GEOID), color= 'black', linewidth = 0.01) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.10) +
  scale_fill_gradient('2018-19\nSyndromic\nFlu Prev.\n', low = "white", high = '#9B59B6', na.value = "grey90",
                      breaks = c(0.04, 0.20, 0.40), 
                      limits = c(0, 0.41),
                      labels = c("0","0.20","0.40")) + 
  theme_void() + map_theme 

map_2018 

ggsave('./figs/map_2018.jpg', plot = map_2018, height = 4, width = 7)


# Filter 
flu_county_2017_perc_filt <- flu_county_2017_perc |> dplyr::filter(all_cause_sum > 20000)

# Merge on prev data
usa_albers_county_2017 <- left_join(usa_albers_county, flu_county_2017_perc_filt, by = c('GEOID' = 'county_fips'))

map_2017 <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_2017), aes(fill = perc_50, group = GEOID), color= 'black', linewidth = 0.01) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.10) +
  scale_fill_gradient('2017-18\nSyndromic\nFlu Burden \n', low = "white", high = '#d7642c', na.value = "grey90",
                      breaks = c(0, 0.20, 0.40), 
                      limits = c(0, 0.41),
                      labels = c("0%","20%","40%")) + 
  theme_void() + map_theme 

map_2017 


# Filter 
flu_county_2016_perc_filt <- flu_county_2016_perc |> dplyr::filter(all_cause_sum > 20000)

# Merge on prev data
usa_albers_county_2016 <- left_join(usa_albers_county, flu_county_2016_perc_filt, by = c('GEOID' = 'county_fips'))

map_2016 <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_2016), aes(fill = perc_50, group = GEOID), color= 'black', linewidth = 0.01) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.10) +
  scale_fill_gradient('2016-17\nSyndromic\nFlu Burden \n', low = "white", high = '#3cbb75ff', na.value = "grey90",
                      breaks = c(0, 0.20, 0.40), 
                      limits = c(0, 0.41),
                      labels = c("0%","20%","40%")) + 
  theme_void() + map_theme 


den_2016 <- ggplot(usa_albers_county_2016) + geom_density(aes(perc_50), color = '#3cbb75ff', fill = '#3cbb75ff', alpha = 0.5,
                                              linewidth = 1.5) + 
  geom_vline(xintercept=all_estimates[4, ]$perc_50, lty=2, linewidth = 1.3, color = '#3cbb75ff') +
  theme_minimal() + theme(axis.title = element_text(size=16),
                          strip.background = element_blank(),
                          axis.text.y = element_blank(),
                          axis.text.x = element_text(size=14),
                          panel.grid.major.y = element_blank(),
                          panel.grid.minor = element_blank()) + xlab('Percent of Population') + ylab('Density') + ggtitle('') +
  scale_x_continuous(labels = label_percent(), limits = c(0, 0.41))

den_2017 <- ggplot(usa_albers_county_2017) + geom_density(aes(perc_50), color = '#d7642c', fill = '#d7642c', alpha = 0.5,
                                                          linewidth = 1.5) + 
  geom_vline(xintercept=all_estimates[5, ]$perc_50, lty=2, linewidth = 1.3, color = '#d7642c') +
  theme_minimal() + theme(axis.title = element_text(size=16),
                          strip.background = element_blank(),
                          axis.text.y = element_blank(),
                          axis.text.x = element_text(size=14),
                          panel.grid.major.y = element_blank(),
                          panel.grid.minor = element_blank()) + xlab('Percent of Population') + ylab('Density') + ggtitle('') +
  scale_x_continuous(labels = label_percent(), limits = c(0, 0.41))


den_2018 <- ggplot(usa_albers_county_2018) + geom_density(aes(perc_50), color = '#9B59B6', fill = '#9B59B6', alpha = 0.5,
                                                          linewidth = 1.5) + 
  geom_vline(xintercept=all_estimates[6, ]$perc_50, lty=2, linewidth = 1.3, color = '#9B59B6') +
  theme_minimal() + theme(axis.title = element_text(size=16),
                          strip.background = element_blank(),
                          axis.text.y = element_blank(),
                          axis.text.x = element_text(size=14),
                          panel.grid.major.y = element_blank(),
                          panel.grid.minor = element_blank()) + xlab('Percent of Population') + ylab('Density') + ggtitle('') +
  scale_x_continuous(labels = label_percent(), limits = c(0, 0.41))


map_2016 





figure_5_map <- cowplot::plot_grid(map_2018, den_2018, map_2017,den_2017,  map_2016, den_2016,
                               nrow = 3,
                               labels = c("", '', "", '', "", ''),
                               label_size = 20)

figure_5 <- cowplot::plot_grid(season_comp_plot, figure_5_map,
                                   rel_widths = c(1, 1),
                                   nrow = 1,
                                   labels = c("a", "b"),
                                   label_size = 20)

ggsave('./figs/figure_5.jpg', plot = figure_5, height = 8, width = 20)


measure <- readRDS('./figs/measure_plot.RData')
peak <- readRDS('./figs/peak_plot_2018.RData')


figure_test_map <- cowplot::plot_grid(map_2018, den_2018,
                                   nrow = 2,
                                   labels = c("", ''),
                                   label_size = 20)



flu_county_perc_filt <- left_join(flu_county_2018_perc_filt, flu_county_2017_perc_filt, by = 'county_fips')

cor(flu_county_perc_filt$perc_50.x, flu_county_perc_filt$perc_50.y, use = 'complete.obs', method = 'spearman')

cdc_overall_age <- cdc_est |>
  dplyr::filter(`Flu Season` == '2016-2017') |>
  dplyr::filter(`Age Group` != 'All Ages') |>
  dplyr::select(c(`Age Group`, Population, `Symptomatic Illnesses`, `Symptomatic Illnesses 95% LL`, `Symptomatic Illnesses 95% UL`)) |>
  mutate(perc_50 = `Symptomatic Illnesses` / Population,
         perc_05 = `Symptomatic Illnesses 95% LL` / Population,
         perc_95 = `Symptomatic Illnesses 95% UL` / Population) |>
  dplyr::select(c(`Age Group`, perc_50, perc_05, perc_95)) |>
  mutate(Source = 'CDC')

flu_county_2018_age <- flu_pred_dat_2018 |>
  group_by(county_fips, year_week_dt, age_grp) |>
  mutate(prediction = sum(pred_count)) |>
  distinct(county_fips, age_grp, year_week_dt, prediction)


# All cause
all_cause_age <- read.csv('./raw/county_week_age_ac_v3_imputed.csv', header = T)
all_cause_age$year_week <- sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", all_cause_age$year_week)
all_cause_age$year_week_dt <- ISOweek2date(all_cause_age$year_week)
all_cause_age$county_fips <- ifelse(nchar(all_cause_age$county_fips) == 4, str_pad(all_cause_age$county_fips, 5, pad = "0"), all_cause_age$county_fips)


flu_county_2018_age <- left_join(flu_county_2018_age, all_cause_age[, c(1, 4, 6, 7)],
                             by = c('county_fips', 'year_week_dt', 'age_grp'))

flu_county_2018_age$age_grp <- ifelse(flu_county_2018_age$age_grp == '5-12', '5-17', flu_county_2018_age$age_grp)
flu_county_2018_age$age_grp <- ifelse(flu_county_2018_age$age_grp == '13-17', '5-17', flu_county_2018_age$age_grp)

flu_county_2018_perc <- flu_county_2018_age |>
  dplyr::filter(year_week_dt > as.Date('2018-10-01') & year_week_dt < as.Date('2019-06-01')) |>
  group_by(county_fips, age_grp) |>
  mutate(prediction = sum(prediction),
         all_cause_sum = sum(all_cause_wtd)) |>
  distinct(county_fips, age_grp, prediction, all_cause_sum) |>
  mutate(prediction_50 = median(prediction * runif(10000, min = 1/0.58, max = 1/0.42)),
         prediction_05 = quantile(prediction * runif(10000, min = 1/0.58, max = 1/0.42), probs = 0.05),
         prediction_95 = quantile(prediction * runif(10000, min = 1/0.58, max = 1/0.42), probs = 0.95)) |>
  mutate(perc_50 = prediction_50 / all_cause_sum,
         perc_05 = prediction_05 / all_cause_sum,
         perc_95 = prediction_95 / all_cause_sum)

flu_2018_perc <- flu_county_2018_perc |> 
  group_by(age_grp) |>
  mutate(perc_50 = sum(prediction_50) / sum(all_cause_sum),
         perc_05 = sum(prediction_05) / sum(all_cause_sum),
         perc_95 = sum(prediction_95) / sum(all_cause_sum)) |>
  distinct(age_grp, perc_50, perc_05, perc_95) 




flu_county_2018_age_grp <- flu_county_2018_age |>
  group_by(county_fips, age_grp) |>
  mutate(prediction = sum(prediction),
         all_cause_sum = sum(all_cause_wtd)) |>
  distinct(county_fips, age_grp, year_week_dt, prediction, all_cause_sum)



flu_2018_perc$`Age Group` <- c('0-4', '4-17', '18-49', '50-64', '65+')
flu_2018_perc$Source <- 'Prediction'
cdc_overall_age$`Age Group` <- c('0-4', '4-17', '18-49', '50-64', '65+')

age_group <- rbind(flu_2018_perc, cdc_overall_age)



age_group$Name <- c('Syndromic Estimate 1', 'Syndromic Estimate 2', 'Syndromic Estimate 3', 'Syndromic Estimate 4', 'Syndromic Estimate 5', 'CDC Estimate', 'CDC Estimate', 'CDC Estimate', 'CDC Estimate', 'CDC Estimate')

all_estimates$Name <- factor(all_estimates$Name, levels = c('Syndromic Estimate 1', 'Syndromic Estimate 2', 'Syndromic Estimate 3', 'Syndromic Estimate 4', 'Syndromic Estimate 5', 'CDC Estimate'))

age_group$`Age Group` <- factor(age_group$`Age Group`, levels = c('0-4', '4-17', '18-49', '50-64', '65+'))

colors <- c('CDC Estimate' = 'darkgray', 'Syndromic Estimate 1' = '#b9d4ec', 'Syndromic Estimate 2' = '#8bb7e0', 'Syndromic Estimate 3' = '#347DC1', 'Syndromic Estimate 4' = '#285f93', 'Syndromic Estimate 5' = '#1b4164')

forrest_theme <- theme(legend.position =  "bottom",
                       legend.position.inside = c(0.873, 0.237),
                       axis.text = element_text(size=14, color = 'black'),
                       axis.title = element_text(size=16),
                       strip.background = element_blank(),
                       legend.title = element_text(size=16),
                       legend.text = element_text(size=14),
                       plot.title = element_text(size=18, hjust = 0),
                       legend.box="vertical",
                       legend.box.background = element_rect(colour = "white", 
                                                            fill = 'white'))

age_comp_plot <- ggplot(data=age_group, aes(x=`Age Group`, y=perc_50, ymin=perc_05, ymax=perc_95)) +
  
  geom_pointrange(aes(group = Name, color = Name, fill = Name),
                  position = position_dodge2(width = 0.5), size = 1.7, alpha = 0.75, linewidth = 1.2) + 
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Age Group") + ylab("Percent of Population") + #scale_y_continuous(breaks = seq(0, 0.35, by = 0.05), limits = c(0, 0.35)) +
  ggtitle("Influenza Burden by Age Group") +
  theme_minimal() +  # use a white background
  forrest_theme +
  scale_fill_manual('', values = colors, labels = c("Syndromic Estimate 1" = "", "Syndromic Estimate 2" = "", "Syndromic Estimate 3" = "", "Syndromic Estimate 4" = "", 'Syndromic Estimate 5' = "Syndromic Estimate", "CDC Estimate" = 'CDC Estimate')) +
  scale_color_manual('', values = colors, labels = c("Syndromic Estimate 1" = "", "Syndromic Estimate 2" = "", "Syndromic Estimate 3" = "", "Syndromic Estimate 4" = "", 'Syndromic Estimate 5' = "Syndromic Estimate", "CDC Estimate" = 'CDC Estimate')) +
  guides(fill=guide_legend(nrow=1,byrow=FALSE))

age_comp_plot





# Filter 
flu_county_2018_perc_age_filt <- left_join(flu_county_2018_perc_filt[, c(1)], flu_county_2018_perc, by = 'county_fips')

# Merge on prev data
usa_albers_county_0_4 <- left_join(usa_albers_county,
                                   flu_county_2018_perc_age_filt[flu_county_2018_perc_age_filt$age_grp == '0-4',], 
                                   by = c('GEOID' = 'county_fips'))

# Set the map theme
map_theme <- theme(legend.position = "right",
                   axis.title = element_text(size=16),
                   strip.background = element_blank(),
                   plot.title = element_blank(),
                   panel.border = element_rect(fill=NA, linewidth = 0.8, color = 'white'),
                   legend.title = element_text(size = 16),
                   legend.text = element_text(size = 12),
                   legend.key.height = unit(0.7, 'cm'),
                   legend.key.width = unit(0.4, "cm"))

map_0_4 <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_0_4), aes(fill = perc_50, group = GEOID), color= 'black', linewidth = 0.01) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.10) +
  scale_fill_gradient('0-4\nSyndromic\nFlu Burden \n', low = "white", high = '#347DC1', na.value = "grey90",
                      breaks = c(0, 0.50, 1), 
                      limits = c(0, 1),
                      labels = c("0%","50%", "100%")) + 
  theme_void() + map_theme 

map_0_4 


# Merge on prev data
usa_albers_county_5_17 <- left_join(usa_albers_county,
                                   flu_county_2018_perc_age_filt[flu_county_2018_perc_age_filt$age_grp == '5-17',], 
                                   by = c('GEOID' = 'county_fips'))

map_5_17 <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_5_17), aes(fill = perc_50, group = GEOID), color= 'black', linewidth = 0.01) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.10) +
  scale_fill_gradient('5-17\nSyndromic\nFlu Burden \n', low = "white", high = '#347DC1', na.value = "grey90",
                      breaks = c(0, 0.50, 1), 
                      limits = c(0, 1),
                      labels = c("0%","50%", "100%")) + 
  theme_void() + map_theme 

map_5_17 

# Merge on prev data
usa_albers_county_18_49 <- left_join(usa_albers_county,
                                    flu_county_2018_perc_age_filt[flu_county_2018_perc_age_filt$age_grp == '18-49',], 
                                    by = c('GEOID' = 'county_fips'))

map_18_49 <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_18_49), aes(fill = perc_50, group = GEOID), color= 'black', linewidth = 0.01) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.10) +
  scale_fill_gradient('18-49\nSyndromic\nFlu Burden \n', low = "white", high = '#347DC1', na.value = "grey90",
                      breaks = c(0, 0.25, 0.50), 
                      limits = c(0, 0.50),
                      labels = c("0%","25%","50%")) + 
  theme_void() + map_theme 

map_18_49

# Merge on prev data
usa_albers_county_50_64 <- left_join(usa_albers_county,
                                     flu_county_2018_perc_age_filt[flu_county_2018_perc_age_filt$age_grp == '50-64',], 
                                     by = c('GEOID' = 'county_fips'))

map_50_64 <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_50_64), aes(fill = perc_50, group = GEOID), color= 'black', linewidth = 0.01) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.10) +
  scale_fill_gradient('50-64\nSyndromic\nFlu Burden \n', low = "white", high = '#347DC1', na.value = "grey90",
                      breaks = c(0, 0.20, 0.40), 
                      limits = c(0, 0.40),
                      labels = c("0%","20%","40%")) + 
  theme_void() + map_theme 

map_50_64

# Merge on prev data
usa_albers_county_65 <- left_join(usa_albers_county,
                                     flu_county_2018_perc_age_filt[flu_county_2018_perc_age_filt$age_grp == '>=65',], 
                                     by = c('GEOID' = 'county_fips'))

map_65 <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_65), aes(fill = perc_50, group = GEOID), color= 'black', linewidth = 0.01) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.10) +
  scale_fill_gradient('65+\nSyndromic\nFlu Burden \n', low = "white", high = '#347DC1', na.value = "grey90",
                      breaks = c(0, 0.20, 0.40), 
                      limits = c(0, 0.40),
                      labels = c("0%","20%","40%")) + 
  theme_void() + map_theme 

map_65



den_0_4 <- ggplot(usa_albers_county_0_4) + geom_density(aes(perc_50), color = '#b9d4ec', fill = '#b9d4ec', alpha = 0.5,
                                                          linewidth = 1.5) + 
  geom_vline(xintercept=age_group[1, ]$perc_50, lty=2, linewidth = 1.3, color = '#b9d4ec') +
  theme_minimal() + theme(axis.title = element_text(size=16),
                          strip.background = element_blank(),
                          axis.text.y = element_blank(),
                          axis.text.x = element_text(size=14),
                          panel.grid.major.y = element_blank(),
                          panel.grid.minor = element_blank()) + xlab('Percent of Population') + ylab('Density') + ggtitle('') +
  scale_x_continuous(labels = label_percent(), limits = c(0, 1))

den_0_4

den_5_17 <- ggplot(usa_albers_county_5_17) + geom_density(aes(perc_50), color = '#8bb7e0', fill = '#8bb7e0', alpha = 0.5,
                                                        linewidth = 1.5) + 
  geom_vline(xintercept=age_group[2, ]$perc_50, lty=2, linewidth = 1.3, color = '#8bb7e0') +
  theme_minimal() + theme(axis.title = element_text(size=16),
                          strip.background = element_blank(),
                          axis.text.y = element_blank(),
                          axis.text.x = element_text(size=14),
                          panel.grid.major.y = element_blank(),
                          panel.grid.minor = element_blank()) + xlab('Percent of Population') + ylab('Density') + ggtitle('') +
  scale_x_continuous(labels = label_percent(), limits = c(0, 1))

den_5_17

den_18_49 <- ggplot(usa_albers_county_18_49) + geom_density(aes(perc_50), color = '#347DC1', fill = '#347DC1', alpha = 0.5,
                                                          linewidth = 1.5) + 
  geom_vline(xintercept=age_group[3, ]$perc_50, lty=2, linewidth = 1.3, color = '#347DC1') +
  theme_minimal() + theme(axis.title = element_text(size=16),
                          strip.background = element_blank(),
                          axis.text.y = element_blank(),
                          axis.text.x = element_text(size=14),
                          panel.grid.major.y = element_blank(),
                          panel.grid.minor = element_blank()) + xlab('Percent of Population') + ylab('Density') + ggtitle('') +
  scale_x_continuous(labels = label_percent(), limits = c(0, 1))

den_18_49

den_50_64 <- ggplot(usa_albers_county_50_64) + geom_density(aes(perc_50), color = '#285f93', fill = '#285f93', alpha = 0.5,
                                                            linewidth = 1.5) + 
  geom_vline(xintercept=age_group[4, ]$perc_50, lty=2, linewidth = 1.3, color = '#285f93') +
  theme_minimal() + theme(axis.title = element_text(size=16),
                          strip.background = element_blank(),
                          axis.text.y = element_blank(),
                          axis.text.x = element_text(size=14),
                          panel.grid.major.y = element_blank(),
                          panel.grid.minor = element_blank()) + xlab('Percent of Population') + ylab('Density') + ggtitle('') +
  scale_x_continuous(labels = label_percent(), limits = c(0, 1))

den_50_64

den_65 <- ggplot(usa_albers_county_65) + geom_density(aes(perc_50), color = '#1b4164', fill = '#1b4164', alpha = 0.5,
                                                            linewidth = 1.5) + 
  geom_vline(xintercept=age_group[5, ]$perc_50, lty=2, linewidth = 1.3, color = '#1b4164') +
  theme_minimal() + theme(axis.title = element_text(size=16),
                          strip.background = element_blank(),
                          axis.text.y = element_blank(),
                          axis.text.x = element_text(size=14),
                          panel.grid.major.y = element_blank(),
                          panel.grid.minor = element_blank()) + xlab('Percent of Population') + ylab('Density') + ggtitle('') +
  scale_x_continuous(labels = label_percent(), limits = c(0, 1))

den_65

figure_5_map_age <- cowplot::plot_grid(map_65, den_65,  map_18_49, den_18_49, map_0_4, den_0_4,
                                   nrow = 3,
                                   labels = c("", "", "", "", "", ""),
                                   label_size = 20)

figure_5 <- cowplot::plot_grid(season_comp_plot, figure_5_map,
                               age_comp_plot, figure_5_map_age,
                               rel_widths = c(1, 1, 1, 1),
                               nrow = 2,
                               labels = c("a", "b", "c", "d"),
                               label_size = 20)

ggsave('./figs/figure_5_age.jpg', plot = figure_5, height = 16, width = 20)

flu_county_perc_filt <- left_join(usa_albers_county_18_49, as.data.frame(usa_albers_county_65), by = 'GEOID')

cor(flu_county_perc_filt$perc_50.x, flu_county_perc_filt$perc_50.y, use = 'complete.obs', method = 'spearman')





figure_5 <- cowplot::plot_grid(age_comp_plot, figure_5_map_age,
                               rel_widths = c(1, 1),
                               nrow = 1,
                               labels = c("a", ""),
                               label_size = 20)

ggsave('./figs/figure_5_age.jpg', plot = figure_5, height = 8, width = 20)










