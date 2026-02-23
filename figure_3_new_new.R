################################################################################
# File Name: figure_3                                                          #
#                                                                              #
# Purpose:   Create figure 3 for the manuscript.                               #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Create figure 3                                                #
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
library(arrow)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

######################
# 2. CREATE FIGURE 3 #
######################

# Load other surveillance data
ili_nat <- readRDS('./tmp/ili_net_national.rds')
ili_state <- readRDS('./tmp/ili_net_state.rds')
flu_surv_nat <- readRDS('./tmp/flu_surv_net_national.rds')
flu_surv_state <- readRDS('./tmp/flu_surv_net_state.rds')
test_pos <- readRDS('./tmp/test_positivity_county.rds')
nys_county <- readRDS('./tmp/nys_county.rds')

# LOAD DATA
# 2016
flu_pred_dat_2016 <- readRDS('./tmp/flu_2016_predictions')

flu_county_2016 <- flu_pred_dat_2016 |> 
  group_by(year_week_dt, county_fips) |>
  mutate(pred_sum = sum(pred_case_count, na.rm = T),
         flu_sum = sum(flu_count)) |>
  distinct(year_week_dt, county_fips, pred_sum, flu_sum)
remove(flu_pred_dat_2016)

# 2017
flu_pred_dat_2017 <- readRDS('./tmp/flu_2017_predictions')

flu_county_2017 <- flu_pred_dat_2017 |> 
  group_by(year_week_dt, county_fips) |>
  mutate(pred_sum = sum(pred_case_count, na.rm = T),
         flu_sum = sum(flu_count)) |>
  distinct(year_week_dt, county_fips, pred_sum, flu_sum)
remove(flu_pred_dat_2017)


# 2018
flu_pred_dat_2018 <- readRDS('./tmp/flu_2018_predictions')

flu_county_2018 <- flu_pred_dat_2018 |> 
  group_by(year_week_dt, county_fips) |>
  mutate(pred_sum = sum(pred_case_count, na.rm = T),
         flu_sum = sum(flu_count)) |>
  distinct(year_week_dt, county_fips, pred_sum, flu_sum)
remove(flu_pred_dat_2018)


# All cause county
# Load all cause data
all_cause <- read.csv('./raw/county_week_ac_v2.csv', header = T)

# Impute all cause counts
all_cause$patient_count_imp <- all_cause$all_cause
all_cause$patient_count_imp[all_cause$patient_count_imp == '<=5'] <- 
  sample(1:5, length(all_cause$patient_count_imp[all_cause$patient_count_imp == '<=5']), 
         replace= TRUE)
# Convert imputed count to numeric
all_cause$patient_count_imp <- as.numeric(all_cause$patient_count_imp)

# Convert year week to week date
all_cause <- all_cause |>
  mutate(week_date = sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", year_week)) |>
  mutate(week_date = as.Date(ISOweek2date(week_date)))

# Collapse to county level
all_cause_county <- all_cause |> group_by(week_date, county_fips) |>
  mutate(all_cause_sum = sum(patient_count_imp)) |>
  distinct(week_date, county_fips, all_cause_sum) |>
  ungroup() |>
  dplyr::filter(week_date < as.Date('2019-09-01')) |>
  dplyr::filter(week_date > as.Date('2016-08-31')) |>
  group_by(county_fips) |>
  mutate(all_cause_avg = mean(all_cause_sum))

# Merge on all cause county
flu_county_2016 <- left_join(flu_county_2016, all_cause_county, by = c('county_fips','year_week_dt' = 'week_date'))
flu_county_2017 <- left_join(flu_county_2017, all_cause_county, by = c('county_fips','year_week_dt' = 'week_date'))
flu_county_2018 <- left_join(flu_county_2018, all_cause_county, by = c('county_fips','year_week_dt' = 'week_date'))


# Adjust count
flu_county_2016$flu_obs_scale <- flu_county_2016$flu_sum *
  (1/(flu_county_2016$all_cause_sum / flu_county_2016$all_cause_avg))
flu_county_2016$flu_pred_scale <- flu_county_2016$pred_sum *
  (1/(flu_county_2016$all_cause_sum / flu_county_2016$all_cause_avg))

flu_county_2017$flu_obs_scale <- flu_county_2017$flu_sum *
  (1/(flu_county_2017$all_cause_sum / flu_county_2017$all_cause_avg))
flu_county_2017$flu_pred_scale <- flu_county_2017$pred_sum *
  (1/(flu_county_2017$all_cause_sum / flu_county_2017$all_cause_avg))

flu_county_2018$flu_obs_scale <- flu_county_2018$flu_sum *
  (1/(flu_county_2018$all_cause_sum / flu_county_2018$all_cause_avg))
flu_county_2018$flu_pred_scale <- flu_county_2018$pred_sum *
  (1/(flu_county_2018$all_cause_sum / flu_county_2018$all_cause_avg))


# Add year and month variables
flu_county_2016$year <- year(flu_county_2016$year_week_dt)
flu_county_2016$month <- month(flu_county_2016$year_week_dt)

flu_county_2017$year <- year(flu_county_2017$year_week_dt)
flu_county_2017$month <- month(flu_county_2017$year_week_dt)

flu_county_2018$year <- year(flu_county_2018$year_week_dt)
flu_county_2018$month <- month(flu_county_2018$year_week_dt)


# Calculate summer means
flu_county_summer_2017 <- flu_county_2016 |> 
  dplyr::filter(month < 8 & month > 5) |> 
  group_by(year, county_fips) |>
  dplyr::filter(year > 2016) |>
  mutate(pred_summ_mean = mean(pred_sum),
         obs_summ_mean = mean(flu_sum)) |>
  distinct(year, county_fips, pred_summ_mean, 
           obs_summ_mean)

flu_county_summer_2018 <- flu_county_2017 |> 
  dplyr::filter(month < 8 & month > 5) |> 
  group_by(year, county_fips) |>
  dplyr::filter(year > 2017) |>
  mutate(pred_summ_mean = mean(pred_sum),
         obs_summ_mean = mean(flu_sum)) |>
  distinct(year, county_fips, pred_summ_mean, 
           obs_summ_mean)


# Merge on summer data to the correct season (summer 2016 should align with the 2016-2017 flu season)
flu_county_2016$year_2 <- ifelse(flu_county_2016$month < 6, 2017, 2017)
flu_county_2017$year_2 <- ifelse(flu_county_2017$month < 6, 2017, 2017)
flu_county_2018$year_2 <- ifelse(flu_county_2018$month < 6, 2018, 2018)

flu_county_2016 <- left_join(flu_county_2016, flu_county_summer_2017[, c(1, 2, 3, 4)], by = c('year_2' = 'year', 'county_fips' = 'county_fips'))
flu_county_2017 <- left_join(flu_county_2017, flu_county_summer_2017[, c(1, 2, 3, 4)], by = c('year_2' = 'year', 'county_fips' = 'county_fips'))
flu_county_2018 <- left_join(flu_county_2018, flu_county_summer_2018[, c(1, 2, 3, 4)], by = c('year_2' = 'year', 'county_fips' = 'county_fips'))

# Calculate Z score
flu_county_2016 <- flu_county_2016 |> group_by(county_fips) |>
  mutate(pred_z = (pred_sum - pred_summ_mean) / sd(pred_sum, na.rm = T),
         obs_z = (flu_sum - obs_summ_mean) / sd(flu_sum, na.rm = T))

flu_county_2017 <- flu_county_2017 |> group_by(county_fips) |>
  mutate(pred_z = (pred_sum - pred_summ_mean) / sd(pred_sum, na.rm = T),
         obs_z = (flu_sum - obs_summ_mean) / sd(flu_sum, na.rm = T))

flu_county_2018 <- flu_county_2018 |> group_by(county_fips) |>
  mutate(pred_z = (pred_sum - pred_summ_mean) / sd(pred_sum, na.rm = T),
         obs_z = (flu_sum - obs_summ_mean) / sd(flu_sum, na.rm = T))

# CALCULATE PEAKS #
flu_county_peak_2016 <- flu_county_2016
flu_county_peak_2017 <- flu_county_2017
flu_county_peak_2018 <- flu_county_2018

flu_county_peak_2016$Season <- '2016-2017'
flu_county_peak_2017$Season <- '2017-2018'
flu_county_peak_2018$Season <- '2018-2019'

flu_county_peak_2016 <- flu_county_peak_2016 |> arrange(year_week_dt) |> group_by(county_fips) |> 
  mutate(flu_sum_roll = rollmean(flu_sum, k = 4, align = 'right', na.pad = TRUE),
         flu_pred_roll = rollmean(pred_sum, k = 4, align = 'right', na.pad = TRUE))

flu_county_peak_2017 <- flu_county_peak_2017 |> arrange(year_week_dt) |> group_by(county_fips) |> 
  mutate(flu_sum_roll = rollmean(flu_sum, k = 4, align = 'right', na.pad = TRUE),
         flu_pred_roll = rollmean(pred_sum, k = 4, align = 'right', na.pad = TRUE))

flu_county_peak_2018 <- flu_county_peak_2018 |> arrange(year_week_dt) |> group_by(county_fips) |> 
  mutate(flu_sum_roll = rollmean(flu_sum, k = 4, align = 'right', na.pad = TRUE),
         flu_pred_roll = rollmean(pred_sum, k = 4, align = 'right', na.pad = TRUE))

# Load filtering data
all_cause_season <- read_csv('raw/county_season_ac_v2.csv')

county_pop <- read.csv('./raw/co-est2019-alldata.csv', stringsAsFactors=FALSE)
county_pop$fips_state <- ifelse(nchar(county_pop$STATE) == 1, paste0("0", as.character(county_pop$STATE)), as.character(county_pop$STATE))
county_pop$fips_county <- ifelse(nchar(county_pop$COUNTY) == 1, paste0("00", as.character(county_pop$COUNTY)), as.character(county_pop$COUNTY))
county_pop$fips_county <- ifelse(nchar(county_pop$COUNTY) == 2, paste0("0", as.character(county_pop$COUNTY)), county_pop$fips_county)
county_pop$fips <- paste(county_pop$fips_state, county_pop$fips_county, sep = '')


# Merge on filtering data
flu_county_peak_2016 <- left_join(flu_county_peak_2016, all_cause_season, by = c('Season' = 'season', 'county_fips'))
flu_county_peak_2017 <- left_join(flu_county_peak_2017, all_cause_season, by = c('Season' = 'season', 'county_fips'))
flu_county_peak_2018 <- left_join(flu_county_peak_2018, all_cause_season, by = c('Season' = 'season', 'county_fips'))

flu_county_peak_2016 <- left_join(flu_county_peak_2016, county_pop, by = c('county_fips' = 'fips'))
flu_county_peak_2017 <- left_join(flu_county_peak_2017, county_pop, by = c('county_fips' = 'fips'))
flu_county_peak_2018 <- left_join(flu_county_peak_2018, county_pop, by = c('county_fips' = 'fips'))

# Limit data based on county pop and number of claims in a season
flu_county_peak_filt_2016 <- flu_county_peak_2016 |>
  dplyr::filter(POPESTIMATE2019 > 10000) |>
  dplyr::filter(all_cause > 5000)

flu_county_peak_filt_2017 <- flu_county_peak_2017 |>
  dplyr::filter(POPESTIMATE2019 > 10000) |>
  dplyr::filter(all_cause > 5000)

flu_county_peak_filt_2018 <- flu_county_peak_2018 |>
  dplyr::filter(POPESTIMATE2019 > 10000) |>
  dplyr::filter(all_cause > 5000)

# FIND PEAKS

# Calculate peaks
county_peaks_usa_2016 <- flu_county_peak_filt_2016 |> 
  ungroup() |>
  dplyr::filter(!is.na(flu_sum_roll)) |>
  dplyr::filter(!is.na(flu_pred_roll)) |>
  # Calculate the row number
  mutate(index = row_number()) |>
  group_by(Season, county_fips) |> 
  # Calculate multiple peaks greater than 10 cases and larger than 10 neighbors
  mutate(peaks_obs = seq_along(index) %in% 
           findpeaks(flu_sum_roll, MinPeakDistance = 10, MinPeakHeight = 20)$loc) |>
  mutate(peaks_pred = seq_along(index) %in% 
           findpeaks(flu_pred_roll, MinPeakDistance = 10, MinPeakHeight = 20)$loc) |>
  # Remove first and last observations, which cannot be peaks but might be max
  dplyr::filter(row_number() != 1 & row_number() != n()) |>
  # Calculate the maximum value among the peaks
  mutate(is_max_obs = flu_sum_roll == max(flu_sum_roll),
         is_max_pred = flu_pred_roll == max(flu_pred_roll)) |>
  select(c(Season, county_fips, year_week_dt, flu_pred_roll, flu_sum_roll, obs_z, pred_z,
           peaks_obs, peaks_pred, is_max_obs, is_max_pred))

county_peaks_usa_2017 <- flu_county_peak_filt_2017 |> 
  ungroup() |>
  dplyr::filter(!is.na(flu_sum_roll)) |>
  dplyr::filter(!is.na(flu_pred_roll)) |>
  # Calculate the row number
  mutate(index = row_number()) |>
  group_by(Season, county_fips) |> 
  # Calculate multiple peaks greater than 10 cases and larger than 10 neighbors
  mutate(peaks_obs = seq_along(index) %in% 
           findpeaks(flu_sum_roll, MinPeakDistance = 10, MinPeakHeight = 20)$loc) |>
  mutate(peaks_pred = seq_along(index) %in% 
           findpeaks(flu_pred_roll, MinPeakDistance = 10, MinPeakHeight = 20)$loc) |>
  # Remove first and last observations, which cannot be peaks but might be max
  dplyr::filter(row_number() != 1 & row_number() != n()) |>
  # Calculate the maximum value among the peaks
  mutate(is_max_obs = flu_sum_roll == max(flu_sum_roll),
         is_max_pred = flu_pred_roll == max(flu_pred_roll)) |>
  select(c(Season, county_fips, year_week_dt, flu_pred_roll, flu_sum_roll, obs_z, pred_z,
           peaks_obs, peaks_pred, is_max_obs, is_max_pred))

county_peaks_usa_2018 <- flu_county_peak_filt_2018 |> 
  ungroup() |>
  dplyr::filter(!is.na(flu_sum_roll)) |>
  dplyr::filter(!is.na(flu_pred_roll)) |>
  # Calculate the row number
  mutate(index = row_number()) |>
  group_by(Season, county_fips) |> 
  # Calculate multiple peaks greater than 10 cases and larger than 10 neighbors
  mutate(peaks_obs = seq_along(index) %in% 
           findpeaks(flu_sum_roll, MinPeakDistance = 10, MinPeakHeight = 20)$loc) |>
  mutate(peaks_pred = seq_along(index) %in% 
           findpeaks(flu_pred_roll, MinPeakDistance = 10, MinPeakHeight = 20)$loc) |>
  # Remove first and last observations, which cannot be peaks but might be max
  dplyr::filter(row_number() != 1 & row_number() != n()) |>
  # Calculate the maximum value among the peaks
  mutate(is_max_obs = flu_sum_roll == max(flu_sum_roll),
         is_max_pred = flu_pred_roll == max(flu_pred_roll)) |>
  select(c(Season, county_fips, year_week_dt, flu_pred_roll, flu_sum_roll, obs_z, pred_z,
           peaks_obs, peaks_pred, is_max_obs, is_max_pred))

county_peaks_obs_2016 <- county_peaks_usa_2016 |> 
  # Select the peaks
  dplyr::filter(peaks_obs == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_obs == T) |>
  select(c(Season, county_fips, year_week_dt, peaks_obs, is_max_obs, obs_z)) |>
  rename('week_obs' = 'year_week_dt') |>
  group_by(Season, county_fips) |>
  # Break any ties by taking the earliest peak
  arrange(week_obs) |>
  slice_head(n = 1) |>
  mutate(week_obs_num = week(week_obs)) #|>
  # Filter for significant peaks w/ z score > 1
  #dplyr::filter(obs_z > 1)

county_peaks_obs_2017 <- county_peaks_usa_2017 |> 
  # Select the peaks
  dplyr::filter(peaks_obs == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_obs == T) |>
  select(c(Season, county_fips, year_week_dt, peaks_obs, is_max_obs, obs_z)) |>
  rename('week_obs' = 'year_week_dt') |>
  group_by(Season, county_fips) |>
  # Break any ties by taking the earliest peak
  arrange(week_obs) |>
  slice_head(n = 1) |>
  mutate(week_obs_num = week(week_obs)) #|>
  # Filter for significant peaks w/ z score > 1
  #dplyr::filter(obs_z > 1)

county_peaks_obs_2018 <- county_peaks_usa_2018 |> 
  # Select the peaks
  dplyr::filter(peaks_obs == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_obs == T) |>
  select(c(Season, county_fips, year_week_dt, peaks_obs, is_max_obs, obs_z)) |>
  rename('week_obs' = 'year_week_dt') |>
  group_by(Season, county_fips) |>
  # Break any ties by taking the earliest peak
  arrange(week_obs) |>
  slice_head(n = 1) |>
  mutate(week_obs_num = week(week_obs)) #|>
  # Filter for significant peaks w/ z score > 1
  #dplyr::filter(obs_z > 1)

county_peaks_pred_2016 <- county_peaks_usa_2016 |> 
  # Select the peaks
  dplyr::filter(peaks_pred == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_pred == T) |>
  select(c(Season, county_fips, year_week_dt, peaks_pred, is_max_pred, pred_z)) |>
  rename('week_pred' = 'year_week_dt') |>
  group_by(Season, county_fips) |>
  # Break any ties by taking the earliest peak
  arrange(week_pred) |>
  slice_head(n = 1) |>
  mutate(week_pred_num = week(week_pred)) #|>
  # Filter for significant peaks w/ z score > 1
  #dplyr::filter(pred_z > 1)

county_peaks_pred_2017 <- county_peaks_usa_2017 |> 
  # Select the peaks
  dplyr::filter(peaks_pred == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_pred == T) |>
  select(c(Season, county_fips, year_week_dt, peaks_pred, is_max_pred, pred_z)) |>
  rename('week_pred' = 'year_week_dt') |>
  group_by(Season, county_fips) |>
  # Break any ties by taking the earliest peak
  arrange(week_pred) |>
  slice_head(n = 1) |>
  mutate(week_pred_num = week(week_pred)) #|>
  # Filter for significant peaks w/ z score > 1
  #dplyr::filter(pred_z > 1)

county_peaks_pred_2018 <- county_peaks_usa_2018 |> 
  # Select the peaks
  dplyr::filter(peaks_pred == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_pred == T) |>
  select(c(Season, county_fips, year_week_dt, peaks_pred, is_max_pred, pred_z)) |>
  rename('week_pred' = 'year_week_dt') |>
  group_by(Season, county_fips) |>
  # Break any ties by taking the earliest peak
  arrange(week_pred) |>
  slice_head(n = 1) |>
  mutate(week_pred_num = week(week_pred)) #|>
  # Filter for significant peaks w/ z score > 1
  #dplyr::filter(pred_z > 1)

# Merge obs and pred peaks
county_peaks_both_2016 <- left_join(county_peaks_pred_2016, county_peaks_obs_2016,
                               by = c('Season' = 'Season',
                                      'county_fips' = 'county_fips'))

county_peaks_both_2017 <- left_join(county_peaks_pred_2017, county_peaks_obs_2017,
                               by = c('Season' = 'Season',
                                      'county_fips' = 'county_fips'))

county_peaks_both_2018 <- left_join(county_peaks_pred_2018, county_peaks_obs_2018,
                               by = c('Season' = 'Season',
                                      'county_fips' = 'county_fips'))

# Filter NAs
county_peaks_both_2016_filt <- county_peaks_both_2016 |>
  dplyr::filter(!is.na(week_pred)) |>
  dplyr::filter(!is.na(week_obs))

county_peaks_both_2017_filt <- county_peaks_both_2017 |>
  dplyr::filter(!is.na(week_pred)) |>
  dplyr::filter(!is.na(week_obs))
  
county_peaks_both_2018_filt <- county_peaks_both_2018 |>
  dplyr::filter(!is.na(week_pred)) |>
  dplyr::filter(!is.na(week_obs))

# Calculate statistics
county_peaks_both_2016_filt$diff <- as.numeric(county_peaks_both_2016_filt$week_pred - county_peaks_both_2016_filt$week_obs)
county_peaks_both_2017_filt$diff <- as.numeric(county_peaks_both_2017_filt$week_pred - county_peaks_both_2017_filt$week_obs)
county_peaks_both_2018_filt$diff <- as.numeric(county_peaks_both_2018_filt$week_pred - county_peaks_both_2018_filt$week_obs)

rmse_2016 <- sqrt(sum((county_peaks_both_2016_filt$diff)^2, na.rm = T) / length(county_peaks_both_2016_filt$diff)) 
rmse_2017 <- sqrt(sum((county_peaks_both_2017_filt$diff)^2, na.rm = T) / length(county_peaks_both_2017_filt$diff)) 
rmse_2018 <- sqrt(sum((county_peaks_both_2018_filt$diff)^2, na.rm = T) / length(county_peaks_both_2018_filt$diff)) 

rmse_2016
rmse_2017
rmse_2018

county_peaks_both_2018_filt <- county_peaks_both_2018_filt |>
  mutate(week_obs_num_adj = ifelse(week_obs_num > 30, week_obs_num - 52, week_obs_num),
         week_pred_num_adj = ifelse(week_pred_num > 30, week_pred_num - 52, week_pred_num))
county_peaks_both_2017_filt <- county_peaks_both_2017_filt |>
  mutate(week_obs_num_adj = ifelse(week_obs_num > 30, week_obs_num - 52, week_obs_num),
         week_pred_num_adj = ifelse(week_pred_num > 30, week_pred_num - 52, week_pred_num))
county_peaks_both_2016_filt <- county_peaks_both_2016_filt |>
  mutate(week_obs_num_adj = ifelse(week_obs_num > 30, week_obs_num - 52, week_obs_num),
         week_pred_num_adj = ifelse(week_pred_num > 30, week_pred_num - 52, week_pred_num))


cor(county_peaks_both_2016_filt$week_pred_num_adj, county_peaks_both_2016_filt$week_obs_num_adj, 
    use = 'complete.obs', method = 'spearman')
cor(county_peaks_both_2017_filt$week_pred_num_adj, county_peaks_both_2017_filt$week_obs_num_adj, 
    use = 'complete.obs', method = 'spearman')
cor(county_peaks_both_2018_filt$week_pred_num_adj, county_peaks_both_2018_filt$week_obs_num_adj, 
    use = 'complete.obs', method = 'spearman')

# Collapse data for graphing purposes
county_peaks_both_col_2016 <- county_peaks_both_2016_filt |>
  mutate(count = 1) |>
  group_by(week_obs_num, week_pred_num) |>
  mutate(count_sum = sum(count)) |>
  distinct(week_pred, week_obs, week_obs_num, week_pred_num, count_sum)

county_peaks_both_col_2017 <- county_peaks_both_2017_filt |>
  mutate(count = 1) |>
  group_by(week_obs_num, week_pred_num) |>
  mutate(count_sum = sum(count)) |>
  distinct(week_pred, week_obs, week_obs_num, week_pred_num, count_sum)

county_peaks_both_col_2018 <- county_peaks_both_2018_filt |>
  mutate(count = 1) |>
  group_by(week_obs_num, week_pred_num) |>
  mutate(count_sum = sum(count)) |>
  distinct(week_pred, week_obs, week_obs_num, week_pred_num, count_sum)

forrest_theme <- theme(legend.position =  "bottom",
                       legend.position.inside = c(0.873, 0.237),
                       axis.text = element_text(size=14, color = 'black'),
                       axis.title = element_text(size=16),
                       strip.background = element_blank(),
                       legend.title = element_text(size=16),
                       legend.text = element_text(size=14),
                       plot.title = element_text(size=18, hjust = 0),
                       legend.box="horizontal",
                       legend.box.background = element_rect(colour = "white", 
                                                            fill = 'white'),
                       legend.key.width = unit(1, "cm"))


peak_comp_plot_2016 <- ggplot(county_peaks_both_col_2016[!is.na(county_peaks_both_col_2016$week_obs),], aes(week_obs, week_pred, fill= count_sum)) + 
  geom_tile(color = 'black', linewidth = 0.25) + 
  scale_fill_gradient('Number of Counties\n', high = '#3cbb75ff', low = 'white') +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black",
              linewidth = 1) +
  theme_minimal() + 
  labs(title = 'Peak Time Comparison\n2016-2017 Season',
       x = "Claims Peak Week",
       y = "Syndromic Peak Week") + forrest_theme +
  scale_x_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2016-11-01'), as.Date('2017-05-27'))) +
  scale_y_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2016-11-01'), as.Date('2017-05-27')))
peak_comp_plot_2016

peak_comp_plot_2017 <- ggplot(county_peaks_both_col_2017[!is.na(county_peaks_both_col_2017$week_obs),], aes(week_obs, week_pred, fill= count_sum)) + 
  geom_tile(color = 'black', linewidth = 0.25) + 
  scale_fill_gradient('Number of Counties\n', low = 'white', high = '#d7642c') +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black",
              linewidth = 1) +
  theme_minimal() + 
  labs(title = 'Peak Time Comparison\n2017-2018 Season',
       x = "Claims Peak Week",
       y = "Syndromic Peak Week") + forrest_theme +
  scale_x_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2017-11-01'), as.Date('2018-05-27'))) +
  scale_y_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2017-11-01'), as.Date('2018-05-27')))
peak_comp_plot_2017

forrest_theme <- theme(legend.position =  "right",
                       legend.position.inside = c(0.873, 0.237),
                       axis.text = element_text(size=14, color = 'black'),
                       axis.title = element_text(size=16),
                       strip.background = element_blank(),
                       legend.title = element_text(size=16),
                       legend.text = element_text(size=14),
                       plot.title = element_text(size=18, hjust = 0),
                       legend.box="horizontal",
                       legend.box.background = element_rect(colour = "white", 
                                                            fill = 'white'),
                       legend.key.width = unit(0.5, "cm"),
                       legend.key.height = unit(0.75, "cm"))


peak_comp_plot_2018 <- ggplot(county_peaks_both_col_2018[!is.na(county_peaks_both_col_2018$week_obs),], aes(week_obs, week_pred, fill= count_sum)) + 
  geom_tile(color = 'black', linewidth = 0.25) + 
  scale_fill_gradient('Number of \nCounties', low = 'white', high = '#9B59B6') +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black",
              linewidth = 1) +
  theme_minimal() + 
  labs(title = 'Peak Time Comparison\n2018-2019 Season',
       x = "Claims Peak Week",
       y = "Syndromic Peak Week") + forrest_theme +
  scale_x_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2018-11-01'), as.Date('2019-05-27'))) +
  scale_y_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2018-11-01'), as.Date('2019-05-27')))
peak_comp_plot_2018

ggsave('./figs/peak_comp.jpg', plot = peak_comp_plot_2018, height = 4, width = 7)


# STATE LEVEL #

flu_state_2016 <- flu_pred_dat_2016 |> ungroup() |>
  group_by(year_week_dt, state_fips) |>
  mutate(pred_sum = sum(pred_count, na.rm = T),
         flu_sum = sum(patient_count * flu)) |>
  distinct(year_week_dt, state_fips, pred_sum, flu_sum)

flu_state_2017 <- flu_pred_dat_2017 |> ungroup() |>
  group_by(year_week_dt, state_fips) |>
  mutate(pred_sum = sum(pred_count, na.rm = T),
         flu_sum = sum(patient_count * flu)) |>
  distinct(year_week_dt, state_fips, pred_sum, flu_sum)

flu_state_2018 <- flu_pred_dat_2018 |> ungroup() |>
  group_by(year_week_dt, state_fips) |>
  mutate(pred_sum = sum(pred_count, na.rm = T),
         flu_sum = sum(patient_count * flu)) |>
  distinct(year_week_dt, state_fips, pred_sum, flu_sum)

# Set year and month variables
flu_state_2016$month <- month(flu_state_2016$year_week_dt)
flu_state_2016$year <- year(flu_state_2016$year_week_dt)

flu_state_2017$month <- month(flu_state_2017$year_week_dt)
flu_state_2017$year <- year(flu_state_2017$year_week_dt)

flu_state_2018$month <- month(flu_state_2018$year_week_dt)
flu_state_2018$year <- year(flu_state_2018$year_week_dt)

# Calculate summer averages
flu_state_summer_2017 <- flu_state_2016 |> 
  dplyr::filter(month < 9 & month > 5) |> group_by(year, state_fips) |>
  dplyr::filter(year > 2016) |>
  mutate(pred_summ_mean = mean(pred_sum),
         obs_summ_mean = mean(flu_sum)) |>
  distinct(year, state_fips, pred_summ_mean, 
           obs_summ_mean)

flu_state_summer_2018 <- flu_state_2017 |> 
  dplyr::filter(month < 9 & month > 5) |> group_by(year, state_fips) |>
  dplyr::filter(year > 2016) |>
  mutate(pred_summ_mean = mean(pred_sum),
         obs_summ_mean = mean(flu_sum)) |>
  distinct(year, state_fips, pred_summ_mean, 
           obs_summ_mean)

# Merge on summer data to the correct season (summer 2016 should align with the 2016-2017 flu season)
flu_state_2016$year_2 <- ifelse(flu_state_2016$month < 6, 2017, 2017)
flu_state_2017$year_2 <- ifelse(flu_state_2017$month < 6, 2017, 2017)
flu_state_2018$year_2 <- ifelse(flu_state_2018$month < 6, 2018, 2018)

flu_state_2016 <- left_join(flu_state_2016, flu_state_summer_2017, by = c('year_2' = 'year', 'state_fips' = 'state_fips'))
flu_state_2017 <- left_join(flu_state_2017, flu_state_summer_2017, by = c('year_2' = 'year', 'state_fips' = 'state_fips'))
flu_state_2018 <- left_join(flu_state_2018, flu_state_summer_2018, by = c('year_2' = 'year', 'state_fips' = 'state_fips'))

flu_state_2016 <- flu_state_2016 |> arrange(year_week_dt) |> group_by(state_fips) |> 
  mutate(flu_sum_roll = rollmean(flu_sum, k = 4, align = 'right', na.pad = TRUE),
         flu_pred_roll = rollmean(pred_sum, k = 4, align = 'right', na.pad = TRUE))
flu_state_2017 <- flu_state_2017 |> arrange(year_week_dt) |> group_by(state_fips) |> 
  mutate(flu_sum_roll = rollmean(flu_sum, k = 4, align = 'right', na.pad = TRUE),
         flu_pred_roll = rollmean(pred_sum, k = 4, align = 'right', na.pad = TRUE))
flu_state_2018 <- flu_state_2018 |> arrange(year_week_dt) |> group_by(state_fips) |> 
  mutate(flu_sum_roll = rollmean(flu_sum, k = 4, align = 'right', na.pad = TRUE),
         flu_pred_roll = rollmean(pred_sum, k = 4, align = 'right', na.pad = TRUE))

# Calculate Z score
flu_state_2016$z_score <- (flu_state_2016$pred_sum - flu_state_2016$pred_summ_mean) / sd(flu_state_2016$pred_sum, na.rm = T)
flu_state_2017$z_score <- (flu_state_2017$pred_sum - flu_state_2017$pred_summ_mean) / sd(flu_state_2017$pred_sum, na.rm = T)
flu_state_2018$z_score <- (flu_state_2018$pred_sum - flu_state_2018$pred_summ_mean) / sd(flu_state_2018$pred_sum, na.rm = T)

flu_state_2016$week_date <- flu_state_2016$year_week_dt
flu_state_2017$week_date <- flu_state_2017$year_week_dt
flu_state_2018$week_date <- flu_state_2018$year_week_dt

state_peaks_2016 <- flu_state_2016 |> 
  dplyr::filter(!is.na(flu_sum_roll)) |>
  ungroup() |>
  # Calculate the row number
  mutate(index = row_number()) |>
  group_by(state_fips) |> 
  # Calculate multiple peaks greater than 10 cases and larger than 10 neighbors
  mutate(peaks_obs = seq_along(index) %in% 
           findpeaks(flu_sum_roll, MinPeakDistance = 10, MinPeakHeight = 20)$loc) |>
  mutate(peaks_pred = seq_along(index) %in% 
           findpeaks(flu_pred_roll, MinPeakDistance = 10, MinPeakHeight = 20)$loc) |>
  # Remove first and last observations, which cannot be peaks but might be max
  dplyr::filter(row_number() != 1 & row_number() != n()) |>
  # Calculate the maximum value among the peaks
  mutate(is_max_obs = flu_sum_roll == max(flu_sum_roll)) |>
  mutate(is_max_pred = flu_pred_roll == max(flu_pred_roll)) |>
  select(c(state_fips, week_date, flu_sum_roll,flu_pred_roll,  peaks_obs, peaks_pred, is_max_obs, is_max_pred))

state_peaks_2017 <- flu_state_2017 |> 
  dplyr::filter(!is.na(flu_sum_roll)) |>
  ungroup() |>
  # Calculate the row number
  mutate(index = row_number()) |>
  group_by(state_fips) |> 
  # Calculate multiple peaks greater than 10 cases and larger than 10 neighbors
  mutate(peaks_obs = seq_along(index) %in% 
           findpeaks(flu_sum_roll, MinPeakDistance = 10, MinPeakHeight = 20)$loc) |>
  mutate(peaks_pred = seq_along(index) %in% 
           findpeaks(flu_pred_roll, MinPeakDistance = 10, MinPeakHeight = 20)$loc) |>
  # Remove first and last observations, which cannot be peaks but might be max
  dplyr::filter(row_number() != 1 & row_number() != n()) |>
  # Calculate the maximum value among the peaks
  mutate(is_max_obs = flu_sum_roll == max(flu_sum_roll)) |>
  mutate(is_max_pred = flu_pred_roll == max(flu_pred_roll)) |>
  select(c(state_fips, week_date, flu_sum_roll,flu_pred_roll,  peaks_obs, peaks_pred, is_max_obs, is_max_pred))

state_peaks_2018 <- flu_state_2018 |> 
  dplyr::filter(!is.na(flu_sum_roll)) |>
  ungroup() |>
  # Calculate the row number
  mutate(index = row_number()) |>
  group_by(state_fips) |> 
  # Calculate multiple peaks greater than 10 cases and larger than 10 neighbors
  mutate(peaks_obs = seq_along(index) %in% 
           findpeaks(flu_sum_roll, MinPeakDistance = 10, MinPeakHeight = 20)$loc) |>
  mutate(peaks_pred = seq_along(index) %in% 
           findpeaks(flu_pred_roll, MinPeakDistance = 10, MinPeakHeight = 20)$loc) |>
  # Remove first and last observations, which cannot be peaks but might be max
  dplyr::filter(row_number() != 1 & row_number() != n()) |>
  # Calculate the maximum value among the peaks
  mutate(is_max_obs = flu_sum_roll == max(flu_sum_roll)) |>
  mutate(is_max_pred = flu_pred_roll == max(flu_pred_roll)) |>
  select(c(state_fips, week_date, flu_sum_roll,flu_pred_roll,  peaks_obs, peaks_pred, is_max_obs, is_max_pred))

state_peaks_obs_2016 <- state_peaks_2016 |> 
  # Select the peaks
  dplyr::filter(peaks_obs == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_obs == T) |>
  select(c(state_fips, week_date, flu_sum_roll, peaks_obs, is_max_obs)) |>
  rename('week_state' = 'week_date') |>
  group_by(state_fips) |>
  # Break any ties by taking the earliest peak
  arrange(week_state) |>
  slice_head(n = 1) |>
  mutate(week_obs_num = week(week_state)) 

state_peaks_obs_2017 <- state_peaks_2017 |> 
  # Select the peaks
  dplyr::filter(peaks_obs == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_obs == T) |>
  select(c(state_fips, week_date, flu_sum_roll, peaks_obs, is_max_obs)) |>
  rename('week_state' = 'week_date') |>
  group_by(state_fips) |>
  # Break any ties by taking the earliest peak
  arrange(week_state) |>
  slice_head(n = 1) |>
  mutate(week_obs_num = week(week_state)) 

state_peaks_obs_2018 <- state_peaks_2018 |> 
  # Select the peaks
  dplyr::filter(peaks_obs == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_obs == T) |>
  select(c(state_fips, week_date, flu_sum_roll, peaks_obs, is_max_obs)) |>
  rename('week_state' = 'week_date') |>
  group_by(state_fips) |>
  # Break any ties by taking the earliest peak
  arrange(week_state) |>
  slice_head(n = 1) |>
  mutate(week_obs_num = week(week_state)) 

state_peaks_pred_2016 <- state_peaks_2016 |> 
  # Select the peaks
  dplyr::filter(peaks_pred == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_pred == T) |>
  select(c(state_fips, week_date, flu_pred_roll, peaks_pred, is_max_pred)) |>
  rename('week_pred' = 'week_date') |>
  group_by(state_fips) |>
  # Break any ties by taking the earliest peak
  arrange(week_pred) |>
  slice_head(n = 1) |>
  mutate(week_pred_num = week(week_pred)) |>
  dplyr::select(state_fips, week_pred, week_pred_num)

state_peaks_pred_2017 <- state_peaks_2017 |> 
  # Select the peaks
  dplyr::filter(peaks_pred == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_pred == T) |>
  select(c(state_fips, week_date, flu_pred_roll, peaks_pred, is_max_pred)) |>
  rename('week_pred' = 'week_date') |>
  group_by(state_fips) |>
  # Break any ties by taking the earliest peak
  arrange(week_pred) |>
  slice_head(n = 1) |>
  mutate(week_pred_num = week(week_pred)) |>
  dplyr::select(state_fips, week_pred, week_pred_num)

state_peaks_pred_2018 <- state_peaks_2018 |> 
  # Select the peaks
  dplyr::filter(peaks_pred == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_pred == T) |>
  select(c(state_fips, week_date, flu_pred_roll, peaks_pred, is_max_pred)) |>
  rename('week_pred' = 'week_date') |>
  group_by(state_fips) |>
  # Break any ties by taking the earliest peak
  arrange(week_pred) |>
  slice_head(n = 1) |>
  mutate(week_pred_num = week(week_pred)) |>
  dplyr::select(state_fips, week_pred, week_pred_num)

county_peaks_pred_2016$state_fips <- substr(county_peaks_pred_2016$county_fips, 1, 2)
early_peaks_2016 <- left_join(county_peaks_pred_2016, state_peaks_obs_2016, by = 'state_fips')
early_peaks_2016$Season <- '2016-2017'

county_peaks_pred_2017$state_fips <- substr(county_peaks_pred_2017$county_fips, 1, 2)
early_peaks_2017 <- left_join(county_peaks_pred_2017, state_peaks_obs_2017, by = 'state_fips')
early_peaks_2017$Season <- '2017-2018'

county_peaks_pred_2018$state_fips <- substr(county_peaks_pred_2018$county_fips, 1, 2)
early_peaks_2018 <- left_join(county_peaks_pred_2018, state_peaks_obs_2018, by = 'state_fips')
early_peaks_2018$Season <- '2018-2019'

early_peaks <- rbind(early_peaks_2016, early_peaks_2017, early_peaks_2018)

early_peaks <- early_peaks |> group_by(state_fips, Season) |>
  mutate(value = as.numeric(week_state - week_pred)) |>
  dplyr::filter(!is.na(value)) |>
  dplyr::filter(value == max(value)) |>
  arrange(week_state) |>
  slice_head(n = 1) |>
  dplyr::select(c(state_fips, Season, value)) |>
  mutate(Measure = 'Peak Time')

forrest_theme <- theme(legend.position =  "none",
                       legend.position.inside = c(0.873, 0.237),
                       axis.text = element_text(size=14, color = 'black'),
                       axis.title = element_text(size=16),
                       strip.background = element_blank(),
                       legend.title = element_text(size=16),
                       legend.text = element_text(size=14),
                       plot.title = element_text(size=18, hjust = 0),
                       legend.box="horizontal",
                       legend.box.background = element_rect(colour = "white", 
                                                            fill = 'white'),
                       legend.key.width = unit(1.5, "cm")) 



measure_plot <- ggplot(early_peaks) + 
  geom_boxplot(aes(y = Season, x = value, fill = Season), width=0.35, linewidth = 0.8,
               outlier.size = 2.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black",
             linewidth = 1) +
  theme_minimal() +
  xlim(-25, 65) +
  ggtitle('Syndromic vs. Observed \nPeak Lead Time') + 
  labs(
       x = "Number of Days",
       y = "Influenza Season") + forrest_theme +
  scale_fill_manual(values = c('2016-2017' = '#3cbb75ff' , '2017-2018' = '#d7642c', '2018-2019' = '#9B59B6'))

measure_plot


ggsave('./figs/measure_plot.jpg', plot = measure_plot, height = 4, width = 7)


early_peaks |> group_by(Season) |>
  summarise(mean(value))


# EXTRERNAL COMPARISON!!!

flu_state <- rbind(flu_state_2016, flu_state_2017, flu_state_2018)

flu_surv_state_avg <- flu_surv_state |> group_by(state_fips, week_date) |>
  mutate(avg_value = mean(value, na.rm = T)) |>
  distinct(state_fips, week_date, avg_value) |>
  dplyr::filter(week_date > as.Date('2016-09-01') & week_date < as.Date('2019-08-31'))

# Merge on external data
flu_state <- left_join(flu_state, ili_state[, c(2, 3, 4)], by = c('year_week_dt' = 'week_date', 'state_fips'))
flu_state <- left_join(flu_state, nrevss_state[, c(2, 3, 4)], by = c('year_week_dt' = 'week_date', 'state_fips'))
flu_state <- left_join(flu_state, flu_surv_state_avg[, c(1, 2, 3)], by = c('year_week_dt' = 'week_date', 'state_fips'))


ggplot(flu_state) + geom_line(aes(x = year_week_dt, y = pred_sum), color = 'blue') +
  facet_wrap(~state_fips, scale = 'free')

ggplot(flu_state) + 
  geom_line(aes(x = year_week_dt, y = value.x), color = 'red') +
  facet_wrap(~state_fips, scale = 'free')

# CORRELATIONS # 
flu_corr <- flu_state |>
  group_by(state_fips) |>
  mutate(ili_corr = cor(pred_sum, value.x, use = 'pairwise.complete.obs'),
         nrevss_corr = cor(pred_sum, value.y, use = 'pairwise.complete.obs'),
         surv_corr = cor(pred_sum, avg_value, use = 'pairwise.complete.obs')) |>
  distinct(state_fips, ili_corr, nrevss_corr, surv_corr)

# Load map
usa_albers_st <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_state.rds'))   # convert to sf

# Merge on data
usa_albers_st_corr <- left_join(usa_albers_st , flu_corr, by = c('STATEFP' = 'state_fips'))

map_theme <- theme(legend.position = "bottom",
                   axis.title = element_text(size=16),
                   strip.background = element_blank(),
                   plot.title = element_text(size=18, hjust = 0.5),
                   panel.border = element_rect(fill=NA, linewidth = 0.8, color = 'white'),
                   legend.title = element_text(size = 16),
                   legend.text = element_text(size = 12),
                   legend.key.height = unit(0.38, 'cm'),
                   legend.key.width = unit(0.7, "cm"))

corr_ili <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_st_corr), aes(fill = ili_corr, group = STATEFP), color= 'black', linewidth = 0.1) +
  scale_fill_gradient('Correlation   \n', low = "white", high = '#347DC1', na.value = "grey",
                      breaks = c(0.25, 0.50, 0.75, 1), 
                      limits = c(0.25, 1.00)) + 
  ggtitle('ILI-NET\n') +
  theme_void() + map_theme 

corr_ili


corr_nrevss <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_st_corr), aes(fill = nrevss_corr, group = STATEFP), color= 'black', linewidth = 0.1) +
  scale_fill_gradient('Correlation   \n', low = "white", high = '#F8766D', na.value = "grey",
                      breaks = c(0.25, 0.50, 0.75, 1), 
                      limits = c(0.25, 1.00)) + 
  ggtitle('NREVSS\n') +
  theme_void() + map_theme 

corr_nrevss

corr_surv <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_st_corr), aes(fill = surv_corr, group = STATEFP), color= 'black', linewidth = 0.1) +
  scale_fill_gradient('Correlation   \n', low = "white", high = '#9B59B6', na.value = "grey",
                      breaks = c(0.25, 0.50, 0.75, 1), 
                      limits = c(0.25, 1.00))+ 
  ggtitle('FluSurv-NET\n') +
  theme_void() + map_theme 

corr_surv

nys_county$county_fips <- as.character(nys_county$county_fips )

flu_county <- rbind(flu_county_2016, flu_county_2017, flu_county_2018)

nys_county <- left_join(nys_county, flu_county[, c(1, 2, 3)], by = c('week_date' = 'year_week_dt',
                                                                     'county_fips'))

nys_corr <- nys_county |> group_by(county_fips) |>
  dplyr::filter(county_fips != '36041') |>
  mutate(nys_corr = cor(pred_sum, value, use = 'pairwise.complete.obs')) |>
  distinct(county_fips, nys_corr)


county <- read_sf(dsn = './raw/cb_2021_us_county_20m/', 
                  layer = 'cb_2021_us_county_20m')

county_nys <- county |> dplyr::filter(STATEFP == '36')

county_tomp <- county_nys |> dplyr::filter(COUNTYFP == '109' | COUNTYFP == '107')

# Merge on prev data
county_nys_corr <- left_join(county_nys, nys_corr, by = c('GEOID' = 'county_fips'))

ggplot() +
  geom_sf(data = st_as_sf(county_tomp), fill = 'white', color= 'black', linewidth = 1)

corr_map <- ggplot() +
  geom_sf(data = st_as_sf(county_nys_corr ), aes(fill = nys_corr, group = GEOID), color= 'black', linewidth = 0.1) +
  scale_fill_gradient('Correlation   \n', low = "white", high = "#e6a532", na.value = "grey",
                      breaks = c(0.25, 0.50, 0.75, 1), 
                      limits = c(0.25, 1.00)) + 
  ggtitle('NYS Lab-Confirmed') +
  theme_void() + map_theme 



figure_3_corr <- cowplot::plot_grid(corr_ili, corr_surv, corr_nrevss, corr_map,
                                    nrow = 2,
                                    labels = c("f"),
                                    label_size = 20)

figure_3_corr


# PEAK TIMES #


flu_surv_state_avg <- flu_surv_state |> group_by(state_fips, week_date) |>
  mutate(avg_value = mean(value_roll, na.rm = T)) |>
  distinct(state_fips, week_date, avg_value) |>
  dplyr::filter(week_date > as.Date('2016-09-01') & week_date < as.Date('2019-08-30'))


flu_surv_state_avg$month <- month(flu_surv_state_avg$week_date)
flu_surv_state_avg$year <- year(flu_surv_state_avg$week_date)
flu_surv_state_avg$Season <- ifelse(flu_surv_state_avg$year == 2016, '2016-2017', '')
flu_surv_state_avg$Season <- ifelse(flu_surv_state_avg$year == 2017 & flu_surv_state_avg$month < 9, '2016-2017', flu_surv_state_avg$Season)
flu_surv_state_avg$Season <- ifelse(flu_surv_state_avg$year == 2017 & flu_surv_state_avg$month > 8, '2017-2018', flu_surv_state_avg$Season)
flu_surv_state_avg$Season <- ifelse(flu_surv_state_avg$year == 2018 & flu_surv_state_avg$month < 9, '2017-2018', flu_surv_state_avg$Season)
flu_surv_state_avg$Season <- ifelse(flu_surv_state_avg$year == 2018 & flu_surv_state_avg$month > 8, '2018-2019', flu_surv_state_avg$Season)
flu_surv_state_avg$Season <- ifelse(flu_surv_state_avg$year == 2019 & flu_surv_state_avg$month < 9, '2018-2019', flu_surv_state_avg$Season)


state_peaks_surv <- flu_surv_state_avg |> 
  ungroup() |>
  dplyr::filter(!is.na(avg_value) & !is.na(state_fips)) |>
  # Calculate the row number
  mutate(index = row_number()) |>
  group_by(state_fips, Season) |> 
  # Calculate multiple peaks greater than 10 cases and larger than 10 neighbors
  mutate(peaks_obs = seq_along(index) %in% 
           findpeaks(avg_value, MinPeakDistance = 10, MinPeakHeight = 1)$loc) |>
  # Remove first and last observations, which cannot be peaks but might be max
  dplyr::filter(row_number() != 1 & row_number() != n()) |>
  # Calculate the maximum value among the peaks
  mutate(is_max_obs = avg_value == max(avg_value)) |>
  select(c(state_fips, Season, week_date, avg_value, peaks_obs, is_max_obs))

state_peaks_obs_surv <- state_peaks_surv |> 
  # Select the peaks
  dplyr::filter(peaks_obs == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_obs == T) |>
  select(c(state_fips, Season, week_date, avg_value, peaks_obs, is_max_obs)) |>
  rename('week_surv' = 'week_date') |>
  group_by(state_fips, Season) |>
  # Break any ties by taking the earliest peak
  arrange(week_surv, Season) |>
  slice_head(n = 1) |>
  mutate(week_surv_num = week(week_surv)) |>
  dplyr::select(c(state_fips, Season, week_surv, week_surv_num))

ili_state$month <- month(ili_state$week_date)
ili_state$year <- year(ili_state$week_date)
ili_state$Season <- ifelse(ili_state$year == 2016, '2016-2017', '')
ili_state$Season <- ifelse(ili_state$year == 2017 & ili_state$month < 9, '2016-2017', ili_state$Season)
ili_state$Season <- ifelse(ili_state$year == 2017 & ili_state$month > 8, '2017-2018', ili_state$Season)
ili_state$Season <- ifelse(ili_state$year == 2018 & ili_state$month < 9, '2017-2018', ili_state$Season)
ili_state$Season <- ifelse(ili_state$year == 2018 & ili_state$month > 8, '2018-2019', ili_state$Season)
ili_state$Season <- ifelse(ili_state$year == 2019 & ili_state$month < 9, '2018-2019', ili_state$Season)

state_peaks_ili <- ili_state |> 
  ungroup() |>
  dplyr::filter(!is.na(value_roll) & !is.na(state_fips)) |>
  mutate(value_roll = ifelse(value_roll < 0, 0, value_roll)) |>
  dplyr::filter(week_date > as.Date('2016-09-01') & week_date < as.Date('2019-08-01')) |>
  # Calculate the row number
  mutate(index = row_number()) |>
  group_by(state_fips, Season) |> 
  # Calculate multiple peaks greater than 10 cases and larger than 10 neighbors
  mutate(peaks_obs = seq_along(index) %in% 
           findpeaks(value_roll, MinPeakDistance = 10, MinPeakHeight = 1)$loc) |>
  # Remove first and last observations, which cannot be peaks but might be max
  dplyr::filter(row_number() != 1 & row_number() != n()) |>
  # Calculate the maximum value among the peaks
  mutate(is_max_obs = value_roll == max(value_roll)) |>
  select(c(state_fips, Season, week_date, value_roll, peaks_obs, is_max_obs))

state_peaks_obs_ili <- state_peaks_ili  |> 
  # Select the peaks
  dplyr::filter(peaks_obs == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_obs == T) |>
  select(c(state_fips, Season, week_date, value_roll, peaks_obs, is_max_obs)) |>
  rename('week_ili' = 'week_date') |>
  group_by(state_fips, Season) |>
  # Break any ties by taking the earliest peak
  arrange(week_ili) |>
  slice_head(n = 1) |>
  mutate(week_ili_num = week(week_ili)) |>
  dplyr::select(c(state_fips, Season, week_ili, week_ili_num))

nrevss_state$month <- month(nrevss_state$week_date)
nrevss_state$year <- year(nrevss_state$week_date)
nrevss_state$Season <- ifelse(nrevss_state$year == 2016, '2016-2017', '')
nrevss_state$Season <- ifelse(nrevss_state$year == 2017 & nrevss_state$month < 9, '2016-2017', nrevss_state$Season)
nrevss_state$Season <- ifelse(nrevss_state$year == 2017 & nrevss_state$month > 8, '2017-2018', nrevss_state$Season)
nrevss_state$Season <- ifelse(nrevss_state$year == 2018 & nrevss_state$month < 9, '2017-2018', nrevss_state$Season)
nrevss_state$Season <- ifelse(nrevss_state$year == 2018 & nrevss_state$month > 8, '2018-2019', nrevss_state$Season)
nrevss_state$Season <- ifelse(nrevss_state$year == 2019 & nrevss_state$month < 9, '2018-2019', nrevss_state$Season)


state_peaks_nrevss <- nrevss_state |> 
  ungroup() |>
  dplyr::filter(!is.na(value_roll) & !is.na(state_fips)) |>
  mutate(value_roll = ifelse(value_roll < 0, 0, value_roll)) |>
  mutate(value_roll = ifelse(is.na(value_roll), 0, value_roll)) |>
  dplyr::filter(week_date > as.Date('2016-09-01') & week_date < as.Date('2019-08-01')) |>
  dplyr::filter(state_fips != '02') |>
  dplyr::filter(state_fips != '38') |>
  # Calculate the row number
  mutate(index = row_number()) |>
  group_by(state_fips, Season) |> 
  # Calculate multiple peaks greater than 10 cases and larger than 10 neighbors
  mutate(peaks_obs = seq_along(index) %in% 
           findpeaks(value_roll, MinPeakDistance = 10, MinPeakHeight = 1)$loc) |>
  # Remove first and last observations, which cannot be peaks but might be max
  dplyr::filter(row_number() != 1 & row_number() != n()) |>
  # Calculate the maximum value among the peaks
  mutate(is_max_obs = value_roll == max(value_roll)) |>
  select(c(state_fips, Season, week_date, value_roll, peaks_obs, is_max_obs))

state_peaks_obs_nrevss <- state_peaks_nrevss  |> 
  # Select the peaks
  dplyr::filter(peaks_obs == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_obs == T) |>
  select(c(state_fips, Season, week_date, value_roll, peaks_obs, is_max_obs)) |>
  rename('week_nrevss' = 'week_date') |>
  group_by(state_fips, Season) |>
  # Break any ties by taking the earliest peak
  arrange(week_nrevss) |>
  slice_head(n = 1) |>
  mutate(week_nrevss_num = week(week_nrevss)) |>
  dplyr::select(c(state_fips, Season, week_nrevss, week_nrevss_num))


nys_county$month <- month(nys_county$week_date)
nys_county$year <- year(nys_county$week_date)
nys_county$Season <- ifelse(nys_county$year == 2016, '2016-2017', '')
nys_county$Season <- ifelse(nys_county$year == 2017 & nys_county$month < 9, '2016-2017', nys_county$Season)
nys_county$Season <- ifelse(nys_county$year == 2017 & nys_county$month > 8, '2017-2018', nys_county$Season)
nys_county$Season <- ifelse(nys_county$year == 2018 & nys_county$month < 9, '2017-2018', nys_county$Season)
nys_county$Season <- ifelse(nys_county$year == 2018 & nys_county$month > 8, '2018-2019', nys_county$Season)
nys_county$Season <- ifelse(nys_county$year == 2019 & nys_county$month < 9, '2018-2019', nys_county$Season)


county_peaks_nys <- nys_county |> 
  ungroup() |>
  dplyr::filter(!is.na(value_roll) & !is.na(county_fips)) |>
  mutate(value_roll = ifelse(value_roll < 0, 0, value_roll)) |>
  dplyr::filter(week_date > as.Date('2016-09-01') & week_date < as.Date('2019-08-01')) |>
  # Calculate the row number
  mutate(index = row_number()) |>
  group_by(county_fips, Season) |> 
  # Calculate multiple peaks greater than 10 cases and larger than 10 neighbors
  mutate(peaks_obs = seq_along(index) %in% 
           findpeaks(value_roll, MinPeakDistance = 10, MinPeakHeight = 1)$loc) |>
  # Remove first and last observations, which cannot be peaks but might be max
  dplyr::filter(row_number() != 1 & row_number() != n()) |>
  # Calculate the maximum value among the peaks
  mutate(is_max_obs = value_roll == max(value_roll)) |>
  select(c(county_fips, Season, week_date, value_roll, peaks_obs, is_max_obs))

county_peaks_obs_nys <- county_peaks_nys  |> 
  # Select the peaks
  dplyr::filter(peaks_obs == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_obs == T) |>
  select(c(county_fips, Season, week_date, value_roll, peaks_obs, is_max_obs)) |>
  rename('week_nys' = 'week_date') |>
  group_by(county_fips, Season) |>
  # Break any ties by taking the earliest peak
  arrange(week_nys) |>
  slice_head(n = 1) |>
  mutate(week_nys_num = week(week_nys)) |>
  dplyr::select(c(county_fips, Season, week_nys, week_nys_num))


state_peaks_pred <- rbind(state_peaks_pred_2018, state_peaks_pred_2017, state_peaks_pred_2016)

state_peaks_pred$month <- month(state_peaks_pred$week_pred)
state_peaks_pred$year <- year(state_peaks_pred$week_pred)
state_peaks_pred$Season <- ifelse(state_peaks_pred$year == 2016, '2016-2017', '')
state_peaks_pred$Season <- ifelse(state_peaks_pred$year == 2017 & state_peaks_pred$month < 9, '2016-2017', state_peaks_pred$Season)
state_peaks_pred$Season <- ifelse(state_peaks_pred$year == 2017 & state_peaks_pred$month > 8, '2017-2018', state_peaks_pred$Season)
state_peaks_pred$Season <- ifelse(state_peaks_pred$year == 2018 & state_peaks_pred$month < 9, '2017-2018', state_peaks_pred$Season)
state_peaks_pred$Season <- ifelse(state_peaks_pred$year == 2018 & state_peaks_pred$month > 8, '2018-2019', state_peaks_pred$Season)
state_peaks_pred$Season <- ifelse(state_peaks_pred$year == 2019 & state_peaks_pred$month < 9, '2018-2019', state_peaks_pred$Season)


state_peaks_corr <- left_join(state_peaks_pred, state_peaks_obs_ili, by = c('state_fips', 'Season'))

state_peaks_corr <- left_join(state_peaks_corr, state_peaks_obs_nrevss, by = c('state_fips', 'Season'))

state_peaks_corr <- left_join(state_peaks_corr, state_peaks_obs_surv, by = c('state_fips', 'Season'))

state_peaks_corr$diff <- as.numeric(state_peaks_corr$week_pred - state_peaks_corr$week_nrevss)
state_peaks_corr_filt <- state_peaks_corr |> dplyr::filter(!is.na(diff))
sqrt(sum((state_peaks_corr_filt$diff)^2, na.rm = T) / length(state_peaks_corr_filt$diff)) 
state_peaks_corr_filt <- state_peaks_corr_filt |>
  mutate(week_nrevss_num_adj = ifelse(week_nrevss_num > 30, week_nrevss_num - 52, week_nrevss_num),
         week_pred_num_adj = ifelse(week_pred_num > 30, week_pred_num - 52, week_pred_num))
cor(state_peaks_corr_filt$week_pred_num_adj, state_peaks_corr_filt$week_nrevss_num_adj, 
    use = 'complete.obs', method = 'spearman')


state_peaks_corr_ili <- state_peaks_corr |>
  mutate(count = 1) |>
  group_by(week_pred, week_ili) |>
  mutate(count_sum = sum(count, na.rm = T)) |>
  distinct(week_pred, week_ili, count_sum)



forrest_theme <- theme(legend.position =  "bottom",
                       legend.position.inside = c(0.873, 0.237),
                       axis.text = element_text(size=14, color = 'black'),
                       axis.title = element_text(size=16),
                       strip.background = element_blank(),
                       legend.title = element_text(size=16),
                       legend.text = element_text(size=12),
                       plot.title = element_text(size=18, hjust = 0.5),
                       legend.box="horizontal",
                       legend.box.background = element_rect(colour = "white", 
                                                            fill = 'white'),
                       legend.key.height = unit(0.38, 'cm'),
                       legend.key.width = unit(0.7, "cm")) 


state_peaks_corr_ili$count_sum <- ifelse(state_peaks_corr_ili$count_sum == 9, 8, state_peaks_corr_ili$count_sum)
corr_ili_scatter <- ggplot(state_peaks_corr_ili[!is.na(state_peaks_corr_ili$week_ili),], aes(as.Date(week_ili), as.Date(week_pred), fill= count_sum)) + 
  geom_tile(color = 'black', linewidth = 0.25) + 
  scale_fill_gradient2('Number of \nStates', low = 'white', high = '#347DC1') +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black",
              linewidth = 1) +
  theme_minimal() + 
  labs(title = 'ILI-NET',
       x = "ILI-NET Peak",
       y = "Syndromic Peak") + forrest_theme +
  scale_x_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2016-12-01'), as.Date('2017-05-27'))) +
  scale_y_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2016-12-01'), as.Date('2017-05-27')))
corr_ili_scatter


state_peaks_corr_surv <- state_peaks_corr |>
  mutate(count = 1) |>
  group_by(week_pred, week_surv) |>
  mutate(count_sum = sum(count, na.rm = T)) |>
  distinct(week_pred, week_surv, count_sum)

state_peaks_corr_surv[!is.na(state_peaks_corr_surv$week_surv),]

corr_surv_scatter <- ggplot(state_peaks_corr_surv[!is.na(state_peaks_corr_surv$week_surv),], aes(as.Date(week_surv), as.Date(week_pred), fill= count_sum)) + 
  geom_tile(color = 'black', linewidth = 0.25) + 
  scale_fill_gradient2('Number of \nStates', low = 'white', high = '#9B59B6') +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black",
              linewidth = 1) +
  theme_minimal() + 
  labs(title = 'FluSurv-NET',
       x = "FluSurv-NET Peak",
       y = "Syndromic Peak") + forrest_theme +
  scale_x_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2016-12-01'), as.Date('2017-05-27'))) +
  scale_y_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2016-12-01'), as.Date('2017-05-27')))
corr_surv_scatter

state_peaks_corr_nrevss <- state_peaks_corr |>
  mutate(count = 1) |>
  group_by(week_pred, week_nrevss) |>
  mutate(count_sum = sum(count, na.rm = T)) |>
  distinct(week_pred, week_nrevss, count_sum)


corr_nrevss_scatter <- ggplot(state_peaks_corr_nrevss[!is.na(state_peaks_corr_nrevss$week_nrevss),], aes(as.Date(week_nrevss), as.Date(week_pred), fill= count_sum)) + 
  geom_tile(color = 'black', linewidth = 0.25) + 
  scale_fill_gradient2('Number of \nStates', low = 'white', high = '#F8766D') +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black",
              linewidth = 1) +
  theme_minimal() + 
  labs(title = 'NREVSS',
       x = "NREVSS Peak",
       y = "Syndromic Peak") + forrest_theme +
  scale_x_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2016-12-01'), as.Date('2017-05-27'))) +
  scale_y_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2016-12-01'), as.Date('2017-05-27')))
corr_nrevss_scatter


county_peaks_pred <- rbind(county_peaks_pred_2018, county_peaks_pred_2017, county_peaks_pred_2016)

county_peaks_pred$month <- month(county_peaks_pred$week_pred)
county_peaks_pred$year <- year(county_peaks_pred$week_pred)
county_peaks_pred$Season <- ifelse(county_peaks_pred$year == 2016, '2016-2017', '')
county_peaks_pred$Season <- ifelse(county_peaks_pred$year == 2017 & county_peaks_pred$month < 9, '2016-2017', county_peaks_pred$Season)
county_peaks_pred$Season <- ifelse(county_peaks_pred$year == 2017 & county_peaks_pred$month > 8, '2017-2018', county_peaks_pred$Season)
county_peaks_pred$Season <- ifelse(county_peaks_pred$year == 2018 & county_peaks_pred$month < 9, '2017-2018', county_peaks_pred$Season)
county_peaks_pred$Season <- ifelse(county_peaks_pred$year == 2018 & county_peaks_pred$month > 8, '2018-2019', county_peaks_pred$Season)
county_peaks_pred$Season <- ifelse(county_peaks_pred$year == 2019 & county_peaks_pred$month < 9, '2018-2019', county_peaks_pred$Season)



county_peaks_corr_nys <- left_join(county_peaks_obs_nys, county_peaks_pred, by = c('county_fips', 'Season'))

county_peaks_corr_nys$diff <- as.numeric(county_peaks_corr_nys$week_pred - county_peaks_corr_nys$week_nys)
county_peaks_corr_nys_filt <- county_peaks_corr_nys |> dplyr::filter(!is.na(diff))
sqrt(sum((county_peaks_corr_nys$diff)^2, na.rm = T) / length(county_peaks_corr_nys$diff)) 
county_peaks_corr_nys <- county_peaks_corr_nys |>
  mutate(week_nys_num_adj = ifelse(week_nys_num > 30, week_nys_num - 52, week_nys_num),
         week_pred_num_adj = ifelse(week_pred_num > 30, week_pred_num - 52, week_pred_num))
cor(county_peaks_corr_nys$week_pred_num_adj, county_peaks_corr_nys$week_nys_num_adj, 
    use = 'complete.obs', method = 'spearman')



county_peaks_corr_nys_sum <- county_peaks_corr_nys |>
  mutate(count = 1) |>
  group_by(week_pred, week_nys) |>
  mutate(count_sum = sum(count, na.rm = T)) |>
  distinct(week_pred, week_nys, count_sum)

county_peaks_corr_nys_sum <- county_peaks_corr_nys_sum |>
  dplyr::filter(week_pred < as.Date('2017-05-27'))

corr_nys_scatter <- ggplot(county_peaks_corr_nys_sum[!is.na(county_peaks_corr_nys_sum$week_pred),], aes(as.Date(week_nys), as.Date(week_pred), fill= count_sum)) + 
  geom_tile(color = 'black', linewidth = 0.25) + 
  scale_fill_gradient2('Number of \nCounties', low = 'white', high = '#e6a532') +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black",
              linewidth = 1) +
  theme_minimal() + 
  labs(title = 'NYS Lab-Confirmed',
       x = "NYS Peak",
       y = "Syndromic Peak") + forrest_theme +
  scale_x_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2016-12-01'), as.Date('2017-05-27'))) +
  scale_y_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2016-12-01'), as.Date('2017-05-27')))
corr_nys_scatter





figure_3_scatter <- cowplot::plot_grid(corr_ili_scatter, corr_surv_scatter , corr_nrevss_scatter, corr_nys_scatter,
                                       nrow = 2,
                                       labels = c("g"),
                                       label_size = 20)

figure_3_scatter



flu_national <- flu_county |> group_by(year_week_dt) |>
  mutate(pred_sum = sum(pred_sum, na.rm = T)) |>
  distinct(year_week_dt, pred_sum) |>
  dplyr::filter(year_week_dt < as.Date('2019-09-01')) |>
  dplyr::filter(year_week_dt > as.Date('2016-08-31')) |>
  ungroup()

# Add month and year variables
flu_national$month <- month(flu_national$year_week_dt)
flu_national$year <- year(flu_national$year_week_dt)

# Load all cause data
all_cause <- read.csv('./raw/county_week_ac_v2.csv', header = T)

# Impute all cause counts
all_cause$patient_count_imp <- all_cause$all_cause
all_cause$patient_count_imp[all_cause$patient_count_imp == '<=5'] <- 
  sample(1:5, length(all_cause$patient_count_imp[all_cause$patient_count_imp == '<=5']), 
         replace= TRUE)
# Convert imputed count to numeric
all_cause$patient_count_imp <- as.numeric(all_cause$patient_count_imp)

# Convert year week to week date
all_cause <- all_cause |>
  mutate(week_date = sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", year_week)) |>
  mutate(week_date = ISOweek2date(week_date))

# Collapse to national level
all_cause_national <- all_cause |> group_by(week_date) |>
  mutate(all_cause_sum = sum(patient_count_imp)) |>
  distinct(week_date, all_cause_sum) |>
  ungroup() |>
  dplyr::filter(week_date < as.Date('2020-03-01')) |>
  dplyr::filter(week_date > as.Date('2016-08-31')) |>
  mutate(all_cause_avg = mean(all_cause_sum))

# Quickly visualize all cause data
ggplot(all_cause_national) + 
  geom_line(aes(x = week_date, y = all_cause_sum / all_cause_avg), color = 'blue')


# Calculate summer averages
flu_national_summer <- flu_national |> 
  dplyr::filter(month < 9 & month > 5) |> group_by(year) |>
  dplyr::filter(year > 2016) |>
  mutate(pred_summ_mean = mean(pred_sum)) |>
  distinct(year, pred_summ_mean)

# Merge on summer data to the correct season (summer 2016 should align with the 2016-2017 flu season)
flu_national_summer[nrow(flu_national_summer) + 1, ] <- flu_national_summer[1, ] 
flu_national_summer$year[4] <- 2016
flu_national$year_2 <- ifelse(flu_national$month < 6, flu_national$year - 1, flu_national$year)
flu_national <- left_join(flu_national, flu_national_summer, by = c('year_2' = 'year'))

# Calculate rolling average
flu_national$roll_pred <- rollmean(flu_national$pred_sum, 3, fill = NA, align = 'center')

# Calculate Z score
flu_national$pred_z_roll <- (flu_national$roll_pred - flu_national$pred_summ_mean) / sd(flu_national$roll_pred, na.rm = T)
flu_national$pred_z <- (flu_national$pred_sum - flu_national$pred_summ_mean) / sd(flu_national$pred_sum, na.rm = T)


# Create Season Variable
flu_national$Season <- ifelse(flu_national$year == 2016, '2016-2017', '')
flu_national$Season <- ifelse(flu_national$year == 2017 & flu_national$month < 9, '2016-2017', flu_national$Season)
flu_national$Season <- ifelse(flu_national$year == 2017 & flu_national$month > 8, '2017-2018', flu_national$Season)
flu_national$Season <- ifelse(flu_national$year == 2018 & flu_national$month < 9, '2017-2018', flu_national$Season)
flu_national$Season <- ifelse(flu_national$year == 2018 & flu_national$month > 8, '2018-2019', flu_national$Season)
flu_national$Season <- ifelse(flu_national$year == 2019 & flu_national$month < 9, '2018-2019', flu_national$Season)
flu_national$Season <- ifelse(flu_national$year == 2019 & flu_national$month > 8, '2019-2020', flu_national$Season)
flu_national$Season <- ifelse(flu_national$year == 2020, '2019-2020', flu_national$Season)

flu_national$pred_above <- ifelse(flu_national$pred_z > 0.1, 1, 0)

flu_national <- left_join(flu_national, ili_nat[, c(3, 6)], by = c('year_week_dt' = 'week_date'))
flu_national <- flu_national |>
  rename('ili_z' = 'z_score')

flu_national <- left_join(flu_national, nrevss_nat[, c(3, 6)], by = c('year_week_dt' = 'week_date'))
flu_national <- flu_national |>
  rename('nrevss_z' = 'z_score')

flu_national <- left_join(flu_national, flu_surv_nat[, c(2, 6)], by = c('year_week_dt' = 'week_date'))
flu_national <- flu_national |>
  rename('surv_z' = 'z_score')


flu_national$ili_above <- ifelse(flu_national$ili_z > 0.1, 1, 0)
flu_national$nrevss_above <- ifelse(flu_national$nrevss_z > 0.1, 1, 0)
flu_national$surv_above <- ifelse(flu_national$surv_z > 0.1, 1, 0)

flu_national_duration <- flu_national |>
  group_by(Season) |>
  mutate(pred_duration = sum(pred_above),
         ili_duration = sum(ili_above),
         nrevss_duration = sum(nrevss_above),
         surv_duration = sum(surv_above, na.rm = T)) |>
  distinct(Season, pred_duration, ili_duration, 
           nrevss_duration, surv_duration)

flu_national_duration_long <- flu_national_duration |>
  pivot_longer(
    cols = pred_duration:surv_duration, 
    names_to = "Source",
    values_to = "value"
  ) |>
  dplyr::filter(Season != '2019-2020')


flu_national_duration_long$Source[flu_national_duration_long$Source == 'pred_duration'] <- 'Predicted'
flu_national_duration_long$Source[flu_national_duration_long$Source == 'obs_duration'] <- 'Claims'
flu_national_duration_long$Source[flu_national_duration_long$Source == 'ili_duration'] <- 'ILI-NET'
flu_national_duration_long$Source[flu_national_duration_long$Source == 'nrevss_duration'] <- 'NREVSS'
flu_national_duration_long$Source[flu_national_duration_long$Source == 'surv_duration'] <- 'FluSurv-NET'

flu_national_duration_long_col <- flu_national_duration_long |>
  group_by(Source) |>
  mutate(min = min(value),
         max = max(value)) |>
  distinct(Source, min, max) |>
  dplyr::filter(Source != 'Claims')

Lines <- c('Claims' = 'solid', 'FluSurv-NET' = 'solid', 
           'ILI-NET' = 'solid', 'NREVSS' = 'solid', 'Predicted' = 'dashed')

Colors <- c('Claims' = '#bf6fa7', 'FluSurv-NET' = '#9B59B6', 
            'ILI-NET' = '#347DC1', 'NREVSS' = '#F8766D', 'Predicted' = 'gray39')

percent_theme <- theme(legend.position = "none",
                       axis.text = element_text(size=14),
                       axis.text.x  = element_text(size=14, angle = 45, hjust = 1),
                       axis.title = element_text(size=16),
                       strip.background = element_blank(),
                       plot.title = element_text(size=18),
                       legend.title = element_text(size = 16),
                       legend.text = element_text(size = 14))

flu_national_duration_plot <- ggplot(flu_national_duration_long_col, 
                                     aes(x = Source, group = Source,
                                         ymin = min, ymax = max)) +
  geom_errorbar(linewidth = 2, width = 0.2, aes(color = Source)) +
  scale_color_manual(values = Colors) + ylab('Number of Weeks') +
  labs(title = "Season Duration by Source") + theme_minimal() + 
  percent_theme

flu_national_duration_plot

figure_4 <- cowplot::plot_grid(flu_national_duration_plot
                               , figure_3_corr, figure_3_scatter,
                               rel_widths = c(0.5, 1, 1),
                               nrow = 1,
                               labels = c("e", "f", "g"),
                               label_size = 20)

ggsave('./figs/figure_4.jpg', plot = figure_4, height = 8, width = 20)

saveRDS(measure_plot, './figs/measure_plot.RData')
saveRDS(peak_comp_plot_2018, './figs/peak_plot_2018.RData')


figure_3 <- cowplot::plot_grid(peak_comp_plot_2016, peak_comp_plot_2017, peak_comp_plot_2018, measure_plot, 
                               nrow = 1,
                               labels = c("a", "b", "c", "d"),
                               label_size = 20)

ggsave('./figs/figure_3.jpg', plot = figure_3, height = 6, width = 20)




figure_3.5 <- cowplot::plot_grid(figure_3, figure_4, 
                               nrow = 2,
                               labels = c("", ""),
                               label_size = 20)

ggsave('./figs/figure_3.5.jpg', plot = figure_3.5, height = 14, width = 20)





# Collapse data for graphing purposes
merge_onset_col <- merge_onset |>
  mutate(count = 1) |>
  group_by(onset_week_obs_adj, onset_week_pred_adj) |>
  mutate(count_sum = sum(count)) |>
  distinct(onset_week_obs_adj, onset_week_pred_adj, count_sum)


flu_state$Source <- 'Predicted'

duration_state <- rbind(flu_state[, c('state_fips', 'week_date', 'z_score', 'Source')],
                        ili_state[, c('state_fips', 'week_date', 'z_score', 'Source')],
                        nrevss_state[, c('state_fips', 'week_date', 'z_score', 'Source')],
                        flu_surv_state[, c('state_fips', 'week_date', 'z_score', 'Source')])

duration_dat <- duration_state |>
  dplyr::filter(week_date < as.Date('2017-09-01')) |>
  group_by(state_fips, Source) |>
  mutate(count = ifelse(z_score > 0.05, 1, 0),
         duration = sum(count, na.rm = T)) |>
  distinct(state_fips, Source, duration) |>
  dplyr::filter(duration > 0) 
  


ggplot(duration_dat) + geom_violin(aes(x = Source, y = duration, fill = Source))


ggplot(flu_state) + 
  geom_line(aes(x = week_date, y = pred_z), color = 'red') +
  facet_wrap(~state_fips, scale = 'free')













