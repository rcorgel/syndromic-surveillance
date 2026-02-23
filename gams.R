################################################################################
# File Name: figure_2                                                          #
#                                                                              #
# Purpose:   Create figure 2 for the manuscript.                               #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Create figure 2                                                #
#               a. Create profile plots                                        #
#               b. Define candidate models                                     #
#               c. Create age forest plots                                     #
#               d. Create time forest plots                                    #
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
library(lme4)
library(car)
library(mgcv)
library(gsignal)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

######################
# 2. CREATE FIGURE 2 #
######################

flu_2016 <- readRDS('./tmp/flu_filt_symp_p2_2016.rds')

test_2016 <- flu_2016 |> group_by(year_week_dt, county_fips, state_fips) |>
  mutate(flu = sum(flu * patient_count),
         fever = sum(fever * patient_count),
         myalgia = sum(myalgia * patient_count),
         cough = sum(cough * patient_count),
         sore_throat = sum(sore_throat * patient_count), 
         short_breath = sum(short_breath * patient_count), 
         hypoxemia = sum(hypoxemia * patient_count), 
         chest_pain = sum(chest_pain * patient_count), 
         bronchitis = sum(bronchitis * patient_count),
         nausea_vom = sum(nausea_vom * patient_count), 
         diarrhea = sum(diarrhea * patient_count), 
         fatigue = sum(fatigue * patient_count), 
         headache = sum(headache * patient_count),
         congestion = sum(congestion * patient_count), 
         sneezing = sum(sneezing * patient_count)) |>
  distinct(year_week_dt, county_fips, state_fips, flu, fever, myalgia, cough, sore_throat, 
           short_breath, hypoxemia, chest_pain, bronchitis,
           nausea_vom, diarrhea, fatigue, headache,
           congestion, sneezing)

test_2016$county_fips <- as.factor(test_2016$county_fips)

all_cause_season <- read_csv('raw/county_season_ac_v3_imputed.csv')
all_cause_season <- all_cause_season |> dplyr::filter(season == '2016-2017')

test_2016 <- left_join(test_2016, all_cause_season, by = c('county_fips'))

test_counties <- test_2016 |>
  mutate(flu = flu / all_cause_wtd,
         fever = fever / all_cause_wtd,
         myalgia = myalgia / all_cause_wtd,
         cough = cough / all_cause_wtd,
         sore_throat = sore_throat / all_cause_wtd, 
         short_breath = short_breath / all_cause_wtd, 
         hypoxemia = hypoxemia / all_cause_wtd, 
         chest_pain = chest_pain / all_cause_wtd, 
         bronchitis = bronchitis / all_cause_wtd,
         nausea_vom = nausea_vom / all_cause_wtd, 
         diarrhea = diarrhea / all_cause_wtd, 
         fatigue = fatigue / all_cause_wtd, 
         headache = headache / all_cause_wtd,
         congestion = congestion / all_cause_wtd, 
         sneezing = sneezing / all_cause_wtd)

# Epi week
test_2016$epi_week <- epiweek(test_2016$year_week_dt)

# Change reference week to week 26
test_2016$epi_week_adj <- ((test_2016$epi_week - 26) %% 52) + 1

model <- gam(flu ~ 
               s(fever, k = 3) +
               s(cough, k = 3) +
               s(sore_throat, k = 3) +
               s(congestion, k = 3) +
               s(bronchitis, k = 3) +
               s(chest_pain, k = 3) +
               s(nausea_vom, k = 3) +
               s(fatigue, k = 3) +
               s(diarrhea, k = 3) +
               s(sneezing, k = 3) +
               s(hypoxemia, k = 3) +
               s(short_breath, k = 3) +
               s(myalgia, k = 3) +
               s(headache, k = 3),
             family = quasipoisson(link = "log"),
             data = test_2016, 
             method = "REML")

test_2016$pred <- predict(model, newdata = test_2016, type = "response") 

library(zoo)
flu_county_peak <- test_2016 |> arrange(year_week_dt) |> group_by(county_fips) |> 
  mutate(flu_sum_roll = rollmean(flu, k = 4, align = 'right', na.pad = TRUE),
         flu_pred_roll = rollmean(pred, k = 4, align = 'right', na.pad = TRUE))

# Calculate peaks
county_peaks_usa <- flu_county_peak |> dplyr::filter(county_fips != '15005') |>
  ungroup() |>
  dplyr::filter(!is.na(flu_sum_roll)) |>
  dplyr::filter(!is.na(flu_pred_roll)) |>
  # Calculate the row number
  mutate(index = row_number()) |>
  group_by(county_fips) |> 
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
  select(c(county_fips, year_week_dt, flu_pred_roll, flu_sum_roll,
           peaks_obs, peaks_pred, is_max_obs, is_max_pred))

county_peaks_obs <- county_peaks_usa |> 
  # Select the peaks
  dplyr::filter(peaks_obs == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_obs == T) |>
  select(c(county_fips, year_week_dt, peaks_obs, is_max_obs)) |>
  rename('week_obs' = 'year_week_dt') |>
  group_by(county_fips) |>
  # Break any ties by taking the earliest peak
  arrange(week_obs) |>
  slice_head(n = 1) |>
  mutate(week_obs_num = week(week_obs)) 

county_peaks_pred <- county_peaks_usa |> 
  # Select the peaks
  dplyr::filter(peaks_pred == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_pred == T) |>
  select(c(county_fips, year_week_dt, peaks_pred, is_max_pred)) |>
  rename('week_pred' = 'year_week_dt') |>
  group_by(county_fips) |>
  # Break any ties by taking the earliest peak
  arrange(week_pred) |>
  slice_head(n = 1) |>
  mutate(week_pred_num = week(week_pred)) 

# Merge obs and pred peaks
county_peaks_both <- left_join(county_peaks_pred, county_peaks_obs,
                               by = c('county_fips' = 'county_fips'))

sd(county_peaks_obs$week_obs_num)

# Change ordering of Epi weeks for calculations and visualizations
county_peaks_both_filt <- county_peaks_both |>
  mutate(week_obs_num_adj = ifelse(week_obs_num > 30, week_obs_num - 52, week_obs_num),
         week_pred_num_adj = ifelse(week_pred_num > 30, week_pred_num - 52, week_pred_num))

# Calculate statistics
county_peaks_both_filt$diff <- as.numeric(as.Date(county_peaks_both_filt$week_pred) - as.Date(county_peaks_both_filt$week_obs))
rmse <- sum(sqrt(sum((county_peaks_both_filt$diff)^2, na.rm = T) / length(county_peaks_both_filt$diff)))
rmse
cor(county_peaks_both_filt$week_pred_num, county_peaks_both_filt$week_obs_num, 
    use = 'complete.obs', method = 'spearman')

# Collapse data for graphing purposes
county_peaks_both_col <- county_peaks_both_filt |>
  mutate(count = 1) |>
  group_by(week_obs_num, week_pred_num) |>
  mutate(count_sum = sum(count)) |>
  distinct(week_obs_num, week_pred_num, week_obs_num_adj, week_pred_num_adj, count_sum)

# Plot
ggplot(county_peaks_both_col, aes(x = week_obs_num, y = week_pred_num)) +
  geom_point(aes(size = count_sum), color = "#e6a532", alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  labs(title = 'Syndromic Prediction vs. Claims\nPeak Timing',
       x = "Obs Peak Week",
       y = "Syndromic Peak Week") + xlim(-10, 20) + ylim(-10, 20) +
  scale_size_continuous(range = c(1, 10), breaks = c(1, 10, 100, 200, 400)) +
  theme_minimal()


library(viridis)
peak_comp_plot <- ggplot(county_peaks_both_col, aes(week_obs_num_adj, week_pred_num_adj, fill= count_sum)) + 
  geom_tile() + xlim(-8, 20) + ylim(-8, 20) +
  scale_fill_viridis(discrete=FALSE, direction = -1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black",
              linewidth = 1) +
  theme_minimal() + 
  labs(title = 'Peak Timing Comparison',
       x = "Claims Peak Epidemiological Week",
       y = "Syndromic Peak Epidemiological Week")
peak_comp_plot












