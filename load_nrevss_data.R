################################################################################
# File Name: load_nrevss_data                                                  #
#                                                                              #
# Purpose:   Load NREVSS data from the CDC.                                    #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load national data                                             #
#            3. Load state data                                                #
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
library(ISOweek)
library(zoo)
library(usmap)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

#########################
# 2. LOAD NATIONAL DATA #
#########################

# Load NREVSS National Data
nrevss_nat <- read.csv('./raw/WHO_NREVSS_Clinical_Labs.csv', skip = 1, header = TRUE, sep = ',')

# Convert year week to week date
nrevss_nat$WEEK_2 <- ifelse(nchar(nrevss_nat$WEEK) == 1, paste0('0', nrevss_nat$WEEK), nrevss_nat$WEEK)
nrevss_nat$year_week <- paste(nrevss_nat$YEAR, nrevss_nat$WEEK_2, sep = '-')
nrevss_nat$week_date = sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", nrevss_nat$year_week)
nrevss_nat$week_date = ISOweek2date(nrevss_nat$week_date)

# Select certain variables of interest
nrevss_nat_clean <- nrevss_nat[, c(2, 8, 13)]

# Rename and fill variables
nrevss_nat_clean <- nrevss_nat_clean |>
  rename('region' = 'REGION',
         'percent_positive_nrevss' = 'PERCENT.POSITIVE') |>
  mutate(region = 'United States',
         state_fips = '00')

# Create month and year variables
nrevss_nat_clean$month <- month(nrevss_nat_clean$week_date)
nrevss_nat_clean$year <- year(nrevss_nat_clean$week_date)

# Subset data to July and calculate the July mean
nrevss_summer <- nrevss_nat_clean |> 
  dplyr::filter(month < 9 & month > 5) |> group_by(year) |>
  dplyr::filter(year >= 2016) |>
  dplyr::filter(year < 2020) |>
  mutate(nrevss_summ_mean = mean(percent_positive_nrevss)) |>
  distinct(year, nrevss_summ_mean)

# Merge on summer data to the correct season (summer 2016 should align with the 2016-2017 flu season)
nrevss_nat_clean$year_2 <- ifelse(nrevss_nat_clean$month < 6, nrevss_nat_clean$year - 1, nrevss_nat_clean$year)
nrevss_nat_clean <- left_join(nrevss_nat_clean, nrevss_summer, by = c('year_2' = 'year'))

# Calculate rolling average
nrevss_nat_clean <- nrevss_nat_clean |>
  mutate(roll_nrevss = rollmean(percent_positive_nrevss, 4, fill = NA, align = 'right'))

# Limit data to 2016 to 2020
nrevss_nat_clean <- nrevss_nat_clean |> 
  dplyr::filter(week_date < as.Date('2020-03-01')) |>
  dplyr::filter(week_date > as.Date('2016-08-31'))

# Calculate Z score
nrevss_nat_clean$z_score <- (nrevss_nat_clean$percent_positive_nrevss - 
                                nrevss_nat_clean$nrevss_summ_mean) / sd(nrevss_nat_clean$percent_positive_nrevss)
nrevss_nat_clean$z_score_roll <- (nrevss_nat_clean$roll_nrevss - 
                                    nrevss_nat_clean$nrevss_summ_mean) / sd(nrevss_nat_clean$roll_nrevss, na.rm = TRUE)

# Finalize national data
nrevss_final <- nrevss_nat_clean[, c(1, 2, 3, 4, 9, 10, 11)]
nrevss_final <- nrevss_final |>
  rename('value' = 'percent_positive_nrevss',
         'value_roll' = 'roll_nrevss')
nrevss_final$Source <- 'NREVSS'

# A quick visual
ggplot(nrevss_final) + 
  geom_line(aes(x = week_date, y = value), color = 'blue') +
  geom_line(aes(x = week_date, y = value_roll), color = 'red')

# Save
saveRDS(nrevss_final, './tmp/nrevss_national.rds')

######################
# 3. LOAD STATE DATA #
######################

# Load NREVSS State Data
nrevss_state <- read.csv('./raw/ICL_NREVSS_Clinical_Labs_State.csv', skip = 1, header = TRUE, sep = ',')

# Convert year week to week date
nrevss_state$WEEK_2 <- ifelse(nchar(nrevss_state$WEEK) == 1, paste0('0', nrevss_state$WEEK), nrevss_state$WEEK)
nrevss_state$year_week <- paste(nrevss_state$YEAR, nrevss_state$WEEK_2, sep = '-')
nrevss_state$week_date = sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", nrevss_state$year_week)
nrevss_state$week_date = ISOweek2date(nrevss_state$week_date)

# Convert state name to state fips
nrevss_state$state_fips <- fips(nrevss_state$REGION)

# Select certain variables of interest
nrevss_state_clean <- nrevss_state[, c(2, 8, 13, 14)]

# Rename and fill variables
nrevss_state_clean <- nrevss_state_clean |>
  rename('region' = 'REGION',
         'percent_positive_nrevss' = 'PERCENT.POSITIVE') 

# Create month and year variables
nrevss_state_clean$month <- month(nrevss_state_clean$week_date)
nrevss_state_clean$year <- year(nrevss_state_clean$week_date)

# Subset data to July and calculate the July mean
# There are some NAs
nrevss_state_summer <- nrevss_state_clean |> 
  dplyr::filter(month < 9 & month > 5) |> group_by(year, state_fips) |>
  dplyr::filter(year >= 2016) |>
  dplyr::filter(year < 2020) |>
  mutate(nrevss_summ_mean = mean(as.numeric(percent_positive_nrevss), na.rm = TRUE)) |>
  distinct(year, state_fips, nrevss_summ_mean)

# Merge on summer data to the correct season (summer 2016 should align with the 2016-2017 flu season)
nrevss_state_clean$year_2 <- ifelse(nrevss_state_clean$month < 6, nrevss_state_clean$year - 1, nrevss_state_clean$year)
nrevss_state_clean <- left_join(nrevss_state_clean, nrevss_state_summer, by = c('year_2' = 'year',
                                                                                'state_fips' = 'state_fips'))

# Calculate rolling average
# There are some NAs
nrevss_state_clean <- nrevss_state_clean |> group_by(state_fips) |>
  mutate(roll_nrevss = rollmean(as.numeric(percent_positive_nrevss), 4, fill = NA, align = 'right'))

# Limit data to 2016 to 2020
nrevss_state_clean <- nrevss_state_clean |> 
  dplyr::filter(week_date < as.Date('2020-03-01')) |>
  dplyr::filter(week_date > as.Date('2016-08-31'))

# Calculate Z score
nrevss_state_clean <- nrevss_state_clean |>
  group_by(state_fips) |>
  mutate(z_score = (as.numeric(percent_positive_nrevss) - nrevss_summ_mean) / 
           sd(as.numeric(percent_positive_nrevss), na.rm = TRUE),
         z_score_roll = (as.numeric(roll_nrevss) - nrevss_summ_mean) / 
           sd(roll_nrevss, na.rm = TRUE))

# Finalize state data
nrevss_state_final <- nrevss_state_clean[, c(1, 2, 3, 4, 9, 10, 11)]
nrevss_state_final <- nrevss_state_final |>
  rename('value' = 'percent_positive_nrevss',
         'value_roll' = 'roll_nrevss') |>
  mutate(value = as.numeric(value)) |>
  dplyr::filter(!is.na(state_fips))
nrevss_state_final$Source <- 'NREVSS'

# A quick visual
# All states plus NYC have data besides Florida
ggplot(nrevss_state_final) + 
  geom_line(aes(x = week_date, y = value), color = 'blue') +
  geom_line(aes(x = week_date, y = value_roll), color = 'red') +
  facet_wrap(~ region, nrow = 10, scales = 'free')

# Save
saveRDS(nrevss_state_final, './tmp/nrevss_state.rds')

################################################################################
################################################################################
