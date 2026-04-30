################################################################################
# File Name: load_ili_net_data                                                 #
#                                                                              #
# Purpose:   Load ILI-NET data from the CDC.                                   #
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

# Load ILI NET National Data
ili_nat <- read.csv('./raw/ILINet.csv', skip = 1, header = TRUE, sep = ',')

# Convert year week to week date
ili_nat$WEEK_2 <- ifelse(nchar(ili_nat$WEEK) == 1, paste0('0', ili_nat$WEEK), ili_nat$WEEK)
ili_nat$year_week <- paste(ili_nat$YEAR, ili_nat$WEEK_2, sep = '-')
ili_nat$week_date = sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", ili_nat$year_week)
ili_nat$week_date = ISOweek2date(ili_nat$week_date)

# Select certain variables of interest
ili_nat_clean <- ili_nat[, c(2, 5, 6, 18)]

# Rename and fill variables
ili_nat_clean <- ili_nat_clean |>
  rename('region' = 'REGION',
         'unweighted_ili_net' = 'X.UNWEIGHTED.ILI',
         'weighted_ili_net' = 'X..WEIGHTED.ILI') |>
  mutate(region = 'United States',
         state_fips = '00')

# Create month and year variables
ili_nat_clean$month <- month(ili_nat_clean$week_date)
ili_nat_clean$year <- year(ili_nat_clean$week_date)

# Subset data to July and calculate the July mean
ili_summer <- ili_nat_clean |> 
  dplyr::filter(month < 9 & month > 5) |> group_by(year) |>
  dplyr::filter(year >= 2016) |>
  dplyr::filter(year < 2025) |>
  mutate(ili_summ_mean = mean(unweighted_ili_net)) |>
  distinct(year, ili_summ_mean)

# Merge on summer data to the correct season (summer 2016 should align with the 2016-2017 flu season)
ili_nat_clean$year_2 <- ifelse(ili_nat_clean$month < 6, ili_nat_clean$year - 1, ili_nat_clean$year)
ili_nat_clean <- left_join(ili_nat_clean, ili_summer, by = c('year_2' = 'year'))

# Calculate rolling average
ili_nat_clean <- ili_nat_clean |>
  mutate(roll_ili = rollmean(unweighted_ili_net, 4, fill = NA, align = 'right'))

# Limit data to 2016 to 2020
ili_nat_clean <- ili_nat_clean |> 
  dplyr::filter(week_date < as.Date('2023-08-30')) |>
  dplyr::filter(week_date > as.Date('2016-08-31'))

# Calculate Z score
ili_nat_clean$z_score <- (ili_nat_clean$unweighted_ili_net - 
                            ili_nat_clean$ili_summ_mean) / sd(ili_nat_clean$unweighted_ili_net)
ili_nat_clean$z_score_roll <- (ili_nat_clean$roll_ili - 
                            ili_nat_clean$ili_summ_mean) / sd(ili_nat_clean$roll_ili)

# Finalize national data
ili_final <- ili_nat_clean[, c(1, 3, 4, 5, 10, 11, 12)]
ili_final <- ili_final |>
  rename('value' = 'unweighted_ili_net',
         'value_roll' = 'roll_ili')
ili_final$Source <- 'ILI-NET'

# A quick visual
ggplot(ili_final) + 
  geom_line(aes(x = week_date, y = value), color = 'blue') +
  geom_line(aes(x = week_date, y = value_roll), color = 'red')

# Save
saveRDS(ili_final, './tmp/ili_net_national.rds')

######################
# 3. LOAD STATE DATA #
######################

# Load ILI NET State Data
ili_state <- read.csv('./raw/ILINet_State.csv', skip = 1, header = TRUE, sep = ',')

# Convert year week to week date
ili_state$WEEK_2 <- ifelse(nchar(ili_state$WEEK) == 1, paste0('0', ili_state$WEEK), ili_state$WEEK)
ili_state$year_week <- paste(ili_state$YEAR, ili_state$WEEK_2, sep = '-')
ili_state$week_date = sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", ili_state$year_week)
ili_state$week_date = ISOweek2date(ili_state$week_date)

# Convert state name to state fips
ili_state$state_fips <- fips(ili_state$REGION)

# Add NYC to state fips
ili_state$state_fips = ifelse(ili_state$REGION == 'New York City', '99', ili_state$state_fips)

# Select certain variables of interest
ili_state_clean <- ili_state[, c(2, 6, 18, 19)]

# Rename variables
ili_state_clean <- ili_state_clean |>
  rename('region' = 'REGION',
         'unweighted_ili_net' = 'X.UNWEIGHTED.ILI')

# Add year and month variables
ili_state_clean$month <- month(ili_state_clean$week_date)
ili_state_clean$year <- year(ili_state_clean$week_date)

# Subset data to July and calculate the July mean
# NAs are produced for Florida
ili_state_summer <- ili_state_clean |> 
  dplyr::filter(month < 9 & month > 5) |> group_by(year, state_fips) |>
  dplyr::filter(year >= 2016) |>
  mutate(ili_summ_mean = mean(as.numeric(unweighted_ili_net), na.rm = TRUE)) |>
  distinct(year, state_fips, ili_summ_mean)

# Merge on summer data to the correct season (summer 2016 should align with the 2016-2017 flu season)
ili_state_clean$year_2 <- ifelse(ili_state_clean$month < 5, ili_state_clean$year - 1, ili_state_clean$year)
ili_state_clean <- left_join(ili_state_clean, ili_state_summer, by = c('year_2' = 'year',
                                                                       'state_fips' = 'state_fips'))

# Calculate rolling average
# NAs are produced for Florida
ili_state_clean <- ili_state_clean |> group_by(state_fips) |>
  mutate(roll_ili = rollmean(as.numeric(unweighted_ili_net), 4, fill = NA, align = 'right'))

# Limit data to 2016 to 2020
ili_state_clean <- ili_state_clean |> 
  dplyr::filter(week_date > as.Date('2016-08-31'))

# Calculate Z score
# NAs are produced for Florida
ili_state_clean <- ili_state_clean |>
  group_by(state_fips) |>
  mutate(z_score = (as.numeric(unweighted_ili_net) - ili_summ_mean) / 
           sd(as.numeric(unweighted_ili_net), na.rm = TRUE),
         z_score_roll = (roll_ili - ili_summ_mean) / 
           sd(roll_ili, na.rm = TRUE))

# Finalize state data
ili_state_final <- ili_state_clean[, c(1, 2, 3, 4, 9, 10, 11)]
ili_state_final <- ili_state_final |>
  rename('value' = 'unweighted_ili_net',
         'value_roll' = 'roll_ili') |>
  mutate(value = as.numeric(value)) |>
  dplyr::filter(!is.na(state_fips))
ili_state_final$Source <- 'ILI-NET'

# A quick visual
# All states plus NYC have data besides Florida
ggplot(ili_state_final) + 
  geom_line(aes(x = week_date, y = z_score), color = 'blue') +
  geom_line(aes(x = week_date, y = z_score_roll), color = 'red') +
  facet_wrap(~ state_fips, nrow = 10, scales = 'free')

# Save
saveRDS(ili_state_final, './tmp/ili_net_state.rds')

################################################################################
################################################################################
