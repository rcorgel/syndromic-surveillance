################################################################################
# File Name: load_flu_surv_net_data                                            #
#                                                                              #
# Purpose:   Load FluSurv-NET data from the CDC.                               #
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

# Load FLU SURV-NET National Data
surv_nat <- read.csv('./raw/FluSurveillance_Custom_Download_Data.csv', skip = 2, header = TRUE, sep = ',')

# Limit to all all individuals and all subtypes
surv_nat <- surv_nat |> dplyr::filter(AGE.CATEGORY == 'Overall') |>
  dplyr::filter(SEX.CATEGORY == 'Overall') |>
  dplyr::filter(RACE.CATEGORY == 'Overall') |>
  dplyr::filter(VIRUS.TYPE.CATEGORY == 'Overall')

# Convert year week to week date
surv_nat$WEEK_2 <- ifelse(nchar(surv_nat$WEEK) == 1, paste0('0', surv_nat$WEEK), surv_nat$WEEK)
surv_nat$year_week <- paste(surv_nat$YEAR.1, surv_nat$WEEK_2, sep = '-')
surv_nat$week_date = sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", surv_nat$year_week)
surv_nat$week_date = ISOweek2date(surv_nat$week_date)

# Select certain variables of interest
surv_nat_clean <- surv_nat[, c(11, 16)]

# Rename and fill variables
surv_nat_clean <- surv_nat_clean |>
  rename('unweighted_flu_surv' = 'WEEKLY.RATE') |>
  mutate(region = 'United States',
         state_fips = '00')

# Create month and year variables
surv_nat_clean$month <- month(surv_nat_clean$week_date)
surv_nat_clean$year <- year(surv_nat_clean$week_date)

# Subset data to July and calculate the July mean
surv_summer <- surv_nat_clean |> 
  dplyr::filter(month < 11 & month > 9) |> group_by(year) |>
  dplyr::filter(year >= 2016) |>
  dplyr::filter(year < 2020) |>
  mutate(surv_summ_mean = mean(unweighted_flu_surv)) |>
  distinct(year, surv_summ_mean)

# Merge on summer data to the correct season (summer 2016 should align with the 2016-2017 flu season)
surv_nat_clean$year_2 <- ifelse(surv_nat_clean$month < 7, surv_nat_clean$year - 1, surv_nat_clean$year)
surv_nat_clean <- left_join(surv_nat_clean, surv_summer, by = c('year_2' = 'year'))

# Fill in missing weeks
ili_weeks <- readRDS('./tmp/ili_net_national.rds')
surv_nat_clean <- left_join(ili_weeks[, c(1, 3)], surv_nat_clean, by = c('week_date' = 'week_date'))
surv_nat_clean <- surv_nat_clean |>
  dplyr::select(-c(region.y)) |>
  rename('region' = 'region.x')

# Calculate rolling average
surv_nat_clean <- surv_nat_clean |>
  mutate(roll_surv = rollmean(unweighted_flu_surv, 3, fill = NA, align = 'right'))

# Limit data to 2016 to 2020
surv_nat_clean <- surv_nat_clean |> 
  dplyr::filter(week_date < as.Date('2020-03-01')) |>
  dplyr::filter(week_date > as.Date('2016-08-31'))

# Calculate Z score
surv_nat_clean$z_score <- (surv_nat_clean$unweighted_flu_surv - 
                            surv_nat_clean$surv_summ_mean) / sd(surv_nat_clean$unweighted_flu_surv, na.rm = TRUE)
surv_nat_clean$z_score_roll <- (surv_nat_clean$roll_surv - 
                                  surv_nat_clean$surv_summ_mean) / sd(surv_nat_clean$roll_surv, na.rm = TRUE)

# Finalize national data
surv_final <- surv_nat_clean[, c(1, 2, 3, 4, 9, 10, 11)]
surv_final <- surv_final |>
  rename('value' = 'unweighted_flu_surv',
         'value_roll' = 'roll_surv')
surv_final$Source <- 'FluSurv-NET'

# A quick visual
ggplot(surv_final) + 
  geom_line(aes(x = week_date, y = value), color = 'blue') +
  geom_line(aes(x = week_date, y = value_roll), color = 'red')

# Save
saveRDS(surv_final, './tmp/flu_surv_net_national.rds')

######################
# 3. LOAD STATE DATA #
######################

# Load FLU SURV-NET State Data
surv_state <- read.csv('./raw/FluSurveillance_Custom_Download_Data (3).csv', skip = 2, header = TRUE, sep = ',')

# Limit to all all individuals and all subtypes
surv_state <- surv_state |> dplyr::filter(AGE.CATEGORY == 'Overall') |>
  dplyr::filter(SEX.CATEGORY == 'Overall') |>
  dplyr::filter(RACE.CATEGORY == 'Overall') |>
  dplyr::filter(VIRUS.TYPE.CATEGORY == 'Overall')

# Convert year week to week date
surv_state$WEEK_2 <- ifelse(nchar(surv_state$WEEK) == 1, paste0('0', surv_state$WEEK), surv_state$WEEK)
surv_state$year_week <- paste(surv_state$YEAR.1, surv_state$WEEK_2, sep = '-')
surv_state$week_date = sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", surv_state$year_week)
surv_state$week_date = ISOweek2date(surv_state$week_date)

# Convert state name to state fips
surv_state$state_fips <- fips(surv_state$CATCHMENT)
# There are two New York sites
surv_state$state_fips <- ifelse(surv_state$CATCHMENT == 'New York - Albany', '36', surv_state$state_fips)
surv_state$state_fips <- ifelse(surv_state$CATCHMENT == 'New York - Rochester', '36', surv_state$state_fips)
surv_state <- surv_state |> dplyr::filter(CATCHMENT != 'Entire Network')

# Select certain variables of interest
surv_state_clean <- surv_state[, c(1, 11, 16, 17)]

# Rename and fill variables
surv_state_clean <- surv_state_clean |>
  rename('unweighted_flu_surv' = 'WEEKLY.RATE',
         'region' = 'CATCHMENT')

# Create month and year variables
surv_state_clean$month <- month(surv_state_clean$week_date)
surv_state_clean$year <- year(surv_state_clean$week_date)

# Subset data to July and calculate the July mean
surv_state_summer <- surv_state_clean |> 
  dplyr::filter(month < 11 & month > 9) |> group_by(year, state_fips, region) |>
  dplyr::filter(year >= 2016) |>
  dplyr::filter(year < 2020) |>
  mutate(surv_summ_mean = mean(as.numeric(unweighted_flu_surv), na.rm = TRUE)) |>
  distinct(year, state_fips, region, surv_summ_mean)

# Merge on summer data to the correct season (summer 2016 should align with the 2016-2017 flu season)
surv_state_clean$year_2 <- ifelse(surv_state_clean$month < 7, surv_state_clean$year - 1, surv_state_clean$year)
surv_state_clean <- left_join(surv_state_clean, surv_state_summer, by = c('year_2' = 'year',
                                                                          'state_fips' = 'state_fips',
                                                                          'region' = 'region'))

# Fill in missing weeks
ili_weeks <- readRDS('./tmp/ili_net_national.rds')
ili_all_state <- NULL
for (i in unique(surv_state_clean$region)) {
  print(i)
  ili_state <- ili_weeks[,c(1, 3)]
  ili_state$region <- i
  ili_all_state <- rbind(ili_all_state, ili_state)
}
surv_state_clean <- left_join(ili_all_state, surv_state_clean, by = c('week_date' = 'week_date',
                                                                      'region' = 'region'))

# Calculate rolling average
surv_state_clean <- surv_state_clean |> group_by(region) |>
  mutate(roll_surv = rollmean(as.numeric(unweighted_flu_surv), 3, fill = NA, align = 'right'))

# Limit data to 2016 to 2020
surv_state_clean <- surv_state_clean |> 
  dplyr::filter(week_date < as.Date('2020-03-01')) |>
  dplyr::filter(week_date > as.Date('2016-08-31'))

# Calculate Z score
surv_state_clean <- surv_state_clean |>
  group_by(state_fips, region) |>
  mutate(z_score = (as.numeric(unweighted_flu_surv) - surv_summ_mean) / 
           sd(unweighted_flu_surv, na.rm = TRUE),
         z_score_roll = (as.numeric(roll_surv) - surv_summ_mean) / 
           sd(roll_surv, na.rm = TRUE))

# Finalize state data
surv_state_final <- surv_state_clean[, c(1, 2, 3, 4, 9, 10, 11)]
surv_state_final <- surv_state_final |>
  rename('value' = 'unweighted_flu_surv',
         'value_roll' = 'roll_surv')
surv_state_final$Source <- 'FluSurv-NET'

# A quick visual
# All states plus NYC have data besides Florida
ggplot(surv_state_final) + 
  geom_line(aes(x = week_date, y = z_score), color = 'blue') +
  geom_line(aes(x = week_date, y = z_score_roll), color = 'red') +
  facet_wrap(~ region, nrow = 4, scales = 'free')

# Save
saveRDS(surv_state_final, './tmp/flu_surv_net_state.rds')

################################################################################
################################################################################
