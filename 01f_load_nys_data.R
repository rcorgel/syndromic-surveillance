################################################################################
# File Name: load_nys_data                                                     #
#                                                                              #
# Purpose:   Load NYS County data from the CDC.                                #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load data                                                      #
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

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

################
# 2. LOAD DATA #
################

# Load NYS Data
nys <- read.csv('./raw/Influenza_Laboratory-Confirmed_Cases_by_County__Beginning_2009-10_Season_20250715.csv', header = TRUE, sep = ',')

# Convert year week to week date
nys$year <- ifelse(nys$CDC.Week < 35, substr(nys$Season, 6, 9), substr(nys$Season, 1, 4))
nys$WEEK_2 <- ifelse(nchar(nys$CDC.Week) == 1, paste0('0', nys$CDC.Week), nys$CDC.Week)
nys$year_week <- paste(nys$year, nys$WEEK_2, sep = '-')
nys$week_date = sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", nys$year_week)
nys$week_date = ISOweek2date(nys$week_date)

# Select certain variables of interest
nys_clean <- nys[, c(3, 7, 9, 13)]

# Rename and fill variables
nys_clean <- nys_clean |>
  rename('county' = 'County',
         'county_fips' = 'FIPS',
         'nys_count' = 'Count')

# Collapse to the county week level
nys_clean <- nys_clean |> group_by(county_fips, week_date) |>
  mutate(case_count = sum(nys_count)) |>
  distinct(county_fips, week_date, county, case_count)

# Create month and year variables
nys_clean$month <- month(nys_clean$week_date)
nys_clean$year <- year(nys_clean$week_date)

# Subset data to July and calculate the July mean
# There are some NAs
nys_summer <- nys_clean |> 
  dplyr::filter(month < 11 & month > 9) |> group_by(year, county_fips) |>
  dplyr::filter(year >= 2016) |>
  dplyr::filter(year < 2020) |>
  mutate(nys_summ_mean = mean(as.numeric(case_count), na.rm = TRUE)) |>
  distinct(year, county_fips, nys_summ_mean)

# Merge on summer data to the correct season (summer 2016 should align with the 2016-2017 flu season)
nys_clean$year_2 <- ifelse(nys_clean$month < 7, nys_clean$year - 1, nys_clean$year)
nys_clean <- left_join(nys_clean, nys_summer, by = c('year_2' = 'year',
                                                     'county_fips' = 'county_fips'))

# Fill in missing weeks
ili_weeks <- readRDS('./tmp/ili_net_national.rds')
ili_all_county <- NULL
for (i in unique(nys_clean$county)) {
  print(i)
  ili_county <- ili_weeks[,c(1, 3)]
  ili_county$county <- i
  ili_all_county <- rbind(ili_all_county, ili_county)
}
nys_clean <- left_join(ili_all_county, nys_clean, by = c('week_date' = 'week_date',
                                                         'county' = 'county'))

# Calculate rolling average
# There are some NAs
nys_clean <- nys_clean |> group_by(county_fips) |>
  mutate(roll_nys = rollmean(as.numeric(case_count), 4, fill = NA, align = 'right'))

# Limit data to 2016 to 2020
nys_clean <- nys_clean |> 
  dplyr::filter(week_date > as.Date('2016-08-31'))

# Calculate Z score
nys_clean <- nys_clean |>
  group_by(county_fips) |>
  mutate(z_score = (as.numeric(case_count) - nys_summ_mean) / 
           sd(as.numeric(roll_nys), na.rm = TRUE),
         z_score_roll = (as.numeric(roll_nys) - nys_summ_mean) / 
           sd(roll_nys, na.rm = TRUE))

# Finalize state data
nys_final <- nys_clean[, c(2, 3, 4, 5, 10, 11, 12)]
nys_final <- nys_final |>
  rename('value' = 'case_count',
         'value_roll' = 'roll_nys') 
nys_final$Source <- 'NYS'

# A quick visual
# All states plus NYC have data besides Florida
ggplot(nys_final) + 
  geom_line(aes(x = week_date, y = value), color = 'blue') +
  geom_line(aes(x = week_date, y = value_roll), color = 'red') +
  facet_wrap(~ county, nrow = 10, scales = 'free')

# Save
saveRDS(nys_final, './tmp/nys_county.rds')  

################################################################################
################################################################################
