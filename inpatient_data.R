################################################################################
# File Name: inpatient_data                                                    #
#                                                                              #
# Purpose:   Load influenza, COVID-19, and RSV inpatient claims data.          #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load, explore, and save data                                   #
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

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1CgyyYhIWFVG2aMITI-HBpKamvwwpQ_ED/CHC_data/inpatient/2024-06-21/')

################
# 2. LOAD DATA #
################

# Load data
flu <- read.csv('county_week_conf_flu_inpatient.csv')
covid <- read.csv('county_week_covid_inpatient.csv')
rsv <- read.csv('county_week_rsv_inpatient.csv')

# Get list of unique counties
unique_county <- data.frame(unique(flu$county_fips)) |> 
  rename('county_fips' = 'unique.flu.county_fips.')
unique_county_1 <- data.frame(unique(covid$county_fips)) |> 
  rename('county_fips' = 'unique.covid.county_fips.')
unique_county_2 <- data.frame(unique(rsv$county_fips)) |> 
  rename('county_fips' = 'unique.rsv.county_fips.')
unique_counties <- rbind(unique_county, unique_county_1, unique_county_2) |>
  distinct(county_fips)

# Get a unique list of weeks
unique_week <- data.frame(unique(flu$year_week)) |> 
  rename('year_week' = 'unique.flu.year_week.')
unique_week_1 <- data.frame(unique(covid$year_week)) |> 
  rename('year_week' = 'unique.covid.year_week.')
unique_week_2 <- data.frame(unique(rsv$year_week)) |> 
  rename('year_week' = 'unique.rsv.year_week.')
unique_weeks <- rbind(unique_week, unique_week_1, unique_week_2) |>
  distinct(year_week)

# Convert week number to date 
unique_weeks <- unique_weeks |>
  mutate(year_week_date = 
           as.Date(paste(year_week, 1, sep = '-'), format = '%Y-%U-%u'))





