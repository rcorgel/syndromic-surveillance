################################################################################
# File Name: load_weather_data                                                 #
#                                                                              #
# Purpose:   Load weather data from: DATA/visual_crossing_weather_bulk/.       #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load weather data                                              #
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
library(assertr)
library(ISOweek)
library(lubridate)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

########################
# 2. LOAD WEATHER DATA #
########################

# Load weather data
temp_16 <- read.csv('raw/tot_2016_crosswalk_county.csv')
temp_17 <- read.csv('raw/tot_2017_crosswalk_county.csv')
temp_18 <- read.csv('raw/tot_2018_crosswalk_county.csv')
temp_19 <- read.csv('raw/tot_2019_crosswalk_county.csv')
temp_20 <- read.csv('raw/tot_2020_crosswalk_county.csv')

# Combine to one dataset
weather <- rbind(temp_16, temp_17, temp_18, temp_19, temp_20)
remove(temp_16, temp_17, temp_18, temp_19, temp_20)

# Add leading 0's to county fips
weather <- weather |> 
  mutate(county_fips = ifelse(nchar(county_fips) == 4, 
                              paste0("0", county_fips), county_fips))

# Drop unknown counties and territories
weather$state_fips <- substr(weather$county_fips, 1, 2)
weather <- weather |> filter(county_fips != '99999') |>
  filter(state_fips != '78') |> filter(state_fips != '72') |>
  filter(state_fips != '60')

# Convert date variable from string to date
weather$date <- as.Date(weather$datetime, format = '%Y-%m-%d')

# Collapse data to county group to align with claims
# Load crosswalk
x_walk <- readRDS('./tmp/county_group_xwalk.rds')
x_walk <- x_walk |> rename('county_group' = 'county_fips')
# Join county to county group
weather <- left_join(weather, x_walk[, c(1, 2)],
                     by = c('county_fips' = 'fips'))
# Assert that all counties are merged on
weather |> assert(not_na, county_group)
# Collapse to county group
weather_group <- weather |> group_by(county_group, date) |>
  mutate(feelslike = mean(feelslike),
         dew = mean(dew),
         precip = mean(precip),
         windspeed = mean(windspeed),
         solarradiation = mean(solarradiation),
         sealevelpressure = mean(sealevelpressure)) |>
  distinct(county_group, date, feelslike, dew, precip, windspeed, 
           solarradiation, sealevelpressure)
  
# Collapse data to month level
weather_group$month <- month(weather_group$date)
weather_group$year <- year(weather_group$date)
weather_month <- weather_group |> group_by(county_group, month, year) |>
  mutate(feelslike = mean(feelslike),
         dew = mean(dew),
         precip = mean(precip),
         windspeed = mean(windspeed),
         solarradiation = mean(solarradiation),
         sealevelpressure = mean(sealevelpressure)) |>
  distinct(county_group, year, month, feelslike, dew, precip, windspeed, 
           solarradiation, sealevelpressure)
  
# Collapse to week level
weather_group$week_date_ISO <- as.character(ISOweek(weather_group$date))
weather_group$week_date <- sub("W", "", weather_group$week_date_ISO)
weather_group$week_date <- sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", weather_group$week_date)
weather_group$week_date <- ISOweek2date(weather_group$week_date)
weather_week <- weather_group |> group_by(county_group, week_date) |>
  mutate(feelslike = mean(feelslike),
         dew = mean(dew),
         precip = mean(precip),
         windspeed = mean(windspeed),
         solarradiation = mean(solarradiation),
         sealevelpressure = mean(sealevelpressure)) |>
  distinct(county_group, week_date, feelslike, dew, precip, windspeed, 
           solarradiation, sealevelpressure)

# Save weather data
saveRDS(weather_week, './tmp/county_week_weather.rds') 
saveRDS(weather_month, './tmp/county_month_weather.rds') 

################################################################################
################################################################################
