################################################################################
# File Name: load_data_p1                                                      #
#                                                                              #
# Purpose:   Load Influenza, COVID-19, and RSV medical claims data.            #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load, append, check, and save data                             #
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

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

################
# 2. LOAD DATA #
################

## INFLUENZA ##

# Influenza symptom data
flu_16_17 <- read.csv('raw/flu_2016-09_2017-08_1b.csv')
flu_17_18 <- read.csv('raw/flu_2017-09_2018-08_1b.csv')
flu_18_19 <- read.csv('raw/flu_2018-09_2019-08_1b.csv')
flu_19_20 <- read.csv('raw/flu_2019-09_2020-08_1b.csv')
# Don't include 2020-2023 due to low flu circulation & covid-19 co-circulation
# Append
flu <- rbind(flu_16_17, flu_17_18, flu_18_19, flu_19_20)
remove(flu_16_17, flu_17_18, flu_18_19, flu_19_20)

# Run sanity checks on the data
flu |>
  # There should be 48 months
  verify(length(unique(year_month)) == 48) |>
  # There should be more than 2500 county groups
  verify(length(unique(county_fips)) > 2500) |>
  # There are unknown age groups
  verify(length(unique(age_grp)) == 7) |>
  # There are unknown genders
  verify(length(unique(patient_gender_code)) == 3) |>
  # Check disease, symptom, and count vars are not missing
  assert(not_na, flu:patient_count)

# Save full data
saveRDS(flu, './tmp/flu_p1_dat_raw.rds')
remove(flu)

## COVID-19 ##

# Covid symptom data
covid_wild <- read.csv('raw/covid_wild_1b.csv')
covid_alpha <- read.csv('raw/covid_alpha_1b.csv')
covid_delta <- read.csv('raw/covid_delta_1b.csv')
covid_omicron <- read.csv('raw/covid_omicron_1b.csv')
# Append
covid <- rbind(covid_wild, covid_alpha, covid_delta, covid_omicron)
remove(covid_wild, covid_alpha, covid_delta, covid_omicron)

# Run sanity checks on the data
covid |>
  # There should be 48 months
  verify(length(unique(year_month)) == 48) |>
  # There should be more than 2500 county groups
  verify(length(unique(county_fips)) > 2500) |>
  # There are unknown age groups
  verify(length(unique(age_grp)) == 7) |>
  # There are unknown genders
  verify(length(unique(patient_gender_code)) == 3) |>
  # Check disease, symptom, and count vars are not missing
  assert(not_na, covid:patient_count)

# Save full data
saveRDS(covid, './tmp/covid_p1_dat_raw.rds')
remove(covid)

## RSV ##

# RSV symptom data
rsv_16_17 <- read.csv('raw/rsv_2016-09_2017-08_1d.csv')
rsv_17_18 <- read.csv('raw/rsv_2017-09_2018-08_1d.csv')
rsv_18_19 <- read.csv('raw/rsv_2018-09_2019-08_1d.csv')
rsv_19_20 <- read.csv('raw/rsv_2019-09_2020-08_1d.csv')
# Don't include 2020-2023 due to low rsv circulation & covid-19 co-circulation
# Append
rsv <- rbind(rsv_16_17, rsv_17_18, rsv_18_19, rsv_19_20)
remove(rsv_16_17, rsv_17_18, rsv_18_19, rsv_19_20)

# Run sanity checks on the data
rsv |>
  # There should be 48 months
  verify(length(unique(year_month)) == 48) |>
  # There should be more than 2500 county groups
  verify(length(unique(county_fips)) > 2500) |>
  # There are unknown age groups, only 6 groups for RSV
  verify(length(unique(age_grp)) == 6) |>
  # There are unknown genders
  verify(length(unique(patient_gender_code)) == 3) |>
  # Check disease, symptom, and count vars are not missing
  assert(not_na, rsv:patient_count)

# Save full data
saveRDS(rsv, './tmp/rsv_p1_dat_raw.rds')
remove(rsv)

################################################################################
################################################################################
