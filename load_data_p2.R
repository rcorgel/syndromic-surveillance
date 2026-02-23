################################################################################
# File Name: load_data                                                         #
#                                                                              #
# Purpose:   Load and process Influenza medical claims data from 2016-2020.    #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load, append, check, and save data                             #
#                                                                              #
# Project:   Syndromic Influenza                                               #
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
library(arrow)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

################
# 2. LOAD DATA #
################

## INFLUENZA ##

# Influenza symptom data
# Convert to parquet
flu_16_17 <- read_csv_arrow('raw/flu_2016-09_2017-08_p2_imputed.csv')
write_parquet(flu_16_17, 'tmp/flu_2016-09_2017-08_p2_imputed.parquet')

flu_17_18 <- read_csv_arrow('raw/flu_2017-09_2018-08_p2_imputed.csv')
write_parquet(flu_17_18, 'tmp/flu_2017-09_2018-08_p2_imputed.parquet')

flu_18_19 <- read_csv_arrow('raw/flu_2018-09_2019-08_p2_imputed.csv')
write_parquet(flu_18_19, 'tmp/flu_2018-09_2019-08_p2_imputed.parquet')

flu_19_20 <- read_csv_arrow('raw/flu_2019-09_2020-08_p2_imputed.csv')
write_parquet(flu_19_20, 'tmp/flu_2019-09_2020-08_p2_imputed.parquet')

# Read in parquet and convert to all strings
flu_16_17 <- read_parquet('tmp/flu_2016-09_2017-08_p2_imputed.parquet')
flu_16_17 <- flu_16_17 |>
  mutate(across(everything(), as.character))
flu_17_18 <- read_parquet('tmp/flu_2017-09_2018-08_p2_imputed.parquet')
flu_17_18 <- flu_17_18 |>
  mutate(across(everything(), as.character))
flu_18_19 <- read_parquet('tmp/flu_2018-09_2019-08_p2_imputed.parquet')
flu_18_19 <- flu_18_19 |>
  mutate(across(everything(), as.character))
flu_19_20 <- read_parquet('tmp/flu_2019-09_2020-08_p2_imputed.parquet')
flu_19_20 <- flu_19_20 |>
  mutate(across(everything(), as.character))

# Don't include 2020-2023 due to low flu circulation & covid-19 co-circulation

# Combine into one data set
table_1 <- arrow_table(flu_16_17)
table_2 <- arrow_table(flu_17_18)
table_3 <- arrow_table(flu_18_19)
table_4 <- arrow_table(flu_19_20)
flu_table <- concat_tables(table_1, table_2, table_3, table_4)

# Remove data
remove(flu_16_17, flu_17_18, flu_18_19, flu_19_20)

# Store checks as you go
checks_passed <- TRUE

# Check unique counts using Arrow
n_weeks <- flu_table |>
  distinct(year_week_dt) |>
  count() |>
  pull(n, as_vector = TRUE) |>
  as.numeric()

checks_passed <- checks_passed && (n_weeks > 204 && n_weeks < 220)

n_counties <- flu_table |> 
  distinct(county_fips) |> 
  count() |> 
  pull(n, as_vector = TRUE) |> 
  as.numeric()

checks_passed <- checks_passed && (n_counties > 2500)

# Check for NAs
na_check <- flu_table |> 
  filter(is.na(flu) | is.na(patient_count)) |> 
  count() |> 
  pull(n) |> 
  as.numeric()

checks_passed <- checks_passed && (na_check == 0)

if (!checks_passed) {
  stop("Data validation failed")
}





# Check data
flu_table |>
  collect()
  # There should be 52-54 weeks
  verify(length(unique(year_week_dt)) > 204) |>
  verify(length(unique(year_week_dt)) < 220) |>
  # There should be more than 2500 county groups
  verify(length(unique(county_fips)) > 2500) |>
  # There are unknown age groups, 6 plus 1
  verify(length(unique(age_grp)) == 7) |>
  # There are unknown genders, 2 plus 1
  verify(length(unique(patient_gender_code)) == 3) |>
  # Check disease, symptom, and count vars are not missing
  assert(not_na, flu:patient_count)



# Loop through flu data and perform sanity checks
for (i in 1:4) {
  flu_check <- as.data.frame(flu[i])
  # Run sanity checks on the data
  flu_check |>
    # There should be 52-54 weeks
    verify(length(unique(year_week)) > 51) |>
    verify(length(unique(year_week)) < 55) |>
    # There should be more than 2500 county groups
    verify(length(unique(county_fips)) > 2500) |>
    # There are unknown age groups, 6 plus 1
    verify(length(unique(age_grp)) == 7) |>
    # There are unknown genders, 2 plus 1
    verify(length(unique(patient_gender_code)) == 3) |>
    # Check disease, symptom, and count vars are not missing
    assert(not_na, flu:patient_count)
  remove(flu_check)
}

# Save full data
saveRDS(flu, './tmp/flu_p2_dat_raw_list.rds')
remove(flu)

## COVID-19 ##

# Covid symptom data
covid_wild <- read.csv('raw/covid_wild_p2.csv')
covid_alpha <- read.csv('raw/covid_alpha_p2.csv')
covid_delta <- read.csv('raw/covid_delta_p2.csv')
covid_omicron <- read.csv('raw/covid_omicron_p2.csv')
# Append
covid <- rbind(covid_wild, covid_alpha, covid_delta, covid_omicron)
remove(covid_wild, covid_alpha, covid_delta, covid_omicron)

# Run sanity checks on the data
covid |>
  # There should be 48 months
  verify(length(unique(year_month)) == 210) |>
  # There should be more than 2500 county groups
  verify(length(unique(county_fips)) > 2500) |>
  # There are unknown age groups
  verify(length(unique(age_grp)) == 7) |>
  # There are unknown genders
  verify(length(unique(patient_gender_code)) == 3) |>
  # Check disease, symptom, and count vars are not missing
  assert(not_na, covid:patient_count)

# Save full data
saveRDS(covid, './tmp/covid_p2_dat_raw.rds')
remove(covid)

## RSV ##

# RSV symptom data
rsv_16_17 <- read.csv('raw/rsv_2016-09_2017-08_p2.csv')
rsv_17_18 <- read.csv('raw/rsv_2017-09_2018-08_p2.csv')
rsv_18_19 <- read.csv('raw/rsv_2018-09_2019-08_p2.csv')
rsv_19_20 <- read.csv('raw/rsv_2019-09_2020-08_p2.csv')
# Don't include 2020-2023 due to low rsv circulation & covid-19 co-circulation
# Append
rsv <- rbind(rsv_16_17, rsv_17_18, rsv_18_19, rsv_19_20)
remove(rsv_16_17, rsv_17_18, rsv_18_19, rsv_19_20)

# Run sanity checks on the data
rsv |>
  # There should be > 200 weeks
  verify(length(unique(year_week)) == 210) |>
  # There should be more than 2500 county groups
  verify(length(unique(county_fips)) > 2500) |>
  # There are unknown age groups, only 6 groups for RSV
  verify(length(unique(age_grp)) == 6) |>
  # There are unknown genders
  verify(length(unique(patient_gender_code)) == 3) |>
  # Check disease, symptom, and count vars are not missing
  assert(not_na, rsv:patient_count)

# Save full data
saveRDS(rsv, './tmp/rsv_p2_dat_raw.rds')
remove(rsv)

################################################################################
################################################################################
