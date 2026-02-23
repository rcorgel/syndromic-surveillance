################################################################################
# File Name: process_data_p1                                                   #
#                                                                              #
# Purpose:   Process Influenza, COVID-19, and RSV medical claims data.         #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Process Influenza data                                         #
#            3. Process COVID-19 data                                          #
#            4. Process RSV data                                               #
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
library(sf)
library(spdep)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

#############################
# 2. PROCESS INFLUENZA DATA #
#############################

# Load weather data to merge onto all disease data
weather <- readRDS('./tmp/county_month_weather.rds')

# Load data
flu <- readRDS('./tmp/flu_p1_dat_raw.rds')

# Impute data
# Check the percent of rows that will be imputed (81.2%)
nrow(flu[flu$patient_count == "<=5",]) / nrow(flu)
# Impute rows where patient count is <=5 to a random number 1 to 5
flu$patient_count_imp <- flu$patient_count
flu$patient_count_imp[flu$patient_count_imp == '<=5'] <- 
  sample(1:5, length(flu$patient_count_imp[flu$patient_count_imp == '<=5']), 
         replace= TRUE)
# Convert imputed count to numeric
flu$patient_count_imp <- as.numeric(flu$patient_count_imp)
# Check the percent of patients that were imputed (8.7%)
sum(flu[flu$patient_count_imp <= 5,]$patient_count_imp) / sum(flu$patient_count_imp)

# Create new variables
flu <- flu |>
  # Convert month date to date
  mutate(month_date = as.Date(paste(year_month, "-01", sep=""))) |>
  # Create state fips
  mutate(state_fips = substr(county_fips, 1, 2)) |>
  # Create symptom count
  mutate(symp_count = fever + myalgia + cough + 
           sore_throat + short_breath + hypoxemia + 
           chest_pain + bronchitis + nausea_vom + 
           diarrhea + fatigue + headache + congestion + sneezing)

# Check the percent of patients that will be dropped for:
# Missing county or not US state & DC (0.4%)
sum(flu[flu$state_fips %in% c("72","99",""),]$patient_count_imp) / sum(flu$patient_count_imp)
# Unknown or missing ages (1.4%)
sum(flu[flu$age_grp %in% c("U",""),]$patient_count_imp) / sum(flu$patient_count_imp)
# Unknown or missing genders (<0.1%)
sum(flu[flu$patient_gender_code %in% c("U",""),]$patient_count_imp) / sum(flu$patient_count_imp)
# Asymptomatic or non-coded individuals (83.4%)
sum(flu[flu$symp_count == 0,]$patient_count_imp) / sum(flu$patient_count_imp)

# Drop based on:
flu_filt <- flu |>
  # Missing county or not US state & DC
  filter(!(state_fips %in% c("72","99",""))) |>
  # Unknown or missing ages
  filter(!(age_grp %in% c("U",""))) |>
  # Unknown or missing genders
  filter(!(patient_gender_code %in% c("U","")))
# Check the number of individuals kept (98.2%)
sum(flu_filt$patient_count_imp) / sum(flu$patient_count_imp)

# Confirm drop was completed
flu_filt |>
  # There should be more than 2500 county groups
  verify(state_fips != "99" | state_fips != "72") |>
  # There are unknown age groups
  verify(age_grp != "U") |>
  # There are unknown genders
  verify(patient_gender_code != "U") |>
  # Check disease, symptom, and count vars are not missing
  assert(not_na, county_fips:symp_count)

# Merge on weather data
# weather$month_date <- paste(weather$year, weather$month, sep="-")
# weather$month_date <- as.Date(paste(weather$month_date, "-01", sep=""), '%Y-%m-%d')
# flu_filt <- left_join(flu_filt, weather, by = c('county_fips' = 'county_group',
#                                                 'month_date' = 'month_date'))

# Create lagged fever, sore throat, and cough percent variables
lagged_vars <- flu_filt |> group_by(county_fips, month_date) |>
  mutate(fever_perc_lag = sum(fever*patient_count_imp) / sum(patient_count_imp),
         cough_perc_lag = sum(cough*patient_count_imp) / sum(patient_count_imp),
         sore_throat_perc_lag = sum(sore_throat*patient_count_imp) / sum(patient_count_imp)) |>
  distinct(county_fips, month_date, fever_perc_lag, cough_perc_lag, sore_throat_perc_lag, .keep_all = FALSE) |>
  mutate(lag_month_date = month_date %m+% months(1)) |> ungroup() |>
  select(-c(month_date))

# Merge on lagged variables
flu_filt <- left_join(flu_filt, lagged_vars, by = c('county_fips' = 'county_fips',
                                                'month_date' = 'lag_month_date'))

usa_albers_state <- st_as_sf(readRDS('./tmp/usa_albers_state.rds'))   # convert to sf
usa_albers_county <- st_as_sf(readRDS('./tmp/usa_albers_county.rds')) # convert to sf

#neighbors <- poly2nb(usa_albers_county , queen = TRUE)
#usa_albers_county$neighbor_list <- sapply(neighbors, function(n){c(usa_albers_county$GEOID[n])})

# Collapse county map to county group
usa_albers_county_x_walk <- left_join(x_walk, usa_albers_county[, c(5, 6)], by = c('fips' = 'GEOID'))
usa_albers_county_group <- usa_albers_county_x_walk %>% 
  group_by(county_fips) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()

# Drop symptomatic individuals
flu_filt_symp <- flu_filt |>
  filter(symp_count > 0)
# Check the number of individuals kept (16.3%)
sum(flu_filt_symp$patient_count_imp) / sum(flu$patient_count_imp)

# Save filtered versions
saveRDS(flu_filt, './tmp/flu_p1_dat_proc.rds')
saveRDS(flu_filt_symp, './tmp/flu_p1_dat_proc_symp.rds')
remove(flu, flu_filt, flu_filt_symp)

#########################
# 3. PROCESS COVID DATA #
#########################

# Load data
covid <- readRDS('./tmp/covid_p1_dat_raw.rds')

# Impute data
# Check the percent of rows that will be imputed (82.9%)
nrow(covid[covid$patient_count == "<=5",]) / nrow(covid)
# Impute rows where patient count is <=5 to a random number 1 to 5
covid$patient_count_imp <- covid$patient_count
covid$patient_count_imp[covid$patient_count_imp == '<=5'] <- 
  sample(1:5, length(covid$patient_count_imp[covid$patient_count_imp == '<=5']), 
         replace= TRUE)
# Convert imputed count to numeric
covid$patient_count_imp <- as.numeric(covid$patient_count_imp)
# Check the percent of patients that were imputed (9.5%)
sum(covid[covid$patient_count_imp <= 5,]$patient_count_imp) / sum(covid$patient_count_imp)

# Create new variables
covid <- covid |>
  # Convert month date to date
  mutate(month_date = as.Date(paste(year_month, "-01", sep=""))) |>
  # Create state fips
  mutate(state_fips = substr(county_fips, 1, 2)) |>
  # Create symptom count
  mutate(symp_count = fever + myalgia + cough + 
           sore_throat + short_breath + hypoxemia + 
           chest_pain + bronchitis + nausea_vom + loss_appetite +
           diarrhea + fatigue + headache + congestion + loss_smell_taste)

# Check the percent of patients that will be dropped for:
# Missing county or not US state & DC (0.2%)
sum(covid[covid$state_fips %in% c("72","99",""),]$patient_count_imp) / sum(covid$patient_count_imp)
# Unknown or missing ages (1.9%)
sum(covid[covid$age_grp %in% c("U",""),]$patient_count_imp) / sum(covid$patient_count_imp)
# Unknown or missing genders (<0.1%)
sum(covid[covid$patient_gender_code %in% c("U",""),]$patient_count_imp) / sum(covid$patient_count_imp)
# Asymptomatic or non-coded individuals (83.9%)
sum(covid[covid$symp_count == 0,]$patient_count_imp) / sum(covid$patient_count_imp)

# Drop based on:
covid_filt <- covid |>
  # Missing county or not US state & DC
  filter(!(state_fips %in% c("72","99",""))) |>
  # Unknown or missing ages
  filter(!(age_grp %in% c("U",""))) |>
  # Unknown or missing genders
  filter(!(patient_gender_code %in% c("U","")))
# Check the number of individuals kept (97.9%)
sum(covid_filt$patient_count_imp) / sum(covid$patient_count_imp)

# Confirm drop was completed
covid_filt |>
  # There should be more than 2500 county groups
  verify(state_fips != "99" | state_fips != "72") |>
  # There are unknown age groups
  verify(age_grp != "U") |>
  # There are unknown genders
  verify(patient_gender_code != "U") |>
  # Check disease, symptom, and count vars are not missing
  assert(not_na, county_fips:symp_count)

# Drop symptomatic individuals
covid_filt_symp <- covid_filt |>
  filter(symp_count > 0)
# Check the number of individuals kept (15.7%)
sum(covid_filt_symp$patient_count_imp) / sum(covid$patient_count_imp)

# Save filtered versions
saveRDS(covid_filt, './tmp/covid_p1_dat_proc.rds')
saveRDS(covid_filt_symp, './tmp/covid_p1_dat_proc_symp.rds')
remove(covid, covid_filt, covid_filt_symp)

#######################
# 4. PROCESS RSV DATA #
#######################

# Load data
rsv <- readRDS('./tmp/rsv_p1_dat_raw.rds')

# Impute data
# Check the percent of rows that will be imputed (76.0%)
nrow(rsv[rsv$patient_count == "<=5",]) / nrow(rsv)
# Impute rows where patient count is <=5 to a random number 1 to 5
rsv$patient_count_imp <- rsv$patient_count
rsv$patient_count_imp[rsv$patient_count_imp == '<=5'] <- 
  sample(1:5, length(rsv$patient_count_imp[rsv$patient_count_imp == '<=5']), 
         replace= TRUE)
# Convert imputed count to numeric
rsv$patient_count_imp <- as.numeric(rsv$patient_count_imp)
# Check the percent of patients that were imputed (4.7%)
sum(rsv[rsv$patient_count_imp <= 5,]$patient_count_imp) / sum(rsv$patient_count_imp)

# Create new variables
rsv <- rsv |>
  # Convert month date to date
  mutate(month_date = as.Date(paste(year_month, "-01", sep=""))) |>
  # Create state fips
  mutate(state_fips = substr(county_fips, 1, 2)) |>
  # Create symptom count
  mutate(symp_count = fever + cough + sore_throat +
           short_breath + hypoxemia + 
           bronchitis + loss_appetite +
           fatigue + headache + congestion + sneezing)

# Check the percent of patients that will be dropped for:
# Missing county or not US state & DC (0.3%)
sum(rsv[rsv$state_fips %in% c("72","99",""),]$patient_count_imp) / sum(rsv$patient_count_imp)
# Unknown or missing ages (1.4%)
sum(rsv[rsv$age_grp %in% c("U",""),]$patient_count_imp) / sum(rsv$patient_count_imp)
# Unknown or missing genders (<0.1%)
sum(rsv[rsv$patient_gender_code %in% c("U",""),]$patient_count_imp) / sum(rsv$patient_count_imp)
# Asymptomatic or non-coded individuals (88.7%)
sum(rsv[rsv$symp_count == 0,]$patient_count_imp) / sum(rsv$patient_count_imp)

# Drop based on:
rsv_filt <- rsv |>
  # Missing county or not US state & DC
  filter(!(state_fips %in% c("72","99",""))) |>
  # Unknown or missing ages
  filter(!(age_grp %in% c("U",""))) |>
  # Unknown or missing genders
  filter(!(patient_gender_code %in% c("U","")))
# Check the number of individuals kept (98.2%)
sum(rsv_filt$patient_count_imp) / sum(rsv$patient_count_imp)

# Confirm drop was completed
rsv_filt |>
  # There should be more than 2500 county groups
  verify(state_fips != "99" | state_fips != "72") |>
  # There are unknown age groups
  verify(age_grp != "U") |>
  # There are unknown genders
  verify(patient_gender_code != "U") |>
  # Check disease, symptom, and count vars are not missing
  assert(not_na, county_fips:symp_count)

# Drop symptomatic individuals
rsv_filt_symp <- rsv_filt |>
  filter(symp_count > 0)
# Check the number of individuals kept (11.1%)
sum(rsv_filt_symp$patient_count_imp) / sum(rsv$patient_count_imp)

# Save filtered versions
saveRDS(rsv_filt, './tmp/rsv_p1_dat_proc.rds')
saveRDS(rsv_filt_symp, './tmp/rsv_p1_dat_proc_symp.rds')
remove(rsv, rsv_filt, rsv_filt_symp)

################################################################################
################################################################################
