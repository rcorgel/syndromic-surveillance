################################################################################
# File Name: sample_part_2_dat                                                 #
#                                                                              #
# Purpose:   Load and process Influenza medical claims data from 2016-2020.    #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load and check data                                            #
#            3. Process data                                                   #
#            4. Subset data                                                    #
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

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

################
# 2. LOAD DATA #
################

# Load Influenza symptom data
flu <- read.csv('./raw/flu_2016-09_2017-08_p2_imputed.csv')
# Don't include 2020-2023 due to low flu circulation & covid-19 co-circulation

# Run sanity checks on the data
flu |>
  # There should be 52-54 weeks
  verify(length(unique(year_week)) > 51) |>
  verify(length(unique(year_week)) < 55) |>
  # There should be more than 3000 counties
  verify(length(unique(county_fips)) > 3000) |>
  # There are unknown age groups, 6 plus 1
  verify(length(unique(age_grp)) == 7) |>
  # There are unknown genders, 2 plus 1
  verify(length(unique(patient_gender_code)) == 3) |>
  # Check disease, symptom, and count vars are not missing
  assert(not_na, flu:patient_count)

###################
# 3. PROCESS DATA #
###################

# CREATE/ADD NEW VARIABLES #

flu <- flu |>
  # Add time variables
  mutate(month = month(year_week_dt)) |>
  mutate(year = year(year_week_dt)) |>
  # Fix county fips to add a leading 0
  mutate(county_fips = ifelse(nchar(county_fips) == 4, 
                              paste0('0', county_fips), county_fips)) |>
  verify(nchar(county_fips) == 5) |>
  # Create state fips
  mutate(state_fips = substr(county_fips, 1, 2)) |>
  # Convert symptom vars from true/false to 1/0
  mutate(flu = ifelse(flu == 'true', 1, 0),
         fever = ifelse(fever == 'true', 1, 0),
         myalgia = ifelse(myalgia == 'true', 1, 0),
         cough = ifelse(cough == 'true', 1, 0),
         sore_throat = ifelse(sore_throat == 'true', 1, 0),
         short_breath = ifelse(short_breath == 'true', 1, 0),
         hypoxemia = ifelse(hypoxemia == 'true', 1, 0),
         chest_pain = ifelse(chest_pain == 'true', 1, 0),
         bronchitis = ifelse(bronchitis == 'true', 1, 0),
         nausea_vom = ifelse(nausea_vom == 'true', 1, 0),
         diarrhea = ifelse(diarrhea == 'true', 1, 0),
         fatigue = ifelse(fatigue == 'true', 1, 0),
         headache = ifelse(headache == 'true', 1, 0),
         congestion = ifelse(congestion == 'true', 1, 0),
         sneezing = ifelse(sneezing == 'true', 1, 0),) |>
  # Create symptom count
  mutate(symp_count = fever + myalgia + cough + 
           sore_throat + short_breath + hypoxemia + 
           chest_pain + bronchitis + nausea_vom + 
           diarrhea + fatigue + headache + congestion + sneezing)

# Merge on geographic data
urban_south <-readRDS("./tmp/urban_south_binary.rds")
flu <- left_join(flu, urban_south, by = 'county_fips')

# FILTER DATA #

# Drop based on:
flu_filt <- flu |>
  # Missing county or not US state & DC
  dplyr::filter(!(state_fips %in% c("72","99",""))) |>
  # Unknown or missing ages
  dplyr::filter(!(age_grp %in% c("U",""))) |>
  # Unknown or missing genders
  dplyr::filter(!(patient_gender_code %in% c("U","")))

# Check the number of individuals kept
print(paste("The percent of patients kept is:",
            sum(flu_filt$patient_count) / sum(flu$patient_count)))

# Confirm drop was completed
flu_filt |>
  # There should be no unknown states or PR
  verify(state_fips != "99" | state_fips != "72") |>
  # There are unknown age groups
  verify(age_grp != "U") |>
  # There are unknown genders
  verify(patient_gender_code != "U") |>
  # Check disease, symptom, and count vars are not missing
  assert(not_na, county_fips:symp_count)

# Drop symptomatic individuals
flu_filt_symp <- flu_filt |>
  dplyr::filter(symp_count > 0)

# Check the number of individuals kept
print(paste("The percent of patients kept and symptomatic is:",
            sum(flu_filt_symp$patient_count) / sum(flu$patient_count)))

# SAVE DATA #
saveRDS(flu_filt_symp, './tmp/flu_filt_symp_p2_2016.rds')

##################
# 4. SUBSET DATA #
##################

# Remove some of the previous data
remove(flu, flu_filt, urban_south)

# Sample 250,000 individuals from the symptomatic data
flu_balanced <- flu_filt_symp |> group_by(flu) |>
  uncount(patient_count) |> sample_n(125000)

# Save balanced data
saveRDS(flu_balanced, './tmp/flu_filt_symp_p2_bal_2016.rds')

################################################################################
################################################################################
