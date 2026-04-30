################################################################################
# File Name: 02b_process_rsv_data                                              #
#                                                                              #
# Purpose:   Process RSV medical claims data from 2016-2023.                   #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Read in and process data                                       #
#                 a. Check data                                                #
#                 b. Format variables                                          #
#                 c. Filter data                                               #
#                 d. Save data                                                 #
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
library(arrow)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

###############################
# 2. READ IN AND PROCESS DATA #
###############################

###############
# PART I DATA #
###############

# Create import file list
filelist <- c('tmp/rsv_2016-09_2017-08_1d_imputed.parquet', 
              'tmp/rsv_2017-09_2018-08_1d_imputed.parquet',
              'tmp/rsv_2018-09_2019-08_1d_imputed.parquet',
              'tmp/rsv_2019-09_2020-08_1d_imputed.parquet',
              'tmp/rsv_2020-09_2021-08_1d_imputed.parquet',
              'tmp/rsv_2021-09_2022-08_1d_imputed.parquet',
              'tmp/rsv_2022-09_2023-08_1d_imputed.parquet')

# Loop through each season of data
for (i in 1:7) {
  
  # Load data
  rsv_dat <- read_parquet(filelist[i])
  cat("Loading", filelist[i], "...")
  
  # Convert to all strings
  rsv_dat <- rsv_dat |>
    mutate(across(everything(), as.character))
  
  rsv <- rsv_dat |>
    # Convert and create time variables
    mutate(year_month = paste0(year_month, '-01'),
           year_month_dt = as.Date(year_month),
           year = year(year_month_dt),
           month = month(year_month_dt)) |>
    # Convert patient county numeric
    mutate(patient_count = as.numeric(patient_count)) |>
    # Replace true/false with 1/0
    mutate(rsv = ifelse(rsv == 'TRUE', 1, 0),
           fever = ifelse(fever == 'TRUE', 1, 0),
           loss_appetite = ifelse(loss_appetite == 'TRUE', 1, 0),
           short_breath = ifelse(short_breath == 'TRUE', 1, 0),
           cough = ifelse(cough == 'TRUE', 1, 0),
           bronchitis = ifelse(bronchitis == 'TRUE', 1, 0),
           sore_throat = ifelse(sore_throat == 'TRUE', 1, 0),
           fatigue = ifelse(fatigue == 'TRUE', 1, 0),
           hypoxemia = ifelse(hypoxemia == 'TRUE', 1, 0),
           headache = ifelse(headache == 'TRUE', 1, 0),
           congestion = ifelse(congestion == 'TRUE', 1, 0),
           sneezing = ifelse(sneezing == 'TRUE', 1, 0)) |>
    # Add a leading 0 to county FIPS code
    mutate(county_fips = ifelse(nchar(county_fips) == 4, paste0("0", county_fips), county_fips)) |>
    # Create state fips
    mutate(state_fips = substr(county_fips, 1, 2)) |>
    # Create symptom count
    mutate(symp_count = fever + loss_appetite + cough + 
             sore_throat + short_breath + hypoxemia + 
             bronchitis + fatigue + headache + 
             congestion + sneezing) |>
    # Filter to only symptomatic patients
    dplyr::filter(symp_count > 0) |>
    # Filter based on:
    # Missing county or not US state & DC
    dplyr::filter(!(state_fips %in% c("72","99",""))) |>
    # Unknown or missing ages
    dplyr::filter(!(age_grp %in% c("U",""))) |>
    # Unknown or missing genders
    dplyr::filter(!(patient_gender_code %in% c("U","")))
  
  # Save to list
  filename <- paste0("tmp/rsv_", i, "_proc_p1.parquet")
  write_parquet(rsv, filename)
  cat("Saved:", filename, "\n")
}

################
# PART II DATA #
################

# Create import file list
filelist <- c('tmp/rsv_2016-09_2017-08_p2_imputed.parquet', 
              'tmp/rsv_2017-09_2018-08_p2_imputed.parquet',
              'tmp/rsv_2018-09_2019-08_p2_imputed.parquet',
              'tmp/rsv_2019-09_2020-08_p2_imputed.parquet',
              'tmp/rsv_2020-09_2021-08_p2_imputed.parquet',
              'tmp/rsv_2021-09_2022-08_p2_imputed.parquet',
              'tmp/rsv_2022-09_2023-08_p2_imputed.parquet')

# Loop through each season of data
for (i in 1:7) {

  # Load data
  rsv_dat <- read_parquet(filelist[i])
  cat("Loading", filelist[i], "...")
  
  # Convert to all strings
  rsv_dat <- rsv_dat |>
    mutate(across(everything(), as.character))
  
  rsv <- rsv_dat |>
    # Convert and create time variables
    mutate(year_week_dt = as.Date(year_week_dt),
           year = year(year_week_dt),
           month = month(year_week_dt)) |>
    # Convert patient county numeric
    mutate(patient_count = as.numeric(patient_count)) |>
    # Replace true/false with 1/0
    mutate(rsv = ifelse(rsv == 'TRUE', 1, 0),
           fever = ifelse(fever == 'TRUE', 1, 0),
           loss_appetite = ifelse(loss_appetite == 'TRUE', 1, 0),
           short_breath = ifelse(short_breath == 'TRUE', 1, 0),
           cough = ifelse(cough == 'TRUE', 1, 0),
           bronchitis = ifelse(bronchitis == 'TRUE', 1, 0),
           sore_throat = ifelse(sore_throat == 'TRUE', 1, 0),
           fatigue = ifelse(fatigue == 'TRUE', 1, 0),
           hypoxemia = ifelse(hypoxemia == 'TRUE', 1, 0),
           headache = ifelse(headache == 'TRUE', 1, 0),
           congestion = ifelse(congestion == 'TRUE', 1, 0),
           sneezing = ifelse(sneezing == 'TRUE', 1, 0)) |>
    # Add a leading 0 to county FIPS code
    mutate(county_fips = ifelse(nchar(county_fips) == 4, paste0("0", county_fips), county_fips)) |>
    # Create state fips
    mutate(state_fips = substr(county_fips, 1, 2)) |>
    # Create symptom count
    mutate(symp_count = fever + loss_appetite + cough + 
             sore_throat + short_breath + hypoxemia + 
             bronchitis + fatigue + headache + 
             congestion + sneezing) |>
    # Filter to only symptomatic patients
    dplyr::filter(symp_count > 0) |>
    # Filter based on:
    # Missing county or not US state & DC
    dplyr::filter(!(state_fips %in% c("72","99",""))) |>
    # Unknown or missing ages
    dplyr::filter(!(age_grp %in% c("U",""))) |>
    # Unknown or missing genders
    dplyr::filter(!(patient_gender_code %in% c("U","")))
  
  # Save to list
  filename <- paste0("tmp/rsv_", i, "_proc_p2.parquet")
  write_parquet(rsv, filename)
  cat("Saved:", filename, "\n")
}

################################################################################
################################################################################
