################################################################################
# File Name: process_data_p2                                                   #
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
library(ISOweek)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

#############################
# 2. PROCESS INFLUENZA DATA #
#############################

# Load data
flu_dat <- readRDS('./tmp/flu_p2_dat_raw_list.rds')

# Create empty lists to fill
flu_proc_symp <- list()

# Loop through each flu season to format data
for (i in 1:4) {
  # Print the progress and set season data
  print(i)
  flu <- as.data.frame(flu_dat[i])
  flu_sum <- sum(flu$patient_count)
  
  # Create new variables
  flu <- flu |>
    # Replace true/false with 1/0
    mutate(flu = ifelse(flu == 'true', 1, 0),
           fever = ifelse(fever == 'true', 1, 0),
           myalgia = ifelse(myalgia == 'true', 1, 0),
           hypoxemia = ifelse(hypoxemia == 'true', 1, 0),
           short_breath = ifelse(short_breath == 'true', 1, 0),
           cough = ifelse(cough == 'true', 1, 0),
           bronchitis = ifelse(bronchitis == 'true', 1, 0),
           chest_pain = ifelse(chest_pain == 'true', 1, 0),
           nausea_vom = ifelse(nausea_vom == 'true', 1, 0),
           sore_throat = ifelse(sore_throat == 'true', 1, 0),
           fatigue = ifelse(fatigue == 'true', 1, 0),
           diarrhea = ifelse(diarrhea == 'true', 1, 0),
           headache = ifelse(headache == 'true', 1, 0),
           congestion = ifelse(congestion == 'true', 1, 0),
           sneezing = ifelse(sneezing == 'true', 1, 0)) |>
    # Add a leading 0 to county FIPS code
    mutate(county_fips = as.character(county_fips),
           county_fips = ifelse(nchar(county_fips) == 4, paste0("0", county_fips), county_fips)) |>
    # Create state fips
    mutate(state_fips = substr(county_fips, 1, 2)) |>
    # Create symptom count
    mutate(symp_count = fever + myalgia + cough + 
             sore_throat + short_breath + hypoxemia + 
             chest_pain + bronchitis + nausea_vom + 
             diarrhea + fatigue + headache + congestion + sneezing) 
  
  # Drop symptomatic individuals
  flu_symp <- flu |>
    filter(symp_count > 0)
  # Check the number of individuals kept
  print(paste("The percent of patients kept and symptomatic is:",
              sum(flu_symp$patient_count) / flu_sum))
  
  # Remove to save space
  remove(flu)
  
  # Drop based on:
  flu_filt_symp <- flu_symp |>
    # Missing county or not US state & DC
    filter(!(state_fips %in% c("72","99",""))) |>
    # Unknown or missing ages
    filter(!(age_grp %in% c("U",""))) |>
    # Unknown or missing genders
    filter(!(patient_gender_code %in% c("U","")))
  # Check the number of individuals kept
  print(paste("The percent of patients kept is:",
              sum(flu_filt_symp$patient_count) / sum(flu_symp$patient_count)))
  
  # Confirm drop was completed
  flu_filt_symp |>
    # There should be more than 2500 county groups
    verify(state_fips != "99" | state_fips != "72") |>
    # There are unknown age groups
    verify(age_grp != "U") |>
    # There are unknown genders
    verify(patient_gender_code != "") |>
    # Check disease, symptom, and count vars are not missing
    assert(not_na, county_fips:symp_count)
  
  # Fill lists
  flu_proc_symp[[i]] <- as.data.frame(flu_filt_symp)
  remove(flu_symp, flu_filt_symp)

}
  
# Save filtered versions
saveRDS(flu_proc_symp, './tmp/flu_p2_dat_proc_symp.rds')
remove(flu_dat)

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
rsv <- readRDS('./tmp/rsv_p2_dat_raw.rds')

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
  # Convert week date to date
  mutate(week_date = sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", year_week)) |>
  mutate(week_date = ISOweek2date(week_date)) |>
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
saveRDS(rsv_filt, './tmp/rsv_p2_dat_proc.rds')
saveRDS(rsv_filt_symp, './tmp/rsv_p2_dat_proc_symp.rds')
remove(rsv, rsv_filt, rsv_filt_symp)

################################################################################
################################################################################
