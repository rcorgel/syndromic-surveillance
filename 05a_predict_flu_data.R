################################################################################
# File Name: 05a_predict_flu_data                                              #
#                                                                              #
# Purpose:   Predict influenza from medical claims symptoms.                   #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load models                                                    #
#            3. Predict data                                                   #
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
library(cowplot)
library(lubridate)
library(arrow)
library(mgcv)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

##################
# 2. LOAD MODELS #
##################

# Load V1 Models
fever_model <- readRDS("./tmp/fever_model_v1_flu.rds")
all_symptom_model  <- readRDS("./tmp/all_symptom_model_v1_flu.rds")
symptom_interaction_model  <- readRDS("./tmp/symptom_interaction_model_v1_flu.rds")
interaction_case_model  <- readRDS("./tmp/interaction_case_model_v1_flu.rds")
interaction_nrevss_model  <- readRDS("./tmp/interaction_nrevss_model_v1_flu.rds")

# Cutpoints
cutpoints_v1 <- readRDS("./tmp/model_cutpoints_v1_flu.rds")
cutpoints_v2 <- readRDS("./tmp/model_cutpoints_v2_flu.rds")

# Case per visit
cases_visits <- readRDS("./tmp/visit_positivity_county_flu.rds")

# Load nrevss data
# Load NREVSS data
nrevss_state <- readRDS('./tmp/nrevss_state.rds')
nrevss_nat <- readRDS('./tmp/nrevss_national.rds')

# Set prediction function
predict_flu <- function(df, cutoff_df) {
  
  # Made predictions based off the four models
  df$pred_fever <- predict(fever_model, newdata = df, type = "response")
  df$pred_symptom <- predict(all_symptom_model, newdata = df, type = "response")
  df$pred_interaction <- predict(symptom_interaction_model, newdata = df, type = "response")
  df$pred_case <- predict(interaction_case_model, newdata = df, type = "response")
  df$pred_nrevss <- predict(interaction_nrevss_model, newdata = df, type = "response")
  
  # Cut predictions into 0/1 based off the cutoff data
  result <- df |>
    mutate(pred_fever_count = ifelse(pred_fever >= cutoff_df[cutoff_df$age_grp == '0-4',]$cutoff_fever 
                                     & age_grp == '0-4', 1, NA),
           pred_fever_count = ifelse(pred_fever >= cutoff_df[cutoff_df$age_grp == '5-12',]$cutoff_fever 
                                     & age_grp == '5-12', 1, pred_fever_count),
           pred_fever_count = ifelse(pred_fever >= cutoff_df[cutoff_df$age_grp == '13-17',]$cutoff_fever 
                                     & age_grp == '13-17', 1, pred_fever_count),
           pred_fever_count = ifelse(pred_fever >= cutoff_df[cutoff_df$age_grp == '18-49',]$cutoff_fever 
                                     & age_grp == '18-49', 1, pred_fever_count),
           pred_fever_count = ifelse(pred_fever >= cutoff_df[cutoff_df$age_grp == '50-64',]$cutoff_fever 
                                     & age_grp == '50-64', 1, pred_fever_count),
           pred_fever_count = ifelse(pred_fever >= cutoff_df[cutoff_df$age_grp == '>=65',]$cutoff_fever 
                                     & age_grp == '>=65', 1, pred_fever_count),
           pred_fever_count = ifelse(is.na(pred_fever_count), 0, pred_fever_count),
           pred_fever_count = pred_fever_count * patient_count) |>
    mutate(pred_symptom_count = ifelse(pred_symptom >= cutoff_df[cutoff_df$age_grp == '0-4',]$cutoff_symptom 
                                       & age_grp == '0-4', 1, NA),
           pred_symptom_count = ifelse(pred_symptom >= cutoff_df[cutoff_df$age_grp == '5-12',]$cutoff_symptom 
                                       & age_grp == '5-12', 1, pred_symptom_count),
           pred_symptom_count = ifelse(pred_symptom >= cutoff_df[cutoff_df$age_grp == '13-17',]$cutoff_symptom 
                                       & age_grp == '13-17', 1, pred_symptom_count),
           pred_symptom_count = ifelse(pred_symptom >= cutoff_df[cutoff_df$age_grp == '18-49',]$cutoff_symptom 
                                       & age_grp == '18-49', 1, pred_symptom_count),
           pred_symptom_count = ifelse(pred_symptom >= cutoff_df[cutoff_df$age_grp == '50-64',]$cutoff_symptom 
                                       & age_grp == '50-64', 1, pred_symptom_count),
           pred_symptom_count = ifelse(pred_symptom >= cutoff_df[cutoff_df$age_grp == '>=65',]$cutoff_symptom 
                                       & age_grp == '>=65', 1, pred_symptom_count),
           pred_symptom_count = ifelse(is.na(pred_symptom_count), 0, pred_symptom_count),
           pred_symptom_count = pred_symptom_count * patient_count) |>
    mutate(pred_interaction_count = ifelse(pred_interaction >= cutoff_df[cutoff_df$age_grp == '0-4',]$cutoff_interaction 
                                           & age_grp == '0-4', 1, NA),
           pred_interaction_count = ifelse(pred_interaction >= cutoff_df[cutoff_df$age_grp == '5-12',]$cutoff_interaction 
                                           & age_grp == '5-12', 1, pred_interaction_count),
           pred_interaction_count = ifelse(pred_interaction >= cutoff_df[cutoff_df$age_grp == '13-17',]$cutoff_interaction 
                                           & age_grp == '13-17', 1, pred_interaction_count),
           pred_interaction_count = ifelse(pred_interaction >= cutoff_df[cutoff_df$age_grp == '18-49',]$cutoff_interaction 
                                           & age_grp == '18-49', 1, pred_interaction_count),
           pred_interaction_count = ifelse(pred_interaction >= cutoff_df[cutoff_df$age_grp == '50-64',]$cutoff_interaction 
                                           & age_grp == '50-64', 1, pred_interaction_count),
           pred_interaction_count = ifelse(pred_interaction >= cutoff_df[cutoff_df$age_grp == '>=65',]$cutoff_interaction 
                                           & age_grp == '>=65', 1, pred_interaction_count),
           pred_interaction_count = ifelse(is.na(pred_interaction_count), 0, pred_interaction_count),
           pred_interaction_count = pred_interaction_count * patient_count) |>
    mutate(pred_case_count = ifelse(pred_case >= cutoff_df[cutoff_df$age_grp == '0-4',]$cutoff_case 
                                    & age_grp == '0-4', 1, NA),
           pred_case_count = ifelse(pred_case >= cutoff_df[cutoff_df$age_grp == '5-12',]$cutoff_case 
                                    & age_grp == '5-12', 1, pred_case_count),
           pred_case_count = ifelse(pred_case >= cutoff_df[cutoff_df$age_grp == '13-17',]$cutoff_case 
                                    & age_grp == '13-17', 1, pred_case_count),
           pred_case_count = ifelse(pred_case >= cutoff_df[cutoff_df$age_grp == '18-49',]$cutoff_case 
                                    & age_grp == '18-49', 1, pred_case_count),
           pred_case_count = ifelse(pred_case >= cutoff_df[cutoff_df$age_grp == '50-64',]$cutoff_case 
                                    & age_grp == '50-64', 1, pred_case_count),
           pred_case_count = ifelse(pred_case >= cutoff_df[cutoff_df$age_grp == '>=65',]$cutoff_case 
                                    & age_grp == '>=65', 1, pred_case_count),
           pred_case_count = ifelse(is.na(pred_case_count), 0, pred_case_count),
           pred_case_count = pred_case_count * patient_count) |>
    mutate(pred_nrevss_count = ifelse(pred_nrevss >= cutoff_df[cutoff_df$age_grp == '0-4',]$cutoff_nrevss 
                                    & age_grp == '0-4', 1, NA),
           pred_nrevss_count = ifelse(pred_nrevss >= cutoff_df[cutoff_df$age_grp == '5-12',]$cutoff_nrevss 
                                    & age_grp == '5-12', 1, pred_nrevss_count),
           pred_nrevss_count = ifelse(pred_nrevss >= cutoff_df[cutoff_df$age_grp == '13-17',]$cutoff_nrevss 
                                    & age_grp == '13-17', 1, pred_nrevss_count),
           pred_nrevss_count = ifelse(pred_nrevss >= cutoff_df[cutoff_df$age_grp == '18-49',]$cutoff_nrevss 
                                    & age_grp == '18-49', 1, pred_nrevss_count),
           pred_nrevss_count = ifelse(pred_nrevss >= cutoff_df[cutoff_df$age_grp == '50-64',]$cutoff_nrevss 
                                    & age_grp == '50-64', 1, pred_nrevss_count),
           pred_nrevss_count = ifelse(pred_nrevss >= cutoff_df[cutoff_df$age_grp == '>=65',]$cutoff_nrevss 
                                    & age_grp == '>=65', 1, pred_nrevss_count),
           pred_nrevss_count = ifelse(is.na(pred_nrevss_count), 0, pred_nrevss_count),
           pred_nrevss_count = pred_nrevss_count * patient_count) |>
    # Collapse observations and predictions by week, county, and age group
    group_by(year_week_dt, county_fips, age_grp) |>
    mutate(flu_count = sum(flu * patient_count),
           pred_fever_count = sum(pred_fever_count),
           pred_symptom_count = sum(pred_symptom_count),
           pred_interaction_count = sum(pred_interaction_count),
           pred_case_count = sum(pred_case_count),
           pred_nrevss_count = sum(pred_nrevss_count)) |>
    distinct(year_week_dt, county_fips, age_grp, flu_count, 
           pred_fever_count, pred_symptom_count,
           pred_interaction_count, pred_case_count, pred_nrevss_count)
  
  # Output result
  return(result)
}

###################
# 3. PREDICT DATA #
###################

########
# 2016 #
########

# Load 2016 data
flu_2016 <- read_parquet('tmp/flu_1_proc_p2.parquet')

# Megge on visit data
flu_2016 <- left_join(flu_2016, cases_visits[, c(1, 2, 6, 7)], by = c('county_fips', 'year_week_dt'))

# Merge on NREVSS data
flu_2016 <- left_join(flu_2016, nrevss_state[, c(3, 4, 5)], by = c('year_week_dt' = 'week_date',
                                                               'state_fips' = 'state_fips'))
flu_2016 <- left_join(flu_2016, nrevss_nat[, c(3, 5)], by = c('year_week_dt' = 'week_date'))

# Replace missing state values with national
flu_2016$nrevss_pos <- ifelse(is.na(flu_2016$value_roll.x), flu_2016$value_roll.y, flu_2016$value_roll.x)

# Split data by week
flu_2016_list <- split(flu_2016, flu_2016$year_week_dt)
remove(flu_2016)

# Predict data
flu_2016_list <- lapply(flu_2016_list, predict_flu, cutoff_df = cutpoints_v1)

# Combine data
flu_2016 <- do.call(rbind, flu_2016_list)
remove(flu_2016_list)

# Save data
saveRDS(flu_2016, 'tmp/flu_2016_predictions_v1.rds')

########
# 2017 #
########

# Load 2017 data
flu_2017 <- read_parquet('tmp/flu_2_proc_p2.parquet')

# Megge on visit data
flu_2017 <- left_join(flu_2017, cases_visits[, c(1, 2, 6, 7)], by = c('county_fips', 'year_week_dt'))

# Merge on NREVSS data
flu_2017 <- left_join(flu_2017, nrevss_state[, c(3, 4, 5)], by = c('year_week_dt' = 'week_date',
                                                                   'state_fips' = 'state_fips'))
flu_2017 <- left_join(flu_2017, nrevss_nat[, c(3, 5)], by = c('year_week_dt' = 'week_date'))

# Replace missing state values with national
flu_2017$nrevss_pos <- ifelse(is.na(flu_2017$value_roll.x), flu_2017$value_roll.y, flu_2017$value_roll.x)

# Split data by week
flu_2017_list <- split(flu_2017, flu_2017$year_week_dt)
remove(flu_2017)

# Predict data
flu_2017_list <- lapply(flu_2017_list, predict_flu, cutoff_df = cutpoints_v1)

# Combine data
flu_2017 <- do.call(rbind, flu_2017_list)
remove(flu_2017_list)

# Save data
saveRDS(flu_2017, 'tmp/flu_2017_predictions_v1.rds')

########
# 2018 #
########

# Load 2018 data
flu_2018 <- read_parquet('tmp/flu_3_proc_p2.parquet')

# Megge on visit data
flu_2018 <- left_join(flu_2018, cases_visits[, c(1, 2, 6, 7)], by = c('county_fips', 'year_week_dt'))

# Merge on NREVSS data
flu_2018 <- left_join(flu_2018, nrevss_state[, c(3, 4, 5)], by = c('year_week_dt' = 'week_date',
                                                                   'state_fips' = 'state_fips'))
flu_2018 <- left_join(flu_2018, nrevss_nat[, c(3, 5)], by = c('year_week_dt' = 'week_date'))

# Replace missing state values with national
flu_2018$nrevss_pos <- ifelse(is.na(flu_2018$value_roll.x), flu_2018$value_roll.y, flu_2018$value_roll.x)

# Split data by week
flu_2018_list <- split(flu_2018, flu_2018$year_week_dt)
remove(flu_2018)

# Predict data
flu_2018_list <- lapply(flu_2018_list, predict_flu, cutoff_df = cutpoints_v1)

# Combine data
flu_2018 <- do.call(rbind, flu_2018_list)
remove(flu_2018_list)

# Save data
saveRDS(flu_2018, 'tmp/flu_2018_predictions_v1.rds')

########
# 2019 #
########

# Load 2019 data
flu_2019 <- read_parquet('tmp/flu_4_proc_p2.parquet')

# Megge on visit data
flu_2019 <- left_join(flu_2019, cases_visits[, c(1, 2, 6, 7)], by = c('county_fips', 'year_week_dt'))

# Merge on NREVSS data
flu_2019 <- left_join(flu_2019, nrevss_state[, c(3, 4, 5)], by = c('year_week_dt' = 'week_date',
                                                                   'state_fips' = 'state_fips'))
flu_2019 <- left_join(flu_2019, nrevss_nat[, c(3, 5)], by = c('year_week_dt' = 'week_date'))

# Replace missing state values with national
flu_2019$nrevss_pos <- ifelse(is.na(flu_2019$value_roll.x), flu_2019$value_roll.y, flu_2019$value_roll.x)

# Split data by week
flu_2019_list <- split(flu_2019, flu_2019$year_week_dt)
remove(flu_2019)

# Predict data
flu_2019_list <- lapply(flu_2019_list, predict_flu, cutoff_df = cutpoints_v1)

# Combine data
flu_2019 <- do.call(rbind, flu_2019_list)
remove(flu_2019_list)

# Save data
saveRDS(flu_2019, 'tmp/flu_2019_predictions_v1.rds')

########
# 2020 #
########

# Load 2020 data
flu_2020 <- read_parquet('tmp/flu_5_proc_p2.parquet')

# Megge on visit data
flu_2020 <- left_join(flu_2020, cases_visits[, c(1, 2, 6, 7)], by = c('county_fips', 'year_week_dt'))

# Merge on NREVSS data
flu_2020 <- left_join(flu_2020, nrevss_state[, c(3, 4, 5)], by = c('year_week_dt' = 'week_date',
                                                                   'state_fips' = 'state_fips'))
flu_2020 <- left_join(flu_2020, nrevss_nat[, c(3, 5)], by = c('year_week_dt' = 'week_date'))

# Replace missing state values with national
flu_2020$nrevss_pos <- ifelse(is.na(flu_2020$value_roll.x), flu_2020$value_roll.y, flu_2020$value_roll.x)

# Split data by week
flu_2020_list <- split(flu_2020, flu_2020$year_week_dt)
remove(flu_2020)

# Predict data
flu_2020_list <- lapply(flu_2020_list, predict_flu, cutoff_df = cutpoints_v1)

# Combine data
flu_2020 <- do.call(rbind, flu_2020_list)
remove(flu_2020_list)

# Save data
saveRDS(flu_2020, 'tmp/flu_2020_predictions_v1.rds')

########
# 2021 #
########

# Load 2021 data
flu_2021 <- read_parquet('tmp/flu_6_proc_p2.parquet')

# Megge on visit data
flu_2021 <- left_join(flu_2021, cases_visits[, c(1, 2, 6, 7)], by = c('county_fips', 'year_week_dt'))

# Merge on NREVSS data
flu_2021 <- left_join(flu_2021, nrevss_state[, c(3, 4, 5)], by = c('year_week_dt' = 'week_date',
                                                                   'state_fips' = 'state_fips'))
flu_2021 <- left_join(flu_2021, nrevss_nat[, c(3, 5)], by = c('year_week_dt' = 'week_date'))

# Replace missing state values with national
flu_2021$nrevss_pos <- ifelse(is.na(flu_2021$value_roll.x), flu_2021$value_roll.y, flu_2021$value_roll.x)

# Split data by week
flu_2021_list <- split(flu_2021, flu_2021$year_week_dt)
remove(flu_2021)

# Predict data
flu_2021_list <- lapply(flu_2021_list, predict_flu, cutoff_df = cutpoints_v1)

# Combine data
flu_2021 <- do.call(rbind, flu_2021_list)
remove(flu_2021_list)

# Save data
saveRDS(flu_2021, 'tmp/flu_2021_predictions_v1.rds')

########
# 2022 #
########

# Load 2022 data
flu_2022 <- read_parquet('tmp/flu_7_proc_p2.parquet')

# Megge on visit data
flu_2022 <- left_join(flu_2022, cases_visits[, c(1, 2, 6, 7)], by = c('county_fips', 'year_week_dt'))

# Merge on NREVSS data
flu_2022 <- left_join(flu_2022, nrevss_state[, c(3, 4, 5)], by = c('year_week_dt' = 'week_date',
                                                                   'state_fips' = 'state_fips'))
flu_2022 <- left_join(flu_2022, nrevss_nat[, c(3, 5)], by = c('year_week_dt' = 'week_date'))

# Replace missing state values with national
flu_2022$nrevss_pos <- ifelse(is.na(flu_2022$value_roll.x), flu_2022$value_roll.y, flu_2022$value_roll.x)

# Split data by week
flu_2022_list <- split(flu_2022, flu_2022$year_week_dt)
remove(flu_2022)

# Predict data
flu_2022_list <- lapply(flu_2022_list, predict_flu, cutoff_df = cutpoints_v1)

# Combine data
flu_2022 <- do.call(rbind, flu_2022_list)
remove(flu_2022_list)

# Save data
saveRDS(flu_2022, 'tmp/flu_2022_predictions_v1.rds')

## Load V2 Models ##
fever_model <- readRDS("./tmp/fever_model_v2_flu.rds")
all_symptom_model  <- readRDS("./tmp/all_symptom_model_v2_flu.rds")
symptom_interaction_model  <- readRDS("./tmp/symptom_interaction_model_v2_flu.rds")
interaction_case_model  <- readRDS("./tmp/interaction_case_model_v2_flu.rds")
interaction_nrevss_model  <- readRDS("./tmp/interaction_nrevss_model_v2_flu.rds")

########
# 2020 #
########

# Load 2020 data
flu_2020 <- read_parquet('tmp/flu_5_proc_p2.parquet')

# Megge on visit data
flu_2020 <- left_join(flu_2020, cases_visits[, c(1, 2, 6, 7)], by = c('county_fips', 'year_week_dt'))

# Merge on NREVSS data
flu_2020 <- left_join(flu_2020, nrevss_state[, c(3, 4, 5)], by = c('year_week_dt' = 'week_date',
                                                                   'state_fips' = 'state_fips'))
flu_2020 <- left_join(flu_2020, nrevss_nat[, c(3, 5)], by = c('year_week_dt' = 'week_date'))

# Replace missing state values with national
flu_2020$nrevss_pos <- ifelse(is.na(flu_2020$value_roll.x), flu_2020$value_roll.y, flu_2020$value_roll.x)

# Split data by week
flu_2020_list <- split(flu_2020, flu_2020$year_week_dt)
remove(flu_2020)

# Predict data
flu_2020_list <- lapply(flu_2020_list, predict_flu, cutoff_df = cutpoints_v2)

# Combine data
flu_2020 <- do.call(rbind, flu_2020_list)
remove(flu_2020_list)

# Save data
saveRDS(flu_2020, 'tmp/flu_2020_predictions_v2.rds')

########
# 2021 #
########

# Load 2021 data
flu_2021 <- read_parquet('tmp/flu_6_proc_p2.parquet')

# Megge on visit data
flu_2021 <- left_join(flu_2021, cases_visits[, c(1, 2, 6, 7)], by = c('county_fips', 'year_week_dt'))

# Merge on NREVSS data
flu_2021 <- left_join(flu_2021, nrevss_state[, c(3, 4, 5)], by = c('year_week_dt' = 'week_date',
                                                                   'state_fips' = 'state_fips'))
flu_2021 <- left_join(flu_2021, nrevss_nat[, c(3, 5)], by = c('year_week_dt' = 'week_date'))

# Replace missing state values with national
flu_2021$nrevss_pos <- ifelse(is.na(flu_2021$value_roll.x), flu_2021$value_roll.y, flu_2021$value_roll.x)

# Split data by week
flu_2021_list <- split(flu_2021, flu_2021$year_week_dt)
remove(flu_2021)

# Predict data
flu_2021_list <- lapply(flu_2021_list, predict_flu, cutoff_df = cutpoints_v2)

# Combine data
flu_2021 <- do.call(rbind, flu_2021_list)
remove(flu_2021_list)

# Save data
saveRDS(flu_2021, 'tmp/flu_2021_predictions_v2.rds')

########
# 2022 #
########

# Load 2022 data
flu_2022 <- read_parquet('tmp/flu_7_proc_p2.parquet')

# Megge on visit data
flu_2022 <- left_join(flu_2022, cases_visits[, c(1, 2, 6, 7)], by = c('county_fips', 'year_week_dt'))

# Merge on NREVSS data
flu_2022 <- left_join(flu_2022, nrevss_state[, c(3, 4, 5)], by = c('year_week_dt' = 'week_date',
                                                                   'state_fips' = 'state_fips'))
flu_2022 <- left_join(flu_2022, nrevss_nat[, c(3, 5)], by = c('year_week_dt' = 'week_date'))

# Replace missing state values with national
flu_2022$nrevss_pos <- ifelse(is.na(flu_2022$value_roll.x), flu_2022$value_roll.y, flu_2022$value_roll.x)

# Split data by week
flu_2022_list <- split(flu_2022, flu_2022$year_week_dt)
remove(flu_2022)

# Predict data
flu_2022_list <- lapply(flu_2022_list, predict_flu, cutoff_df = cutpoints_v2)

# Combine data
flu_2022 <- do.call(rbind, flu_2022_list)
remove(flu_2022_list)

# Save data
saveRDS(flu_2022, 'tmp/flu_2022_predictions_v2.rds')

################################################################################
################################################################################
