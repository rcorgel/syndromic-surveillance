################################################################################
# File Name: 05b_predict_rsv_data                                              #
#                                                                              #
# Purpose:   Predict RSV from medical claims symptoms.                         #
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

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

##################
# 2. LOAD MODELS #
##################

# Load V1 Models
fever_model <- readRDS("./tmp/fever_model_v1_rsv.rds")
all_symptom_model  <- readRDS("./tmp/all_symptom_model_v1_rsv.rds")
symptom_interaction_model  <- readRDS("./tmp/symptom_interaction_model_v1_rsv.rds")
interaction_case_model  <- readRDS("./tmp/interaction_case_model_v1_rsv.rds")

# Cutpoints
cutpoints_v1 <- readRDS("./tmp/model_cutpoints_v1_rsv.rds")
cutpoints_v2 <- readRDS("./tmp/model_cutpoints_v2_rsv.rds")

# Case per visit
cases_visits <- readRDS("./tmp/visit_positivity_county_rsv.rds")

# Set prediction function
predict_rsv <- function(df, cutoff_df) {
  
  # Made predictions based off the four models
  df$pred_fever <- predict(fever_model, newdata = df, type = "response")
  df$pred_symptom <- predict(all_symptom_model, newdata = df, type = "response")
  df$pred_interaction <- predict(symptom_interaction_model, newdata = df, type = "response")
  df$pred_case <- predict(interaction_case_model, newdata = df, type = "response")

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
           pred_case_count = ifelse(is.na(pred_case_count), 0, pred_case_count),
           pred_case_count = pred_case_count * patient_count) |>
    # Collapse observations and predictions by week, county, and age group
    group_by(year_week_dt, county_fips, age_grp) |>
    mutate(rsv_count = sum(rsv * patient_count),
           pred_fever_count = sum(pred_fever_count),
           pred_symptom_count = sum(pred_symptom_count),
           pred_interaction_count = sum(pred_interaction_count),
           pred_case_count = sum(pred_case_count)) |>
    distinct(year_week_dt, county_fips, age_grp, rsv_count, 
             pred_fever_count, pred_symptom_count,
             pred_interaction_count, pred_case_count)
  
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
rsv_2016 <- read_parquet('tmp/rsv_1_proc_p2.parquet')

# Megge on visit data
rsv_2016 <- left_join(rsv_2016, cases_visits[, c(1, 2, 5)], by = c('county_fips', 'year_week_dt'))

# Split data by week
rsv_2016_list <- split(rsv_2016, rsv_2016$year_week_dt)
remove(rsv_2016)

# Predict data
rsv_2016_list <- lapply(rsv_2016_list, predict_rsv, cutoff_df = cutpoints_v1)

# Combine data
rsv_2016 <- do.call(rbind, rsv_2016_list)
remove(rsv_2016_list)

# Save data
saveRDS(rsv_2016, 'tmp/rsv_2016_predictions_v1.rds')

########
# 2017 #
########

# Load 2017 data
rsv_2017 <- read_parquet('tmp/rsv_2_proc_p2.parquet')

# Megge on visit data
rsv_2017 <- left_join(rsv_2017, cases_visits[, c(1, 2, 5)], by = c('county_fips', 'year_week_dt'))

# Split data by week
rsv_2017_list <- split(rsv_2017, rsv_2017$year_week_dt)
remove(rsv_2017)

# Predict data
rsv_2017_list <- lapply(rsv_2017_list, predict_rsv, cutoff_df = cutpoints_v1)

# Combine data
rsv_2017 <- do.call(rbind, rsv_2017_list)
remove(rsv_2017_list)

# Save data
saveRDS(rsv_2017, 'tmp/rsv_2017_predictions_v1.rds')

########
# 2018 #
########

# Load 2018 data
rsv_2018 <- read_parquet('tmp/rsv_3_proc_p2.parquet')

# Megge on visit data
rsv_2018 <- left_join(rsv_2018, cases_visits[, c(1, 2, 5)], by = c('county_fips', 'year_week_dt'))

# Split data by week
rsv_2018_list <- split(rsv_2018, rsv_2018$year_week_dt)
remove(rsv_2018)

# Predict data
rsv_2018_list <- lapply(rsv_2018_list, predict_rsv, cutoff_df = cutpoints_v1)

# Combine data
rsv_2018 <- do.call(rbind, rsv_2018_list)
remove(rsv_2018_list)

# Save data
saveRDS(rsv_2018, 'tmp/rsv_2018_predictions_v1.rds')

########
# 2019 #
########

# Load 2019 data
rsv_2019 <- read_parquet('tmp/rsv_4_proc_p2.parquet')

# Megge on visit data
rsv_2019 <- left_join(rsv_2019, cases_visits[, c(1, 2, 5)], by = c('county_fips', 'year_week_dt'))

# Split data by week
rsv_2019_list <- split(rsv_2019, rsv_2019$year_week_dt)
remove(rsv_2019)

# Predict data
rsv_2019_list <- lapply(rsv_2019_list, predict_rsv, cutoff_df = cutpoints_v1)

# Combine data
rsv_2019 <- do.call(rbind, rsv_2019_list)
remove(rsv_2019_list)

# Save data
saveRDS(rsv_2019, 'tmp/rsv_2019_predictions_v1.rds')

########
# 2020 #
########

# Load 2020 data
rsv_2020 <- read_parquet('tmp/rsv_5_proc_p2.parquet')

# Megge on visit data
rsv_2020 <- left_join(rsv_2020, cases_visits[, c(1, 2, 5)], by = c('county_fips', 'year_week_dt'))

# Split data by week
rsv_2020_list <- split(rsv_2020, rsv_2020$year_week_dt)
remove(rsv_2020)

# Predict data
rsv_2020_list <- lapply(rsv_2020_list, predict_rsv, cutoff_df = cutpoints_v1)

# Combine data
rsv_2020 <- do.call(rbind, rsv_2020_list)
remove(rsv_2020_list)

# Save data
saveRDS(rsv_2020, 'tmp/rsv_2020_predictions_v1.rds')

########
# 2021 #
########

# Load 2021 data
rsv_2021 <- read_parquet('tmp/rsv_6_proc_p2.parquet')

# Megge on visit data
rsv_2021 <- left_join(rsv_2021, cases_visits[, c(1, 2, 5)], by = c('county_fips', 'year_week_dt'))

# Split data by week
rsv_2021_list <- split(rsv_2021, rsv_2021$year_week_dt)
remove(rsv_2021)

# Predict data
rsv_2021_list <- lapply(rsv_2021_list, predict_rsv, cutoff_df = cutpoints_v1)

# Combine data
rsv_2021 <- do.call(rbind, rsv_2021_list)
remove(rsv_2021_list)

# Save data
saveRDS(rsv_2021, 'tmp/rsv_2021_predictions_v1.rds')

########
# 2022 #
########

# Load 2022 data
rsv_2022 <- read_parquet('tmp/rsv_7_proc_p2.parquet')

# Megge on visit data
rsv_2022 <- left_join(rsv_2022, cases_visits[, c(1, 2, 5)], by = c('county_fips', 'year_week_dt'))

# Split data by week
rsv_2022_list <- split(rsv_2022, rsv_2022$year_week_dt)
remove(rsv_2022)

# Predict data
rsv_2022_list <- lapply(rsv_2022_list, predict_rsv, cutoff_df = cutpoints_v1)

# Combine data
rsv_2022 <- do.call(rbind, rsv_2022_list)
remove(rsv_2022_list)

# Save data
saveRDS(rsv_2022, 'tmp/rsv_2022_predictions_v1.rds')

## Load V2 Models ##
fever_model <- readRDS("./tmp/fever_model_v2_rsv.rds")
all_symptom_model  <- readRDS("./tmp/all_symptom_model_v2_rsv.rds")
symptom_interaction_model  <- readRDS("./tmp/symptom_interaction_model_v2_rsv.rds")
interaction_case_model  <- readRDS("./tmp/interaction_case_model_v2_rsv.rds")

########
# 2020 #
########

# Load 2020 data
rsv_2020 <- read_parquet('tmp/rsv_5_proc_p2.parquet')

# Megge on visit data
rsv_2020 <- left_join(rsv_2020, cases_visits[, c(1, 2, 5)], by = c('county_fips', 'year_week_dt'))

# Split data by week
rsv_2020_list <- split(rsv_2020, rsv_2020$year_week_dt)
remove(rsv_2020)

# Predict data
rsv_2020_list <- lapply(rsv_2020_list, predict_rsv, cutoff_df = cutpoints_v2)

# Combine data
rsv_2020 <- do.call(rbind, rsv_2020_list)
remove(rsv_2020_list)

# Save data
saveRDS(rsv_2020, 'tmp/rsv_2020_predictions_v2.rds')

########
# 2021 #
########

# Load 2021 data
rsv_2021 <- read_parquet('tmp/rsv_6_proc_p2.parquet')

# Megge on visit data
rsv_2021 <- left_join(rsv_2021, cases_visits[, c(1, 2, 5)], by = c('county_fips', 'year_week_dt'))

# Split data by week
rsv_2021_list <- split(rsv_2021, rsv_2021$year_week_dt)
remove(rsv_2021)

# Predict data
rsv_2021_list <- lapply(rsv_2021_list, predict_rsv, cutoff_df = cutpoints_v2)

# Combine data
rsv_2021 <- do.call(rbind, rsv_2021_list)
remove(rsv_2021_list)

# Save data
saveRDS(rsv_2021, 'tmp/rsv_2021_predictions_v2.rds')

########
# 2022 #
########

# Load 2022 data
rsv_2022 <- read_parquet('tmp/rsv_7_proc_p2.parquet')

# Megge on visit data
rsv_2022 <- left_join(rsv_2022, cases_visits[, c(1, 2, 5)], by = c('county_fips', 'year_week_dt'))

# Split data by week
rsv_2022_list <- split(rsv_2022, rsv_2022$year_week_dt)
remove(rsv_2022)

# Predict data
rsv_2022_list <- lapply(rsv_2022_list, predict_rsv, cutoff_df = cutpoints_v12)

# Combine data
rsv_2022 <- do.call(rbind, rsv_2022_list)
remove(rsv_2022_list)

# Save data
saveRDS(rsv_2022, 'tmp/rsv_2022_predictions_v2.rds')

################################################################################
################################################################################