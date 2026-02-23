################################################################################
# File Name: 03_predict_flu_data                                               #
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

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

##################
# 2. LOAD MODELS #
##################

# Models
fever_model <- readRDS("./tmp/fever_model.rds")
all_symptom_model <- readRDS("./tmp/all_symptom_model.rds")
symptom_interaction_model <- readRDS("./tmp/symptom_interaction_model.rds")
interaction_case_model <- readRDS("./tmp/interaction_case_model.rds")

# Cutpoints
cutpoints <- readRDS("./tmp/model_cutpoints.rds")

# Case per visit
cases_visits <- readRDS("./tmp/visit_positivity_county.rds")

# Set prediction function
predict_flu <- function(df, cutoff_df) {
  
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
    # Collapse observations and predictions by week, county, and age group
    group_by(year_week_dt, county_fips, age_grp) |>
    mutate(flu_count = sum(flu * patient_count),
           pred_fever_count = sum(pred_fever_count),
           pred_symptom_count = sum(pred_symptom_count),
           pred_interaction_count = sum(pred_interaction_count),
           pred_case_count = sum(pred_case_count)) |>
    distinct(year_week_dt, county_fips, age_grp, flu_count, 
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
flu_2016 <- read_parquet('tmp/flu_1_proc.parquet')
flu_2016 <- as.data.frame(flu_2016)
flu_2016 <- left_join(flu_2016, cases_visits[, c(1, 4, 10)], by = c('county_fips', 'year_week_dt'))

# Split data by week
flu_2016_list <- split(flu_2016, flu_2016$year_week_dt)
remove(flu_2016)

# Predict data
flu_2016_list <- lapply(flu_2016_list, predict_flu, cutoff_df = cutpoints)

# Combine data
flu_2016 <- do.call(rbind, flu_2016_list)
remove(flu_2016_list)

# Save data
saveRDS(flu_2016, 'tmp/flu_2016_predictions')

########
# 2017 #
########

# Load 2017 data
flu_2017 <- read_parquet('tmp/flu_2_proc.parquet')
flu_2017 <- as.data.frame(flu_2017)
flu_2017 <- left_join(flu_2017, cases_visits[, c(1, 4, 10)], by = c('county_fips', 'year_week_dt'))

# Split data by week
flu_2017_list <- split(flu_2017, flu_2017$year_week_dt)
remove(flu_2017)

# Predict data
flu_2017_list <- lapply(flu_2017_list, predict_flu, cutoff_df = cutpoints)

# Combine data
flu_2017 <- do.call(rbind, flu_2017_list)
remove(flu_2017_list)

# Save data
saveRDS(flu_2017, 'tmp/flu_2017_predictions')

########
# 2018 #
########

# Load 2018 data
flu_2018 <- read_parquet('tmp/flu_3_proc.parquet')
flu_2018 <- as.data.frame(flu_2018)
flu_2018 <- left_join(flu_2018, cases_visits[, c(1, 4, 10)], by = c('county_fips', 'year_week_dt'))

# Split data by week
flu_2018_list <- split(flu_2018, flu_2018$year_week_dt)
remove(flu_2018)

# Predict data
flu_2018_list <- lapply(flu_2018_list, predict_flu, cutoff_df = cutpoints)

# Combine data
flu_2018 <- do.call(rbind, flu_2018_list)
remove(flu_2018_list)

# Save data
saveRDS(flu_2018, 'tmp/flu_2018_predictions')

########
# 2019 #
########

# Load 2019 data
flu_2019 <- read_parquet('tmp/flu_4_proc.parquet')
flu_2019 <- as.data.frame(flu_2019)
flu_2019 <- left_join(flu_2019, cases_visits[, c(1, 4, 10)], by = c('county_fips', 'year_week_dt'))

# Split data by week
flu_2019_list <- split(flu_2019, flu_2019$year_week_dt)
remove(flu_2019)

# Predict data
flu_2019_list <- lapply(flu_2019_list, predict_flu, cutoff_df = cutpoints)

# Combine data
flu_2019 <- do.call(rbind, flu_2019_list)
remove(flu_2019_list)

# Save data
saveRDS(flu_2019, 'tmp/flu_2019_predictions')

################################################################################
################################################################################
