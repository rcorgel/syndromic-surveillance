################################################################################
# File Name: 02_regress_flu_data                                               #
#                                                                              #
# Purpose:   Subset processed flu data and conduct regression analysis.        #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load & subset data                                             #
#            3. Calculate cases per visit                                      #
#            4. Test different models                                          #
#            5. Run full models & save                                         #
#            6. Determine cutoff values                                        #
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
library(ISOweek)
library(pROC)
library(mgcv)
library(cutpointr)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

#########################
# 2. LOAD & SUBSET DATA #
#########################

# Load and subset data
# 2016
flu_1 <- read_parquet('tmp/flu_1_proc.parquet')
flu_balanced_16 <- flu_1 |> 
  group_by(flu) |> 
  uncount(patient_count) |> 
  sample_n(125000)
remove(flu_1)
# 2017
flu_2 <- read_parquet('tmp/flu_2_proc.parquet')
flu_balanced_17 <- flu_2 |> 
  group_by(flu) |> 
  uncount(patient_count) |> 
  sample_n(125000)
remove(flu_2)
# 2018
flu_3 <- read_parquet('tmp/flu_3_proc.parquet')
flu_balanced_18 <- flu_3 |> 
  group_by(flu) |> 
  uncount(patient_count) |> 
  sample_n(125000)
remove(flu_3)
# 2019
flu_4 <- read_parquet('tmp/flu_4_proc.parquet')
flu_balanced_19 <- flu_4 |> 
  group_by(flu) |> 
  uncount(patient_count) |> 
  sample_n(125000)
remove(flu_4)

# Combine balanced data
flu <- rbind(flu_balanced_16, flu_balanced_17, flu_balanced_18, flu_balanced_19)
saveRDS(flu, 'tmp/flu_sample.rds')

################################
# 3. CALCULATE CASES PER VISIT #
################################

# Load flu cases and visit
cases <- read.csv('raw/county_week_flu_v3_imputed (1).csv')
visits <- read.csv('raw/county_week_ac_v3_imputed.csv')

# Convert year week to date
visits$year_week_dt = sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", visits$year_week)
visits$year_week_dt = as.character(ISOweek2date(visits$year_week_dt))

# Merge cases and visits
cases_visits <- left_join(cases, visits, by = c('county_fips', 'year_week_dt'))

# Calculate cases / visits
cases_visits$visit_pos <- cases_visits$conf_flu / cases_visits$all_cause_wtd

# Add leading 0 to county fips
cases_visits$county_fips = ifelse(nchar(cases_visits$county_fips) == 4, 
                                 paste0("0", cases_visits$county_fips), 
                                 cases_visits$county_fips)

# Convert to character data to date
cases_visits$year_week_dt <- as.Date(cases_visits$year_week_dt)

# Add 7 days to the week date such that date is lagged when merged
cases_visits$year_week_dt <- cases_visits$year_week_dt + 7

# Save
saveRDS(cases_visits, "./tmp/visit_positivity_county.rds")

############################
# 4. TEST DIFFERENT MODELS #
############################

# Re-load sampled flu data
flu <- readRDS('tmp/flu_sample.rds')

# Merge to sampled data
flu <- left_join(flu, cases_visits[, c(1, 4, 10)], by = c('county_fips', 'year_week_dt'))

# Conduct 10-fold cross validation to test multiple models

# Set the number of folds
folds <- sample(1:10, nrow(flu), replace = TRUE)

# Create empty vector to fill
roc_values <- rep(NA, 10) 

####################
# FEVER ONLY MODEL #
####################

# Run the cross validation
for (i in 1:10) {
  # Set training data
  train_data <- flu[folds != i, ]
  
  # Set testing data
  test_data <- flu[folds == i, ]
  
  # Fit the model on the training data
  model_fit <- glm(flu ~ fever, 
                   data = train_data,
                   family = binomial(link = "logit"))
  
  # Make predictions on the test data
  test_data$predictions <- predict(model_fit, newdata = test_data, type = "response")
  
  # Calculate the ROC
  roc <- pROC::auc(pROC::roc(test_data$flu, test_data$predictions))[1]
  roc_values[i] <- roc
}

# Calculate the average ROC across all folds
print(paste("Average K-fold ROC for the fever only model:", round(mean(roc_values), 3)))

#####################
# ALL SYMPTOM MODEL #
#####################

# Reset roc values
roc_values <- rep(NA, 10) 

# Run the cross validation
for (i in 1:10) {
  # Set training data
  train_data <- flu[folds != i, ]
  
  # Set testing data
  test_data <- flu[folds == i, ]
  
  # Fit the model on the training data
  model_fit <- glm(flu ~ fever + 
                     cough + 
                     sore_throat +
                     myalgia + 
                     short_breath + 
                     hypoxemia + 
                     nausea_vom + 
                     bronchitis + 
                     chest_pain + 
                     diarrhea + 
                     fatigue + 
                     headache + 
                     congestion + 
                     sneezing, 
                   data = train_data,
                   family = binomial(link = "logit"))
  
  # Make predictions on the test data
  test_data$predictions <- predict(model_fit, newdata = test_data, type = "response")
  
  # Calculate the ROC
  roc <- pROC::auc(pROC::roc(test_data$flu, test_data$predictions))[1]
  roc_values[i] <- roc
}

# Calculate the average ROC across all folds
print(paste("Average K-fold ROC for the all symptom model:", round(mean(roc_values), 3)))

#############################
# SYMPTOM INTERACTION MODEL #
#############################

# Reset roc values
roc_values <- rep(NA, 10) 

# Run the cross validation
for (i in 1:10) {
  # Set training data
  train_data <- flu[folds != i, ]
  
  # Set testing data
  test_data <- flu[folds == i, ]
  
  # Fit the model on the training data
  model_fit <- glm(flu ~ fever*as.factor(age_grp) + 
                     cough*as.factor(age_grp) + 
                     sore_throat*as.factor(age_grp) +
                     myalgia*as.factor(age_grp) + 
                     short_breath*as.factor(age_grp) + 
                     hypoxemia*as.factor(age_grp) + 
                     nausea_vom*as.factor(age_grp) + 
                     bronchitis*as.factor(age_grp) + 
                     chest_pain*as.factor(age_grp) + 
                     diarrhea*as.factor(age_grp) + 
                     fatigue*as.factor(age_grp) + 
                     headache*as.factor(age_grp) + 
                     congestion*as.factor(age_grp) + 
                     sneezing*as.factor(age_grp) +
                     fever*as.factor(patient_gender_code) + 
                     cough*as.factor(patient_gender_code) + 
                     sore_throat*as.factor(patient_gender_code) +
                     myalgia*as.factor(patient_gender_code) + 
                     short_breath*as.factor(patient_gender_code) + 
                     hypoxemia*as.factor(patient_gender_code) + 
                     nausea_vom*as.factor(patient_gender_code) + 
                     bronchitis*as.factor(patient_gender_code) + 
                     chest_pain*as.factor(patient_gender_code) + 
                     diarrhea*as.factor(patient_gender_code) + 
                     fatigue*as.factor(patient_gender_code) + 
                     headache*as.factor(patient_gender_code) + 
                     congestion*as.factor(patient_gender_code) + 
                     sneezing*as.factor(patient_gender_code), 
                   data = train_data,
                   family = binomial(link = "logit"))
  
  # Make predictions on the test data
  test_data$predictions <- predict(model_fit, newdata = test_data, type = "response")
  
  # Calculate the ROC
  roc <- pROC::auc(pROC::roc(test_data$flu, test_data$predictions))[1]
  roc_values[i] <- roc
}

# Calculate the average ROC across all folds
print(paste("Average K-fold ROC for the symptom interaction model:", round(mean(roc_values), 3)))

##########################
# INTERACTION CASE MODEL #
##########################

# Reset roc values
roc_values <- rep(NA, 10) 

# Run the cross validation
for (i in 1:10) {
  # Set training data
  train_data <- flu[folds != i, ]
  
  # Set testing data
  test_data <- flu[folds == i, ]
  
  # Fit the model on the training data
  model_fit <- glm(flu ~ fever*as.factor(age_grp) + 
                     cough*as.factor(age_grp) + 
                     sore_throat*as.factor(age_grp) +
                     myalgia*as.factor(age_grp) + 
                     short_breath*as.factor(age_grp) + 
                     hypoxemia*as.factor(age_grp) + 
                     nausea_vom*as.factor(age_grp) + 
                     bronchitis*as.factor(age_grp) + 
                     chest_pain*as.factor(age_grp) + 
                     diarrhea*as.factor(age_grp) + 
                     fatigue*as.factor(age_grp) + 
                     headache*as.factor(age_grp) + 
                     congestion*as.factor(age_grp) + 
                     sneezing*as.factor(age_grp) +
                     fever*as.factor(patient_gender_code) + 
                     cough*as.factor(patient_gender_code) + 
                     sore_throat*as.factor(patient_gender_code) +
                     myalgia*as.factor(patient_gender_code) + 
                     short_breath*as.factor(patient_gender_code) + 
                     hypoxemia*as.factor(patient_gender_code) + 
                     nausea_vom*as.factor(patient_gender_code) + 
                     bronchitis*as.factor(patient_gender_code) + 
                     chest_pain*as.factor(patient_gender_code) + 
                     diarrhea*as.factor(patient_gender_code) + 
                     fatigue*as.factor(patient_gender_code) + 
                     headache*as.factor(patient_gender_code) + 
                     congestion*as.factor(patient_gender_code) + 
                     sneezing*as.factor(patient_gender_code) +
                     visit_pos,
                     #s(visit_pos, bs = "cc", k = 20), 
                   data = train_data,
                   family = binomial(link = "logit"))
  
  # Make predictions on the test data
  test_data$predictions <- predict(model_fit, newdata = test_data, type = "response")
  
  # Calculate the ROC
  roc <- pROC::auc(pROC::roc(test_data$flu, test_data$predictions))[1]
  roc_values[i] <- roc
}

# Calculate the average ROC across all folds
print(paste("Average K-fold ROC for the interaction case model:", round(mean(roc_values), 3)))

#############################
# 5. RUN FULL MODELS & SAVE #
#############################

# Fever model
fever_model <- glm(flu ~ fever, 
                   data = flu,
                   family = binomial(link = "logit"))

# All symptom model
all_symptom_model <- glm(flu ~ fever + 
                           cough + 
                           sore_throat +
                           myalgia + 
                           short_breath + 
                           hypoxemia + 
                           nausea_vom + 
                           bronchitis + 
                           chest_pain + 
                           diarrhea + 
                           fatigue + 
                           headache + 
                           congestion + 
                           sneezing, 
                         data = flu,
                         family = binomial(link = "logit"))

# Symptom interaction model
symptom_interaction_model <- glm(flu ~ fever*as.factor(age_grp) + 
                                   cough*as.factor(age_grp) + 
                                   sore_throat*as.factor(age_grp) +
                                   myalgia*as.factor(age_grp) + 
                                   short_breath*as.factor(age_grp) + 
                                   hypoxemia*as.factor(age_grp) + 
                                   nausea_vom*as.factor(age_grp) + 
                                   bronchitis*as.factor(age_grp) + 
                                   chest_pain*as.factor(age_grp) + 
                                   diarrhea*as.factor(age_grp) + 
                                   fatigue*as.factor(age_grp) + 
                                   headache*as.factor(age_grp) + 
                                   congestion*as.factor(age_grp) + 
                                   sneezing*as.factor(age_grp) +
                                   fever*as.factor(patient_gender_code) + 
                                   cough*as.factor(patient_gender_code) + 
                                   sore_throat*as.factor(patient_gender_code) +
                                   myalgia*as.factor(patient_gender_code) + 
                                   short_breath*as.factor(patient_gender_code) + 
                                   hypoxemia*as.factor(patient_gender_code) + 
                                   nausea_vom*as.factor(patient_gender_code) + 
                                   bronchitis*as.factor(patient_gender_code) + 
                                   chest_pain*as.factor(patient_gender_code) + 
                                   diarrhea*as.factor(patient_gender_code) + 
                                   fatigue*as.factor(patient_gender_code) + 
                                   headache*as.factor(patient_gender_code) + 
                                   congestion*as.factor(patient_gender_code) + 
                                   sneezing*as.factor(patient_gender_code), 
                                 data = flu,
                                 family = binomial(link = "logit"))

# Interaction case model
interaction_case_model <- glm(flu ~ fever*as.factor(age_grp) + 
                                cough*as.factor(age_grp) + 
                                sore_throat*as.factor(age_grp) +
                                myalgia*as.factor(age_grp) + 
                                short_breath*as.factor(age_grp) + 
                                hypoxemia*as.factor(age_grp) + 
                                nausea_vom*as.factor(age_grp) + 
                                bronchitis*as.factor(age_grp) + 
                                chest_pain*as.factor(age_grp) + 
                                diarrhea*as.factor(age_grp) + 
                                fatigue*as.factor(age_grp) + 
                                headache*as.factor(age_grp) + 
                                congestion*as.factor(age_grp) + 
                                sneezing*as.factor(age_grp) +
                                fever*as.factor(patient_gender_code) + 
                                cough*as.factor(patient_gender_code) + 
                                sore_throat*as.factor(patient_gender_code) +
                                myalgia*as.factor(patient_gender_code) + 
                                short_breath*as.factor(patient_gender_code) + 
                                hypoxemia*as.factor(patient_gender_code) + 
                                nausea_vom*as.factor(patient_gender_code) + 
                                bronchitis*as.factor(patient_gender_code) + 
                                chest_pain*as.factor(patient_gender_code) + 
                                diarrhea*as.factor(patient_gender_code) + 
                                fatigue*as.factor(patient_gender_code) + 
                                headache*as.factor(patient_gender_code) + 
                                congestion*as.factor(patient_gender_code) + 
                                sneezing*as.factor(patient_gender_code) +
                                visit_pos,
                                #s(visit_pos, bs = "cc", k = 20), 
                              data = flu,
                              family = binomial(link = "logit"))

# Save models
saveRDS(fever_model, "./tmp/fever_model.rds")
saveRDS(all_symptom_model, "./tmp/all_symptom_model.rds")
saveRDS(symptom_interaction_model, "./tmp/symptom_interaction_model.rds")
saveRDS(interaction_case_model, "./tmp/interaction_case_model.rds")

##############################
# 6. DETERMINE CUTOFF VALUES #
##############################

# Predict based on models
flu$pred_fever <- predict(fever_model, newdata = flu, type = "response")
flu$pred_symptom <- predict(all_symptom_model, newdata = flu, type = "response")
flu$pred_interaction <- predict(symptom_interaction_model, newdata = flu, type = "response")
flu$pred_case <- predict(interaction_case_model, newdata = flu, type = "response")

# Determine cutpoints by age group
cutpoints <- flu |> group_by(age_grp) |>
  mutate(cutoff_fever = cutpointr(pred_fever, flu, 
                                  method = maximize_metric, 
                                  metric = youden, 
                                  na.rm = TRUE)$optimal_cutpoint) |>
  mutate(cutoff_symptom = cutpointr(pred_symptom, flu, 
                                  method = maximize_metric, 
                                  metric = youden, 
                                  na.rm = TRUE)$optimal_cutpoint) |>
  mutate(cutoff_interaction = cutpointr(pred_interaction, flu, 
                                     method = maximize_metric, 
                                     metric = youden, 
                                     na.rm = TRUE)$optimal_cutpoint) |>
  mutate(cutoff_case = cutpointr(pred_case, flu, 
                                     method = maximize_metric, 
                                     metric = youden, 
                                     na.rm = TRUE)$optimal_cutpoint) |>
  distinct(age_grp, cutoff_fever, cutoff_symptom, cutoff_interaction, cutoff_case)

# Save
saveRDS(cutpoints, "./tmp/model_cutpoints.rds")

################################################################################
################################################################################
