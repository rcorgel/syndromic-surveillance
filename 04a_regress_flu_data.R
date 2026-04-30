################################################################################
# File Name: 04a_regress_flu_data                                              #
#                                                                              #
# Purpose:   Conduct regression analysis on sampled flu data.                  #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load and merge auxiliary data                                  #
#            4. Test different models                                          #
#            5. Run full models & save                                         #
#            6. Determine cutoff values                                        #
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
library(ISOweek)
library(pROC)
library(mgcv)
library(cutpointr)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

####################################
# 2. LOAD AND MERGE AUXILIARY DATA #
####################################

# Load flu samples data
flu_v1 <- readRDS('tmp/flu_sample_v1.rds')
flu_v2 <- readRDS('tmp/flu_sample_v2.rds')

# Load flu cases, tests, and all cause visits
cases <- read.csv('raw/county_week_flu_v3_imputed.csv')
visits <- read.csv('raw/county_week_ac_v3_imputed.csv')
tests <- read_csv('raw/county_week_flu_tests_v2_imputed.csv')

# Convert year week to date
visits$year_week_dt = sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", visits$year_week)
visits$year_week_dt = as.character(ISOweek2date(visits$year_week_dt))

# Merge cases and visits
cases_visits <- left_join(cases, visits, by = c('county_fips', 'year_week_dt'))

# Add leading 0 to county fips
cases_visits$county_fips = ifelse(nchar(cases_visits$county_fips) == 4, 
                                  paste0("0", cases_visits$county_fips), 
                                  cases_visits$county_fips)

# Convert character date to date
cases_visits$year_week_dt <- as.Date(cases_visits$year_week_dt)

# Merge on tests
cases_visits <- left_join(cases_visits, tests, by = c('county_fips', 'year_week_dt'))

# Calculate cases / visits
cases_visits$visit_pos <- cases_visits$conf_flu / cases_visits$all_cause_wtd

# Calculate tests / visits
cases_visits$tests_per <- cases_visits$flu_tests / cases_visits$all_cause_wtd

# Add 7 days to the week date such that date is lagged when merged
cases_visits$year_week_dt <- cases_visits$year_week_dt + 7

# Select variables of interest
cases_visits <- cases_visits |> dplyr::select(c(county_fips, year_week_dt, conf_flu, 
                                                all_cause_wtd, flu_tests, visit_pos,
                                                tests_per))

# Save
saveRDS(cases_visits, "./tmp/visit_positivity_county_flu.rds")

# Merge on to sample data
flu_v1 <- left_join(flu_v1, cases_visits, by = c('county_fips', 'year_week_dt'))
flu_v2 <- left_join(flu_v2, cases_visits, by = c('county_fips', 'year_week_dt'))

# Load NREVSS data
nrevss_state <- readRDS('./tmp/nrevss_state.rds')
nrevss_nat <- readRDS('./tmp/nrevss_national.rds')

# Lag by a week
nrevss_state$week_date <- nrevss_state$week_date + 7
nrevss_nat$week_date <- nrevss_nat$week_date + 7

# Merge on NREVSS
flu_v1 <- left_join(flu_v1, nrevss_state[, c(3, 4, 5)], by = c('year_week_dt' = 'week_date',
                                                         'state_fips' = 'state_fips'))
flu_v1 <- left_join(flu_v1, nrevss_nat[, c(3, 5)], by = c('year_week_dt' = 'week_date'))
flu_v2 <- left_join(flu_v2, nrevss_state[, c(3, 4, 5)], by = c('year_week_dt' = 'week_date',
                                                               'state_fips' = 'state_fips'))
flu_v2 <- left_join(flu_v2, nrevss_nat[, c(3, 5)], by = c('year_week_dt' = 'week_date'))

# Replace missing state values with national
flu_v1$nrevss_pos <- ifelse(is.na(flu_v1$value_roll.x), flu_v1$value_roll.y, flu_v1$value_roll.x)
flu_v2$nrevss_pos <- ifelse(is.na(flu_v2$value_roll.x), flu_v2$value_roll.y, flu_v2$value_roll.x)

############################
# 3. TEST DIFFERENT MODELS #
############################

# Conduct 10-fold cross validation to test multiple models

# Set the number of folds
folds_v1 <- sample(1:10, nrow(flu_v1), replace = TRUE)
folds_v2 <- sample(1:10, nrow(flu_v2), replace = TRUE)

# Create empty vector to fill
auc_values_v1 <- rep(NA, 10) 
auc_values_v2 <- rep(NA, 10) 
roc_curve_v1 <- list(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
roc_curve_v2 <- list(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

####################
# FEVER ONLY MODEL #
####################

# Run the cross validation
for (i in 1:10) {
  # Set training data
  train_data_v1 <- flu_v1[folds_v1 != i, ]
  train_data_v2 <- flu_v2[folds_v2 != i, ]
  
  # Set testing data
  test_data_v1 <- flu_v1[folds_v1 == i, ]
  test_data_v2 <- flu_v2[folds_v2 == i, ]
  
  # Fit the model on the training data
  model_fit_v1 <- glm(flu ~ fever, 
                   data = train_data_v1,
                   family = binomial(link = "logit"))
  model_fit_v2 <- glm(flu ~ fever, 
                   data = train_data_v2,
                   family = binomial(link = "logit"))
  
  # Make predictions on the test data
  test_data_v1$predictions <- stats::predict(model_fit_v1, newdata = test_data_v1, type = "response")
  test_data_v2$predictions <- stats::predict(model_fit_v2, newdata = test_data_v2, type = "response")
  
  # Calculate the AUC & ROC
  roc_v1 <- pROC::roc(test_data_v1$flu, test_data_v1$predictions)
  roc_curve_v1[[i]] <- data.frame(sensitivity = roc_v1$sensitivities, specificity = roc_v1$specificities, run = i)
  auc_values_v1[i] <- pROC::auc(roc_v1)[1]
  roc_v2 <- pROC::roc(test_data_v2$flu, test_data_v2$predictions)
  roc_curve_v2[[i]] <- data.frame(sensitivity = roc_v2$sensitivities, specificity = roc_v2$specificities, run = i)
  auc_values_v2[i] <- pROC::auc(roc_v2)[1]
}

# Calculate the average ROC across all folds
print(paste("Average K-fold ROC for the fever only model (2016-2020):", round(mean(auc_values_v1), 3)))
print(paste("Average K-fold ROC for the fever only model (2020-2023):", round(mean(auc_values_v2), 3)))

# Append ROC data for plotting
roc_data_v1 <- do.call(rbind, roc_curve_v1)
roc_data_v2 <- do.call(rbind, roc_curve_v2)

# Plot
line_theme <- theme(legend.position = "none",
                    axis.text = element_text(size=14),
                    axis.title = element_text(size=16),
                    strip.background = element_blank(),
                    plot.title = element_text(size=16))

auc_v1 <- round(mean(auc_values_v1), 3)
fever_roc_v1 <- ggplot(roc_data_v1, aes(x = 1 - specificity, y = sensitivity, group = as.factor(run))) +
  geom_line(linewidth = 2, alpha = 0.5, color = '#F8766D') + # Use geom_line or geom_step
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1.2) +
  labs(
    title = "Fever Model",
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme_bw() + line_theme +
  coord_equal() +
  annotate("text", x = 0.75, y = 0.25, label = glue::glue("AUC = {auc_v1}"), color = "black", size = 5)
fever_roc_v1

auc_v2 <- round(mean(auc_values_v2), 3)
fever_roc_v2 <- ggplot(roc_data_v2, aes(x = 1 - specificity, y = sensitivity, group = as.factor(run))) +
  geom_line(linewidth = 2, alpha = 0.5, color = '#F8766D') + # Use geom_line or geom_step
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1.2) +
  labs(
    title = "Fever Model",
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme_bw() + line_theme +
  coord_equal() +
  annotate("text", x = 0.75, y = 0.25, label = glue::glue("AUC = {auc_v2}"), color = "black", size = 5)
fever_roc_v2

#####################
# ALL SYMPTOM MODEL #
#####################

# Reset roc values
auc_values_v1 <- rep(NA, 10) 
auc_values_v2 <- rep(NA, 10) 
roc_curve_v1 <- list(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
roc_curve_v2 <- list(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

# Run the cross validation
for (i in 1:10) {
  # Set training data
  train_data_v1 <- flu_v1[folds_v1 != i, ]
  train_data_v2 <- flu_v2[folds_v2 != i, ]
  
  # Set testing data
  test_data_v1 <- flu_v1[folds_v1 == i, ]
  test_data_v2 <- flu_v2[folds_v2 == i, ]
  
  # Fit the model on the training data
  model_fit_v1 <- glm(flu ~ fever + 
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
                      data = train_data_v1,
                      family = binomial(link = "logit"))
  model_fit_v2 <- glm(flu ~ fever + 
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
                      data = train_data_v2,
                      family = binomial(link = "logit"))
  
  # Make predictions on the test data
  test_data_v1$predictions <- stats::predict(model_fit_v1, newdata = test_data_v1, type = "response")
  test_data_v2$predictions <- stats::predict(model_fit_v2, newdata = test_data_v2, type = "response")
  
  # Calculate the AUC & ROC
  roc_v1 <- pROC::roc(test_data_v1$flu, test_data_v1$predictions)
  roc_curve_v1[[i]] <- data.frame(sensitivity = roc_v1$sensitivities, specificity = roc_v1$specificities, run = i)
  auc_values_v1[i] <- pROC::auc(roc_v1)[1]
  roc_v2 <- pROC::roc(test_data_v2$flu, test_data_v2$predictions)
  roc_curve_v2[[i]] <- data.frame(sensitivity = roc_v2$sensitivities, specificity = roc_v2$specificities, run = i)
  auc_values_v2[i] <- pROC::auc(roc_v2)[1]
}

# Calculate the average ROC across all folds
print(paste("Average K-fold ROC for the symptom model (2016-2020):", round(mean(auc_values_v1), 3)))
print(paste("Average K-fold ROC for the symptom model (2020-2023):", round(mean(auc_values_v2), 3)))

# Append ROC data for plotting
roc_data_v1 <- do.call(rbind, roc_curve_v1)
roc_data_v2 <- do.call(rbind, roc_curve_v2)

# Plot
auc_v1 <- round(mean(auc_values_v1), 3)
symptom_roc_v1 <- ggplot(roc_data_v1, aes(x = 1 - specificity, y = sensitivity, group = as.factor(run))) +
  geom_line(linewidth = 2, alpha = 0.5, color = '#00B0F6') + # Use geom_line or geom_step
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1.2) +
  labs(
    title = "Symptom Model",
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme_bw() + line_theme +
  coord_equal() +
  annotate("text", x = 0.75, y = 0.25, label = glue::glue("AUC = {auc_v1}"), color = "black", size = 5)
symptom_roc_v1

auc_v2 <- round(mean(auc_values_v2), 3)
symptom_roc_v2 <- ggplot(roc_data_v2, aes(x = 1 - specificity, y = sensitivity, group = as.factor(run))) +
  geom_line(linewidth = 2, alpha = 0.5, color = '#00B0F6') + # Use geom_line or geom_step
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1.2) +
  labs(
    title = "Symptom Model",
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme_bw() + line_theme +
  coord_equal() +
  annotate("text", x = 0.75, y = 0.25, label = glue::glue("AUC = {auc_v2}"), color = "black", size = 5)
symptom_roc_v2

#############################
# SYMPTOM INTERACTION MODEL #
#############################

# Reset roc values
auc_values_v1 <- rep(NA, 10) 
auc_values_v2 <- rep(NA, 10) 
roc_curve_v1 <- list(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
roc_curve_v2 <- list(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

# Run the cross validation
for (i in 1:10) {
  # Set training data
  train_data_v1 <- flu_v1[folds_v1 != i, ]
  train_data_v2 <- flu_v2[folds_v2 != i, ]
  
  # Set testing data
  test_data_v1 <- flu_v1[folds_v1 == i, ]
  test_data_v2 <- flu_v2[folds_v2 == i, ]
  
  # Fit the model on the training data
  model_fit_v1 <- glm(flu ~ fever*as.factor(age_grp) + 
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
                      data = train_data_v1,
                      family = binomial(link = "logit"))
  model_fit_v2 <- glm(flu ~ fever*as.factor(age_grp) + 
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
                      data = train_data_v2,
                      family = binomial(link = "logit"))
  
  # Make predictions on the test data
  test_data_v1$predictions <- stats::predict(model_fit_v1, newdata = test_data_v1, type = "response")
  test_data_v2$predictions <- stats::predict(model_fit_v2, newdata = test_data_v2, type = "response")
  
  # Calculate the AUC & ROC
  roc_v1 <- pROC::roc(test_data_v1$flu, test_data_v1$predictions)
  roc_curve_v1[[i]] <- data.frame(sensitivity = roc_v1$sensitivities, specificity = roc_v1$specificities, run = i)
  auc_values_v1[i] <- pROC::auc(roc_v1)[1]
  roc_v2 <- pROC::roc(test_data_v2$flu, test_data_v2$predictions)
  roc_curve_v2[[i]] <- data.frame(sensitivity = roc_v2$sensitivities, specificity = roc_v2$specificities, run = i)
  auc_values_v2[i] <- pROC::auc(roc_v2)[1]
}

# Calculate the average ROC across all folds
print(paste("Average K-fold ROC for the interaction model (2016-2020):", round(mean(auc_values_v1), 3)))
print(paste("Average K-fold ROC for the interaction model (2020-2023):", round(mean(auc_values_v2), 3)))

# Append ROC data for plotting
roc_data_v1 <- do.call(rbind, roc_curve_v1)
roc_data_v2 <- do.call(rbind, roc_curve_v2)

# Plot
auc_v1 <- round(mean(auc_values_v1), 3)
interaction_roc_v1 <- ggplot(roc_data_v1, aes(x = 1 - specificity, y = sensitivity, group = as.factor(run))) +
  geom_line(linewidth = 2, alpha = 0.5, color = '#00BF7D') + # Use geom_line or geom_step
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1.2) +
  labs(
    title = "Interaction Model",
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme_bw() + line_theme +
  coord_equal() +
  annotate("text", x = 0.75, y = 0.25, label = glue::glue("AUC = {auc_v1}"), color = "black", size = 5)
interaction_roc_v1

auc_v2 <- round(mean(auc_values_v2), 3)
interaction_roc_v2 <- ggplot(roc_data_v2, aes(x = 1 - specificity, y = sensitivity, group = as.factor(run))) +
  geom_line(linewidth = 2, alpha = 0.5, color = '#00BF7D') + # Use geom_line or geom_step
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1.2) +
  labs(
    title = "Interaction Model",
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme_bw() + line_theme +
  coord_equal() +
  annotate("text", x = 0.75, y = 0.25, label = glue::glue("AUC = {auc_v2}"), color = "black", size = 5)
interaction_roc_v2

##########################
# INTERACTION CASE MODEL #
##########################

# Reset roc values
auc_values_v1 <- rep(NA, 10) 
auc_values_v2 <- rep(NA, 10) 
roc_curve_v1 <- list(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
roc_curve_v2 <- list(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

# Run the cross validation
for (i in 1:10) {
  # Set training data
  train_data_v1 <- flu_v1[folds_v1 != i, ]
  train_data_v2 <- flu_v2[folds_v2 != i, ]
  
  # Set testing data
  test_data_v1 <- flu_v1[folds_v1 == i, ]
  test_data_v2 <- flu_v2[folds_v2 == i, ]
  
  # Fit the model on the training data
  model_fit_v1 <- gam(flu ~ fever*as.factor(age_grp) + 
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
                        s(visit_pos, bs = "cc", k = 20), 
                      data = train_data_v1,
                      family = binomial(link = "logit"))
  model_fit_v2 <- gam(flu ~ fever*as.factor(age_grp) + 
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
                        s(visit_pos, bs = "cc", k = 20), 
                      data = train_data_v2,
                      family = binomial(link = "logit"))
  
  # Make predictions on the test data
  test_data_v1$predictions <- stats::predict(model_fit_v1, newdata = test_data_v1, type = "response")
  test_data_v2$predictions <- stats::predict(model_fit_v2, newdata = test_data_v2, type = "response")
  
  # Calculate the AUC & ROC
  roc_v1 <- pROC::roc(test_data_v1$flu, test_data_v1$predictions)
  roc_curve_v1[[i]] <- data.frame(sensitivity = roc_v1$sensitivities, specificity = roc_v1$specificities, run = i)
  auc_values_v1[i] <- pROC::auc(roc_v1)[1]
  roc_v2 <- pROC::roc(test_data_v2$flu, test_data_v2$predictions)
  roc_curve_v2[[i]] <- data.frame(sensitivity = roc_v2$sensitivities, specificity = roc_v2$specificities, run = i)
  auc_values_v2[i] <- pROC::auc(roc_v2)[1]
}

# Calculate the average ROC across all folds
print(paste("Average K-fold ROC for the interaction case model:", round(mean(auc_values_v1), 3)))
print(paste("Average K-fold ROC for the interaction case model:", round(mean(auc_values_v2), 3)))

# Append ROC data for plotting
roc_data_v1 <- do.call(rbind, roc_curve_v1)
roc_data_v2 <- do.call(rbind, roc_curve_v2)

# Plot
auc_v1 <- round(mean(auc_values_v1), 3)
case_roc_v1 <- ggplot(roc_data_v1, aes(x = 1 - specificity, y = sensitivity, group = as.factor(run))) +
  geom_line(linewidth = 2, alpha = 0.5, color = '#A3A500') + # Use geom_line or geom_step
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1.2) +
  labs(
    title = "Case Model",
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme_bw() + line_theme +
  coord_equal() +
  annotate("text", x = 0.75, y = 0.25, label = glue::glue("AUC = {auc_v1}"), color = "black", size = 5)
case_roc_v1

auc_v2 <- round(mean(auc_values_v2), 3)
case_roc_v2 <- ggplot(roc_data_v2, aes(x = 1 - specificity, y = sensitivity, group = as.factor(run))) +
  geom_line(linewidth = 2, alpha = 0.5, color = '#A3A500') + # Use geom_line or geom_step
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1.2) +
  labs(
    title = "Case Model",
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme_bw() + line_theme +
  coord_equal() +
  annotate("text", x = 0.75, y = 0.25, label = glue::glue("AUC = {auc_v2}"), color = "black", size = 5)
case_roc_v2

############################
# INTERACTION NREVSS MODEL #
############################

# Reset roc values
auc_values_v1 <- rep(NA, 10) 
auc_values_v2 <- rep(NA, 10) 
roc_curve_v1 <- list(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
roc_curve_v2 <- list(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

# Run the cross validation
for (i in 1:10) {
  # Set training data
  train_data_v1 <- flu_v1[folds_v1 != i, ]
  train_data_v2 <- flu_v2[folds_v2 != i, ]
  
  # Set testing data
  test_data_v1 <- flu_v1[folds_v1 == i, ]
  test_data_v2 <- flu_v2[folds_v2 == i, ]
  
  # Fit the model on the training data
  model_fit_v1 <- gam(flu ~ fever*as.factor(age_grp) + 
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
                        s(nrevss_pos, bs = "cc", k = 20), 
                      data = train_data_v1,
                      family = binomial(link = "logit"))
  model_fit_v2 <- gam(flu ~ fever*as.factor(age_grp) + 
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
                        s(nrevss_pos, bs = "cc", k = 20), 
                      data = train_data_v2,
                      family = binomial(link = "logit"))
  
  # Make predictions on the test data
  test_data_v1$predictions <- stats::predict(model_fit_v1, newdata = test_data_v1, type = "response")
  test_data_v2$predictions <- stats::predict(model_fit_v2, newdata = test_data_v2, type = "response")
  
  # Calculate the AUC & ROC
  roc_v1 <- pROC::roc(test_data_v1$flu, test_data_v1$predictions)
  roc_curve_v1[[i]] <- data.frame(sensitivity = roc_v1$sensitivities, specificity = roc_v1$specificities, run = i)
  auc_values_v1[i] <- pROC::auc(roc_v1)[1]
  roc_v2 <- pROC::roc(test_data_v2$flu, test_data_v2$predictions)
  roc_curve_v2[[i]] <- data.frame(sensitivity = roc_v2$sensitivities, specificity = roc_v2$specificities, run = i)
  auc_values_v2[i] <- pROC::auc(roc_v2)[1]
}

# Calculate the average ROC across all folds
print(paste("Average K-fold ROC for the interaction nrevss model:", round(mean(auc_values_v1), 3)))
print(paste("Average K-fold ROC for the interaction nrevss model:", round(mean(auc_values_v2), 3)))

# Append ROC data for plotting
roc_data_v1 <- do.call(rbind, roc_curve_v1)
roc_data_v2 <- do.call(rbind, roc_curve_v2)

# Plot
auc_v1 <- round(mean(auc_values_v1), 3)
nrevss_roc_v1 <- ggplot(roc_data_v1, aes(x = 1 - specificity, y = sensitivity, group = as.factor(run))) +
  geom_line(linewidth = 2, alpha = 0.5, color = '#E76BF3') + # Use geom_line or geom_step
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1.2) +
  labs(
    title = "NREVSS Model",
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme_bw() + line_theme +
  coord_equal() +
  annotate("text", x = 0.75, y = 0.25, label = glue::glue("AUC = {auc_v1}"), color = "black", size = 5)
nrevss_roc_v1

auc_v2 <- round(mean(auc_values_v2), 3)
nrevss_roc_v2 <- ggplot(roc_data_v2, aes(x = 1 - specificity, y = sensitivity, group = as.factor(run))) +
  geom_line(linewidth = 2, alpha = 0.5, color = '#E76BF3') + # Use geom_line or geom_step
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1.2) +
  labs(
    title = "NREVSS Model",
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme_bw() + line_theme +
  coord_equal() +
  annotate("text", x = 0.75, y = 0.25, label = glue::glue("AUC = {auc_v2}"), color = "black", size = 5)
nrevss_roc_v2

# Save fig 
library(cowplot)

figure_1 <- cowplot::plot_grid(fever_roc_v1, symptom_roc_v1, interaction_roc_v1, case_roc_v1, nrevss_roc_v1,
                               fever_roc_v2, symptom_roc_v2, interaction_roc_v2, case_roc_v2, nrevss_roc_v2,
                               nrow = 2,
                               labels = c("a", "", "", "", "",
                                          "b", "", "", "", ""),
                               label_size = 20)

saveRDS(figure_1, './tmp/roc_by_model_object_flu.rds')

# Save figure 
ggsave('./figs/roc_by_model_flu.jpg', plot = figure_1, height = 8, width = 25)

#############################
# 5. RUN FULL MODELS & SAVE #
#############################

library(lme4)


# Fever model
fever_model_v1 <- glm(flu ~ fever, 
                   data = flu_v1,
                   family = binomial(link = "logit"))

fever_model_v2 <- glm(flu ~ fever, 
                   data = flu_v2,
                   family = binomial(link = "logit"))

# All symptom model
all_symptom_model_v1 <- glm(flu ~ fever*as.factor(state_fips) + 
                              cough*as.factor(state_fips) + 
                              sore_throat*as.factor(state_fips) +
                              myalgia*as.factor(state_fips) + 
                              short_breath*as.factor(state_fips) + 
                              hypoxemia*as.factor(state_fips) + 
                              nausea_vom*as.factor(state_fips) + 
                              bronchitis*as.factor(state_fips) + 
                              chest_pain*as.factor(state_fips) + 
                              diarrhea*as.factor(state_fips) + 
                              fatigue*as.factor(state_fips) + 
                              headache*as.factor(state_fips) + 
                              congestion*as.factor(state_fips) + 
                              sneezing*as.factor(state_fips), 
                            data = flu_v1,
                            family = binomial(link = "logit"))
all_symptom_model_v2 <- glm(flu ~ fever + 
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
                            data = flu_v2,
                            family = binomial(link = "logit"))

# Symptom interaction model
symptom_interaction_model_v1 <- glm(flu ~ fever*as.factor(age_grp) + 
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
                                    data = flu_v1,
                                    family = binomial(link = "logit"))
symptom_interaction_model_v2 <- glm(flu ~ fever*as.factor(age_grp) + 
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
                                    data = flu_v2,
                                    family = binomial(link = "logit"))

# Interaction case model
interaction_case_model_v1 <- gam(flu ~ fever*as.factor(age_grp) + 
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
                                   s(visit_pos, bs = "cc", k = 20), 
                                 data = flu_v1,
                                 family = binomial(link = "logit"))
interaction_case_model_v2 <- gam(flu ~ fever*as.factor(age_grp) + 
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
                                   s(visit_pos, bs = "cc", k = 20), 
                                 data = flu_v2,
                                 family = binomial(link = "logit"))

# Interaction nrevss model
interaction_nrevss_model_v1 <- gam(flu ~ fever*as.factor(age_grp) + 
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
                                     s(nrevss_pos, bs = "cc", k = 20), 
                                   data = flu_v1,
                                   family = binomial(link = "logit"))
interaction_nrevss_model_v2 <- gam(flu ~ fever*as.factor(age_grp) + 
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
                                     s(nrevss_pos, bs = "cc", k = 20), 
                                   data = flu_v2,
                                   family = binomial(link = "logit"))

# Save models
saveRDS(fever_model_v1, "./tmp/fever_model_v1_flu.rds")
saveRDS(all_symptom_model_v1, "./tmp/all_symptom_model_v1_flu.rds")
saveRDS(symptom_interaction_model_v1, "./tmp/symptom_interaction_model_v1_flu.rds")
saveRDS(interaction_case_model_v1, "./tmp/interaction_case_model_v1_flu.rds")
saveRDS(interaction_nrevss_model_v1, "./tmp/interaction_nrevss_model_v1_flu.rds")
saveRDS(fever_model_v2, "./tmp/fever_model_v2_flu.rds")
saveRDS(all_symptom_model_v2, "./tmp/all_symptom_model_v2_flu.rds")
saveRDS(symptom_interaction_model_v2, "./tmp/symptom_interaction_model_v2_flu.rds")
saveRDS(interaction_case_model_v2, "./tmp/interaction_case_model_v2_flu.rds")
saveRDS(interaction_nrevss_model_v2, "./tmp/interaction_nrevss_model_v2_flu.rds")

##############################
# 6. DETERMINE CUTOFF VALUES #
##############################

# Predict based on models
flu_v1$pred_fever <- predict(fever_model_v1, newdata = flu_v1, type = "response")
flu_v1$pred_symptom <- predict(all_symptom_model_v1, newdata = flu_v1, type = "response")
flu_v1$pred_interaction <- predict(symptom_interaction_model_v1, newdata = flu_v1, type = "response")
flu_v1$pred_case <- predict(interaction_case_model_v1, newdata = flu_v1, type = "response")
flu_v1$pred_nrevss <- predict(interaction_nrevss_model_v1, newdata = flu_v1, type = "response")
flu_v2$pred_fever <- predict(fever_model_v2, newdata = flu_v2, type = "response")
flu_v2$pred_symptom <- predict(all_symptom_model_v2, newdata = flu_v2, type = "response")
flu_v2$pred_interaction <- predict(symptom_interaction_model_v2, newdata = flu_v2, type = "response")
flu_v2$pred_case <- predict(interaction_case_model_v2, newdata = flu_v2, type = "response")
flu_v2$pred_nrevss <- predict(interaction_nrevss_model_v2, newdata = flu_v2, type = "response")

# Determine cutpoints by age group
cutpoints_v1 <- flu_v1 |> group_by(age_grp) |>
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
  mutate(cutoff_nrevss = cutpointr(pred_nrevss, flu, 
                                 method = maximize_metric, 
                                 metric = youden, 
                                 na.rm = TRUE)$optimal_cutpoint) |>
  distinct(age_grp, cutoff_fever, cutoff_symptom, cutoff_interaction, cutoff_case, cutoff_nrevss)

cutpoints_v2 <- flu_v2 |> group_by(age_grp) |>
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
  mutate(cutoff_nrevss = cutpointr(pred_nrevss, flu, 
                                   method = maximize_metric, 
                                   metric = youden, 
                                   na.rm = TRUE)$optimal_cutpoint) |>
  distinct(age_grp, cutoff_fever, cutoff_symptom, cutoff_interaction, cutoff_case, cutoff_nrevss)


# Save
saveRDS(cutpoints_v1, "./tmp/model_cutpoints_v1_flu.rds")
saveRDS(cutpoints_v2, "./tmp/model_cutpoints_v2_flu.rds")

################################################################################
################################################################################
