################################################################################
# File Name: predict_rsv_data                                                  #
#                                                                              #
# Purpose:   Predict RSV from medical claims symptoms.                         #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load data                                                      #
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
library(plotROC)
library(pROC)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

###############
# 2. RSV DATA #
###############

# Combine flu cases and non-cases
rsv_cases_ind <- readRDS('./tmp/cases_exp_sub_18_19_rsv.rds')
rsv_cases_ind <- rsv_cases_ind |>  mutate(symp_count = fever + cough + 
                                            sore_throat + short_breath + hypoxemia + 
                                            loss_appetite + bronchitis + fatigue + headache + 
                                            congestion + sneezing,
                                          no_symptoms = ifelse(symp_count == 0, 1, 0))
rsv_non_cases_ind <- readRDS('./tmp/non_cases_exp_sub_18_19_rsv.rds')
rsv_non_cases_ind <- rsv_non_cases_ind |>  mutate(symp_count = fever + cough + 
                                                    sore_throat + short_breath + hypoxemia + 
                                                    loss_appetite + bronchitis + fatigue + headache + 
                                                    congestion + sneezing,
                                                  no_symptoms = ifelse(symp_count == 0, 1, 0))
rsv_subset_ind <- rbind(rsv_cases_ind, rsv_non_cases_ind)
remove(rsv_cases_ind, rsv_non_cases_ind)

# Create training and testing data sets
sample <- sample(c(TRUE, FALSE), nrow(rsv_subset_ind), replace=TRUE, prob=c(0.5,0.5))
rsv_train  <- rsv_subset_ind[sample, ]
rsv_test   <- rsv_subset_ind[!sample, ]
remove(rsv_subset_ind)

# Collapse train data to the group level
rsv_train$num <- 1
rsv_train_count <- rsv_train |> group_by(rsv, fever, short_breath, hypoxemia, 
                                         cough, bronchitis, loss_appetite,
                                         sore_throat, fatigue, headache, 
                                         congestion, sneezing, year, month, state_fips, 
                                         county_fips, no_symptoms, symp_count, age_grp, 
                                         patient_gender_code, year_month) |>
  mutate(count = sum(num)) |>
  distinct(rsv, fever, short_breath, hypoxemia, 
           cough, bronchitis, loss_appetite,
           sore_throat, fatigue, headache, 
           congestion, sneezing, year, month, state_fips, 
           county_fips, no_symptoms, symp_count, age_grp, 
           patient_gender_code,count, year_month, .keep_all = FALSE)

# Collapse test data to the group level
rsv_test$num <- 1
rsv_test_count <- rsv_test |> group_by(rsv, fever, short_breath, hypoxemia, 
                                       cough, bronchitis, loss_appetite,
                                       sore_throat, fatigue, headache, 
                                       congestion, sneezing, year, month, state_fips, 
                                       county_fips, no_symptoms, symp_count, age_grp, 
                                       patient_gender_code, year_month) |>
  mutate(count = sum(num)) |>
  distinct(rsv, fever, short_breath, hypoxemia, 
           cough, bronchitis, loss_appetite,
           sore_throat, fatigue, headache, 
           congestion, sneezing, year, month, state_fips, 
           county_fips, no_symptoms, symp_count, age_grp, 
           patient_gender_code,count, year_month, .keep_all = FALSE)

# Perform simple logistic regression
train_model <- glm(rsv ~ fever + short_breath + hypoxemia +
                     cough + bronchitis + loss_appetite + 
                     sore_throat + fatigue + headache + 
                     congestion + sneezing + as.factor(state_fips) +
                     as.factor(age_grp) + as.factor(patient_gender_code),
                   family = binomial(link = "logit"),
                   data = rsv_train_count,
                   weights = count)
summary(train_model)


library(car)
vif(train_model)
#exp(coef(train_model))

# Predict testing data
rsv_test_count$pred <- predict(train_model, newdata = rsv_test_count, type = "response")

# Display ROC curve
basicplot <- ggplot(rsv_test_count, aes(d = rsv, m = pred)) + 
  geom_roc(n.cuts=20,labels=FALSE) + 
  style_roc(theme = theme_grey)
basicplot + theme_minimal() +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed')

roc_object <- pROC::roc(rsv_test_count$rsv, rsv_test_count$pred)

# Calculate area under curve
pROC::auc(roc_object)
