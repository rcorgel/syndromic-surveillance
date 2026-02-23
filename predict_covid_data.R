################################################################################
# File Name: predict_covid_data                                                #
#                                                                              #
# Purpose:   Predict COVID-19 from medical claims symptoms.                    #
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

#################
# 2. COVID DATA #
#################

# Combine covid cases and non-cases
covid_cases_ind <- readRDS('./tmp/cases_exp_sub_delta_covid.rds')
covid_cases_ind <- covid_cases_ind |>  mutate(symp_count = fever + cough + 
                                            sore_throat + short_breath + hypoxemia + loss_smell_taste +
                                            loss_appetite + bronchitis + fatigue + headache + nausea_vom + 
                                            diarrhea + congestion + chest_pain,
                                          no_symptoms = ifelse(symp_count == 0, 1, 0))
covid_non_cases_ind <- readRDS('./tmp/non_cases_exp_sub_delta_covid.rds')
covid_non_cases_ind <- covid_non_cases_ind |>  mutate(symp_count = fever + cough + 
                                                        sore_throat + short_breath + hypoxemia + loss_smell_taste +
                                                        loss_appetite + bronchitis + fatigue + headache + nausea_vom + 
                                                        diarrhea + congestion + chest_pain,
                                                  no_symptoms = ifelse(symp_count == 0, 1, 0))
covid_subset_ind <- rbind(covid_cases_ind, covid_non_cases_ind)
remove(covid_cases_ind, covid_non_cases_ind)

# Create training and testing data sets
sample <- sample(c(TRUE, FALSE), nrow(covid_subset_ind), replace=TRUE, prob=c(0.5,0.5))
covid_train  <- covid_subset_ind[sample, ]
covid_test   <- covid_subset_ind[!sample, ]
remove(covid_subset_ind)

# Collapse train data to the group level
covid_train$num <- 1
covid_train_count <- covid_train |> group_by(covid, fever, cough, sore_throat, short_breath, 
                                             hypoxemia, loss_smell_taste, loss_appetite, bronchitis,
                                             fatigue, headache, nausea_vom, diarrhea, congestion, chest_pain, 
                                             year, month, state_fips, 
                                         county_fips, no_symptoms, symp_count, age_grp, 
                                         patient_gender_code, year_month) |>
  mutate(count = sum(num)) |>
  distinct(covid, fever, cough, sore_throat, short_breath, 
           hypoxemia, loss_smell_taste, loss_appetite, bronchitis,
           fatigue, headache, nausea_vom, diarrhea, congestion, chest_pain, 
           year, month, state_fips,
           county_fips, no_symptoms, symp_count, count, age_grp, 
           patient_gender_code, year_month, .keep_all = FALSE)

# Collapse test data to the group level
covid_test$num <- 1
covid_test_count <- covid_test |> group_by(covid, fever, cough, sore_throat, short_breath, 
                                           hypoxemia, loss_smell_taste, loss_appetite, bronchitis,
                                           fatigue, headache, nausea_vom, diarrhea, congestion, chest_pain, 
                                           year, month, state_fips, 
                                           county_fips, no_symptoms, symp_count, age_grp, 
                                           patient_gender_code, year_month) |>
  mutate(count = sum(num)) |>
  distinct(covid, fever, cough, sore_throat, short_breath, 
           hypoxemia, loss_smell_taste, loss_appetite, bronchitis,
           fatigue, headache, nausea_vom, diarrhea, congestion, chest_pain, 
           year, month, state_fips, 
           county_fips, no_symptoms, symp_count, count, age_grp, 
           patient_gender_code, year_month, .keep_all = FALSE)

# Perform simple logistic regression
train_model <- glm(covid ~ fever + cough + 
                     sore_throat + short_breath + hypoxemia + loss_smell_taste +
                     loss_appetite + bronchitis + fatigue + headache + nausea_vom + 
                     diarrhea + congestion + chest_pain + as.factor(state_fips) +
                     as.factor(age_grp) + as.factor(patient_gender_code),
                   family = binomial(link = "logit"),
                   data = covid_train_count,
                   weights = count)
summary(train_model)


library(car)
vif(train_model)
#exp(coef(train_model))

# Predict testing data
covid_test_count$pred <- predict(train_model, newdata = covid_test_count, type = "response")

# Display ROC curve
basicplot <- ggplot(covid_test_count, aes(d = covid, m = pred)) + 
  geom_roc(n.cuts=20,labels=FALSE) + 
  style_roc(theme = theme_grey)
basicplot + theme_minimal() +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed')

roc_object <- pROC::roc(covid_test_count$covid, covid_test_count$pred)

# Calculate area under curve
pROC::auc(roc_object)






























# Create covid data format function
covid_data_format <- function(data) {
  # Replace True/False with 1 and 0
  data$covid <- ifelse(data$covid == 'True', 1, 0)
  data$fever <- ifelse(data$fever == 'True', 1, 0)
  data$loss_smell <- ifelse(data$loss_smell == 'True', 1, 0)
  data$loss_app <- ifelse(data$loss_app == 'True', 1, 0)
  data$myalgia <- ifelse(data$myalgia == 'True', 1, 0)
  data$hypoxemia <- ifelse(data$hypoxemia == 'True', 1, 0)
  data$short_breath <- ifelse(data$short_breath == 'True', 1, 0)
  data$cough <- ifelse(data$cough == 'True', 1, 0)
  data$chest_pain <- ifelse(data$chest_pain == 'True', 1, 0)
  data$nausea <- ifelse(data$nausea == 'True', 1, 0)
  data$sore_throat <- ifelse(data$sore_throat == 'True', 1, 0)
  data$fatigue <- ifelse(data$fatigue == 'True', 1, 0)
  
  # Make state and year variables
  data$state_fips <- substr(data$county_fips, 1, 2)
  data$year <- year(as.Date(paste(data$year_week, 1, sep = '-'), format = '%Y-%U-%u'))
  data$year <- ifelse(is.na(data$year), 2020, data$year) # throws errors for week 53
  data$month <- month(as.Date(paste(data$year_week, 1, sep = '-'), format = '%Y-%U-%u'))
  data$month <- ifelse(is.na(data$month), 12, data$month) # throws errors for week 53
  
  # Create testing and training data
  sample <- sample(c(TRUE, FALSE), nrow(data), replace=TRUE, prob=c(0.5,0.5))
  covid_train  <- data[sample, ]
  covid_test   <- data[!sample, ]
  
  # Remove data
  remove(data)
  
  # Collapse data
  # Collapse train data to the group level
  covid_train$num <- 1
  covid_train_count <- covid_train |> group_by(covid, fever, loss_app, loss_smell, 
                                               myalgia, hypoxemia, short_breath, cough, 
                                               chest_pain, nausea, sore_throat, fatigue,
                                               year, month, state_fips, 
                                               county_fips, year_week) |>
    mutate(count = sum(num)) |>
    distinct(covid, fever, loss_app, loss_smell, 
             myalgia, hypoxemia, short_breath, cough, 
             chest_pain, nausea, sore_throat, fatigue, 
             year, month, state_fips, 
             county_fips, year_week, count)
  
  # Collapse test data to the group level
  covid_test$num <- 1
  covid_test_count <- covid_test |> group_by(covid, fever, loss_app, loss_smell, 
                                             myalgia, hypoxemia, short_breath, cough, 
                                             chest_pain, nausea, sore_throat, fatigue,
                                             year, month, state_fips, 
                                             county_fips, year_week) |>
    mutate(count = sum(num)) |>
    distinct(covid, fever, loss_app, loss_smell, 
             myalgia, hypoxemia, short_breath, cough, 
             chest_pain, nausea, sore_throat, fatigue, 
             year, month, state_fips, 
             county_fips, year_week, count)
  
  # Remove data
  remove(covid_train, covid_test)
  
  # Return collapsed train and test data
  covid_data <- list(covid_train_count, covid_test_count)
  return(covid_data)
}

# Loop through data to format, split into testing/training, and collapsing
covid_data_list <- NULL
for (i in c('./tmp/covid_cases_full.rds', 
            './tmp/covid_non_cases_sample_2020.rds',
            './tmp/covid_non_cases_sample_2021.rds',
            './tmp/covid_non_cases_sample_2022.rds')) {
  # Display and load data
  print(i)
  data <- readRDS(i)
  
  # Format data
  data_format <- covid_data_format(data)
  
  # Add to data list
  covid_data_list <- append(covid_data_list, data_format)
}

# Append data together
covid_train_count <- rbind(covid_data_list[[1]],
                           covid_data_list[[3]],
                           covid_data_list[[5]],
                           covid_data_list[[7]])
covid_test_count <- rbind(covid_data_list[[2]],
                          covid_data_list[[4]],
                          covid_data_list[[6]],
                          covid_data_list[[8]])
remove(covid_data_list)
# Simple logistic regression
train_model <- glm(covid ~ fever + loss_app + loss_smell + myalgia + hypoxemia + 
                     short_breath + cough + chest_pain + nausea + sore_throat + 
                     fatigue + as.factor(state_fips) + as.factor(year) + as.factor(month),
                   family = binomial(link = "logit"),
                   data = covid_train_count,
                   weights = count)
summary(train_model)

# Predict testing data
covid_test_count$pred <- predict(train_model, newdata = covid_test_count, type = "response")

# Display ROC curve
basicplot <- ggplot(covid_test_count, aes(d = covid, m = pred)) + 
  geom_roc(n.cuts=20,labels=FALSE) + 
  style_roc(theme = theme_grey)
basicplot

roc_object <- pROC::roc(covid_test_count$covid, covid_test_count$pred)

# Calculate area under curve
pROC::auc(roc_object)

# Save model
saveRDS(train_model, file = './tmp/covid_19_train_model.rds')





remove(roc_object, basicplot)

cutoff <- NULL
value <- NULL 
value_2 <- NULL
for (i in seq(0, 1, 0.01)) {
  print(i)
  covid_test_count$pos <- ifelse(covid_test_count$pred >= i, 1, 0)
  covid_test_count$correct <- 
    ifelse(covid_test_count$pos == covid_test_count$covid, 1, 0) * covid_test_count$count
  accuracy <- sum(covid_test_count$correct) / sum(covid_test_count$count)
  all_covid <- sum(covid_test_count$pos * covid_test_count$count) / 
    sum(covid_test_count$covid * covid_test_count$count)
  cutoff <- rbind(cutoff, i)
  value <- rbind(value, accuracy)
  value_2 <- rbind(value_2, all_covid)
}

data <- data.frame(cutoff, value, value_2)
ggplot(data) + geom_line(aes(x = cutoff, y = value))

library(cutpointr)
opt_cut <- cutpointr(covid_test_count, pred, covid, direction = ">=", pos_class = 1,
                     neg_class = 0, method = maximize_metric, metric = youden)
plot(opt_cut)
summary(opt_cut)


rm(list = ls())
train_model <- readRDS('./tmp/covid_train_model.rds')
covid <- readRDS('./tmp/covid_full.rds')

covid$state_fips <- substr(covid$county_fips, 1, 2)
covid$year <- year(as.Date(paste(covid$year_week, 1, sep = '-'), format = '%Y-%U-%u'))
covid$month <- month(as.Date(paste(covid$year_week, 1, sep = '-'), format = '%Y-%U-%u'))

covid <- covid %>% filter(state_fips != '99')

covid$covid <- ifelse(covid$covid == 'True', 1, 0)
covid$fever <- ifelse(covid$fever == 'True', 1, 0)
covid$myalgia <- ifelse(covid$myalgia == 'True', 1, 0)
covid$short_breath <- ifelse(covid$short_breath == 'True', 1, 0)
covid$cough <- ifelse(covid$cough == 'True', 1, 0)
covid$nausea <- ifelse(covid$nausea == 'True', 1, 0)
covid$sore_throat <- ifelse(covid$sore_throat == 'True', 1, 0)
covid$loss_app <- ifelse(covid$loss_app == 'True', 1, 0)
covid$fatigue <- ifelse(covid$fatigue == 'True', 1, 0)
covid$loss_smell <- ifelse(covid$loss_smell == 'True', 1, 0)
covid$hypoxemia <- ifelse(covid$hypoxemia == 'True', 1, 0)
covid$chest_pain <- ifelse(covid$chest_pain == 'True', 1, 0)

covid <- covid |> filter(year != 2023)

covid$pred <- predict(train_model, newdata = covid, type = "response")
covid$covid_pred <- ifelse(covid$pred >= 0.7147, 1, 0)

covid_time <- covid %>% group_by(year_week) %>%
  mutate(sum_covid = sum(covid * patient_count_imp),
         sum_pred = sum(covid_pred * patient_count_imp),
         year_week_date = as.Date(paste(year_week, 1, sep = '-'), format = '%Y-%U-%u')) %>%
  distinct(year_week_date, sum_covid, sum_pred)
ggplot(covid_time) + geom_line(aes(x = year_week_date, y = sum_covid), linewidth = 1.2, color = 'blue') +
  geom_line(aes(x = year_week_date, y = sum_pred), linewidth = 1.2, color = 'green') + theme_minimal()




