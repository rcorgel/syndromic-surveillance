################################################################################
# File Name: predict_flu_data                                                  #
#                                                                              #
# Purpose:   Predict influenza from medical claims symptoms.                   #
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
library(car)
library(cutpointr)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

###############
# 2. FLU DATA #
###############

# Combine flu cases and non-cases for each season
# 2016-17
flu_cases_ind <- readRDS('./tmp/cases_exp_sub_16_17_flu_symp.rds')
flu_non_cases_ind <- readRDS('./tmp/non_cases_exp_sub_16_17_flu_symp.rds')
flu_subset_ind_1 <- rbind(flu_cases_ind, flu_non_cases_ind) #readRDS('./tmp/flu_exp_sub_16_17_flu_symp.rds')
# 2017-18
flu_cases_ind <- readRDS('./tmp/cases_exp_sub_17_18_flu_symp.rds')
flu_non_cases_ind <- readRDS('./tmp/non_cases_exp_sub_17_18_flu_symp.rds')
flu_subset_ind_2 <- rbind(flu_cases_ind, flu_non_cases_ind) #readRDS('./tmp/flu_exp_sub_17_18_flu_symp.rds')
# 2018-19
flu_cases_ind <- readRDS('./tmp/cases_exp_sub_18_19_flu_symp.rds')
flu_non_cases_ind <- readRDS('./tmp/non_cases_exp_sub_18_19_flu_symp.rds')
flu_subset_ind_3 <- rbind(flu_cases_ind, flu_non_cases_ind) #readRDS('./tmp/flu_exp_sub_18_19_flu_symp.rds')
# 2019-20
flu_cases_ind <- readRDS('./tmp/cases_exp_sub_19_20_flu_symp.rds')
flu_non_cases_ind <- readRDS('./tmp/non_cases_exp_sub_19_20_flu_symp.rds')
flu_subset_ind_4 <- rbind(flu_cases_ind, flu_non_cases_ind) #readRDS('./tmp/flu_exp_sub_19_20_flu_symp.rds')
remove(flu_cases_ind, flu_non_cases_ind)
# Combine all dat ainto a list
flu_subset_ind <- list(flu_subset_ind_1, flu_subset_ind_2, flu_subset_ind_3, flu_subset_ind_4)

# Set empty objects of all of the training models and testing data
training_models <- NULL
testing_data <- NULL
test_dat <- readRDS('./tmp/flu_p2_dat_proc_select_counties.rds')
test_dat$count <- test_dat$patient_count_imp

# Loop through all seasons
for (i in seq(1, 4, 1)) {
  print(i)
  # Create symptom count variable
  flu_subset_ind_lvl <- flu_subset_ind[[i]] |>  mutate(ili = ifelse((fever == 1 & cough ==1) |  (fever == 1 & sore_throat == 1), 1, 0),
                                                       fever_sore_throat = ifelse((fever == 1 & sore_throat == 1), 1, 0),
                                                       fever_sore_throat_cough = ifelse((fever == 1 & sore_throat == 1 & cough == 1), 1, 0),
                                                       fever_cough = ifelse((fever == 1 & cough == 1), 1, 0),
                                                       fever_nausea = ifelse((fever == 1 & nausea_vom == 1), 1, 0),
                                                       fever_nausea_cough = ifelse((fever == 1 & cough == 1 & nausea_vom == 1), 1, 0),
                                                       fever_short_breath = ifelse((fever == 1 & short_breath == 1), 1, 0),
                                                       cough_short_breath = ifelse((cough == 1 & short_breath == 1), 1, 0),
                                                       cough_congestion = ifelse((cough == 1 & congestion == 1), 1, 0),
                                                       fever_congestion = ifelse((fever == 1 & congestion == 1), 1, 0),
                                                       fever_cough_congestion = ifelse((fever == 1 & cough == 1 & congestion == 1), 1, 0),
                                                       cough_myalgia = ifelse((cough == 1 & myalgia == 1), 1, 0),
                                                       fever_myalgia = ifelse((fever == 1 & myalgia == 1), 1, 0),
                                                       fever_cough_myalgia = ifelse((fever == 1 & cough == 1 & myalgia == 1), 1, 0),
                                                       summer = ifelse((month == 6 | month == 7 | month == 8), 1, 0))

  # Create training and testing data sets
  sample <- sample(c(TRUE, FALSE), nrow(flu_subset_ind_lvl), replace=TRUE, prob=c(0.5, 0.5))
  flu_train  <- flu_subset_ind_lvl[sample, ]
  flu_test   <- flu_subset_ind_lvl[!sample, ]
  remove(flu_subset_ind_lvl)
  
  # Collapse train data to the group level
  flu_train$num <- 1
  flu_train_count <- flu_train |> group_by(flu, ili, fever, myalgia, short_breath, hypoxemia, 
                                           cough, nausea_vom, bronchitis, chest_pain,
                                           sore_throat, diarrhea, fatigue, headache, 
                                           congestion, sneezing, year, month, state_fips, 
                                           county_fips, symp_count, age_grp, 
                                           patient_gender_code, month_date, fever_sore_throat, fever_cough, fever_sore_throat_cough,
                                           fever_nausea, fever_nausea_cough, fever_short_breath, cough_short_breath, 
                                           cough_congestion, fever_congestion, fever_cough_congestion,
                                           cough_myalgia, fever_myalgia, fever_cough_myalgia, 
                                           cough_perc_lag, fever_perc_lag, sore_throat_perc_lag, summer) |>
    mutate(count = sum(num)) |>
    distinct(flu, ili, fever, myalgia, short_breath, hypoxemia, 
             cough, nausea_vom, bronchitis, chest_pain,
             sore_throat, diarrhea, fatigue, headache, 
             congestion, sneezing, year, month, state_fips, 
             county_fips, symp_count, age_grp, 
             patient_gender_code,count, month_date, fever_sore_throat, fever_cough, fever_sore_throat_cough,
             fever_nausea, fever_nausea_cough, fever_short_breath, cough_short_breath, 
             cough_congestion, fever_congestion, fever_cough_congestion,
             cough_myalgia, fever_myalgia, fever_cough_myalgia, cough_perc_lag, fever_perc_lag, sore_throat_perc_lag, summer, .keep_all = FALSE)

  # Collapse test data to the group level
  flu_test$num <- 1
  flu_test_count <- flu_test |> group_by(flu, ili, fever, myalgia, short_breath, hypoxemia, 
                                         cough, nausea_vom, bronchitis, chest_pain,
                                         sore_throat, diarrhea, fatigue, headache, 
                                         congestion, sneezing, year, month, state_fips, 
                                         county_fips, symp_count, age_grp, 
                                         patient_gender_code, month_date, fever_sore_throat, fever_cough, fever_sore_throat_cough,
                                         fever_nausea, fever_nausea_cough, fever_short_breath, cough_short_breath, 
                                         cough_congestion, fever_congestion, fever_cough_congestion,
                                         cough_myalgia, fever_myalgia, fever_cough_myalgia, cough_perc_lag, fever_perc_lag, sore_throat_perc_lag, summer) |>
    mutate(count = sum(num)) |>
    distinct(flu, ili, fever, myalgia, short_breath, hypoxemia, 
             cough, nausea_vom, bronchitis, chest_pain,
             sore_throat, diarrhea, fatigue, headache, 
             congestion, sneezing, year, month, state_fips, 
             county_fips, symp_count, age_grp, 
             patient_gender_code,count, month_date, fever_sore_throat, fever_cough, fever_sore_throat_cough,
             fever_nausea, fever_nausea_cough, fever_short_breath, cough_short_breath, 
             cough_congestion, fever_congestion, fever_cough_congestion,
             cough_myalgia, fever_myalgia, fever_cough_myalgia, cough_perc_lag, fever_perc_lag, sore_throat_perc_lag, summer, .keep_all = FALSE)

  # Perform simple logistic regression
  train_model <- glm(flu ~ fever*cough*sore_throat + myalgia + short_breath + hypoxemia +
                     nausea_vom + bronchitis + chest_pain + 
                     diarrhea + fatigue + headache + 
                     congestion + sneezing - 0 + as.factor(month),
                     family = binomial(link = "logit"),
                     data = flu_train_count,
                     weights = count)
  
  # Append model to object
  training_models[[i]] <- train_model
  
  # Predict testing data
  test_dat$pred <- predict(train_model, newdata = test_dat, type = "response")
  
  # Append data to object
  testing_data[[i]] <- test_dat
}

summary(training_models[[1]])
vif(training_models[[1]])

basicplot <- ggplot(testing_data[[1]], aes(d = flu, m = pred)) + 
  geom_roc(n.cuts=40,labels=FALSE) + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  theme_minimal() + ylab('True Positive Rate') + xlab('False Positive Rate') +
  theme_minimal() + theme(
                          axis.text.x = element_text(size = 12),
                          axis.text.y = element_text(size = 12),
                          axis.title.x = element_text(size = 14),
                          axis.title.y = element_text(size = 14),
                          strip.text = element_text(size=14),
                          plot.title = element_text(size=18),
                          legend.text=element_text(size=12),
                          legend.title=element_text(size=14))
basicplot 

basicplot <- ggplot(testing_data[[2]], aes(d = flu, m = pred)) + 
  geom_roc(n.cuts=40,labels=FALSE) + 
  style_roc(theme = theme_grey) + theme_minimal() +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed')
basicplot 

basicplot <- ggplot(testing_data[[3]], aes(d = flu, m = pred)) + 
  geom_roc(n.cuts=40,labels=FALSE) + 
  style_roc(theme = theme_grey) + theme_minimal() +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed')
basicplot 



# 
# summary(training_models[[1]])
# 
# exp(coef(training_models[[1]]))
# 
# 1 / (1 + exp(-1 * (-1.357075 + 0.476540 + 0.702647)))
# 
# exp(coef(training_models[[1]])) / (1 + exp(coef(training_models[[1]])))

# Calculate AUC
# 2016-17
pROC::auc(pROC::roc(testing_data[[1]]$flu, testing_data[[1]]$pred))
# 2017-18
pROC::auc(pROC::roc(testing_data[[2]]$flu, testing_data[[2]]$pred))
# 2018-19
pROC::auc(pROC::roc(testing_data[[3]]$flu, testing_data[[3]]$pred))
# 2019-20
pROC::auc(pROC::roc(testing_data[[4]]$flu, testing_data[[4]]$pred))

# Evaluate ideal cut points
# 2016-17
# Youden Index
opt_cut_1617_youden <- cutpointr(testing_data[[1]], pred, flu, direction = ">=", 
                                 pos_class = 1, neg_class = 0, 
                                 method = maximize_metric, metric = youden, na.rm = TRUE)
summary(opt_cut_1617_youden)
# Accuracy
opt_cut_1617_accuracy <- cutpointr(testing_data[[1]], pred, flu, direction = ">=", 
                                 pos_class = 1, neg_class = 0, 
                                 method = maximize_metric, metric = accuracy, na.rm = TRUE)
summary(opt_cut_1617_accuracy)
# 2017-18
# Youden Index
opt_cut_1718_youden <- cutpointr(testing_data[[2]], pred, flu, direction = ">=", 
                                 pos_class = 1, neg_class = 0, 
                                 method = maximize_metric, metric = youden, na.rm = TRUE)
summary(opt_cut_1718_youden)
# Accuracy
opt_cut_1718_accuracy <- cutpointr(testing_data[[2]], pred, flu, direction = ">=", 
                                   pos_class = 1, neg_class = 0, 
                                   method = maximize_metric, metric = accuracy, na.rm = TRUE)
summary(opt_cut_1718_accuracy)
# 2018-19
# Youden Index
opt_cut_1819_youden <- cutpointr(testing_data[[3]], pred, flu, direction = ">=", 
                                 pos_class = 1, neg_class = 0, 
                                 method = maximize_metric, metric = youden, na.rm = TRUE)
summary(opt_cut_1819_youden)
# Accuracy
opt_cut_1819_accuracy <- cutpointr(testing_data[[3]], pred, flu, direction = ">=", 
                                   pos_class = 1, neg_class = 0, 
                                   method = maximize_metric, metric = accuracy, na.rm = TRUE)
summary(opt_cut_1819_accuracy)
# 2019-20
# Youden Index
opt_cut_1920_youden <- cutpointr(testing_data[[4]], pred, flu, direction = ">=", 
                                 pos_class = 1, neg_class = 0, 
                                 method = maximize_metric, metric = youden, na.rm = TRUE)
summary(opt_cut_1920_youden)
# Accuracy
opt_cut_1920_accuracy <- cutpointr(testing_data[[4]], pred, flu, direction = ">=", 
                                   pos_class = 1, neg_class = 0, 
                                   method = maximize_metric, metric = accuracy, na.rm = TRUE)
summary(opt_cut_1920_accuracy)

# Save models and cut off values
save(list = c('training_models',
              'opt_cut_1617_youden', 
              'opt_cut_1718_youden', 
              'opt_cut_1819_youden', 
              'opt_cut_1920_youden',
              'opt_cut_1617_accuracy',
              'opt_cut_1718_accuracy', 
              'opt_cut_1819_accuracy', 
              'opt_cut_1920_accuracy'), 
     file = './tmp/flu_model_cut_off_results.RData')


saveRDS(training_models, './tmp/training_models_complex.rds') 









# Display ROC curve
ggplot(testing_data[[2]], aes(d = flu, m = pred)) + 
  geom_roc(n.cuts=20,labels=FALSE) + 
  style_roc(theme = theme_grey) + theme_minimal() +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed')

# Calculate area under curve
pROC::auc(pROC::roc(testing_data[[4]]$flu, testing_data[[4]]$pred))

opt_cut <- cutpointr(flu_test_count, pred, flu, direction = ">=", pos_class = 1,
                     neg_class = 0, method = maximize_metric, metric = prod_sens_spec)
plot(opt_cut)
summary(opt_cut)









# Prediction
flu_train_count$pred <- predict(train_model, type = "response")
flu_train_uncount <- testing_data[[1]] |> uncount(count)
ggplot(flu_train_uncount, aes(pred, color = as.factor(flu))) +
  geom_density(size = 1) + theme_minimal()


cutoff <- NULL
value <- NULL 
value_2 <- NULL
for (i in seq(0, 1, 0.01)) {
  print(i)
  flu_test_count$pos <- ifelse(flu_test_count$pred >= i, 1, 0)
  flu_test_count$correct <- 
    ifelse(flu_test_count$pos == flu_test_count$flu, 1, 0) * flu_test_count$count
  accuracy <- sum(flu_test_count$correct) / sum(flu_test_count$count)
  all_flu <- sum(flu_test_count$pos * flu_test_count$count) / 
    sum(flu_test_count$flu * flu_test_count$count)
  cutoff <- rbind(cutoff, i)
  value <- rbind(value, accuracy)
  value_2 <- rbind(value_2, all_flu)
}

data <- data.frame(cutoff, value, value_2)
ggplot(data) + geom_line(aes(x = cutoff, y = value))










#0.7147 
rm(list = ls())
train_model <- readRDS('./tmp/flu_train_model.rds')
summary(train_model)
flu <- readRDS('./tmp/flu_full.rds')

flu$state_fips <- substr(flu$county_fips, 1, 2)
flu$year <- year(as.Date(paste(flu$year_week, 1, sep = '-'), format = '%Y-%U-%u'))
flu$month <- month(as.Date(paste(flu$year_week, 1, sep = '-'), format = '%Y-%U-%u'))

flu <- flu %>% filter(state_fips != '99')

flu$flu <- ifelse(flu$flu == 'True', 1, 0)
flu$fever <- ifelse(flu$fever == 'True', 1, 0)
flu$myalgia <- ifelse(flu$myalgia == 'True', 1, 0)
flu$short_breath <- ifelse(flu$short_breath == 'True', 1, 0)
flu$cough <- ifelse(flu$cough == 'True', 1, 0)
flu$nausea <- ifelse(flu$nausea == 'True', 1, 0)
flu$sore_throat <- ifelse(flu$sore_throat == 'True', 1, 0)
flu$diarrhea <- ifelse(flu$diarrhea == 'True', 1, 0)
flu$fatigue <- ifelse(flu$fatigue == 'True', 1, 0)
flu$headache <- ifelse(flu$headache == 'True', 1, 0)
flu$congestion <- ifelse(flu$congestion == 'True', 1, 0)
flu$sneezing <- ifelse(flu$sneezing == 'True', 1, 0)

flu_1 <- flu |> filter(year < 2018)
flu_2 <- flu |> filter(year >= 2018)
remove(flu)

flu_1$pred <- predict(train_model, newdata = flu_1, type = "response")
flu_1$flu_pred <- ifelse(flu_1$pred >= 0.90, 1, 0)

flu_2$pred <- predict(train_model, newdata = flu_2, type = "response")
flu_2$flu_pred <- ifelse(flu_2$pred >= 0.90, 1, 0)

flu <- rbind(flu_1, flu_2)
remove(flu_1, flu_2)

flu$ili <- ifelse(flu$fever == 1 & (flu$sore_throat == 1 | flu$cough == 1), 1, 0)

flu_time <- flu %>% group_by(year_week) %>%
  mutate(sum_flu = sum(flu * patient_count_imp),
         sum_pred = sum(flu_pred * patient_count_imp),
         sum_ili = sum(ili * patient_count_imp),
         year_week_date = as.Date(paste(year_week, 1, sep = '-'), format = '%Y-%U-%u')) %>%
  distinct(year_week_date, sum_flu, sum_pred, sum_ili)


flu_month <- flu %>% group_by(year, month) %>%
  mutate(sum_flu = sum(flu * patient_count_imp),
         sum_pred = sum(flu_pred * patient_count_imp),
         sum_ili = sum(ili * patient_count_imp),
         year_month_date = as.Date(paste(year, month, 1, sep = '-'), format = '%Y-%m-%d')) %>%
  distinct(year_month_date, sum_flu, sum_pred, sum_ili)

flu_year <- flu %>% group_by(year) %>%
  mutate(sum_flu = sum(flu * patient_count_imp),
         sum_pred = sum(flu_pred * patient_count_imp),
         sum_ili = sum(ili * patient_count_imp)) %>%
  distinct(year, sum_flu, sum_pred, sum_ili)

coeff <- max(flu_time$sum_flu) / max(flu_time$sum_pred)

ggplot(flu_time, aes(x=year_week_date)) +
  geom_line(aes(y = sum_pred), linewidth = 1.1, color = '#00ABFD') +
  geom_line(aes(y = sum_ili / coeff), linewidth = 1.1, color = '#00BA38') +
  geom_line(aes(y = sum_flu / coeff), linewidth = 0.6, color = 'black', alpha = 0.65, linetype = 2) + 
  scale_y_continuous(
    # Features of the first axis
    name = "Predicted Flu Cases",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name= "Observed Flu/ILI Cases")
  ) + theme_minimal() + xlab('Date (weekly)') + 
  theme(axis.title.y.left = element_text(color = "#00ABFD"))

ggplot(flu_month, aes(x=year_month_date)) +
  geom_line(aes(y = sum_pred), linewidth = 1.1, color = '#00ABFD') +
  geom_line(aes(y = sum_ili), linewidth = 1.1, color = '#00BA38') +
  geom_line(aes(y = sum_flu), linewidth = 0.6, color = 'black', alpha = 0.65, linetype = 2) + 
  theme_minimal() + xlab('Date (monthly)') + ylab('Number of Cases')



all_cause <- read.csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1CgyyYhIWFVG2aMITI-HBpKamvwwpQ_ED/CHC_data/all_cause/2023_08_22_county_weekly_flu_ac_export.csv')
all_cause$all_cause <- as.numeric(all_cause$all_cause)
all_cause$all_cause <- ifelse(is.na(all_cause$all_cause), 0, all_cause$all_cause)
all_cause$year_week_date = as.Date(paste(all_cause$year_week, 1, sep = '-'), format = '%Y-%U-%u')

all_cause_time <- all_cause %>% group_by(year_week_date) %>%
  mutate(sum_cause = sum(all_cause)) |>
  distinct(year_week_date, sum_cause)

flu_time <- left_join(flu_time, all_cause_time, by = c('year_week_date' = 'year_week_date'))
flu_time$obs_prop <- flu_time$sum_flu / flu_time$sum_cause
flu_time$pred_prop <- flu_time$sum_pred / flu_time$sum_cause

coeff <- max(flu_time$obs_prop) / max(flu_time$pred_prop)

ggplot(flu_time, aes(x=year_week_date)) +
  geom_line(aes(y = pred_prop), linewidth = 1.1, color = '#00BA38') +
  geom_line(aes(y = obs_prop / coeff), linewidth = 0.6, color = 'black', alpha = 0.65, linetype = 2) + 
  scale_y_continuous(
    # Features of the first axis
    name = "Predicted Flu Proportion",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name= "Observed Flu Proportion")
  ) + theme_minimal() + xlab('Date (weekly)') + 
  theme(axis.title.y.left = element_text(color = "#00BA38"))




flu_time$ratio <- flu_time$sum_pred/flu_time$sum_flu


rm(list = ls())
covid_subset <- readRDS('./tmp/covid_subset.rds')

# Replace True/False with 1 and 0
covid_subset$covid <- ifelse(covid_subset$covid == 'True', 1, 0)
covid_subset$fever <- ifelse(covid_subset$fever == 'True', 1, 0)
covid_subset$loss_smell <- ifelse(covid_subset$loss_smell == 'True', 1, 0)
covid_subset$loss_app <- ifelse(covid_subset$loss_app == 'True', 1, 0)
covid_subset$myalgia <- ifelse(covid_subset$myalgia == 'True', 1, 0)
covid_subset$hypoxemia <- ifelse(covid_subset$hypoxemia == 'True', 1, 0)
covid_subset$short_breath <- ifelse(covid_subset$short_breath == 'True', 1, 0)
covid_subset$cough <- ifelse(covid_subset$cough == 'True', 1, 0)
covid_subset$chest_pain <- ifelse(covid_subset$chest_pain == 'True', 1, 0)
covid_subset$nausea <- ifelse(covid_subset$nausea == 'True', 1, 0)
covid_subset$sore_throat <- ifelse(covid_subset$sore_throat == 'True', 1, 0)
covid_subset$fatigue <- ifelse(covid_subset$fatigue == 'True', 1, 0)

# Make state and year variables
covid_subset$state_fips <- substr(covid_subset$county_fips, 1, 2)
covid_subset$year <- year(as.Date(paste(covid_subset$year_week, 1, sep = '-'), format = '%Y-%U-%u'))
covid_subset$year <- ifelse(is.na(covid_subset$year), 2020, covid_subset$year)

# Drop asymptomatic
covid_subset <- covid_subset %>% 
  mutate(symp_sum = fever + loss_app + loss_smell + myalgia + hypoxemia + 
           short_breath + cough + chest_pain + nausea + sore_throat)

covid_subset$patient_count_imp <- as.numeric(covid_subset$patient_count_imp)
sum(covid_subset$patient_count_imp[covid_subset$covid == 1]) / sum(covid_subset$patient_count_imp)

covid_subset <- covid_subset %>% filter(symp_sum != 0)
sum(covid_subset$patient_count_imp[covid_subset$covid == 1]) / sum(covid_subset$patient_count_imp)


# Create training and testing data sets
sample <- sample(c(TRUE, FALSE), nrow(covid_subset), replace=TRUE, prob=c(0.4,0.6))
covid_train  <- covid_subset[sample, ]
covid_test   <- covid_subset[!sample, ]

# Simple logistic regression
train_model <- glm(covid ~ fever + loss_app + loss_smell + myalgia + hypoxemia + 
                     short_breath + cough + chest_pain + nausea + sore_throat + 
                     fatigue + state_fips + as.factor(year), 
                   family = binomial(link = "logit"),
                   data = covid_train,
                   weights = patient_count_imp)
#fever + myalgia + short_breath + cough + nausea + sore_throat +
#  diarrhea + fatigue + headache + congestion + sneezing
summary(train_model)

covid_test$pred <- predict(train_model, newdata = covid_test, type = "response")

basicplot <- ggplot(covid_test, aes(d = covid, m = pred)) + 
  geom_roc(n.cuts=20,labels=FALSE) + 
  style_roc(theme = theme_grey)
basicplot

roc_object <- pROC::roc(covid_test$covid, covid_test$pred)

# calculate area under curve
pROC::auc(roc_object)







library(randomForest)

data <- flu_subset_ind_1 |> ungroup() |>
  select(c(flu, fever, cough, sore_throat,
             myalgia, short_breath, hypoxemia, 
             nausea_vom, bronchitis, chest_pain,
             diarrhea, fatigue, headache,
             congestion, sneezing)) |> 
  sample_n(1000)
set.seed(222)
ind <- sample(2, nrow(data), replace = TRUE, prob = c(0.7, 0.3))
train <- data[ind==1,]
test <- data[ind==2,]

rf <- randomForest(flu~., data=train, proximity=TRUE) 


print(rf)
library(caret)

p1 <- predict(rf, train)
confusionMatrix(p1, train$flu)











