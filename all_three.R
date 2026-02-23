################################################################################
# File Name: process_data_p1                                                   #
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
library(sf)
library(spdep)
library(ISOweek)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

#############################
# 2. PROCESS INFLUENZA DATA #
#############################

flu <- read.csv('raw/flu_2022-09_2023-08_p3.csv')

# Impute data
# Check the percent of rows that will be imputed (81.2%)
nrow(flu[flu$patient_count == "<=5",]) / nrow(flu)
# Impute rows where patient count is <=5 to a random number 1 to 5
flu$patient_count_imp <- flu$patient_count
flu$patient_count_imp[flu$patient_count_imp == '<=5'] <- 
  sample(1:5, length(flu$patient_count_imp[flu$patient_count_imp == '<=5']), 
         replace= TRUE)
# Convert imputed count to numeric
flu$patient_count_imp <- as.numeric(flu$patient_count_imp)
# Check the percent of patients that were imputed (8.7%)
sum(flu[flu$patient_count_imp <= 5,]$patient_count_imp) / sum(flu$patient_count_imp)

# Create new variables
# Create new variables
flu <- flu |>
  # Convert week date to date
  mutate(week_date = sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", year_week)) |>
  mutate(week_date = ISOweek2date(week_date)) |>
  # Create state fips
  mutate(state_fips = substr(county_fips, 1, 2)) |>
  # Create symptom count
  mutate(symp_count = fever + myalgia + cough + 
           sore_throat + short_breath + hypoxemia + 
           chest_pain + bronchitis + nausea_vom + 
           diarrhea + fatigue + headache + congestion + sneezing)

# Check the percent of patients that will be dropped for:
# Missing county or not US state & DC (0.4%)
sum(flu[flu$state_fips %in% c("72","99",""),]$patient_count_imp) / sum(flu$patient_count_imp)
# Asymptomatic or non-coded individuals (83.4%)
sum(flu[flu$symp_count == 0,]$patient_count_imp) / sum(flu$patient_count_imp)

# Drop based on:
flu_filt <- flu |>
  # Missing county or not US state & DC
  filter(!(state_fips %in% c("72","99","")))
# Check the number of individuals kept (98.2%)
sum(flu_filt$patient_count_imp) / sum(flu$patient_count_imp)

# Confirm drop was completed
flu_filt |>
  # There should be more than 2500 county groups
  verify(state_fips != "99" | state_fips != "72") |>
  # Check disease, symptom, and count vars are not missing
  assert(not_na, county_fips:symp_count)

flu_filt_17 <- flu_filt |> 
  group_by(week_date, flu, covid, rsv) |>
  mutate(patient_sum = sum(patient_count_imp)) |>
  distinct(week_date, flu, covid, rsv, patient_sum, .keep_all = FALSE) |>
  filter((flu != 0 | covid != 0 | rsv != 0)) |>
  mutate(disease = ifelse((flu == 1 & covid == 0 & rsv == 0), 'flu', ''),
         disease = ifelse((flu == 0 & covid == 0 & rsv == 1), 'rsv', disease),
         disease = ifelse((flu == 0 & covid == 1 & rsv == 0), 'cov', disease),
         disease = ifelse((flu == 1 & covid == 1 & rsv == 0), 'flu + cov', disease),
         disease = ifelse((flu == 1 & covid == 0 & rsv == 1), 'flu + rsv', disease),
         disease = ifelse((flu == 0 & covid == 1 & rsv == 1), 'cov + rsv', disease),
         disease = ifelse((flu == 1 & covid == 1 & rsv == 1), 'all', disease)) |>
  distinct(week_date, disease, patient_sum, .keep_all = FALSE) |>
  ungroup() |> group_by(week_date) |>
  mutate(percent = patient_sum / sum(patient_sum))
  
ggplot() + 
  geom_line(data = flu_filt_17, 
            aes(x = week_date, y = patient_sum, color = disease), size = 0.9) + 
  theme_minimal() + ylab('Count') + xlab('Week')

ggplot() + 
  geom_line(data = flu_filt_17, 
            aes(x = week_date, y = patient_sum, color = disease), size = 0.9) + 
  theme_minimal() + ylab('Count') + xlab('Week') + ylim(0, 8000)

ggplot() + 
  geom_line(data = flu_filt_17, 
            aes(x = week_date, y = percent, color = disease), size = 0.9) + 
  theme_minimal() + ylab('Percent') + xlab('Week')

ggplot() + 
  geom_line(data = flu_filt_17, 
            aes(x = week_date, y = percent, color = disease), size = 0.9) + 
  theme_minimal() + ylab('Percent') + xlab('Week') + ylim(0, 0.02)





# Create lagged fever, sore throat, and cough percent variables
lagged_vars <- flu_filt |> group_by(county_fips, week_date) |>
  mutate(fever_perc_lag = sum(fever*patient_count_imp) / sum(patient_count_imp),
         cough_perc_lag = sum(cough*patient_count_imp) / sum(patient_count_imp),
         sore_throat_perc_lag = sum(sore_throat*patient_count_imp) / sum(patient_count_imp)) |>
  distinct(county_fips, week_date, fever_perc_lag, cough_perc_lag, sore_throat_perc_lag, .keep_all = FALSE) |>
  mutate(lag_week_date = week_date + 7) |> ungroup() |>
  select(-c(week_date))

# Merge on lagged variables
flu_filt <- left_join(flu_filt, lagged_vars, by = c('county_fips' = 'county_fips',
                                                    'week_date' = 'lag_week_date'))

# Drop symptomatic individuals
flu_filt_symp <- flu_filt |>
  filter(symp_count > 0)
# Check the number of individuals kept (16.3%)
sum(flu_filt_symp$patient_count_imp) / sum(flu$patient_count_imp)

flu_filt_3 <- flu_filt_symp |> 
  mutate(flu_only = ifelse((flu == 1 & covid == 0 & rsv == 0), 1, 0),
         other_resp = ifelse(((flu == 0 & covid == 1) | (flu == 0 & rsv == 1)), 1, 0),
         no_dis = ifelse((flu == 0 & covid == 0 & rsv == 0), 1, 0)) #|>
  #filter(flu_only == 1 | no_dis != 1)




cases_16_17 <- flu_filt_3[flu_filt_3$flu == 1, ]
non_cases_16_17 <- flu_filt_3[flu_filt_3$flu != 1, ]


cases_exp_sub_16_17 <- cases_16_17 |> uncount(patient_count_imp) |> 
  group_by(year_week) |> sample_n(800)

non_cases_exp_sub_16_17 <- non_cases_16_17 |> uncount(patient_count_imp) |> 
  group_by(year_week) |> filter(week_date <= as.Date('2023-08-07')) |>
  sample_n(800)

flu_test <- rbind(cases_exp_sub_16_17, non_cases_exp_sub_16_17)

sample <- sample(c(TRUE, FALSE), nrow(flu_test), replace=TRUE, prob=c(0.5, 0.5))
flu_train  <- flu_test[sample, ]
flu_test   <- flu_test[!sample, ]



model <- glm(flu ~ fever + cough + sore_throat +
               myalgia + short_breath + hypoxemia + 
               nausea_vom + bronchitis + chest_pain +
               diarrhea + fatigue + headache +
               congestion + sneezing,
    family = binomial(link = "logit"),
    data = flu_train)

summary(model)


flu_test$pred <- predict(model, newdata = flu_test, type = "response")

basicplot <- ggplot(flu_test, aes(d = flu, m = pred)) + 
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

pROC::auc(pROC::roc(flu_test$flu, flu_test$pred))

opt_cut_1617_youden <- cutpointr(flu_test, pred, flu, direction = ">=", 
                                 pos_class = 1, neg_class = 0, 
                                 method = maximize_metric, metric = youden, na.rm = TRUE)
summary(opt_cut_1617_youden)


saveRDS(model, './tmp/training_models_all_three_the_rest.rds') 








flu <- read.csv('raw/flu_2016-09_2017-08_p3.csv')

collapse_1 <- 

