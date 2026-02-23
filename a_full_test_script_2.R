# Start with a clear environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(lubridate)
library(car)
library(plotROC)
library(pROC)
library(cutpointr)
library(zoo)

# Set seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

# load model
train_model <- readRDS('./tmp/train_model_bal.rds')

# Part 2 flu data
flu_p2 <- readRDS('./tmp/flu_p2_dat_proc.rds')

# Add time variables to Part 2 data
flu_p2$month <- month(flu_p2$week_date)
flu_p2$year <- year(flu_p2$week_date)

# Part 2
flu_symp_p2 <- flu_p2 |> dplyr::filter(symp_count > 0)
remove(flu_p2)

# Split data in 5
flu_symp_p2_2016 <- flu_symp_p2 |> filter(year == 2016)
flu_symp_p2_2017 <- flu_symp_p2 |> filter(year == 2017)
flu_symp_p2_2018 <- flu_symp_p2 |> filter(year == 2018)
flu_symp_p2_2019 <- flu_symp_p2 |> filter(year == 2019)
flu_symp_p2_2020 <- flu_symp_p2 |> filter(year == 2020)
remove(flu_symp_p2)

# Predict the data
flu_symp_p2_2016$pred <- predict(train_model, 
                            newdata = flu_symp_p2_2016, 
                            type = "response")
flu_symp_p2_2017$pred <- predict(train_model, 
                              newdata = flu_symp_p2_2017, 
                              type = "response")
flu_symp_p2_2018$pred <- predict(train_model, 
                            newdata = flu_symp_p2_2018, 
                            type = "response")
flu_symp_p2_2019$pred <- predict(train_model, 
                              newdata = flu_symp_p2_2019, 
                              type = "response")
flu_symp_p2_2020$pred <- predict(train_model, 
                                 newdata = flu_symp_p2_2020, 
                                 type = "response")

# Assign flu positive/negative
flu_symp_p2_2016$pred_flu <- ifelse(flu_symp_p2_2016$pred > 0.7034, 1, 0)
flu_symp_p2_2017$pred_flu <- ifelse(flu_symp_p2_2017$pred > 0.7034, 1, 0)
flu_symp_p2_2018$pred_flu <- ifelse(flu_symp_p2_2018$pred > 0.7034, 1, 0)
flu_symp_p2_2019$pred_flu <- ifelse(flu_symp_p2_2019$pred > 0.7034, 1, 0)
flu_symp_p2_2020$pred_flu <- ifelse(flu_symp_p2_2020$pred > 0.7034, 1, 0)

# Append the data
flu_symp_p2 <- rbind(flu_symp_p2_2016, flu_symp_p2_2017, flu_symp_p2_2018,
                     flu_symp_p2_2019, flu_symp_p2_2020)
remove(flu_symp_p2_2016, flu_symp_p2_2017, flu_symp_p2_2018,
       flu_symp_p2_2019, flu_symp_p2_2020)

# Collapse to age, week, county level
flu_age_county_week <- flu_symp_p2 |> 
  group_by(age_grp, county_fips, week_date) |>
  mutate(flu_cases = sum(flu*patient_count_imp),
         pred_flu_cases = sum(pred_flu*patient_count_imp)) |>
  distinct(county_fips, week_date, age_grp, flu_cases, pred_flu_cases)

saveRDS(flu_age_county_week, './tmp/flu_age_county_week.rds') 

flu_county_week <- flu_age_county_week |> 
  group_by(county_fips, week_date) |>
  mutate(flu_cases = sum(flu_cases),
         pred_flu_cases = sum(pred_flu_cases)) |>
  distinct(county_fips, week_date, flu_cases, pred_flu_cases)

saveRDS(flu_county_week, './tmp/flu_county_week.rds') 


