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

# Load processed data
# Part 1 full data
flu_p1 <- readRDS('./tmp/flu_p1_dat_proc.rds')

# Add time variables to Part 1 data
flu_p1$month <- month(flu_p1$month_date)
flu_p1$year <- year(flu_p1$month_date)

# Add time variables to Part 2 data
flu_p2$month <- month(flu_p2$week_date)
flu_p2$year <- year(flu_p2$week_date)

# Create age group variable
# Part 1
flu_p1$age_group <- factor(flu_p1$age_grp,
                           levels = c('0', '1', '2', '3', '4', '5'),
                           labels = c('0-4', '5-12', '13-17', '18-49', 
                                      '50-64', '65+'))

# Part 2
flu_p2$age_group <- factor(flu_p2$age_grp,
                           levels = c('0', '1', '2', '3', '4', '5'),
                           labels = c('0-4', '5-12', '13-17', '18-49', 
                                      '50-64', '65+'))

# Restrict data to symptomatic cases only
# Part 1
flu_symp_p1 <- flu_p1 |> dplyr::filter(symp_count > 0)

# Part 2
flu_symp_p2 <- flu_p2 |> dplyr::filter(symp_count > 0)

# Split part 1 data into training and testing
counties <- unique(flu_symp_p1$county_fips)
# Take a random sample of 50% of counties
train <- as.data.frame(sample(counties, size = length(counties)*0.50, replace = F))
colnames(train) <- 'county_fips'
train$train <- 1
flu_symp_p1 <- left_join(flu_symp_p1, train, by = c('county_fips' = 'county_fips'))
flu_symp_p1_train <- flu_symp_p1 |> filter(train == 1)
flu_symp_p1_test <- flu_symp_p1 |> filter(is.na(train))
remove(flu_symp_p1)

# Balance the training data
flu_symp_p1_train_bal <- flu_symp_p1_train |> group_by(flu) |>
  uncount(patient_count_imp) |> sample_n(250000)

# Run a simple model
train_model_bal <- glm(flu ~ fever + cough + sore_throat +
                         myalgia + short_breath + hypoxemia +
                         nausea_vom + bronchitis + chest_pain + 
                         diarrhea + fatigue + headache + 
                         congestion + sneezing + as.factor(month),
                       family = binomial(link = "logit"),
                       data = flu_symp_p1_train_bal)

summary(train_model_bal)

# Run prediction on test set
flu_symp_p1_test$pred <- predict(train_model_bal, 
                                 newdata = flu_symp_p1_test, 
                                 type = "response")
pROC::auc(pROC::roc(flu_symp_p1_test$flu, flu_symp_p1_test$pred))

# Examine cutpoint metrics
cp_bal <- cutpointr(flu_symp_p1_test, pred, flu, 
                    method = maximize_metric, metric = F1_score)
summary(cp_bal)
plot(cp_bal)
plot_metric(cp_bal)

# Save model
saveRDS(train_model_bal, './tmp/train_model_bal.rds') 
