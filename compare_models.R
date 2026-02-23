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
library(ggcorrplot)

# Set seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

#########################
# LOAD AND PREPARE DATA #
#########################

# Load processed data
# Part 1 full data
flu_p1 <- readRDS('./tmp/flu_p1_dat_proc.rds')

# Part 2 flu data (major counties)
flu_p2 <- readRDS('./tmp/flu_p2_dat_proc_select_counties.rds')

# Add time variables to Part 1 data
flu_p1$month <- month(flu_p1$month_date)
flu_p1$year <- year(flu_p1$month_date)
flu_p1$jan <- ifelse(flu_p1$month == 1, 1, 0)
flu_p1$feb <- ifelse(flu_p1$month == 2, 1, 0)
flu_p1$mar <- ifelse(flu_p1$month == 3, 1, 0)
flu_p1$apr <- ifelse(flu_p1$month == 4, 1, 0)
flu_p1$may <- ifelse(flu_p1$month == 5, 1, 0)
flu_p1$jun <- ifelse(flu_p1$month == 6, 1, 0)
flu_p1$jul <- ifelse(flu_p1$month == 7, 1, 0)
flu_p1$aug <- ifelse(flu_p1$month == 8, 1, 0)
flu_p1$sep <- ifelse(flu_p1$month == 9, 1, 0)
flu_p1$oct <- ifelse(flu_p1$month == 10, 1, 0)
flu_p1$nov <- ifelse(flu_p1$month == 11, 1, 0)
flu_p1$dec <- ifelse(flu_p1$month == 12, 1, 0)

# Add time variables to Part 2 data
flu_p2$month <- month(flu_p2$week_date)
flu_p2$year <- year(flu_p2$week_date)
flu_p2$month_date <- as.Date(paste(flu_p2$month, '01', flu_p2$year, sep = '-'), "%m-%d-%Y")
flu_p2$jan <- ifelse(flu_p2$month == 1, 1, 0)
flu_p2$feb <- ifelse(flu_p2$month == 2, 1, 0)
flu_p2$mar <- ifelse(flu_p2$month == 3, 1, 0)
flu_p2$apr <- ifelse(flu_p2$month == 4, 1, 0)
flu_p2$may <- ifelse(flu_p2$month == 5, 1, 0)
flu_p2$jun <- ifelse(flu_p2$month == 6, 1, 0)
flu_p2$jul <- ifelse(flu_p2$month == 7, 1, 0)
flu_p2$aug <- ifelse(flu_p2$month == 8, 1, 0)
flu_p2$sep <- ifelse(flu_p2$month == 9, 1, 0)
flu_p2$oct <- ifelse(flu_p2$month == 10, 1, 0)
flu_p2$nov <- ifelse(flu_p2$month == 11, 1, 0)
flu_p2$dec <- ifelse(flu_p2$month == 12, 1, 0)

# Create age group variable
# Part 1
flu_p1$age_group <- factor(flu_p1$age_grp,
                           levels = c('0', '1', '2', '3', '4', '5'),
                           labels = c('0-4', '5-12', '13-17', '18-49', '50-64', '65+'))

# Part 2
flu_p2$age_group <- factor(flu_p2$age_grp,
                           levels = c('0', '1', '2', '3', '4', '5'),
                           labels = c('0-4', '5-12', '13-17', '18-49', '50-64', '65+'))

# Restrict data to symptomatic cases only
# Part 1
flu_symp_p1 <- flu_p1 |> dplyr::filter(symp_count > 0)
#remove(flu_p1)
# Part 2
flu_symp_p2 <- flu_p2 |> dplyr::filter(symp_count > 0)
remove(flu_p2)

# TRAIN AND TEST MODELS
# Split up counties into test and train
counties <- flu_symp_p1 |> group_by(county_fips) |>
  distinct(county_fips)
sample_counties <- sample(c(TRUE, FALSE), nrow(counties), replace=TRUE, prob=c(0.5, 0.5))
counties_train  <- counties[sample_counties, ]
counties_train$train <- 1

# Split the data by merging on train counties
flu_symp_p1 <- left_join(flu_symp_p1, counties_train, by = c('county_fips' = 'county_fips'))
flu_symp_p1_train <- flu_symp_p1 |> filter(train == 1)
flu_symp_p1_test <- flu_symp_p1 |> filter(is.na(train))

# Sample from the full Part 1 data
# Unbalanced data
unbalanced <- flu_symp_p1_train |> group_by(age_group) |> 
  uncount(patient_count_imp) |> sample_n(10000)

# Balanced Data
cases <- flu_symp_p1_train |> filter(flu == 1) |> group_by(age_group) |>
  uncount(patient_count_imp) |> sample_n(5000)
non_cases <- flu_symp_p1_train |> filter(flu == 0) |> group_by(age_group) |>
  uncount(patient_count_imp) |> sample_n(5000)
balanced <- rbind(cases, non_cases)

# Train the model by age
train_model_fever <- glm(flu ~ fever + 0,
                         family = binomial(link = "logit"),
                         data = unbalanced)

train_model_symptoms <- glm(flu ~ fever + cough + sore_throat +
                                   + myalgia + short_breath + hypoxemia +
                                   nausea_vom + bronchitis + chest_pain + 
                                   diarrhea + fatigue + headache + 
                                   congestion + sneezing + 0,
                                 family = binomial(link = "logit"),
                                 data = unbalanced)

train_model_symptoms_age <- glm(flu ~ fever*as.factor(age_grp) + 
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
                                  sneezing*as.factor(age_grp) + 0,
                                family = binomial(link = "logit"),
                                data = unbalanced)

train_model_fever_time <- glm(flu ~ fever + 0 + jan + feb + mar + 
                                apr + may + jun + jul + aug + sep + 
                                oct + nov + dec,
                              family = binomial(link = "logit"),
                              data = unbalanced)

train_model_symptoms_time <- glm(flu ~ fever + cough + sore_throat +
                                   + myalgia + short_breath + hypoxemia +
                                   nausea_vom + bronchitis + chest_pain + 
                                   diarrhea + fatigue + headache + 
                                   congestion + sneezing + 0 + jan + feb + mar + 
                                   apr + may + jun + jul + aug + sep + oct + 
                                   nov + dec,
                                 family = binomial(link = "logit"),
                                 data = unbalanced)

train_model_symptoms_time_age <- glm(flu ~ fever*as.factor(age_grp) + 
                                       cough*as.factor(age_grp) + 
                                       sore_throat*as.factor(age_grp) +
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
                                       sneezing + 0 + jan + feb + mar + 
                                   apr + may + jun + jul + aug + sep + oct + 
                                     nov + dec,
                                 family = binomial(link = "logit"),
                                 data = unbalanced)



# Predict testing data
# Expand the testing data
flu_symp_p1_test_expand <- flu_symp_p1_test |> group_by(age_group) |>
  uncount(patient_count_imp) |> sample_n(10000)

# Predict
flu_symp_p1_test$pred_fever <- predict(train_model_fever, 
                                     newdata = flu_symp_p1_test, 
                                     type = "response")
flu_symp_p1_test$pred_symptoms <- predict(train_model_symptoms, 
                                      newdata = flu_symp_p1_test, 
                                      type = "response")
flu_symp_p1_test$pred_symptoms_age <- predict(train_model_symptoms_age, 
                                          newdata = flu_symp_p1_test, 
                                          type = "response")
flu_symp_p1_test$pred_fever_time <- predict(train_model_fever_time, 
                                              newdata = flu_symp_p1_test, 
                                              type = "response")
flu_symp_p1_test$pred_symptoms_time <- predict(train_model_symptoms_time, 
                                            newdata = flu_symp_p1_test, 
                                            type = "response")
flu_symp_p1_test$pred_symptoms_time_age <- predict(train_model_symptoms_time_age, 
                                               newdata = flu_symp_p1_test, 
                                               type = "response")

pROC::auc(pROC::roc(flu_symp_p1_test$flu, flu_symp_p1_test$pred_fever))
pROC::auc(pROC::roc(flu_symp_p1_test$flu, flu_symp_p1_test$pred_symptoms))
pROC::auc(pROC::roc(flu_symp_p1_test$flu, flu_symp_p1_test$pred_symptoms_age))
pROC::auc(pROC::roc(flu_symp_p1_test$flu, flu_symp_p1_test$pred_fever_time))
pROC::auc(pROC::roc(flu_symp_p1_test$flu, flu_symp_p1_test$pred_symptoms_time))
pROC::auc(pROC::roc(flu_symp_p1_test$flu, flu_symp_p1_test$pred_symptoms_time_age))

saveRDS(train_model_fever, './tmp/training_model_fever.rds') 
saveRDS(train_model_symptoms, './tmp/training_model_symptoms.rds') 
saveRDS(train_model_symptoms_age, './tmp/training_model_symptoms_age.rds') 
saveRDS(train_model_fever_time, './tmp/training_model_fever_time.rds') 
saveRDS(train_model_symptoms_time, './tmp/training_model_symptoms_time.rds') 
saveRDS(train_model_symptoms_time_age, './tmp/training_model_symptoms_time_age.rds') 





# Reproduce on Part 2 data
flu_symp_p2$pred_fever <- predict(train_model_fever, newdata = flu_symp_p2, type = "response")
#flu_p2_symp$se <- predict(all, newdata = flu_p2_symp, type = "response", se.fit = TRUE)$se.fit
#flu_p2_symp$upr <- flu_p2_symp$pred + (1.96 * flu_p2_symp$se)
#flu_p2_symp$lwr <- flu_p2_symp$pred - (1.96 * flu_p2_symp$se)

flu_symp_p2$pred_symptoms <- predict(train_model_symptoms, newdata = flu_symp_p2, type = "response")
#flu_p2_symp$se_f <- predict(fever, newdata = flu_p2_symp, type = "response", se.fit = TRUE)$se.fit
#flu_p2_symp$upr_f <- flu_p2_symp$pred_f + (1.96 * flu_p2_symp$se_f)
#flu_p2_symp$lwr_f <- flu_p2_symp$pred_f - (1.96 * flu_p2_symp$se_f)

flu_symp_p2$pred_symptoms_age <- predict(train_model_symptoms_age, newdata = flu_symp_p2, type = "response")
flu_symp_p2$pred_fever_time <- predict(train_model_fever_time, newdata = flu_symp_p2, type = "response")
flu_symp_p2$pred_symptoms_time <- predict(train_model_symptoms_time, newdata = flu_symp_p2, type = "response")
flu_symp_p2$pred_symptoms_time_age <- predict(train_model_symptoms_time_age, newdata = flu_symp_p2, type = "response")



p2_week_county_sum <- flu_symp_p2 |> group_by(week_date, county_fips) |> 
  mutate(sum_flu = sum(flu),
         sum_fever = sum(pred_fever),
         sum_symptoms = sum(pred_symptoms),
         sum_fever_time = sum(pred_fever_time),
         sum_symptoms_time = sum(pred_symptoms_time),
         sum_symptoms_age = sum(pred_symptoms_age),
         sum_symptoms_time_age = sum(pred_symptoms_time_age)) |>
  distinct(week_date, county_fips, sum_flu, sum_fever, sum_fever_time,
           sum_symptoms, sum_symptoms_time, sum_symptoms_age, 
           sum_symptoms_time_age, .keep_all = FALSE) |>
  ungroup() |>
  group_by(county_fips) |>
  mutate(roll_flu = rollmean(x = sum_flu, k = 3, fill = NA, align = "right"),
         roll_fever = rollmean(x = sum_fever, k = 3, fill = NA, align = "right"),
         roll_fever_time = rollmean(x = sum_fever_time, k = 3, fill = NA, align = "right"),
         roll_symptoms = rollmean(x = sum_symptoms, k = 3, fill = NA, align = "right"),
         roll_symptoms_age = rollmean(x = sum_symptoms_age, k = 3, fill = NA, align = "right"),
         roll_symptoms_time = rollmean(x = sum_symptoms_time, k = 3, fill = NA, align = "right"),
         roll_symptoms_time_age = rollmean(x = sum_symptoms_time_age, k = 3, fill = NA, align = "right")) 


library(scales)

hue_pal()(6)

ggplot(p2_week_county_sum, aes(x=week_date, y=roll_flu)) + 
  geom_line(aes(y = roll_flu), color = 'black') + 
  geom_line(aes(y = roll_fever), color = "#F8766D") +
  geom_line(aes(y = roll_fever_time), color = "#B79F00") +
  geom_line(aes(y = roll_symptoms), color = "#00BA38") + 
  geom_line(aes(y = roll_symptoms_age), color = "#00BFC4") +
  geom_line(aes(y = roll_symptoms_time), color = "#619CFF") +
  geom_line(aes(y = roll_symptoms_time_age), color = "#F564E3") +
  facet_wrap(vars(county_fips), scales = "free") +
  theme_minimal() + ylab('Number of Cases') + xlab('Time (weeks)')


