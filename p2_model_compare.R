
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

train_model_fever <- readRDS('./tmp/training_model_fever.rds') 
train_model_symptoms <- readRDS('./tmp/training_model_symptoms.rds') 
train_model_symptoms_age <- readRDS('./tmp/training_model_symptoms_age.rds') 
train_model_fever_time <- readRDS('./tmp/training_model_fever_time.rds') 
train_model_symptoms_time <- readRDS('./tmp/training_model_symptoms_time.rds') 
train_model_symptoms_time_age <- readRDS('./tmp/training_model_symptoms_time_age.rds') 

flu_p2 <- readRDS('./tmp/flu_p2_dat_proc.rds')
county_pop_group <- readRDS('./tmp/county_group_pop.rds')

flu_p2 <- flu_p2 |> filter(county_fips == '36061')

# flu_p2 <- left_join(flu_p2, county_pop_group, by = c('county_fips' = 'county_fips'))
# flu_p2_lim <- flu_p2 |> filter(county_group_pop > 50000)
# remove(flu_p2)
# flu_p2 <- flu_p2_lim
# remove(flu_p2_lim)

flu_symp_p2_16_17 <- flu_p2 |> 
  filter(week_date > as.Date('2016-08-31', "%Y-%m-%d")) |>
  filter(week_date < as.Date('2017-09-01', "%Y-%m-%d"))

#flu_p2 <- flu_symp_p2_16_17
#remove(flu_symp_p2_16_17)

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

flu_symp_p2 <- flu_p2 |> dplyr::filter(symp_count > 0)
remove(flu_p2)

flu_symp_p2$ili <- ifelse((flu_symp_p2$fever == 1 & flu_symp_p2$cough == 1), 1, 0)
flu_symp_p2$ili <- ifelse((flu_symp_p2$fever == 1 & flu_symp_p2$sore_throat == 1), 1, flu_symp_p2$ili)
flu_symp_p2$pred_fever <- predict(train_model_fever, newdata = flu_symp_p2, type = "response")
flu_symp_p2$pred_symptoms <- predict(train_model_symptoms, newdata = flu_symp_p2, type = "response")
flu_symp_p2$pred_symptoms_age <- predict(train_model_symptoms_age, newdata = flu_symp_p2, type = "response")
flu_symp_p2$pred_fever_time <- predict(train_model_fever_time, newdata = flu_symp_p2, type = "response")
flu_symp_p2$pred_symptoms_time <- predict(train_model_symptoms_time, newdata = flu_symp_p2, type = "response")
flu_symp_p2$pred_symptoms_time_age <- predict(train_model_symptoms_time_age, newdata = flu_symp_p2, type = "response")

all_cause <- read.csv('raw/county_season_ac_v2.csv')
all_cause <- all_cause |> filter(season == '2016-2017')

flu_symp_p2 <- left_join(flu_symp_p2, all_cause, by = c('county_fips' = 'county_fips'))

p2_week_county_sum <- flu_symp_p2 |> group_by(week_date, county_fips) |> 
  mutate(sum_flu = sum(flu*patient_count_imp) / all_cause,
         sum_ili = sum(ili*patient_count_imp) / all_cause,
         sum_fever = sum(pred_fever*patient_count_imp) / all_cause,
         sum_symptoms = sum(pred_symptoms*patient_count_imp) / all_cause,
         sum_fever_time = sum(pred_fever_time*patient_count_imp) / all_cause,
         sum_symptoms_time = sum(pred_symptoms_time*patient_count_imp) / all_cause,
         sum_symptoms_age = sum(pred_symptoms_age*patient_count_imp) / all_cause,
         sum_symptoms_time_age = sum(pred_symptoms_time_age*patient_count_imp) / all_cause) |>
  distinct(week_date, county_fips, sum_flu, sum_ili, sum_fever, sum_fever_time,
           sum_symptoms, sum_symptoms_time, sum_symptoms_age, 
           sum_symptoms_time_age, .keep_all = FALSE) |>
  ungroup() |>
  group_by(county_fips) |>
  mutate(roll_flu = rollmean(x = sum_flu, k = 4, fill = NA, align = "right"),
         roll_ili = rollmean(x = sum_ili, k = 4, fill = NA, align = "right"),
         roll_fever = rollmean(x = sum_fever, k = 4, fill = NA, align = "right"),
         roll_fever_time = rollmean(x = sum_fever_time, k = 4, fill = NA, align = "right"),
         roll_symptoms = rollmean(x = sum_symptoms, k = 4, fill = NA, align = "right"),
         roll_symptoms_age = rollmean(x = sum_symptoms_age, k = 4, fill = NA, align = "right"),
         roll_symptoms_time = rollmean(x = sum_symptoms_time, k = 4, fill = NA, align = "right"),
         roll_symptoms_time_age = rollmean(x = sum_symptoms_time_age, k = 4, fill = NA, align = "right")) |>
  mutate(flu_scale = ((roll_flu - mean(roll_flu, na.rm = TRUE)) / sd(roll_flu, na.rm = TRUE)),
         ili_scale = ((roll_ili - mean(roll_ili, na.rm = TRUE)) / sd(roll_ili, na.rm = TRUE)),
         fever_scale = ((roll_fever - mean(roll_fever, na.rm = TRUE)) / sd(roll_fever, na.rm = TRUE)),
         fever_time_scale = ((roll_fever_time - mean(roll_fever_time, na.rm = TRUE)) / sd(roll_fever_time, na.rm = TRUE)),
         symptoms_scale = ((roll_symptoms - mean(roll_symptoms, na.rm = TRUE)) / sd(roll_symptoms, na.rm = TRUE)),
         symptoms_time_scale = ((roll_symptoms_time - mean(roll_symptoms_time, na.rm = TRUE)) / sd(roll_symptoms_time, na.rm = TRUE)),
         symptoms_age_scale = ((roll_symptoms_age - mean(roll_symptoms_age, na.rm = TRUE)) / sd(roll_symptoms_age, na.rm = TRUE)),
         symptoms_time_age_scale = ((roll_symptoms_time_age - mean(roll_symptoms_time_age, na.rm = TRUE)) / sd(roll_symptoms_time_age, na.rm = TRUE)))

ggplot(p2_week_county_sum[p2_week_county_sum$county_fips == '12086',], aes(x=week_date, y=roll_flu)) + 
  geom_line(aes(y = flu_scale), color = 'black') + 
  geom_line(aes(y = ili_scale), color = 'black', linetype = "dashed") + 
  geom_line(aes(y = fever_time_scale), color = "#B79F00") +
  geom_line(aes(y = symptoms_scale), color = "#00BA38") + 
  geom_line(aes(y = symptoms_age_scale), color = "#00BFC4") +
  geom_line(aes(y = symptoms_time_scale), color = "#619CFF") +
  geom_line(aes(y = symptoms_time_age_scale), color = "#F564E3") +
  theme_minimal() + ylab('Number of Cases') + xlab('Time (weeks)')


p2_week_county_max <- p2_week_county_sum |> select(c(county_fips, week_date, flu_scale, ili_scale,
                               fever_time_scale, fever_scale, symptoms_scale, symptoms_age_scale,
                               symptoms_time_scale, symptoms_time_age_scale)) |>
  pivot_longer(!c(county_fips, week_date), names_to = "measure", values_to = "z_score") |>
  group_by(county_fips, measure) |>
  slice(which.max(z_score)) |> ungroup() |>
  group_by(county_fips) |> mutate(reference = week_date[3]) |>
  mutate(difference = week_date - reference)

ggplot(p2_week_county_max, aes(x=measure, y=difference)) + 
  geom_boxplot()

p2_week_county_max |> group_by(measure) |>
  mutate(rmse = sqrt(mean(as.numeric(difference)^2))) |>
  distinct(measure, rmse)


p2_week_county_start <- p2_week_county_sum |> select(c(county_fips, week_date, flu_scale, ili_scale,
                                                     fever_time_scale, fever_scale, symptoms_scale, symptoms_age_scale,
                                                     symptoms_time_scale, symptoms_time_age_scale)) |>
  pivot_longer(!c(county_fips, week_date), names_to = "measure", values_to = "z_score") |>
  arrange(county_fips, week_date) |>
  group_by(county_fips, measure) |>
  filter(z_score > 0) |>
  slice(1) |> ungroup() |>
  group_by(county_fips) |> mutate(reference = week_date[3]) |>
  mutate(difference = week_date - reference)

ggplot(p2_week_county_start, aes(x=measure, y=difference)) + 
  geom_boxplot()

p2_week_county_start |> group_by(measure) |>
  mutate(rmse = sqrt(mean(as.numeric(difference)^2))) |>
  distinct(measure, rmse)


