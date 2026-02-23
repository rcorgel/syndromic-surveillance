####################
# 1. SET-UP SCRIPT #
####################

# Start with a clear environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(lubridate)
library(ISOweek)
library(sf)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

###############
# 2. FLU DATA #
###############

flu <- read_csv('raw/flu_2016-09_2017-08_p2.csv')

# Impute <=5 to a random number between 1 and 5
flu$patient_count_imp <- flu$patient_count
flu$patient_count_imp[flu$patient_count_imp == '<=5'] <- 
  sample(1:5, length(flu$patient_count_imp[flu$patient_count_imp == '<=5']), replace= TRUE)
flu$patient_count_imp <- as.numeric(flu$patient_count_imp)

# Count flu symptoms
flu_16_17_group <- flu |> 
  mutate(symp_count = fever + myalgia + 
           cough + sore_throat + 
           short_breath + hypoxemia + 
           chest_pain + bronchitis + 
           nausea_vom + diarrhea + fatigue + 
           headache + congestion + sneezing) |>
  dplyr::filter(symp_count > 0) |>
  group_by(year_week) |>
  mutate(symptom_claims = sum(patient_count_imp),
         flu_claims = sum(flu*patient_count_imp)) |>
  distinct(year_week, symptom_claims, flu_claims)

remove(flu)
  
flu <- read_csv('raw/flu_2017-09_2018-08_p2.csv')

# Impute <=5 to a random number between 1 and 5
flu$patient_count_imp <- flu$patient_count
flu$patient_count_imp[flu$patient_count_imp == '<=5'] <- 
  sample(1:5, length(flu$patient_count_imp[flu$patient_count_imp == '<=5']), replace= TRUE)
flu$patient_count_imp <- as.numeric(flu$patient_count_imp)

# Count flu symptoms
flu_17_18_group <- flu |> 
  mutate(symp_count = fever + myalgia + 
           cough + sore_throat + 
           short_breath + hypoxemia + 
           chest_pain + bronchitis + 
           nausea_vom + diarrhea + fatigue + 
           headache + congestion + sneezing) |>
  dplyr::filter(symp_count > 0) |>
  group_by(year_week) |>
  mutate(symptom_claims = sum(patient_count_imp),
         flu_claims = sum(flu*patient_count_imp)) |>
  distinct(year_week, symptom_claims, flu_claims)

remove(flu)

flu <- read_csv('raw/flu_2018-09_2019-08_p2.csv')

# Impute <=5 to a random number between 1 and 5
flu$patient_count_imp <- flu$patient_count
flu$patient_count_imp[flu$patient_count_imp == '<=5'] <- 
  sample(1:5, length(flu$patient_count_imp[flu$patient_count_imp == '<=5']), replace= TRUE)
flu$patient_count_imp <- as.numeric(flu$patient_count_imp)

# Count flu symptoms
flu_18_19_group <- flu |> 
  mutate(symp_count = fever + myalgia + 
           cough + sore_throat + 
           short_breath + hypoxemia + 
           chest_pain + bronchitis + 
           nausea_vom + diarrhea + fatigue + 
           headache + congestion + sneezing) |>
  dplyr::filter(symp_count > 0) |>
  group_by(year_week) |>
  mutate(symptom_claims = sum(patient_count_imp),
         flu_claims = sum(flu*patient_count_imp)) |>
  distinct(year_week, symptom_claims, flu_claims)

remove(flu)

flu <- read_csv('raw/flu_2019-09_2020-08_p2.csv')

# Impute <=5 to a random number between 1 and 5
flu$patient_count_imp <- flu$patient_count
flu$patient_count_imp[flu$patient_count_imp == '<=5'] <- 
  sample(1:5, length(flu$patient_count_imp[flu$patient_count_imp == '<=5']), replace= TRUE)
flu$patient_count_imp <- as.numeric(flu$patient_count_imp)

# Count flu symptoms
flu_19_20_group <- flu |> 
  mutate(symp_count = fever + myalgia + 
           cough + sore_throat + 
           short_breath + hypoxemia + 
           chest_pain + bronchitis + 
           nausea_vom + diarrhea + fatigue + 
           headache + congestion + sneezing) |>
  dplyr::filter(symp_count > 0) |>
  group_by(year_week) |>
  mutate(symptom_claims = sum(patient_count_imp),
         flu_claims = sum(flu*patient_count_imp)) |>
  distinct(year_week, symptom_claims, flu_claims)

remove(flu)

flu <- read_csv('raw/flu_2020-09_2021-08_p2.csv')

# Impute <=5 to a random number between 1 and 5
flu$patient_count_imp <- flu$patient_count
flu$patient_count_imp[flu$patient_count_imp == '<=5'] <- 
  sample(1:5, length(flu$patient_count_imp[flu$patient_count_imp == '<=5']), replace= TRUE)
flu$patient_count_imp <- as.numeric(flu$patient_count_imp)

# Count flu symptoms
flu_20_21_group <- flu |> 
  mutate(symp_count = fever + myalgia + 
           cough + sore_throat + 
           short_breath + hypoxemia + 
           chest_pain + bronchitis + 
           nausea_vom + diarrhea + fatigue + 
           headache + congestion + sneezing) |>
  dplyr::filter(symp_count > 0) |>
  group_by(year_week) |>
  mutate(symptom_claims = sum(patient_count_imp),
         flu_claims = sum(flu*patient_count_imp)) |>
  distinct(year_week, symptom_claims, flu_claims)

remove(flu)

flu <- read_csv('raw/flu_2021-09_2022-08_p2.csv')

# Impute <=5 to a random number between 1 and 5
flu$patient_count_imp <- flu$patient_count
flu$patient_count_imp[flu$patient_count_imp == '<=5'] <- 
  sample(1:5, length(flu$patient_count_imp[flu$patient_count_imp == '<=5']), replace= TRUE)
flu$patient_count_imp <- as.numeric(flu$patient_count_imp)

# Count flu symptoms
flu_21_22_group <- flu |> 
  mutate(symp_count = fever + myalgia + 
           cough + sore_throat + 
           short_breath + hypoxemia + 
           chest_pain + bronchitis + 
           nausea_vom + diarrhea + fatigue + 
           headache + congestion + sneezing) |>
  dplyr::filter(symp_count > 0) |>
  group_by(year_week) |>
  mutate(symptom_claims = sum(patient_count_imp),
         flu_claims = sum(flu*patient_count_imp)) |>
  distinct(year_week, symptom_claims, flu_claims)

remove(flu)

flu <- read_csv('raw/flu_2022-09_2023-08_p2.csv')

# Impute <=5 to a random number between 1 and 5
flu$patient_count_imp <- flu$patient_count
flu$patient_count_imp[flu$patient_count_imp == '<=5'] <- 
  sample(1:5, length(flu$patient_count_imp[flu$patient_count_imp == '<=5']), replace= TRUE)
flu$patient_count_imp <- as.numeric(flu$patient_count_imp)

# Count flu symptoms
flu_22_23_group <- flu |> 
  mutate(symp_count = fever + myalgia + 
           cough + sore_throat + 
           short_breath + hypoxemia + 
           chest_pain + bronchitis + 
           nausea_vom + diarrhea + fatigue + 
           headache + congestion + sneezing) |>
  dplyr::filter(symp_count > 0) |>
  group_by(year_week) |>
  mutate(symptom_claims = sum(patient_count_imp),
         flu_claims = sum(flu*patient_count_imp)) |>
  distinct(year_week, symptom_claims, flu_claims)

remove(flu)

# Combine data
flu <- rbind(flu_16_17_group, flu_17_18_group, flu_18_19_group, 
             flu_19_20_group, flu_20_21_group, flu_21_22_group, 
             flu_22_23_group)

# Convert to week date
flu$week_date <- sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", flu$year_week)
flu$week_date <- ISOweek2date(flu$week_date)

# Re-collapse data to the week level
symptoms_week <- flu |> group_by(week_date) |>
  mutate(symptoms_sum = sum(symptom_claims),
         flu_sum = sum(flu_claims)) |>
  distinct(week_date, symptoms_sum, flu_sum)

# Load test data
tests <- read_csv('raw/county_week_flu_tests_v2.csv')

# Impute test data
tests$flu_tests_imp <- tests$flu_tests
tests$flu_tests_imp[tests$flu_tests_imp == '<=5'] <- 
  sample(1:5, length(tests$flu_tests_imp[tests$flu_tests_imp == '<=5']), 
         replace= TRUE)
tests$flu_tests_imp <- as.numeric(tests$flu_tests_imp)

# Convert to week date
tests$week_date <- sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", tests$year_week)
tests$week_date <- ISOweek2date(tests$week_date)

# Collapse data to week level
tests_week <- tests |> group_by(week_date) |>
  mutate(tests_sum = sum(flu_tests_imp)) |>
  distinct(week_date, tests_sum)

# Plot tests over time
ggplot(tests_week, aes(x = week_date, y = tests_sum)) + 
  geom_line() + theme_minimal() + ylab('Number of Flu Tests') + 
  xlab('Week') + ggtitle('Influenza Tests, 2016 - 2024')

# Merge on symptom counts
tests_week <- left_join(tests_week, symptoms_week, 
                        by = c('week_date' = 'week_date'))

# Plot symptoms over time
ggplot(symptoms_week, aes(x = week_date, y = symptoms_sum)) + 
  geom_line() + theme_minimal() + ylab('Number of Symptom Visits') + 
  xlab('Week') + ggtitle('Symptom Visits, 2016 - 2023')

# Plot normalized tests
tests_week_norm <- tests_week
tests_week_norm <- tests_week_norm |> dplyr::filter(week_date < '2023-08-07')

tests_week_norm$tests_norm <- tests_week_norm$tests_sum / tests_week_norm$symptoms_sum
tests_week_norm$flu_norm <- tests_week_norm$flu_sum / tests_week_norm$tests_sum

ggplot(tests_week_norm, aes(x = week_date, y = flu_sum)) + 
  geom_line() + theme_minimal() + ylab('Number of Flu Tests / Symptom Visits') + 
  xlab('Week') + ggtitle('Influenza Tests / Symptom Visits, 2016 - 2023')

ggplot(tests_week_norm, aes(x = week_date, y = tests_norm)) + 
  geom_line() + theme_minimal() + ylab('Number of Flu Tests / Symptom Visits') + 
  xlab('Week') + ggtitle('Influenza Tests / Symptom Visits, 2016 - 2023')


ggplot(tests_week_norm, aes(x = week_date, y = flu_norm)) + 
  geom_line() + theme_minimal() + ylab('Number of Flu / Tests') + 
  xlab('Week') + ggtitle('Flu / Tests, 2016 - 2023')

tests_week_norm$month <- month(tests_week_norm$week_date)
tests_week_norm$year <- year(tests_week_norm$week_date)

ggplot(tests_week_norm, aes(x = flu_sum, y = tests_sum, color = as.factor(month))) + 
  geom_point() + theme_minimal() + ylab('tests') + 
  xlab('flu') + ggtitle('Flu & Tests Scatter, 2016 - 2023')

tests_week_norm$season <- ifelse(tests_week_norm$year == 2016 & tests_week_norm$month > 8, '2016-2017', '')
tests_week_norm$season <- ifelse(tests_week_norm$year == 2017 & tests_week_norm$month < 9, '2016-2017', tests_week_norm$season)
tests_week_norm$season <- ifelse(tests_week_norm$year == 2017 & tests_week_norm$month > 8, '2017-2018', tests_week_norm$season)
tests_week_norm$season <- ifelse(tests_week_norm$year == 2018 & tests_week_norm$month < 9, '2017-2018', tests_week_norm$season)
tests_week_norm$season <- ifelse(tests_week_norm$year == 2018 & tests_week_norm$month > 8, '2018-2019', tests_week_norm$season)
tests_week_norm$season <- ifelse(tests_week_norm$year == 2019 & tests_week_norm$month < 9, '2018-2019', tests_week_norm$season)
tests_week_norm$season <- ifelse(tests_week_norm$year == 2019 & tests_week_norm$month > 8, '2019-2020', tests_week_norm$season)
tests_week_norm$season <- ifelse(tests_week_norm$year == 2020 & tests_week_norm$month < 9, '2019-2020', tests_week_norm$season)
tests_week_norm$season <- ifelse(tests_week_norm$year == 2020 & tests_week_norm$month > 8, '2020-2021', tests_week_norm$season)
tests_week_norm$season <- ifelse(tests_week_norm$year == 2021 & tests_week_norm$month < 9, '2020-2021', tests_week_norm$season)
tests_week_norm$season <- ifelse(tests_week_norm$year == 2021 & tests_week_norm$month > 8, '2021-2022', tests_week_norm$season)
tests_week_norm$season <- ifelse(tests_week_norm$year == 2022 & tests_week_norm$month < 9, '2021-2022', tests_week_norm$season)
tests_week_norm$season <- ifelse(tests_week_norm$year == 2022 & tests_week_norm$month > 8, '2022-2023', tests_week_norm$season)
tests_week_norm$season <- ifelse(tests_week_norm$year == 2023 & tests_week_norm$month < 9, '2022-2023', tests_week_norm$season)

ggplot(tests_week_norm, aes(x = flu_sum, y = tests_sum, color = as.factor(season))) + 
  geom_point() + theme_minimal() + ylab('tests') + 
  xlab('flu') + ggtitle('Flu & Tests Scatter, 2016 - 2023')

# Value used to transform the data
coeff <- 10

ggplot(tests_week_norm, aes(x=week_date)) +
  
  geom_line( aes(y=tests_sum)) + 
  geom_line( aes(y=flu_sum / 0.5), color = 'red') + # Divide by 10 to get the same range than the temperature
  
  scale_y_continuous(
    
    # Features of the first axis
    name = "Tests",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*0.5, name="Flu", )
  ) + 
  theme(
    axis.title.y.right = element_text(color = 'red', size=13)
  )

################################################################################
################################################################################
