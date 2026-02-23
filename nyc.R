
# Load libraries
library(tidyverse)
library(lubridate)
library(car)
library(plotROC)
library(pROC)
library(cutpointr)
library(zoo)
library('olsrr')

# Set seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

# Load processed data
# Part 1 full data
flu_p1 <- readRDS('./tmp/flu_p1_dat_proc.rds')

# Part 2 flu data (major counties)
flu_p2 <- readRDS('./tmp/flu_p2_dat_proc_select_counties.rds')
#flu_p2 <- flu_p2 |> filter(county_fips == '36061' | county_fips == '17031')

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
flu_p2$epi_week <- epiweek(flu_p2$week_date)
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

# # Create lagged fever, sore throat, and cough percent variables for Part 2
# lagged_vars <- flu_p2 |> group_by(county_fips, week_date) |>
#   mutate(fever_perc_lag = sum(fever*patient_count_imp) / sum(patient_count_imp),
#          cough_perc_lag = sum(cough*patient_count_imp) / sum(patient_count_imp),
#          sore_throat_perc_lag = sum(sore_throat*patient_count_imp) / sum(patient_count_imp)) |>
#   distinct(county_fips, week_date, fever_perc_lag, cough_perc_lag, sore_throat_perc_lag, .keep_all = FALSE) |>
#   mutate(lag_week_date = week_date + 14) |> ungroup() |>
#   select(-c(week_date))
# 
# # Merge on lagged variables for Part 2
# flu_p2 <- left_join(flu_p2, lagged_vars, by = c('county_fips' = 'county_fips',
#                                                 'week_date' = 'lag_week_date'))

# Restrict data to symptomatic cases only
# Part 1
flu_symp_p1 <- flu_p1 |> dplyr::filter(symp_count > 0)
#remove(flu_p1)
# Part 2
flu_symp_p2 <- flu_p2 |> dplyr::filter(symp_count > 0)

# Load population data
county_pop <- readRDS('tmp/county_group_pop.rds')

# Restrict to counties with 10,000 or more individuals
flu_symp_p1 <- left_join(flu_symp_p1, county_pop, 
                         by = c('county_fips' = 'county_fips'))
flu_symp_p1 <- flu_symp_p1 |> filter(county_group_pop > 10000)

# Make data individual level (sample 1 million individuals)
# Balance flu casea and non-flu cases
balanced <- flu_symp_p1 |> group_by(flu) |>
  uncount(patient_count_imp) |> sample_n(500000)
# Observed data
unbalanced <- flu_symp_p1 |>
  uncount(patient_count_imp) |> sample_n(1000000)

# Split the data randomly into testing and training
sample <- sample(c(TRUE, FALSE), nrow(balanced), replace=TRUE, prob=c(0.5, 0.5))
flu_train  <- balanced[sample, ]
flu_test   <- balanced[!sample, ]

## Define candidate models ##

# Model with only fever
fev_model <- glm(flu ~ fever,
                 family = binomial(link = "logit"),
                 data = flu_train)
summary(fev_model)

# Model with all symptoms
all_model <- glm(flu ~ fever + cough + sore_throat +
                   myalgia + short_breath + hypoxemia + nausea_vom + 
                   bronchitis + chest_pain + diarrhea + 
                   fatigue + headache + congestion + sneezing,
                 family = binomial(link = "logit"),
                 data = flu_train)
summary(all_model)

# Model with all symptoms modified by age
age_model <- glm(flu ~ as.factor(age_group) + fever*as.factor(age_group) + 
                   cough*as.factor(age_group) + sore_throat*as.factor(age_group) +
                    myalgia*as.factor(age_group) + short_breath*as.factor(age_group) + 
                    hypoxemia*as.factor(age_group) +
                    nausea_vom*as.factor(age_group) + bronchitis*as.factor(age_group) + 
                    chest_pain*as.factor(age_group) + 
                    diarrhea*as.factor(age_group) + fatigue*as.factor(age_group) + 
                    headache*as.factor(age_group) + 
                    congestion*as.factor(age_group) + sneezing*as.factor(age_group),
                  family = binomial(link = "logit"),
                  data = flu_train)
summary(age_model)

# Model with all symptoms modified by age in addition to time
time_model <- glm(flu ~ as.factor(age_group) + fever*as.factor(age_group) + cough*as.factor(age_group) + 
                     sore_throat*as.factor(age_group) +
                     myalgia*as.factor(age_group) + short_breath*as.factor(age_group) + 
                     hypoxemia*as.factor(age_group) +
                     nausea_vom*as.factor(age_group) + bronchitis*as.factor(age_group) + 
                     chest_pain*as.factor(age_group) + 
                     diarrhea*as.factor(age_group) + fatigue*as.factor(age_group) + 
                     headache*as.factor(age_group) + 
                     congestion*as.factor(age_group) + sneezing*as.factor(age_group) + as.factor(month),
                   family = binomial(link = "logit"),
                   data = flu_train)
summary(time_model)

# Add predictions to the testing data
# Fever model
flu_test$pred_fev <- predict(fev_model, 
                          newdata = flu_test, 
                          type = "response")
unbalanced$pred_fev <- predict(fev_model, 
                         newdata = unbalanced, 
                         type = "response")
# All symptom model
flu_test$pred_all <- predict(all_model, 
                             newdata = flu_test, 
                             type = "response")
unbalanced$pred_all <- predict(all_model, 
                               newdata = unbalanced, 
                               type = "response")

# Age model
flu_test$pred_age <- predict(age_model, 
                             newdata = flu_test, 
                             type = "response")

unbalanced$pred_age <- predict(age_model, 
                               newdata = unbalanced, 
                               type = "response")

# Time model
flu_test$pred_time <- predict(time_model, 
                             newdata = flu_test, 
                             type = "response")

unbalanced$pred_time <- predict(time_model, 
                               newdata = unbalanced, 
                               type = "response")

# Check AUC values for all models
pROC::auc(pROC::roc(flu_test$flu, flu_test$pred_fev))
pROC::auc(pROC::roc(unbalanced$flu, unbalanced$pred_fev))
pROC::auc(pROC::roc(flu_test$flu, flu_test$pred_all))
pROC::auc(pROC::roc(unbalanced$flu, unbalanced$pred_all))
pROC::auc(pROC::roc(flu_test$flu, flu_test$pred_age))
pROC::auc(pROC::roc(unbalanced$flu, unbalanced$pred_age))
pROC::auc(pROC::roc(flu_test$flu, flu_test$pred_time))
pROC::auc(pROC::roc(unbalanced$flu, unbalanced$pred_time))

# Check densities of predictions
ggplot(unbalanced) + geom_density(aes(x = pred_fev, group = as.factor(flu), color = as.factor(flu)))
ggplot(unbalanced) + geom_density(aes(x = pred_all, group = as.factor(flu), color = as.factor(flu)))
ggplot(unbalanced) + geom_density(aes(x = pred_age, group = as.factor(flu), color = as.factor(flu)))
ggplot(unbalanced) + geom_density(aes(x = pred_time, group = as.factor(flu), color = as.factor(flu)))

# Save model results
saveRDS(fev_model, "./tmp/fever_model.rds")
saveRDS(all_model, "./tmp/symptoms_model.rds")
saveRDS(age_model, "./tmp/age_model.rds")
saveRDS(time_model, "./tmp/time_model.rds")

cp <- cutpointr(unbalanced, pred_fev, flu, 
                method = maximize_metric,
                metric = youden)
summary(cp)

cp <- cutpointr(unbalanced, pred_all, flu, 
                method = maximize_metric,
                metric = youden)
summary(cp)

cp <- cutpointr(unbalanced, pred_age, flu, 
                method = maximize_metric,
                metric = youden)
summary(cp)

cp <- cutpointr(unbalanced, pred_time, flu, 
                method = maximize_metric,
                metric = youden)
summary(cp)

# Load part 2 data
flu_proc_symp <- readRDS('./tmp/flu_p2_dat_proc_symp.rds')

# Create empty lists to fill
flu_nat_filt <- list()

for (i in 1:4) {
  # Print the progress and set season data
  print(i)
  flu <- as.data.frame(flu_proc_symp[i])
  
  # add month and age variables
  flu$month <- month(flu$week_date)
  flu$age_group <- factor(flu$age_grp,
                          levels = c('0', '1', '2', '3', '4', '5'),
                          labels = c('0-4', '5-12', '13-17', '18-49', '50-64', '65+'))
  
  # Restrict data to the influenza season
  flu <- flu |> dplyr::filter(month > 9 | month < 6)
  
  # Collapse to the national level
  flu_nat <- flu |> group_by(flu, fever, cough, sore_throat, myalgia, 
                             short_breath, hypoxemia, nausea_vom, 
                             bronchitis, chest_pain, diarrhea,
                             fatigue, headache, congestion, sneezing,
                             age_group, month) |>
    mutate(patient_sum = sum(patient_count_imp)) |>
    distinct(flu, fever, cough, sore_throat, myalgia, 
             short_breath, hypoxemia, nausea_vom, 
             bronchitis, chest_pain, diarrhea,
             fatigue, headache, congestion, sneezing,
             age_group, month, patient_sum)
  
  # Fill lists
  flu_nat_filt[[i]] <- as.data.frame(flu_nat)
  remove(flu, flu_nat)
}

flu_16_17 <- flu_nat_filt[[1]]
flu_16_17$pred_fev <- predict(fev_model, 
                             newdata = flu_16_17, 
                             type = "response")
flu_16_17$pred_all <- predict(all_model, 
                              newdata = flu_16_17, 
                              type = "response")
flu_16_17$pred_age <- predict(age_model, 
                              newdata = flu_16_17, 
                              type = "response")
flu_16_17$pred_time <- predict(time_model, 
                              newdata = flu_16_17, 
                              type = "response")
flu_16_17$pred_fev_count <- flu_16_17$pred_fev * flu_16_17$patient_sum
flu_16_17$pred_all_count <- flu_16_17$pred_all * flu_16_17$patient_sum
flu_16_17$pred_age_count <- flu_16_17$pred_age * flu_16_17$patient_sum
flu_16_17$pred_time_count <- flu_16_17$pred_time * flu_16_17$patient_sum

sum(flu_16_17$pred_time_count) * runif(10000, min = 1/0.58, max = 1/0.42)


cp <- cutpointr(unbalanced, pred, flu, 
                method = maximize_metric,
                metric = youden)

summary(cp)



flu_symp_p2$season <- ifelse(flu_symp_p2$month > 8 | flu_symp_p2$month < 6, 1, 0)

flu_symp_p2$pred <-  predict(train_model, 
                                 newdata = flu_symp_p2, 
                                 type = "response")

flu_symp_p2$ili_flu <- ifelse(flu_symp_p2$fever == 1 & flu_symp_p2$cough == 1, 1, 0)
flu_symp_p2$ili_flu <- ifelse(flu_symp_p2$fever == 1 & flu_symp_p2$sore_throat == 1, 1, flu_symp_p2$ili_flu)

flu_symp_p2$pred_flu <- ifelse(flu_symp_p2$pred >  0.4907, 1, 0)


flu_symp_p2$systemic <- ifelse(flu_symp_p2$fever == 1 | flu_symp_p2$fatigue == 1 | flu_symp_p2$myalgia == 1, 1, 0)
flu_symp_p2$resp <- ifelse(flu_symp_p2$cough == 1 | flu_symp_p2$sore_throat == 1 | flu_symp_p2$short_breath == 1, 1, 0)

flu_symp_p2$ecdc <- ifelse(flu_symp_p2$systemic == 1 & flu_symp_p2$resp == 1, 1, 0)

flu_symp_p2_filt <- flu_symp_p2 |> filter(season == 1)

p2_county_week <- flu_symp_p2_filt |> group_by(county_fips, week_date) |>
  mutate(flu = sum(flu*patient_count_imp),
         fever = sum(fever*patient_count_imp),
         cough = sum(cough*patient_count_imp),
         sore_throat = sum(sore_throat*patient_count_imp),
         fatigue = sum(fatigue*patient_count_imp),
         myalgia = sum(myalgia*patient_count_imp),
         hypoxemia = sum(hypoxemia*patient_count_imp),
         short_breath = sum(short_breath*patient_count_imp),
         bronchitis = sum(bronchitis*patient_count_imp),
         chest_pain = sum(chest_pain*patient_count_imp),
         nausea_vom = sum(nausea_vom*patient_count_imp),
         headache = sum(headache*patient_count_imp),
         diarrhea = sum(diarrhea*patient_count_imp),
         congestion = sum(congestion*patient_count_imp),
         sneezing = sum(sneezing*patient_count_imp),
         flu_pred = sum(pred_flu * patient_count_imp),
         all_cause = sum(patient_count_imp)) |>
  distinct(county_fips, week_date, month, flu, fever,
           myalgia, hypoxemia, short_breath, cough, bronchitis, chest_pain, 
           nausea_vom, sore_throat, fatigue, diarrhea, headache, congestion, 
           sneezing, flu_pred, all_cause, .keep_all = FALSE)

ggplot(p2_county_week) + geom_line(aes(y = flu / flu_pred , x = week_date), color = 'black') +
# geom_line(aes(y = flu , x = week_date), color = 'red') +
  facet_wrap(vars(county_fips), scales = "free")


  
p2_county_week$year <- lubridate::year(p2_county_week$week_date)

p2_county_week$season <- ifelse(p2_county_week$year == 2016 & p2_county_week$month > 8, '2016-2017', '')
p2_county_week$season <- ifelse(p2_county_week$year == 2017 & p2_county_week$month < 9, '2016-2017', p2_county_week$season)
p2_county_week$season <- ifelse(p2_county_week$year == 2017 & p2_county_week$month > 8, '2017-2018', p2_county_week$season)
p2_county_week$season <- ifelse(p2_county_week$year == 2018 & p2_county_week$month < 9, '2017-2018', p2_county_week$season)
p2_county_week$season <- ifelse(p2_county_week$year == 2018 & p2_county_week$month > 8, '2018-2019', p2_county_week$season)
p2_county_week$season <- ifelse(p2_county_week$year == 2019 & p2_county_week$month < 9, '2018-2019', p2_county_week$season)
p2_county_week$season <- ifelse(p2_county_week$year == 2019 & p2_county_week$month > 8, '2019-2020', p2_county_week$season)
p2_county_week$season <- ifelse(p2_county_week$year == 2020 & p2_county_week$month < 9, '2019-2020', p2_county_week$season)


p2_county_perc_all <- p2_county_week |> dplyr::filter(month > 9 | month < 5) |>
  group_by(county_fips, season) |>
  mutate(sum_pred = sum(flu_pred),
         sum_all_cause = sum(all_cause)) |>
    distinct(county_fips, season, sum_pred, sum_all_cause, .keep_all = FALSE)

p2_county_perc_all$perc <- p2_county_perc_all$sum_pred / p2_county_perc_all$sum_all_cause
p2_county_perc_all$Cat <- 'All Symptoms'


p2_county_perc <- rbind(p2_county_perc_obs, p2_county_perc_null, p2_county_perc_fever,
                        p2_county_perc_all, p2_county_perc_all_age, p2_county_perc_all_age_time)


library(ggridges)
ggplot(p2_county_perc, aes(x = perc, y = Cat, fill = Cat)) + geom_density_ridges(scale = 1)



# Convert all count data to log(count)
p2_county_week$flu_log <- log(p2_county_week$flu + 1)
p2_county_week$fever_log <- log(p2_county_week$fever + 1)
p2_county_week$cough_log <- log(p2_county_week$cough + 1)
p2_county_week$sore_throat_log <- log(p2_county_week$sore_throat + 1)
p2_county_week$fatigue_log <- log(p2_county_week$fatigue + 1)
p2_county_week$myalgia_log <- log(p2_county_week$myalgia + 1)
p2_county_week$hypoxemia_log <- log(p2_county_week$hypoxemia + 1)
p2_county_week$short_breath_log <- log(p2_county_week$short_breath + 1)
p2_county_week$bronchitis_log <- log(p2_county_week$bronchitis + 1)
p2_county_week$chest_pain_log <- log(p2_county_week$chest_pain + 1)
p2_county_week$nausea_vom_log <- log(p2_county_week$nausea_vom + 1)
p2_county_week$headache_log <- log(p2_county_week$headache + 1)
p2_county_week$diarrhea_log <- log(p2_county_week$diarrhea + 1)
p2_county_week$congestion_log <- log(p2_county_week$congestion + 1)
p2_county_week$sneezing_log <- log(p2_county_week$sneezing + 1)


model <- lm(flu_log ~ fever_log + cough_log + sore_throat_log + 
      myalgia_log + short_breath_log + hypoxemia_log + nausea_vom_log + 
      bronchitis_log + chest_pain_log + diarrhea_log + fatigue_log + headache_log + 
      congestion_log + sneezing_log,
    data = p2_county_week)

summary(model)

p2_county_week$pred <-  predict(model, 
                             newdata = p2_county_week)
p2_county_week$flu_lin <- exp(p2_county_week$pred) - 1


ggplot(p2_county_week) + geom_line(aes(y = flu, x = week_date), color = 'red') +
  geom_line(aes(y = flu_lin, x = week_date), color = 'blue')





flu_p2 <- readRDS('./tmp/flu_p2_dat_proc.rds')
p2_county_week <- flu_p2 |> dplyr::filter(symp_count > 0) |>
  group_by(county_fips, week_date) |>
  mutate(flu = sum(flu*patient_count_imp),
         fever = sum(fever*patient_count_imp),
         cough = sum(cough*patient_count_imp),
         sore_throat = sum(sore_throat*patient_count_imp),
         fatigue = sum(fatigue*patient_count_imp),
         myalgia = sum(myalgia*patient_count_imp),
         hypoxemia = sum(hypoxemia*patient_count_imp),
         short_breath = sum(short_breath*patient_count_imp),
         bronchitis = sum(bronchitis*patient_count_imp),
         chest_pain = sum(chest_pain*patient_count_imp),
         nausea_vom = sum(nausea_vom*patient_count_imp),
         headache = sum(headache*patient_count_imp),
         diarrhea = sum(diarrhea*patient_count_imp),
         congestion = sum(congestion*patient_count_imp),
         sneezing = sum(sneezing*patient_count_imp),
         flu_pred = sum(pred_flu * patient_count_imp)) |>
  distinct(county_fips, week_date, month, flu, fever, 
           myalgia, hypoxemia, short_breath, cough, bronchitis, chest_pain, 
           nausea_vom, sore_throat, fatigue, diarrhea, headache, congestion, 
           sneezing, flu_pred, .keep_all = FALSE)

saveRDS(p2_county_week, './tmp/p2_county_week_full.rds')








p2_county_week <- flu_symp_p2 |> group_by(county_fips, week_date) |>
  mutate(flu_obs_tot = sum(flu * patient_count_imp),
         flu_pred_tot = sum(pred_flu * patient_count_imp),
         flu_pred_tot_TEST = sum(pred * patient_count_imp),
         average_pred = mean(pred)) |>
  distinct(county_fips, week_date, flu_obs_tot, flu_pred_tot, flu_pred_tot_TEST,
           average_pred)

long <- p2_county_week %>% 
  pivot_longer(
    cols = flu_obs_tot:average_pred, 
    names_to = "Cat",
    values_to = "value"
  )

all_week <- read.csv('./raw/county_week_ac_v2.csv')
all_season <- read.csv('./raw/county_season_ac_v2.csv')
all_season_avg <- all_season |> group_by(county_fips) |>
  mutate(avg_all = mean(all_cause)) |>
  distinct(county_fips, avg_all)
  

library(ISOweek)
all_week$week_date <- sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", all_week$year_week)
all_week$week_date <- ISOweek2date(all_week$week_date)

long <- left_join(long, all_week, by = c('county_fips' = 'county_fips',
                                         'week_date' = 'week_date'))
long <- left_join(long, all_season_avg, by = c('county_fips' = 'county_fips'))

long$percent <- long$value / as.numeric(long$all_cause)
long$percent_2 <- long$value / as.numeric(long$avg_all)

ggplot(long) + geom_line(aes(y = value, x = week_date, group = Cat, color = Cat)) +
  facet_wrap(vars(county_fips), scales = "free") 


p2_county_week$percent <- p2_county_week$flu_obs_tot / p2_county_week$flu_pred_tot
p2_county_week_filt <- p2_county_week |> filter(week_date > '2016-10-24') |>
  filter(week_date < '2017-04-03')
ggplot(p2_county_week_filt) + geom_line(aes(y = flu_obs_tot / flu_pred_tot, x = week_date)) +
  facet_wrap(vars(county_fips), scales = "free") 


p2_county_week_avg <- flu_symp_p2 |> group_by(county_fips, week_date) |>
  mutate(average_pred = mean(pred),
         average_pred_bal = mean(pred_bal),
         pred_25 = quantile(pred, probs = c(0.25)),
         pred_25_bal = quantile(pred_bal, probs = c(0.25)),
         pred_75 = quantile(pred, probs = c(0.75)),
         pred_75_bal = quantile(pred_bal, probs = c(0.75))) |>
  distinct(county_fips, week_date,
           average_pred, average_pred_bal, pred_25,
           pred_25_bal, pred_75, pred_75_bal)

long_avg <- p2_county_week_avg %>% dplyr::select(c(county_fips, week_date,
                                        average_pred, average_pred_bal)) |>
  pivot_longer(
    cols = average_pred:average_pred_bal, 
    names_to = "Cat",
    values_to = "value"
  )

ggplot(long_avg) + 
  facet_wrap(vars(county_fips), scales = "free") +
  geom_ribbon(data = p2_county_week_avg, aes(x = week_date,
                                             ymin = pred_25,
                                             ymax = pred_75),
              fill = '#F8766D') + 
  geom_ribbon(data = p2_county_week_avg, aes(x = week_date,
                                             ymin = pred_25_bal,
                                             ymax = pred_75_bal),
              fill = '#00BFC4')




flu_symp_p1_test$pred_test <- ifelse(flu_symp_p1_test$pred > i, 1, 0)


for (i in seq(0, 1, 0.01)) {
  flu_symp_p1_test$pred_test <- ifelse(flu_symp_p1_test$pred > i, 1, 0)
  flu_symp_p1_test$correct <- ifelse(flu_symp_p1_test$pred_test == 
                                       flu_symp_p1_test$flu, 1, 0)
  accuracy
}



library(caret)
conf_matrix <- confusionMatrix(flu_symp_p1_test$flu, flu_symp_p1_test$pred)


unbalanced <- flu_symp_p1 |> filter(county_fips == "36061") |>
  uncount(patient_count_imp) |> sample_n(100000)

train_model_nyc <- glm(flu ~ fever + cough + sore_throat +
                             myalgia + short_breath + hypoxemia +
                             nausea_vom + bronchitis + chest_pain + 
                             diarrhea + fatigue + headache + 
                             congestion + sneezing + as.factor(month),
                           family = binomial(link = "logit"),
                           data = unbalanced)

unbalanced <- flu_symp_p1 |> filter(county_fips == "06037") |>
  uncount(patient_count_imp) |> sample_n(100000)

train_model_la <- glm(flu ~ fever + cough + sore_throat +
                         myalgia + short_breath + hypoxemia +
                         nausea_vom + bronchitis + chest_pain + 
                         diarrhea + fatigue + headache + 
                         congestion + sneezing + as.factor(month),
                       family = binomial(link = "logit"),
                       data = unbalanced)

unbalanced <- flu_symp_p1 |> filter(county_fips == "11001") |>
  uncount(patient_count_imp) |> sample_n(100000)

train_model_dc <- glm(flu ~ fever + cough + sore_throat +
                        myalgia + short_breath + hypoxemia +
                        nausea_vom + bronchitis + chest_pain + 
                        diarrhea + fatigue + headache + 
                        congestion + sneezing + as.factor(month),
                      family = binomial(link = "logit"),
                      data = unbalanced)

flu_symp_p2$pred_nyc <-  predict(train_model_nyc, 
                                 newdata = flu_symp_p2, 
                                 type = "response")

flu_symp_p2$pred_la <-  predict(train_model_la, 
                                 newdata = flu_symp_p2, 
                                 type = "response")

flu_symp_p2$pred_dc <-  predict(train_model_dc, 
                                 newdata = flu_symp_p2, 
                                 type = "response")


p2_county_week <- flu_symp_p2 |> group_by(county_fips, week_date) |>
  mutate(flu_obs_tot = sum(flu * patient_count_imp),
         flu_pred_nyc = sum(pred_nyc * patient_count_imp),
         flu_pred_la = sum(pred_la * patient_count_imp),
         flu_pred_dc = sum(pred_dc * patient_count_imp)) |>
  distinct(county_fips, week_date, flu_obs_tot, flu_pred_nyc, flu_pred_la, flu_pred_dc)

long <- p2_county_week %>% 
  pivot_longer(
    cols = flu_obs_tot:flu_pred_dc, 
    names_to = "Cat",
    values_to = "value"
  )

ggplot(long) + geom_line(aes(y = value, x = week_date, group = Cat, color = Cat)) +
  facet_wrap(vars(county_fips), scales = "free") 








num <- 10000
loop <- NULL
list <- c( "17031", "53033", "06037", "04013", "36061", 
           "48113", "11001", "13121", "29510")
for (i in 1:100) {
  
  unbalanced <- flu_symp_p1_train |> 
    uncount(patient_count_imp) |> sample_n(num)
  
  train_model_symptom <- glm(flu ~ fever + cough + sore_throat +
                               myalgia + short_breath + hypoxemia +
                               nausea_vom + bronchitis + chest_pain + 
                               diarrhea + fatigue + headache + 
                               congestion + sneezing,
                             family = binomial(link = "logit"),
                             data = unbalanced)
  test <- data.frame(t(train_model_symptom$coefficients))
  test$num <- num
  loop <- rbind(loop, test)
  num <- num + 10000
}

long <- loop %>% 
  pivot_longer(
    cols = `X.Intercept.`:`sneezing`, 
    names_to = "coeff",
    values_to = "value"
  )

ggplot(long) + geom_line(aes(y = value, x = num, group = coeff)) +
  facet_wrap(vars(coeff), scales = "free")



loop <- NULL
list <- c( "17031", "53033", "06037", "04013", "36061", 
           "48113", "11001", "13121")


flu_symp_p1 |> group_by(county_fips) |>
  mutate(sum = sum(patient_count_imp)) |>
  distinct(county_fips, sum)

for (i in list) {
  unbalanced <- flu_symp_p1 |> filter(county_fips == i) |>
    uncount(patient_count_imp) |> sample_n(100000)
  
  train_model_symptom <- glm(flu ~ fever + cough + sore_throat +
                               myalgia + short_breath + hypoxemia +
                               nausea_vom + bronchitis + chest_pain + 
                               diarrhea + fatigue + headache + 
                               congestion + sneezing,
                             family = binomial(link = "logit"),
                             data = unbalanced)
  
  test <- data.frame(t(train_model_symptom$coefficients))
  test$county <- i
  loop <- rbind(loop, test)
}

long <- loop %>% 
  pivot_longer(
    cols = `X.Intercept.`:`sneezing`, 
    names_to = "coeff",
    values_to = "value"
  )

ggplot(long) + geom_point(aes(y = value, x = county, group = coeff, color = coeff), size = 1.5) +
  facet_wrap(vars(coeff), scales = "free") + theme(legend.position = 'none')


# Sample from the full Part 1 data
# Unbalanced data
unbalanced <- flu_symp_p2_train |> 
  uncount(patient_count_imp) |> sample_n(100000)

cases <- flu_symp_p2_train |> 
  uncount(patient_count_imp) |> sample_n(50000)

non_cases <- flu_symp_p2_train |> 
  uncount(patient_count_imp) |> sample_n(50000)

balanced <- rbind(cases, non_cases)
  
train_model_symptom <- glm(flu ~ 0 + fever + cough + sore_throat +
                             myalgia + short_breath + hypoxemia +
                             nausea_vom + bronchitis + chest_pain + 
                             diarrhea + fatigue + headache + 
                             congestion + sneezing + fever_perc_lag + 
                             cough_perc_lag + sore_throat_perc_lag,
                           family = binomial(link = "logit"),
                           data = balanced)
summary(train_model_symptom)
train_model_fever <- glm(flu ~ fever,
                           family = binomial(link = "logit"),
                           data = unbalanced)

train_model_symptom_age <- glm(flu ~ fever*as.factor(age_grp) + cough*as.factor(age_grp) + 
                                 sore_throat*as.factor(age_grp) +
                                 + myalgia*as.factor(age_grp) + short_breath*as.factor(age_grp) + 
                                 hypoxemia*as.factor(age_grp) +
                                 nausea_vom *as.factor(age_grp)+ bronchitis*as.factor(age_grp) + 
                                 chest_pain*as.factor(age_grp) + 
                                 diarrhea*as.factor(age_grp) + fatigue*as.factor(age_grp) + 
                                 headache*as.factor(age_grp) + congestion*as.factor(age_grp) + 
                                 sneezing*as.factor(age_grp),
                               family = binomial(link = "logit"),
                               data = unbalanced)

train_model_time <- glm(flu ~ fever + cough + sore_throat +
                             myalgia + short_breath + hypoxemia +
                             nausea_vom + bronchitis + chest_pain + 
                             diarrhea + fatigue + headache + 
                             congestion + sneezing + jan + feb + mar + 
                          apr + may + jun + jul + aug + sep + oct + nov + dec,
                           family = binomial(link = "logit"),
                           data = unbalanced)

summary(train_model_symptom)
summary(train_model_fever)
summary(train_model_symptom_age)
summary(train_model_time)

test <- flu_symp_p2_test |> 
  uncount(patient_count_imp) |> sample_n(100000)

test$pred_age <- predict(train_model_symptom_age, 
                          newdata = test, 
                          type = "response")

test$pred_fever <- predict(train_model_fever, 
                           newdata = test, 
                           type = "response")

test$pred_symp <- predict(train_model_symptom, 
                          newdata = test, 
                          type = "response")

test$pred_time <- predict(train_model_time, 
                          newdata = test, 
                          type = "response")
library(plotROC)
library(pROC)

pROC::auc(pROC::roc(test$flu, test$pred_age))
pROC::auc(pROC::roc(test$flu, test$pred_fever))
pROC::auc(pROC::roc(test$flu, test$pred_symp))
pROC::auc(pROC::roc(test$flu, test$pred_time))


test$diff <- test$flu - test$pred_symp

ggplot(test) + geom_density(aes(x = pred_time, group = flu, color = flu))
mean(test[test$flu == 1, ]$pred_symp)
mean(test[test$flu == 0, ]$pred_symp, na.rm = TRUE)

library(cutpointr)
opt_cut_1617_youden <- cutpointr(test, pred_symp, flu, direction = ">=", 
                                 pos_class = 1, neg_class = 0, 
                                 method = maximize_metric, metric = youden, 
                                 na.rm = TRUE)
opt_cut_1617_youden


library(subsampling)

flu_p1_symp_expand <- flu_symp_p1 |> uncount(patient_count_imp)


ssp.relogit(
  formula = flu ~ fever + cough + sore_throat +
    myalgia + short_breath + hypoxemia +
    nausea_vom + bronchitis + chest_pain + 
    diarrhea + fatigue + headache + 
    congestion + sneezing,
  data = flu_p1_symp_expand,
  n.plt = 200,
  n.ssp = 100000,
  criterion = "optL",
  likelihood = "logOddsCorrection"
)









