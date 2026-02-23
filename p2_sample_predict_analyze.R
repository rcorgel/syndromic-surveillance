
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

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

# Load processed data
flu <- readRDS('./tmp/flu_p1_dat_proc.rds')
flu_p2 <- readRDS('./tmp/flu_p2_dat_proc_select_counties.rds')
flu_p2$month <- month(flu_p2$week_date)
flu_p2$year <- year(flu_p2$week_date)
flu_p2$month_date <- as.Date(paste(flu_p2$month, '01', flu_p2$year, sep = '-'), "%m-%d-%Y")

# Create lagged fever, sore throat, and cough percent variables
lagged_vars <- flu_p2 |> group_by(county_fips, month_date) |>
  mutate(fever_perc_lag = sum(fever*patient_count_imp) / sum(patient_count_imp),
         cough_perc_lag = sum(cough*patient_count_imp) / sum(patient_count_imp),
         sore_throat_perc_lag = sum(sore_throat*patient_count_imp) / sum(patient_count_imp)) |>
  distinct(county_fips, month_date, fever_perc_lag, cough_perc_lag, sore_throat_perc_lag, .keep_all = FALSE) |>
  mutate(lag_month_date = month_date %m+% months(1)) |> ungroup() |>
  select(-c(month_date))

# Merge on lagged variables
flu_p2 <- left_join(flu_p2, lagged_vars, by = c('county_fips' = 'county_fips',
                                                    'month_date' = 'lag_month_date'))

# Remove asymptomatic cases for disease and non-diseased
flu_symp <- flu |> dplyr::filter(symp_count > 0)

flu_county_select <- flu_symp |> 
  dplyr::filter(county_fips == "17031" | county_fips == "53033" | 
           county_fips == "06037" | county_fips == "04013" |
           county_fips == "48113" | county_fips == "36061" |
           county_fips == "11001" | county_fips == "13121" |
           county_fips == "29510")

# Bring data back to season level
flu_16_17 <- flu_county_select |>
  filter(month_date > as.Date("2016-08-31")) |>
  filter(month_date < as.Date("2017-09-01"))

flu_16_17$month <- month(flu_16_17$month_date)

nyc <- flu_p2 |> dplyr::filter(county_fips == "36061") |>
  filter(symp_count > 0)
flu_16_17 <- flu_16_17 |> dplyr::filter(!county_fips == "36061")

# Restrict data to flu season months (September - April
# Split into cases and non-cases
# 16-17
cases_16_17 <- flu_16_17[flu_16_17$flu == 1, ]
non_cases_16_17 <- flu_16_17[flu_16_17$flu == 0, ]

month <- cases_16_17 |> group_by(month) |>
  mutate(sum = sum(patient_count_imp)) |>
  distinct(month, sum)

# Expand case data and save
cases_exp_sub_16_17 <- cases_16_17 |> uncount(patient_count_imp) |> 
  sample_n(10000)

# Expand and subset non-case data
non_cases_exp_sub_16_17 <- non_cases_16_17 |> uncount(patient_count_imp) |> 
  sample_n(10000)

as_is <- flu_16_17 |> uncount(patient_count_imp) |> 
  sample_n(30000)

flu_subset_ind_1 <- rbind(cases_exp_sub_16_17 , non_cases_exp_sub_16_17) 

sample <- sample(c(TRUE, FALSE), nrow(as_is), replace=TRUE, prob=c(0.5, 0.5))
flu_train  <- as_is[sample, ]
flu_test   <- as_is[!sample, ]

train_model <- glm(flu ~ fever + cough + sore_throat +
                     + myalgia + short_breath + hypoxemia +
                     nausea_vom + bronchitis + chest_pain + 
                     diarrhea + fatigue + headache + 
                     congestion + sneezing + cough_perc_lag,
                   family = binomial(link = "logit"),
                   data = flu_train)

summary(train_model)
vif(train_model)

# Predict testing data
flu_test$pred <- predict(train_model, newdata = flu_test, type = "response")


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
# coeffitients 
pROC::auc(pROC::roc(flu_test$flu, flu_test$pred))

nyc_uncount <- nyc |> uncount(patient_count_imp)
  
nyc_uncount$pred <- predict(train_model, newdata = nyc_uncount, type = "response", se.fit = TRUE)$fit
nyc_uncount$se <- predict(train_model, newdata = nyc_uncount, type = "response", se.fit = TRUE)$se.fit
nyc_uncount$upr <- nyc_uncount$pred + (1.96 * nyc_uncount$se)
nyc_uncount$lwr <- nyc_uncount$pred - (1.96 * nyc_uncount$se)

#test <- predict(train_model, newdata = nyc_uncount, type = "response", se.fit = TRUE)

nyc_test <- nyc_uncount |> group_by(week_date) |> 
  mutate(sum_flu = sum(flu),
         sum_pred = sum(pred),
         sum_up = sum(upr),
         sum_low = sum(lwr)) |>
  distinct(week_date, sum_flu, sum_pred, sum_up, sum_low, .keep_all = FALSE)

ggplot(nyc_test) +
  geom_line(aes(x = week_date, y = sum_flu), color = 'red') +
  geom_line(aes(x = week_date, y = sum_pred), color = 'blue') +
  geom_line(aes(x = week_date, y = sum_up), color = 'blue') +
  geom_line(aes(x = week_date, y = sum_low), color = 'blue')





nyc_symp_count <- nyc |> 
  group_by(week_date) |>
  mutate(fever = sum(fever*patient_count_imp),
         myalgia = sum(myalgia*patient_count_imp),
         cough = sum(cough*patient_count_imp),
         sore_throat = sum(sore_throat*patient_count_imp), 
         short_breath = sum(short_breath*patient_count_imp), 
         hypoxemia = sum(hypoxemia*patient_count_imp), 
         chest_pain = sum(chest_pain*patient_count_imp), 
         bronchitis = sum(bronchitis*patient_count_imp),
         nausea_vom = sum(nausea_vom*patient_count_imp), 
         diarrhea = sum(diarrhea*patient_count_imp), 
         fatigue = sum(fatigue*patient_count_imp), 
         headache = sum(headache*patient_count_imp),
         congestion = sum(congestion*patient_count_imp), 
         sneezing = sum(sneezing*patient_count_imp)) |>
  distinct(week_date, fever, myalgia, cough, sore_throat, 
           short_breath, hypoxemia, chest_pain, bronchitis,
           nausea_vom, diarrhea, fatigue, headache,
           congestion, sneezing) |>
  pivot_longer(cols = fever:sneezing,
               names_to = "symptom",
               values_to = "value")


ggplot(nyc_symp_count) +
  geom_line(aes(x = week_date, y = value, color = symptom, group = symptom))



nyc_symp <- nyc_uncount |> group_by(week_date) |> 
  mutate(sum_fever = sum(fever),
         sum_st = sum(sore_throat),
         sum_cough = sum(cough)) |>
  distinct(week_date, sum_fever, sum_st, sum_cough, .keep_all = FALSE)


long <- wide %>% 
  pivot_longer(
    cols = `1950`:`1954`, 
    names_to = "year",
    values_to = "value"
  )



ggplot(nyc_symp) +
  geom_line(aes(x = week_date, y = sum_fever), color = 'red') +
  geom_line(aes(x = week_date, y = sum_st), color = 'blue') +
  geom_line(aes(x = week_date, y = sum_cough), color = 'green') +
  ylim(0, 3000)








