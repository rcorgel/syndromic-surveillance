################################################################################
# File Name: analyze_flu_data                                                  #
#                                                                              #
# Purpose:   Analyze influenza from medical claims symptoms.                   #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load data                                                      #
#            3. Analyze data                                                   #
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
library(cutpointr)
library(ISOweek)
library(zoo)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

################
# 2. LOAD DATA #
################

# Load models and cut offs
load('./tmp/flu_model_cut_off_results.RData')

mod_fever_combos <- readRDS('./tmp/training_models_all_three_other_resp.rds')
mod_fever <- readRDS('./tmp/training_model_time.rds')
mod_fever_cough_sore <- readRDS('./tmp/training_models_all_three_the_rest.rds')


# Load p2 data
flu_16_17 <- readRDS('./tmp/flu_p2_dat_proc.rds')

flu_16_17_county_select <- flu_16_17 |> 
  filter(county_fips == "17031" | county_fips == "53033" | 
           county_fips == "06037" | county_fips == "04013" |
           county_fips == "36061" | county_fips == "48113" |
           county_fips == "11001" | county_fips == "13121" |
           county_fips == "29510")

saveRDS(flu_16_17_county_select, './tmp/flu_p2_dat_proc_select_counties.rds')

flu_16_17 <- flu_16_17 |>
  filter(week_date > as.Date("2016-08-31")) |>
  filter(week_date < as.Date("2017-09-01"))



# Make ili var
flu_16_17 <- flu_16_17 |> mutate(ili = ifelse((fever == 1 & cough ==1) |  (fever == 1 & sore_throat == 1) |  (fever == 1 & sore_throat == 1 & cough ==1), 1, 0),
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
                                 fever_cough_myalgia = ifelse((fever == 1 & cough == 1 & myalgia == 1), 1, 0))

flu_16_17$month <- month(flu_16_17$week_date)
flu_16_17$jan <- ifelse(flu_16_17$month == 1, 1, 0)
flu_16_17$feb <- ifelse(flu_16_17$month == 2, 1, 0)
flu_16_17$mar <- ifelse(flu_16_17$month == 3, 1, 0)
flu_16_17$apr <- ifelse(flu_16_17$month == 4, 1, 0)
flu_16_17$may <- ifelse(flu_16_17$month == 5, 1, 0)
flu_16_17$jun <- ifelse(flu_16_17$month == 6, 1, 0)
flu_16_17$jul <- ifelse(flu_16_17$month == 7, 1, 0)
flu_16_17$aug <- ifelse(flu_16_17$month == 8, 1, 0)
flu_16_17$sep <- ifelse(flu_16_17$month == 9, 1, 0)
flu_16_17$oct <- ifelse(flu_16_17$month == 10, 1, 0)
flu_16_17$nov <- ifelse(flu_16_17$month == 11, 1, 0)
flu_16_17$dec <- ifelse(flu_16_17$month == 12, 1, 0)

# Merge on weather data
#weather <- readRDS('./tmp/county_week_weather.rds')
#flu_16_17 <- left_join(flu_16_17, weather, by = c('county_fips' = 'county_group',
                                                  #'week_date' = 'week_date'))

# Create lagged fever, sore throat, and cough percent variables
lagged_vars <- flu_16_17 |> group_by(county_fips, week_date) |>
  mutate(fever_perc_lag = sum(fever*patient_count_imp) / sum(patient_count_imp),
         cough_perc_lag = sum(cough*patient_count_imp) / sum(patient_count_imp),
         sore_throat_perc_lag = sum(sore_throat*patient_count_imp) / sum(patient_count_imp)) |>
  distinct(county_fips, week_date, fever_perc_lag, cough_perc_lag, sore_throat_perc_lag, .keep_all = FALSE) |>
  mutate(lag_week_date = week_date - 7) |> ungroup() |>
  select(-c(week_date))

# Merge on lagged variables
flu_16_17 <- left_join(flu_16_17, lagged_vars, by = c('county_fips' = 'county_fips',
                                                    'week_date' = 'lag_week_date'))

# Break up data into two sections
sample <- sample(c(TRUE, FALSE), nrow(flu_16_17), replace=TRUE, prob=c(0.5, 0.5))
flu_16_17_1  <- flu_16_17[sample, ]
flu_16_17_2   <- flu_16_17[!sample, ]
remove(flu_16_17)

# Predict using the model
flu_16_17_1$pred_com <- predict(mod_fever_combos, newdata = flu_16_17_1, type = "response")
flu_16_17_2$pred_com <- predict(mod_fever_combos, newdata = flu_16_17_2, type = "response")

flu_16_17_1$pred_fev <- predict(mod_fever, newdata = flu_16_17_1, type = "response")
flu_16_17_2$pred_fev <- predict(mod_fever, newdata = flu_16_17_2, type = "response")

flu_16_17_1$pred_fcst <- predict(mod_fever_cough_sore, newdata = flu_16_17_1, type = "response")
flu_16_17_2$pred_fcst <- predict(mod_fever_cough_sore, newdata = flu_16_17_2, type = "response")

# Recombine data
flu_16_17 <- rbind(flu_16_17_1, flu_16_17_2)
remove(flu_16_17_1, flu_16_17_2)

# Examine different cut-offs
flu_16_17$flu_combos <- ifelse(flu_16_17$pred_com > 0.50, 1, 0)
flu_16_17$flu_fever <- ifelse(flu_16_17$pred_fev >  0.50, 1, 0)
flu_16_17$flu_complex <- ifelse(flu_16_17$pred_fcst > 0.50, 1, 0)
#flu_16_17$flu_high_high <- ifelse(flu_16_17$pred > 0.85, 1, 0)





###
# count_dat <- flu_16_17 |> mutate(symp_count = fever + myalgia + cough + 
#                                    sore_throat + short_breath + hypoxemia + 
#                                    chest_pain + bronchitis + nausea_vom + 
#                                    diarrhea + fatigue + headache + congestion + sneezing,
#                                  no_symptoms = ifelse(symp_count == 0, patient_count_imp, 0),
#                                  fever = ifelse(fever == 1, patient_count_imp, 0),
#                                  myalgia = ifelse( myalgia == 1, patient_count_imp, 0),
#                                  short_breath = ifelse(short_breath == 1, patient_count_imp, 0),
#                                  hypoxemia = ifelse(hypoxemia == 1, patient_count_imp, 0),
#                                  chest_pain = ifelse(chest_pain == 1, patient_count_imp, 0),
#                                  cough = ifelse(cough == 1, patient_count_imp, 0),
#                                  bronchitis = ifelse(bronchitis == 1, patient_count_imp, 0),
#                                  nausea_vom = ifelse(nausea_vom == 1, patient_count_imp, 0),
#                                  sore_throat = ifelse(sore_throat == 1, patient_count_imp, 0),
#                                  diarrhea = ifelse(diarrhea == 1, patient_count_imp, 0),
#                                  fatigue = ifelse(fatigue == 1, patient_count_imp, 0),
#                                  headache = ifelse(headache == 1, patient_count_imp, 0),
#                                  congestion = ifelse(congestion == 1, patient_count_imp, 0),
#                                  sneezing = ifelse(sneezing == 1, patient_count_imp, 0))
# 
# by_flu <- count_dat |> group_by(flu) |> mutate(
#   fever_count_week = sum(as.numeric(fever)),
#   myalgia_count_week = sum(as.numeric(myalgia)),
#   short_breath_count_week = sum(as.numeric(short_breath)),
#   hypoxemia_count_week = sum(as.numeric(hypoxemia)),
#   chest_pain_count_week = sum(as.numeric(chest_pain)),
#   cough_count_week = sum(as.numeric(cough)),
#   bronchitis_count_week = sum(as.numeric(bronchitis)),
#   nausea_count_week = sum(as.numeric(nausea_vom)),
#   sore_throat_count_week = sum(as.numeric(sore_throat)),
#   diarrhea_count_week = sum(as.numeric(diarrhea)),
#   fatigue_count_week = sum(as.numeric(fatigue)),
#   headache_count_week = sum(as.numeric(headache)),
#   congestion_count_week = sum(as.numeric(congestion)),
#   sneezing_count_week = sum(as.numeric(sneezing)),
#   patient_count_week = sum(as.numeric(patient_count_imp)),
#   asymptomatic = sum(as.numeric(no_symptoms)),
#   symptom_count = mean(as.numeric(symp_count))) |>
#   distinct(flu, fever_count_week, myalgia_count_week, hypoxemia_count_week,
#            short_breath_count_week, cough_count_week, nausea_count_week, chest_pain_count_week,
#            sore_throat_count_week, diarrhea_count_week, fatigue_count_week, bronchitis_count_week,
#            headache_count_week, congestion_count_week, sneezing_count_week, patient_count_week, asymptomatic, symptom_count,
#            .keep_all = FALSE) |>
#   mutate(
#     fever_perc = fever_count_week / patient_count_week,
#     myalgia_perc = myalgia_count_week / patient_count_week,
#     short_breath_perc = short_breath_count_week / patient_count_week,
#     hypoxemia_perc = hypoxemia_count_week / patient_count_week,
#     chest_pain_perc = chest_pain_count_week / patient_count_week,
#     cough_perc = cough_count_week / patient_count_week,
#     bronchitis_perc = bronchitis_count_week / patient_count_week,
#     nausea_perc = nausea_count_week / patient_count_week,
#     sore_throat_perc = sore_throat_count_week / patient_count_week,
#     diarrhea_perc = diarrhea_count_week / patient_count_week,
#     fatigue_perc = fatigue_count_week / patient_count_week,
#     headache_perc = headache_count_week / patient_count_week,
#     congestion_perc = congestion_count_week / patient_count_week,
#     sneezing_perc = sneezing_count_week / patient_count_week,
#     asymptomatic_perc = asymptomatic / patient_count_week)
# 
# flu_16_17$all_symp <- ifelse(flu_16_17$fever == 1, 1, 0)
# flu_16_17$all_symp <- ifelse(flu_16_17$cough == 1, 1, flu_16_17$all_symp)
# flu_16_17$all_symp <- ifelse(flu_16_17$sore_throat == 1, 1, flu_16_17$all_symp)
# flu_16_17$all_symp <- ifelse(flu_16_17$hypoxemia == 1, 1, flu_16_17$all_symp)
# flu_16_17$all_symp <- ifelse(flu_16_17$congestion == 1, 1, flu_16_17$all_symp)
# flu_16_17$all_symp <- ifelse(flu_16_17$short_breath == 1, 1, flu_16_17$all_symp)
# flu_16_17$all_symp <- ifelse(flu_16_17$fatigue == 1, 1, flu_16_17$all_symp)
# flu_16_17$all_symp <- ifelse(flu_16_17$nausea_vom == 1, 1, flu_16_17$all_symp)
# flu_16_17$all_symp <- ifelse(flu_16_17$myalgia == 1, 1, flu_16_17$all_symp)
# 
# flu_symptoms <- flu_16_17 |> group_by(week_date) |>
#   mutate(flu = sum(flu*patient_count_imp),
#          fever = sum(fever*patient_count_imp),
#          myalgia = sum(myalgia*patient_count_imp),
#          short_breath = sum(short_breath*patient_count_imp),
#          hypoxemia = sum(hypoxemia*patient_count_imp),
#          chest_pain = sum(chest_pain*patient_count_imp),
#          cough = sum(cough*patient_count_imp),
#          bronchitis = sum(bronchitis*patient_count_imp),
#          nausea_vom = sum(nausea_vom*patient_count_imp),
#          sore_throat = sum(sore_throat*patient_count_imp),
#          diarrhea = sum(diarrhea*patient_count_imp),
#          fatigue = sum(fatigue*patient_count_imp),
#          headache = sum(headache*patient_count_imp),
#          congestion = sum(congestion*patient_count_imp),
#          sneezing = sum(sneezing*patient_count_imp)) |>
#   distinct(week_date, flu, fever, myalgia, cough, sore_throat, 
#            short_breath, hypoxemia, chest_pain, bronchitis,
#            nausea_vom, diarrhea, fatigue, headache,
#            congestion, sneezing, .keep_all = FALSE)
# 
# # Reshape data from wide to long format
# flu_symptoms_long <- melt(flu_symptoms, id.vars = c("week_date"), variable.name = "symptom", 
#                         value.name = "count")
# 
# ggplot() + 
#   geom_line(data = flu_symptoms_long, 
#             aes(x = week_date, y = count, color = symptom)) + 
#   theme_minimal()
# 
# 
# 
# 
# all_cause_week <- read.csv('raw/county_week_ac_v2.csv')
# # Impute <=5 to a random number between 1 and 5
# all_cause_week$patient_count_imp <- all_cause_week$all_cause
# all_cause_week$patient_count_imp[all_cause_week$patient_count_imp == '<=5'] <- 
#   sample(1:5, length(all_cause_week$patient_count_imp[all_cause_week$patient_count_imp == '<=5']), replace= TRUE)
# all_cause_week$patient_count_imp <- as.numeric(all_cause_week$patient_count_imp)
# 
# all_cause_week$week_date <- sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", all_cause_week$year_week)
# all_cause_week$week_date <- ISOweek2date(all_cause_week$week_date)
# 
# all_cause_week_county <- all_cause_week |> group_by(week_date) |>
#   mutate(all_cause_sum = sum(patient_count_imp)) |>
#   distinct(week_date, all_cause_sum) |>
#   filter(week_date > as.Date('2016-01-01')) |>
#   filter(week_date < as.Date('2023-12-24'))
# 
# ggplot() + 
#   geom_line(data = all_cause_week_county, 
#             aes(x = week_date, y = all_cause_sum)) + 
#   theme_minimal() + ylim(0, 16000000) + ggtitle('Claims over Time') + 
#   xlab('Time (Weekly)') + ylab('Number of Claims')
###

# Load all cause claims
all_cause <- read.csv('raw/county_season_ac_v2.csv')
all_cause <- all_cause |> filter(season == '2017-2018')
flu_16_17 <- left_join(flu_16_17, all_cause, by = c('county_fips' = 'county_fips'))

# Collapse to county
flu_16_17_county <- flu_16_17 |> group_by(week_date, county_fips) |> 
  mutate(flu_sum = sum(flu*patient_count_imp, na.rm = TRUE),
         ili_sum = sum(ili*patient_count_imp, na.rm = TRUE),
         flu_combos_sum = sum(flu_combos*patient_count_imp, na.rm = TRUE),
         flu_fever_sum = sum(pred_fev*patient_count_imp, na.rm = TRUE),
         flu_complex_sum = sum(flu_complex*patient_count_imp, na.rm = TRUE)) |>
  distinct(week_date, county_fips, flu_sum, ili_sum, flu_combos_sum, flu_fever_sum, flu_complex_sum, all_cause,
           .keep_all = FALSE)

#flu_16_17_county <- left_join(flu_16_17_county, all_cause, 
                              #by= c('week_date' = 'week_date',
                                    #'county_fips' = 'county_fips'))


# Collapse to country
flu_16_17_country <- flu_16_17_county |> group_by(week_date) |> 
  mutate(all_cause_sum = sum(all_cause),
         flu = sum(flu_sum) / all_cause_sum,
         ili = sum(ili_sum) / all_cause_sum,
         flu_combos = sum(flu_combos_sum) / all_cause_sum,
         flu_fever = sum(flu_fever_sum) / all_cause_sum,
         flu_complex = sum(flu_complex_sum) / all_cause_sum) |>
  distinct(week_date, flu, ili, flu_combos, flu_fever, flu_complex, all_cause_sum, .keep_all = FALSE)

flu_16_17_country$flu_scale <- (flu_16_17_country$flu - 
                                  mean(flu_16_17_country$flu)) / sd(flu_16_17_country$flu)

flu_16_17_country$flu_combos_scale <- (flu_16_17_country$flu_combos - 
                                  mean(flu_16_17_country$flu_combos)) / sd(flu_16_17_country$flu_combos)

flu_16_17_country$ili_scale <- (flu_16_17_country$ili - 
                                         mean(flu_16_17_country$ili)) / sd(flu_16_17_country$ili)

flu_16_17_country$flu_fever_scale <- (flu_16_17_country$flu_fever - 
                                  mean(flu_16_17_country$flu_fever)) / sd(flu_16_17_country$flu_fever)

flu_16_17_country$flu_complex_scale <- (flu_16_17_country$flu_complex - 
                                  mean(flu_16_17_country$flu_complex)) / sd(flu_16_17_country$flu_complex)

flu_16_17_country_scale <- flu_16_17_country |> 
  select(c(week_date, flu_scale, flu_combos_scale, flu_fever_scale, flu_complex_scale))

flu_16_17_long <- flu_16_17_country |> select(-c(all_cause_sum, flu_scale, ili_scale, ili, flu_combos_scale, flu_fever_scale, flu_complex_scale)) |>
  pivot_longer(!week_date, names_to = "flu_version", values_to = "count")

ggplot(flu_16_17_long) +
  geom_line(aes(x = week_date, y = as.numeric(count), color = flu_version)) +
  theme_minimal() + ylab('proportion') + theme(legend.position = 'bottom') 

flu_16_17_long_scale <- flu_16_17_country_scale |>
  pivot_longer(!week_date, names_to = "flu_version", values_to = "Z")

ggplot(flu_16_17_long_scale) +
  geom_line(aes(x = week_date, y = as.numeric(Z), color = flu_version)) +
  theme_minimal() + ylab('z') + theme(legend.position = 'bottom') 



flu_16_17_county_select <- flu_16_17_county |> 
  filter(county_fips == "17031" | county_fips == "53033" | 
           county_fips == "06037" | county_fips == "04013" |
           county_fips == "36061" | county_fips == "48113" |
           county_fips == "11001" | county_fips == "13121" |
           county_fips == "29510")

flu_16_17_long <- flu_16_17_county_select |> select(c(county_fips, flu_sum, flu_fever_sum, flu_complex_sum, flu_combos_sum, week_date)) |>
  pivot_longer(!c(week_date, county_fips), names_to = "flu_version", values_to = "count")


ggplot(flu_16_17_long) +
  geom_line(aes(x = week_date, y = as.numeric(count), color = flu_version)) +
  theme_minimal() + ylab('proportion') + theme(legend.position = 'bottom') +
  facet_wrap(vars(county_fips), scales = "free", nrow = 5) 


# cook, king, la, maricopa, nyc, dallas, dc, atl, stl
county_pop_groups <- readRDS('./tmp/county_group_pop.rds') 

flu_16_17_county <- left_join(flu_16_17_county, county_pop_groups, by = c('county_fips' = 'county_fips'))
flu_16_17_county_lim <- flu_16_17_county |> filter(county_group_pop > 100000)

flu_16_17_long_county <- flu_16_17_county_lim |> select(-c(all_cause, flu_combos_sum, county_group_pop)) |>
  pivot_longer(!c(week_date, county_fips), names_to = "flu_version", values_to = "count")



flu_16_17_long_county <- flu_16_17_long_county |> 
  filter(week_date >= as.Date('2017-09-01')) |>
  filter(week_date <= as.Date('2018-08-31')) |>
  group_by(county_fips, flu_version) |>
  mutate(roll_mean = rollmean(count, k = 4, align = 'right', fill = NA)) |>
  mutate(roll_mean_scale = (roll_mean - 
                              mean(roll_mean, na.rm = TRUE)) / sd(roll_mean, na.rm = TRUE)) |>
  mutate(max_date = week_date[which.max(roll_mean_scale)],
         start_date = week_date[which(roll_mean_scale > 0)[1]]) 

flu_16_17_long_county_dist <- flu_16_17_long_county |>
  distinct(county_fips, flu_version, max_date, start_date)
  
ggplot(flu_16_17_long_county_dist, aes(max_date, factor(flu_version))) +
  geom_violin()

ggplot(flu_16_17_long_county_dist, aes(start_date, factor(flu_version))) +
  geom_violin()

flu_16_17_long_county_dist_wide_start <- flu_16_17_long_county_dist |> 
  select(-c(max_date))

flu_16_17_long_county_dist_wide_start <- pivot_wider(flu_16_17_long_county_dist_wide_start, names_from = flu_version, values_from = start_date)

flu_16_17_long_county_dist_wide_start$ili_diff = 
  as.numeric(flu_16_17_long_county_dist_wide_start$flu_sum -
  flu_16_17_long_county_dist_wide_start$ili_sum)

flu_16_17_long_county_dist_wide_start$fever_diff = 
  as.numeric(flu_16_17_long_county_dist_wide_start$flu_sum -
  flu_16_17_long_county_dist_wide_start$flu_fever_sum)
  
flu_16_17_long_county_dist_wide_max <- flu_16_17_long_county_dist |> 
  select(-c(start_date))

flu_16_17_long_county_dist_wide_max <- pivot_wider(flu_16_17_long_county_dist_wide_max, names_from = flu_version, values_from = max_date)

flu_16_17_long_county_dist_wide_max$ili_diff = 
  as.numeric(flu_16_17_long_county_dist_wide_max$flu_sum -
               flu_16_17_long_county_dist_wide_max$ili_sum)

flu_16_17_long_county_dist_wide_max$fever_diff = 
  as.numeric(flu_16_17_long_county_dist_wide_max$flu_sum -
               flu_16_17_long_county_dist_wide_max$flu_fever_sum)



flu_16_17_long_county_dist_wide_max_long <- flu_16_17_long_county_dist_wide_max |>
  select(c(county_fips, ili_diff, fever_diff)) |>
  pivot_longer(!county_fips, names_to = "version", values_to = "difference")

ggplot(flu_16_17_long_county_dist_wide_max_long, aes(difference, version)) +
  geom_boxplot()

saveRDS(flu_16_17_long_county_dist_wide_max_long, './tmp/boxplot_max_17.rds') 

flu_16_17_long_county_dist_wide_start_long <- flu_16_17_long_county_dist_wide_start |>
  select(c(county_fips, ili_diff, fever_diff)) |>
  pivot_longer(!county_fips, names_to = "version", values_to = "difference")

ggplot(flu_16_17_long_county_dist_wide_start_long, aes(difference, version)) +
  geom_boxplot()

saveRDS(flu_16_17_long_county_dist_wide_start_long, './tmp/boxplot_start_17.rds') 

ggplot(flu_16_17_long_county_dist_wide_max, aes(flu_sum, flu_fever_sum)) +
  geom_point()



county_select_16 <- readRDS('./tmp/boxplot_start_16.rds')
county_select_17 <- readRDS('./tmp/boxplot_start_17.rds')
county_select_18 <- readRDS('./tmp/boxplot_start_18.rds')
county_select <- flu_16_17_long_county_dist

county_select$`Disease Estimate` <- county_select$flu_version
county_select$`Disease Estimate`[county_select$`Disease Estimate` == 'flu_fever_sum'] <- 'Fever Model'
county_select$`Disease Estimate`[county_select$`Disease Estimate` == 'ili_sum'] <- 'ILI Definition'
county_select$`Disease Estimate`[county_select$`Disease Estimate` == 'flu_sum'] <- 'Confirmed Influenza'
county_select <- county_select |> filter(flu_version != 'flu_complex_sum')

ggplot(county_select, aes(`Disease Estimate`, start_date, fill = `Disease Estimate`)) +
  geom_jitter(width = 0.25, alpha = 0.1) + geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.9) +
  theme_minimal() + theme(legend.position = 'none') + ylab('Season Onset Date') + xlab('Disease Estimate') +
  scale_fill_manual(values=c("#969696", "#FD8D3C", "#756BB1")) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12))


ggplot(county_select, aes(`Disease Estimate`, max_date, fill = `Disease Estimate`)) +
  geom_jitter(width = 0.25, alpha = 0.1) + geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.9) +
  theme_minimal() + theme(legend.position = 'none') + ylab('Season Peak Date') + xlab('Disease Estimate') +
  scale_fill_manual(values=c("#969696", "#FD8D3C", "#756BB1")) +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12))


county_select_16 <- readRDS('./tmp/boxplot_max_16.rds')
county_select_17 <- readRDS('./tmp/boxplot_max_17.rds')
county_select_18 <- readRDS('./tmp/boxplot_max_18.rds')
county_select <- rbind(county_select_16, county_select_17, county_select_18)

county_select$`Disease Estimate` <- county_select$version
county_select$`Disease Estimate`[county_select$`Disease Estimate` == 'fever_diff'] <- 'Fever Model'
county_select$`Disease Estimate`[county_select$`Disease Estimate` == 'ili_diff'] <- 'ILI Definition'

ggplot(county_select, aes(`Disease Estimate`, difference, fill = `Disease Estimate`), alpha = 1) +
  geom_jitter(width = 0.25, alpha = 0.02) + geom_boxplot(width = 0.5, outlier.shape = NA) +
  theme_minimal() + theme(legend.position = 'none') + ylab('Difference from Observed in Days') + xlab('Disease Estimate Difference') +
  scale_fill_manual(values=c("#00BA38", '#619CFF')) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        strip.text = element_text(size=12),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12))





#mean((as.numeric(flu_16_17_long_county_dist_wide$flu_sum - flu_16_17_long_county_dist_wide$flu_complex_sum)))
#mean((as.numeric(flu_16_17_long_county_dist_wide$flu_sum - flu_16_17_long_county_dist_wide$flu_fever_sum)))
#mean((as.numeric(flu_16_17_long_county_dist_wide$flu_sum - flu_16_17_long_county_dist_wide$ili_sum)))


#cor(flu_16_17_long_county_dist_wide$flu_sum, flu_16_17_long_county_dist_wide$ili_sum)



flu_16_17_county_select <- flu_16_17_long_county |> 
  filter(county_fips == "17031" | county_fips == "53033" | 
           county_fips == "06037" | #county_fips == "04013" |
           county_fips == "36061" | county_fips == "13121") |> #| #county_fips == "48113" |
           #county_fips == "11001" | county_fips == "13121" | county_fips == "29510" county_fips == "13121") 
           filter(flu_version != 'flu_complex_sum') |>
  group_by(county_fips) |>
  mutate(roll_count = rollmean(count, k = 4, align = 'right', fill = NA))


ggplot(flu_16_17_county_select) +
  geom_line(aes(x = week_date, y = as.numeric(roll_count), color = flu_version)) +
  theme_minimal() + ylab('proportion') + theme(legend.position = 'bottom') +
  facet_wrap(vars(county_fips), scales = "free", nrow = 1)

saveRDS(flu_16_17_county_select, './tmp/county_select_17.rds') 


county_select_16 <- readRDS('./tmp/county_select_16.rds')
county_select_17 <- readRDS('./tmp/county_select_17.rds')
county_select_18 <- readRDS('./tmp/county_select_18.rds')
county_select_19 <- readRDS('./tmp/county_select_19.rds')
county_select <- county_select_17

county_select$Name <- county_select$county_fips
county_select$Name[county_select$Name == '06037'] <- 'Los Angeles County, CA'
county_select$Name[county_select$Name == '13121'] <- 'Fulton County, GA'
county_select$Name[county_select$Name == '17031'] <- 'Cook County, IL'
county_select$Name[county_select$Name == '36061'] <- 'New York County, NY'
county_select$Name[county_select$Name == '53033'] <- 'King County, WA'


county_select$`Disease Estimate` <- county_select$flu_version
county_select$`Disease Estimate`[county_select$`Disease Estimate` == 'flu_fever_sum'] <- 'Fever Model'
county_select$`Disease Estimate`[county_select$`Disease Estimate` == 'flu_sum'] <- 'Confirmed Influenza'
county_select$`Disease Estimate`[county_select$`Disease Estimate` == 'ili_sum'] <- 'ILI Definition'

ggplot(county_select) +
  geom_line(aes(x = week_date, y = as.numeric(roll_mean_scale), color = `Disease Estimate`, linetype = `Disease Estimate`), size = 1, alpha = 0.9) +
  theme_minimal() + ylab('Z Score') + theme(legend.position = 'bottom') + xlab('Date') +
  facet_wrap(vars(Name), scales = "free", nrow = 5) +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        strip.text = element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14)) +
  scale_color_manual(values=c("dimgrey", "#FD8D3C", "#756BB1")) +
  scale_linetype_manual(values=c('dashed', 'solid', 'solid'))
        

county_select_filt <- county_select |> filter(`Disease Estimate` == 'ILI Definition') |>
  filter(Name == 'Los Angeles County, CA')
ggplot(county_select_filt) +
  geom_line(aes(x = week_date, y = as.numeric(roll_mean_scale), color = `Disease Estimate`, linetype = `Disease Estimate`), size = 2, alpha = 0.9) +
  theme_classic() + ylab('Z Score') + theme(legend.position = 'none') + xlab('Date') +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_blank(),
        legend.text=element_blank(),
        legend.title=element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_line(colour = 'black', size = 2)) +
  scale_color_manual(values=c("#756BB1")) +
  scale_linetype_manual(values=c('solid'))

                                   
library(RColorBrewer)
par(mar=c(3,4,2,2))
display.brewer.all()

brewer.pal(n=5,"Greys")

 "#F2F0F7" "#CBC9E2" "#9E9AC8" "#756BB1" "#54278F"
"#FEEDDE" "#FDBE85" "#FD8D3C" "#E6550D" "#A63603"
[1] "#F7F7F7" "#CCCCCC" "#969696" "#636363" "#252525"
[1] "#EFF3FF" "#BDD7E7" "#6BAED6" "#3182BD" "#08519C"

ggplot(flu_16_17_country) +
  geom_line(aes(x = week_date, y = as.numeric(flu_sum))) +
  geom_smooth(aes(x = week_date, y = as.numeric(flu_high_sum))) 

sum(flu_16_17_country$flu_high_sum)


chi_flu <- flu_16_17 |> filter(county_fips == '17031') |> group_by(week_date) |>
  mutate(flu_sum = sum(flu*patient_count_imp),
         flu_accuracy_sum = sum(flu_accuracy*patient_count_imp),
         flu_youden_sum = sum(flu_youden*patient_count_imp),
         fever_sum = sum(fever*patient_count_imp),
         sore_throat_sum = sum(sore_throat*patient_count_imp),
         cough_sum = sum(cough*patient_count_imp)) |>
  distinct(week_date, flu_sum, flu_accuracy_sum, flu_youden_sum, fever_sum, sore_throat_sum, cough_sum,
           .keep_all = FALSE)

flu_16_17_long_chi <- chi_flu |>
  pivot_longer(!week_date, names_to = "flu_version", values_to = "count")

ggplot(flu_16_17_long_chi) +
  geom_line(aes(x = week_date, y = as.numeric(count), color = flu_version)) +
  theme_minimal() + ylab('count') + theme(legend.position = 'bottom') 
summary(training_models[[1]])
