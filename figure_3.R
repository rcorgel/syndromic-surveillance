################################################################################
# File Name: figure_3                                                          #
#                                                                              #
# Purpose:   Create figure 3 for the manuscript.                               #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Create figure 3                                                #
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
library(ISOweek)
library(splus2R)
library(scales)
library(sf)
library(geosphere)
library(cowplot)
library(zoo)
library(gsignal)
library(mgcv)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

######################
# 2. CREATE FIGURE 3 #
######################

# Load other surveillance data
ili_nat <- readRDS('./tmp/ili_net_national.rds')
ili_state <- readRDS('./tmp/ili_net_state.rds')
flu_surv_nat <- readRDS('./tmp/flu_surv_net_national.rds')
flu_surv_state <- readRDS('./tmp/flu_surv_net_state.rds')
nrevss_nat <- readRDS('./tmp/nrevss_national.rds')
nrevss_state <- readRDS('./tmp/nrevss_state.rds')
nys_county <- readRDS('./tmp/nys_county.rds')

# Load primary model
time_model <- readRDS("./tmp/nrevss_model_state.rds")

# Load part 2 data
flu_2016 <- readRDS('./tmp/flu_filt_symp_p2_2016.rds')

#flu_2017 <- readRDS('./tmp/flu_filt_symp_p2_2017.rds')
#flu_2018 <- readRDS('./tmp/flu_filt_symp_p2_2018.rds')
#flu_2019 <- readRDS('./tmp/flu_filt_symp_p2_2019.rds')

#flu <- rbind(flu_2016, flu_2017, flu_2018)
#remove(flu_2016, flu_2017, flu_2018)

# Epi week
flu_2016$epi_week <- epiweek(flu_2016$year_week_dt)

# Change reference week to week 26
flu_2016$epi_week_adj <- ((flu_2016$epi_week - 26) %% 52) + 1

flu_2016$year_week_dt <- as.Date(flu_2016$year_week_dt)

# Create empty lists to fill
flu_nat_filt_county_week <- list()

# Set count
count <- 1

nrevss_nat <- readRDS('./tmp/nrevss_national.rds')
nrevss_state <- readRDS('./tmp/nrevss_state.rds')
flu_2016 <- left_join(flu_2016, nrevss_nat[, c(2, 3)], by = c('year_week_dt' = 'week_date'))
flu_2016 <- flu_2016 |> rename('nrevss_nat' = 'value')
flu_2016 <- left_join(flu_2016, nrevss_state[, c(2, 3, 4)], 
                 by = c('year_week_dt' = 'week_date',
                        'state_fips' = 'state_fips'))
flu_2016 <- flu_2016 |> rename('nrevss_state' = 'value')


for (i in unique(flu_2016$year_week_dt)) {
  # Print the progress and set season data
  print(i)
  flu_week <- flu_2016 |> dplyr::filter(year_week_dt == i)
  
  flu_week$pred <- predict(time_model, 
                          newdata = flu_week, 
                          type = "response",
                          allow.new.levels = TRUE)
  
  # Fill lists
  flu_nat_filt_county_week[[count]] <- as.data.frame(flu_week)
  remove(flu_week)
  
  # Update count
  count <- count + 1
}

remove(flu_2016)
remove(time_model)

#saveRDS(flu_nat_filt_county_week, './tmp/flu_pred_dat_county_p2_time_geo_week_list.rds')

#flu <- readRDS('./tmp/flu_pred_dat_county_p2_time_geo_week_list.rds')

# Append all week data frames
flu_pred_dat <- do.call(rbind, flu_nat_filt_county_week)

# Save
saveRDS(flu_pred_dat, './tmp/flu_pred_dat_county_p2_time_dis_week_nrevss_state.rds')

# Remove data frame list
remove(flu_nat_filt_county_week)









#flu_pred_dat <- readRDS('./tmp/flu_pred_dat_county_p2_time_geo_week.rds')



library(mgcv)


test_2019 <- flu_2019 |> group_by(year_week_dt, county_fips, state_fips) |>
  mutate(flu = sum(flu * patient_count),
         fever = sum(fever * patient_count),
         myalgia = sum(myalgia * patient_count),
         cough = sum(cough * patient_count),
         sore_throat = sum(sore_throat * patient_count), 
         short_breath = sum(short_breath * patient_count), 
         hypoxemia = sum(hypoxemia * patient_count), 
         chest_pain = sum(chest_pain * patient_count), 
         bronchitis = sum(bronchitis * patient_count),
         nausea_vom = sum(nausea_vom * patient_count), 
         diarrhea = sum(diarrhea * patient_count), 
         fatigue = sum(fatigue * patient_count), 
         headache = sum(headache * patient_count),
         congestion = sum(congestion * patient_count), 
         sneezing = sum(sneezing * patient_count)) |>
  distinct(year_week_dt, county_fips, state_fips, flu, fever, myalgia, cough, sore_throat, 
           short_breath, hypoxemia, chest_pain, bronchitis,
           nausea_vom, diarrhea, fatigue, headache,
           congestion, sneezing)


testing_nyc <- rbind(test_2017[test_2017$county_fips == '36061',], test_2018[test_2018$county_fips == '36061',], test_2019[test_2019$county_fips == '36061',])
testing_tompkins <- rbind(test_2017[test_2017$county_fips == '36109',], test_2018[test_2018$county_fips == '36109',], test_2019[test_2019$county_fips == '36109',])

testing <- rbind(test_2017, test_2018)

test <- left_join(test, all_cause_season, by = c('county_fips'))
test_test <- test |> dplyr::filter(all_cause_wtd > 50000)

# Get unique counties
counties <- unique(test_test$county_fips)

# Create empty list to store models
county_models <- list()

# Fit a model for each county
for (county_name in counties) {
  cat("Fitting model for:", county_name, "\n")
  
  # Subset data for this county
  county_data <- test_test[test_test$county_fips == county_name, ]
  
  # Fit GAM for this county
  model <- gam(
    flu ~ 
      s(fever, k = 3) +
      s(cough, k = 3) +
      s(sore_throat, k = 3) +
      s(congestion, k = 3) +
      s(bronchitis, k = 3) +
      s(chest_pain, k = 3) +
      s(nausea_vom, k = 3) +
      s(fatigue, k = 3) +
      s(diarrhea, k = 3) +
      s(sneezing, k = 3) +
      s(hypoxemia, k = 3) +
      s(short_breath, k = 3) +
      s(myalgia, k = 3) +
      s(headache, k = 3),
    data = county_data,
    family = quasipoisson(link = "log"),
    method = "REML"
  )
  
  # Store model in list with county name
  county_models[[county_name]] <- model
}






# Epi week
testing_tompkins$epi_week <- epiweek(testing_tompkins$year_week_dt)

# Change reference week to week 26
testing_tompkins$epi_week_adj <- ((testing_tompkins $epi_week - 26) %% 52) + 1

# Epi week
test$epi_week <- epiweek(test$year_week_dt)

# Change reference week to week 26
test$epi_week_adj <- ((test$epi_week - 26) %% 52) + 1
test$state_fips <- as.factor(test$state_fips)
test$county_fips <- as.factor(test$county_fips)

test_counties <- test |> 
  dplyr::filter(county_fips %in% c('17031', '11001', '12086', '36061', '36109', '55025', '53063'))

all_cause_season <- read_csv('raw/county_season_ac_v3_imputed.csv')
all_cause_season <- all_cause_season |> dplyr::filter(season == '2016-2017')

test_counties <- left_join(test_counties, all_cause_season, by = c('county_fips'))

test_counties <- test_counties |>
  mutate(flu = flu / all_cause_wtd,
         fever = fever / all_cause_wtd,
         myalgia = myalgia / all_cause_wtd,
         cough = cough / all_cause_wtd,
         sore_throat = sore_throat / all_cause_wtd, 
         short_breath = short_breath / all_cause_wtd, 
         hypoxemia = hypoxemia / all_cause_wtd, 
         chest_pain = chest_pain / all_cause_wtd, 
         bronchitis = bronchitis / all_cause_wtd,
         nausea_vom = nausea_vom / all_cause_wtd, 
         diarrhea = diarrhea / all_cause_wtd, 
         fatigue = fatigue / all_cause_wtd, 
         headache = headache / all_cause_wtd,
         congestion = congestion / all_cause_wtd, 
         sneezing = sneezing / all_cause_wtd)

test_counties <- test_counties |>
  group_by(county_fips) |>
  mutate(fev_1 = lag(fever, n = 1),
         fev_2 = lag(fever, n = 2),
         fev_3 = lag(fever, n = 3),
         fev_4 = lag(fever, n = 4))

testing_tompkins <- testing_tompkins |>
  group_by(county_fips) |>
  mutate(fev_1 = lag(fever, n = 1),
         fev_2 = lag(fever, n = 2),
         fev_3 = lag(fever, n = 3),
         fev_4 = lag(fever, n = 4))

model <- gam(flu ~ 
               s(fever, k = 3) +
               s(cough, k = 3) +
               s(sore_throat, k = 3) +
               s(congestion, k = 3) +
               s(bronchitis, k = 3) +
               s(chest_pain, k = 3) +
               s(nausea_vom, k = 3) +
               s(fatigue, k = 3) +
               s(diarrhea, k = 3) +
               s(sneezing, k = 3) +
               s(hypoxemia, k = 3) +
               s(short_breath, k = 3) +
               s(myalgia, k = 3) +
               s(headache, k = 3),
             data = test[test$county_fips == '01005 ',], 
             family = quasipoisson(link = "log"),
             method = "REML")

testing_nyc$pred <- predict(model,
                              newdata = testing_nyc, type = "response") 

ggplot(testing_nyc) + 
  geom_line(aes(x = as.Date(year_week_dt), y = pred), color = 'blue') +
  geom_line(aes(x = as.Date(year_week_dt), y = flu), color = 'red') +
  ylim(0, 5000)

test$county_fips <- as.factor(test$county_fips)

model <- gam(flu ~ 
               s(fever, k = 3) +
               s(fever, county_fips, bs = 're') +
               s(cough, k = 3) +
               s(cough, county_fips, bs = 're') +
               s(sore_throat, k = 3) +
               s(sore_throat, county_fips, bs = 're') +
               s(congestion, k = 3) +
               s(congestion, county_fips, bs = 're') +
               s(bronchitis, k = 3) +
               s(bronchitis, county_fips, bs = 're') +
               s(chest_pain, k = 3) +
               s(chest_pain, county_fips, bs = 're') +
               s(nausea_vom, k = 3) +
               s(nausea_vom, county_fips, bs = 're') +
               s(fatigue, k = 3) +
               s(fatigue, county_fips, bs = 're') +
               s(diarrhea, k = 3) +
               s(diarrhea, county_fips, bs = 're') +
               s(sneezing, k = 3) +
               s(sneezing, county_fips, bs = 're') +
               s(hypoxemia, k = 3) +
               s(hypoxemia, county_fips, bs = 're') +
               s(short_breath, k = 3) +
               s(short_breath, county_fips, bs = 're') +
               s(myalgia, k = 3) +
               s(myalgia, county_fips, bs = 're') +
               s(headache, k = 3) +
               s(headache, county_fips, bs = 're') +
               s(county_fips, bs = 're') +
               s(epi_week_adj, bs = "cc", k = 20),
    data = test, 
    family = quasipoisson(link = "log"),
    method = "REML")

summary(model)

test_counties$pred <- predict(model,
        newdata = test_counties, type = "response") 

ggplot(test_counties[test_counties$county_fips == '36109',]) + 
  geom_line(aes(x = as.Date(year_week_dt), y = pred), color = 'blue') +
  geom_line(aes(x = as.Date(year_week_dt), y = flu), color = 'red')


# Collapse data to the national scale
flu_national <- flu_pred_dat |> group_by(year_week_dt) |>
  mutate(pred_sum = sum(pred_count, na.rm = T),
         flu_sum = sum(flu * patient_count)) |>
  distinct(year_week_dt, pred_sum, flu_sum) |>
  dplyr::filter(year_week_dt < as.Date('2020-03-01')) |>
  dplyr::filter(year_week_dt > as.Date('2016-08-31')) |>
  ungroup()

# Quickly plot results
ggplot(flu_national) + 
  geom_line(aes(x = as.Date(year_week_dt), y = pred_sum), color = 'blue') +
  geom_line(aes(x = as.Date(year_week_dt), y = flu_sum), color = 'red')

# Quickly plot results
ggplot(flu_county[flu_county$county_fips == '36109',]) + 
  geom_line(aes(x = as.Date(year_week_dt), y = pred_sum), color = 'blue') +
  geom_line(aes(x = as.Date(year_week_dt), y = flu_sum), color = 'red')


# Set cutoff value for prediction
flu_pred_dat$pred_count <- ifelse(flu_pred_dat$pred > 0.50, 1, 0) * flu_pred_dat$patient_count







# COUNTY LEVEL
flu_county <- flu_pred_dat |> 
  group_by(year_week_dt, county_fips) |>
  mutate(pred_sum = sum(pred_count, na.rm = T),
         flu_sum = sum(patient_count * flu)) |>
  distinct(year_week_dt, county_fips, pred_sum, flu_sum)

#flu_county <- flu_county |> rename('year_week_dt' = 'week_date')
# All cause county
# Load all cause data
all_cause <- read.csv('./raw/county_week_ac_v2.csv', header = T)

# Impute all cause counts
all_cause$patient_count_imp <- all_cause$all_cause
all_cause$patient_count_imp[all_cause$patient_count_imp == '<=5'] <- 
  sample(1:5, length(all_cause$patient_count_imp[all_cause$patient_count_imp == '<=5']), 
         replace= TRUE)
# Convert imputed count to numeric
all_cause$patient_count_imp <- as.numeric(all_cause$patient_count_imp)

# Convert year week to week date
all_cause <- all_cause |>
  mutate(week_date = sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", year_week)) |>
  mutate(week_date = as.Date(ISOweek2date(week_date)))

# Collapse to county level
all_cause_county <- all_cause |> group_by(week_date, county_fips) |>
  mutate(all_cause_sum = sum(patient_count_imp)) |>
  distinct(week_date, county_fips, all_cause_sum) |>
  ungroup() |>
  dplyr::filter(week_date < as.Date('2019-09-01')) |>
  dplyr::filter(week_date > as.Date('2016-08-31')) |>
  group_by(county_fips) |>
  mutate(all_cause_avg = mean(all_cause_sum))

# Merge on all cause county
flu_county$year_week_dt <- as.Date(flu_county$year_week_dt)
flu_county <- left_join(flu_county, all_cause_county, by = c('county_fips','year_week_dt' = 'week_date'))

# Adjust count
flu_county$flu_obs_scale <- flu_county$flu_sum *
  (1/(flu_county$all_cause_sum / flu_county$all_cause_avg))
flu_county$flu_pred_scale <- flu_county$pred_sum *
  (1/(flu_county$all_cause_sum / flu_county$all_cause_avg))

# Add year and month variables
flu_county$year <- year(flu_county$year_week_dt)
flu_county$month <- month(flu_county$year_week_dt)

# Calculate summer means
flu_county_summer <- flu_county |> 
  dplyr::filter(month < 8 & month > 5) |> 
  group_by(year, county_fips) |>
  dplyr::filter(year > 2016) |>
  mutate(pred_summ_mean = mean(pred_sum),
         obs_summ_mean = mean(flu_sum)) |>
  distinct(year, county_fips, pred_summ_mean, 
           obs_summ_mean)

flu_county_summer$count <- 1
flu_county_summer[flu_county_summer$year == 2017,]$count <- 2
flu_county_summer <- flu_county_summer |> uncount(count)
flu_county_summer <- flu_county_summer |> group_by(county_fips) |>
  mutate(index = row_number(),
         year = ifelse(index == 1, 2016, year))

# Merge on summer data to the correct season (summer 2016 should align with the 2016-2017 flu season)
flu_county$year_2 <- ifelse(flu_county$month < 6, flu_county$year - 1, flu_county$year)
flu_county <- left_join(flu_county, flu_county_summer[, c(1, 2, 3, 4)], by = c('year_2' = 'year', 'county_fips' = 'county_fips'))

# Calculate Z score
flu_county <- flu_county |> group_by(county_fips) |>
  mutate(pred_z = (pred_sum - pred_summ_mean) / sd(pred_sum, na.rm = T),
         obs_z = (flu_sum - obs_summ_mean) / sd(flu_sum, na.rm = T))


# CALCULATE PEAKS #
flu_county_peak <- flu_county





county_acf <- flu_county |> group_by(county_fips) |>
  mutate(case_sum = sum(flu_sum)) |>
  dplyr::filter(case_sum > 0) |>
  mutate(max_cor = max(ccf(flu_sum, pred_sum)[["acf"]]),
         lag = ccf(flu_sum, pred_sum)[["lag"]][which.max(ccf(flu_sum, pred_sum)[["acf"]])]) |>
  distinct(county_fips, max_cor, lag)


ggplot(data = county_acf) + geom_density(aes(max_cor))
ggplot(data = county_acf) + geom_density(aes(lag))

usa_albers_county <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_county.rds'))   # convert to sf

usa_albers_county <- left_join(usa_albers_county, county_acf, 
                               by = c('GEOID' = 'county_fips'))

map_theme <- theme(legend.position = "right",
                   axis.title = element_text(size=16),
                   strip.background = element_blank(),
                   plot.title = element_blank(),
                   panel.border = element_rect(fill=NA, linewidth = 0.6, color = '#FFFFFF00'),
                   legend.title = element_text(size = 16),
                   legend.text = element_text(size = 12),
                   legend.key.height = unit(0.7, 'cm'),
                   legend.key.width = unit(0.4, "cm"))

cont_map <- ggplot() +
  geom_sf(data = usa_albers_county[abs(usa_albers_county$lag) < 3,]  , aes(fill = lag), , color= 'black', linewidth = 0.2) +
  theme_void() + map_theme 

county_z <- flu_county |> group_by(county_fips, year) |>
  mutate(max = max(obs_z),
         max_pred = max(pred_z)) |>
  distinct(county_fips, year, max, max_pred)

length(county_z[county_z$max > 1,]$max) 

ggplot(flu_county) + geom_line(aes(x = as.Date(year_week_dt), y = pred_z, group = county_fips))


# Create Season Variable
flu_county_peak$Season <- ifelse(flu_county_peak$year == 2016, '2016-2017', '')
flu_county_peak$Season <- ifelse(flu_county_peak$year == 2017 & flu_county_peak$month < 9, '2016-2017', flu_county_peak$Season)
flu_county_peak$Season <- ifelse(flu_county_peak$year == 2017 & flu_county_peak$month > 8, '2017-2018', flu_county_peak$Season)
flu_county_peak$Season <- ifelse(flu_county_peak$year == 2018 & flu_county_peak$month < 9, '2017-2018', flu_county_peak$Season)
flu_county_peak$Season <- ifelse(flu_county_peak$year == 2018 & flu_county_peak$month > 8, '2018-2019', flu_county_peak$Season)
flu_county_peak$Season <- ifelse(flu_county_peak$year == 2019 & flu_county_peak$month < 9, '2018-2019', flu_county_peak$Season)
flu_county_peak$Season <- ifelse(flu_county_peak$year == 2019 & flu_county_peak$month > 8, '2019-2020', flu_county_peak$Season)
flu_county_peak$Season <- ifelse(flu_county_peak$year == 2020, '2019-2020', flu_county_peak$Season)

library(zoo)
flu_county_peak <- flu_county_peak |> arrange(year_week_dt) |> group_by(county_fips) |> 
  mutate(flu_sum_roll = rollmean(flu_sum, k = 4, align = 'right', na.pad = TRUE),
         flu_pred_roll = rollmean(pred_sum, k = 4, align = 'right', na.pad = TRUE))

ggplot(flu_county_peak[flu_county_peak$county_fips == '36109',]) + 
  geom_line(aes(x = year_week_dt, y = flu_sum), color = 'blue') +
  geom_line(aes(x = year_week_dt, y = pred_sum), color = 'red')

# Load filtering data
all_cause_season <- read_csv('raw/county_season_ac_v2.csv')
county_group_pop <- readRDS('tmp/county_group_pop.rds')

# Merge on filtering data
flu_county_peak <- left_join(flu_county_peak, all_cause_season, by = c('Season' = 'season', 'county_fips'))
flu_county_peak <- left_join(flu_county_peak, county_group_pop, by = c('county_fips'))

# Limit data based on county pop and number of claims in a season
flu_county_peak_filt <- flu_county_peak |>
  dplyr::filter(county_group_pop > 10000) |>
  dplyr::filter(all_cause > 5000)

# Calculate peaks
county_peaks_usa <- flu_county_peak_filt |> 
  ungroup() |>
  dplyr::filter(!is.na(flu_sum_roll)) |>
  dplyr::filter(!is.na(flu_pred_roll)) |>
  # Calculate the row number
  mutate(index = row_number()) |>
  group_by(Season, county_fips) |> 
  # Calculate multiple peaks greater than 10 cases and larger than 10 neighbors
  mutate(peaks_obs = seq_along(index) %in% 
           findpeaks(flu_sum_roll, MinPeakDistance = 10, MinPeakHeight = 20)$loc) |>
  mutate(peaks_pred = seq_along(index) %in% 
           findpeaks(flu_pred_roll, MinPeakDistance = 10, MinPeakHeight = 20)$loc) |>
  # Remove first and last observations, which cannot be peaks but might be max
  dplyr::filter(row_number() != 1 & row_number() != n()) |>
  # Calculate the maximum value among the peaks
  mutate(is_max_obs = flu_sum_roll == max(flu_sum_roll),
         is_max_pred = flu_pred_roll == max(flu_pred_roll)) |>
  select(c(Season, county_fips, year_week_dt, flu_pred_roll, flu_sum_roll, obs_z, pred_z,
           peaks_obs, peaks_pred, is_max_obs, is_max_pred))

county_peaks_obs <- county_peaks_usa |> 
  # Select the peaks
  dplyr::filter(peaks_obs == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_obs == T) |>
  select(c(Season, county_fips, year_week_dt, peaks_obs, is_max_obs, obs_z)) |>
  rename('week_obs' = 'year_week_dt') |>
  group_by(Season, county_fips) |>
  # Break any ties by taking the earliest peak
  arrange(week_obs) |>
  slice_head(n = 1) |>
  mutate(week_obs_num = week(week_obs)) |>
  # Filter for significant peaks w/ z score > 1
  dplyr::filter(obs_z > 1)

county_peaks_pred <- county_peaks_usa |> 
  # Select the peaks
  dplyr::filter(peaks_pred == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_pred == T) |>
  select(c(Season, county_fips, year_week_dt, peaks_pred, is_max_pred, pred_z)) |>
  rename('week_pred' = 'year_week_dt') |>
  group_by(Season, county_fips) |>
  # Break any ties by taking the earliest peak
  arrange(week_pred) |>
  slice_head(n = 1) |>
  mutate(week_pred_num = week(week_pred)) |>
  # Filter for significant peaks w/ z score > 1
  dplyr::filter(pred_z > 1)

# Merge obs and pred peaks
county_peaks_both <- left_join(county_peaks_pred, county_peaks_obs,
                               by = c('Season' = 'Season',
                                      'county_fips' = 'county_fips'))

sd(county_peaks_obs$week_obs_num)

# Change ordering of Epi weeks for calculations and visualizations
county_peaks_both_filt <- county_peaks_both |>
  mutate(week_obs_num_adj = ifelse(week_obs_num > 30, week_obs_num - 52, week_obs_num),
         week_pred_num_adj = ifelse(week_pred_num > 30, week_pred_num - 52, week_pred_num))

# Calculate statistics
county_peaks_both_filt$diff <- as.numeric(county_peaks_both_filt$week_pred - county_peaks_both_filt$week_obs)
rmse <- sqrt(sum((county_peaks_both_filt$diff)^2, na.rm = T) / length(county_peaks_both_filt$diff)) 
rmse
cor(county_peaks_both_filt$week_pred_num, county_peaks_both_filt$week_obs_num, 
    use = 'complete.obs', method = 'spearman')

# Collapse data for graphing purposes
county_peaks_both_col <- county_peaks_both_filt |>
  mutate(count = 1) |>
  group_by(week_obs_num, week_pred_num) |>
  mutate(count_sum = sum(count)) |>
  distinct(week_obs_num, week_pred_num, week_obs_num_adj, week_pred_num_adj, count_sum)

# Plot
ggplot(county_peaks_both_col, aes(x = week_obs_num, y = week_pred_num)) +
  geom_point(aes(size = count_sum), color = "#e6a532", alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  labs(title = 'Syndromic Prediction vs. Claims\nPeak Timing',
       x = "Obs Peak Week",
       y = "Syndromic Peak Week") + xlim(-10, 20) + ylim(-10, 20) +
  scale_size_continuous(range = c(1, 10), breaks = c(1, 10, 100, 200, 400)) +
  theme_minimal()


library(viridis)
peak_comp_plot <- ggplot(county_peaks_both_col, aes(week_obs_num_adj, week_pred_num_adj, fill= count_sum)) + 
  geom_tile() + xlim(-5, 22) + ylim(-5, 22) +
  scale_fill_viridis(discrete=FALSE, direction = -1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black",
              linewidth = 1) +
  theme_minimal() + 
  labs(title = 'Peak Timing Comparison',
       x = "Claims Peak Epidemiological Week",
       y = "Syndromic Peak Epidemiological Week")


peak_comp_plot 

usa_albers_county <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_county.rds'))   # convert to sf

usa_albers_county_2018 <- left_join(usa_albers_county, 
                               county_peaks_both_filt, 
                               by = c('GEOID' = 'county_fips'))

usa_albers_county_2018$diff <- usa_albers_county_2018$week_obs - usa_albers_county_2018$week_pred

ggplot() +
  geom_sf(data = usa_albers_county_2018, aes(fill = as.numeric(diff)), color= 'black', linewidth = 0.1) +
  scale_fill_gradient2(low = 'red', high = 'blue', midpoint = 0, mid = 'white') +
  #scale_fill_viridis(discrete=FALSE, direction = -1) +
  ggtitle('2018-2019, Difference') +
  theme_void()

ggplot() +
  geom_density(data = usa_albers_county_2018, aes(x = as.numeric(diff))) 

flu_surv_state$year <- year(flu_surv_state$week_date)
flu_surv_state$month <- month(flu_surv_state$week_date)


# Create Season Variable
flu_surv_state$Season <- ifelse(flu_surv_state$year == 2016, '2016-2017', '')
flu_surv_state$Season <- ifelse(flu_surv_state$year == 2017 & flu_surv_state$month < 9, '2016-2017', flu_surv_state$Season)
flu_surv_state$Season <- ifelse(flu_surv_state$year == 2017 & flu_surv_state$month > 8, '2017-2018', flu_surv_state$Season)
flu_surv_state$Season <- ifelse(flu_surv_state$year == 2018 & flu_surv_state$month < 9, '2017-2018', flu_surv_state$Season)
flu_surv_state$Season <- ifelse(flu_surv_state$year == 2018 & flu_surv_state$month > 8, '2018-2019', flu_surv_state$Season)
flu_surv_state$Season <- ifelse(flu_surv_state$year == 2019 & flu_surv_state$month < 9, '2018-2019', flu_surv_state$Season)
flu_surv_state$Season <- ifelse(flu_surv_state$year == 2019 & flu_surv_state$month > 8, '2019-2020', flu_surv_state$Season)
flu_surv_state$Season <- ifelse(flu_surv_state$year == 2020, '2019-2020', flu_surv_state$Season)

library(zoo)
flu_surv_state_peak <- flu_surv_state |> arrange(week_date) |> group_by(state_fips) |> 
  mutate(flu_sum_roll = rollmean(value, k = 4, align = 'right', na.pad = TRUE))


# Calculate peaks
flu_surv_state_peaks <- flu_surv_state_peak |> 
  ungroup() |>
  dplyr::filter(!is.na(value_roll)) |>
  # Calculate the row number
  mutate(index = row_number()) |>
  group_by(Season, region) |> 
  # Calculate multiple peaks greater than 10 cases and larger than 10 neighbors
  mutate(peaks_obs = seq_along(index) %in% 
           findpeaks(value_roll, MinPeakDistance = 10, MinPeakHeight = 1)$loc) |>
  # Remove first and last observations, which cannot be peaks but might be max
  dplyr::filter(row_number() != 1 & row_number() != n()) |>
  # Calculate the maximum value among the peaks
  mutate(is_max_obs = flu_sum_roll == max(flu_sum_roll)) |>
  select(c(Season, region, state_fips, week_date, value_roll,
           peaks_obs,  is_max_obs))

county_peaks_obs <- county_peaks_usa |> 
  # Select the peaks
  dplyr::filter(peaks_obs == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_obs == T) |>
  select(c(Season, county_fips, year_week_dt, peaks_obs, is_max_obs, obs_z)) |>
  rename('week_obs' = 'year_week_dt') |>
  group_by(Season, county_fips) |>
  # Break any ties by taking the earliest peak
  arrange(week_obs) |>
  slice_head(n = 1) |>
  mutate(week_obs_num = week(week_obs)) |>
  # Filter for significant peaks w/ z score > 1
  dplyr::filter(obs_z > 1)








nyc <- flu_county_peak |> 
  ungroup() |>
  mutate(index = row_number()) |>
  group_by(Season, county_fips) |>
  mutate(peaks = seq_along(index) %in% 
           findpeaks(pred_sum, MinPeakDistance = 5, MinPeakHeight = 5)$loc)
  



flu_county_peak_test <- flu_county_peak |>
  dplyr::filter(Season == '2016-2017' & county_fips == '05059')


flag_peaks(flu_county_peak_test$pred_sum, mpd = 10, mph = 5)


nyc <- flu_county_peak[flu_county_peak$Season == '2016-2017' & flu_county_peak$county_fips == '05059', ]

ggplot(nyc) + 
  geom_line(aes(y = pred_sum, x = week_date)) 

res <- findpeaks(nyc$pred_sum, MinPeakDistance = 52, MinPeakHeight = 5)


nyc <- flu_county_peak |> 
  ungroup() |>
  mutate(index = row_number()) |>
  group_by(Season, county_fips) |>
  mutate(peaks = seq_along(index) %in% 
           findpeaks(pred_sum, MinPeakDistance = 48, MinPeakHeight = 5)$loc)


nyc_peaks <- nyc |> dplyr::filter(peaks == TRUE) 
  
ggplot(nyc) + 
  geom_line(aes(y = flu_sum, x = week_date)) +
  geom_point(data = nyc_peaks, aes(y = flu_sum, x = week_date))




# STATE LEVEL
flu_state <- flu_pred_dat |> 
  mutate(state_fips = substr(county_fips, 1, 2)) |>
  group_by(week_date, state_fips) |>
  mutate(pred_sum = pred_count),
         flu_sum = patient_sum * flu)) |>
  distinct(week_date, pred_sum, flu_sum) |>
  filter(week_date < as.Date('2020-03-01')) |>
  filter(week_date > as.Date('2016-08-31'))


# Merge on all cause county
flu_county <- left_join(flu_county, all_cause_county, by = c('county_fips', 'week_date'))

flu_county_peak <- flu_county

flu_county_peak$month <- month(flu_county_peak$week_date)
flu_county_peak$year <- year(flu_county_peak$week_date)


# Create Season Variable
flu_county_peak$Season <- ifelse(flu_county_peak$year == 2016, '2016-2017', '')
flu_county_peak$Season <- ifelse(flu_county_peak$year == 2017 & flu_county_peak$month < 9, '2016-2017', flu_county_peak$Season)
flu_county_peak$Season <- ifelse(flu_county_peak$year == 2017 & flu_county_peak$month > 8, '2017-2018', flu_county_peak$Season)
flu_county_peak$Season <- ifelse(flu_county_peak$year == 2018 & flu_county_peak$month < 9, '2017-2018', flu_county_peak$Season)
flu_county_peak$Season <- ifelse(flu_county_peak$year == 2018 & flu_county_peak$month > 8, '2018-2019', flu_county_peak$Season)
flu_county_peak$Season <- ifelse(flu_county_peak$year == 2019 & flu_county_peak$month < 9, '2018-2019', flu_county_peak$Season)
flu_county_peak$Season <- ifelse(flu_county_peak$year == 2019 & flu_county_peak$month > 8, '2019-2020', flu_county_peak$Season)
flu_county_peak$Season <- ifelse(flu_county_peak$year == 2020, '2019-2020', flu_county_peak$Season)

# Calculate peaks
county_peaks_usa <- flu_county_peak |> filter(Season != '2019-2020') |>
  group_by(Season, county_fips) |>
  mutate(pred_peak = peaks(pred_sum, span=101, strict=FALSE, endbehavior = 1)) |>
  mutate(obs_peak = peaks(flu_sum, span=101, strict=FALSE, endbehavior = 1)) |>
  select(c(Season, county_fips, week_date, pred_sum, flu_sum,
           pred_peak, obs_peak, urban_code_binary, state_fips))


county_peaks_pred <- county_peaks_usa |> filter(pred_peak == T) |>
  select(c(Season, county_fips, week_date, pred_peak, urban_code_binary, state_fips)) |>
  rename('county_fips' = 'county_fips',
         'week_pred' = 'week_date') |>
  group_by(Season, county_fips) |>
  arrange(week_pred) |>
  slice_head(n = 1) |>
  mutate(week_pred_num = week(week_pred))

county_peaks_obs <- county_peaks_usa |> filter(obs_peak == T) |>
  select(c(Season, county_fips, week_date, obs_peak)) |>
  rename('county_fips' = 'county_fips',
         'week_obs' = 'week_date') |>
  group_by(Season, county_fips) |>
  arrange(week_obs) |>
  slice_head(n = 1) |>
  mutate(week_obs_num = week(week_obs))


county_peaks_both <- left_join(county_peaks_pred, county_peaks_obs,
                               by = c('Season' = 'Season',
                                      'county_fips' = 'county_fips'))

all_cause_season <- read_csv('raw/county_season_ac_v2.csv')
county_group <- readRDS('tmp/county_group_pop.rds')


county_peaks_both <- left_join(county_peaks_both, all_cause_season, by = c('Season', 'county_fips'))
county_peaks_both <- left_join(county_peaks_both, county_group, by = c('county_fips'))

county_peaks_both_filt <- county_peaks_both |>
  filter(county_group_pop > 10000) |>
  filter(all_cause > 10000)

county_peaks_both_col <- county_peaks_both |>
  filter(county_group_pop > 10000) |>
  filter(all_cause > 10000) |>
  mutate(count = 1) |>
  group_by(week_obs_num, week_pred_num) |>
  mutate(count_sum = count)) |>
  distinct(week_obs_num, week_pred_num, count_sum)





county_peaks_both_col <- county_peaks_both_col |>
  mutate(week_obs_num = ifelse(week_obs_num > 30, week_obs_num - 52, week_obs_num),
         week_pred_num = ifelse(week_pred_num > 30, week_pred_num - 52, week_pred_num))

ggplot(county_peaks_both_col, aes(x = week_obs_num, y = week_pred_num)) +
  geom_point(aes(size = count_sum), color = "#e6a532", alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  labs(title = 'Syndromic Prediction vs. Claims\nPeak Timing',
       x = "Obs Peak Week",
       y = "Syndromic Peak Week") +
  scale_size_continuous(range = c(1, 10), breaks = c(1, 10, 100, 200, 400)) +
  theme_minimal()



mean(county_peaks_both_filt$week_pred_num - county_peaks_both_filt$week_obs_num)
cor(county_peaks_both_filt$week_pred_num, county_peaks_both_filt$week_obs_num)

















# Add month and year variables
flu_national$month <- month(flu_national$week_date)
flu_national$year <- year(flu_national$week_date)

# Load all cause data
all_cause <- read.csv('./raw/county_week_ac_v2.csv', header = T)

# Impute all cause counts
all_cause$patient_count_imp <- all_cause$all_cause
all_cause$patient_count_imp[all_cause$patient_count_imp == '<=5'] <- 
  sample(1:5, length(all_cause$patient_count_imp[all_cause$patient_count_imp == '<=5']), 
         replace= TRUE)
# Convert imputed count to numeric
all_cause$patient_count_imp <- as.numeric(all_cause$patient_count_imp)

# Convert year week to week date
all_cause <- all_cause |>
  mutate(week_date = sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", year_week)) |>
  mutate(week_date = ISOweek2date(week_date))

# Collapse to national level
all_cause_national <- all_cause |> group_by(week_date) |>
  mutate(all_cause_sum = patient_count_imp)) |>
  distinct(week_date, all_cause_sum) |>
  ungroup() |>
  filter(week_date < as.Date('2020-03-01')) |>
  filter(week_date > as.Date('2016-08-31')) |>
  mutate(all_cause_avg = mean(all_cause_sum))
  
# Quickly visualize all cause data
ggplot(all_cause_national) + 
  geom_line(aes(x = week_date, y = all_cause_sum / all_cause_avg), color = 'blue')

# Scale observed and predicted data by all cause
flu_national <- left_join(flu_national, all_cause_national, by = 'week_date')

# Calculate adjusted values scales by all cause data
flu_national$flu_obs_scale <- flu_national$flu_sum *
  (1/(flu_national$all_cause_sum / flu_national$all_cause_avg))
flu_national$flu_pred_scale <- flu_national$pred_sum *
  (1/(flu_national$all_cause_sum / flu_national$all_cause_avg))

# Quick visualization of the scale variables
ggplot(flu_national) + 
  geom_line(aes(x = week_date, y = flu_obs_scale), color = 'blue') +
  geom_line(aes(x = week_date, y = flu_pred_scale), color = 'red')

# Calculate summer averages
flu_national_summer <- flu_national |> 
  filter(month < 9 & month > 5) |> group_by(year) |>
  filter(year > 2016) |>
  mutate(pred_summ_mean = mean(flu_pred_scale),
         obs_summ_mean = mean(flu_obs_scale)) |>
  distinct(year, pred_summ_mean, 
           obs_summ_mean)

# Merge on summer data to the correct season (summer 2016 should align with the 2016-2017 flu season)
flu_national_summer[nrow(flu_national_summer) + 1, ] <- flu_national_summer[1, ] 
flu_national_summer$year[4] <- 2016
flu_national$year_2 <- ifelse(flu_national$month < 6, flu_national$year - 1, flu_national$year)
flu_national <- left_join(flu_national, flu_national_summer, by = c('year_2' = 'year'))

# Calculate rolling average
flu_national$roll_obs <- rollmean(flu_national$flu_obs_scale, 3, fill = NA, align = 'center')
flu_national$roll_pred <- rollmean(flu_national$flu_pred_scale, 3, fill = NA, align = 'center')

# Calculate Z score
flu_national$obs_z_roll <- (flu_national$roll_obs - flu_national$obs_summ_mean) / sd(flu_national$roll_obs, na.rm = T)
flu_national$obs_z <- (flu_national$flu_obs_scale - flu_national$obs_summ_mean) / sd(flu_national$flu_obs_scale, na.rm = T)
flu_national$pred_z_roll <- (flu_national$roll_pred - flu_national$pred_summ_mean) / sd(flu_national$roll_pred, na.rm = T)
flu_national$pred_z <- (flu_national$flu_pred_scale - flu_national$pred_summ_mean) / sd(flu_national$flu_pred_scale, na.rm = T)


# Create Season Variable
flu_national$Season <- ifelse(flu_national$year == 2016, '2016-2017', '')
flu_national$Season <- ifelse(flu_national$year == 2017 & flu_national$month < 9, '2016-2017', flu_national$Season)
flu_national$Season <- ifelse(flu_national$year == 2017 & flu_national$month > 8, '2017-2018', flu_national$Season)
flu_national$Season <- ifelse(flu_national$year == 2018 & flu_national$month < 9, '2017-2018', flu_national$Season)
flu_national$Season <- ifelse(flu_national$year == 2018 & flu_national$month > 8, '2018-2019', flu_national$Season)
flu_national$Season <- ifelse(flu_national$year == 2019 & flu_national$month < 9, '2018-2019', flu_national$Season)
flu_national$Season <- ifelse(flu_national$year == 2019 & flu_national$month > 8, '2019-2020', flu_national$Season)
flu_national$Season <- ifelse(flu_national$year == 2020, '2019-2020', flu_national$Season)

flu_national$pred_above <- ifelse(flu_national$pred_z > 0.1, 1, 0)
flu_national$obs_above <- ifelse(flu_national$obs_z > 0.1, 1, 0)

flu_national <- left_join(flu_national, ili_nat[, c(3, 6)], by = 'week_date')
flu_national <- flu_national |>
  rename('ili_z' = 'z_score')

flu_national <- left_join(flu_national, nrevss_nat[, c(3, 6)], by = 'week_date')
flu_national <- flu_national |>
  rename('nrevss_z' = 'z_score')

flu_national <- left_join(flu_national, flu_surv_nat[, c(2, 6)], by = 'week_date')
flu_national <- flu_national |>
  rename('surv_z' = 'z_score')


flu_national$ili_above <- ifelse(flu_national$ili_z > 0.1, 1, 0)
flu_national$nrevss_above <- ifelse(flu_national$nrevss_z > 0.1, 1, 0)
flu_national$surv_above <- ifelse(flu_national$surv_z > 0.1, 1, 0)

flu_national_duration <- flu_national |>
  group_by(Season) |>
  mutate(pred_duration = pred_above),
         obs_duration = obs_above),
         ili_duration = ili_above),
         nrevss_duration = nrevss_above),
         surv_duration = surv_above, na.rm = T)) |>
  distinct(Season, pred_duration, obs_duration, ili_duration, 
           nrevss_duration, surv_duration)

flu_national_duration_long <- flu_national_duration |>
  pivot_longer(
    cols = pred_duration:surv_duration, 
    names_to = "Source",
    values_to = "value"
  ) |>
  filter(Season != '2019-2020')


flu_national_duration_long$Source[flu_national_duration_long$Source == 'pred_duration'] <- 'Predicted'
flu_national_duration_long$Source[flu_national_duration_long$Source == 'obs_duration'] <- 'Claims'
flu_national_duration_long$Source[flu_national_duration_long$Source == 'ili_duration'] <- 'ILI-NET'
flu_national_duration_long$Source[flu_national_duration_long$Source == 'nrevss_duration'] <- 'NREVSS'
flu_national_duration_long$Source[flu_national_duration_long$Source == 'surv_duration'] <- 'FluSurv-NET'

flu_national_duration_long_col <- flu_national_duration_long |>
  group_by(Source) |>
  mutate(min = min(value),
         max = max(value)) |>
  distinct(Source, min, max) |>
  filter(Source != 'Claims')

Lines <- c('Claims' = 'solid', 'FluSurv-NET' = 'solid', 
           'ILI-NET' = 'solid', 'NREVSS' = 'solid', 'Predicted' = 'dashed')

Colors <- c('Claims' = '#bf6fa7', 'FluSurv-NET' = '#029456', 
            'ILI-NET' = '#4363d8', 'NREVSS' = '#bf6fa7', 'Predicted' = 'gray39')

percent_theme <- theme(legend.position = "none",
                       axis.text = element_text(size=14),
                       axis.title = element_text(size=16),
                       strip.background = element_blank(),
                       plot.title = element_text(size=18),
                       legend.title = element_text(size = 16),
                       legend.text = element_text(size = 14))

flu_national_duration_plot <- ggplot(flu_national_duration_long_col, 
                                     aes(x = Source, group = Source,
                                         ymin = min, ymax = max)) +
  geom_errorbar(linewidth = 2.5, width = 0.25, aes(color = Source)) +
  scale_color_manual(values = Colors) + ylab('Number of Weeks') +
  labs(title = "National-Level Season Duration by Source") + theme_minimal() + 
  percent_theme

flu_national_duration_plot

ggplot(flu_national) + 
  geom_line(aes(x = week_date, y = obs_z_roll), color = 'blue') +
  geom_line(aes(x = week_date, y = pred_z_roll), color = 'red')

# Create dataset for predicted flu
flu_pred_z <- flu_national[, c(1, 17)]
flu_pred_z$Source <- 'Predicted'
flu_pred_z$region <- 'United States'
flu_pred_z$state_fips <- '00'
flu_pred_z <- flu_pred_z |> rename('z_score_roll' = 'pred_z_roll')

# Create dataset for observed flu
flu_obs_z <- flu_national[, c(1, 15)]
flu_obs_z$Source <- 'Claims'
flu_obs_z$region <- 'United States'
flu_obs_z$state_fips <- '00'
flu_obs_z <- flu_obs_z |> rename('z_score_roll' = 'obs_z_roll')

national_trends <- rbind(flu_pred_z, ili_nat[, c(1, 3, 4, 7, 8)], 
                         nrevss_nat[, c(1, 3, 4, 7, 8)], flu_surv_nat[, c(1, 2, 4, 7, 8)])

national_trends <- national_trends |> filter(week_date > as.Date('2016-09-05'))

percent_theme <- theme(legend.position = "bottom",
                       axis.text = element_text(size=14),
                       axis.title = element_text(size=16),
                       strip.background = element_blank(),
                       plot.title = element_text(size=18),
                       legend.title = element_text(size = 16),
                       legend.text = element_text(size = 14))

Lines <- c('FluSurv-NET' = 'solid', 
           'ILI-NET' = 'solid', 'NREVSS' = 'solid', 'Predicted' = 'dotted')
Alpha <- c('FluSurv-NET' = 0.65, 
           'ILI-NET' = 0.65, 'NREVSS' = 0.65, 'Predicted' = 1)
Colors <- c('FluSurv-NET' = '#029456', 
           'ILI-NET' = '#4363d8', 'NREVSS' = '#bf6fa7', 'Predicted' = 'gray39')
Width <- c('FluSurv-NET' = 1.5, 
           'ILI-NET' = 1.5, 'NREVSS' = 1.5, 'Predicted' = 4.5)

nat_trend <- ggplot(national_trends) + 
  geom_line(aes(y = z_score_roll , x = week_date, color = Source, linetype = Source,
                alpha = Source, linewidth = Source)) +
  scale_y_continuous(limits = c(-0.5, 5)) +
  scale_linetype_manual(values = Lines) +
  scale_alpha_manual(values = Alpha) + 
  scale_color_manual(values = Colors) + 
  scale_linewidth_manual(values = Width) + 
  guides(colour = guide_legend(nrow = 1), alpha = guide_legend(nrow = 1),
         linetype = guide_legend(nrow = 1)) +
  xlab('Week') + ylab('Standardized Value') +
  theme_minimal() + ggtitle('National Influenza Activity by Source') +
  percent_theme 

nat_trend

usa_albers_state <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_state.rds'))   # convert to sf

usa_albers_state$cont <- ifelse(usa_albers_state$STATEFP == '02' | 
                                  usa_albers_state$STATEFP == '15', 0, 1)
usa_albers_cont <- usa_albers_state %>% 
  group_by(cont) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()



# Set the map theme
map_theme <- theme(legend.position = "right",
                   axis.title = element_text(size=16),
                   strip.background = element_blank(),
                   plot.title = element_blank(),
                   panel.border = element_rect(fill=NA, linewidth = 0.6, color = '#FFFFFF00'),
                   legend.title = element_text(size = 16),
                   legend.text = element_text(size = 12),
                   legend.key.height = unit(0.7, 'cm'),
                   legend.key.width = unit(0.4, "cm"))

cont_map <- ggplot() +
  geom_sf(data = usa_albers_cont , aes(), fill = 'white', color= 'black', linewidth = 0.2) +
  theme_void() + map_theme 

nat_trend |> 
  ggdraw() +
  draw_plot(
    {
      cont_map
    },
    x = 0.78, 
    y = 0.78,
    width = 0.25, 
    height = 0.15)


flu_nat_corr <- left_join(flu_national[c(1, 18)], flu_national[c(1, 16)],
                            by = 'week_date')
flu_nat_corr <- left_join(flu_nat_corr, ili_nat[,c(3, 6)] |> rename('ili' = 'z_score'),
                          by = 'week_date')
flu_nat_corr <- left_join(flu_nat_corr, nrevss_nat[,c(3, 6)] |> rename('nrevss' = 'z_score'),
                          by = 'week_date')
flu_nat_corr <- left_join(flu_nat_corr, flu_surv_nat[,c(2, 6)] |> rename('flu_surv' = 'z_score'),
                          by = 'week_date')

corr <- c(
cor(flu_nat_corr$pred_z, flu_nat_corr$obs_z, use = "pairwise.complete.obs",
     method = 'pearson'),
cor(flu_nat_corr$pred_z, flu_nat_corr$ili, use = "pairwise.complete.obs",
    method = 'pearson'),
cor(flu_nat_corr$pred_z, flu_nat_corr$nrevss, use = "pairwise.complete.obs",
    method = 'pearson'),
cor(flu_nat_corr$pred_z, flu_nat_corr$flu_surv, use = "pairwise.complete.obs",
    method = 'pearson'))

Correlation <- c('0.91', '0.90', '0.94', '0.87')
Source <- c('Claims', 'ILI-NET', 'NREVSS', 'FluSurv-NET')

corr_nat <- data.frame(Source, corr, Correlation)

Colors <- c('Claims' = '#bf6fa7', 'FluSurv-NET' = '#029456', 
            'ILI-NET' = '#4363d8', 'NREVSS' = '#f86a38')

corr_nat_plot <- ggplot(corr_nat, aes(x = Source, y = corr, fill = Source,
                     color = Source)) +
  geom_col(alpha = 0.75, width = 0.65) +
  geom_text(aes(label = Correlation), 
            vjust = 1.5,       
            color = "white",    
            size = 6,
            fontface = "bold") +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1.0)) +
  scale_color_manual(values = Colors) + 
  scale_fill_manual(values = Colors) + 
  labs(title = "National Correlation by Source", y = "Correlation") + 
  theme_minimal() + percent_theme + theme(legend.position = 'none')

# STATE LEVEL
flu_state <- flu_pred_dat |> 
  mutate(state_fips = substr(county_fips, 1, 2)) |>
  group_by(year_week_dt, state_fips) |>
  mutate(pred_sum = pred_count),
         flu_sum = patient_sum * flu)) |>
  distinct(year_week_dt, pred_sum, flu_sum) |>
  dplyr::filter(year_week_dt < as.Date('2020-03-01')) |>
  dplyr::filter(year_week_dt > as.Date('2016-08-31'))

# Restrict data to NYS
flu_nys <- flu_state |> dplyr::filter(state_fips == '36')

# Set year and month variables
flu_nys$month <- month(flu_nys$year_week_dt)
flu_nys$year <- year(flu_nys$year_week_dt)
flu_nys$year_week_dt <- as.Date(flu_nys$year_week_dt)

# Collapse to state level
all_cause_nys <- all_cause |> 
  mutate(state_fips = substr(county_fips, 1, 2)) |>
  group_by(week_date, state_fips) |>
  mutate(all_cause_sum = patient_count_imp)) |>
  distinct(week_date, all_cause_sum) |>
  dplyr::filter(state_fips == '36') |>
  ungroup() |>
  dplyr::filter(week_date < as.Date('2020-03-01')) |>
  dplyr::filter(week_date > as.Date('2016-08-31')) |>
  mutate(all_cause_avg = mean(all_cause_sum))

# Quickly visualize all cause data
ggplot(all_cause_nys) + 
  geom_line(aes(x = week_date, y = all_cause_sum / all_cause_avg), color = 'blue')

# Scale observed and predicted data by all cause
flu_nys <- left_join(flu_nys, all_cause_nys, by = c('year_week_dt' = 'week_date', 'state_fips'))

# Calculate adjusted values scales by all cause data
flu_nys$flu_obs_scale <- flu_nys$flu_sum * 
  (1/(flu_nys$all_cause_sum / flu_nys$all_cause_avg))
flu_nys$flu_pred_scale <- flu_nys$pred_sum *
  (1/(flu_nys$all_cause_sum / flu_nys$all_cause_avg))

# Quick visualization of the scale variables
ggplot(flu_nys) + 
  geom_line(aes(x = year_week_dt, y = flu_obs_scale), color = 'blue') +
  geom_line(aes(x =  year_week_dt, y = flu_pred_scale), color = 'red')

# Calculate summer averages
flu_state_summer <- flu_nys |> 
  dplyr::filter(month < 9 & month > 5) |> group_by(year, state_fips) |>
  dplyr::filter(year > 2016) |>
  mutate(pred_summ_mean = mean(flu_pred_scale),
         obs_summ_mean = mean(flu_obs_scale)) |>
  distinct(year, state_fips, pred_summ_mean, 
           obs_summ_mean)

# Merge on summer data to the correct season (summer 2016 should align with the 2016-2017 flu season)
flu_state_summer[nrow(flu_state_summer) + 1, ] <- flu_state_summer[1, ] 
flu_state_summer$year[4] <- 2016
flu_nys$year_2 <- ifelse(flu_nys$month < 6, flu_nys$year - 1, flu_nys$year)
flu_nys <- left_join(flu_nys, flu_state_summer, by = c('year_2' = 'year', 'state_fips' = 'state_fips'))

# Calculate rolling average
flu_nys$roll_obs <- rollmean(flu_nys$flu_obs_scale, 4, fill = NA, align = 'right')
flu_nys$roll_pred <- rollmean(flu_nys$flu_pred_scale, 4, fill = NA, align = 'right')

# Calculate Z score
flu_nys$obs_z <- (flu_nys$roll_obs - flu_nys$obs_summ_mean) / sd(flu_nys$roll_obs, na.rm = T)
flu_nys$pred_z <- (flu_nys$roll_pred - flu_nys$pred_summ_mean) / sd(flu_nys$roll_pred, na.rm = T)

ggplot(flu_nys) + 
  geom_line(aes(x =  year_week_dt, y = obs_z), color = 'blue') +
  geom_line(aes(x =  year_week_dt, y = pred_z), color = 'red')

# Create dataset for predicted flu
flu_pred_z <- flu_nys[, c(1, 2, 17)]
flu_pred_z$Source <- 'Predicted'
flu_pred_z$region <- 'New York'
flu_pred_z <- flu_pred_z |> rename('z_score_roll' = 'pred_z')

# Create dataset for observed flu
flu_obs_z <- flu_nys[, c(1, 2, 16)]
flu_obs_z$Source <- 'Claims'
flu_obs_z$region <- 'New York'
flu_obs_z <- flu_obs_z |> rename('z_score_roll' = 'obs_z')


ili_nys <- ili_state |> dplyr::filter(state_fips == '36')
flu_surv_nys <- flu_surv_state |> dplyr::filter(region == 'New York - Albany' | region == 'New York - Rochester') |>
  mutate(Source = ifelse(region == 'New York - Albany', 'FluSurv-NET (Albany)', 'FluSurv-NET (Rochester)'))
flu_surv_nys_mean <- flu_surv_nys |> group_by(week_date) |>
  mutate(z_score_roll = mean(z_score_roll, na.rm = T),
         Source = 'FluSurv-NET',
         Region = 'New York') |>
  distinct(week_date, region, state_fips, z_score_roll, Source)

nrevss_nys <- nrevss_state |> dplyr::filter(state_fips == '36')

flu_pred_z <- flu_pred_z |> rename('week_date' = 'year_week_dt')

nys_trends <- rbind(flu_pred_z, ili_nys[, c(1, 3, 4, 7, 8)], 
                         nrevss_nys[, c(1, 3, 4, 7, 8)], flu_surv_nys_mean)


nys_trends <- nys_trends |> dplyr::filter(week_date > as.Date('2016-09-05'))

percent_theme <- theme(legend.position = "bottom",
                       axis.text = element_text(size=14),
                       axis.title = element_text(size=16),
                       strip.background = element_blank(),
                       plot.title = element_text(size=18),
                       legend.title = element_text(size = 16),
                       legend.text = element_text(size = 14))

Lines <- c('FluSurv-NET' = 'solid', 
           'ILI-NET' = 'solid', 'NREVSS' = 'solid', 'Predicted' = 'dashed')
Alpha <- c('FluSurv-NET' = 0.65, 
           'ILI-NET' = 0.65, 'NREVSS' = 0.65, 'Predicted' = 1)
Colors <- c('FluSurv-NET' = '#029456', 
            'ILI-NET' = '#4363d8', 'NREVSS' = '#bf6fa7', 'Predicted' = 'gray39')
Width <- c('FluSurv-NET' = 1.5, 
           'ILI-NET' = 1.5, 'NREVSS' = 1.5, 'Predicted' = 2.5)


nys_trend <- ggplot(nys_trends[nys_trends$week_date < as.Date('2019-09-01'),]) + 
  geom_line(aes(y = z_score_roll , x = week_date, color = Source, linetype = Source,
                alpha = Source, linewidth = Source)) +
  scale_y_continuous(limits = c(-0.2, 5)) +
  scale_linetype_manual(values = Lines) +
  scale_alpha_manual(values = Alpha) + 
  scale_color_manual(values = Colors) + 
  scale_linewidth_manual(values = Width) + 
  guides(colour = guide_legend(nrow = 1), alpha = guide_legend(nrow = 1),
         linetype = guide_legend(nrow = 1)) +
  xlab('Week') + ylab('Standardized Value') +
  theme_minimal() + ggtitle('New York State Influenza Activity by Source') +
  percent_theme 

nys_trend

# Load state-level shape file
state <- read_sf(dsn = './raw/cb_2021_us_state_20m/', 
                 layer = 'cb_2021_us_state_20m')
state_ny <- state |> filter(STATEFP == '36')

nys_map <- ggplot() +
  geom_sf(data = state_ny, aes(), fill = 'white', color= 'black', linewidth = 0.4) +
  theme_void() + map_theme 

nys_trend |> 
  ggdraw() +
  draw_plot(
    {
      nys_map
    },
    x = 0.78, 
    y = 0.78,
    width = 0.25, 
    height = 0.15)


# STATE CORR

# Set year and month variables
flu_state$month <- month(flu_state$week_date)
flu_state$year <- year(flu_state$week_date)

# Collapse to state level
all_cause_state <- all_cause |> 
  mutate(state_fips = substr(county_fips, 1, 2)) |>
  group_by(week_date, state_fips) |>
  mutate(all_cause_sum = patient_count_imp)) |>
  distinct(week_date, all_cause_sum) |>
  ungroup() |>
  group_by(state_fips) |>
  filter(week_date < as.Date('2020-03-01')) |>
  filter(week_date > as.Date('2016-08-31')) |>
  mutate(all_cause_avg = mean(all_cause_sum))

# Quickly visualize all cause data
ggplot(all_cause_state) + 
  geom_line(aes(x = week_date, y = all_cause_sum / all_cause_avg), color = 'blue') +
  facet_wrap(~state_fips)

# Scale observed and predicted data by all cause
flu_state <- left_join(flu_state, all_cause_state, by = c('week_date', 'state_fips'))

# Calculate adjusted values scales by all cause data
flu_state$flu_obs_scale <- flu_state$flu_sum / (flu_state$all_cause_sum)
flu_state$flu_pred_scale <- flu_state$pred_sum / (flu_state$all_cause_sum)

# Quick visualization of the scale variables
ggplot(flu_state) + 
  geom_line(aes(x = week_date, y = flu_obs_scale), color = 'blue') +
  geom_line(aes(x = week_date, y = flu_pred_scale), color = 'red') +
  facet_wrap(~state_fips, scale = 'free')

# Calculate summer averages
flu_state_summer <- flu_state |> 
  filter(month < 9 & month > 5) |> group_by(year, state_fips) |>
  filter(year > 2016) |>
  mutate(pred_summ_mean = mean(flu_pred_scale),
         obs_summ_mean = mean(flu_obs_scale)) |>
  distinct(year, state_fips, pred_summ_mean, 
           obs_summ_mean)


flu_state_summer$count <- 1
flu_state_summer[flu_state_summer$year == 2017,]$count <- 2
flu_state_summer <- flu_state_summer |> uncount(count)
flu_state_summer <- flu_state_summer |> group_by(state_fips) |>
  mutate(index = row_number(),
         year = ifelse(index == 1, 2016, year))

# Merge on summer data to the correct season (summer 2016 should align with the 2016-2017 flu season)
flu_state$year_2 <- ifelse(flu_state$month < 6, flu_state$year - 1, flu_state$year)
flu_state <- left_join(flu_state, flu_state_summer, by = c('year_2' = 'year', 'state_fips' = 'state_fips'))

# Calculate rolling average
flu_state$roll_obs <- rollmean(flu_state$flu_obs_scale, 3, fill = NA, align = 'center')
flu_state$roll_pred <- rollmean(flu_state$flu_pred_scale, 3, fill = NA, align = 'center')

# Calculate Z score
flu_state$pred_z <- (flu_state$flu_pred_scale - flu_state$pred_summ_mean) / sd(flu_state$flu_pred_scale, na.rm = T)


ggplot(flu_state) + 
  geom_line(aes(x = week_date, y = pred_z), color = 'red') +
  facet_wrap(~state_fips, scale = 'free')


flu_state <- left_join(flu_state, ili_state[, c(3, 4, 6)], by = c('week_date', 'state_fips'))

state_corr <- flu_state |> group_by(state_fips) |>
  mutate(correlation = cor(pred_z, z_score, use = "pairwise.complete.obs",
                           method = 'pearson')) |>
  distinct(state_fips, correlation)


usa_albers_state <- left_join(usa_albers_state, state_corr, by = c('STATEFP' = 'state_fips'))

map_theme <- theme(legend.position = "right",
                   axis.title = element_text(size=16),
                   strip.background = element_blank(),
                   plot.title = element_text(size=18),
                   panel.border = element_rect(fill=NA, linewidth = 0.8, color = 'white'),
                   legend.title = element_text(size = 16),
                   legend.text = element_text(size = 12),
                   legend.key.height = unit(1.0, 'cm'),
                   legend.key.width = unit(0.4, "cm"))


state_map <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_state), aes(fill = correlation, group = STATEFP), color= 'black', linewidth = 0.10) +
  geom_sf(data = usa_albers_cont, aes(), fill = '#FFFFFF00', color= 'black', linewidth = 0.5) +
  scale_fill_gradient('Correlation   \n', low = "white", high = "#4363d8", na.value = "grey",
                      breaks = c(0.40, 0.60, 0.80, 1.00), 
                      limits = c(0.6, 1)) + 
  ggtitle('ILI-NET Correlation') +
  theme_void() + map_theme 

state_map

# STATE PEAKS

# Create Season Variable
flu_state$Season <- ifelse(flu_state$year == 2016, '2016-2017', '')
flu_state$Season <- ifelse(flu_state$year == 2017 & flu_state$month < 9, '2016-2017', flu_state$Season)
flu_state$Season <- ifelse(flu_state$year == 2017 & flu_state$month > 8, '2017-2018', flu_state$Season)
flu_state$Season <- ifelse(flu_state$year == 2018 & flu_state$month < 9, '2017-2018', flu_state$Season)
flu_state$Season <- ifelse(flu_state$year == 2018 & flu_state$month > 8, '2018-2019', flu_state$Season)
flu_state$Season <- ifelse(flu_state$year == 2019 & flu_state$month < 9, '2018-2019', flu_state$Season)
flu_state$Season <- ifelse(flu_state$year == 2019 & flu_state$month > 8, '2019-2020', flu_state$Season)
flu_state$Season <- ifelse(flu_state$year == 2020, '2019-2020', flu_state$Season)

# Calculate peaks
flu_state <- left_join(flu_state, nrevss_state[, c(3, 4, 6)], by = c('week_date', 'state_fips'))


ggplot(flu_state[flu_state$state_fips == '02', ]) + 
  geom_line(aes(x = week_date, y = pred_z), color = 'red') +
  facet_wrap(~state_fips, scale = 'free')


state_peaks <- flu_state |> filter(Season != '2019-2020') |>
  group_by(Season, state_fips) |>
  filter(state_fips != '12') |> # Florida is missing ILI data 
  filter(state_fips != '11') |> # DC has very small sample size
  mutate(pred_peak = peaks(flu_pred_scale, span=41, strict=FALSE, endbehavior = 1)) |>
  mutate(obs_peak = peaks(flu_sum, span=41, strict=FALSE, endbehavior = 1)) |>
  mutate(nrevss_peak = peaks(z_score.y, span=41, strict=FALSE, endbehavior = 1)) |>
  select(c(Season, state_fips, week_date, flu_pred_scale,
           pred_peak, nrevss_peak, obs_peak))


state_peaks_pred <- state_peaks |> filter(pred_peak == T) |>
  select(c(Season, state_fips, week_date, pred_peak)) |>
  rename('week_pred' = 'week_date') |>
  group_by(Season, state_fips) |>
  arrange(week_pred) |>
  slice_head(n = 1) |>
  mutate(week_pred_num = week(week_pred))

state_peaks_nrevss <- state_peaks |> filter(nrevss_peak == T) |>
  select(c(Season, state_fips, week_date, nrevss_peak)) |>
  rename('week_nrevss' = 'week_date') |>
  group_by(Season, state_fips) |>
  arrange(week_nrevss) |>
  slice_head(n = 1) |>
  mutate(week_nrevss_num = week(week_nrevss))


state_peaks_both <- left_join(state_peaks_pred, state_peaks_nrevss,
                               by = c('Season' = 'Season',
                                      'state_fips' = 'state_fips'))

state_peaks_both <- state_peaks_both |>
  mutate(week_nrevss_num = ifelse(week_nrevss_num > 26, week_nrevss_num - 52, week_nrevss_num),
         week_pred_num = ifelse(week_pred_num > 26, week_pred_num - 52, week_pred_num))

cor(state_peaks_both$week_pred_num, state_peaks_both$week_nrevss_num)

state_peaks_both_col <- state_peaks_both |>
  mutate(count = 1) |>
  group_by(week_nrevss_num, week_pred_num) |>
  mutate(`Number of Peaks` = count)) |>
  distinct(week_nrevss_num, week_pred_num, `Number of Peaks`)

state_peaks_plot <- ggplot(state_peaks_both_col, aes(x = week_nrevss_num, y = week_pred_num)) +
  geom_point(aes(size = `Number of Peaks`), color = "#4363d8", alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  labs(title = 'State-Level ILI-NET Peak Time Comparison',
       x = "ILI-NET Peak Epi Week",
       y = "Syndromic Peak Epi Week") +
  scale_size_continuous(range = c(3, 10), breaks = c(1, 5, 10, 20)) +
  theme_minimal() + percent_theme 

state_peaks_plot
mean(state_peaks_both$week_pred_num - state_peaks_both$week_ili_num)

# COUNTY LEVEL
flu_county <- flu_pred_dat |> 
  group_by(week_date, county_fips) |>
  mutate(pred_sum = pred_count),
         flu_sum = patient_sum * flu)) |>
  distinct(week_date, pred_sum, flu_sum, urban_code_binary, state_fips) |>
  filter(week_date < as.Date('2020-03-01')) |>
  filter(week_date > as.Date('2016-08-31'))

# Restrict data to NYS
flu_nyc <- flu_county |> filter(county_fips == '36061')

# Set year and month variables
flu_nyc$month <- month(flu_nyc$week_date)
flu_nyc$year <- year(flu_nyc$week_date)

# Collapse to national level
all_cause_nyc <- all_cause |> 
  group_by(week_date, county_fips) |>
  mutate(all_cause_sum = patient_count_imp)) |>
  distinct(week_date, all_cause_sum) |>
  filter(county_fips == '36061') |>
  ungroup() |>
  filter(week_date < as.Date('2020-03-01')) |>
  filter(week_date > as.Date('2016-08-31')) |>
  mutate(all_cause_avg = mean(all_cause_sum))

# Quickly visualize all cause data
ggplot(all_cause_nyc) + 
  geom_line(aes(x = week_date, y = all_cause_sum / all_cause_avg), color = 'blue')

# Scale observed and predicted data by all cause
flu_nyc <- left_join(flu_nyc, all_cause_nyc, by = c('week_date', 'county_fips'))

# Calculate adjusted values scales by all cause data
flu_nyc$flu_obs_scale <- flu_nyc$flu_sum *
  (1/(flu_nyc$all_cause_sum / flu_nyc$all_cause_avg))
flu_nyc$flu_pred_scale <- flu_nyc$pred_sum * 
  (1/(flu_nyc$all_cause_sum / flu_nyc$all_cause_avg))

# Quick visualization of the scale variables
ggplot(flu_nyc) + 
  geom_line(aes(x = week_date, y = flu_obs_scale), color = 'blue') +
  geom_line(aes(x = week_date, y = flu_pred_scale), color = 'red')

# Calculate summer averages
flu_county_summer <- flu_nyc |> 
  filter(month < 8 & month > 5) |> group_by(year, county_fips) |>
  filter(year > 2016) |>
  mutate(pred_summ_mean = mean(flu_pred_scale),
         obs_summ_mean = mean(flu_obs_scale)) |>
  distinct(year, county_fips, pred_summ_mean, 
           obs_summ_mean)

# Merge on summer data to the correct season (summer 2016 should align with the 2016-2017 flu season)
flu_county_summer[nrow(flu_county_summer) + 1, ] <- flu_county_summer[1, ] 
flu_county_summer$year[4] <- 2016
flu_nyc$year_2 <- ifelse(flu_nyc$month < 6, flu_nyc$year - 1, flu_nyc$year)
flu_nyc <- left_join(flu_nyc, flu_county_summer, by = c('year_2' = 'year', 'county_fips' = 'county_fips'))

# Calculate rolling average
flu_nyc$roll_obs <- rollmean(flu_nyc$flu_obs_scale, 3, fill = NA, align = 'center')
flu_nyc$roll_pred <- rollmean(flu_nyc$flu_pred_scale, 3, fill = NA, align = 'center')

# Calculate Z score
flu_nyc$obs_z <- (flu_nyc$roll_obs - flu_nyc$obs_summ_mean) / sd(flu_nyc$roll_obs, na.rm = T)
flu_nyc$pred_z <- (flu_nyc$roll_pred - flu_nyc$pred_summ_mean) / sd(flu_nyc$roll_pred, na.rm = T)

# Create dataset for predicted flu
flu_pred_z <- flu_nyc[, c(1, 2, 19)]
flu_pred_z$Source <- 'Predicted'
flu_pred_z <- flu_pred_z |> rename('z_score_roll' = 'pred_z')

# Create dataset for observed flu
flu_obs_z <- flu_nyc[, c(1, 2, 18)]
flu_obs_z$Source <- 'Claims'
flu_obs_z <- flu_obs_z |> rename('z_score_roll' = 'obs_z')


ili_nyc <- ili_state |> filter(state_fips == '99')
nyc_county <- nys_county |> filter(county_fips == '36061') |>
  mutate(county_fips = as.character(county_fips))

nyc_trends <- rbind(flu_pred_z, ili_nyc[, c(3, 4, 7, 8)], 
                    nyc_county[, c(1, 3, 7, 8)])


nyc_trends <- nyc_trends |> filter(week_date > as.Date('2016-09-05'))

percent_theme <- theme(legend.position = "bottom",
                       axis.text = element_text(size=14),
                       axis.title = element_text(size=16),
                       strip.background = element_blank(),
                       plot.title = element_text(size=18),
                       legend.title = element_text(size = 16),
                       legend.text = element_text(size = 14))

Lines <- c('NYS' = 'solid', 
           'ILI-NET' = 'solid', 'Predicted' = 'dashed')
Alpha <- c( 'NYS' = 0.70, 
           'ILI-NET' = 0.70, 'Predicted' = 1)
Colors <- c('NYS' = '#e6a532', 
            'ILI-NET' = '#4363d8', 'Predicted' = 'gray39')
Width <- c('NYS' = 1.5, 
           'ILI-NET' = 1.5, 'Predicted' = 2.5)

nyc_trend <- ggplot(nyc_trends) + 
  geom_line(aes(y = z_score_roll , x = week_date, color = Source,
                alpha = Source, linetype = Source, linewidth = Source)) +
  scale_y_continuous(limits = c(-0.2, 5)) +
  scale_linetype_manual(values = Lines) +
  scale_alpha_manual(values = Alpha) + 
  scale_color_manual(values = Colors) + 
  scale_linewidth_manual(values = Width) + 
  guides(colour = guide_legend(nrow = 1), alpha = guide_legend(nrow = 1),
         linetype = guide_legend(nrow = 1)) +
  xlab('Week') + ylab('Standardized Value') +
  theme_minimal() + ggtitle('New York County Influenza Activity by Source') +
  percent_theme 


nyc_trend

# Load county-level shape files
county <- read_sf(dsn = './raw/NYS_SWIS_Codes_1007335001743884489/', 
                  layer = 'NYS_SWIS_Codes')

county_nyc <- county |> filter(COUNTY_NAM == 'NewYork')

nyc_map <- ggplot() +
  geom_sf(data = county_nyc, aes(), fill = 'white', color= 'black', linewidth = 0.4) +
  theme_void() + map_theme 

nyc_trend |> 
  ggdraw() +
  draw_plot(
    {
      nyc_map
    },
    x = 0.78, 
    y = 0.78,
    width = 0.25, 
    height = 0.15)


# National CORRELATION


# COUNTY CORRELATION

flu_county_nys <- flu_county |>
  mutate(state_fips = substr(county_fips, 1, 2)) |>
  filter(state_fips == '36')

# Set year and month variables
flu_county_nys$month <- month(flu_county_nys$week_date)
flu_county_nys$year <- year(flu_county_nys$week_date)

# Collapse to national level
all_cause_county <- all_cause |> 
  group_by(week_date, county_fips) |>
  mutate(all_cause_sum = patient_count_imp)) |>
  distinct(week_date, all_cause_sum) |>
  ungroup() |>
  group_by(county_fips) |>
  filter(week_date < as.Date('2020-03-01')) |>
  filter(week_date > as.Date('2016-08-31')) |>
  mutate(all_cause_avg = mean(all_cause_sum))

# Quickly visualize all cause data
ggplot(all_cause_county[all_cause_county$county_fips == '36061',]) + 
  geom_line(aes(x = week_date, y = all_cause_sum / all_cause_avg), color = 'blue')

# Scale observed and predicted data by all cause
flu_county_nys <- left_join(flu_county_nys, all_cause_county, by = c('week_date', 'county_fips'))

# Calculate adjusted values scales by all cause data
flu_county_nys$flu_obs_scale <- flu_county_nys$flu_sum *
  (1/(flu_county_nys$all_cause_sum / flu_county_nys$all_cause_avg))
flu_county_nys$flu_pred_scale <- flu_county_nys$pred_sum * 
  (1/(flu_county_nys$all_cause_sum / flu_county_nys$all_cause_avg))

# Quick visualization of the scale variables
ggplot(flu_county_nys) + 
  geom_line(aes(x = week_date, y = flu_obs_scale), color = 'blue') +
  geom_line(aes(x = week_date, y = flu_pred_scale), color = 'red')

# Separate out county groups

# Load county-group cross walk
county_group <- read.csv('./raw/county_fips_grp_pop_wts.csv')

# Merge on rsv/symptom data
county_group_nys <- left_join(flu_county_nys, county_group,
                         by = c('county_fips' = 'county_fips_grp'))

# Disaggregate counties
county_group_nys <- county_group_nys |>
  mutate(flu_pred_scale = flu_pred_scale * pop_wt,
         flu_obs_scale = flu_obs_scale * pop_wt) 

# Calculate summer averages
flu_county_summer <- county_group_nys |> 
  filter(month < 8 & month > 5) |> group_by(year, county_fips.y) |>
  filter(year > 2016) |>
  mutate(pred_summ_mean = mean(flu_pred_scale),
         obs_summ_mean = mean(flu_obs_scale)) |>
  distinct(year, county_fips.y, pred_summ_mean, 
           obs_summ_mean)

flu_county_summer$count <- 1
flu_county_summer[flu_county_summer$year == 2017,]$count <- 2
flu_county_summer <- flu_county_summer |> uncount(count)
flu_county_summer <- flu_county_summer |> group_by(county_fips.y) |>
  mutate(index = row_number(),
         year = ifelse(index == 1, 2016, year))

# Merge on summer data to the correct season (summer 2016 should align with the 2016-2017 flu season)
county_group_nys$year_2 <- ifelse(county_group_nys$month < 6, county_group_nys$year - 1, county_group_nys$year)
county_group_nys <- left_join(county_group_nys, flu_county_summer, by = c('year_2' = 'year', 'county_fips.y' = 'county_fips.y'))

# Calculate rolling average
county_group_nys$roll_obs <- rollmean(county_group_nys$flu_obs_scale, 3, fill = NA, align = 'center')
county_group_nys$roll_pred <- rollmean(county_group_nys$flu_pred_scale, 3, fill = NA, align = 'center')

# Calculate Z score
county_group_nys$pred_z <- (county_group_nys$flu_pred_scale - county_group_nys$pred_summ_mean) / sd(county_group_nys$flu_pred_scale, na.rm = T)




nys_county$county_fips <- as.numeric(nys_county$county_fips)
county_group_nys <- left_join(county_group_nys, nys_county, by = c('week_date',
                                                               'county_fips.y' = 'county_fips'))

nys_corr <- county_group_nys |> group_by(county_fips.y) |>
  mutate(correlation = cor(pred_z, z_score, use = "pairwise.complete.obs",
                           method = 'pearson')) |>
  distinct(county_fips, correlation) |>
  mutate(county_fips.y = as.character(county_fips.y))


county <- read_sf(dsn = './raw/cb_2021_us_county_20m/', 
                  layer = 'cb_2021_us_county_20m')

county_nys <- county |> filter(STATEFP == '36')

# Merge on prev data
usa_albers_county_nys <- left_join(county_nys, nys_corr, by = c('GEOID' = 'county_fips.y'))

# Set the map theme
map_theme <- theme(legend.position = "bottom",
                   axis.title = element_text(size=16),
                   strip.background = element_blank(),
                   plot.title = element_text(size=18),
                   panel.border = element_rect(fill=NA, linewidth = 0.8, color = 'white'),
                   legend.title = element_text(size = 16),
                   legend.text = element_text(size = 12),
                   legend.key.height = unit(0.5, 'cm'),
                   legend.key.width = unit(1.0, "cm"))


corr_map <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_nys), aes(fill = correlation, group = GEOID), color= 'black', linewidth = 0.1) +
  geom_sf(data = state_ny, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.5) +
  scale_fill_gradient('Correlation   \n', low = "white", high = "#e6a532", na.value = "grey",
                      breaks = c(0.50, 0.60, 0.70, 0.80, 0.90), 
                      limits = c(0.50, 0.90)) + 
  ggtitle('      County-Level NYS Lab-Confirmed Correlation') +
  theme_void() + map_theme 

corr_map 

library(splus2R)
county_group_nys

# Create Season Variable
county_group_nys$Season <- ifelse(county_group_nys$year == 2016, '2016-2017', '')
county_group_nys$Season <- ifelse(county_group_nys$year == 2017 & county_group_nys$month < 9, '2016-2017', county_group_nys$Season)
county_group_nys$Season <- ifelse(county_group_nys$year == 2017 & county_group_nys$month > 8, '2017-2018', county_group_nys$Season)
county_group_nys$Season <- ifelse(county_group_nys$year == 2018 & county_group_nys$month < 9, '2017-2018', county_group_nys$Season)
county_group_nys$Season <- ifelse(county_group_nys$year == 2018 & county_group_nys$month > 8, '2018-2019', county_group_nys$Season)
county_group_nys$Season <- ifelse(county_group_nys$year == 2019 & county_group_nys$month < 9, '2018-2019', county_group_nys$Season)
county_group_nys$Season <- ifelse(county_group_nys$year == 2019 & county_group_nys$month > 8, '2019-2020', county_group_nys$Season)
county_group_nys$Season <- ifelse(county_group_nys$year == 2020, '2019-2020', county_group_nys$Season)

# Calculate peaks
county_peaks <- county_group_nys |> filter(Season != '2019-2020') |>
  group_by(Season, county_fips.y) |>
  mutate(pred_peak = peaks(pred_sum, span=41, strict=FALSE, endbehavior = 1)) |>
  mutate(nys_peak = peaks(value, span=41, strict=FALSE, endbehavior = 1)) |>
  select(c(Season, county_fips.y, week_date, flu_pred_scale, value,
           pred_peak, nys_peak))

county_peaks_pred <- county_peaks |> filter(pred_peak == T) |>
  select(c(Season, county_fips.y, week_date, pred_peak)) |>
  rename('county_fips' = 'county_fips.y',
         'week_pred' = 'week_date') |>
  group_by(Season, county_fips) |>
  arrange(week_pred) |>
  slice_head(n = 1) |>
  mutate(week_pred_num = week(week_pred))

county_peaks_nys <- county_peaks |> filter(nys_peak == T) |>
  select(c(Season, county_fips.y, week_date, nys_peak)) |>
  rename('county_fips' = 'county_fips.y',
         'week_nys' = 'week_date') |>
  group_by(Season, county_fips) |>
  arrange(week_nys) |>
  slice_head(n = 1) |>
  mutate(week_nys_num = week(week_nys))
  

county_peaks_both <- left_join(county_peaks_pred, county_peaks_nys,
                               by = c('Season' = 'Season',
                                    'county_fips' = 'county_fips'))


county_peaks_both_col <- county_peaks_both |>
  mutate(count = 1) |>
  group_by(week_nys_num, week_pred_num) |>
  mutate(count_sum = count)) |>
  distinct(week_nys_num, week_pred_num, count_sum)
 
ggplot(county_peaks_both_col, aes(x = week_nys_num, y = week_pred_num)) +
  geom_point(aes(size = count_sum), color = "#e6a532", alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  labs(title = 'Syndromic Prediction vs. NYS Lab Confirmed\nPeak Timing',
       x = "NYS Peak Week",
       y = "Syndromic Peak Week") +
  scale_size_continuous(range = c(3, 9), breaks = c(1, 5, 10, 20)) +
  scale_y_continuous(limits = c(0, 15)) +
  scale_x_continuous(limits = c(0, 15)) +
  theme_minimal()




flu_county_peak <- flu_county

flu_county_peak$month <- month(flu_county_peak$week_date)
flu_county_peak$year <- year(flu_county_peak$week_date)


# Create Season Variable
flu_county_peak$Season <- ifelse(flu_county_peak$year == 2016, '2016-2017', '')
flu_county_peak$Season <- ifelse(flu_county_peak$year == 2017 & flu_county_peak$month < 9, '2016-2017', flu_county_peak$Season)
flu_county_peak$Season <- ifelse(flu_county_peak$year == 2017 & flu_county_peak$month > 8, '2017-2018', flu_county_peak$Season)
flu_county_peak$Season <- ifelse(flu_county_peak$year == 2018 & flu_county_peak$month < 9, '2017-2018', flu_county_peak$Season)
flu_county_peak$Season <- ifelse(flu_county_peak$year == 2018 & flu_county_peak$month > 8, '2018-2019', flu_county_peak$Season)
flu_county_peak$Season <- ifelse(flu_county_peak$year == 2019 & flu_county_peak$month < 9, '2018-2019', flu_county_peak$Season)
flu_county_peak$Season <- ifelse(flu_county_peak$year == 2019 & flu_county_peak$month > 8, '2019-2020', flu_county_peak$Season)
flu_county_peak$Season <- ifelse(flu_county_peak$year == 2020, '2019-2020', flu_county_peak$Season)

# Calculate peaks
county_peaks_usa <- flu_county_peak |> filter(Season != '2019-2020') |>
  group_by(Season, county_fips) |>
  mutate(pred_peak = peaks(pred_sum, span=101, strict=FALSE, endbehavior = 1)) |>
  mutate(obs_peak = peaks(flu_sum, span=101, strict=FALSE, endbehavior = 1)) |>
  select(c(Season, county_fips, week_date, pred_sum, flu_sum,
           pred_peak, obs_peak, urban_code_binary, state_fips))


county_peaks_pred <- county_peaks_usa |> filter(pred_peak == T) |>
  select(c(Season, county_fips, week_date, pred_peak, urban_code_binary, state_fips)) |>
  rename('county_fips' = 'county_fips',
         'week_pred' = 'week_date') |>
  group_by(Season, county_fips) |>
  arrange(week_pred) |>
  slice_head(n = 1) |>
  mutate(week_pred_num = week(week_pred))

county_peaks_obs <- county_peaks_usa |> filter(obs_peak == T) |>
  select(c(Season, county_fips, week_date, obs_peak)) |>
  rename('county_fips' = 'county_fips',
         'week_obs' = 'week_date') |>
  group_by(Season, county_fips) |>
  arrange(week_obs) |>
  slice_head(n = 1) |>
  mutate(week_obs_num = week(week_obs))


county_peaks_both <- left_join(county_peaks_pred, county_peaks_obs,
                               by = c('Season' = 'Season',
                                      'county_fips' = 'county_fips'))

all_cause_season <- read_csv('raw/county_season_ac_v2.csv')
county_group <- readRDS('tmp/county_group_pop.rds')


county_peaks_both <- left_join(county_peaks_both, all_cause_season, by = c('Season' = 'season', 'county_fips'))
county_peaks_both <- left_join(county_peaks_both, county_group, by = c('county_fips'))

county_peaks_both_filt <- county_peaks_both |>
  filter(county_group_pop > 10000) |>
  filter(all_cause > 10000)

county_peaks_both_col <- county_peaks_both |>
  filter(county_group_pop > 10000) |>
  filter(all_cause > 10000) |>
  mutate(count = 1) |>
  group_by(week_obs_num, week_pred_num) |>
  mutate(count_sum = count)) |>
  distinct(week_obs_num, week_pred_num, count_sum)





county_peaks_both_col <- county_peaks_both_col |>
  mutate(week_obs_num = ifelse(week_obs_num > 30, week_obs_num - 52, week_obs_num),
         week_pred_num = ifelse(week_pred_num > 30, week_pred_num - 52, week_pred_num))

ggplot(county_peaks_both_col, aes(x = week_obs_num, y = week_pred_num)) +
  geom_point(aes(size = count_sum), color = "#e6a532", alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
  labs(title = 'Syndromic Prediction vs. Claims\nPeak Timing',
       x = "Obs Peak Week",
       y = "Syndromic Peak Week") +
  scale_size_continuous(range = c(1, 10), breaks = c(1, 10, 100, 200, 400)) +
  theme_minimal()



mean(county_peaks_both_filt$week_pred_num - county_peaks_both_filt$week_obs_num)
cor(county_peaks_both_filt$week_pred_num, county_peaks_both_filt$week_obs_num)


county_peaks_both_filt$diff <- county_peaks_both_filt$week_pred - county_peaks_both_filt$week_obs



usa_albers_county <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_county.rds'))   # convert to sf

x_walk <- readRDS('./tmp/county_group_xwalk.rds')

# Collapse county map to county group
usa_albers_county_x_walk <- left_join(x_walk, usa_albers_county[, c(5, 6)], by = c('fips' = 'GEOID'))
usa_albers_county_group <- usa_albers_county_x_walk %>% 
  group_by(county_fips) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()

county_peaks_both_filt_filt <- county_peaks_both_filt |>
  filter(diff < 50 & diff > -50)

usa_albers_county_group <- left_join(usa_albers_county_group, 
                               county_peaks_both_filt_filt[county_peaks_both_filt_filt$Season == '2016-2017', c(2, 13)],
                               by = c('county_fips' = 'county_fips'))


ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_group), aes(fill = as.numeric(diff), group = county_fips), color= 'black', linewidth = 0.01) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.10) +
  scale_fill_gradient2(' \n', low = "blue", mid = 'white', high = "red", na.value = "grey") + 
  theme_void() + map_theme 



peak_time_comp <- county_peaks_both_filt |> group_by(state_fips) |>
  mutate(mean_diff = mean(week_pred_num - week_obs_num)) |>
  distinct(state_fips, mean_diff)

ggplot(county_peaks_both_filt) + 
  geom_point(aes(y = abs(week_pred - week_obs), 
                x = all_cause), color = 'blue')

cor(as.numeric(abs(county_peaks_both_filt$week_pred - county_peaks_both_filt$week_obs)), 
    county_peaks_both_filt$county_group_pop)

ggplot(flu_county[flu_county$county_fips == '13283_13107',]) + 
  geom_line(aes(x = week_date, y = pred_sum), color = 'blue') +
  geom_line(aes(x = week_date, y = flu_sum), color = 'red')





figure_3 <- cowplot::plot_grid(nat_trend |> 
                                 ggdraw() +
                                 draw_plot(
                                   {
                                     cont_map
                                   },
                                   x = 0.78, 
                                   y = 0.78,
                                   width = 0.25, 
                                   height = 0.15), 
                               nys_trend |> 
                                 ggdraw() +
                                 draw_plot(
                                   {
                                     nys_map
                                   },
                                   x = 0.78, 
                                   y = 0.78,
                                   width = 0.25, 
                                   height = 0.15),
                               nyc_trend |> 
                                 ggdraw() +
                                 draw_plot(
                                   {
                                     nyc_map
                                   },
                                   x = 0.78, 
                                   y = 0.78,
                                   width = 0.25, 
                                   height = 0.15),
                               flu_national_duration_plot, state_peaks_plot, corr_map,
                               nrow = 2,
                               rel_heights = c(0.7, 1),
                               labels = c("a", "b", "c", "d", "e", "f"),
                               label_size = 20)

# Save figure 
ggsave('./figs/figure_3_state.jpg', plot = figure_3, height = 10, width = 20)






flu_state_summer$count <- 1
flu_state_summer[flu_state_summer$year == 2017,]$count <- 2

flu_state_summer <- flu_state_summer |> uncount(count)

flu_state_summer <- flu_state_summer |> group_by(state_fips) |>
  mutate(index = row_number(),
         year = ifelse(index == 1, 2016, year))

flu_state$year_2 <- ifelse(flu_state$month < 7, flu_state$year - 1, flu_state$year)
flu_state <- left_join(flu_state, flu_state_summer, by = c('year_2' = 'year', 'state_fips'))

flu_state <- flu_state |> filter(week_date < as.Date('2020-03-01'))

flu_state$obs_z <- (flu_state$flu_sum - flu_state$obs_summ_mean) / sd(flu_state$flu_sum)
flu_state$pred_z <- (flu_state$pred_sum - flu_state$pred_summ_mean) / sd(flu_state$pred_sum)

flu_pred_z_state <- flu_state[, c(1, 2, 12)]
flu_pred_z_state$Source <- 'Predicted'
flu_pred_z_state <- flu_pred_z_state |> rename('value' = 'pred_z')

flu_obs_z_state <- flu_state[, c(1, 2, 11)]
flu_obs_z_state$Source <- 'Claims'
flu_obs_z_state <- flu_obs_z_state |> rename('value' = 'obs_z')

state_trends <- rbind(flu_pred_z_state[flu_pred_z_state$state_fips == '06',], 
                      flu_obs_z_state[flu_obs_z_state$state_fips == '06',], 
                      ili_state_z[ili_state_z$state_fips == '06',], 
                      nrevss_state_z[nrevss_state_z$state_fips == '06',], 
                      surv_state_z_fill[surv_state_z_fill$state_fips == '06',])

state_trends <- state_trends |> filter(week_date > as.Date('2016-09-05'))

thick <- c('Claims' = 0.5, 'FluSurv-NET' = 0.5, 'ILI-NET' = 0.5, 
           'NREVSS' = 0.5, 'Predicted' = 1.1
)

percent_theme <- theme(legend.position = "bottom",
                       axis.text = element_text(size=14),
                       axis.title = element_text(size=16),
                       strip.background = element_blank(),
                       plot.title = element_text(size=18),
                       legend.title = element_text(size = 16),
                       legend.text = element_text(size = 14))

cal_trend <- ggplot(state_trends) + geom_line(aes(y = value , x = week_date, color = Source, linewidth = Source)) +
  scale_linewidth_manual(values = thick) + xlab('Week') + ylab('Standardized Value') +
  theme_minimal() + ggtitle('California Influenza Activity by Source') +
  percent_theme 
cal_trend 

flu_obs_z <- flu_obs_z |> rename('obs_z' = 'value')
flu_pred_z <- flu_pred_z |> rename('pred_z' = 'value')
ili_z <- ili_z |> rename('ili_z' = 'value')
nrevss_z <- nrevss_z |> rename('nrevss_z' = 'value')
surv_z <- surv_z |> rename('surv_z' = 'value')


national_corr <- left_join(flu_pred_z, flu_obs_z[, c(1, 2)], by = c('week_date'))
national_corr <- left_join(national_corr, ili_z[, c(2, 4)], by = c('week_date'))
national_corr <- left_join(national_corr, nrevss_z[, c(2, 4)], by = c('week_date'))
national_corr <- left_join(national_corr, surv_z[, c(1, 4)], by = c('week_date'))

national_corr <- national_corr |> filter(week_date > as.Date('2016-09-05'))

national_corr$year <- year(national_corr$week_date)
national_corr$month <- month(national_corr$week_date)

national_corr$Season <- ifelse(national_corr$year == 2016, '2016-2017', '')
national_corr$Season <- ifelse(national_corr$year == 2017 & national_corr$month < 9, '2016-2017', national_corr$Season)
national_corr$Season <- ifelse(national_corr$year == 2017 & national_corr$month > 8, '2017-2018', national_corr$Season)
national_corr$Season <- ifelse(national_corr$year == 2018 & national_corr$month < 9, '2017-2018', national_corr$Season)
national_corr$Season <- ifelse(national_corr$year == 2018 & national_corr$month > 8, '2018-2019', national_corr$Season)
national_corr$Season <- ifelse(national_corr$year == 2019 & national_corr$month < 9, '2018-2019', national_corr$Season)
national_corr$Season <- ifelse(national_corr$year == 2019 & national_corr$month > 8, '2019-2020', national_corr$Season)
national_corr$Season <- ifelse(national_corr$year == 2020, '2019-2020', national_corr$Season)

national_obs <- national_corr |> group_by(region, Season) |>
  mutate(correlation = cor(pred_z, obs_z),
         Source = 'Claims') |> distinct(region, Season, Source, correlation)

national_ili <- national_corr |> group_by(region, Season) |>
  mutate(correlation = cor(pred_z, ili_z),
         Source = 'ILI-NET') |> distinct(region, Season, Source, correlation)

national_nrevss <- national_corr |> group_by(region, Season) |>
  mutate(correlation = cor(pred_z, nrevss_z),
         Source = 'NREVSS') |> distinct(region, Season, Source, correlation)

national_surv <- national_corr |> group_by(region, Season) |>
  mutate(correlation = cor(pred_z, surv_z, use = 'complete.obs'),
         Source = 'FluSurv-NET') |> distinct(region, Season, Source, correlation)

nat_corr <- rbind(national_obs, national_ili, national_nrevss, national_surv)



flu_obs_z_state <- flu_obs_z_state |> rename('obs_z' = 'value')
flu_pred_z_state <- flu_pred_z_state |> rename('pred_z' = 'value')
ili_state_z <- ili_state_z |> rename('ili_z' = 'value')
nrevss_state_z <- nrevss_state_z |> rename('nrevss_z' = 'value')
surv_state_z <- surv_state_z |> rename('surv_z' = 'value')


state_corr <- left_join(flu_pred_z_state, flu_obs_z_state[, c(1, 2, 3)], by = c('week_date', 'state_fips'))
state_corr <- left_join(state_corr, ili_state_z[, c(1, 2, 3, 4)], by = c('week_date', 'state_fips'))
state_corr <- left_join(state_corr, nrevss_state_z[, c(2, 3, 4)], by = c('week_date', 'state_fips'))
state_corr <- left_join(state_corr, surv_state_z[, c(2, 3, 4)], by = c('week_date', 'state_fips'))

state_corr <- state_corr |> filter(week_date > as.Date('2016-09-05'))

state_corr$year <- year(state_corr$week_date)
state_corr$month <- month(state_corr$week_date)

state_corr$Season <- ifelse(state_corr$year == 2016, '2016-2017', '')
state_corr$Season <- ifelse(state_corr$year == 2017 & state_corr$month < 9, '2016-2017', state_corr$Season)
state_corr$Season <- ifelse(state_corr$year == 2017 & state_corr$month > 8, '2017-2018', state_corr$Season)
state_corr$Season <- ifelse(state_corr$year == 2018 & state_corr$month < 9, '2017-2018', state_corr$Season)
state_corr$Season <- ifelse(state_corr$year == 2018 & state_corr$month > 8, '2018-2019', state_corr$Season)
state_corr$Season <- ifelse(state_corr$year == 2019 & state_corr$month < 9, '2018-2019', state_corr$Season)
state_corr$Season <- ifelse(state_corr$year == 2019 & state_corr$month > 8, '2019-2020', state_corr$Season)
state_corr$Season <- ifelse(state_corr$year == 2020, '2019-2020', state_corr$Season)

state_obs <- state_corr |> group_by(region, Season) |>
  mutate(correlation = cor(pred_z, obs_z),
         Source = 'Claims') |> distinct(region, Season, Source, correlation)

state_ili <- state_corr |> group_by(region, Season) |>
  mutate(correlation = cor(pred_z, ili_z),
         Source = 'ILI-NET') |> distinct(region, Season, Source, correlation)

state_nrevss <- state_corr |> group_by(region, Season) |>
  mutate(correlation = cor(pred_z, nrevss_z),
         Source = 'NREVSS') |> distinct(region, Season, Source, correlation)

state_surv <- state_corr |> group_by(region, Season) |>
  mutate(correlation = cor(pred_z, surv_z, use = 'na.or.complete'),
         Source = 'FluSurv-NET') |> distinct(region, Season, Source, correlation)

sta_corr <- rbind(state_obs, state_ili, state_nrevss, state_surv)


corr <- rbind(sta_corr, nat_corr)

merge <- ili_state_z[, c(1, 3)] |> distinct(region, state_fips)
corr <- left_join(corr, merge, by = c('region'))

corr[corr$region == 'United States',]$state_fips <- '00'

library(viridis)
library(forcats)
p <-ggplot(corr,aes(Source,fct_reorder(region, -as.numeric(state_fips)),fill=correlation))+
  geom_tile(color= "white",linewidth =0.1) + 
  scale_fill_viridis(name="Correlation",option ="D", limits = c(-1,1)) +
  facet_grid(cols = vars(Season)) + ylab('') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.position = 'bottom')
p




ggplot(flu_national[2:nrow(flu_national),]) + geom_line(aes(y = obs_z , x = week_date), color = 'black') +
  geom_line(aes(y = pred_z , x = week_date), color = 'red')



flu_state <- flu |> group_by(state_fips, week_date) |>
  mutate(pred_sum = pred_age_count),
         flu_sum = patient_sum * flu)) |>
  distinct(state_fips, week_date, pred_sum, flu_sum)





# flu_county_nys <- flu |>   filter(state_fips == '36') |>
# group_by(county_fips, week_date) |>
#   mutate(pred_sum = pred_age_count),
#          flu_sum = patient_sum * flu)) |>
#   distinct(county_fips, week_date, pred_sum, flu_sum)


all_cause_16_17_week_state <- all_cause_16_17_week |> group_by(state_fips, week_date) |>
  mutate(cause_sum = patient_count_imp)) |>
  distinct(state_fips, week_date, cause_sum)

flu_state <- left_join(flu_state, all_cause_16_17_week_state, by = c('state_fips', 'week_date'))

flu_state$inci <- flu_state$pred_sum / flu_state$cause_sum

ggplot(flu_state) + geom_line(aes(y = pred_sum , x = week_date), color = 'black') +
  # geom_line(aes(y = flu , x = week_date), color = 'red') +
  facet_wrap(vars(state_fips), scales = "free")


flu_state <- left_join(flu_state, ili_state[, c(6, 18, 19)], by = c('state_fips', 'week_date'))
flu_state$ili_rate <- as.numeric(flu_state$X.UNWEIGHTED.ILI)
flu_state <- flu_state |> group_by(state_fips) |>
  mutate(ili_norm = (ili_rate - mean(ili_rate)) / sd(ili_rate),
         pred_norm = (pred_sum - mean(pred_sum)) / sd(pred_sum)) 

ggplot(flu_state) + geom_line(aes(y = ili_norm , x = week_date), color = 'black') +
  geom_line(aes(y = pred_norm , x = week_date), color = 'red') +
  # geom_line(aes(y = flu , x = week_date), color = 'red') +
  facet_wrap(vars(state_fips), scales = "free")