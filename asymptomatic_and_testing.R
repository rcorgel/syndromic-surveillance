################################################################################
# File Name: predict_flu_data                                                  #
#                                                                              #
# Purpose:   Predict influenza from medical claims symptoms.                   #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load data                                                      #
#            3. Predict data                                                   #
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
library(ISOweek)
library(sf)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

###############
# 2. FLU DATA #
###############

flu_16_17 <- read_csv('raw/flu_2016-09_2017-08_1b.csv')
flu_17_18 <- read_csv('raw/flu_2017-09_2018-08_1b.csv')
flu_18_19 <- read_csv('raw/flu_2018-09_2019-08_1b.csv')
flu_19_20 <- read_csv('raw/flu_2019-09_2020-08_1b.csv')
flu_20_21 <- read_csv('raw/flu_2020-09_2021-08_1b.csv')
flu_21_22 <- read_csv('raw/flu_2021-09_2022-08_1b.csv')
flu_22_23 <- read_csv('raw/flu_2022-09_2023-08_1b.csv')

flu <- rbind(flu_16_17, flu_17_18, flu_18_19, flu_19_20, flu_20_21, flu_21_22, flu_22_23)
remove(flu_16_17, flu_17_18, flu_18_19, flu_19_20, flu_20_21, flu_21_22, flu_22_23)

# Filter to flu cases only
flu <- flu |> filter(flu == 1)

# Create state and time variables
flu$state_fips <- substr(flu$county_fips, 1, 2)
flu$year <- year(as.Date(paste(flu$year_month, 1, sep = '-'), format = '%Y-%m-%d'))
flu$month <- month(as.Date(paste(flu$year_month, 1, sep = '-'), format = '%Y-%m-%d'))
flu$month_date <- as.Date(paste(flu$year, flu$month, 1, sep = '-'), format = '%Y-%m-%d')

# Impute <=5 to a random number between 1 and 5
flu$patient_count_imp <- flu$patient_count
flu$patient_count_imp[flu$patient_count_imp == '<=5'] <- 
  sample(1:5, length(flu$patient_count_imp[flu$patient_count_imp == '<=5']), replace= TRUE)
flu$patient_count_imp <- as.numeric(flu$patient_count_imp)

# Count flu symptoms
flu <- flu |> mutate(symp_count = fever + myalgia + 
                       cough + sore_throat + 
                       short_breath + hypoxemia + 
                       chest_pain + bronchitis + 
                       nausea_vom + diarrhea + fatigue + 
                       headache + congestion + sneezing,
                     no_symptoms = ifelse(symp_count == 0, 1, 0))

asymp_month <- flu |> group_by(month_date) |>
  mutate(all_flu = sum(patient_count_imp)) |>
  filter(no_symptoms == 1) |>
  mutate(all_asymp = sum(patient_count_imp)) |>
  distinct(month_date, all_asymp, all_flu) |>
  mutate(asymp_perc = all_asymp / all_flu)

ggplot(asymp_month, aes(x = month_date, y = asymp_perc)) + 
  geom_line() + theme_minimal() + ylab('Asymptomatic Percent') + 
  xlab('Month') + ggtitle('Asymptomatic Flu Percent, United States') + ylim(0, 1)

asymp_county <- flu |> group_by(county_fips) |>
  mutate(all_flu = sum(patient_count_imp)) |>
  filter(no_symptoms == 1) |>
  mutate(all_asymp = sum(patient_count_imp)) |>
  distinct(county_fips, all_asymp, all_flu) |>
  mutate(asymp_perc = all_asymp / all_flu)

# Create county to county group x-walk
x_walk <- flu |> separate(county_fips, sep = "_", remove = FALSE, 
                                into = c('fips_1', 'fips_2', 'fips_3', 'fips_4',
                                         'fips_5', 'fips_6', 'fips_7', 'fips_8')) |>
  pivot_longer(cols=c('fips_1', 'fips_2', 'fips_3', 'fips_4', 'fips_5', 'fips_6', 
                      'fips_7', 'fips_8'),
               names_to='fips_num',
               values_to='fips') |>
  filter(!is.na(fips)) |> select(c(county_fips, fips))
x_walk$state_fips <- substr(x_walk$county_fips, 1, 2)
# Remove PR and unknown (99999)
x_walk <- x_walk |> filter(state_fips != '72') |>
  filter(county_fips != '99999')

# Load maps
usa_albers_state <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_state.rds'))   # convert to sf
usa_albers_county <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_county.rds')) # convert to sf

# Collapse county map to county group
usa_albers_county_x_walk <- left_join(x_walk, usa_albers_county[, c(5, 6)], by = c('fips' = 'GEOID'))
usa_albers_county_group <- usa_albers_county_x_walk %>% 
  group_by(county_fips) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()

usa_albers_county_group <- left_join(usa_albers_county_group, asymp_county,
                                     by = c('county_fips' = 'county_fips'))
# Create maps
map <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_group), aes(fill = asymp_perc, group = county_fips), color= 'black', linewidth = 0.10) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.3) +
  scale_fill_gradient('Asymptomatic Percent   ', low = "white", high = "black", limits = c(0, 1)) + ggtitle('Asymptomatic Flu Percent, 2016 - 2023') + 
  theme_void() + theme(legend.position = 'bottom',
                       plot.title = element_text(size = 16, hjust = 0.5),
                       panel.border = element_rect(fill=NA, linewidth = 0.8, color = 'white'),
                       legend.title = element_text(size = 12),
                       legend.text = element_text(size = 12),
                       legend.key.height = unit(0.4, 'cm'),
                       legend.key.width = unit(1.5, "cm")) 
map

ggplot(asymp_county, aes(x = asymp_perc)) + 
  geom_density() + theme_minimal() + xlab('Asymptomatic Percent') + 
  ggtitle('Asymptomatic Flu Percent, US Counties')

tests <- read_csv('raw/county_week_flu_tests_v2.csv')

tests$flu_tests_imp <- tests$flu_tests
tests$flu_tests_imp[tests$flu_tests_imp == '<=5'] <- 
  sample(1:5, length(tests$flu_tests_imp[tests$flu_tests_imp == '<=5']), replace= TRUE)
tests$flu_tests_imp <- as.numeric(tests$flu_tests_imp)

library(ISOweek)
tests$week_date <- sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", tests$year_week)
tests$week_date <- ISOweek2date(tests$week_date)

tests_week <- tests |> group_by(week_date) |>
  mutate(tests_sum = sum(flu_tests_imp)) |>
  distinct(week_date, tests_sum)

ggplot(tests_week, aes(x = week_date, y = tests_sum)) + 
  geom_line() + theme_minimal() + ylab('Number of Tests') + 
  xlab('Week') + ggtitle('Influenza Tests, 2016 - 2024')
