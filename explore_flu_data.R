################################################################################
# File Name: explore_flu_data                                                  #
#                                                                              #
# Purpose:   Explore influenza medical claims data.                            #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load data                                                      #
#            3. Explore data                                                   #
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
library(splitstackshape)
library(sf)
library(cowplot)
library(maptools)
library(sp)
library(mapproj)
library(rgeos)
library(lubridate)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1CgyyYhIWFVG2aMITI-HBpKamvwwpQ_ED/CHC_data/syndromic_surveillance/2024-05-06/')

################
# 2. LOAD DATA #
################

# Influenza symptom data
flu_16_17 <- read.csv('flu_2016-09_2017-08_1b_imputed.csv')

# Re-set directory 
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

###################
# 3. EXPLORE DATA #
###################

flu_descriptives_gender <- flu_16_17 |> group_by(flu, patient_gender_code) |> 
  mutate(group_sum = sum(patient_count)) |> distinct(flu, patient_gender_code, group_sum) |>
  ungroup() |> group_by(flu) |> 
  mutate(flu_sum = sum(group_sum),
         group_prop = group_sum / flu_sum)

flu_descriptives_age <- flu_16_17 |> group_by(flu, age_grp) |> 
  mutate(group_sum = sum(patient_count)) |> distinct(flu, age_grp, group_sum) |>
  ungroup() |> group_by(flu) |> 
  mutate(flu_sum = sum(group_sum),
         group_prop = group_sum / flu_sum)

flu_descriptives_symp <- flu_16_17 |> mutate(symp_count = fever + myalgia + cough + 
                                         sore_throat + short_breath + hypoxemia + 
                                         chest_pain + bronchitis + nausea_vom + 
                                         diarrhea + fatigue + headache + congestion + sneezing,
                                       no_symptoms = ifelse(symp_count == 0, patient_count, 0),
                                       fever = ifelse(fever == 1, patient_count, 0),
                                       myalgia = ifelse( myalgia == 1, patient_count, 0),
                                       short_breath = ifelse(short_breath == 1, patient_count, 0),
                                       hypoxemia = ifelse(hypoxemia == 1, patient_count, 0),
                                       chest_pain = ifelse(chest_pain == 1, patient_count, 0),
                                       cough = ifelse(cough == 1, patient_count, 0),
                                       bronchitis = ifelse(bronchitis == 1, patient_count, 0),
                                       nausea_vom = ifelse(nausea_vom == 1, patient_count, 0),
                                       sore_throat = ifelse(sore_throat == 1, patient_count, 0),
                                       diarrhea = ifelse(diarrhea == 1, patient_count, 0),
                                       fatigue = ifelse(fatigue == 1, patient_count, 0),
                                       headache = ifelse(headache == 1, patient_count, 0),
                                       congestion = ifelse(congestion == 1, patient_count, 0),
                                       sneezing = ifelse(sneezing == 1, patient_count, 0)) |>
  group_by(flu) |> mutate(fever_sum = sum(fever),
                          myalgia_sum = sum(myalgia),
                          short_breath_sum = sum(short_breath),
                          hypoxemia_sum = sum(hypoxemia),
                          chest_pain_sum = sum(chest_pain),
                          cough_sum = sum(cough),
                          bronchitis_sum = sum(bronchitis),
                          nausea_vom_sum = sum(nausea_vom),
                          sore_throat_sum = sum(sore_throat),
                          diarrhea_sum = sum(diarrhea),
                          fatigue_sum = sum(fatigue),
                          headache_sum = sum(headache),
                          congestion_sum = sum(congestion),
                          sneezing_sum = sum(sneezing),
                          mean_symp = mean(symp_count),
                          patient_sum = sum(patient_count)) |>
  distinct(flu, fever_sum, myalgia_sum, short_breath_sum, hypoxemia_sum, chest_pain_sum, cough_sum,
           bronchitis_sum, nausea_vom_sum, sore_throat_sum, diarrhea_sum, fatigue_sum, headache_sum,
           congestion_sum, sneezing_sum, mean_symp, patient_sum, .keep_all = FALSE) |>
  mutate(fever_prop = fever_sum / patient_sum,
         myalgia_prop = myalgia_sum / patient_sum,
         short_breath_prop = short_breath_sum / patient_sum,
         hypoxemia_prop = hypoxemia_sum / patient_sum,
         chest_pain_prop = chest_pain_sum / patient_sum,
         cough_prop = cough_sum / patient_sum,
         bronchitis_prop = bronchitis_sum / patient_sum,
         nausea_vom_prop = nausea_vom_sum / patient_sum,
         sore_throat_prop = sore_throat_sum / patient_sum,
         diarrhea_prop = diarrhea_sum / patient_sum,
         fatigue_prop = fatigue_sum / patient_sum,
         headache_prop = headache_sum / patient_sum,
         congestion_prop = congestion_sum / patient_sum,
         sneezing_prop = sneezing_sum / patient_sum)


# Full Flu
# Re-set directory 
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

flu <- readRDS('./tmp/flu_full_05-06_1b.rds')


# Collapse by week
flu_week_count_func <- function(data) {
  flu_week_count <- data |> 
    group_by(year_month) |>
    mutate(flu_count_week = sum(flu),
           fever_count_week = sum(fever),
           myalgia_count_week = sum(myalgia),
           short_breath_count_week = sum(short_breath),
           hypoxemia_count_week = sum(hypoxemia),
           chest_pain_count_week = sum(chest_pain),
           cough_count_week = sum(cough),
           bronchitis_count_week = sum(bronchitis),
           nausea_count_week = sum(nausea_vom),
           sore_throat_count_week = sum(sore_throat),
           diarrhea_count_week = sum(diarrhea),
           fatigue_count_week = sum(fatigue),
           headache_count_week = sum(headache),
           congestion_count_week = sum(congestion),
           sneezing_count_week = sum(sneezing)) |>
    distinct(year_month, flu_count_week, fever_count_week, myalgia_count_week, hypoxemia_count_week,
             short_breath_count_week, cough_count_week, nausea_count_week, chest_pain_count_week,
             sore_throat_count_week, diarrhea_count_week, fatigue_count_week, bronchitis_count_week,
             headache_count_week, congestion_count_week, sneezing_count_week,
             .keep_all = FALSE) |>
    mutate(year_month_date = as.Date(paste(year_month, '1', sep = '-'), format = '%Y-%m-%d'))
  return(flu_week_count)
}

# Collapse all flu data to the week level across counties
flu_week <- flu_week_count_func(flu)

flu_week_count <- read_csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1CgyyYhIWFVG2aMITI-HBpKamvwwpQ_ED/CHC_data/flu_cases/2023_Dec/county_week_conf_flu.csv')
flu_week_count$conf_flu_imp <- flu_week_count$conf_flu
flu_week_count$conf_flu_imp[flu_week_count$conf_flu_imp == '<=5'] <- 
  sample(1:5, length(flu_week_count$conf_flu_imp[flu_week_count$conf_flu_imp == '<=5']), replace= TRUE)
flu_week_count_sum <- flu_week_count %>% group_by(year_week) %>%
  mutate(flu_week_sum = sum(as.numeric(conf_flu_imp)),
         year_week_date = as.Date(paste(year_week, 1, sep = '-'), format = '%Y-%U-%u')) %>%
  distinct(year_week_date, flu_week_sum)


flu_week <- left_join(flu_week, flu_week_count_sum, by = c('year_week_date' = 'year_week_date'))


fever_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_month_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_month_date, y = fever_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Fever') +
  theme(plot.title = element_text(hjust = 0.5))

cough_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_month_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_month_date, y = cough_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Cough') +
  theme(plot.title = element_text(hjust = 0.5))

myalgia_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_month_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_month_date, y = myalgia_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Myalgia') +
  theme(plot.title = element_text(hjust = 0.5))

short_breath_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_month_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_month_date, y = short_breath_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Short Breath') +
  theme(plot.title = element_text(hjust = 0.5))

chest_pain_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_month_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_month_date, y = chest_pain_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Chest Pain') +
  theme(plot.title = element_text(hjust = 0.5))

hypoxemia_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_month_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_month_date, y = hypoxemia_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Hypoxemia') +
  theme(plot.title = element_text(hjust = 0.5))

bronchitis_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_month_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_month_date, y = bronchitis_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Bronchitis') +
  theme(plot.title = element_text(hjust = 0.5))

nausea_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_month_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_month_date, y = nausea_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Nausea') +
  theme(plot.title = element_text(hjust = 0.5))

sore_throat_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_month_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_month_date, y = sore_throat_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Sore Throat') +
  theme(plot.title = element_text(hjust = 0.5))

diarrhea_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_month_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_month_date, y = diarrhea_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Diarrhea') +
  theme(plot.title = element_text(hjust = 0.5))

fatigue_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_month_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_month_date, y = fatigue_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Fatigue') +
  theme(plot.title = element_text(hjust = 0.5))

headache_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_month_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_month_date, y = headache_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Headache') +
  theme(plot.title = element_text(hjust = 0.5))

congestion_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_month_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_month_date, y = congestion_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Congestion') +
  theme(plot.title = element_text(hjust = 0.5))

sneezing_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_month_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_month_date, y = sneezing_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Sneezing') +
  theme(plot.title = element_text(hjust = 0.5))

plot <- plot_grid(fever_plot, cough_plot, bronchitis_plot, myalgia_plot, short_breath_plot, hypoxemia_plot, chest_pain_plot,
                  nausea_plot, sore_throat_plot, diarrhea_plot, fatigue_plot,
                  headache_plot, congestion_plot, sneezing_plot)

ggsave('./figs/symptom_week_flu_p2_month.jpg', plot = plot, height = 8, width = 14)

# Extra code 


flu_16_17_p2 <- read.csv('flu_2016-09_2017-08_p2_imputed.csv')
flu_17_18_p2 <- read.csv('flu_2017-09_2018-08_p2_imputed.csv')
flu_18_19 <- read.csv('flu_2018-09_2019-08_p2_imputed.csv')
flu_19_20 <- read.csv('flu_2019-09_2020-08_p2_imputed.csv')
flu_new_2 <- rbind(flu_16_17_p2, flu_17_18_p2)
remove(flu_16_17_p2, flu_17_18_p2)

flu_16_17 <- read.csv('flu_2016-09_2017-08_imputed.csv')
flu_17_18 <- read.csv('flu_2017-09_2018-08_imputed.csv')
flu_new_1 <- rbind(flu_16_17, flu_17_18)
remove(flu_16_17, flu_17_18)

flu_17_old <- read.csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1CgyyYhIWFVG2aMITI-HBpKamvwwpQ_ED/CHC_data/syndromic_surveillance/old_version/flu_2017_imputed.csv')

flu_new_2_sub <- flu_new_2 %>% filter(county_fips == '36061') %>% 
  mutate(year_week_date = as.Date(paste(year_week, 1, sep = '-'), format = '%Y-%U-%u')) %>%
  filter(year(year_week_date) == '2017')

flu_new_1_sub <- flu_new_1 %>% filter(county_fips == '36061') %>% 
  mutate(year_week_date = as.Date(paste(year_week, 1, sep = '-'), format = '%Y-%U-%u')) %>%
  filter(year(year_week_date) == '2017')

flu_17_old_sub <- flu_17_old %>% filter(county_fips == '36061') %>% 
  mutate(year_week_date = as.Date(paste(year_week, 1, sep = '-'), format = '%Y-%U-%u')) %>%
  filter(year(year_week_date) == '2017')

flu_new_2_sub_count <- flu_new_2 |>
  mutate(year_week_date = as.Date(paste(year_week, 1, sep = '-'), format = '%Y-%U-%u')) |>
  group_by(year_week_date) |>
  mutate(flu = ifelse(flu == 1, patient_count, flu),
         fever = ifelse(fever == 1, patient_count, flu),
         sneezing = ifelse(sneezing == 1, patient_count, flu),
         flu_count_week_new_2 = sum(flu),
         fever_count_week_new_2 = sum(fever),
         sneezing_count_week_new_2 = sum(sneezing)) |>
  distinct(year_week_date, flu_count_week_new_2, fever_count_week_new_2, sneezing_count_week_new_2, 
           .keep_all = FALSE) %>%
  filter(year(year_week_date) == '2017')

flu_new_1_sub_count <- flu_new_1 |>
  mutate(year_week_date = as.Date(paste(year_week, 1, sep = '-'), format = '%Y-%U-%u')) |>
  group_by(year_week_date) |>
  mutate(flu = ifelse(flu == 1, patient_count, flu),
         fever = ifelse(fever == 1, patient_count, flu),
         sneezing = ifelse(sneezing == 1, patient_count, flu),
         flu_count_week_new_1 = sum(flu),
         fever_count_week_new_1 = sum(fever),
         sneezing_count_week_new_1 = sum(sneezing)) |>
  distinct(year_week_date, flu_count_week_new_1, fever_count_week_new_1, sneezing_count_week_new_1,
           .keep_all = FALSE) %>%
  filter(year(year_week_date) == '2017')

flu_17_old_sub_count <- flu_17_old |>
  mutate(year_week_date = as.Date(paste(year_week, 1, sep = '-'), format = '%Y-%U-%u')) |> 
  group_by(year_week_date) |>
  mutate(flu = ifelse(flu == 'True', patient_count_imp, 0),
         fever = ifelse(fever == 'True', patient_count_imp, 0),
         sneezing = ifelse(sneezing == 'True', patient_count_imp, 0),
         flu_count_week_old = sum(flu),
         fever_count_week_old = sum(fever),
         sneezing_count_week_old = sum(sneezing)) |>
  distinct(year_week_date, flu_count_week_old, fever_count_week_old, sneezing_count_week_old,
           .keep_all = FALSE)

flu_count <- left_join(flu_new_2_sub_count, flu_new_1_sub_count,
                       by = c('year_week_date' = 'year_week_date'))
flu_count_all <- left_join(flu_count, flu_17_old_sub_count,
                           by = c('year_week_date' = 'year_week_date'))




flu_comp <- ggplot(flu_count_all) + 
  geom_line(aes(x = year_week_date, y = flu_count_week_old), color = 'black') +
  geom_line(aes(x = year_week_date, y = flu_count_week_new_1), color = '#E7861B') + 
  geom_line(aes(x = year_week_date, y = flu_count_week_new_2), color = '#00BA42') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Flu Comparison') +
  theme(plot.title = element_text(hjust = 0.5))
flu_comp

fever_comp <- ggplot(flu_count_all) + 
  geom_line(aes(x = year_week_date, y = fever_count_week_old), color = 'black') +
  geom_line(aes(x = year_week_date, y = fever_count_week_new_1), color = '#E7861B') + 
  geom_line(aes(x = year_week_date, y = fever_count_week_new_2), color = '#00BA42') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Fever Comparison') +
  theme(plot.title = element_text(hjust = 0.5))
fever_comp

sneeze_comp <- ggplot(flu_count_all) + 
  geom_line(aes(x = year_week_date, y = sneezing_count_week_old), color = 'black') +
  geom_line(aes(x = year_week_date, y = sneezing_count_week_new_1), color = '#E7861B') + 
  geom_line(aes(x = year_week_date, y = sneezing_count_week_new_2), color = '#00BA42') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Sneeze Comparison') +
  theme(plot.title = element_text(hjust = 0.5))
sneeze_comp




setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1CgyyYhIWFVG2aMITI-HBpKamvwwpQ_ED/CHC_data/')


flu <- read.csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1CgyyYhIWFVG2aMITI-HBpKamvwwpQ_ED/CHC_data/flu_cases/2023_Oct/county_season_conf_flu.csv')
all_cause <- read.csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1CgyyYhIWFVG2aMITI-HBpKamvwwpQ_ED/CHC_data/all_cause/2023-10-27/county_season_ac_Oct2023.csv')

map_dat <- left_join(all_cause, flu, by = c('county_fips' = 'county_fips', 'season' = 'season'))
map_dat$conf_flu <- ifelse(map_dat$conf_flu == '<=5', 
                           sample(1:5, length(map_dat[map_dat$conf_flu == '<=5', 1]), replace = TRUE), map_dat$conf_flu)
map_dat$prev <- as.numeric(map_dat$conf_flu) / as.numeric(map_dat$all_cause)


