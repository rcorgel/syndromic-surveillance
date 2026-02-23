################################################################################
# File Name: explore_rsv_data                                                  #
#                                                                              #
# Purpose:   Explore RSV medical claims data.                                  #
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

rsv_16_17 <- read.csv('rsv_2016-09_2017-08_1d_imputed.csv')
rsv_17_18 <- read.csv('rsv_2017-09_2018-08_1d_imputed.csv')
rsv_18_19 <- read.csv('rsv_2018-09_2019-08_1d_imputed.csv')
rsv_19_20 <- read.csv('rsv_2019-09_2020-08_1d_imputed.csv')
# Don't include 2020-2023 due to low rsv circulation & covid-19 co-circulation
# Append
rsv <- rbind(rsv_16_17, rsv_17_18, rsv_18_19, rsv_19_20)
remove(rsv_16_17, rsv_17_18, rsv_18_19, rsv_19_20)

rsv_descriptives_gender <- rsv |> group_by(rsv, patient_gender_code) |> 
  mutate(group_sum = sum(patient_count)) |> distinct(rsv, patient_gender_code, group_sum) |>
  ungroup() |> group_by(rsv) |> 
  mutate(rsv_sum = sum(group_sum),
         group_prop = group_sum / rsv_sum)

rsv_descriptives_age <- rsv |> group_by(rsv, age_grp) |> 
  mutate(group_sum = sum(patient_count)) |> distinct(rsv, age_grp, group_sum) |>
  ungroup() |> group_by(rsv) |> 
  mutate(rsv_sum = sum(group_sum),
         group_prop = group_sum / rsv_sum)

rsv_descriptives_symp <- rsv |> mutate(symp_count = fever + loss_appetite + cough + 
                                               sore_throat + short_breath + hypoxemia + 
                                               bronchitis +
                                               fatigue + headache + congestion + sneezing,
                                             no_symptoms = ifelse(symp_count == 0, patient_count, 0),
                                             fever = ifelse(fever == 1, patient_count, 0),
                                             loss_appetite = ifelse(loss_appetite == 1, patient_count, 0),
                                             short_breath = ifelse(short_breath == 1, patient_count, 0),
                                             hypoxemia = ifelse(hypoxemia == 1, patient_count, 0),
                                             cough = ifelse(cough == 1, patient_count, 0),
                                             bronchitis = ifelse(bronchitis == 1, patient_count, 0),
                                             sore_throat = ifelse(sore_throat == 1, patient_count, 0),
                                             fatigue = ifelse(fatigue == 1, patient_count, 0),
                                             headache = ifelse(headache == 1, patient_count, 0),
                                             congestion = ifelse(congestion == 1, patient_count, 0),
                                             sneezing = ifelse(sneezing == 1, patient_count, 0)) |>
  group_by(rsv) |> mutate(fever_sum = sum(fever),
                          loss_appetite_sum = sum(loss_appetite),
                          short_breath_sum = sum(short_breath),
                          hypoxemia_sum = sum(hypoxemia),
                          cough_sum = sum(cough),
                          bronchitis_sum = sum(bronchitis),
                          sore_throat_sum = sum(sore_throat),
                          fatigue_sum = sum(fatigue),
                          headache_sum = sum(headache),
                          congestion_sum = sum(congestion),
                          sneezing_sum = sum(sneezing),
                          mean_symp = mean(symp_count),
                          patient_sum = sum(patient_count)) |>
  distinct(rsv, fever_sum, loss_appetite_sum, short_breath_sum, hypoxemia_sum, cough_sum,
           bronchitis_sum, sore_throat_sum, fatigue_sum, headache_sum,
           congestion_sum, sneezing_sum, mean_symp, patient_sum, .keep_all = FALSE) |>
  mutate(fever_prop = fever_sum / patient_sum,
         loss_appetite_prop = loss_appetite_sum / patient_sum,
         short_breath_prop = short_breath_sum / patient_sum,
         hypoxemia_prop = hypoxemia_sum / patient_sum,
         cough_prop = cough_sum / patient_sum,
         bronchitis_prop = bronchitis_sum / patient_sum,
         sore_throat_prop = sore_throat_sum / patient_sum,
         fatigue_prop = fatigue_sum / patient_sum,
         headache_prop = headache_sum / patient_sum,
         congestion_prop = congestion_sum / patient_sum,
         sneezing_prop = sneezing_sum / patient_sum)


# Collapse by week
rsv_week_count_func <- function(data) {
  rsv_week_count <- data |> 
    group_by(year_month) |>
    mutate(rsv_count_week = sum(rsv),
           fever_count_week = sum(fever),
           loss_appetite_count_week = sum(loss_appetite),
           short_breath_count_week = sum(short_breath),
           hypoxemia_count_week = sum(hypoxemia),
           cough_count_week = sum(cough),
           bronchitis_count_week = sum(bronchitis),
           sore_throat_count_week = sum(sore_throat),
           fatigue_count_week = sum(fatigue),
           headache_count_week = sum(headache),
           congestion_count_week = sum(congestion),
           sneezing_count_week = sum(sneezing)) |>
    distinct(year_month, rsv_count_week, fever_count_week, loss_appetite_count_week, hypoxemia_count_week,
             short_breath_count_week, cough_count_week,
             sore_throat_count_week, fatigue_count_week, bronchitis_count_week,
             headache_count_week, congestion_count_week, sneezing_count_week,
             .keep_all = FALSE) |>
    mutate(year_month_date = as.Date(paste(year_month, '1', sep = '-'), format = '%Y-%m-%d'))
  return(rsv_week_count)
}

# Collapse all flu data to the week level across counties
rsv_week <- rsv_week_count_func(rsv)

flu_week_count <- read_csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1CgyyYhIWFVG2aMITI-HBpKamvwwpQ_ED/CHC_data/flu_cases/2023_Dec/county_week_conf_flu.csv')
flu_week_count$conf_flu_imp <- flu_week_count$conf_flu
flu_week_count$conf_flu_imp[flu_week_count$conf_flu_imp == '<=5'] <- 
  sample(1:5, length(flu_week_count$conf_flu_imp[flu_week_count$conf_flu_imp == '<=5']), replace= TRUE)
flu_week_count_sum <- flu_week_count %>% group_by(year_week) %>%
  mutate(flu_week_sum = sum(as.numeric(conf_flu_imp)),
         year_week_date = as.Date(paste(year_week, 1, sep = '-'), format = '%Y-%U-%u')) %>%
  distinct(year_week_date, flu_week_sum)


flu_week <- left_join(flu_week, flu_week_count_sum, by = c('year_week_date' = 'year_week_date'))


fever_plot <- ggplot(rsv_week) + 
  geom_line(aes(x = year_month_date, y = rsv_count_week), color = '#00BA42') +
  geom_line(aes(x = year_month_date, y = fever_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Fever') +
  theme(plot.title = element_text(hjust = 0.5))

cough_plot <- ggplot(rsv_week) + 
  geom_line(aes(x = year_month_date, y = rsv_count_week), color = '#00BA42') +
  geom_line(aes(x = year_month_date, y = cough_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Cough') +
  theme(plot.title = element_text(hjust = 0.5))

loss_appetite_plot <- ggplot(rsv_week) + 
  geom_line(aes(x = year_month_date, y = rsv_count_week), color = '#00BA42') +
  geom_line(aes(x = year_month_date, y = loss_appetite_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Loss Appetite') +
  theme(plot.title = element_text(hjust = 0.5))

short_breath_plot <- ggplot(rsv_week) + 
  geom_line(aes(x = year_month_date, y = rsv_count_week), color = '#00BA42') +
  geom_line(aes(x = year_month_date, y = short_breath_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Short Breath') +
  theme(plot.title = element_text(hjust = 0.5))

hypoxemia_plot <- ggplot(rsv_week) + 
  geom_line(aes(x = year_month_date, y = rsv_count_week), color = '#00BA42') +
  geom_line(aes(x = year_month_date, y = hypoxemia_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Hypoxemia') +
  theme(plot.title = element_text(hjust = 0.5))

bronchitis_plot <- ggplot(rsv_week) + 
  geom_line(aes(x = year_month_date, y = rsv_count_week), color = '#00BA42') +
  geom_line(aes(x = year_month_date, y = bronchitis_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Bronchitis') +
  theme(plot.title = element_text(hjust = 0.5))

sore_throat_plot <- ggplot(rsv_week) + 
  geom_line(aes(x = year_month_date, y = rsv_count_week), color = '#00BA42') +
  geom_line(aes(x = year_month_date, y = sore_throat_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Sore Throat') +
  theme(plot.title = element_text(hjust = 0.5))

fatigue_plot <- ggplot(rsv_week) + 
  geom_line(aes(x = year_month_date, y = rsv_count_week), color = '#00BA42') +
  geom_line(aes(x = year_month_date, y = fatigue_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Fatigue') +
  theme(plot.title = element_text(hjust = 0.5))

headache_plot <- ggplot(rsv_week) + 
  geom_line(aes(x = year_month_date, y = rsv_count_week), color = '#00BA42') +
  geom_line(aes(x = year_month_date, y = headache_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Headache') +
  theme(plot.title = element_text(hjust = 0.5))

congestion_plot <- ggplot(rsv_week) + 
  geom_line(aes(x = year_month_date, y = rsv_count_week), color = '#00BA42') +
  geom_line(aes(x = year_month_date, y = congestion_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Congestion') +
  theme(plot.title = element_text(hjust = 0.5))

sneezing_plot <- ggplot(rsv_week) + 
  geom_line(aes(x = year_month_date, y = rsv_count_week), color = '#00BA42') +
  geom_line(aes(x = year_month_date, y = sneezing_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Sneezing') +
  theme(plot.title = element_text(hjust = 0.5))

plot <- plot_grid(fever_plot, cough_plot, bronchitis_plot, loss_appetite_plot, short_breath_plot, hypoxemia_plot,
                  sore_throat_plot, fatigue_plot,
                  headache_plot, congestion_plot, sneezing_plot)

plot



