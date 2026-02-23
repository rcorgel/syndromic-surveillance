################################################################################
# File Name: explore_data                                                      #
#                                                                              #
# Purpose:   Explore influenza and COVID-19 medical claims data.               #
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
library(sf)
library(cowplot)
library(lubridate)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')


flu_16_17 <- read.csv('raw/flu_2016-09_2017-08_p2.csv')
flu_17_18 <- read.csv('raw/flu_2017-09_2018-08_p2.csv')
flu_18_19 <- read.csv('raw/flu_2018-09_2019-08_p2.csv')
flu_19_20 <- read.csv('raw/flu_2019-09_2020-08_p2.csv')
flu_20_21 <- read.csv('raw/flu_2020-09_2021-08_p2.csv')
flu_21_22 <- read.csv('raw/flu_2021-09_2022-08_p2.csv')
flu_22_23 <- read.csv('raw/flu_2022-09_2023-08_p2.csv')
flu_23_24 <- read.csv('raw/flu_2023-09_2024-08_p2.csv')

flu <- rbind(flu_16_17, flu_17_18, flu_18_19, flu_19_20, flu_20_21, flu_21_22, flu_22_23, flu_23_24)
remove(flu_16_17, flu_17_18, flu_18_19, flu_19_20, flu_20_21, flu_21_22, flu_22_23, flu_23_24)


flu$patient_count_imp <- flu$patient_count
flu$patient_count_imp[flu$patient_count_imp == '<=5'] <- 
  sample(1:5, length(flu$patient_count_imp[flu$patient_count_imp == '<=5']), replace= TRUE)

flu$state_fips <- substr(flu$county_fips, 1, 2)
# Remove PR and unknown (99999)
flu_filt <- flu |> filter(state_fips != '72') |>
  filter(county_fips != '99999')
library(ISOweek)
flu_filt$week_date <- sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", flu_filt$year_week)
flu_filt$week_date <- ISOweek2date(flu_filt$week_date)

flu_week <- flu_filt |> mutate(fever = ifelse(fever == 1, patient_count, 0),
                                       cough = ifelse(cough == 1, patient_count, 0),
                                       sore_throat = ifelse(sore_throat == 1, patient_count, 0)) |>
  group_by(week_date) |> mutate(fever_sum = sum(fever),
                          cough_sum = sum(cough),
                          sore_throat_sum = sum(sore_throat)) |>
  distinct(week_date, fever_sum, cough_sum,
           sore_throat_sum) 


viruses_us <- viruses |> group_by(week_date, Virus) |>
  mutate(count_sum = sum(as.numeric(count))) |>
  distinct(week_date, Virus, count_sum) |> ungroup() |>
  group_by(Virus) |>
  mutate(average = mean(count_sum),
         sd = sd(count_sum),
         z = (count_sum - average) / sd)

################
# 2. LOAD DATA #
################

# Influenza symptom data
flu <- readRDS('./tmp/flu_full_05-06_1b.rds')

###################
# 3. EXPLORE DATA #
###################

## Descriptive Statistics ##

flu_descriptives_gender <- flu |> group_by(flu, patient_gender_code) |> 
  mutate(group_sum = sum(patient_count)) |> distinct(flu, patient_gender_code, group_sum) |>
  ungroup() |> group_by(flu) |> 
  mutate(flu_sum = sum(group_sum),
         group_prop = group_sum / flu_sum)

flu_descriptives_age <- flu |> group_by(flu, age_grp) |> 
  mutate(group_sum = sum(patient_count)) |> distinct(flu, age_grp, group_sum) |>
  ungroup() |> group_by(flu) |> 
  mutate(flu_sum = sum(group_sum),
         group_prop = group_sum / flu_sum)

flu_descriptives_symp <- flu |> mutate(symp_count = fever + myalgia + cough + 
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
           congestion_sum, sneezing_sum, mean_symp, patient_sum) |> 
  mutate(fever_sum = sum(fever),
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
         patient_sum = sum(patient_count))
  








##############################
# Spatial Patterns by Season #
##############################

# Load all cause data by season
all_cause_season <- read_csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1CgyyYhIWFVG2aMITI-HBpKamvwwpQ_ED/CHC_data/all_cause/2023-10-27/county_season_ac_Oct2023.csv')

# Load flu data by season
flu_season <- read_csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1CgyyYhIWFVG2aMITI-HBpKamvwwpQ_ED/CHC_data/flu_cases/2023_Dec/county_season_conf_flu.csv')
flu_season$conf_flu_imp <- flu_season$conf_flu
# Randomly impute small numbers of flu
flu_season$conf_flu_imp[flu_season$conf_flu_imp == '<=5'] <- 
  sample(1:5, length(flu_season$conf_flu_imp[flu_season$conf_flu_imp == '<=5']), replace= TRUE)

flu_season <- left_join(flu_season, all_cause_season, by = c('county_fips' = 'county_fips',
                                                             'season' = 'season'))
flu_season$flu_prop <- as.numeric(flu_season$conf_flu_imp) / 
  as.numeric(flu_season$all_cause)

# Break up data by season
flu_16_17 <- flu_season[flu_season$season == '2016-2017',]
flu_16_17 <- flu_16_17 %>% rename('flu_prop_16_17' = 'flu_prop')
flu_17_18 <- flu_season[flu_season$season == '2017-2018',]
flu_17_18 <- flu_17_18 %>% rename('flu_prop_17_18' = 'flu_prop')
flu_18_19 <- flu_season[flu_season$season == '2018-2019',]
flu_18_19 <- flu_18_19 %>% rename('flu_prop_18_19' = 'flu_prop')
flu_19_20 <- flu_season[flu_season$season == '2019-2020',]
flu_19_20 <- flu_19_20 %>% rename('flu_prop_19_20' = 'flu_prop')

#################
# Create x-walk #
#################

# Create county to county group x-walk
x_walk <- flu_19_20 |> separate(county_fips, sep = "_", remove = FALSE, 
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
test <- left_join(usa_albers_county[, c(5, 6)], x_walk, by = c('GEOID' = 'fips'))


usa_albers_county_x_walk <- left_join(x_walk, usa_albers_county[, c(5, 6)], by = c('fips' = 'GEOID'))
usa_albers_county_group <- usa_albers_county_x_walk %>% 
  group_by(county_fips) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()

# Merge on data
usa_albers_county_group <- left_join(usa_albers_county_group, flu_16_17[, c(1, 6)],
                               by = c('county_fips' = 'county_fips'))
usa_albers_county_group <- left_join(usa_albers_county_group, flu_17_18[, c(1, 6)],
                                     by = c('county_fips' = 'county_fips'))
usa_albers_county_group <- left_join(usa_albers_county_group, flu_18_19[, c(1, 6)],
                                     by = c('county_fips' = 'county_fips'))
usa_albers_county_group <- left_join(usa_albers_county_group, flu_19_20[, c(1, 6)],
                                     by = c('county_fips' = 'county_fips'))

# Create maps
map <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_group), aes(fill = flu_prop_16_17, group = county_fips), color= 'black', linewidth = 0.15) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.4) +
  scale_fill_viridis_c('Prevalence', direction = -1) + ggtitle('Influenza Prevalence, 2016-17 Season') + 
  theme_void() + theme(legend.position = 'right',
                                     plot.title = element_text(size = 20, hjust = 0.5),
                                     panel.border = element_rect(fill=NA, linewidth = 0.8, color = 'white'),
                                     legend.title = element_text(size = 12),
                                     legend.text = element_text(size = 12),
                                     legend.key.size = unit(0.6, 'cm')) 
map

############
# Flu data #
############

# Collapse by week
flu_week_count_func <- function(data) {
  flu_week_count <- data |> 
    group_by(year_week) |>
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
    distinct(year_week, flu_count_week, fever_count_week, myalgia_count_week, hypoxemia_count_week,
             short_breath_count_week, cough_count_week, nausea_count_week, chest_pain_count_week,
             sore_throat_count_week, diarrhea_count_week, fatigue_count_week, bronchitis_count_week,
             headache_count_week, congestion_count_week, sneezing_count_week,
             .keep_all = FALSE) |>
    mutate(year_week_date = as.Date(paste(year_week, 1, sep = '-'), format = '%Y-%U-%u'))
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
  geom_line(aes(x = year_week_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_week_date, y = fever_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Fever') +
  theme(plot.title = element_text(hjust = 0.5))

cough_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_week_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_week_date, y = cough_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Cough') +
  theme(plot.title = element_text(hjust = 0.5))

myalgia_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_week_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_week_date, y = myalgia_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Myalgia') +
  theme(plot.title = element_text(hjust = 0.5))

short_breath_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_week_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_week_date, y = short_breath_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Short Breath') +
  theme(plot.title = element_text(hjust = 0.5))

chest_pain_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_week_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_week_date, y = chest_pain_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Chest Pain') +
  theme(plot.title = element_text(hjust = 0.5))

hypoxemia_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_week_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_week_date, y = hypoxemia_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Hypoxemia') +
  theme(plot.title = element_text(hjust = 0.5))

bronchitis_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_week_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_week_date, y = bronchitis_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Bronchitis') +
  theme(plot.title = element_text(hjust = 0.5))

nausea_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_week_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_week_date, y = nausea_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Nausea') +
  theme(plot.title = element_text(hjust = 0.5))

sore_throat_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_week_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_week_date, y = sore_throat_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Sore Throat') +
  theme(plot.title = element_text(hjust = 0.5))

diarrhea_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_week_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_week_date, y = diarrhea_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Diarrhea') +
  theme(plot.title = element_text(hjust = 0.5))

fatigue_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_week_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_week_date, y = fatigue_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Fatigue') +
  theme(plot.title = element_text(hjust = 0.5))

headache_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_week_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_week_date, y = headache_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Headache') +
  theme(plot.title = element_text(hjust = 0.5))

congestion_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_week_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_week_date, y = congestion_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Congestion') +
  theme(plot.title = element_text(hjust = 0.5))

sneezing_plot <- ggplot(flu_week) + 
  geom_line(aes(x = year_week_date, y = flu_count_week), color = '#00BA42') +
  geom_line(aes(x = year_week_date, y = sneezing_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('Sneezing') +
  theme(plot.title = element_text(hjust = 0.5))

plot <- plot_grid(fever_plot, cough_plot, bronchitis_plot, myalgia_plot, short_breath_plot, hypoxemia_plot, chest_pain_plot,
                  nausea_plot, sore_throat_plot, diarrhea_plot, fatigue_plot,
                  headache_plot, congestion_plot, sneezing_plot)

ggsave('./figs/symptom_week_flu_p2.jpg', plot = plot, height = 8, width = 14)

cor(flu_week$fever_count_week, flu_week$flu_count_week)
cor(flu_week$cough_count_week, flu_week$flu_count_week)
cor(flu_week$myalgia_count_week, flu_week$flu_count_week)
cor(flu_week$short_breath_count_week, flu_week$flu_count_week)
cor(flu_week$nausea_count_week, flu_week$flu_count_week)
cor(flu_week$sore_throat_count_week, flu_week$flu_count_week)
cor(flu_week$diarrhea_count_week, flu_week$flu_count_week)
cor(flu_week$fatigue_count_week, flu_week$flu_count_week)
cor(flu_week$headache_count_week, flu_week$flu_count_week)
cor(flu_week$congestion_count_week, flu_week$flu_count_week)
cor(flu_week$sneezing_count_week, flu_week$flu_count_week)

time_series <- left_join(all_cause_week, flu_week, by = c('year_week_date' = 
                                                            'year_week_date'))

plot <- ggplot(time_series) + 
  geom_line(aes(x = year_week_date, y = flu_count_week), color = '#E7861B') +
  geom_line(aes(x = year_week_date, y = all_claims), color = '#00A9FF') + 
  theme_minimal() + xlab('Date') + ylab('Claim Count') + ggtitle('All Claims vs. Flu') +
  theme(plot.title = element_text(hjust = 0.5))

ggsave('./figs/all_claims_week.jpg', plot = plot, height = 4, width = 7)

plot <- ggplot(time_series) + 
  geom_line(aes(x = year_week_date, y = flu_count_week/all_claims), color = '#00BA42') +
  theme_minimal() + xlab('Date') + ylab('Claim Proportion') + ggtitle('Flu / All Claims') +
  theme(plot.title = element_text(hjust = 0.5))

ggsave('./figs/all_claims_week_prop.jpg', plot = plot, height = 4, width = 7)


flu_week <- flu_week |> arrange(year_week)
flu_time_cor_func <- function(data, col_num) {
  date_list <- NULL
  corr_list <- NULL
  for (i in 10:208) {
    flu_week_sub <- data[i-8:i,]
    corr <- cor(flu_week_sub[, col_num], flu_week_sub$flu_count_week)
    date_list <- rbind(date_list, flu_week$year_week[i])
    corr_list <- rbind(corr_list, corr)
  }
  corr_over_time <- data.frame(cbind(date_list, corr_list))
  corr_over_time$X2 <- as.numeric(corr_over_time$X2)
  corr_over_time$X1 <- as.Date(paste(corr_over_time$X1, 1, sep = '-'), format = '%Y-%U-%u')
  return(corr_over_time)
}

fever_corr_time_plot <- ggplot(flu_time_cor_func(flu_week, 3)) + 
  geom_line(aes(x = X1, y = X2), color = '#00A9FF') +
  theme_minimal() + xlab('Date') + ylab('Correlation') + ggtitle('Fever vs. Flu') +
  theme(plot.title = element_text(hjust = 0.5)) + scale_y_continuous(limits = c(-1, 1))

myalgia_corr_time_plot <- ggplot(flu_time_cor_func(flu_week, 4)) + 
  geom_line(aes(x = X1, y = X2), color = '#00A9FF') +
  theme_minimal() + xlab('Date') + ylab('Correlation') + ggtitle('Myalgia vs. Flu') +
  theme(plot.title = element_text(hjust = 0.5)) + scale_y_continuous(limits = c(-1, 1))

short_breath_corr_time_plot <- ggplot(flu_time_cor_func(flu_week, 5)) + 
  geom_line(aes(x = X1, y = X2), color = '#00A9FF') +
  theme_minimal() + xlab('Date') + ylab('Correlation') + ggtitle('Short Breath vs. Flu') +
  theme(plot.title = element_text(hjust = 0.5)) + scale_y_continuous(limits = c(-1, 1))

cough_corr_time_plot <- ggplot(flu_time_cor_func(flu_week, 6)) + 
  geom_line(aes(x = X1, y = X2), color = '#00A9FF') +
  theme_minimal() + xlab('Date') + ylab('Correlation') + ggtitle('Cough vs. Flu') +
  theme(plot.title = element_text(hjust = 0.5)) + scale_y_continuous(limits = c(-1, 1))

nausea_corr_time_plot <- ggplot(flu_time_cor_func(flu_week, 7)) + 
  geom_line(aes(x = X1, y = X2), color = '#00A9FF') +
  theme_minimal() + xlab('Date') + ylab('Correlation') + ggtitle('Nausea vs. Flu') +
  theme(plot.title = element_text(hjust = 0.5)) + scale_y_continuous(limits = c(-1, 1))

sore_throat_corr_time_plot <- ggplot(flu_time_cor_func(flu_week, 8)) + 
  geom_line(aes(x = X1, y = X2), color = '#00A9FF') +
  theme_minimal() + xlab('Date') + ylab('Correlation') + ggtitle('Sore Throat vs. Flu') +
  theme(plot.title = element_text(hjust = 0.5)) + scale_y_continuous(limits = c(-1, 1))

diarrhea_corr_time_plot <- ggplot(flu_time_cor_func(flu_week, 9)) + 
  geom_line(aes(x = X1, y = X2), color = '#00A9FF') +
  theme_minimal() + xlab('Date') + ylab('Correlation') + ggtitle('Diarrhea vs. Flu') +
  theme(plot.title = element_text(hjust = 0.5)) + scale_y_continuous(limits = c(-1, 1))

fatigue_corr_time_plot <- ggplot(flu_time_cor_func(flu_week, 10)) + 
  geom_line(aes(x = X1, y = X2), color = '#00A9FF') +
  theme_minimal() + xlab('Date') + ylab('Correlation') + ggtitle('Fatigue vs. Flu') +
  theme(plot.title = element_text(hjust = 0.5)) + scale_y_continuous(limits = c(-1, 1))

headache_corr_time_plot <- ggplot(flu_time_cor_func(flu_week, 11)) + 
  geom_line(aes(x = X1, y = X2), color = '#00A9FF') +
  theme_minimal() + xlab('Date') + ylab('Correlation') + ggtitle('Headache vs. Flu') +
  theme(plot.title = element_text(hjust = 0.5)) + scale_y_continuous(limits = c(-1, 1))

congestion_corr_time_plot <- ggplot(flu_time_cor_func(flu_week, 12)) + 
  geom_line(aes(x = X1, y = X2), color = '#00A9FF') +
  theme_minimal() + xlab('Date') + ylab('Correlation') + ggtitle('Congestion vs. Flu') +
  theme(plot.title = element_text(hjust = 0.5)) + scale_y_continuous(limits = c(-1, 1))

sneezing_corr_time_plot <- ggplot(flu_time_cor_func(flu_week, 13)) + 
  geom_line(aes(x = X1, y = X2), color = '#00A9FF') +
  theme_minimal() + xlab('Date') + ylab('Correlation') + ggtitle('Sneezing vs. Flu') +
  theme(plot.title = element_text(hjust = 0.5)) + scale_y_continuous(limits = c(-1, 1))

plot <- plot_grid(fever_corr_time_plot, cough_corr_time_plot, myalgia_corr_time_plot, 
                  short_breath_corr_time_plot, 
                  nausea_corr_time_plot, sore_throat_corr_time_plot, diarrhea_corr_time_plot, 
                  fatigue_corr_time_plot,
                  headache_corr_time_plot, congestion_corr_time_plot, sneezing_corr_time_plot)

ggsave('./figs/symptom_week_cor.jpg', plot = plot, height = 8, width = 14)

# Collapse by county corr
flu_week_fips_count_func <- function(data) {
  flu_week_count <- data |> 
    mutate(flu_count = ifelse(flu == 'True', patient_count_imp, 0),
           fever_count = ifelse(fever == 'True', patient_count_imp, 0),
           myalgia_count = ifelse(myalgia == 'True', patient_count_imp, 0),
           short_breath_count = ifelse(short_breath == 'True', patient_count_imp, 0),
           cough_count = ifelse(cough == 'True', patient_count_imp, 0),
           nausea_count = ifelse(nausea == 'True', patient_count_imp, 0),
           sore_throat_count = ifelse(sore_throat == 'True', patient_count_imp, 0),
           diarrhea_count = ifelse(diarrhea == 'True', patient_count_imp, 0),
           fatigue_count = ifelse(fatigue == 'True', patient_count_imp, 0),
           headache_count = ifelse(headache == 'True', patient_count_imp, 0),
           congestion_count = ifelse(congestion == 'True', patient_count_imp, 0),
           sneezing_count = ifelse(sneezing == 'True', patient_count_imp, 0)) |>
    group_by(county_fips, year_week) |>
    mutate(flu_count_week = sum(flu_count),
           fever_count_week = sum(fever_count),
           myalgia_count_week = sum(myalgia_count),
           short_breath_count_week = sum(short_breath_count),
           cough_count_week = sum(cough_count),
           nausea_count_week = sum(nausea_count),
           sore_throat_count_week = sum(sore_throat_count),
           diarrhea_count_week = sum(diarrhea_count),
           fatigue_count_week = sum(fatigue_count),
           headache_count_week = sum(headache_count),
           congestion_count_week = sum(congestion_count),
           sneezing_count_week = sum(sneezing_count)) |>
    distinct(county_fips, year_week, flu_count_week, fever_count_week, myalgia_count_week,
             short_breath_count_week, cough_count_week, nausea_count_week,
             sore_throat_count_week, diarrhea_count_week, fatigue_count_week,
             headache_count_week, congestion_count_week, sneezing_count_week,
             .keep_all = FALSE) |>
    mutate(year_week_date = as.Date(paste(year_week, 1, sep = '-'), format = '%Y-%U-%u'))
  return(flu_week_count)
}

flu_2016_week_fips <- flu_week_fips_count_func(flu_2016)
flu_2017_week_fips <- flu_week_fips_count_func(flu_2017)
flu_2018_week_fips <- flu_week_fips_count_func(flu_2018)
flu_2019_week_fips <- flu_week_fips_count_func(flu_2019)
flu_week_fips <- rbind(flu_2016_week_fips, flu_2017_week_fips, flu_2018_week_fips, flu_2019_week_fips)

cor(flu_week$fever_count_week, flu_week$flu_count_week)
cor(flu_week$cough_count_week, flu_week$flu_count_week)
cor(flu_week$myalgia_count_week, flu_week$flu_count_week)
cor(flu_week$short_breath_count_week, flu_week$flu_count_week)
cor(flu_week$nausea_count_week, flu_week$flu_count_week)
cor(flu_week$sore_throat_count_week, flu_week$flu_count_week)
cor(flu_week$diarrhea_count_week, flu_week$flu_count_week)
cor(flu_week$fatigue_count_week, flu_week$flu_count_week)
cor(flu_week$headache_count_week, flu_week$flu_count_week)
cor(flu_week$congestion_count_week, flu_week$flu_count_week)
cor(flu_week$sneezing_count_week, flu_week$flu_count_week)

flu_week_fips_cor <- flu_week_fips |> group_by(county_fips) |> 
  mutate(fever_cor = cor(fever_count_week, flu_count_week),
         cough_cor = cor(cough_count_week, flu_count_week),
         myalgia_cor = cor(myalgia_count_week, flu_count_week),
         short_breath_cor = cor(short_breath_count_week, flu_count_week),
         nausea_cor = cor(nausea_count_week, flu_count_week),
         sore_throat_cor = cor(sore_throat_count_week, flu_count_week),
         diarrhea_cor = cor(diarrhea_count_week, flu_count_week),
         fatigue_cor = cor(fatigue_count_week, flu_count_week),
         headache_cor = cor(headache_count_week, flu_count_week),
         congestion_cor = cor(congestion_count_week, flu_count_week),
         sneezing_cor = cor(sneezing_count_week, flu_count_week)) |>
  distinct(county_fips, .keep_all = TRUE)

#################
# Create x-walk #
#################

# Create county to county group x-walk
x_walk <- flu_2018_county |> separate(county_fips, sep = "_", remove = FALSE, 
                                      into = c('fips_1', 'fips_2', 'fips_3', 'fips_4',
                                               'fips_5', 'fips_6', 'fips_7', 'fips_8')) |>
  pivot_longer(cols=c('fips_1', 'fips_2', 'fips_3', 'fips_4', 'fips_5', 'fips_6', 
                      'fips_7', 'fips_8'),
               names_to='fips_num',
               values_to='fips') |>
  filter(!is.na(fips)) |> select(-c(flu_symp_claims))
x_walk$state_fips <- substr(x_walk$county_fips, 1, 2)
# Remove PR and unknown (99999)
x_walk <- x_walk |> filter(state_fips != '72') |>
  filter(county_fips != '99999')

###################
# 3. EXPLORE DATA #
###################

# Load maps
usa_albers_state <- st_as_sf(readRDS('./tmp/usa_albers_state.rds'))   # convert to sf
usa_albers_county <- st_as_sf(readRDS('./tmp/usa_albers_county.rds')) # convert to sf

# Collapse county map to county group
usa_albers_county_x_walk <- left_join(x_walk, usa_albers_county[, c(5, 6)], by = c('fips' = 'GEOID'))
usa_albers_county_group <- usa_albers_county_x_walk %>% 
  group_by(county_fips) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()

# Collapse population data to county group
pop_2018_x_walk <- left_join(x_walk, pop_2018_format, by = c('fips' = 'fips'))
pop_2018_county_group <- pop_2018_x_walk |> group_by(county_fips) %>%
  mutate(pop_2018 = sum(POPESTIMATE2018)) |>
  distinct(county_fips, pop_2018, .keep_all = FALSE)

# Combine data
usa_albers_county_group <- left_join(usa_albers_county_group, pop_2018_county_group,
                                     by = c('county_fips' = 'county_fips'))
usa_albers_county_group <- left_join(usa_albers_county_group, flu_2018_county,
                                     by = c('county_fips' = 'county_fips'))
usa_albers_county_group <- left_join(usa_albers_county_group, all_cause_2018_county,
                                     by = c('county_fips' = 'county_fips'))
usa_albers_county_group <- left_join(usa_albers_county_group, flu_week_fips_cor,
                                     by = c('county_fips' = 'county_fips'))

# Map
map <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_group), aes(fill = all_claims/pop_2018, group = county_fips), color= 'black', linewidth = 0.15) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.4) +
  scale_fill_viridis_c('Claims per \nPerson', direction = -1) +
  theme_void() + ggtitle('') + theme(legend.position = 'right',
                                     plot.title = element_text(size = 30, hjust = 0.5),
                                     panel.border = element_rect(fill=NA, linewidth = 0.8, color = 'white'),
                                     legend.title = element_text(size = 24),
                                     legend.text = element_text(size = 24),
                                     legend.key.size = unit(1.2, 'cm')) 

map

ggsave('./figs/map_claims_pop.jpg', plot = map, height = 16, width = 28)

map <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_group), aes(fill = flu_symp_claims/all_claims, group = county_fips), color= 'black', linewidth = 0.15) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.4) +
  scale_fill_viridis_c('Flu Claims \nper All Claims', direction = -1) +
  theme_void() + ggtitle('') + theme(legend.position = 'right',
                                     plot.title = element_text(size = 30, hjust = 0.5),
                                     panel.border = element_rect(fill=NA, linewidth = 0.8, color = 'white'),
                                     legend.title = element_text(size = 24),
                                     legend.text = element_text(size = 24),
                                     legend.key.size = unit(1.2, 'cm')) 

map

ggsave('./figs/map_flu_claims.jpg', plot = map, height = 16, width = 28)


map <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_group), aes(fill = sneezing_cor, group = county_fips), color= 'black', linewidth = 0.15) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.4) +
  scale_fill_viridis_c('Sneezing Corr.', direction = -1) +
  theme_void() + ggtitle('') + theme(legend.position = 'right',
                                     plot.title = element_text(size = 30, hjust = 0.5),
                                     panel.border = element_rect(fill=NA, linewidth = 0.8, color = 'white'),
                                     legend.title = element_text(size = 24),
                                     legend.text = element_text(size = 24),
                                     legend.key.size = unit(1.2, 'cm')) 

map

ggsave('./figs/map_sneezing_cor.jpg', plot = map, height = 16, width = 28)





########
grouped_county <- x_walk_2 %>% 
  group_by(county_fips) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()

x_walk_grouped <- x_walk |> group_by(county_fips) |>
  distinct(county_fips, .keep_all = TRUE)

grouped_county <- left_join(grouped_county, x_walk_grouped, by = c('county_fips'))


ggplot(data = county_trans_sf) +
  geom_sf() + theme_void()

ggplot(data = st_as_sf(grouped_county)) +
  geom_sf() + theme_void()

ggsave('./figs/map.pdf', plot = last_plot(), height = 25, width = 30)


flu_2016_county

flu_2017 <- read.csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1GMtRVMDqj_Hg-871S7XP6TsrJn9WlgX2/syndromic_surveillance/flu_2017_imputed.csv')
flu_2018 <- read.csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1GMtRVMDqj_Hg-871S7XP6TsrJn9WlgX2/syndromic_surveillance/flu_2018_imputed.csv')
flu_2019 <- read.csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1GMtRVMDqj_Hg-871S7XP6TsrJn9WlgX2/syndromic_surveillance/flu_2019_imputed.csv')
flu_2020 <- read.csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1GMtRVMDqj_Hg-871S7XP6TsrJn9WlgX2/syndromic_surveillance/flu_2020_imputed.csv')
flu_2021 <- read.csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1GMtRVMDqj_Hg-871S7XP6TsrJn9WlgX2/syndromic_surveillance/flu_2021_imputed.csv')
flu_2022 <- read.csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1GMtRVMDqj_Hg-871S7XP6TsrJn9WlgX2/syndromic_surveillance/flu_2022_imputed.csv')

# COVID-19
covid_2020 <- read.csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1GMtRVMDqj_Hg-871S7XP6TsrJn9WlgX2/syndromic_surveillance/covid_2020_imputed.csv')
covid_2021 <- read.csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1GMtRVMDqj_Hg-871S7XP6TsrJn9WlgX2/syndromic_surveillance/covid_2021_imputed.csv')
covid_2022 <- read.csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1GMtRVMDqj_Hg-871S7XP6TsrJn9WlgX2/syndromic_surveillance/covid_2022_imputed.csv')


# Expand rows?










covid_week_count_func <- function(data) {
  covid_week_count <- data |> 
    mutate(covid_count = ifelse(covid == 'True', patient_count_imp, 0),
           fever_count = ifelse(fever == 'True', patient_count_imp, 0),
           myalgia_count = ifelse(myalgia == 'True', patient_count_imp, 0),
           short_breath_count = ifelse(short_breath == 'True', patient_count_imp, 0),
           cough_count = ifelse(cough == 'True', patient_count_imp, 0),
           nausea_count = ifelse(nausea == 'True', patient_count_imp, 0),
           sore_throat_count = ifelse(sore_throat == 'True', patient_count_imp, 0),
           loss_app_count = ifelse(loss_app == 'True', patient_count_imp, 0),
           fatigue_count = ifelse(fatigue == 'True', patient_count_imp, 0),
           loss_smell_count = ifelse(loss_smell == 'True', patient_count_imp, 0),
           hypoxemia_count = ifelse(hypoxemia == 'True', patient_count_imp, 0),
           chest_pain_count = ifelse(chest_pain == 'True', patient_count_imp, 0)) |>
    group_by(year_week) |>
    mutate(covid_count_week = sum(covid_count),
           fever_count_week = sum(fever_count),
           myalgia_count_week = sum(myalgia_count),
           short_breath_count_week = sum(short_breath_count),
           cough_count_week = sum(cough_count),
           nausea_count_week = sum(nausea_count),
           sore_throat_count_week = sum(sore_throat_count),
           loss_app_count_week = sum(loss_app_count),
           fatigue_count_week = sum(fatigue_count),
           loss_smell_count_week = sum(loss_smell_count),
           hypoxemia_count_week = sum(hypoxemia_count),
           chest_pain_count_week = sum(chest_pain_count)) |>
    distinct(year_week, covid_count_week, fever_count_week, myalgia_count_week,
             short_breath_count_week, cough_count_week, nausea_count_week,
             sore_throat_count_week, loss_app_count_week, fatigue_count_week,
             loss_smell_count_week, hypoxemia_count_week, chest_pain_count_week,
             .keep_all = FALSE) |>
    mutate(year_week_date = as.Date(paste(year_week, 1, sep = '-'), format = '%Y-%U-%u'))
  return(covid_week_count)
}

covid_2020_week <- covid_week_count_func(covid_2020)
covid_2021_week <- covid_week_count_func(covid_2021)
covid_2022_week <- covid_week_count_func(covid_2022)
covid_2023_week <- covid_week_count_func(covid_2023)
covid_week <- rbind(covid_2020_week, covid_2021_week, covid_2022_week, covid_2023_week)

fever_plot <- ggplot(covid_week) + 
  geom_line(aes(x = year_week_date, y = covid_count_week), color = '#00BA42') +
  geom_line(aes(x = year_week_date, y = fever_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Claim Count') + ggtitle('Fever') +
  theme(plot.title = element_text(hjust = 0.5))

cough_plot <- ggplot(covid_week) + 
  geom_line(aes(x = year_week_date, y = covid_count_week), color = '#00BA42') +
  geom_line(aes(x = year_week_date, y = cough_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Claim Count') + ggtitle('Cough') +
  theme(plot.title = element_text(hjust = 0.5))

myalgia_plot <- ggplot(covid_week) + 
  geom_line(aes(x = year_week_date, y = covid_count_week), color = '#00BA42') +
  geom_line(aes(x = year_week_date, y = myalgia_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Claim Count') + ggtitle('Myalgia') +
  theme(plot.title = element_text(hjust = 0.5))

short_breath_plot <- ggplot(covid_week) + 
  geom_line(aes(x = year_week_date, y = covid_count_week), color = '#00BA42') +
  geom_line(aes(x = year_week_date, y = short_breath_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Claim Count') + ggtitle('Short Breath') +
  theme(plot.title = element_text(hjust = 0.5))

nausea_plot <- ggplot(covid_week) + 
  geom_line(aes(x = year_week_date, y = covid_count_week), color = '#00BA42') +
  geom_line(aes(x = year_week_date, y = nausea_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Claim Count') + ggtitle('Nausea') +
  theme(plot.title = element_text(hjust = 0.5))

sore_throat_plot <- ggplot(covid_week) + 
  geom_line(aes(x = year_week_date, y = covid_count_week), color = '#00BA42') +
  geom_line(aes(x = year_week_date, y = sore_throat_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Claim Count') + ggtitle('Sore Throat') +
  theme(plot.title = element_text(hjust = 0.5))

loss_app_plot <- ggplot(covid_week) + 
  geom_line(aes(x = year_week_date, y = covid_count_week), color = '#00BA42') +
  geom_line(aes(x = year_week_date, y = loss_app_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Claim Count') + ggtitle('Lost App.') +
  theme(plot.title = element_text(hjust = 0.5))

fatigue_plot <- ggplot(covid_week) + 
  geom_line(aes(x = year_week_date, y = covid_count_week), color = '#00BA42') +
  geom_line(aes(x = year_week_date, y = fatigue_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Claim Count') + ggtitle('Fatigue') +
  theme(plot.title = element_text(hjust = 0.5))

loss_smell_plot <- ggplot(covid_week) + 
  geom_line(aes(x = year_week_date, y = covid_count_week), color = '#00BA42') +
  geom_line(aes(x = year_week_date, y = loss_smell_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Claim Count') + ggtitle('Lost Smell') +
  theme(plot.title = element_text(hjust = 0.5))

congestion_plot <- ggplot(covid_week) + 
  geom_line(aes(x = year_week_date, y = covid_count_week), color = '#00BA42') +
  geom_line(aes(x = year_week_date, y = hypoxemia_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Claim Count') + ggtitle('Hypoxemia') +
  theme(plot.title = element_text(hjust = 0.5))

sneezing_plot <- ggplot(covid_week) + 
  geom_line(aes(x = year_week_date, y = covid_count_week), color = '#00BA42') +
  geom_line(aes(x = year_week_date, y = chest_pain_count_week), color = '#E7861B') + 
  theme_minimal() + xlab('Date') + ylab('Claim Count') + ggtitle('Chest Pain') +
  theme(plot.title = element_text(hjust = 0.5))

plot <- plot_grid(fever_plot, cough_plot, myalgia_plot, short_breath_plot, 
                  nausea_plot, sore_throat_plot, loss_app_plot, fatigue_plot,
                  loss_smell_plot, congestion_plot, sneezing_plot)

ggsave('./figs/symptom_week_covid.jpg', plot = plot, height = 8, width = 14)

cor(covid_week$fever_count_week, covid_week$covid_count_week)
cor(covid_week$cough_count_week, covid_week$covid_count_week)
cor(covid_week$myalgia_count_week, covid_week$covid_count_week)
cor(covid_week$short_breath_count_week, covid_week$covid_count_week)
cor(covid_week$nausea_count_week, covid_week$covid_count_week)
cor(covid_week$sore_throat_count_week, covid_week$covid_count_week)
cor(covid_week$loss_app_count_week, covid_week$covid_count_week)
cor(covid_week$fatigue_count_week, covid_week$covid_count_week)
cor(covid_week$loss_smell_count_week, covid_week$covid_count_week)
cor(covid_week$hypoxemia_count_week, covid_week$covid_count_week)
cor(covid_week$chest_pain_count_week, covid_week$covid_count_week)




# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')


delay_date <- read.csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1CgyyYhIWFVG2aMITI-HBpKamvwwpQ_ED/CHC_data/date_lag/2024-06-13/date_lag_all_clm.csv')
# Randomly impute small numbers
delay_date$patient_count_imp <- delay_date$patient_count
delay_date$patient_count_imp[delay_date$patient_count == '<=3'] <- 
  sample(1:3, length(delay_date$patient_count_imp[delay_date$patient_count == '<=3']), replace= TRUE)
# Create plot
delay_date$lag <- as.numeric(as.Date(delay_date$chc_received_date) - as.Date(delay_date$claim_date))
delay_date$lag_cat <- cut(delay_date$lag, breaks = c(min(delay_date$lag), -1, 2, 7, 14, 28, 50, 70, max(delay_date$lag)),
                          labels = c("Negative","0-2 days","3-7 days","8-14 days","15-28 days", "29-50 days", "51-70 days", "70+ days"), 
                          include.lowest = TRUE)

delay_date_group <- delay_date %>% group_by(claim_date, lag_cat) %>% 
  mutate(lag_sum = sum(as.numeric(patient_count_imp))) %>% 
  distinct(claim_date, lag_cat, lag_sum) %>% 
  ungroup() %>% 
  group_by(claim_date) %>% 
  mutate(total_claims = sum(lag_sum),
         proportion = lag_sum / total_claims)

delay_date_group_dow <- delay_date_group %>% 
  mutate(day_of_week = wday(claim_date, label=TRUE)) %>%
  group_by(day_of_week, lag_cat) %>% 
  mutate(lag_sum_dow = sum(as.numeric(lag_sum))) %>% 
  distinct(day_of_week, lag_cat, lag_sum_dow) %>% 
  ungroup() %>% 
  group_by(day_of_week) %>% 
  mutate(total_claims = sum(lag_sum_dow),
         proportion = lag_sum_dow / total_claims)
  


cat <- data.frame(unique(delay_date$lag_cat))
names(cat) <- 'lag_cat'
cat$merge <- 1
claim_date <- data.frame(unique(delay_date$claim_date))
names(claim_date) <- 'claim_date'
claim_date$merge <- 1
shell <- dplyr::left_join(claim_date, cat, by = c('merge' = 'merge'))

graph_data <- left_join(shell, delay_date_group, by = c('claim_date' = 'claim_date', 
                                                        'lag_cat' = 'lag_cat'))
library(forcats)
ggplot(graph_data, aes(as.Date(claim_date), proportion, fill= forcats::fct_rev(lag_cat))) + 
  geom_bar(stat="summary", width=1) + scale_fill_brewer(palette = "RdYlGn", direction = 1) +
  theme_minimal() +
  labs(title = "Claim Delay by Date",
       x = "Date",
       y = "Proportion",
       fill = "Delay") 

ggplot(delay_date_group_dow, aes(day_of_week, proportion, fill= forcats::fct_rev(lag_cat))) + 
  geom_bar(stat="summary") + scale_fill_brewer(palette = "RdYlGn", direction = 1) +
  theme_minimal() +
  labs(title = "Claim Delay by Day of Week",
       x = "Date",
       y = "Proportion",
       fill = "Delay") 


delay_age <- read.csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1CgyyYhIWFVG2aMITI-HBpKamvwwpQ_ED/CHC_data/date_lag/2024-06-13/date_lag_age_all_clm.csv')

delay_age$patient_count_imp <- delay_age$patient_count
delay_age$patient_count_imp[delay_age$patient_count == '<=3'] <- 
  sample(1:3, length(delay_age$patient_count_imp[delay_age$patient_count == '<=3']), replace= TRUE)

delay_age$lag <- delay_age$chc_rec_claim_date_lag
delay_age$lag_cat <- cut(delay_age$lag, breaks = c(min(delay_age$lag), -1, 2, 7, 14, 28, 50, 70, max(delay_age$lag)),
                          labels = c("Negative","0-2 days","3-7 days","8-14 days","15-28 days", "29-50 days", "51-70 days", "70+ days"), 
                          include.lowest = TRUE)

delay_age_group <- delay_age %>% group_by(age_grp, lag_cat) %>% 
  mutate(lag_sum = sum(as.numeric(patient_count_imp))) %>% 
  distinct(age_grp, lag_cat, lag_sum) %>% 
  ungroup() %>% 
  group_by(age_grp) %>% 
  mutate(total_claims = sum(lag_sum),
         proportion = lag_sum / total_claims)

delay_age_group$age_cat <- ifelse(delay_age_group$age_grp == 'U', 'Unknown', '')
delay_age_group$age_cat <- ifelse(delay_age_group$age_grp == '0', '0-4', delay_age_group$age_cat)
delay_age_group$age_cat <- ifelse(delay_age_group$age_grp == '1', '5-12', delay_age_group$age_cat)
delay_age_group$age_cat <- ifelse(delay_age_group$age_grp == '2', '13-17', delay_age_group$age_cat)
delay_age_group$age_cat <- ifelse(delay_age_group$age_grp == '3', '18-49', delay_age_group$age_cat)
delay_age_group$age_cat <- ifelse(delay_age_group$age_grp == '4', '50-64', delay_age_group$age_cat)
delay_age_group$age_cat <- ifelse(delay_age_group$age_grp == '5', '65+', delay_age_group$age_cat)


ggplot(delay_age_group, aes(age_cat, proportion, fill= forcats::fct_rev(lag_cat))) + 
  geom_bar(stat="summary") + scale_fill_brewer(palette = "RdYlGn", direction = 1) +
  theme_minimal() +
  labs(title = "Claim Delay by Age Group",
       x = "Age Group",
       y = "Proportion",
       fill = "Delay") 

delay_state <- read.csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1CgyyYhIWFVG2aMITI-HBpKamvwwpQ_ED/CHC_data/date_lag/2024-06-13/date_lag_state_all_clm.csv')

delay_state$patient_count_imp <- delay_state$patient_count
delay_state$patient_count_imp[delay_state$patient_count == '<=3'] <- 
  sample(1:3, length(delay_state$patient_count_imp[delay_state$patient_count == '<=3']), replace= TRUE)

delay_state$lag <- delay_state$chc_rec_claim_date_lag
delay_state$lag_cat <- cut(delay_state$lag, breaks = c(min(delay_state$lag), -1, 2, 7, 14, 28, 50, 70, max(delay_state$lag)),
                         labels = c("Negative","0-2 days","3-7 days","8-14 days","15-28 days", "29-50 days", "51-70 days", "70+ days"), 
                         include.lowest = TRUE)


delay_state_avg <- delay_state %>% group_by(state_fips) %>%
  mutate(state_total = sum(as.numeric(patient_count_imp)),
         weight = as.numeric(patient_count_imp) / state_total,
         mean = weighted.mean(lag, w = weight)) %>%
  distinct(state_fips, mean)


delay_state_group <- delay_state %>% group_by(state_fips, lag_cat) %>% 
  mutate(lag_sum = sum(as.numeric(patient_count_imp))) %>% 
  distinct(state_fips, lag_cat, lag_sum) %>% 
  ungroup() %>% 
  group_by(state_fips) %>% 
  mutate(total_claims = sum(lag_sum),
         proportion = lag_sum / total_claims)

delay_state_group$state_fips_char <- as.character(delay_state_group$state_fips)
delay_state_group <- delay_state_group %>% mutate(state_fips_char = ifelse(nchar(state_fips_char) == 1, 
                                                      paste0('0', state_fips_char, ''), 
                                                      state_fips_char))

delay_state_group$state_fips_char <- ifelse(length(delay_state_group$state_fips_char) == 1, 
                                            paste0('0', delay_state_group$state_fips_char), 
                                            delay_state_group$state_fips_char)
  
ggplot(delay_state_group, aes(as.character(state_fips_char), proportion, fill= forcats::fct_rev(lag_cat))) + 
  geom_bar(stat="summary") + scale_fill_brewer(palette = "RdYlGn", direction = 1) +
  theme_minimal() +
  labs(title = "Claim Delay by State",
       x = "State",
       y = "Proportion",
       fill = "Delay") 




delay_county <- read.csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1CgyyYhIWFVG2aMITI-HBpKamvwwpQ_ED/CHC_data/date_lag/2024-06-13/date_lag_county_all_clm.csv')

delay_county$patient_count_imp <- delay_county$patient_count
delay_county$patient_count_imp[delay_county$patient_count == '<=5'] <- 
  sample(1:5, length(delay_county$patient_count_imp[delay_county$patient_count == '<=5']), replace= TRUE)

delay_county$lag <- delay_county$chc_rec_claim_date_lag

delay_county_avg <- delay_county %>% group_by(county_fips) %>%
  mutate(county_total = sum(as.numeric(patient_count_imp)),
         weight = as.numeric(patient_count_imp) / county_total,
         mean = weighted.mean(lag, w = weight)) %>%
  distinct(county_fips, mean)



usa_albers_state <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_state.rds'))   # convert to sf




usa_albers_county <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_county.rds')) # convert to sf


usa_albers_county_group_map <- left_join(usa_albers_county_group, delay_county_avg, by = c('county_fips' = 'county_fips'))
map <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_group_map), aes(fill = mean, group = county_fips), color= 'black', linewidth = 0.15) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.4) +
  scale_fill_distiller(palette = "RdYlGn", direction = -1) + ggtitle('Average County Delay') + 
  theme_void() + theme(legend.position = 'right',
                       plot.title = element_text(size = 20, hjust = 0.5),
                       panel.border = element_rect(fill=NA, linewidth = 0.8, color = 'white'),
                       legend.title = element_text(size = 12),
                       legend.text = element_text(size = 12),
                       legend.key.size = unit(0.6, 'cm')) 
map







flu_16_17 <- read.csv('flu_2016-09_2017-08_1b_imputed.csv')

count_dat <- flu_16_17 |> mutate(symp_count = fever + myalgia + cough + 
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
                                 sneezing = ifelse(sneezing == 1, patient_count, 0))


symptom_count <- count_dat |> group_by(symp_count, flu) |> mutate(symptom_sum = sum(patient_count)) |>
  distinct(symp_count, flu, symptom_sum) |> ungroup() |> group_by(flu) |> mutate(sum = sum(symptom_sum),
                                                                                 perc = symptom_sum/sum)


ggplot(symptom_count, aes(fill=as.factor(flu), y=perc, x=symp_count)) + 
  geom_bar(position="dodge", stat="identity") + theme_minimal() + 
  labs(fill = "Influenza?", fill = "Influenza?") + xlab('Number of Symptoms') + 
  ylab('Percent')


ggplot(symptom_count, aes(symp_count, weight = symptom_sum, color = as.factor(flu), fill = as.factor(flu))) +
  geom_density(alpha = 0.1) + theme_minimal() + xlab('Symptom Count') + labs(fill = "Influenza?", color = "Influenza?")

ggplot(diamonds, aes(carat)) +
  geom_density()

overall <- count_dat |> mutate(
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
  sneezing_count_week = sum(sneezing),
  patient_count_week = sum(patient_count)) |>
  distinct(fever_count_week, myalgia_count_week, hypoxemia_count_week,
           short_breath_count_week, cough_count_week, nausea_count_week, chest_pain_count_week,
           sore_throat_count_week, diarrhea_count_week, fatigue_count_week, bronchitis_count_week,
           headache_count_week, congestion_count_week, sneezing_count_week, patient_count_week,
           .keep_all = FALSE) |>
  mutate(
    fever_perc = fever_count_week / patient_count_week,
    myalgia_perc = myalgia_count_week / patient_count_week,
    short_breath_perc = short_breath_count_week / patient_count_week,
    hypoxemia_perc = hypoxemia_count_week / patient_count_week,
    chest_pain_perc = chest_pain_count_week / patient_count_week,
    cough_perc = cough_count_week / patient_count_week,
    bronchitis_perc = bronchitis_count_week / patient_count_week,
    nausea_perc = nausea_count_week / patient_count_week,
    sore_throat_perc = sore_throat_count_week / patient_count_week,
    diarrhea_perc = diarrhea_count_week / patient_count_week,
    fatigue_perc = fatigue_count_week / patient_count_week,
    headache_perc = headache_count_week / patient_count_week,
    congestion_perc = congestion_count_week / patient_count_week,
    sneezing_perc = sneezing_count_week / patient_count_week)

by_flu <- count_dat |> group_by(flu) |> mutate(
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
  sneezing_count_week = sum(sneezing),
  patient_count_week = sum(patient_count),
  asymptomatic = sum(no_symptoms),
  symptom_count = mean(symp_count)) |>
  distinct(flu, fever_count_week, myalgia_count_week, hypoxemia_count_week,
           short_breath_count_week, cough_count_week, nausea_count_week, chest_pain_count_week,
           sore_throat_count_week, diarrhea_count_week, fatigue_count_week, bronchitis_count_week,
           headache_count_week, congestion_count_week, sneezing_count_week, patient_count_week, asymptomatic, symptom_count,
           .keep_all = FALSE) |>
  mutate(
    fever_perc = fever_count_week / patient_count_week,
    myalgia_perc = myalgia_count_week / patient_count_week,
    short_breath_perc = short_breath_count_week / patient_count_week,
    hypoxemia_perc = hypoxemia_count_week / patient_count_week,
    chest_pain_perc = chest_pain_count_week / patient_count_week,
    cough_perc = cough_count_week / patient_count_week,
    bronchitis_perc = bronchitis_count_week / patient_count_week,
    nausea_perc = nausea_count_week / patient_count_week,
    sore_throat_perc = sore_throat_count_week / patient_count_week,
    diarrhea_perc = diarrhea_count_week / patient_count_week,
    fatigue_perc = fatigue_count_week / patient_count_week,
    headache_perc = headache_count_week / patient_count_week,
    congestion_perc = congestion_count_week / patient_count_week,
    sneezing_perc = sneezing_count_week / patient_count_week,
    asymptomatic_perc = asymptomatic / patient_count_week)


by_symptom <- flu_16_17 |> filter(flu == 1) |> group_by(fever, myalgia, cough, sore_throat, short_breath,
                                                        hypoxemia, chest_pain, bronchitis, nausea_vom, diarrhea,
                                                        fatigue, headache, congestion, sneezing) |>
  mutate(sum = sum(patient_count)) |>
  distinct(sum, fever, myalgia, cough, sore_throat, short_breath,
           hypoxemia, chest_pain, bronchitis, nausea_vom, diarrhea,
           fatigue, headache, congestion, sneezing)


expand <- flu_16_17 |> filter(flu == 1) |> uncount(patient_count)


corr <- cor(expand[ , c(4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)])

library(corrplot)
corrplot(corr)

co_occur <- flu_16_17 |> filter(flu == 1) |> uncount(patient_count) |> select(-c(flu, county_fips, year_week))


dat <- co_occur |> filter(i == 1) |> mutate(symp_count = fever + myalgia + cough + 
                                              sore_throat + short_breath + hypoxemia + 
                                              chest_pain + bronchitis + nausea_vom + 
                                              diarrhea + fatigue + headache + congestion + sneezing,
                                            only = ifelse(symp_count == 1, 1, 0))




occur <- function(i) {
  dat <- co_occur |> filter(co_occur[, i] == 1) |> mutate(symp_count = fever + myalgia + cough + 
                                                            sore_throat + short_breath + hypoxemia + 
                                                            chest_pain + bronchitis + nausea_vom + 
                                                            diarrhea + fatigue + headache + congestion + sneezing,
                                                          only = ifelse(symp_count == 1, 1, 0)) |> 
    summarise_all(mean)
  return(dat)
  
}

symptom_co_occurance <- data.frame(sapply(1:14, occur))
colnames(symptom_co_occurance) <- colnames(co_occur)

mean_dat <- data.frame(colMeans(dat))

for (i in colnames(co_occur)) {
  print(i)
}



flu_16_17 <- read.csv('flu_2016-09_2017-08_1b_imputed.csv')
flu_17_18 <- read.csv('flu_2017-09_2018-08_1b_imputed.csv')
flu_18_19 <- read.csv('flu_2018-09_2019-08_1b_imputed.csv')
flu_19_20 <- read.csv('flu_2019-09_2020-08_1b_imputed.csv')

flu_append <- rbind(flu_16_17, flu_17_18, flu_18_19, flu_19_20)

count_dat <- flu_append |> mutate(symp_count = fever + myalgia + cough + 
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
                                  sneezing = ifelse(sneezing == 1, patient_count, 0))

flu_no_symptoms_count <- count_dat |> filter(flu == 1) |> group_by(year_week) |>
  mutate(no_symptoms_sum = sum(no_symptoms),
         year_week_date = as.Date(paste(year_week, 1, sep = '-'), format = '%Y-%U-%u')) |>
  distinct(year_week_date, no_symptoms_sum, .keep_all = FALSE)

flu_symptoms_count <- count_dat |> filter(flu == 1) |> group_by(county_fips, symp_count) |>
  mutate(symptoms_sum = sum(patient_count)) |>
  distinct(county_fips, symp_count, symptoms_sum, .keep_all = FALSE) |>
  ungroup() |> group_by(county_fips) |> mutate(week_tot = sum(symptoms_sum),
                                               percent = symptoms_sum / week_tot) |>
  filter(symp_count == 0)

ggplot(flu_no_symptoms_count) + 
  geom_line(aes(x = as.Date(year_week_date), y = no_symptoms_sum), color = '#00BA42') +
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('No Symptoms') +
  theme(plot.title = element_text(hjust = 0.5)) + theme_minimal()

ggplot(flu_symptoms_count) + 
  geom_line(aes(x = as.Date(year_week_date), y = percent), color = '#00BA42') +
  theme_minimal() + xlab('Date') + ylab('Percent') + ggtitle('No Symptoms Percent') +
  theme(plot.title = element_text(hjust = 0.5)) + theme_minimal() + ylim(0, 1)

ggplot(flu_symptoms_count, aes(x = percent)) +
  geom_density(fill = '#00BA42', color = '#00BA42', alpha = 0.5) + theme_minimal() +
  xlab('Percent') + ylab('Density') + ggtitle('No Symptoms Percent Spatial Distribution')



usa_albers_county_group <- left_join(usa_albers_county_group, flu_symptoms_count[, c(1, 5)],
                                     by = c('county_fips' = 'county_fips'))
map <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_group), aes(fill = percent, group = county_fips), color= 'black', linewidth = 0.15) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.4) +
  scale_fill_viridis_c('Percent', direction = -1) + ggtitle('Percent of Asymptomatic Flu Cases') + 
  theme_void() + theme(legend.position = 'right',
                       plot.title = element_text(size = 20, hjust = 0.5),
                       panel.border = element_rect(fill=NA, linewidth = 0.8, color = 'white'),
                       legend.title = element_text(size = 12),
                       legend.text = element_text(size = 12),
                       legend.key.size = unit(0.6, 'cm')) 
map

