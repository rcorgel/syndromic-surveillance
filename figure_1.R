################################################################################
# File Name: figure_1                                                          #
#                                                                              #
# Purpose:   Create figure 1 for the manuscript.                               #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Create figure 1                                                #
#               a. Create time series plots                                    #
#               b. Create map plots                                            #
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

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

######################
# 2. CREATE FIGURE 1 #
######################

###############################
# A. CREATE TIME SERIES PLOTS #
###############################

## FLU DATA ##

# Load flu data for time series
flu_2016 <- readRDS('./tmp/flu_filt_symp_p2_2016.rds')
flu_2017 <- readRDS('./tmp/flu_filt_symp_p2_2017.rds')
flu_2018 <- readRDS('./tmp/flu_filt_symp_p2_2018.rds')
flu_2019 <- readRDS('./tmp/flu_filt_symp_p2_2019.rds')

# Combine data
flu <- rbind(flu_2016, flu_2017, flu_2018, flu_2019)
remove(flu_2016, flu_2017, flu_2018, flu_2019)

# Collapse to the national level
flu_nat <- flu |> group_by(week_date) |>
  mutate(flu_sum = sum(flu * patient_count_imp),
         fever_sum = sum(fever * patient_count_imp),
         cough_sum = sum(cough * patient_count_imp),
         sore_throat_sum = sum(sore_throat * patient_count_imp),
         myalgia_sum = sum(myalgia * patient_count_imp),
         hypoxemia_sum = sum(hypoxemia * patient_count_imp),
         short_breath_sum = sum(short_breath * patient_count_imp),
         bronchitis_sum = sum(bronchitis * patient_count_imp),
         chest_pain_sum = sum(chest_pain * patient_count_imp),
         nausea_vom_sum = sum(nausea_vom * patient_count_imp),
         fatigue_sum = sum(fatigue * patient_count_imp),
         diarrhea_sum = sum(diarrhea * patient_count_imp),
         headache_sum = sum(headache * patient_count_imp),
         congestion_sum = sum(congestion * patient_count_imp),
         sneezing_sum = sum(sneezing * patient_count_imp)) |>
  distinct(week_date, flu_sum, fever_sum, cough_sum, sore_throat_sum,
           myalgia_sum, hypoxemia_sum, short_breath_sum, bronchitis_sum,
           chest_pain_sum, nausea_vom_sum, fatigue_sum, diarrhea_sum,
           headache_sum, congestion_sum, sneezing_sum)

# Re-collapse data to the national-week level
flu_week <- flu_week |>
  filter(week_date > as.Date('2016-08-31')) |>
  filter(week_date < as.Date('2020-08-30')) |>
  group_by(week_date) |>
  mutate(flu_sum = sum(flu_sum),
         fever_sum = sum(fever_sum),
         cough_sum = sum(cough_sum),
         sore_throat_sum = sum(sore_throat_sum),
         myalgia_sum = sum(myalgia_sum),
         hypoxemia_sum = sum(hypoxemia_sum),
         short_breath_sum = sum(short_breath_sum),
         bronchitis_sum = sum(bronchitis_sum),
         chest_pain_sum = sum(chest_pain_sum),
         nausea_vom_sum = sum(nausea_vom_sum),
         fatigue_sum = sum(fatigue_sum),
         diarrhea_sum = sum(diarrhea_sum),
         headache_sum = sum(headache_sum),
         congestion_sum = sum(congestion_sum),
         sneezing_sum = sum(sneezing_sum)) |>
  distinct(week_date, flu_sum, fever_sum, cough_sum, sore_throat_sum,
           myalgia_sum, hypoxemia_sum, short_breath_sum, bronchitis_sum,
           chest_pain_sum, nausea_vom_sum, fatigue_sum, diarrhea_sum,
           headache_sum, congestion_sum, sneezing_sum)

## RSV DATA ##

# Load rsv data for time series
rsv_dat <- readRDS('./tmp/rsv_p2_dat_proc_symp.rds')

# Collapse data to national-week level
rsv_week <- rsv_dat |> 
  filter(week_date > as.Date('2016-08-31')) |>
  filter(week_date < as.Date('2020-08-30')) |>
  group_by(week_date) |>
  mutate(rsv_sum = sum(rsv * patient_count_imp),
         loss_appetite_sum = sum(loss_appetite * patient_count_imp)) |>
  distinct(week_date, rsv_sum, loss_appetite_sum)

## PEAK CALCULATION ##

# Calculate season peaks for each disease and symptom
flu_week$flu_peak <- peaks(flu_week$flu_sum, span=31, strict=TRUE)
flu_peaks <- flu_week |> filter(flu_peak == TRUE)

flu_week$fever_peak <- peaks(flu_week$fever_sum, span=31, strict=TRUE)
fever_peaks <- flu_week |> filter(fever_peak == TRUE)

flu_week$cough_peak <- peaks(flu_week$cough_sum, span=31, strict=TRUE)
cough_peaks <- flu_week |> filter(cough_peak == TRUE)

flu_week$sore_peak <- peaks(flu_week$sore_throat_sum, span=31, strict=TRUE)
sore_peaks <- flu_week |> filter(sore_peak == TRUE)

rsv_week$rsv_peak <- peaks(rsv_week$rsv_sum, span=31, strict=TRUE)
rsv_peaks <- rsv_week |> filter(rsv_peak == TRUE)

## CREATE TIME SERIES PLOTS ##

# Set the plot theme
line_theme <- theme(legend.position = "none",
                    axis.text = element_text(size=14),
                    axis.title = element_text(size=16),
                    strip.background = element_blank(),
                    plot.title = element_text(size=16))

# Make basic line plots
line_flu <- ggplot(flu_week, 
       aes(x = week_date, y = flu_sum)) +
  geom_rect(aes(xmin=ymd('2018-09-01'),
                xmax = ymd('2019-08-31'),
                ymin = 0,
                ymax = 225000), , fill="#FFFFFF00", color = 'black', linewidth = 0.6) +
  geom_area(fill = "#d7642c", alpha = 0.4) +  
  geom_line(color = "#d7642c", linewidth = 0.5) +
  geom_segment(aes(x = flu_peaks$week_date[1], y = 0, xend = flu_peaks$week_date[1], yend = 225000), 
             color = "#d7642c", linetype = 2, linewidth = 1.1) +
  geom_segment(aes(x = flu_peaks$week_date[2], y = 0, xend = flu_peaks$week_date[2], yend = 225000), 
               color = "#d7642c", linetype = 2, linewidth = 1.1) +
  geom_segment(aes(x = flu_peaks$week_date[3], y = 0, xend = flu_peaks$week_date[3], yend = 225000), 
               color = "#d7642c", linetype = 2, linewidth = 1.1) +
  geom_segment(aes(x = flu_peaks$week_date[4], y = 0, xend = flu_peaks$week_date[4], yend = 225000), 
               color = "#d7642c", linetype = 2, linewidth = 1.1) +
  theme_minimal() + line_theme +
  labs(title = "Influenza National Healthcare Visits") +
  xlab('Week') + ylab('Visits') +
  scale_y_continuous(labels = comma)

line_rsv <- ggplot(rsv_week, 
       aes(x = week_date, y = rsv_sum)) +
  geom_rect(aes(xmin=ymd('2018-09-01'),
                xmax = ymd('2019-08-31'),
                ymin = 0,
                ymax = 25000), , fill="#FFFFFF00", color = 'black', linewidth = 0.6) +
  geom_area(fill = "#41afaa", alpha = 0.4) +  
  geom_line(color = "#41afaa", linewidth = 0.5) +
  geom_segment(aes(x = rsv_peaks$week_date[1], y = 0, xend = rsv_peaks$week_date[1], yend = 25000), 
               color = "#41afaa", linetype = 2, linewidth = 1.1) +
  geom_segment(aes(x = rsv_peaks$week_date[2], y = 0, xend = rsv_peaks$week_date[2], yend = 25000), 
               color = "#41afaa", linetype = 2, linewidth = 1.1) +
  geom_segment(aes(x = rsv_peaks$week_date[3], y = 0, xend = rsv_peaks$week_date[3], yend = 25000), 
               color = "#41afaa", linetype = 2, linewidth = 1.1) +
  geom_segment(aes(x = rsv_peaks$week_date[4], y = 0, xend = rsv_peaks$week_date[4], yend = 25000), 
               color = "#41afaa", linetype = 2, linewidth = 1.1) +
  theme_minimal() + line_theme +
  labs(title = "RSV National Healthcare Visits") +
  xlab('Week') + ylab('Visits') +
  scale_y_continuous(labels = comma)

line_fever <- ggplot(flu_week, 
       aes(x = week_date, y = fever_sum)) +
  geom_rect(aes(xmin=ymd('2018-09-01'),
                xmax = ymd('2019-08-31'),
                ymin = 0,
                ymax = 360000), , fill="#FFFFFF00", color = 'black', linewidth = 0.6) +
  geom_area(fill = "#af4b91", alpha = 0.4) +  
  geom_line(color = "#af4b91", linewidth = 0.5) +
  geom_segment(aes(x = fever_peaks$week_date[1], y = 0, xend = fever_peaks$week_date[1], yend = 360000), 
               color = "#af4b91", linetype = 2, linewidth = 1.1) +
  geom_segment(aes(x = fever_peaks$week_date[2], y = 0, xend = fever_peaks$week_date[2], yend = 360000), 
               color = "#af4b91", linetype = 2, linewidth = 1.1) +
  geom_segment(aes(x = fever_peaks$week_date[3], y = 0, xend = fever_peaks$week_date[3], yend = 360000), 
               color = "#af4b91", linetype = 2, linewidth = 1.1) +
  geom_segment(aes(x = fever_peaks$week_date[4], y = 0, xend = fever_peaks$week_date[4], yend = 360000), 
               color = "#af4b91", linetype = 2, linewidth = 1.1) +
  theme_minimal() + line_theme +
  labs(title = "Fever National Healthcare Visits") +
  xlab('Week') + ylab('Visits') +
  scale_y_continuous(labels = comma)

line_cough <- ggplot(flu_week, 
       aes(x = week_date, y = cough_sum)) +
  geom_rect(aes(xmin=ymd('2018-09-01'),
                xmax = ymd('2019-08-31'),
                ymin = 0,
                ymax = 500000), fill="#FFFFFF00", color = 'black', linewidth = 0.6) +
  geom_area(fill = "#e6a532", alpha = 0.4) +  
  geom_line(color = "#e6a532", linewidth = 0.5) +
  geom_segment(aes(x = cough_peaks$week_date[1], y = 0, xend = cough_peaks$week_date[1], yend = 500000), 
               color = "#e6a532", linetype = 2, linewidth = 1.1) +
  geom_segment(aes(x = cough_peaks$week_date[2], y = 0, xend = cough_peaks$week_date[2], yend = 500000), 
               color = "#e6a532", linetype = 2, linewidth = 1.1) +
  geom_segment(aes(x = cough_peaks$week_date[3], y = 0, xend = cough_peaks$week_date[3], yend = 500000), 
               color = "#e6a532", linetype = 2, linewidth = 1.1) +
  geom_segment(aes(x = cough_peaks$week_date[4], y = 0, xend = cough_peaks$week_date[4], yend = 500000), 
               color = "#e6a532", linetype = 2, linewidth = 1.1) +
  theme_minimal() + line_theme +
  labs(title = "Cough National Healthcare Visits") +
  xlab('Week') + ylab('Visits') +
  scale_y_continuous(labels = comma)

line_sore <- ggplot(flu_week, 
       aes(x = week_date, y = sore_throat_sum)) +
  geom_rect(aes(xmin=ymd('2018-09-01'),
                xmax = ymd('2019-08-31'),
                ymin = 0,
                ymax = 350000), , fill="#FFFFFF00", color = 'gray', linewidth = 0.7) +
  geom_area(fill = "#e6a532", alpha = 0.4) +  
  geom_line(color = "#e6a532", linewidth = 0.5) +
  geom_segment(aes(x = sore_peaks$week_date[1], y = 0, xend = sore_peaks$week_date[1], yend = 350000), 
               color = "#e6a532", linetype = 2, linewidth = 1.1) +
  geom_segment(aes(x = sore_peaks$week_date[2], y = 0, xend = sore_peaks$week_date[2], yend = 350000), 
               color = "#e6a532", linetype = 2, linewidth = 1.1) +
  geom_segment(aes(x = sore_peaks$week_date[3], y = 0, xend = sore_peaks$week_date[3], yend = 350000), 
               color = "#e6a532", linetype = 2, linewidth = 1.1) +
  geom_segment(aes(x = sore_peaks$week_date[4], y = 0, xend = sore_peaks$week_date[4], yend = 350000), 
               color = "#e6a532", linetype = 2, linewidth = 1.1) +
  theme_minimal() + line_theme +
  labs(title = "Sore Throat National Healthcare Visits") +
  xlab('Week') + ylab('Visits') +
  scale_y_continuous(labels = comma)

#######################
# B. CREATE MAP PLOTS #
#######################

## FLU ##

# Load flu data
flu_dat_p1 <- readRDS('./tmp/flu_p1_dat_proc_symp.rds')

# For a single season, collapse to county level
flu_county <- flu_dat_p1 |> 
  filter(month_date > as.Date('2018-08-31')) |>
  filter(month_date < as.Date('2019-09-01')) |>
  group_by(county_fips) |>
  mutate(flu_sum = sum(flu * patient_count_imp),
         fever_sum = sum(fever * patient_count_imp),
         cough_sum = sum(cough * patient_count_imp),
         sore_throat_sum = sum(sore_throat * patient_count_imp),
         myalgia_sum = sum(myalgia * patient_count_imp),
         hypoxemia_sum = sum(hypoxemia * patient_count_imp),
         short_breath_sum = sum(short_breath * patient_count_imp),
         bronchitis_sum = sum(bronchitis * patient_count_imp),
         chest_pain_sum = sum(chest_pain * patient_count_imp),
         nausea_vom_sum = sum(nausea_vom * patient_count_imp),
         fatigue_sum = sum(fatigue * patient_count_imp),
         diarrhea_sum = sum(diarrhea * patient_count_imp),
         headache_sum = sum(headache * patient_count_imp),
         congestion_sum = sum(congestion * patient_count_imp),
         sneezing_sum = sum(sneezing * patient_count_imp)) |>
  distinct(county_fips, flu_sum, fever_sum, cough_sum, sore_throat_sum,
           myalgia_sum, hypoxemia_sum, short_breath_sum, bronchitis_sum,
           chest_pain_sum, nausea_vom_sum, fatigue_sum, diarrhea_sum,
           headache_sum, congestion_sum, sneezing_sum)

## RSV ##

# Load rsv data
rsv_dat_p1 <- readRDS('./tmp/rsv_p1_dat_proc_symp.rds')

# For a single season, collapse to county level
rsv_county <- rsv_dat_p1 |> 
  filter(month_date > as.Date('2018-08-31')) |>
  filter(month_date < as.Date('2019-09-01')) |>
  group_by(county_fips) |>
  mutate(rsv_sum = sum(rsv * patient_count_imp),
         loss_appetite_sum = sum(loss_appetite * patient_count_imp)) |>
  distinct(county_fips, rsv_sum, loss_appetite_sum)

## ALL CAUSE ##

# Load all cause data
all_cause_season <- read_csv('raw/county_season_ac_v2.csv')

# Impute all cause data
all_cause_season$patient_count_imp <- all_cause_season$all_cause
all_cause_season$patient_count_imp[all_cause_season$patient_count_imp == '<=5'] <- 
  sample(1:5, length(all_cause_season$patient_count_imp[all_cause_season$patient_count_imp == '<=5']), 
         replace= TRUE)
# Convert imputed count to numeric
all_cause_season$patient_count_imp <- as.numeric(all_cause_season$patient_count_imp)

# Convert week date to date and add state
all_cause_season <- all_cause_season |> 
  mutate(state_fips = substr(county_fips, 1, 2)) |>
  filter(!(state_fips %in% c("72","99",""))) |>
  filter(season == '2018-2019')

## DISAGGREGATE COUNTY GROUPS ##

# Load county-group cross walk
county_group <- read.csv('./raw/county_fips_grp_pop_wts.csv')

# Merge on rsv/symptom data
county_prev <- left_join(county_group, rsv_county, 
                         by = c('county_fips_grp' = 'county_fips'))

# Merge on flu/symptom data
county_prev <- left_join(county_prev, flu_county, 
                         by = c('county_fips_grp' = 'county_fips'))

# Merge on all cause data
county_prev <- left_join(county_prev, all_cause_season, 
                         by = c('county_fips_grp' = 'county_fips'))

# Disaggregate counties
county_prev <- county_prev |>
  mutate(flu = flu_sum * pop_wt,
         rsv = rsv_sum * pop_wt,
         fever = fever_sum * pop_wt,
         cough = cough_sum * pop_wt,
         sore_throat = sore_throat_sum * pop_wt,
         myalgia  = myalgia_sum * pop_wt,
         hypoxemia = hypoxemia_sum * pop_wt,
         short_breath = short_breath_sum * pop_wt,
         bronchitis = bronchitis_sum * pop_wt,
         chest_pain = chest_pain_sum * pop_wt,
         nausea_vom = nausea_vom_sum * pop_wt,
         fatigue = fatigue_sum * pop_wt,
         diarrhea = diarrhea_sum * pop_wt,
         headache = headache_sum * pop_wt,
         congestion = congestion_sum * pop_wt,
         sneezing = sneezing_sum * pop_wt,
         loss_appetite = loss_appetite_sum * pop_wt,
         all = patient_count_imp * pop_wt) |>
  mutate(flu = flu / all,
         rsv = rsv / all,
         fever = fever / all,
         cough = cough / all,
         sore_throat = sore_throat / all,
         myalgia  = myalgia / all,
         hypoxemia = hypoxemia / all,
         short_breath = short_breath / all,
         bronchitis = bronchitis / all,
         chest_pain = chest_pain / all,
         nausea_vom = nausea_vom / all,
         fatigue = fatigue / all,
         diarrhea = diarrhea / all,
         headache = headache / all,
         congestion = congestion / all,
         sneezing = sneezing / all,
         loss_appetite = loss_appetite / all)

## FILTER OUT BASED ON POPULATION, CLAIM VOLUME, AND OUTLIERS ##

# Filter out low population and low claim volume counties
county_prev_filt <- county_prev |> 
  filter(pop_2020 > 10000) |> 
  filter(all > 5000) |>
  mutate(county_fips_str = ifelse(nchar(as.character(county_fips)) == 4,
                                  paste0('0', as.character(county_fips)),
                                  as.character(county_fips))) |>
  dplyr::select(c(county_fips_str, flu, fever, cough, sore_throat,
                  myalgia, hypoxemia, short_breath, bronchitis,
                  chest_pain, nausea_vom, fatigue, diarrhea,
                  headache, congestion, sneezing, rsv, loss_appetite))

# Calculate outliers

# Reshape the data to long format
county_prev_filt_long <- county_prev_filt |>
  pivot_longer(!c(county_fips_str), names_to = "condition", values_to = "prev")

# Hide outliers
county_prev_filt_long <- county_prev_filt_long |>
  group_by(condition) |>
  mutate(mean = mean(prev),
         sd = sd(prev),
         z = (prev - mean) / sd,
         prev_filt = ifelse(z > 3 | z < -3, NA, prev))

# Reshape to wide data
county_prev_filt_wide <- county_prev_filt_long |>
  select(-c(mean, sd, z, prev)) |>
  pivot_wider(names_from = condition, values_from = prev_filt)
  
## CREATE MAPS ##

# Load map data
usa_albers_state <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_state.rds'))   # convert to sf
usa_albers_county <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_county.rds'))   # convert to sf

# Merge on prev data
usa_albers_county <- left_join(usa_albers_county, county_prev_filt_wide, by = c('GEOID' = 'county_fips_str'))

# Set the map theme
map_theme <- theme(legend.position = "right",
                   axis.title = element_text(size=16),
                   strip.background = element_blank(),
                   plot.title = element_blank(),
                   panel.border = element_rect(fill=NA, linewidth = 0.8, color = 'white'),
                   legend.title = element_text(size = 16),
                   legend.text = element_text(size = 12),
                   legend.key.height = unit(0.7, 'cm'),
                   legend.key.width = unit(0.4, "cm"))

flu_map <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county), aes(fill = flu, group = GEOID), color= 'black', linewidth = 0.01) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.10) +
  scale_fill_gradient('2018-19\nInfluenza\nPrevalence \n', low = "white", high = "#d7642c", na.value = "grey",
                      breaks = c(0, 0.015, 0.03, 0.045), 
                      limits = c(0, max(usa_albers_county$flu)),
                      labels = c("0.0%","1.5%","3.0%", "4.5%")) + 
  theme_void() + map_theme 

rsv_map <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county), aes(fill = rsv, group = GEOID), color= 'black', linewidth = 0.01) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.10) +
  scale_fill_gradient('2018-19\nRSV\nPrevalence \n', low = "white", high = "#41afaa", na.value = "grey",
                      breaks = c(0, 0.002, 0.004, 0.006, 0.008), 
                      limits = c(0, max(usa_albers_county$rsv)),
                      labels = c("0.0%","0.2%","0.4%", "0.6%", "0.8%")) + 
  theme_void() + map_theme 

fever_map <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county), aes(fill = fever, group = GEOID), color= 'black', linewidth = 0.01) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.10) +  
  scale_fill_gradient('2018-19\nFever\nPrevalence \n', low = "white", high = "#af4b91", na.value = "grey",
                      breaks = c(0, 0.02, 0.04, 0.06, 0.08), 
                      limits = c(0, max(usa_albers_county$fever)),
                      labels = c("0.0%","2.0%","4.0%", "6.0%", "8.0%")) + 
  theme_void() + map_theme 

cough_map <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county), aes(fill = cough, group = GEOID), color= 'black', linewidth = 0.01) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.10) +
  scale_fill_gradient('2018-19\nCough\nPrevalence \n', low = "white", high = "#e6a532", na.value = "grey",
                      breaks = c(0, 0.03, 0.06, 0.09), 
                      limits = c(0, max(usa_albers_county$cough)),
                      labels = c("0.0%","3.0%","6.0%", "9.0%")) + 
  theme_void() + map_theme 

sore_map <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county), aes(fill = sore_throat, group = GEOID), color= 'black', linewidth = 0.01) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.10) +
  scale_fill_gradient('Sore Throat \nPrevalence \n', low = "white", high = "#e6a532",
                      breaks = c(0, 0.021, 0.041, 0.062), limits = c(0, 0.063), na.value = "grey",
                      labels = scales::label_number(accuracy = 0.001)) + 
  theme_void() + map_theme 

## CONBINE SUBFIGURES ##

# Create figure, with labels
figure_1 <- cowplot::plot_grid(line_flu, line_rsv, line_fever, line_cough,
                               flu_map, rsv_map, fever_map, cough_map,
                               nrow = 2,
                               rel_heights = c(1, 1),
                               labels = c("a", "b", "c", "d",
                                          "", "", "", ""),
                               label_size = 20)

# Save figure 
ggsave('./figs/figure_1.jpg', plot = figure_1, height = 5.5, width = 20)

################################################################################
################################################################################
