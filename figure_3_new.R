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
time_model <- readRDS("./tmp/test_pos_model.rds")

# Load part 2 data
flu_2016 <- readRDS('./tmp/flu_filt_symp_p2_2018.rds')

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


# CHange var to urban_south
flu_2016$urban_south <- ifelse(flu_2016$urban_code_binary == 1 & flu_2016$south == 1, 'Urban South', '')
flu_2016$urban_south <- ifelse(flu_2016$urban_code_binary == 0 & flu_2016$south == 1, 'Rural South', flu_2016$urban_south)
flu_2016$urban_south <- ifelse(flu_2016$urban_code_binary == 1 & flu_2016$south == 0, 'Urban North', flu_2016$urban_south)
flu_2016$urban_south <- ifelse(flu_2016$urban_code_binary == 0 & flu_2016$south == 0, 'Rural North', flu_2016$urban_south)





# Add nrevss flu data
nrevss_nat <- readRDS('./tmp/nrevss_national.rds')
nrevss_nat$week_date_ALT <- nrevss_nat$week_date + 0
nrevss_state <- readRDS('./tmp/nrevss_state.rds')
nrevss_state$week_date_ALT <- nrevss_state$week_date + 0
flu_2016 <- left_join(flu_2016, nrevss_nat[, c(2, 9)], by = c('year_week_dt' = 'week_date_ALT'))
flu_2016 <- flu_2016 |> rename('nrevss_nat' = 'value')
flu_2016 <- left_join(flu_2016, nrevss_state[, c(2, 9, 4)], 
                 by = c('year_week_dt' = 'week_date_ALT',
                        'state_fips' = 'state_fips'))
flu_2016 <- flu_2016 |> rename('nrevss_state' = 'value')

flu_2016 <- flu_2016 |> mutate(nrevss_state = ifelse(is.na(nrevss_state), nrevss_nat, nrevss_state))

flu_2016_filt <- flu_2016 |> dplyr::filter(!is.na(nrevss_state)) |>
  dplyr::filter(!is.na(nrevss_nat))



for (i in unique(flu_2016_filt$year_week_dt)) {
  # Print the progress and set season data
  print(i)
  flu_week <- flu_2016_filt |> dplyr::filter(year_week_dt == i)
  
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
saveRDS(flu_pred_dat, './tmp/flu_pred_dat_county_p2_time_dis_week_nrevss_state_urban_2018.rds')

# Remove data frame list
remove(flu_nat_filt_county_week)
  
  
  flu_pred_dat <- readRDS('./tmp/flu_pred_dat_county_p2_time_dis_week_nrevss_state_urban.rds')
  
  
  
  
  
  # Set cutoff value for prediction
  flu_pred_dat$pred_count <- ifelse(flu_pred_dat$pred > 0.50, 1, 0) * flu_pred_dat$patient_count
  
  
  
  
flu_pred_dat_pos <- flu_pred_dat |> dplyr::filter(flu == 1)
  

flu_pred_dat_pos$not <- ifelse(flu_pred_dat_pos$fever == 0 & flu_pred_dat_pos$cough == 0 & flu_pred_dat_pos$sore_throat == 0, 1, 0)

flu_test <- flu_pred_dat_pos |> 
  group_by(year_week_dt, county_fips) |>
  mutate(fever_sum = sum(patient_count * fever),
         flu_sum = sum(patient_count * flu),
         cough_sum = sum(patient_count * cough),
         sore_throat_sum = sum(patient_count * sore_throat),
         not_sum = sum(patient_count * not)) |>
  distinct(year_week_dt, county_fips, fever_sum, cough_sum, sore_throat_sum, not_sum, flu_sum)


ggplot(flu_test[flu_test$county_fips == '39035',]) + 
  geom_line(aes(x = year_week_dt, y = fever_sum), color = 'blue') +
  geom_line(aes(x = year_week_dt, y = cough_sum), color = 'red') +
  geom_line(aes(x = year_week_dt, y = sore_throat_sum), color = 'green') +
  geom_line(aes(x = year_week_dt, y = not_sum), color = 'orange') +
  geom_line(aes(x = year_week_dt, y = flu_sum), color = 'black') 
  
  
  
  # COUNTY LEVEL
  flu_county <- flu_pred_dat |> 
    group_by(year_week_dt, county_fips) |>
    mutate(pred_sum = sum(pred_count, na.rm = T),
           flu_sum = sum(patient_count * flu)) |>
    distinct(year_week_dt, county_fips, pred_sum, flu_sum)
  
  
  percent_theme <- theme(legend.position = "inside",
                         legend.position.inside = c(0.78, 0.88),
                         axis.text = element_text(size=16),
                         axis.title.y = element_text(size=18),
                         axis.title.x = element_blank(),
                         strip.background = element_blank(),
                         plot.title = element_text(size=20),
                         axis.text.x = element_text(angle = 35, vjust = 1, 
                                                    hjust=1, color = 'black'),
                         legend.title = element_text(size = 16),
                         legend.text = element_text(size = 14),
                         legend.box.background = element_rect(colour = "black", 
                                                              fill = 'white'))
  
  # Set colors for influenza
  flu_colors <- c(
    "Tested Cases" = "#CC6594",
    "Syndromic Predicted Cases" = "#3cbb75ff")
  
  
pred <- flu_county[flu_county$county_fips == '36061', c(1, 3)]
pred <- pred |> rename('flu_sum' = 'pred_sum')
pred$Signal <- "Syndromic Predicted Cases"
obs <- flu_county[flu_county$county_fips == '36061', c(1, 4)]
obs$Signal <- "Tested Cases"

nyc_dat <- rbind(pred, obs)

  
  
  # Create flight
ggplot(nyc_dat, aes(x = year_week_dt, y = flu_sum)) +
    geom_line(aes(group = Signal, color = Signal), linewidth = 1.75, alpha = 0.9) +  
    theme_minimal() + percent_theme +
    labs(title = "New York County, NY") +
    ylab('Number of Visits') + xlab('') + expand_limits(y = c(0, 0.5)) +
    scale_color_manual('Signal', values = flu_colors)
  
  

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
  
  ggplot(flu_county_peak[flu_county_peak$county_fips == '36061',]) + 
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
    distinct(week_pred, week_obs, week_obs_num, week_pred_num, week_obs_num_adj, week_pred_num_adj, count_sum)
  
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
  
  forrest_theme <- theme(legend.position =  "bottom",
                         legend.position.inside = c(0.873, 0.237),
                         axis.text = element_text(size=14, color = 'black'),
                         axis.title = element_text(size=16),
                         strip.background = element_blank(),
                         legend.title = element_text(size=16),
                         legend.text = element_text(size=14),
                         plot.title = element_text(size=18, hjust = 0),
                         legend.box="horizontal",
                         legend.box.background = element_rect(colour = "white", 
                                                              fill = 'white'),
                         legend.key.width = unit(1.5, "cm"))
  
  
  peak_comp_plot <- ggplot(county_peaks_both_col[!is.na(county_peaks_both_col$week_obs),], aes(week_obs, week_pred, fill= count_sum)) + 
    geom_tile(color = 'black') + 
    scale_fill_gradient2('Number of Counties\n', low = 'white', high = '#3cbb75ff') +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black",
                linewidth = 1) +
    theme_minimal() + 
    labs(title = 'Peak Time Comparison',
         x = "Claims Peak Week",
         y = "Syndromic Peak Week") + forrest_theme +
    scale_x_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2016-09-01'), as.Date('2017-05-27'))) +
    scale_y_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2016-09-01'), as.Date('2017-05-27')))
  


peak_comp_plot 


usa_albers_county <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_county.rds'))   # convert to sf


county_tomp <- usa_albers_county  |> dplyr::filter(COUNTYFP == '109' | COUNTYFP == '107') |>
  dplyr::filter(STATEFP == '36')


ggplot() +
  geom_sf(data = st_as_sf(county_tomp), fill = 'white', color= 'black', linewidth = 1) + theme_void()


#county_peaks_pred_filt <- county_peaks_pred |> dplyr::filter(week_pred_num < 16)
usa_albers_county_2018 <- left_join(usa_albers_county, 
                                    county_peaks_pred, 
                                    by = c('GEOID' = 'county_fips'))


county_peaks_both_filt$diff <- county_peaks_both_filt$week_pred - county_peaks_both_filt$week_obs
usa_albers_county_map <- left_join(usa_albers_county, 
                                   county_peaks_both_filt, 
                                    by = c('GEOID' = 'county_fips'))



colors <- c('0-4' = '#b9d4ec', '5-12' = '#8bb7e0', '13-17' = '#5d9ad3', 
            '18-49' =  '#347DC1', '50-64' = '#285f93', '65+' = '#1b4164')


usa_albers_county_map$month <- month(usa_albers_county_map$week_pred)
usa_albers_county_map$Month <- factor(usa_albers_county_map$month, levels = c(1, 2, 3, 4, 5), labels = c('Jan', 'Feb', 'Mar', 'Apr', 'May'))

custom_colors <- c("Jan" = "#285f93", "Feb" = "#5d9ad3", "Mar" = "#b9d4ec", "Apr" = '#b3e2b5', "May" = '#45B649')

ggplot() +
  geom_sf(data = usa_albers_county_map, aes(fill = Month), color= 'black', linewidth = 0.1) +
  #scale_fill_viridis(discrete=FALSE, direction = -1) +
  ggtitle('Syndromic Peak Month') +
  theme_void() +
  scale_fill_manual(values = custom_colors)


flu_surv_state_avg <- flu_surv_state |> group_by(state_fips, week_date) |>
  mutate(avg_value = mean(value, na.rm = T)) |>
  distinct(state_fips, week_date, avg_value) |>
  dplyr::filter(week_date > as.Date('2016-09-01') & week_date < as.Date('2017-08-01'))

  
county_peaks_surv <- flu_surv_state_avg |> 
  ungroup() |>
  dplyr::filter(!is.na(avg_value) & !is.na(state_fips)) |>
  # Calculate the row number
  mutate(index = row_number()) |>
  group_by(state_fips) |> 
  # Calculate multiple peaks greater than 10 cases and larger than 10 neighbors
  mutate(peaks_obs = seq_along(index) %in% 
           findpeaks(avg_value, MinPeakDistance = 10, MinPeakHeight = 1)$loc) |>
  # Remove first and last observations, which cannot be peaks but might be max
  dplyr::filter(row_number() != 1 & row_number() != n()) |>
  # Calculate the maximum value among the peaks
  mutate(is_max_obs = avg_value == max(avg_value)) |>
  select(c(state_fips, week_date, avg_value, peaks_obs, is_max_obs))

county_peaks_obs_surv <- county_peaks_surv |> 
  # Select the peaks
  dplyr::filter(peaks_obs == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_obs == T) |>
  select(c(state_fips, week_date, avg_value, peaks_obs, is_max_obs)) |>
  rename('week_surv' = 'week_date') |>
  group_by(state_fips) |>
  # Break any ties by taking the earliest peak
  arrange(week_surv) |>
  slice_head(n = 1) |>
  mutate(week_obs_num = week(week_surv)) 

usa_albers_st <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_state.rds'))   # convert to sf


usa_albers_county_map_surv <- left_join(usa_albers_st, 
                                        county_peaks_obs_surv, 
                                   by = c('STATEFP' = 'state_fips'))




usa_albers_county_map_surv$month <- month(usa_albers_county_map_surv$week_surv)
usa_albers_county_map_surv$Month <- factor(usa_albers_county_map_surv$month, levels = c(1, 2, 3, 4, 5), labels = c('Jan', 'Feb', 'Mar', 'Apr', 'May'))

custom_colors <- c("Jan" = "#285f93", "Feb" = "#5d9ad3", "Mar" = "#b9d4ec", "Apr" = '#b3e2b5', "May" = '#45B649')

ggplot() +
  geom_sf(data = usa_albers_county_map_surv, aes(fill = Month), color= 'black', linewidth = 0.1) +
  #scale_fill_viridis(discrete=FALSE, direction = -1) +
  ggtitle('FluSurv-NET Peak Month') +
  theme_void() +
  scale_fill_manual(values = custom_colors)

usa_albers_county_map <- left_join(usa_albers_county_map, 
                                   county_peaks_obs_surv , 
                                   by = c('STATEFP' = 'state_fips'))

usa_albers_county_map$diff <- usa_albers_county_map$week_pred - usa_albers_county_map$week_surv

ggplot() +
  geom_sf(data = usa_albers_county_map, aes(fill = as.numeric(diff)), color= 'black', linewidth = 0.1) +
  #scale_fill_viridis(discrete=FALSE, direction = -1) +
  ggtitle('Syndromic Peak Month') +
  theme_void() +
  scale_fill_gradient2(low = 'red', high = 'blue', midpoint = 0, mid = 'white')
  

library(lubridate)
date_example <- as.Date("2023-01-15")
iso_week <- isoweek(flu_county$year_week_dt)
iso_year <- isoyear(flu_county$year_week_dt)
flu_county$year_week <- paste0(iso_year, "-", sprintf("%02d", iso_week)) # Format week with leading zero

export <- flu_county |> ungroup() |>
  mutate(node = county_fips,
         year_week = year_week,
         incidence = pred_sum) |>
  dplyr::select(node, year_week, incidence)


write.csv(export, "./tmp/incidence_national_2016_2017.csv")



pred <- read.csv('./tmp/growth_rate_based_arrival_times_pred_10.csv', header = TRUE, sep = ',')
obs <- read.csv('./tmp/growth_rate_based_arrival_times_obs.csv', header = TRUE, sep = ',')

obs$onset_week_obs_adj = lubridate::parse_date_time(paste0(obs$onset_week, "-1"), "Y-W-u")

merge_onset <- left_join(pred, obs, by = 'node')
merge_onset <- merge_onset |>
  rename('onset_week_pred' = 'onset_week.x',
         'onset_week_obs' = 'onset_week.y') |>
  mutate(onset_week_pred_adj = lubridate::parse_date_time(paste0(onset_week_pred, "-1"), "Y-W-u"))

# Collapse data for graphing purposes
merge_onset_col <- merge_onset |>
  mutate(count = 1) |>
  group_by(onset_week_obs_adj, onset_week_pred_adj) |>
  mutate(count_sum = sum(count)) |>
  distinct(onset_week_obs_adj, onset_week_pred_adj, count_sum)

library(viridis)

merge_onset$diff <- as.numeric(as.Date(merge_onset$onset_week_pred_adj) - as.Date(merge_onset$onset_week_obs_adj))
rmse <- sqrt(sum((merge_onset$diff)^2, na.rm = T) / length(merge_onset$diff)) 
rmse
cor(as.numeric(merge_onset$onset_week_pred_adj), as.numeric(merge_onset$onset_week_obs_adj), 
    use = 'complete.obs', method = 'spearman')

forrest_theme <- theme(legend.position =  "bottom",
                       legend.position.inside = c(0.873, 0.237),
                       axis.text = element_text(size=14, color = 'black'),
                       axis.title = element_text(size=16),
                       strip.background = element_blank(),
                       legend.title = element_text(size=16),
                       legend.text = element_text(size=14),
                       plot.title = element_text(size=18, hjust = 0),
                       legend.box="horizontal",
                       legend.box.background = element_rect(colour = "white", 
                                                            fill = 'white'),
                       legend.key.width = unit(1.5, "cm")) 


onset_comp_plot <- ggplot(merge_onset_col[!is.na(merge_onset_col$onset_week_obs_adj),], aes(as.Date(onset_week_obs_adj), as.Date(onset_week_pred_adj), fill= count_sum)) + 
  geom_tile(color = 'white') + 
  scale_fill_viridis('Number of Counties\n', discrete=FALSE, option = 'C') +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black",
              linewidth = 1) +
  theme_minimal() + 
  labs(title = 'Season Onset Time Comparison',
       x = "Claims Onset Week",
       y = "Syndromic Onset Week") + forrest_theme +
  scale_x_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2016-09-01'), as.Date('2017-05-27'))) +
  scale_y_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2016-09-01'), as.Date('2017-05-27')))
onset_comp_plot




merge_onset$county_fips <- ifelse(nchar(merge_onset$node) == 4, str_pad(merge_onset$node, 5, pad = "0"), merge_onset$node)


merge_onset_peak <- left_join(merge_onset, county_peaks_both_filt, by = 'county_fips')

merge_onset_peak$diff_pred <- merge_onset_peak$week_pred - as.Date(merge_onset_peak$onset_week_pred_adj)
merge_onset_peak$diff_obs <- merge_onset_peak$week_obs - as.Date(merge_onset_peak$onset_week_obs_adj)


ggplot(merge_onset_peak) + geom_density(aes(diff_pred))

usa_albers_county_onset <- left_join(usa_albers_county, merge_onset_peak, by = c('GEOID' = 'county_fips'))


ggplot() +
  geom_sf(data = usa_albers_county_onset, aes(fill = as.numeric(diff_obs)), color= 'black', linewidth = 0.1) +
  #scale_fill_viridis(discrete=FALSE, direction = -1) +
  ggtitle('Syndromic Peak Month') +
  theme_void() 



growth <- read.csv('./tmp/weekly_growth_rate_2016_2017.csv', header = TRUE, sep = ',')

growth$week = lubridate::parse_date_time(paste0(growth$year_week, "-1"), "Y-W-u")


ggplot(growth[growth$node == '36061',]) + geom_line(aes(x = week, y = slope_3week))

ggplot(growth[growth$node == '36061',]) + geom_line(aes(x = week, y = incidence))




# state level
flu_state <- flu_pred_dat |> ungroup() |>
  group_by(year_week_dt, state_fips) |>
  mutate(pred_sum = sum(pred_count, na.rm = T),
         flu_sum = sum(patient_count * flu)) |>
  distinct(year_week_dt, state_fips, pred_sum, flu_sum)


# Set year and month variables
flu_state$month <- month(flu_state$year_week_dt)
flu_state$year <- year(flu_state$year_week_dt)

# Collapse to state level
all_cause_state <- all_cause |> ungroup() |>
  mutate(state_fips = substr(county_fips, 1, 2)) |>
  group_by(week_date, state_fips) |>
  mutate(all_cause_sum = sum(patient_count_imp)) |>
  distinct(week_date, state_fips, all_cause_sum) |>
  ungroup() |>
  group_by(state_fips) |>
  dplyr::filter(week_date < as.Date('2020-03-01')) |>
  dplyr::filter(week_date > as.Date('2016-08-31')) |>
  mutate(all_cause_avg = mean(all_cause_sum))

# Quickly visualize all cause data
ggplot(all_cause_state) + 
  geom_line(aes(x = week_date, y = all_cause_sum / all_cause_avg), color = 'blue') +
  facet_wrap(~state_fips)

# Scale observed and predicted data by all cause
flu_state <- left_join(flu_state, all_cause_state, by = c('year_week_dt' = 'week_date', 'state_fips'))

# Calculate adjusted values scales by all cause data
flu_state$flu_obs_scale <- flu_state$flu_sum / (flu_state$all_cause_sum)
flu_state$flu_pred_scale <- flu_state$pred_sum / (flu_state$all_cause_sum)

# Quick visualization of the scale variables
ggplot(flu_state) + 
  geom_line(aes(x = year_week_dt, y = flu_obs_scale), color = 'blue') +
  geom_line(aes(x = year_week_dt, y = flu_pred_scale), color = 'red') +
  facet_wrap(~state_fips, scale = 'free')

# Calculate summer averages
flu_state_summer <- flu_state |> 
  dplyr::filter(month < 9 & month > 5) |> group_by(year, state_fips) |>
  dplyr::filter(year > 2016) |>
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


flu_state <- flu_state |> arrange(year_week_dt) |> group_by(state_fips) |> 
  mutate(flu_sum_roll = rollmean(flu_sum, k = 4, align = 'right', na.pad = TRUE),
         flu_pred_roll = rollmean(pred_sum, k = 4, align = 'right', na.pad = TRUE))

# Calculate Z score
flu_state$z_score <- (flu_state$pred_sum - flu_state$pred_summ_mean) / sd(flu_state$pred_sum, na.rm = T)
#flu_state$obs_z <- (flu_state$flu_sum - flu_state$obs_summ_mean) / sd(flu_state$flu_sum, na.rm = T)



flu_state$Source <- 'Predicted'
flu_state$week_date <- flu_state$year_week_dt




state_peaks <- flu_state |> 
  dplyr::filter(!is.na(flu_sum_roll)) |>
  ungroup() |>
  # Calculate the row number
  mutate(index = row_number()) |>
  group_by(state_fips) |> 
  # Calculate multiple peaks greater than 10 cases and larger than 10 neighbors
  mutate(peaks_obs = seq_along(index) %in% 
           findpeaks(flu_sum_roll, MinPeakDistance = 10, MinPeakHeight = 20)$loc) |>
  mutate(peaks_pred = seq_along(index) %in% 
           findpeaks(flu_pred_roll, MinPeakDistance = 10, MinPeakHeight = 20)$loc) |>
  # Remove first and last observations, which cannot be peaks but might be max
  dplyr::filter(row_number() != 1 & row_number() != n()) |>
  # Calculate the maximum value among the peaks
  mutate(is_max_obs = flu_sum_roll == max(flu_sum_roll)) |>
  mutate(is_max_pred = flu_pred_roll == max(flu_pred_roll)) |>
  select(c(state_fips, week_date, flu_sum_roll,flu_pred_roll,  peaks_obs, peaks_pred, is_max_obs, is_max_pred))

state_peaks_obs <- state_peaks |> 
  # Select the peaks
  dplyr::filter(peaks_obs == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_obs == T) |>
  select(c(state_fips, week_date, flu_sum_roll, peaks_obs, is_max_obs)) |>
  rename('week_state' = 'week_date') |>
  group_by(state_fips) |>
  # Break any ties by taking the earliest peak
  arrange(week_state) |>
  slice_head(n = 1) |>
  mutate(week_obs_num = week(week_state)) 

state_peaks_pred <- state_peaks |> 
  # Select the peaks
  dplyr::filter(peaks_pred == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_pred == T) |>
  select(c(state_fips, week_date, flu_pred_roll, peaks_pred, is_max_pred)) |>
  rename('week_pred' = 'week_date') |>
  group_by(state_fips) |>
  # Break any ties by taking the earliest peak
  arrange(week_pred) |>
  slice_head(n = 1) |>
  mutate(week_obs_num = week(week_pred)) |>
  dplyr::select(state_fips, week_pred)


merge_onset_peak$state_fips <- substr(merge_onset_peak$county_fips, 1, 2)

early_peaks <- left_join(merge_onset_peak, state_peaks_obs, by = 'state_fips')

early_peaks <- early_peaks |> group_by(state_fips) |>
  mutate(value = as.numeric(week_state - week_pred)) |>
  dplyr::filter(!is.na(value)) |>
  dplyr::filter(value == max(value)) |>
  arrange(week_state) |>
  slice_head(n = 1) |>
  dplyr::select(c(state_fips, value)) |>
  mutate(Measure = 'Peak Time',
         order = 2)



ggplot(early_peaks) + geom_density(aes(min_peak_diff))


iso_week <- isoweek(flu_state$year_week_dt)
iso_year <- isoyear(flu_state$year_week_dt)
flu_state$year_week <- paste0(iso_year, "-", sprintf("%02d", iso_week)) # Format week with leading zero


export <- flu_state |> ungroup() |>
  mutate(node = state_fips,
         year_week = year_week,
         incidence = flu_sum) |>
  dplyr::select(node, year_week, incidence)


write.csv(export, "./tmp/incidence_national_2016_2017_state.csv")




state_onset <- read.csv('./tmp/growth_rate_based_arrival_times_obs_state.csv', header = TRUE, sep = ',')

state_onset$onset_week_state_adj = lubridate::parse_date_time(paste0(state_onset$onset_week, "-1"), "Y-W-u")

state_onset$state_fips <- ifelse(nchar(state_onset$node) == 1, str_pad(state_onset$node, 2, pad = "0"), state_onset$node)


merge_onset_state <- left_join(merge_onset_peak, state_onset, by = 'state_fips')

early_onsets <- merge_onset_state |> group_by(state_fips) |>
  mutate(value = as.numeric(as.Date(onset_week_state_adj) - as.Date(onset_week_pred_adj)))|>
  dplyr::filter(!is.na(value)) |>
  dplyr::filter(value == max(value)) |>
  arrange(onset_week_state_adj) |>
  slice_head(n = 1) |>
  dplyr::select(c(state_fips, value)) |>
  mutate(Measure = 'Onset Time',
         order = 1)


early_dat <- rbind(early_onsets, early_peaks)

forrest_theme <- theme(legend.position =  "none",
                       legend.position.inside = c(0.873, 0.237),
                       axis.text = element_text(size=14, color = 'black'),
                       axis.title = element_text(size=16),
                       strip.background = element_blank(),
                       legend.title = element_text(size=16),
                       legend.text = element_text(size=14),
                       plot.title = element_text(size=18, hjust = 0),
                       legend.box="horizontal",
                       legend.box.background = element_rect(colour = "white", 
                                                            fill = 'white'),
                       legend.key.width = unit(1.5, "cm")) 




measure_plot <- ggplot(early_dat) + 
  geom_boxplot(aes(y = reorder(Measure, -order), x = value, fill = Measure), width=0.35, linewidth = 0.8,
               outlier.size = 2.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black",
              linewidth = 1) +
  theme_minimal() + 
  labs(title = 'Syndromic and Spatial Lead Time',
       x = "Number of Days",
       y = "Measure") + forrest_theme +
  scale_fill_manual(values = c('Peak Time' = '#3cbb75ff' , 'Onset Time' = '#A92395ff'))

measure_plot


# Create figure, with labels
figure_3 <- cowplot::plot_grid(onset_comp_plot, peak_comp_plot, measure_plot,
                               nrow = 1,
                               labels = c("a","b", "c"),
                               label_size = 20)



# Save figure 
ggsave('./figs/figure_3.jpg', plot = figure_3, height = 6.5, width = 20)




flu_state <- left_join(flu_state, ili_state[, c(2, 3, 4)], by = c('year_week_dt' = 'week_date', 'state_fips'))
flu_state <- left_join(flu_state, nrevss_state[, c(2, 3, 4)], by = c('year_week_dt' = 'week_date', 'state_fips'))
flu_state <- left_join(flu_state, flu_surv_state_avg[, c(1, 2, 3)], by = c('year_week_dt' = 'week_date', 'state_fips'))


# CORRELATIONS # 

flu_corr <- flu_state |>
  group_by(state_fips) |>
  mutate(ili_corr = cor(pred_sum, value.x, use = 'pairwise.complete.obs'),
         nrevss_corr = cor(pred_sum, value.y, use = 'pairwise.complete.obs'),
         surv_corr = cor(pred_sum, avg_value, use = 'pairwise.complete.obs')) |>
  distinct(state_fips, ili_corr, nrevss_corr, surv_corr)



usa_albers_st <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_state.rds'))   # convert to sf


usa_albers_st_corr <- left_join(usa_albers_st , flu_corr, by = c('STATEFP' = 'state_fips'))


map_theme <- theme(legend.position = "bottom",
                   axis.title = element_text(size=16),
                   strip.background = element_blank(),
                   plot.title = element_text(size=18, hjust = 0.5),
                   panel.border = element_rect(fill=NA, linewidth = 0.8, color = 'white'),
                   legend.title = element_text(size = 16),
                   legend.text = element_text(size = 12),
                   legend.key.height = unit(0.38, 'cm'),
                   legend.key.width = unit(0.7, "cm"))



corr_ili <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_st_corr), aes(fill = ili_corr, group = STATEFP), color= 'black', linewidth = 0.1) +
  scale_fill_gradient('Correlation   \n', low = "white", high = '#347DC1', na.value = "grey",
                      breaks = c(0.25, 0.50, 0.75, 1), 
                      limits = c(0.25, 1.00)) + 
  ggtitle('ILI-NET\n') +
  theme_void() + map_theme 

corr_ili


corr_nrevss <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_st_corr), aes(fill = nrevss_corr, group = STATEFP), color= 'black', linewidth = 0.1) +
  scale_fill_gradient('Correlation   \n', low = "white", high = '#3cbb75ff', na.value = "grey",
                      breaks = c(0.25, 0.50, 0.75, 1), 
                      limits = c(0.25, 1.00)) + 
  ggtitle('NREVSS\n') +
  theme_void() + map_theme 

corr_nrevss

corr_surv <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_st_corr), aes(fill = surv_corr, group = STATEFP), color= 'black', linewidth = 0.1) +
  scale_fill_gradient('Correlation   \n', low = "white", high = '#9B59B6', na.value = "grey",
                      breaks = c(0.25, 0.50, 0.75, 1), 
                      limits = c(0.25, 1.00))+ 
  ggtitle('FluSurv-NET\n') +
  theme_void() + map_theme 

corr_surv

'#3cbb75ff'

'#9B59B6'

'#347DC1'






'#FF6B6B'




nys_county$county_fips <- as.character(nys_county$county_fips )

nys_county <- left_join(nys_county, flu_county[, c(1, 2, 3)], by = c('week_date' = 'year_week_dt',
                                                                     'county_fips'))

nys_corr <- nys_county |> group_by(county_fips) |>
  dplyr::filter(county_fips != '36041') |>
  mutate(nys_corr = cor(pred_sum, value, use = 'pairwise.complete.obs')) |>
  distinct(county_fips, nys_corr)


county <- read_sf(dsn = './raw/cb_2021_us_county_20m/', 
                  layer = 'cb_2021_us_county_20m')

county_nys <- county |> dplyr::filter(STATEFP == '36')

county_tomp <- county_nys |> dplyr::filter(COUNTYFP == '109' | COUNTYFP == '107')

# Merge on prev data
county_nys_corr <- left_join(county_nys, nys_corr, by = c('GEOID' = 'county_fips'))

ggplot() +
  geom_sf(data = st_as_sf(county_tomp), fill = 'white', color= 'black', linewidth = 1)
  
corr_map <- ggplot() +
  geom_sf(data = st_as_sf(county_nys_corr ), aes(fill = nys_corr, group = GEOID), color= 'black', linewidth = 0.1) +
  scale_fill_gradient('Correlation   \n', low = "white", high = "#e6a532", na.value = "grey",
                      breaks = c(0.25, 0.50, 0.75, 1), 
                      limits = c(0.25, 1.00)) + 
  ggtitle('NYS Lab-Confirmed') +
  theme_void() + map_theme 


corr_map

state_peaks_obs


figure_3_corr <- cowplot::plot_grid(corr_ili, corr_surv,corr_nrevss, corr_map,
                               nrow = 2,
                               labels = c("b"),
                               label_size = 20)

figure_3_corr


# PEAK TIMES #


flu_surv_state_avg <- flu_surv_state |> group_by(state_fips, week_date) |>
  mutate(avg_value = mean(value_roll, na.rm = T)) |>
  distinct(state_fips, week_date, avg_value) |>
  dplyr::filter(week_date > as.Date('2016-09-01') & week_date < as.Date('2017-08-01'))


state_peaks_surv <- flu_surv_state_avg |> 
  ungroup() |>
  dplyr::filter(!is.na(avg_value) & !is.na(state_fips)) |>
  # Calculate the row number
  mutate(index = row_number()) |>
  group_by(state_fips) |> 
  # Calculate multiple peaks greater than 10 cases and larger than 10 neighbors
  mutate(peaks_obs = seq_along(index) %in% 
           findpeaks(avg_value, MinPeakDistance = 10, MinPeakHeight = 1)$loc) |>
  # Remove first and last observations, which cannot be peaks but might be max
  dplyr::filter(row_number() != 1 & row_number() != n()) |>
  # Calculate the maximum value among the peaks
  mutate(is_max_obs = avg_value == max(avg_value)) |>
  select(c(state_fips, week_date, avg_value, peaks_obs, is_max_obs))

state_peaks_obs_surv <- county_peaks_surv |> 
  # Select the peaks
  dplyr::filter(peaks_obs == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_obs == T) |>
  select(c(state_fips, week_date, avg_value, peaks_obs, is_max_obs)) |>
  rename('week_surv' = 'week_date') |>
  group_by(state_fips) |>
  # Break any ties by taking the earliest peak
  arrange(week_surv) |>
  slice_head(n = 1) |>
  mutate(week_obs_num = week(week_surv)) |>
  dplyr::select(c(state_fips, week_surv))


state_peaks_ili <- ili_state |> 
  ungroup() |>
  dplyr::filter(!is.na(value_roll) & !is.na(state_fips)) |>
  mutate(value_roll = ifelse(value_roll < 0, 0, value_roll)) |>
  dplyr::filter(week_date > as.Date('2016-09-01') & week_date < as.Date('2017-08-01')) |>
  # Calculate the row number
  mutate(index = row_number()) |>
  group_by(state_fips) |> 
  # Calculate multiple peaks greater than 10 cases and larger than 10 neighbors
  mutate(peaks_obs = seq_along(index) %in% 
           findpeaks(value_roll, MinPeakDistance = 10, MinPeakHeight = 1)$loc) |>
  # Remove first and last observations, which cannot be peaks but might be max
  dplyr::filter(row_number() != 1 & row_number() != n()) |>
  # Calculate the maximum value among the peaks
  mutate(is_max_obs = value_roll == max(value_roll)) |>
  select(c(state_fips, week_date, value_roll, peaks_obs, is_max_obs))

state_peaks_obs_ili <- state_peaks_ili  |> 
  # Select the peaks
  dplyr::filter(peaks_obs == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_obs == T) |>
  select(c(state_fips, week_date, value_roll, peaks_obs, is_max_obs)) |>
  rename('week_ili' = 'week_date') |>
  group_by(state_fips) |>
  # Break any ties by taking the earliest peak
  arrange(week_ili) |>
  slice_head(n = 1) |>
  mutate(week_obs_num = week(week_ili)) |>
  dplyr::select(c(state_fips, week_ili))

state_peaks_nrevss <- nrevss_state |> 
  ungroup() |>
  dplyr::filter(!is.na(value_roll) & !is.na(state_fips)) |>
  mutate(value_roll = ifelse(value_roll < 0, 0, value_roll)) |>
  dplyr::filter(week_date > as.Date('2016-09-01') & week_date < as.Date('2017-08-01')) |>
  # Calculate the row number
  mutate(index = row_number()) |>
  group_by(state_fips) |> 
  # Calculate multiple peaks greater than 10 cases and larger than 10 neighbors
  mutate(peaks_obs = seq_along(index) %in% 
           findpeaks(value_roll, MinPeakDistance = 10, MinPeakHeight = 1)$loc) |>
  # Remove first and last observations, which cannot be peaks but might be max
  dplyr::filter(row_number() != 1 & row_number() != n()) |>
  # Calculate the maximum value among the peaks
  mutate(is_max_obs = value_roll == max(value_roll)) |>
  select(c(state_fips, week_date, value_roll, peaks_obs, is_max_obs))

state_peaks_obs_nrevss <- state_peaks_nrevss  |> 
  # Select the peaks
  dplyr::filter(peaks_obs == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_obs == T) |>
  select(c(state_fips, week_date, value_roll, peaks_obs, is_max_obs)) |>
  rename('week_nrevss' = 'week_date') |>
  group_by(state_fips) |>
  # Break any ties by taking the earliest peak
  arrange(week_nrevss) |>
  slice_head(n = 1) |>
  mutate(week_obs_num = week(week_nrevss)) |>
  dplyr::select(c(state_fips, week_nrevss))

county_peaks_nys <- nys_county |> 
  ungroup() |>
  dplyr::filter(!is.na(value_roll) & !is.na(county_fips)) |>
  mutate(value_roll = ifelse(value_roll < 0, 0, value_roll)) |>
  dplyr::filter(week_date > as.Date('2016-09-01') & week_date < as.Date('2017-08-01')) |>
  # Calculate the row number
  mutate(index = row_number()) |>
  group_by(county_fips) |> 
  # Calculate multiple peaks greater than 10 cases and larger than 10 neighbors
  mutate(peaks_obs = seq_along(index) %in% 
           findpeaks(value_roll, MinPeakDistance = 10, MinPeakHeight = 1)$loc) |>
  # Remove first and last observations, which cannot be peaks but might be max
  dplyr::filter(row_number() != 1 & row_number() != n()) |>
  # Calculate the maximum value among the peaks
  mutate(is_max_obs = value_roll == max(value_roll)) |>
  select(c(county_fips, week_date, value_roll, peaks_obs, is_max_obs))

county_peaks_obs_nys <- county_peaks_nys  |> 
  # Select the peaks
  dplyr::filter(peaks_obs == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_obs == T) |>
  select(c(county_fips, week_date, value_roll, peaks_obs, is_max_obs)) |>
  rename('week_nys' = 'week_date') |>
  group_by(county_fips) |>
  # Break any ties by taking the earliest peak
  arrange(week_nys) |>
  slice_head(n = 1) |>
  mutate(week_obs_num = week(week_nys)) |>
  dplyr::select(c(county_fips, week_nys))



state_peaks_corr <- left_join(state_peaks_pred, state_peaks_obs_ili, by = 'state_fips')

state_peaks_corr <- left_join(state_peaks_corr, state_peaks_obs_nrevss, by = 'state_fips')

state_peaks_corr <- left_join(state_peaks_corr, state_peaks_obs_surv, by = 'state_fips')



state_peaks_corr_ili <- state_peaks_corr |>
  mutate(count = 1) |>
  group_by(week_pred, week_ili) |>
  mutate(count_sum = sum(count, na.rm = T)) |>
  distinct(week_pred, week_ili, count_sum)



forrest_theme <- theme(legend.position =  "bottom",
                       legend.position.inside = c(0.873, 0.237),
                       axis.text = element_text(size=14, color = 'black'),
                       axis.title = element_text(size=16),
                       strip.background = element_blank(),
                       legend.title = element_text(size=16),
                       legend.text = element_text(size=12),
                       plot.title = element_text(size=18, hjust = 0.5),
                       legend.box="horizontal",
                       legend.box.background = element_rect(colour = "white", 
                                                            fill = 'white'),
                       legend.key.height = unit(0.38, 'cm'),
                       legend.key.width = unit(0.7, "cm")) 



corr_ili_scatter <- ggplot(state_peaks_corr_ili[!is.na(state_peaks_corr_ili$week_ili),], aes(as.Date(week_ili), as.Date(week_pred), fill= count_sum)) + 
  geom_tile(color = 'black') + 
  scale_fill_gradient2('Number of \nStates', low = 'white', high = '#347DC1') +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black",
              linewidth = 1) +
  theme_minimal() + 
  labs(title = 'ILI-NET',
       x = "ILI-NET Peak Week",
       y = "Syndromic Peak Week") + forrest_theme +
  scale_x_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2016-12-01'), as.Date('2017-05-27'))) +
  scale_y_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2016-12-01'), as.Date('2017-05-27')))
corr_ili_scatter


state_peaks_corr_surv <- state_peaks_corr |>
  mutate(count = 1) |>
  group_by(week_pred, week_surv) |>
  mutate(count_sum = sum(count, na.rm = T)) |>
  distinct(week_pred, week_surv, count_sum)


corr_surv_scatter <- ggplot(state_peaks_corr_surv[!is.na(state_peaks_corr_surv$week_surv),], aes(as.Date(week_surv), as.Date(week_pred), fill= count_sum)) + 
  geom_tile(color = 'black') + 
  scale_fill_gradient2('Number of \nStates', low = 'white', high = '#9B59B6') +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black",
              linewidth = 1) +
  theme_minimal() + 
  labs(title = 'FluSurv-NET',
       x = "FluSurv-NET Peak Week",
       y = "Syndromic Peak Week") + forrest_theme +
  scale_x_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2016-12-01'), as.Date('2017-05-27'))) +
  scale_y_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2016-12-01'), as.Date('2017-05-27')))
corr_surv_scatter

state_peaks_corr_nrevss <- state_peaks_corr |>
  mutate(count = 1) |>
  group_by(week_pred, week_nrevss) |>
  mutate(count_sum = sum(count, na.rm = T)) |>
  distinct(week_pred, week_nrevss, count_sum)


corr_nrevss_scatter <- ggplot(state_peaks_corr_nrevss[!is.na(state_peaks_corr_nrevss$week_nrevss),], aes(as.Date(week_nrevss), as.Date(week_pred), fill= count_sum)) + 
  geom_tile(color = 'black') + 
  scale_fill_gradient2('Number of \nStates', low = 'white', high = '#3cbb75ff') +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black",
              linewidth = 1) +
  theme_minimal() + 
  labs(title = 'NREVSS',
       x = "NREVSS Peak Week",
       y = "Syndromic Peak Week") + forrest_theme +
  scale_x_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2016-12-01'), as.Date('2017-05-27'))) +
  scale_y_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2016-12-01'), as.Date('2017-05-27')))
corr_nrevss_scatter

county_peaks_corr_nys <- left_join(county_peaks_obs_nys, county_peaks_pred, by = 'county_fips')


county_peaks_corr_nys_sum <- county_peaks_corr_nys |>
  mutate(count = 1) |>
  group_by(week_pred, week_nys) |>
  mutate(count_sum = sum(count, na.rm = T)) |>
  distinct(week_pred, week_nys, count_sum)


corr_nys_scatter <- ggplot(county_peaks_corr_nys_sum[!is.na(county_peaks_corr_nys_sum$week_pred),], aes(as.Date(week_nys), as.Date(week_pred), fill= count_sum)) + 
  geom_tile(color = 'black') + 
  scale_fill_gradient2('Number of \nCounties', low = 'white', high = '#e6a532') +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black",
              linewidth = 1) +
  theme_minimal() + 
  labs(title = 'NYS Lab-Confirmed',
       x = "NYS Peak Week",
       y = "Syndromic Peak Week") + forrest_theme +
  scale_x_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2016-12-01'), as.Date('2017-05-27'))) +
  scale_y_date(date_breaks = "2 months", date_labels = "%B", limits = c(as.Date('2016-12-01'), as.Date('2017-05-27')))
corr_nys_scatter





figure_3_scatter <- cowplot::plot_grid(corr_ili_scatter, corr_surv_scatter , corr_nrevss_scatter, corr_nys_scatter,
                                    nrow = 2,
                                    labels = c("c"),
                                    label_size = 20)

figure_3_scatter





# Collapse data for graphing purposes
merge_onset_col <- merge_onset |>
  mutate(count = 1) |>
  group_by(onset_week_obs_adj, onset_week_pred_adj) |>
  mutate(count_sum = sum(count)) |>
  distinct(onset_week_obs_adj, onset_week_pred_adj, count_sum)

duration_state <- rbind(flu_state[, c('state_fips', 'week_date', 'z_score', 'Source')],
                        ili_state[, c('state_fips', 'week_date', 'z_score', 'Source')],
                        nrevss_state[, c('state_fips', 'week_date', 'z_score', 'Source')],
                        flu_surv_state[, c('state_fips', 'week_date', 'z_score', 'Source')])

duration_dat <- duration_state |>
  dplyr::filter(week_date < as.Date('2017-09-01')) |>
  group_by(state_fips, Source) |>
  mutate(count = ifelse(z_score > 0.1, 1, 0),
         duration = sum(count, na.rm = T)) |>
  distinct(state_fips, Source, duration)


ggplot(duration_dat) + geom_violin(aes(x = Source, y = duration, fill = Source))


ggplot(flu_state) + 
  geom_line(aes(x = week_date, y = pred_z), color = 'red') +
  facet_wrap(~state_fips, scale = 'free')







flu_national <- flu_pred_dat |> group_by(year_week_dt) |>
  mutate(pred_sum = sum(pred_count, na.rm = T),
         flu_sum = sum(flu * patient_count)) |>
  distinct(year_week_dt, pred_sum, flu_sum) |>
  dplyr::filter(year_week_dt < as.Date('2020-03-01')) |>
  dplyr::filter(year_week_dt > as.Date('2016-08-31')) |>
  ungroup()

# Add month and year variables
flu_national$month <- month(flu_national$year_week_dt)
flu_national$year <- year(flu_national$year_week_dt)

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
  mutate(all_cause_sum = sum(patient_count_imp)) |>
  distinct(week_date, all_cause_sum) |>
  ungroup() |>
  dplyr::filter(week_date < as.Date('2020-03-01')) |>
  dplyr::filter(week_date > as.Date('2016-08-31')) |>
  mutate(all_cause_avg = mean(all_cause_sum))

# Quickly visualize all cause data
ggplot(all_cause_national) + 
  geom_line(aes(x = week_date, y = all_cause_sum / all_cause_avg), color = 'blue')

# Scale observed and predicted data by all cause
flu_national <- left_join(flu_national, all_cause_national, by = c('year_week_dt' = 'week_date'))

# Calculate adjusted values scales by all cause data
flu_national$flu_obs_scale <- flu_national$flu_sum *
  (1/(flu_national$all_cause_sum / flu_national$all_cause_avg))
flu_national$flu_pred_scale <- flu_national$pred_sum *
  (1/(flu_national$all_cause_sum / flu_national$all_cause_avg))

# Quick visualization of the scale variables
ggplot(flu_national) + 
  geom_line(aes(x = year_week_dt, y = flu_obs_scale), color = 'blue') +
  geom_line(aes(x = year_week_dt, y = flu_pred_scale), color = 'red')

# Calculate summer averages
flu_national_summer <- flu_national |> 
  dplyr::filter(month < 9 & month > 5) |> group_by(year) |>
  dplyr::filter(year > 2016) |>
  mutate(pred_summ_mean = mean(flu_pred_scale),
         obs_summ_mean = mean(flu_obs_scale)) |>
  distinct(year, pred_summ_mean, 
           obs_summ_mean)

# Merge on summer data to the correct season (summer 2016 should align with the 2016-2017 flu season)
flu_national_summer[nrow(flu_national_summer) + 1, ] <- flu_national_summer[1, ] 
flu_national_summer$year[2] <- 2016
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

flu_national <- left_join(flu_national, ili_nat[, c(3, 6)], by = c('year_week_dt' = 'week_date'))
flu_national <- flu_national |>
  rename('ili_z' = 'z_score')

flu_national <- left_join(flu_national, nrevss_nat[, c(3, 6)], by = c('year_week_dt' = 'week_date'))
flu_national <- flu_national |>
  rename('nrevss_z' = 'z_score')

flu_national <- left_join(flu_national, flu_surv_nat[, c(2, 6)], by = c('year_week_dt' = 'week_date'))
flu_national <- flu_national |>
  rename('surv_z' = 'z_score')


flu_national$ili_above <- ifelse(flu_national$ili_z > 0.1, 1, 0)
flu_national$nrevss_above <- ifelse(flu_national$nrevss_z > 0.1, 1, 0)
flu_national$surv_above <- ifelse(flu_national$surv_z > 0.1, 1, 0)

flu_national_duration <- flu_national |>
  group_by(Season) |>
  mutate(pred_duration = sum(pred_above),
         ili_duration = sum(ili_above),
         nrevss_duration = sum(nrevss_above),
         surv_duration = sum(surv_above, na.rm = T)) |>
  distinct(Season, pred_duration, ili_duration, 
           nrevss_duration, surv_duration)

flu_national_duration_long <- flu_national_duration |>
  pivot_longer(
    cols = pred_duration:surv_duration, 
    names_to = "Source",
    values_to = "value"
  ) |>
  dplyr::filter(Season != '2019-2020')


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
  dplyr::filter(Source != 'Claims')

Lines <- c('Claims' = 'solid', 'FluSurv-NET' = 'solid', 
           'ILI-NET' = 'solid', 'NREVSS' = 'solid', 'Predicted' = 'dashed')

Colors <- c('Claims' = '#bf6fa7', 'FluSurv-NET' = '#9B59B6', 
            'ILI-NET' = '#347DC1', 'NREVSS' = '#3cbb75ff', 'Predicted' = 'gray39')

percent_theme <- theme(legend.position = "none",
                       axis.text = element_text(size=14),
                       axis.text.x  = element_text(size=14, angle = 45, hjust = 1),
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
  labs(title = "Season Duration by Source") + theme_minimal() + 
  percent_theme

flu_national_duration_plot

figure_4 <- cowplot::plot_grid(flu_national_duration_plot
                               , figure_3_corr, figure_3_scatter,
                               rel_widths = c(0.5, 1, 1),
                               nrow = 1,
                               labels = c("a", "", ""),
                               label_size = 20)

ggsave('./figs/figure_4.jpg', plot = figure_4, height = 8, width = 20)


