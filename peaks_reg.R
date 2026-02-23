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
library(cutpointr)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

# Load full data
flu_2016 <- readRDS('./tmp/flu_filt_symp_p2_2017.rds')

# 
# prev_2016 <- flu_2016 |> 
#   group_by(county_fips) |>
#   mutate(flu_sum = sum(patient_count_imp * flu)) |>
#   dplyr::filter(state_fips != '02' & state_fips != '15') |>
#   distinct(county_fips, flu_sum)
# 
# prev_2016 <- left_join(prev_2016,
#                        all_cause_season[all_cause_season$season == '2017-2018',], 
#                        by = c('county_fips'))
# prev_2016$prev <- prev_2016$flu_sum / prev_2016$all_cause
# 
# saveRDS(prev_2016, './tmp/flu_prev_2017.rds')

flu_2016$week <- epiweek(flu_2016$year_week_dt)

# Collapse to county level
# COUNTY LEVEL
flu_county <- flu_2016 |> 
  group_by(year_week_dt, county_fips) |>
  mutate(flu_sum = sum(patient_count * flu)) |>
  dplyr::filter(state_fips != '02' & state_fips != '15') |>
  distinct(year_week_dt, county_fips, flu_sum)
  
# Add year and month variables
flu_county$year <- year(flu_county$year_week_dt)
flu_county$month <- month(flu_county$year_week_dt)

# Calculate peaks
county_peaks_usa <- flu_county |> 
  ungroup() |>
  # Calculate the row number
  mutate(index = row_number()) |>
  group_by(county_fips) |> 
  # Calculate multiple peaks greater than 10 cases and larger than 10 neighbors
  mutate(peaks_obs = seq_along(index) %in% 
           findpeaks(flu_sum, MinPeakDistance = 10, MinPeakHeight = 10)$loc) |>
  # Remove first and last observations, which cannot be peaks but might be max
  dplyr::filter(row_number() != 1 & row_number() != n()) |>
  # Calculate the maximum value among the peaks
  mutate(is_max_obs = flu_sum == max(flu_sum)) |>
  select(c(county_fips, year_week_dt, flu_sum, peaks_obs, is_max_obs))

map_theme <- theme(legend.position = "right",
                   axis.title = element_text(size=16),
                   axis.text = element_text(size=16),
                   title = element_text(size=16),
                   strip.background = element_blank(),
                   panel.border = element_rect(fill=NA, linewidth = 0.8, color = 'white'),
                   legend.title = element_text(size = 16),
                   legend.text = element_text(size = 12),
                   legend.key.height = unit(0.7, 'cm'),
                   legend.key.width = unit(0.5, "cm"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())

library(scales)
ggplot(county_peaks_usa[county_peaks_usa$county_fips == '36061',]) + 
  geom_line(aes(x = as.Date(year_week_dt) - 60, y = flu_sum, group = county_fips), color = 'black', linewidth = 2.5) +
  theme_minimal() + map_theme +
  scale_x_date(
    date_breaks = "3 months", # Place a tick mark every 6 months
    date_labels = "%b %Y"     # Format as "Month Abbreviation Year" (e.g., "Jan 2023")
  ) + ylab('Influenza Cases') + xlab('Time') 

county_peaks_usa_all <- county_peaks_usa |> group_by(year_week_dt) |>
  mutate(sum = sum(flu_sum)) |>
  distinct(year_week_dt, sum)

ggplot(county_peaks_usa_all) + 
  geom_line(aes(x = as.Date(year_week_dt), y = sum), color = 'black', linewidth = 2.5) +
  theme_classic() + map_theme +
  scale_x_date(
    date_breaks = "3 months", # Place a tick mark every 6 months
    date_labels = "%b %Y"     # Format as "Month Abbreviation Year" (e.g., "Jan 2023")
  ) + ylab('Influenza Cases') + xlab('Time')

county_peaks_obs <- county_peaks_usa |> 
  # Select the peaks
  dplyr::filter(peaks_obs == T) |>
  # Limit the peaks to max peak
  dplyr::filter(is_max_obs == T) |>
  select(c(county_fips, year_week_dt, peaks_obs, is_max_obs)) |>
  rename('week_obs' = 'year_week_dt') |>
  group_by(county_fips) |>
  # Break any ties by taking the earliest peak
  arrange(week_obs) |>
  slice_head(n = 1) |>
  mutate(week_obs_num = week(week_obs)) #|>
  # Filter for significant peaks w/ z score > 1
  #dplyr::filter(obs_z > 1)

# Convert peak weeks
county_peaks_obs <- county_peaks_obs |>
  mutate(week_obs_num = ifelse(week_obs_num > 30, week_obs_num - 52, week_obs_num))

# Load maps
usa_albers_state <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_state.rds'))   # convert to sf
usa_albers_county <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_county.rds')) # convert to sf

# Load county to county group x walk
x_walk <- readRDS('./tmp/county_group_xwalk.rds')

# Collapse county map to county group
usa_albers_county_x_walk <- left_join(x_walk, usa_albers_county[, c(1, 5, 6)], by = c('fips' = 'GEOID'))
usa_albers_county_group <- usa_albers_county_x_walk %>% 
  group_by(county_fips, state_fips) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()

# Merge on counties
usa_albers_county_group <- left_join(usa_albers_county_group, county_peaks_obs, by = 'county_fips')

# Map
ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_group), aes(fill = week_obs_num, group = county_fips), color= 'black', linewidth = 0.01) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.10) +
  scale_fill_gradient2('Peak Week\n(Jan. 1 = 0)', low = '#d7642c', mid = 'white', high ="#347DC1", midpoint = 1,
                       na.value = "grey85") + 
  theme_void()

# Add on geograpgy variables
usa_albers_county_group$centroid <- st_centroid(usa_albers_county_group$geometry)
usa_albers_county_group$latlon <- st_transform(usa_albers_county_group$centroid, crs = 4326) # EPSG:4326 for WGS84
usa_albers_county_group <- usa_albers_county_group |>
  group_by(county_fips) |>
  mutate(lat = latlon[[1]][2],
         long = latlon[[1]][1])

# Add on population variables
urban <- readRDS("./tmp/urban_south_binary.rds")
usa_albers_county_group <- left_join(usa_albers_county_group, urban, by = 'county_fips')
all_cause_season <- read_csv('raw/county_season_ac_v2.csv')
usa_albers_county_group <- left_join(usa_albers_county_group, 
                                     all_cause_season[all_cause_season$season == '2017-2018',], 
                                     by = c('county_fips'))
county_group_pop <- readRDS('tmp/county_group_pop.rds')
usa_albers_county_group <- left_join(usa_albers_county_group, county_group_pop, by = c('county_fips'))

# Previous seasons prev
prev <- readRDS('tmp/flu_prev_2016.rds')
usa_albers_county_group <- left_join(usa_albers_county_group, prev, by = c('county_fips'))


# Add school variables
schools <- read_csv('raw/school_calendar_by_county_23_24.csv')
schools$county_fips <- ifelse(nchar(schools$fips) == 4, paste0('0', schools$fips), schools$fips)
usa_albers_county_group <- left_join(usa_albers_county_group, 
                                     schools[schools$Event == 'First Day of School', c(2, 3, 6)], 
                                     by = c('county_fips'))
usa_albers_county_group <- left_join(usa_albers_county_group, 
                                     schools[schools$Event == 'Christmas Break', c(2, 3, 6)], 
                                     by = c('county_fips'))


usa_albers_county_group$week_obs_num <- ifelse(usa_albers_county_group$county_group_pop < 10000, NA, usa_albers_county_group$week_obs_num)
usa_albers_county_group$week_obs_num <- ifelse(usa_albers_county_group$all_cause.x < 10000, NA, usa_albers_county_group$week_obs_num)

usa_albers_county_group$pop <- usa_albers_county_group$county_group_pop / 100000
model <- glm(week_obs_num ~ prev*urban_code_binary + urban_code_binary + lat*urban_code_binary,
                  data = usa_albers_county_group)
summary(model)

ggplot(usa_albers_county_group, aes(x = week_obs_num, y = prev)) + 
  geom_point() +
  geom_smooth(method = "lm")





