
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

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

# Load data
flu_pred_dat <- readRDS('./tmp/flu_pred_dat_county_p2_state.rds')

# Collapse to county level
flu_county <- flu_pred_dat |> 
  group_by(week_date, county_fips) |>
  mutate(flu_sum = sum(patient_sum * flu)) |>
  distinct(week_date, flu_sum, urban_code_binary, state_fips, south) |>
  filter(week_date < as.Date('2020-03-01')) |>
  filter(week_date > as.Date('2016-08-31'))

# Make separate data for peaks analysis
flu_county_peak <- flu_county

# Add month and year
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
  mutate(obs_peak = peaks(flu_sum, span=101, strict=FALSE, endbehavior = 1)) |>
  select(c(Season, county_fips, week_date, flu_sum,
           obs_peak, urban_code_binary, state_fips, south))

# Identify peaks
county_peaks_obs <- county_peaks_usa |> filter(obs_peak == T) |>
  select(c(Season, county_fips, week_date, obs_peak, south, urban_code_binary)) |>
  rename('county_fips' = 'county_fips',
         'week_obs' = 'week_date') |>
  group_by(Season, county_fips) |>
  arrange(week_obs) |>
  slice_head(n = 1) |>
  mutate(week_obs_num = week(week_obs))

# Merge on population and all cause data
all_cause_season <- read_csv('raw/county_season_ac_v2.csv')
county_group <- readRDS('tmp/county_group_pop.rds')

# Merge
county_peaks_obs <- left_join(county_peaks_obs, all_cause_season, by = c('Season' = 'season', 'county_fips'))
county_peaks_obs <- left_join(county_peaks_obs, county_group, by = c('county_fips'))

# Filter
county_peaks_filt <- county_peaks_obs |>
  filter(county_group_pop > 10000) |>
  filter(all_cause > 10000)

# Map
usa_albers_state <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_state.rds'))   # convert to sf
usa_albers_county <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_county.rds'))   # convert to sf
x_walk <- readRDS('./tmp/county_group_xwalk.rds')

# Collapse county map to county group
usa_albers_county_x_walk <- left_join(x_walk, usa_albers_county[, c(5, 6)], by = c('fips' = 'GEOID'))
usa_albers_county_group <- usa_albers_county_x_walk %>% 
  group_by(county_fips) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()

# Merge on peak data
usa_albers_county_group_plot <- left_join(usa_albers_county_group, 
                                     county_peaks_filt[county_peaks_filt$Season == '2016-2017',],
                                     by = c('county_fips' = 'county_fips'))


usa_albers_county_group_plot$week_obs_num_new  <- ifelse(usa_albers_county_group_plot$week_obs_num > 24, usa_albers_county_group_plot$week_obs_num - 52, usa_albers_county_group_plot$week_obs_num)

# Set theme
map_theme <- theme(legend.position = "right",
                   axis.title = element_text(size=16),
                   strip.background = element_blank(),
                   panel.border = element_rect(fill=NA, linewidth = 0.6, color = '#FFFFFF00'),
                   legend.title = element_text(size = 16),
                   legend.text = element_text(size = 12),
                   legend.key.height = unit(0.7, 'cm'),
                   legend.key.width = unit(0.4, "cm"))

# Display Map
ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_group_plot), aes(fill = week_obs_num_new, group = county_fips), color= 'black', linewidth = 0.05) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.10) +
  scale_fill_gradient2('Peak', low = "blue", mid = 'white', high = "red", na.value = "grey", midpoint = 0) + 
  theme_void() + map_theme + ggtitle('2016-17') +
  scale_x_date(date_labels = "%b %Y")

urban <- readRDS("./tmp/urban_south_binary.rds")

usa_albers_county_group_plot <- left_join(usa_albers_county_group_plot, 
                                          urban[, c(1, 3)],
                                          by = c('county_fips' = 'county_fips'))

# Calculate centroid for each county group
usa_albers_county_group_plot$centroid <- st_centroid(usa_albers_county_group_plot$geometry)

# Convert Albers projection to traditional Lat/Long
usa_albers_county_group_plot$latlon <- st_transform(usa_albers_county_group_plot$centroid, crs = 4326) # EPSG:4326 for WGS84



usa_albers_county_group_plot <- 
  usa_albers_county_group_plot |>
  group_by(county_fips) |>
  mutate(lat = latlon[[1]][2],
         lon = latlon[[1]][1])

usa_albers_county_group_plot$latlon[[1]]



model <- lm(data = usa_albers_county_group_plot, week_obs_num_new ~ lat + lon  + all_cause + county_group_pop)

summary(model)


ggplot() +
  geom_density(data = usa_albers_county_group_plot, aes(x = week_obs_num_new))

