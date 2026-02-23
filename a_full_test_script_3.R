# Start with a clear environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(lubridate)
library(car)
library(plotROC)
library(pROC)
library(cutpointr)
library(zoo)
library(sf)

# Set seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

# Load data
flu_county_week <- readRDS('./tmp/flu_county_week.rds') 

# Filter data
flu_county_week_filt <- flu_county_week |> 
  filter(county_fips == "17031" | county_fips == "53033" | 
         county_fips == "06037" | county_fips == "04013" |
         county_fips == "36061" | county_fips == "48113" |
         county_fips == "11001" | county_fips == "13121" |
         county_fips == "29510") |>
  pivot_longer(!c(week_date, county_fips), names_to = "flu_type", values_to = "count")

# Test plot
ggplot(flu_county_week_filt) + geom_line(aes(y = count, x = week_date, 
                                             group = flu_type, color = flu_type)) +
  facet_wrap(vars(county_fips), scales = "free") 

# Calculate a monthly rolling average
flu_county_week_long <- flu_county_week |> 
  # group_by(county_fips) |>
  # mutate(pred_flu_cases_roll = rollmean(pred_flu_cases, k = 3, 
  #                                       align = 'right', fill = NA),
  #        flu_cases = rollmean(flu_cases, k = 3,
  #                             align = 'right', fill = NA)) |>
  # select(-c(pred_flu_cases)) |>
  pivot_longer(!c(week_date, county_fips), names_to = "flu_type", values_to = "count")

# Add time variables
flu_county_week_long$month <- month(flu_county_week_long$week_date)
flu_county_week_long$year <- year(flu_county_week_long$week_date)
# Alter year for season merging
flu_county_week_long$year <- ifelse(flu_county_week_long$month < 9,
                                    flu_county_week_long$year - 1,
                                    flu_county_week_long$year)

# Calculate a summer mean
flu_county_week_long_summer <- flu_county_week_long |>
  filter(month < 9 & month > 5) |> filter(year > 2016) |>
  group_by(county_fips, year, flu_type) |>
  mutate(summer_mean = mean(count, na.rm = TRUE)) |>
  distinct(county_fips, year, flu_type, summer_mean)
flu_county_week_long <- left_join(flu_county_week_long, flu_county_week_long_summer,
                                  by = c('year' = 'year',
                                         'county_fips' = 'county_fips',
                                         'flu_type' = 'flu_type'))

flu_county_week_long <- flu_county_week_long |>
  filter(year > 2016) |> group_by(flu_type) |>
  mutate(sd = sd(count, na.rm = TRUE)) |> ungroup() |>
  mutate(count_scale = (count - summer_mean) / sd)


# Filter data
flu_county_week_long_filt <- flu_county_week_long |> 
  filter(county_fips == "17031" | county_fips == "53033" | 
           county_fips == "06037" | county_fips == "04013" |
           county_fips == "36061" | county_fips == "48113" |
           county_fips == "11001" | county_fips == "13121" |
           county_fips == "29510")

# Test plot
ggplot(flu_county_week_long_filt) + geom_line(aes(y = count_scale, x = week_date, 
                                             group = flu_type, color = flu_type)) +
  facet_wrap(vars(county_fips), scales = "free") 


all_week <- read.csv('./raw/county_week_ac_v2.csv')


library(ISOweek)
all_week$week_date <- sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", all_week$year_week)
all_week$week_date <- ISOweek2date(all_week$week_date)

long <- left_join(flu_county_week_long_filt, all_week, by = c('county_fips' = 'county_fips',
                                         'week_date' = 'week_date'))

long$percent <- long$count / as.numeric(long$all_cause)

ggplot(long) + geom_line(aes(y = percent, x = week_date, group = flu_type, color = flu_type)) +
  facet_wrap(vars(county_fips), scales = "free") 





# Box plot of one season (obs vs pred)
county_pop_groups <- readRDS('./tmp/county_group_pop.rds') 

flu_county_week_long <- left_join(flu_county_week_long, county_pop_groups, 
                                  by = c('county_fips' = 'county_fips'))
flu_county_week_long_lim <- flu_county_week_long |> filter(county_group_pop > 20000)

flu_county_week_long_2017_max <- flu_county_week_long_lim |>
  filter(year == 2019) |> group_by(county_fips, flu_type) |> 
  arrange(county_fips, week_date) |>
  slice_max(count, n=1) |>
  slice(1)

ggplot(flu_county_week_long_2017_max, aes(x=flu_type, y=week_date, fill=flu_type)) + 
  geom_boxplot()

# Map of one season (obs vs pred)
# Load maps
usa_albers_state <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_state.rds'))   # convert to sf
usa_albers_county <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_county.rds')) # convert to sf
county_xwalk <- readRDS('./tmp/county_group_xwalk.rds') 

# Collapse county map to county group
usa_albers_county_x_walk <- left_join(county_xwalk, usa_albers_county, by = c('fips' = 'GEOID'))
usa_albers_county_group <- usa_albers_county_x_walk %>% 
  group_by(county_fips) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()

# Merge data
max_obs_2017 <- flu_county_week_long_2017_max |> 
  filter(flu_type == 'flu_cases') |> select(c(week_date, county_fips)) |>
  rename('max_week_date_obs' = 'week_date')

max_pred_2017 <- flu_county_week_long_2017_max |> 
  filter(flu_type == 'pred_flu_cases') |> select(c(week_date, county_fips)) |>
  rename('max_week_date_pred' = 'week_date')

# Merge on data
usa_albers_county_group <- left_join(usa_albers_county_group, max_obs_2017,
                                     by = c('county_fips' = 'county_fips'))
usa_albers_county_group <- left_join(usa_albers_county_group, max_pred_2017,
                                     by = c('county_fips' = 'county_fips'))

# Create maps
library(scales)
map <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_group), aes(fill = as.Date(max_week_date_obs), 
                                                        group = county_fips), color= 'black', linewidth = 0.15) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.4) +
  ggtitle('Influenza Season Max, 2017-18 Season') + 
  theme_void() + theme(legend.position = 'right',
                       plot.title = element_text(size = 20, hjust = 0.5),
                       panel.border = element_rect(fill=NA, linewidth = 0.8, color = 'white'),
                       legend.title = element_text(size = 12),
                       legend.text = element_text(size = 12),
                       legend.key.size = unit(0.6, 'cm')) +
  scale_fill_date(labels = date_format("%Y-%B-%d"), low = 'white', high = 'blue')
map

map <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_group), aes(fill = as.Date(max_week_date_pred), 
                                                        group = county_fips), color= 'black', linewidth = 0.15) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.4) +
  ggtitle('Influenza Season Max, 2017-18 Season') + 
  theme_void() + theme(legend.position = 'right',
                       plot.title = element_text(size = 20, hjust = 0.5),
                       panel.border = element_rect(fill=NA, linewidth = 0.8, color = 'white'),
                       legend.title = element_text(size = 12),
                       legend.text = element_text(size = 12),
                       legend.key.size = unit(0.6, 'cm')) +
  scale_fill_date(labels = date_format("%Y-%B-%d"), low = 'white', high = 'blue')
map

usa_albers_county_group$diff <- usa_albers_county_group$max_week_date_obs - usa_albers_county_group$max_week_date_pred

map <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_group), aes(fill =  as.numeric(diff), 
                                                        group = county_fips), color= 'black', linewidth = 0.15) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.4) +
  ggtitle('Influenza Season Max, 2019-20 Season Diff') + 
  theme_void() + theme(legend.position = 'right',
                       plot.title = element_text(size = 20, hjust = 0.5),
                       panel.border = element_rect(fill=NA, linewidth = 0.8, color = 'white'),
                       legend.title = element_text(size = 12),
                       legend.text = element_text(size = 12),
                       legend.key.size = unit(0.6, 'cm')) +
  scale_fill_gradient2(low = 'red', mid = 'white', high = 'blue', midpoint = 0)
map


# Comparison of all 3 seasons (difference)

# Correlations
flu_county_week_wide_lim_2017 <- flu_county_week_long_lim |>
  filter(year == 2017) |> select(c(county_fips, week_date, flu_type, count)) |>
  pivot_wider(names_from = flu_type, values_from = count) |>
  group_by(county_fips) |>
  mutate(correlation = cor(flu_cases, pred_flu_cases)) |>
  distinct(county_fips, correlation)

usa_albers_county_group <- left_join(usa_albers_county_group, flu_county_week_wide_lim_2017,
                                     by = c('county_fips' = 'county_fips'))
map <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_group), aes(fill =  as.numeric(correlation), 
                                                        group = county_fips), color= 'black', linewidth = 0.15) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.4) +
  ggtitle('Influenza Season Max, 2019-20 Season Corr') + 
  theme_void() + theme(legend.position = 'right',
                       plot.title = element_text(size = 20, hjust = 0.5),
                       panel.border = element_rect(fill=NA, linewidth = 0.8, color = 'white'),
                       legend.title = element_text(size = 12),
                       legend.text = element_text(size = 12),
                       legend.key.size = unit(0.6, 'cm')) +
  scale_fill_gradient2(low = 'white', high = 'blue', midpoint = 0)
map


# Raw numbers divided by all cause (magnitude)



# Clinical relevance???
# Try one more model







