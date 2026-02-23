library(zoo)
library(tidyverse)
library(usdata)

# Start with a clear environment
rm(list = ls())

# Set working directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

# Load data
eyes <- read.csv('raw/county_week_conjunctivitis.csv')

# Impute rows where patient count is <=5 to a random number 1 to 5
eyes$patient_count_imp <- eyes$conjunctivitis
eyes$patient_count_imp[eyes$conjunctivitis == '<=5'] <- 
  sample(1:5, length(eyes$patient_count_imp[eyes$conjunctivitis == '<=5']), 
         replace= TRUE)
# Convert imputed count to numeric
eyes$patient_count_imp <- as.numeric(eyes$patient_count_imp)

library(ISOweek)
# Convert week to date
eyes$week_date <- sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", eyes$year_week)
eyes$week_date <- ISOweek2date(eyes$week_date)

eyes_filt <- eyes |> dplyr::filter(county_fips == '06067')
eyes_nyc <- eyes |> dplyr::filter(county_fips == '36061')
eyes_nyc <- eyes_nyc |> dplyr::mutate(roll_mean = rollmean(patient_count_imp, k = 8, align = 'right', fill = NA))
eyes_filt <- eyes_filt |> dplyr::mutate(roll_mean = rollmean(patient_count_imp, k = 8, align = 'right', fill = NA))
ggplot(eyes_nyc) + geom_line(aes(y = roll_mean, x = week_date))

# other <- read.csv('raw/other_resp_pathogens.csv')
# other <- other |> dplyr::select(c(year_week, county_fips, adenovirus))
# other$adenovirus_imp <- other$adenovirus
# other$adenovirus_imp[other$adenovirus == '<=5'] <- 
#   sample(1:5, length(other$adenovirus_imp[other$adenovirus == '<=5']), 
#          replace= TRUE)
# 
# # Convert imputed count to numeric
# other$adenovirus_imp <- as.numeric(other$adenovirus_imp)
# 
# # Convert week to date
# other$week_date <- sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", other$year_week)
# other$week_date <- ISOweek2date(other$week_date)
# 
# other_nyc <- other |> filter(county_fips == '06037')
# other_nyc <- other_nyc |> mutate(roll_mean = rollmean(adenovirus_imp, k = 7, align = 'right', fill = NA))
# ggplot(other_nyc) + geom_line(aes(y = roll_mean, x = week_date))

weeks <- eyes_nyc$week_date
county <- unique(eyes$county_fips)
full_eyes <- expand.grid(weeks, county)
full_eyes <- full_eyes |> 
  rename('county_fips' = 'Var2',
         'week_date' = 'Var1')
full_eyes <- left_join(full_eyes, eyes[, c(1, 4, 5)], by = c('county_fips' = 'county_fips',
                                                'week_date' = 'week_date'))
full_eyes$patient_count_imp <- ifelse(is.na(full_eyes$patient_count_imp), 0, 
                                      full_eyes$patient_count_imp)

full_eyes$year <- year(full_eyes$week_date)
eyes_2016 <- full_eyes |> filter(year == 2016) |> 
  group_by(county_fips) |> 
  mutate(mean_2016 = mean(patient_count_imp),
         sd_2016 = sd(patient_count_imp)) |>
  distinct(county_fips, mean_2016, sd_2016)

full_eyes <- left_join(full_eyes, eyes_2016, by = c('county_fips' = 'county_fips'))
full_eyes <- full_eyes |> filter(sd_2016 != 0)


full_eyes_state <- full_eyes |> 
  mutate(state_fips = as.numeric(substr(county_fips, 1, 2))) |>
  group_by(state_fips, week_date) |>
  mutate(conjunc = sum(patient_count_imp)) |>
  distinct(state_fips, week_date, conjunc) |>
  filter(state_fips < 57)

full_eyes_state$year <- year(full_eyes_state$week_date)

eyes_state_2016 <- full_eyes_state |> filter(year == 2016) |> 
  group_by(state_fips) |> 
  mutate(mean_2016 = mean(conjunc),
         sd_2016 = sd(conjunc)) |>
  distinct(state_fips, mean_2016, sd_2016)

full_eyes_state <- left_join(full_eyes_state, eyes_state_2016, by = c('state_fips' = 'state_fips'))

eyes_roll_state <- full_eyes_state |> 
  group_by(state_fips) |>
  mutate(roll_mean = rollmean(conjunc, k = 4, align = 'right', fill = NA)) |>
  mutate(roll_mean_scale = (roll_mean - mean_2016) / sd_2016,
         percent_max = roll_mean / max(roll_mean, na.rm = TRUE))

ggplot(eyes_roll_state) + geom_line(aes(y = roll_mean_scale, x = week_date, 
                                  color = state_fips, group = state_fips),
                              alpha = 0.2) +
  theme(legend.position = 'none')


conjunc_late_2023_county <- full_eyes |>
  mutate(month = month(week_date)) |> 
  filter(year == 2023 & month > 11)

conjunc_late_2019_county <- full_eyes |>
  mutate(month = month(week_date)) |> 
  filter(year == 2019 & month > 11)


conjunc_late_2023 <- eyes_roll_state |>
  mutate(month = month(week_date)) |> 
  filter(year == 2023 & month > 8)

conjunc_late_2019 <- eyes_roll_state |>
  mutate(month = month(week_date)) |> 
  filter(year == 2019 & month > 8)

all_week <- read.csv('./raw/county_week_ac_v2.csv')

library(ISOweek)
all_week$week_date <- sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", all_week$year_week)
all_week$week_date <- ISOweek2date(all_week$week_date)
all_week$state_fips = as.numeric(substr(all_week$county_fips, 1, 2))

all_week$patient_count_imp <- all_week$all_cause
all_week$patient_count_imp[all_week$all_cause == '<=5'] <- 
  sample(1:5, length(all_week$patient_count_imp[all_week$all_cause == '<=5']), 
         replace= TRUE)
# Convert imputed count to numeric
all_week$patient_count_imp <- as.numeric(all_week$patient_count_imp)
all_week$all_cause <- all_week$patient_count_imp
all_week$patient_count_imp <- NULL
all_week_state <- all_week |> group_by(state_fips, week_date) |>
  mutate(all_cause = sum(patient_count_imp)) |>
  distinct(state_fips, week_date, all_cause)

long <- left_join(conjunc_late_2023 , all_week_state, by = c('state_fips' = 'state_fips',
                                         'week_date' = 'week_date'))
long_2 <- left_join(conjunc_late_2019 , all_week_state, by = c('state_fips' = 'state_fips',
                                                             'week_date' = 'week_date'))

conjunc_late_2023 <- long |> group_by(state_fips) |>
  mutate(conjunc_per = sum(conjunc) / sum(all_cause)) |>
  distinct(state_fips, conjunc_per)

conjunc_late_2019 <- long_2 |> group_by(state_fips) |>
  mutate(conjunc_per_2019 = sum(conjunc) / sum(all_cause)) |>
  distinct(state_fips, conjunc_per_2019) 

long_county <- left_join(conjunc_late_2023_county , all_week, by = c('county_fips' = 'county_fips',
                                                             'week_date' = 'week_date'))
long_2_county <- left_join(conjunc_late_2019_county , all_week, by = c('county_fips' = 'county_fips',
                                                               'week_date' = 'week_date'))

conjunc_late_2023 <- long |> group_by(state_fips) |>
  mutate(conjunc_per = sum(conjunc) / sum(all_cause)) |>
  distinct(state_fips, conjunc_per)

conjunc_late_2019 <- long_2 |> group_by(state_fips) |>
  mutate(conjunc_per_2019 = sum(conjunc) / sum(all_cause)) |>
  distinct(state_fips, conjunc_per_2019) 

conjunc_late_2023 <- left_join(conjunc_late_2023, conjunc_late_2019,
                               by = c('state_fips' = 'state_fips'))
conjunc_late_2023$ratio <- conjunc_late_2023$conjunc_per / conjunc_late_2023$conjunc_per_2019

conjunc_late_2023_county <- long_county |> group_by(county_fips) |>
  mutate(conjunc_per = sum(patient_count_imp) / sum(all_cause)) |>
  distinct(county_fips, conjunc_per)

conjunc_late_2019_county <- long_2_county |> group_by(county_fips) |>
  mutate(conjunc_per_2019 = sum(patient_count_imp) / sum(all_cause)) |>
  distinct(county_fips, conjunc_per_2019) 

conjunc_late_2023_county <- left_join(conjunc_late_2023_county, conjunc_late_2019_county,
                               by = c('county_fips' = 'county_fips'))
conjunc_late_2023_county$ratio <- conjunc_late_2023_county$conjunc_per / conjunc_late_2023_county$conjunc_per_2019
conjunc_late_2023_county$diff <- conjunc_late_2023_county$conjunc_per - conjunc_late_2023_county$conjunc_per_2019
#fips_info(conjunc_late_2023_county$county_fips, sortAndRemoveDuplicates = FALSE)


# State Level Cows
library(readxl)
cows <- read_excel('./raw/milkcowsandprod.xlsx', sheet = 'Milk cows', range = "B2:BD74")
cows <- cows |> filter(!is.na(`Back to content page.`)) |>
  rename('state' = 'Back to content page.',
         'cows_2023' = '2023')
library(usmap)
cows$state_fips <- as.numeric(fips(state = cows$state))
cows <- cows |> select(c(state, state_fips, cows_2023)) |>
  filter(!is.na(cows_2023)) |> mutate(cows_2023 = as.numeric(cows_2023) * 1000)

# County Groups
county_groups <- readRDS('./tmp/county_group_xwalk.rds')
county_names <- read_csv("./raw/ssa_fips_state_county_2023.csv")
county_names <- county_names |>
  select(c(countyname_fips, state, fipscounty)) |>
  mutate(state_name = abbr2state(state),
       county_name = gsub('COUNTY','', countyname_fips),
       county_name = gsub('PARISH','', county_name),
       county_name = gsub('CENSUS AREA','', county_name),
       county_name = gsub('BOROUGH','', county_name),
       county_name = trimws(county_name),
       county_name = gsub( " ", "", tolower(county_name)),
       state_name = gsub( " ", "", tolower(state_name)))

county_groups <- left_join(county_groups, county_names, 
                           by = c('fips' = 'fipscounty'))
# County Level Cows
cows_county <- read.csv('./raw/Dairy cattle - Number of heads (4).csv')
cows_county <- cows_county |> filter(!is.na(OBJECTID)) |>
  select(c(State, County, Number.of.heads.in.2017)) |>
  mutate(num_cows = ifelse(Number.of.heads.in.2017 == -999, 0, Number.of.heads.in.2017)) |>
  mutate(county_name = gsub( " ", "", tolower(County)),
         state_name = gsub( " ", "", tolower(State)))

county_groups <- left_join(county_groups, cows_county, 
                           by = c('county_name' = 'county_name',
                                  'state_name' = 'state_name'))
collapse_cows <- county_groups |> group_by(county_fips) |>
  mutate(cows = sum(num_cows, na.rm = TRUE)/1000) |>
  distinct(county_fips, cows)

county_pop <- read_csv("./raw/uscounties.csv")

conjunc_late_2023_county <- left_join(conjunc_late_2023_county, collapse_cows, 
                                      by = c('county_fips' = 'county_fips'))
conjunc_late_2023_county <- left_join(conjunc_late_2023_county, county_pop[, c(4, 9)], 
                                      by = c('county_fips' = 'county_fips'))

conjunc_late_2023_county$cows_per <- (conjunc_late_2023_county$cows*1000) / conjunc_late_2023_county$population

conjunc_late_2023_county_filt <- conjunc_late_2023_county |>
  filter(conjunc_per != 0 & cows_per != 0)

summary(lm(data = conjunc_late_2023_county_filt, diff ~ cows_per))

ggplot(conjunc_late_2023_county_filt) +
  geom_point(aes(x = log(cows_per), y = diff)) +
  geom_smooth(method='lm', aes(x = log(cows_per), y = diff), 
              formula= y~x)

ggplot(conjunc_late_2023_county_filt) +
  geom_point(aes(x = log(cows_per), y = log(conjunc_per))) +
  geom_smooth(method='lm', aes(x = log(cows_per), y = log(conjunc_per)), 
              formula= y~x)

ggplot(conjunc_late_2023_county_filt) +
  geom_point(aes(x = log(cows_per), y = diff)) +
  geom_smooth(method='lm', aes(x = log(cows_per), y = diff), 
              formula= y~x)

library(sf)
# Load maps
usa_albers_state <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_state.rds'))   # convert to sf
usa_albers_county <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_county.rds')) # convert to sf

# Collapse county map to county group
usa_albers_county_x_walk <- left_join(county_groups[, c(1,2)], usa_albers_county[, c(5, 6)], by = c('fips' = 'GEOID'))
usa_albers_county_group <- usa_albers_county_x_walk %>% 
  group_by(county_fips) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()

usa_albers_county_group <- left_join(usa_albers_county_group, conjunc_late_2023_county,
                                     by = c('county_fips' = 'county_fips'))
# Create maps
map <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_group), aes(fill = conjunc_per, group = county_fips), color= 'black', linewidth = 0.10) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.3) +
  scale_fill_gradient('Conjunctivitis', low = "white", high = "black", limits = c(0, 1)) + ggtitle('Conjunctivitis Percent, December 2023') + 
  theme_void() + theme(legend.position = 'bottom',
                       plot.title = element_text(size = 16, hjust = 0.5),
                       panel.border = element_rect(fill=NA, linewidth = 0.8, color = 'white'),
                       legend.title = element_text(size = 12),
                       legend.text = element_text(size = 12),
                       legend.key.height = unit(0.4, 'cm'),
                       legend.key.width = unit(1.5, "cm")) 
map

library(stringr)
pop <- read_excel('./raw/NST-EST2024-POP.xlsx', range = "A4:G60")
pop <- pop |> rename('state' = '...1',
                     'pop_2023' = '2023') |>
  select(c(state, pop_2023)) |>
  mutate(state = substr(state, 2, 30)) |>
  mutate(state_fips = as.numeric(fips(state = state))) |>
  filter(!is.na(state_fips))
  

conjunc_late_2023 <- left_join(conjunc_late_2023, cows, 
                               by = c('state_fips' = 'state_fips'))

conjunc_late_2023 <- left_join(conjunc_late_2023, pop, 
                               by = c('state_fips' = 'state_fips'))

conjunc_late_2023_combo <- conjunc_late_2023 |> filter(!is.na(cows_2023))
conjunc_late_2023_combo$cows_per <- conjunc_late_2023_combo$cows_2023 / conjunc_late_2023_combo$pop_2023

summary(lm(data = conjunc_late_2023_combo, ratio ~ cows_per))



eyes_roll <- full_eyes |> 
  group_by(county_fips) |>
  mutate(roll_mean = rollmean(patient_count_imp, k = 7, align = 'right', fill = NA)) |>
  mutate(roll_mean_scale = (roll_mean - mean_2016) / sd_2016,
         percent_max = roll_mean / max(roll_mean, na.rm = TRUE))











ggplot(eyes_roll) + geom_line(aes(y = percent_max, x = week_date, 
                                  color = county_fips, group = county_fips),
                              alpha = 0.2) +
  theme(legend.position = 'none')

eyes_roll_filt <- eyes_roll |> filter(percent_max == 1 & year == 2023)

test <- eyes_roll |> filter(percent_max == 1 & year == 2023) |>
  distinct(county_fips)
test$keep <-  1

# examine trends in high 2023 counties
# Examine where high 2023 counties are located
# See if counties are more likely to be high in 2023 if they had more cows, control for population
# Need to bring in all cause numbers - conjunctivitus in 2023 / all cause

