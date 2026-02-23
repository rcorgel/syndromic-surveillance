

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

################
# 2. LOAD DATA #
################


rsv_v1 <- read.csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1CgyyYhIWFVG2aMITI-HBpKamvwwpQ_ED/CHC_data/rsv_cases/2024-02-08/county_week_rsv_v1.csv')
rsv_v2 <- read.csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1CgyyYhIWFVG2aMITI-HBpKamvwwpQ_ED/CHC_data/rsv_cases/2024-02-08/county_week_rsv_v2.csv')

rsv_v2$rsv_imp <- rsv_v2$rsv
rsv_v2$rsv_imp[rsv_v2$rsv == '<=5'] <- 
  sample(1:5, length(rsv_v2$rsv_imp[rsv_v2$rsv_imp == '<=5']), replace= TRUE)

rsv_check <- rsv_v2 |>
  mutate(under_5 = ifelse(rsv == '<=5', 'under 5', 'over')) |>
  group_by(under_5) |>
  mutate(sum = sum(as.numeric(rsv_imp))) |>
  ungroup() |>
  distinct(under_5, sum) |>
  mutate(total = sum(sum),
         prop = sum / total)


rsv_time <- rsv_v2 |> group_by(year_week) |>
  mutate(cases = sum(as.numeric(rsv_imp)),
         year_week_date = as.Date(paste(year_week, 1, sep = '-'), format = '%Y-%U-%u')) |>
  distinct(year_week_date, cases)

rsv_time_na <- na.omit(rsv_time)

ggplot(rsv_time_na) + 
  geom_line(aes(x = year_week_date, y = cases), color = '#00BA42') +
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('RSV') +
  theme(plot.title = element_text(hjust = 0.5))


rsv_age_v1 <- read.csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1CgyyYhIWFVG2aMITI-HBpKamvwwpQ_ED/CHC_data/rsv_cases/2024-02-08/county_week_age_rsv_v1.csv')
rsv_age_v2 <- read.csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1CgyyYhIWFVG2aMITI-HBpKamvwwpQ_ED/CHC_data/rsv_cases/2024-02-08/county_week_age_rsv_v2.csv')

rsv_age_v2$rsv_imp <- rsv_age_v2$rsv
rsv_age_v2$rsv_imp[rsv_age_v2$rsv == '<=5'] <- 
  sample(1:5, length(rsv_age_v2$rsv_imp[rsv_age_v2$rsv_imp == '<=5']), replace= TRUE)

rsv_check_age <- rsv_age_v2 |>
  mutate(under_5 = ifelse(rsv == '<=5', 'under 5', 'over')) |>
  group_by(under_5) |>
  mutate(sum = sum(as.numeric(rsv_imp))) |>
  ungroup() |>
  distinct(under_5, sum) |>
  mutate(total = sum(sum),
         prop = sum / total)


rsv_time_age <- rsv_age_v2 |> group_by(year_week, age_grp) |>
  mutate(cases = sum(as.numeric(rsv_imp)),
         year_week_date = as.Date(paste(year_week, 1, sep = '-'), format = '%Y-%U-%u')) |>
  distinct(year_week_date, age_grp, cases)

rsv_time_age_na <- na.omit(rsv_time_age)

ggplot(rsv_time_age_na) + 
  geom_line(aes(x = year_week_date, y = cases, color = age_grp)) +
  theme_minimal() + xlab('Date') + ylab('Count') + ggtitle('RSV') +
  theme(plot.title = element_text(hjust = 0.5))

