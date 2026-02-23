################################################################################
# File Name: analyze_virus_data                                                #
#                                                                              #
# Purpose:   Analyze other viruses from medical claims.                        #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load data                                                      #
#            3. Analyze data                                                   #
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
library(cowplot)
library(lubridate)
library(ISOweek)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1CgyyYhIWFVG2aMITI-HBpKamvwwpQ_ED/CHC_data/other_resp/2024-06-24/')

################
# 2. LOAD DATA #
################

# Load data
other_virus <- read.csv('other_resp_pathogens.csv')

# Convert week to date
other_virus$week_date <- sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", other_virus$year_week)
other_virus$week_date <- ISOweek2date(other_virus$week_date)

# Make unique lists
week_list <- unique(other_virus$week_date)
county_list <- unique(other_virus$county_fips)

# Make a shell of all county/week combinations
shell <- NULL
for (i in county_list) {
  shell <- rbind(shell, data.frame(rep(i, 418), week_list))
}
colnames(shell) <- c('county_fips', 'week_date')

# Merge data to county/week shell to examine weeks with 0's
resp_virus <- left_join(shell, other_virus, by = c('week_date' = 'week_date',
                                                   'county_fips' = 'county_fips'))

# Impute the low-count data, NOTE: WHAT ABOUT EXPLICIT 0's?
resp_virus  <- mutate(resp_virus, adenovirus = ifelse(is.na(adenovirus), 0, adenovirus))
resp_virus$adenovirus[resp_virus$adenovirus == '<=5'] <- 
  sample(0:5, length(resp_virus$adenovirus[resp_virus$adenovirus == '<=5']), replace= TRUE)

resp_virus  <- mutate(resp_virus, enterovirus = ifelse(is.na(enterovirus), 0, enterovirus))
resp_virus$enterovirus[resp_virus$enterovirus == '<=5'] <- 
  sample(0:5, length(resp_virus$enterovirus[resp_virus$enterovirus == '<=5']), replace= TRUE)

resp_virus  <- mutate(resp_virus, rhinovirus = ifelse(is.na(rhinovirus), 0, rhinovirus))
resp_virus$rhinovirus[resp_virus$rhinovirus == '<=5'] <- 
  sample(0:5, length(resp_virus$rhinovirus[resp_virus$rhinovirus == '<=5']), replace= TRUE)

resp_virus  <- mutate(resp_virus, parainfluenza = ifelse(is.na(parainfluenza), 0, parainfluenza))
resp_virus$parainfluenza[resp_virus$parainfluenza == '<=5'] <- 
  sample(0:5, length(resp_virus$parainfluenza[resp_virus$parainfluenza == '<=5']), replace= TRUE)

resp_virus  <- mutate(resp_virus, metapneumovirus = ifelse(is.na(metapneumovirus), 0, metapneumovirus))
resp_virus$metapneumovirus[resp_virus$metapneumovirus == '<=5'] <- 
  sample(0:5, length(resp_virus$metapneumovirus[resp_virus$metapneumovirus == '<=5']), replace= TRUE)

# Remove year_week variable
resp_virus <- resp_virus |> select(-c(year_week)) 

# Exploratory analysis
resp_virus_county_week <- resp_virus |> group_by(week_date) |>
  mutate(adenovirus = sum(as.numeric(adenovirus)),
         enterovirus = sum(as.numeric(enterovirus)),
         rhinovirus = sum(as.numeric(rhinovirus)),
         parainfluenza = sum(as.numeric(parainfluenza)),
         metapneumovirus = sum(as.numeric(metapneumovirus))) |>
  distinct(week_date, adenovirus, enterovirus, rhinovirus,
           parainfluenza, metapneumovirus)

ggplot(resp_virus_county_week) +
  geom_line(aes(x = week_date, y = adenovirus), color = '#F8766D') +
  geom_line(aes(x = week_date, y = enterovirus), color = '#A3A500') +
  geom_line(aes(x = week_date, y = rhinovirus), color = '#00BF7D') +
  geom_line(aes(x = week_date, y = parainfluenza), color = '#00B0F6') +
  geom_line(aes(x = week_date, y = metapneumovirus), color = '#E76BF3') +
  theme_minimal() + ylab('count')

nyc_virus <- resp_virus |> filter(county_fips == '17031')

nyc_virus_long <- nyc_virus |> select(-c(county_fips)) |>
  pivot_longer(!week_date, names_to = "virus", values_to = "count")

ggplot(nyc_virus_long) +
  geom_line(aes(x = week_date, y = as.numeric(count), color = virus)) +
  theme_minimal() + ylab('count') + theme(legend.position = 'bottom') + 
  ggtitle('Chicago, 2016-2023')




