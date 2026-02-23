################################################################################
# File Name: create_county_group_xwalk                                         #
#                                                                              #
# Purpose:   Create a county to county group crosswalk. Due to privacy         #
#            concerns, small counties are grouped together in the medical      #
#            claims data. To align data, a crosswalk must be created.          #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Create crosswalk                                               #
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
library(assertr)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

#######################
# 2. CREATE CROSSWALK #
#######################

# Load processed claims data
flu <- readRDS('./tmp/flu_p1_dat_proc.rds')
covid <- readRDS('./tmp/covid_p1_dat_proc.rds')
rsv <- readRDS('./tmp/rsv_p1_dat_proc.rds')

# Combine unique counties and county groups
flu_county <- flu |> select(county_fips)
covid_county <- covid |> select(county_fips)
rsv_county <- rsv |> select(county_fips)
remove(flu, covid, rsv)
# Combine all counties and county groups
counties <- rbind(flu_county, covid_county, rsv_county)
# Remove repeated counties and county groups
counties_uniq <- unique(counties)

# Create county to county group x-walk
x_walk <- counties_uniq |> separate(county_fips, sep = "_", remove = FALSE, 
                                    into = c('fips_1', 'fips_2', 'fips_3', 'fips_4',
                                             'fips_5', 'fips_6', 'fips_7', 'fips_8')) |>
  pivot_longer(cols=c('fips_1', 'fips_2', 'fips_3', 'fips_4', 'fips_5', 'fips_6', 
                      'fips_7', 'fips_8'),
               names_to='fips_num',
               values_to='fips') |>
  filter(!is.na(fips))
x_walk$state_fips <- substr(x_walk$county_fips, 1, 2)

# Remove PR and unknown (99999)
# These should already be removed from data processing
x_walk <- x_walk |> filter(state_fips != '72') |>
  filter(county_fips != '99999')

# Conduct some sanity checks on the data
# Counties should be unique
verify(x_walk, length(fips) == length(unique(fips)))
# Number of county groups should be the same as in the data
county_group_num <- nrow(counties_uniq)
verify(x_walk, length(unique(county_fips)) == county_group_num)

# Remove number of counties in a group variable
x_walk <- x_walk |> select(-c(fips_num))

# Save xwalk
saveRDS(x_walk, './tmp/county_group_xwalk.rds') 

################################################################################
################################################################################

