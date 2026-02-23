################################################################################
# File Name: urban_south_xwalk                                                 #
#                                                                              #
# Purpose:   Create a crosswalk for county to urban and south binaries.        #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Create south binary                                            #
#            3. Create urban binary                                            #
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

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

############################
# 2. CREATE SOUTH VARIABLE #
############################

# Load county map
usa_albers_county <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_county.rds')) # convert to sf

# Calculate centroid of each county
usa_albers_county$centroid <- st_centroid(usa_albers_county$geometry)

# Convert to lat/long from Albers projection
usa_albers_county$latlon <- st_transform(usa_albers_county$centroid, crs = 4326)

# Create south binary variable
usa_albers_county <- usa_albers_county |>
  group_by(GEOID) |>
  mutate(county_fips = GEOID,
         south = ifelse(latlon[[1]][2] <= 30, 1, 0),
         state_fips = substr(county_fips, 1, 2),
         south = ifelse(state_fips == '02', 0, south))

############################
# 3. CREATE URBAN VARIABLE #
############################

# Load Urban/Rural classification
urban <- read.csv('raw/urban_rural_cat.csv', header = T)

# Clean up data
urban <- urban |> 
  mutate(county_fips = ifelse(nchar(as.character(FIPS.code)) == 4,
                              paste0('0', as.character(FIPS.code)),
                              as.character(FIPS.code)),
         urban_code = as.numeric(substr(X2023.Code, 1, 1)))

# Merge urban data to south data
urban_merge <- left_join(usa_albers_county, urban, by = c('county_fips' = 'county_fips'))

# Fill in missing CT data (since CT redrew county lines in 2023)
urban_merge[urban_merge$county_fips == '09001',]$urban_code <- 3
urban_merge[urban_merge$county_fips == '09003',]$urban_code <- 1
urban_merge[urban_merge$county_fips == '09005',]$urban_code <- 5
urban_merge[urban_merge$county_fips == '09007',]$urban_code <- 2
urban_merge[urban_merge$county_fips == '09009',]$urban_code <- 3
urban_merge[urban_merge$county_fips == '09011',]$urban_code <- 3
urban_merge[urban_merge$county_fips == '09013',]$urban_code <- 2
urban_merge[urban_merge$county_fips == '09015',]$urban_code <- 3

# Create urban binary variable
urban_merge <- urban_merge |> group_by(county_fips) |>
  mutate(urban_code_binary = ifelse(urban_code < 3, 1, 0)) |>
  distinct(county_fips, south, urban_code, urban_code_binary)

# Save data
saveRDS(urban_merge, "./tmp/urban_south_binary.rds", )

################################################################################
################################################################################
