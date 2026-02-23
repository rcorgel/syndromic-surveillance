####################
# 1. SET-UP SCRIPT #
####################

# Start with a clear environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(cowplot)
library(lubridate)
library(plotROC)
library(pROC)
library(car)
library(cutpointr)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')



xwalk <- readRDS('./tmp/county_group_xwalk.rds') 
county_pop <- read.csv('./raw/co-est2019-alldata.csv', stringsAsFactors=FALSE)

county_pop$fips_state <- ifelse(nchar(county_pop$STATE) == 1, paste0("0", as.character(county_pop$STATE)), as.character(county_pop$STATE))
county_pop$fips_county <- ifelse(nchar(county_pop$COUNTY) == 1, paste0("00", as.character(county_pop$COUNTY)), as.character(county_pop$COUNTY))
county_pop$fips_county <- ifelse(nchar(county_pop$COUNTY) == 2, paste0("0", as.character(county_pop$COUNTY)), county_pop$fips_county)

county_pop$fips <- paste(county_pop$fips_state, county_pop$fips_county, sep = '')

county_pop_group <- left_join(xwalk, county_pop[, c('fips', 'POPESTIMATE2019')], by = c('fips' = 'fips'))

county_pop_groups <- county_pop_group |> group_by(county_fips) |>
  mutate(county_group_pop = sum(POPESTIMATE2019, na.rm = TRUE)) |>
  distinct(county_fips, county_group_pop)


saveRDS(county_pop_groups, './tmp/county_group_pop.rds') 
