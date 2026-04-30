################################################################################
# File Name: 03a_sample_flu_data                                               #
#                                                                              #
# Purpose:   Take random samples of visits from influenza medical claims.      #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Sample by season                                               #
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
library(arrow)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

#######################
# 2. SAMPLE BY SEASON #
#######################

########################
# VERSION 1: 2016-2019 #
########################

# Even sample of cases and non-cases 
# Sample equally across seasons, geography
# End with 500,000 cases and 500,000 non-cases
# 2016
flu_1 <- read_parquet('tmp/flu_1_proc_p2.parquet')
flu_balanced_16 <- flu_1 |> 
  group_by(flu) |> 
  uncount(patient_count) |> 
  sample_n(125000)
remove(flu_1)
# 2017
flu_2 <- read_parquet('tmp/flu_2_proc_p2.parquet')
flu_balanced_17 <- flu_2 |> 
  group_by(flu) |> 
  uncount(patient_count) |> 
  sample_n(125000)
remove(flu_2)
# 2018
flu_3 <- read_parquet('tmp/flu_3_proc_p2.parquet')
flu_balanced_18 <- flu_3 |> 
  group_by(flu) |> 
  uncount(patient_count) |> 
  sample_n(125000)
remove(flu_3)
# 2019
flu_4 <- read_parquet('tmp/flu_4_proc_p2.parquet')
flu_balanced_19 <- flu_4 |> 
  group_by(flu) |> 
  uncount(patient_count) |> 
  sample_n(125000)
remove(flu_4)

# Combine balanced data
flu_v1 <- rbind(flu_balanced_16, flu_balanced_17, flu_balanced_18, flu_balanced_19)
saveRDS(flu_v1, 'tmp/flu_sample_v1.rds')

########################
# VERSION 2: 2020-2023 #
########################

# Even sample of cases and non-cases 
# Sample equally across seasons, geography
# End with 500,000 cases and 500,000 non-cases

# 2020
flu_5 <- read_parquet('tmp/flu_5_proc_p2.parquet')
flu_balanced_20 <- flu_5 |> 
  group_by(flu) |> 
  uncount(patient_count) |> 
  sample_n(166666)
remove(flu_5)
# 2021
flu_6 <- read_parquet('tmp/flu_6_proc_p2.parquet')
flu_balanced_21 <- flu_6 |> 
  group_by(flu) |> 
  uncount(patient_count) |> 
  sample_n(166667)
remove(flu_6)
# 2022
flu_7 <- read_parquet('tmp/flu_7_proc_p2.parquet')
flu_balanced_22 <- flu_7 |> 
  group_by(flu) |> 
  uncount(patient_count) |> 
  sample_n(166667)
remove(flu_7)

# Combine balanced data
flu_v2 <- rbind(flu_balanced_20, flu_balanced_21, flu_balanced_22)
saveRDS(flu_v2, 'tmp/flu_sample_v2.rds')

###########################
# VERSION 3: DEDUPLICATED #
###########################

# Set the seed
set.seed(123)


# 2016
flu_1 <- read_parquet('tmp/flu_1_proc_p1.parquet')

flu_balanced_16 <- flu_1 |> 
  group_by(flu) |> 
  uncount(patient_count) |> 
  sample_n(125000)
remove(flu_1)
# 2017
flu_2 <- read_parquet('tmp/flu_2_proc_p1.parquet')
flu_balanced_17 <- flu_2 |> 
  group_by(flu) |> 
  uncount(patient_count) |> 
  sample_n(125000)
remove(flu_2)
# 2018
flu_3 <- read_parquet('tmp/flu_3_proc_p1.parquet')
flu_balanced_18 <- flu_3 |> 
  group_by(flu) |> 
  uncount(patient_count) |> 
  sample_n(125000)
remove(flu_3)
# 2019
flu_4 <- read_parquet('tmp/flu_4_proc_p1.parquet')
flu_balanced_19 <- flu_4 |> 
  group_by(flu) |> 
  uncount(patient_count) |> 
  sample_n(125000)
remove(flu_4)

# Combine balanced data
flu_v3 <- rbind(flu_balanced_16, flu_balanced_17, flu_balanced_18, flu_balanced_19)
saveRDS(flu_v3, 'tmp/flu_sample_v3.rds')

########################
# VERSION 3: GEOGRAPHY #
########################


# 2016
flu_1 <- read_parquet('tmp/flu_1_proc_p2.parquet')
flu_1 <- flu_1 |>
  dplyr::filter(state_fips == '41' | state_fips == '12' |
                  state_fips == '25' | state_fips == '17')
flu_balanced_16 <- flu_1 |> 
  group_by(flu, state_fips) |> 
  uncount(patient_count) |> 
  sample_n(10000)
remove(flu_1)
# 2017
flu_2 <- read_parquet('tmp/flu_2_proc_p2.parquet')
flu_2 <- flu_2 |>
  dplyr::filter(state_fips == '41' | state_fips == '12' |
                  state_fips == '25' | state_fips == '17')
flu_balanced_17 <- flu_2 |> 
  group_by(flu, state_fips) |> 
  uncount(patient_count) |> 
  sample_n(10000)
remove(flu_2)
# 2018
flu_3 <- read_parquet('tmp/flu_3_proc_p2.parquet')
flu_3 <- flu_3 |>
  dplyr::filter(state_fips == '41' | state_fips == '12' |
                  state_fips == '25' | state_fips == '17')
flu_balanced_18 <- flu_3 |> 
  group_by(flu, state_fips) |> 
  uncount(patient_count) |> 
  sample_n(10000)
remove(flu_3)
# 2019
flu_4 <- read_parquet('tmp/flu_4_proc_p2.parquet')
flu_4 <- flu_4 |>
  dplyr::filter(state_fips == '41' | state_fips == '12' |
                  state_fips == '25' | state_fips == '17')
flu_balanced_19 <- flu_4 |> 
  group_by(flu, state_fips) |> 
  uncount(patient_count) |> 
  sample_n(10000)
remove(flu_4)

# Combine balanced data
flu_v4 <- rbind(flu_balanced_16, flu_balanced_17, flu_balanced_18, flu_balanced_19)
saveRDS(flu_v4, 'tmp/flu_sample_v4.rds')

#####################
# VERSION 3: SEASON #
#####################

# 2016
flu_1 <- read_parquet('tmp/flu_1_proc_p2.parquet')
flu_1 <- flu_1 |>
  dplyr::filter(month == 1 | month == 4 |
                  month == 7 | month == 10)
flu_balanced_16 <- flu_1 |> 
  group_by(flu, month) |> 
  uncount(patient_count) |> 
  sample_n(10000)
remove(flu_1)
# 2017
flu_2 <- read_parquet('tmp/flu_2_proc_p2.parquet')
flu_2 <- flu_2 |>
  dplyr::filter(month == 1 | month == 4 |
                  month == 7 | month == 10)
flu_balanced_17 <- flu_2 |> 
  group_by(flu, month) |> 
  uncount(patient_count) |> 
  sample_n(10000)
remove(flu_2)
# 2018
flu_3 <- read_parquet('tmp/flu_3_proc_p2.parquet')
flu_3 <- flu_3 |>
  dplyr::filter(month == 1 | month == 4 |
                  month == 7 | month == 10)
flu_balanced_18 <- flu_3 |> 
  group_by(flu, month) |> 
  uncount(patient_count) |> 
  sample_n(10000)
remove(flu_3)
# 2019
flu_4 <- read_parquet('tmp/flu_4_proc_p2.parquet')
flu_4 <- flu_4 |>
  dplyr::filter(month == 1 | month == 4 |
                  month == 7 | month == 10)
flu_balanced_19 <- flu_4 |> 
  group_by(flu, month) |> 
  uncount(patient_count) |> 
  sample_n(10000)
remove(flu_4)

# Combine balanced data
flu_v5 <- rbind(flu_balanced_16, flu_balanced_17, flu_balanced_18, flu_balanced_19)
saveRDS(flu_v5, 'tmp/flu_sample_v5.rds')

################################################################################
################################################################################
