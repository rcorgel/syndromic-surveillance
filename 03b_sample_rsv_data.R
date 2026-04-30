################################################################################
# File Name: 03b_sample_rsv_data                                               #
#                                                                              #
# Purpose:   Take random samples of visits from RSV medical claims.            #
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
rsv_1 <- read_parquet('tmp/rsv_1_proc_p2.parquet')
rsv_balanced_16 <- rsv_1 |> 
  group_by(rsv) |> 
  uncount(patient_count) |> 
  sample_n(166666)
remove(rsv_1)
# 2017
rsv_2 <- read_parquet('tmp/rsv_2_proc_p2.parquet')
rsv_balanced_17 <- rsv_2 |> 
  group_by(rsv) |> 
  uncount(patient_count) |> 
  sample_n(166667)
remove(rsv_2)
# 2018
rsv_3 <- read_parquet('tmp/rsv_3_proc_p2.parquet')
rsv_balanced_18 <- rsv_3 |> 
  group_by(rsv) |> 
  uncount(patient_count) |> 
  sample_n(166667)
remove(rsv_3)

# Combine balanced data
rsv_v1 <- rbind(rsv_balanced_16, rsv_balanced_17, rsv_balanced_18)
saveRDS(rsv_v1, 'tmp/rsv_sample_v1.rds')

########################
# VERSION 2: 2020-2023 #
########################

# Even sample of cases and non-cases 
# Sample equally across seasons, geography
# End with 500,000 cases and 500,000 non-cases

# 2020
rsv_5 <- read_parquet('tmp/rsv_5_proc_p2.parquet')
rsv_balanced_20 <- rsv_5 |> 
  group_by(rsv) |> 
  uncount(patient_count) |> 
  sample_n(166666)
remove(rsv_5)
# 2021
rsv_6 <- read_parquet('tmp/rsv_6_proc_p2.parquet')
rsv_balanced_21 <- rsv_6 |> 
  group_by(rsv) |> 
  uncount(patient_count) |> 
  sample_n(166667)
remove(rsv_6)
# 2022
rsv_7 <- read_parquet('tmp/rsv_7_proc_p2.parquet')
rsv_balanced_22 <- rsv_7 |> 
  group_by(rsv) |> 
  uncount(patient_count) |> 
  sample_n(166667)
remove(rsv_7)

# Combine balanced data
rsv_v2 <- rbind(rsv_balanced_20, rsv_balanced_21, rsv_balanced_22)
saveRDS(rsv_v2, 'tmp/rsv_sample_v2.rds')

################################################################################
################################################################################
