################################################################################
# File Name: 01a_load_flu_data                                                 #
#                                                                              #
# Purpose:   Load influenza medical claims data from 2016-2023.                #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load data and save to parquet                                  #
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

####################################
# 2. LOAD DATA AND SAVE TO PARQUET #
####################################

# Influenza symptom data

###############
# PART I DATA #
###############

# Part I data is at the month level and deduplicated (patients)

# Convert to parquet
flu_16_17 <- read_csv_arrow('raw/flu_2016-09_2017-08_1b_imputed.csv')
write_parquet(flu_16_17, 'tmp/flu_2016-09_2017-08_1b_imputed.parquet')
remove(flu_16_17)

flu_17_18 <- read_csv_arrow('raw/flu_2017-09_2018-08_1b_imputed.csv')
write_parquet(flu_17_18, 'tmp/flu_2017-09_2018-08_1b_imputed.parquet')
remove(flu_17_18)

flu_18_19 <- read_csv_arrow('raw/flu_2018-09_2019-08_1b_imputed.csv')
write_parquet(flu_18_19, 'tmp/flu_2018-09_2019-08_1b_imputed.parquet')
remove(flu_18_19)

flu_19_20 <- read_csv_arrow('raw/flu_2019-09_2020-08_1b_imputed.csv')
write_parquet(flu_19_20, 'tmp/flu_2019-09_2020-08_1b_imputed.parquet')
remove(flu_19_20)

flu_20_21 <- read_csv_arrow('raw/flu_2020-09_2021-08_1b_imputed.csv')
write_parquet(flu_20_21, 'tmp/flu_2020-09_2021-08_1b_imputed.parquet')
remove(flu_20_21)

flu_21_22 <- read_csv_arrow('raw/flu_2021-09_2022-08_1b_imputed.csv')
write_parquet(flu_21_22, 'tmp/flu_2021-09_2022-08_1b_imputed.parquet')
remove(flu_21_22)

flu_22_23 <- read_csv_arrow('raw/flu_2022-09_2023-08_1b_imputed.csv')
write_parquet(flu_22_23, 'tmp/flu_2022-09_2023-08_1b_imputed.parquet')
remove(flu_22_23)

################
# PART II DATA #
################

# Part II data is at the week level and not deduplicated (visits)

# Convert to parquet
flu_16_17 <- read_csv_arrow('raw/flu_2016-09_2017-08_p2_imputed.csv')
write_parquet(flu_16_17, 'tmp/flu_2016-09_2017-08_p2_imputed.parquet')
remove(flu_16_17)

flu_17_18 <- read_csv_arrow('raw/flu_2017-09_2018-08_p2_imputed.csv')
write_parquet(flu_17_18, 'tmp/flu_2017-09_2018-08_p2_imputed.parquet')
remove(flu_17_18)

flu_18_19 <- read_csv_arrow('raw/flu_2018-09_2019-08_p2_imputed.csv')
write_parquet(flu_18_19, 'tmp/flu_2018-09_2019-08_p2_imputed.parquet')
remove(flu_18_19)

flu_19_20 <- read_csv_arrow('raw/flu_2019-09_2020-08_p2_imputed.csv')
write_parquet(flu_19_20, 'tmp/flu_2019-09_2020-08_p2_imputed.parquet')
remove(flu_19_20)

flu_20_21 <- read_csv_arrow('raw/flu_2020-09_2021-08_p2_imputed.csv')
write_parquet(flu_20_21, 'tmp/flu_2020-09_2021-08_p2_imputed.parquet')
remove(flu_20_21)

flu_21_22 <- read_csv_arrow('raw/flu_2021-09_2022-08_p2_imputed.csv')
write_parquet(flu_21_22, 'tmp/flu_2021-09_2022-08_p2_imputed.parquet')
remove(flu_21_22)

flu_22_23 <- read_csv_arrow('raw/flu_2022-09_2023-08_p2_imputed.csv')
write_parquet(flu_22_23, 'tmp/flu_2022-09_2023-08_p2_imputed.parquet')
remove(flu_22_23)

################################################################################
################################################################################
