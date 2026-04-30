################################################################################
# File Name: 01b_load_rsv_data                                                 #
#                                                                              #
# Purpose:   Load RSV medical claims data from 2016-2023.                      #
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

# RSV symptom data

###############
# PART I DATA #
###############

# Part I data is at the month level and deduplicated (patients)

# Convert to parquet
rsv_16_17 <- read_csv_arrow('raw/rsv_2016-09_2017-08_1d_imputed.csv')
write_parquet(rsv_16_17, 'tmp/rsv_2016-09_2017-08_1d_imputed.parquet')
remove(rsv_16_17)

rsv_17_18 <- read_csv_arrow('raw/rsv_2017-09_2018-08_1d_imputed.csv')
write_parquet(rsv_17_18, 'tmp/rsv_2017-09_2018-08_1d_imputed.parquet')
remove(rsv_17_18)

rsv_18_19 <- read_csv_arrow('raw/rsv_2018-09_2019-08_1d_imputed.csv')
write_parquet(rsv_18_19, 'tmp/rsv_2018-09_2019-08_1d_imputed.parquet')
remove(rsv_18_19)

rsv_19_20 <- read_csv_arrow('raw/rsv_2019-09_2020-08_1d_imputed.csv')
write_parquet(rsv_19_20, 'tmp/rsv_2019-09_2020-08_1d_imputed.parquet')
remove(rsv_19_20)

rsv_20_21 <- read_csv_arrow('raw/rsv_2020-09_2021-08_1d_imputed.csv')
write_parquet(rsv_20_21, 'tmp/rsv_2020-09_2021-08_1d_imputed.parquet')
remove(rsv_20_21)

rsv_21_22 <- read_csv_arrow('raw/rsv_2021-09_2022-08_1d_imputed.csv')
write_parquet(rsv_21_22, 'tmp/rsv_2021-09_2022-08_1d_imputed.parquet')
remove(rsv_21_22)

rsv_22_23 <- read_csv_arrow('raw/rsv_2022-09_2023-08_1d_imputed.csv')
write_parquet(rsv_22_23, 'tmp/rsv_2022-09_2023-08_1d_imputed.parquet')
remove(rsv_22_23)

################
# PART II DATA #
################

# Part II data is at the week level and not deduplicated (visits)

# Convert to parquet
rsv_16_17 <- read_csv_arrow('raw/rsv_2016-09_2017-08_p2_imputed.csv')
write_parquet(rsv_16_17, 'tmp/rsv_2016-09_2017-08_p2_imputed.parquet')
remove(rsv_16_17)

rsv_17_18 <- read_csv_arrow('raw/rsv_2017-09_2018-08_p2_imputed.csv')
write_parquet(rsv_17_18, 'tmp/rsv_2017-09_2018-08_p2_imputed.parquet')
remove(rsv_17_18)

rsv_18_19 <- read_csv_arrow('raw/rsv_2018-09_2019-08_p2_imputed.csv')
write_parquet(rsv_18_19, 'tmp/rsv_2018-09_2019-08_p2_imputed.parquet')
remove(rsv_18_19)

rsv_19_20 <- read_csv_arrow('raw/rsv_2019-09_2020-08_p2_imputed.csv')
write_parquet(rsv_19_20, 'tmp/rsv_2019-09_2020-08_p2_imputed.parquet')
remove(rsv_19_20)

rsv_20_21 <- read_csv_arrow('raw/rsv_2020-09_2021-08_p2_imputed.csv')
write_parquet(rsv_20_21, 'tmp/rsv_2020-09_2021-08_p2_imputed.parquet')
remove(rsv_20_21)

rsv_21_22 <- read_csv_arrow('raw/rsv_2021-09_2022-08_p2_imputed.csv')
write_parquet(rsv_21_22, 'tmp/rsv_2021-09_2022-08_p2_imputed.parquet')
remove(rsv_21_22)

rsv_22_23 <- read_csv_arrow('raw/rsv_2022-09_2023-08_p2_imputed.csv')
write_parquet(rsv_22_23, 'tmp/rsv_2022-09_2023-08_p2_imputed.parquet')
remove(rsv_22_23)

################################################################################
################################################################################
