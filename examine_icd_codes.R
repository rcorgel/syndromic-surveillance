################################################################################
# File Name: examine_icd_codes                                                 #
#                                                                              #
# Purpose:   Examine most popular icd diagnostic codes for flu and COVID-19.   #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load flu data                                                  #
#            3. Load COVID-19 data                                             #
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
library(readxl)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

###############
# 2. FLU DATA #
###############

# Load flu codes
flu_codes <- read.csv('./raw/flu_diag_codes_freq.csv')

# Load all codes with descriptions and append
# ICD-10
icd_10_codes <- read_xlsx('./raw/Section111ValidICD10-Jan2024.xlsx')
icd_10_codes <- icd_10_codes[, c(1, 3)]
icd_10_codes$type <- 'ICD-10'
icd_10_codes <- icd_10_codes |> rename('code' = 'CODE',
                                       'desc' = 'LONG DESCRIPTION (VALID ICD-10 FY2024)')
# ICD-9
icd_9_codes <- read_xlsx('./raw/Section111ValidICD9-Jan2024.xlsx')
icd_9_codes <- icd_9_codes[, c(1, 2)]
icd_9_codes$type <- 'ICD-9'
icd_9_codes <- icd_9_codes |> rename('code' = 'CODE',
                                       'desc' = 'LONG DESCRIPTION (VALID ICD-9 FY2024)')
# Append
all_icd_codes <- rbind(icd_10_codes, icd_9_codes)
all_icd_codes <- all_icd_codes[!duplicated(all_icd_codes$code), ] # deduplicate

# Merge on descriptions
flu_codes <- left_join(flu_codes, all_icd_codes, by = c('diagnosis_codes' = 'code'))

# Save
write.csv(flu_codes, './tmp/flu_diag_codes.csv')

####################
# 3. COVID-19 DATA #
####################

# Load COVID-19 codes
covid_codes <- read.csv('./raw/covid_diag_codes_freq.csv')

# Merge on descriptions
covid_codes <- left_join(covid_codes, all_icd_codes, by = c('diagnosis_codes' = 'code'))

# Save
write.csv(covid_codes, './tmp/covid_diag_codes.csv')




