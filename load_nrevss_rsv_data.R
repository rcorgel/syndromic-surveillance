################################################################################
# File Name: load_nrevss)rsv_data                                              #
#                                                                              #
# Purpose:   Load NREVSS RSV data from the CDC.                                #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load data                                                      #
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
library(ISOweek)
library(zoo)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

#########################
# 2. LOAD NATIONAL DATA #
#########################

# Load RSV NREVSS Data
rsv_nrevss <- read.csv('./raw/Respiratory_Syncytial_Virus_Laboratory_Data__NREVSS__20250721.csv', header = TRUE, sep = ',')

# Convert to week date
rsv_nrevss$week_date <- as.Date(rsv_nrevss$Week.ending.Date, format = "%d %b %Y")
rsv_nrevss$week_date <- rsv_nrevss$date - 5
