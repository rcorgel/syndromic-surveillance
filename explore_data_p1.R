################################################################################
# File Name: explore_data_p1                                                   #
#                                                                              #
# Purpose:   Explore Influenza, COVID-19, and RSV medical claims data.         #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Explore Influenza data                                         #
#            3. Explore COVID-19 data                                          #
#            3. Explore RSV data                                               #
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
library(reshape2)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

#######################
# 2. EXPLORE FLU DATA #
#######################

# Load data
flu <- readRDS('./tmp/flu_p1_dat_proc.rds')

month_symp <- flu |> 
  filter(flu == 1) |>
  group_by(month_date) |>
  mutate(fever = sum(fever*patient_count_imp) / sum(patient_count_imp),
         myalgia = sum(myalgia*patient_count_imp) / sum(patient_count_imp),
         cough = sum(cough*patient_count_imp) / sum(patient_count_imp),
         sore_throat = sum(sore_throat*patient_count_imp) / sum(patient_count_imp), 
         short_breath = sum(short_breath*patient_count_imp) / sum(patient_count_imp), 
         hypoxemia = sum(hypoxemia*patient_count_imp) / sum(patient_count_imp), 
         chest_pain = sum(chest_pain*patient_count_imp) / sum(patient_count_imp), 
         bronchitis = sum(bronchitis*patient_count_imp) / sum(patient_count_imp),
         nausea_vom = sum(nausea_vom*patient_count_imp) / sum(patient_count_imp), 
         diarrhea = sum(diarrhea*patient_count_imp) / sum(patient_count_imp), 
         fatigue = sum(fatigue*patient_count_imp) / sum(patient_count_imp), 
         headache = sum(headache*patient_count_imp) / sum(patient_count_imp),
         congestion = sum(congestion*patient_count_imp) / sum(patient_count_imp), 
         sneezing = sum(sneezing*patient_count_imp) / sum(patient_count_imp)) |>
  distinct(fever, myalgia, cough, sore_throat, 
           short_breath, hypoxemia, chest_pain, bronchitis,
           nausea_vom, diarrhea, fatigue, headache,
           congestion, sneezing)

# Reshape data from wide to long format
month_symp_long <- melt(month_symp, id.vars = c("month_date"), variable.name = "symptom", 
                value.name = "proportion")

ggplot() + 
  geom_line(data = month_symp_long, 
            aes(x = month_date, y = proportion, color = symptom)) + 
  theme_minimal()

month_symp_count <- flu |> 
  #filter(flu == 1) |>
  group_by(month_date) |>
  mutate(fever = sum(fever*patient_count_imp),
         myalgia = sum(myalgia*patient_count_imp),
         cough = sum(cough*patient_count_imp),
         sore_throat = sum(sore_throat*patient_count_imp), 
         short_breath = sum(short_breath*patient_count_imp), 
         hypoxemia = sum(hypoxemia*patient_count_imp), 
         chest_pain = sum(chest_pain*patient_count_imp), 
         bronchitis = sum(bronchitis*patient_count_imp),
         nausea_vom = sum(nausea_vom*patient_count_imp), 
         diarrhea = sum(diarrhea*patient_count_imp), 
         fatigue = sum(fatigue*patient_count_imp), 
         headache = sum(headache*patient_count_imp),
         congestion = sum(congestion*patient_count_imp), 
         sneezing = sum(sneezing*patient_count_imp)) |>
  distinct(fever, myalgia, cough, sore_throat, 
           short_breath, hypoxemia, chest_pain, bronchitis,
           nausea_vom, diarrhea, fatigue, headache,
           congestion, sneezing)

# Reshape data from wide to long format
month_symp_long_count <- melt(month_symp_count, id.vars = c("month_date"), variable.name = "symptom", 
                        value.name = "count")

ggplot() + 
  geom_line(data = month_symp_long_count, 
            aes(x = month_date, y = count, color = symptom)) + 
  theme_minimal()
