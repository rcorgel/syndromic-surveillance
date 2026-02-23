################################################################################
# File Name: explore_data_p2                                                   #
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
library(ISOweek)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

#######################
# 2. EXPLORE FLU DATA #
#######################

# Load data
flu_16_17 <- readRDS('tmp/flu_p2_dat_raw.rds')

# Create state and time variables
# 16-17
flu_16_17$state_fips <- substr(flu_16_17$county_fips, 1, 2)
# Convert week to date
flu_16_17$week_date <- sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", flu_16_17$year_week)
flu_16_17$week_date <- ISOweek2date(flu_16_17$week_date)
flu_16_17$month <- month(as.Date(flu_16_17$week_date))

saveRDS(flu_16_17, './tmp/flu_p2_dat_proc.rds')

flu_16_17 <- readRDS('./tmp/flu_p2_dat_proc.rds')

# Remove unknown states and PR
flu_16_17 <- flu_16_17 |> filter(state_fips != '99') |> filter(state_fips != '72')

saveRDS(flu_16_17, './tmp/flu_p2_dat_proc.rds')



# Remove unknown ages and genders
flu_16_17 <- flu_16_17 |> filter(age_grp != 'U')
flu_16_17 <- flu_16_17 |> filter(patient_gender_code != '')

saveRDS(flu_16_17, './tmp/flu_p2_dat_proc.rds')

# Impute <=5 to a random number between 1 and 5
# 16-17
flu_16_17$patient_count_imp <- flu_16_17$patient_count
flu_16_17$patient_count_imp[flu_16_17$patient_count_imp == '<=5'] <- 
  sample(1:5, length(flu_16_17$patient_count_imp[flu_16_17$patient_count_imp == '<=5']), replace= TRUE)
flu_16_17$patient_count_imp <- as.numeric(flu_16_17$patient_count_imp)

# Symptoms count
flu_16_17 <- flu_16_17 |> mutate(symp_count = fever + myalgia + cough + 
                                   sore_throat + short_breath + hypoxemia + 
                                   chest_pain + bronchitis + nausea_vom + 
                                   diarrhea + fatigue + headache + congestion + 
                                   sneezing)

saveRDS(flu_16_17, './tmp/flu_p2_dat_proc.rds')


flu_16_17 <- readRDS('./tmp/flu_p2_dat_proc.rds')
flu_16_17 <- flu_16_17 |> mutate(year = year(week_date)) |>
  filter(year < 2020)

week_symp <- flu_16_17 |> 
  filter(flu == 1) |>
  filter(symp_count > 0) |>
  group_by(week_date) |>
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
  distinct(week_date, fever, myalgia, cough, sore_throat, 
           short_breath, hypoxemia, chest_pain, bronchitis,
           nausea_vom, diarrhea, fatigue, headache,
           congestion, sneezing)

# Reshape data from wide to long format
week_symp_long <- melt(week_symp, id.vars = c("week_date"), variable.name = "symptom", 
                        value.name = "proportion")

ggplot() + 
  geom_line(data = week_symp_long, 
            aes(x = week_date, y = proportion, color = symptom)) + 
  theme_minimal()


# Symp Combos
symp_combos <- flu_16_17 |> 
  filter(flu == 1) |>
  filter(symp_count > 0) |>
  group_by(fever, myalgia, cough, sore_throat, 
           short_breath, hypoxemia, chest_pain, bronchitis,
           nausea_vom, diarrhea, fatigue, headache,
           congestion, sneezing) |>
  mutate(patient_count_imp = sum(patient_count_imp)) |>
  distinct(fever, myalgia, cough, sore_throat, 
           short_breath, hypoxemia, chest_pain, bronchitis,
           nausea_vom, diarrhea, fatigue, headache,
           congestion, sneezing, patient_count_imp) |>
  ungroup() |> mutate(percent = patient_count_imp / sum(patient_count_imp)) |>
  arrange(-percent) |>
  mutate(cumsum=cumsum(percent)) |>
  filter(cumsum < 0.80) |> 
  mutate(fever = ifelse(fever == 1, "fever", ""),
         myalgia = ifelse(myalgia == 1, "myalgia", ""),
         short_breath = ifelse(short_breath == 1, "short_breath", ""),
         hypoxemia = ifelse(hypoxemia == 1, "hypoxemia", ""),
         chest_pain = ifelse(chest_pain == 1, "chest_pain", ""),
         cough = ifelse(cough == 1, "cough", ""),
         bronchitis = ifelse(bronchitis == 1, "bronchitis", ""),
         nausea_vom = ifelse(nausea_vom == 1, "nausea_vom", ""),
         sore_throat = ifelse(sore_throat == 1, "sore_throat", ""),
         diarrhea = ifelse(diarrhea == 1, "diarrhea", ""),
         fatigue = ifelse(fatigue == 1, "fatigue", ""),
         headache = ifelse(headache == 1, "headache", ""),
         congestion = ifelse(congestion == 1, "congestion", ""),
         sneezing = ifelse(sneezing == 1, "sneezing", ""),
         combo = paste0(fever, ' ', myalgia, ' ',cough, ' ',sore_throat,' ', 
                          short_breath, ' ',hypoxemia, ' ',chest_pain, ' ',bronchitis,' ',
                          nausea_vom, ' ',diarrhea, ' ',fatigue, ' ',headache,' ',
                          congestion, ' ',sneezing)) |>
  mutate(combo = gsub("\\s+" ," + ", str_trim(combo))) |>
  select(c(combo, percent))
  
  
symp_combos_week <- flu_16_17 |> 
  filter(flu == 1) |>
  filter(symp_count > 0) |>
  group_by(week_date, fever, myalgia, cough, sore_throat, 
           short_breath, hypoxemia, chest_pain, bronchitis,
           nausea_vom, diarrhea, fatigue, headache,
           congestion, sneezing) |>
  mutate(patient_count_imp = sum(patient_count_imp)) |>
  distinct(week_date, fever, myalgia, cough, sore_throat, 
           short_breath, hypoxemia, chest_pain, bronchitis,
           nausea_vom, diarrhea, fatigue, headache,
           congestion, sneezing, patient_count_imp) |>
  ungroup() |> mutate(percent = patient_count_imp / sum(patient_count_imp)) |>
  mutate(fever = ifelse(fever == 1, "fever", ""),
         myalgia = ifelse(myalgia == 1, "myalgia", ""),
         short_breath = ifelse(short_breath == 1, "short_breath", ""),
         hypoxemia = ifelse(hypoxemia == 1, "hypoxemia", ""),
         chest_pain = ifelse(chest_pain == 1, "chest_pain", ""),
         cough = ifelse(cough == 1, "cough", ""),
         bronchitis = ifelse(bronchitis == 1, "bronchitis", ""),
         nausea_vom = ifelse(nausea_vom == 1, "nausea_vom", ""),
         sore_throat = ifelse(sore_throat == 1, "sore_throat", ""),
         diarrhea = ifelse(diarrhea == 1, "diarrhea", ""),
         fatigue = ifelse(fatigue == 1, "fatigue", ""),
         headache = ifelse(headache == 1, "headache", ""),
         congestion = ifelse(congestion == 1, "congestion", ""),
         sneezing = ifelse(sneezing == 1, "sneezing", ""),
         combo = paste0(fever, ' ', myalgia, ' ',cough, ' ',sore_throat,' ', 
                        short_breath, ' ',hypoxemia, ' ',chest_pain, ' ',bronchitis,' ',
                        nausea_vom, ' ',diarrhea, ' ',fatigue, ' ',headache,' ',
                        congestion, ' ',sneezing)) |>
  mutate(combo = gsub("\\s+" ," + ", str_trim(combo))) |>
  select(c(week_date, patient_count_imp, combo)) |>
  group_by(week_date) |>
  mutate(percent = patient_count_imp / sum(patient_count_imp)) |> ungroup() |> 
  group_by(combo) |>
  mutate(Z = (percent - mean(percent)) / sd(percent),
         Z_count = (patient_count_imp - mean(patient_count_imp)) / sd(patient_count_imp))

combo_week_flu <- right_join(symp_combos_week, symp_combos, by = c('combo' = 'combo'))
#library("modelbased")


combo_week_flu <- combo_week_flu |> filter(week_date <= as.Date('2019-08-01'))
ggplot() + 
  geom_line(data = combo_week, 
            aes(x = week_date, y = Z_count, color = combo)) + 
  theme_minimal()

ggplot() + 
  geom_line(data = combo_week[combo_week$percent.y > 0.001,], 
            aes(x = week_date, y = patient_count_imp, color = combo)) + 
  theme_minimal()

ggplot() + 
  geom_line(data = combo_week_flu, 
            aes(x = week_date, y = Z_count, color = combo)) + 
  theme_minimal()

ggplot() + 
  geom_line(data = combo_week_flu[combo_week_flu$percent.y > 0.001,], 
            aes(x = week_date, y = patient_count_imp, color = combo)) + 
  theme_minimal()






eeaflu_descriptives_symp <- flu |> mutate(
                                       fever = ifelse(fever == 1, "fever", ""),
                                       myalgia = ifelse(myalgia == 1, "myalgia", ""),
                                       short_breath = ifelse(short_breath == 1, "short_breath", ""),
                                       hypoxemia = ifelse(hypoxemia == 1, "hypoxemia", ""),
                                       chest_pain = ifelse(chest_pain == 1, "chest_pain", ""),
                                       cough = ifelse(cough == 1, "cough", ""),
                                       bronchitis = ifelse(bronchitis == 1, "bronchitis", ""),
                                       nausea_vom = ifelse(nausea_vom == 1, "nausea_vom", ""),
                                       sore_throat = ifelse(sore_throat == 1, "sore_throat", ""),
                                       diarrhea = ifelse(diarrhea == 1, "diarrhea", ""),
                                       fatigue = ifelse(fatigue == 1, "fatigue", ""),
                                       headache = ifelse(headache == 1, "headache", ""),
                                       congestion = ifelse(congestion == 1, "congestion", ""),
                                       sneezing = ifelse(sneezing == 1, "sneezing", "")) |>



symp_combos <- flu_16_17 |> 
  filter(flu == 1) |>
  filter(symp_count > 0) |>
  group_by(fever, myalgia, cough, sore_throat, 
           short_breath, hypoxemia, chest_pain, bronchitis,
           nausea_vom, diarrhea, fatigue, headache,
           congestion, sneezing) |>
  mutate(patient_count_imp = sum(patient_count_imp)) |>
  distinct(fever, myalgia, cough, sore_throat, 
           short_breath, hypoxemia, chest_pain, bronchitis,
           nausea_vom, diarrhea, fatigue, headache,
           congestion, sneezing, patient_count_imp) |>
  ungroup() |> mutate(percent = patient_count_imp / sum(patient_count_imp))















root <- "https://www.google.com/search?q="
reprex_df %>% 
  mutate(new_col = paste0(root, var1, "+", var2, "+", var3))





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