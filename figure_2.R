################################################################################
# File Name: figure_2                                                          #
#                                                                              #
# Purpose:   Create figure 2 for the manuscript.                               #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Create figure 2                                                #
#               a. Create profile plots                                        #
#               b. Create age forest plots                                     #
#               c. Create time forest plots                                    #
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
library(ISOweek)
library(splus2R)
library(scales)
library(sf)
library(geosphere)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

######################
# 2. CREATE FIGURE 2 #
######################

###########################
# A. CREATE PROFILE PLOTS #
###########################

# Load flu data
flu_dat_p1 <- readRDS('./tmp/flu_p1_dat_proc_symp.rds')

# Add month
flu_dat_p1$month <- month(flu_dat_p1$month_date)

# Calculate Percents
flu_symp <- flu_dat_p1 |> 
  group_by(flu) |>
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
  distinct(flu, fever, myalgia, cough, sore_throat, 
           short_breath, hypoxemia, chest_pain, bronchitis,
           nausea_vom, diarrhea, fatigue, headache,
           congestion, sneezing)

# Convert to long
flu_symp_long <- flu_symp |>
  pivot_longer(!c(flu), names_to = "symptom", values_to = "percent") |>
  mutate(symptom = str_to_title(symptom))

# Change symptom names
flu_symp_long$symptom <- ifelse(flu_symp_long$symptom == 'Sore_throat', 
                                'Sore Throat', flu_symp_long$symptom)
flu_symp_long$symptom <- ifelse(flu_symp_long$symptom == 'Nausea_vom', 
                                'Nausea', flu_symp_long$symptom)
flu_symp_long$symptom <- ifelse(flu_symp_long$symptom == 'Chest_pain', 
                                'Chest Pain', flu_symp_long$symptom)
flu_symp_long$symptom <- ifelse(flu_symp_long$symptom == 'Short_breath', 
                                'Shortness of Breath', flu_symp_long$symptom)

# Change case names
flu_symp_long$flu <- ifelse(flu_symp_long$flu == 0, 
                            'Non-Influenza Cases', flu_symp_long$flu)
flu_symp_long$flu <- ifelse(flu_symp_long$flu == 1, 
                            'Influenza Cases', flu_symp_long$flu)

# Merge on rank
flu_symp_long_flu <- flu_symp_long[flu_symp_long$flu == 'Influenza Cases', c(2, 3)]
flu_symp_long_flu$rank <- rank(flu_symp_long_flu$percent, na.last = TRUE)
flu_symp_long <- left_join(flu_symp_long, flu_symp_long_flu[, c(1, 3)], by = 'symptom')

# Load rsv data
rsv_dat_p1 <- readRDS('./tmp/rsv_p1_dat_proc_symp.rds')

# Calculate Percents
rsv_symp <- rsv_dat_p1 |> 
  group_by(rsv) |>
  mutate(fever = sum(fever*patient_count_imp) / sum(patient_count_imp),
         cough = sum(cough*patient_count_imp) / sum(patient_count_imp),
         sore_throat = sum(sore_throat*patient_count_imp) / sum(patient_count_imp), 
         short_breath = sum(short_breath*patient_count_imp) / sum(patient_count_imp), 
         hypoxemia = sum(hypoxemia*patient_count_imp) / sum(patient_count_imp), 
         bronchitis = sum(bronchitis*patient_count_imp) / sum(patient_count_imp),
         loss_appetite = sum(loss_appetite*patient_count_imp) / sum(patient_count_imp), 
         fatigue = sum(fatigue*patient_count_imp) / sum(patient_count_imp), 
         headache = sum(headache*patient_count_imp) / sum(patient_count_imp),
         congestion = sum(congestion*patient_count_imp) / sum(patient_count_imp), 
         sneezing = sum(sneezing*patient_count_imp) / sum(patient_count_imp)) |>
  distinct(rsv, fever, cough, sore_throat, 
           short_breath, hypoxemia, bronchitis,
           loss_appetite, fatigue, headache,
           congestion, sneezing)

# Convert to long
rsv_symp_long <- rsv_symp |>
  pivot_longer(!c(rsv), names_to = "symptom", values_to = "percent") |>
  mutate(symptom = str_to_title(symptom))

# Change symptom names
rsv_symp_long$symptom <- ifelse(rsv_symp_long$symptom == 'Sore_throat', 
                                'Sore Throat', rsv_symp_long$symptom)
rsv_symp_long$symptom <- ifelse(rsv_symp_long$symptom == 'Loss_appetite', 
                                'Loss of Appetite', rsv_symp_long$symptom)
rsv_symp_long$symptom <- ifelse(rsv_symp_long$symptom == 'Short_breath', 
                                'Shortness of Breath', rsv_symp_long$symptom)

# Change case names
rsv_symp_long$rsv <- ifelse(rsv_symp_long$rsv == 0, 
                            'Non-RSV Cases', rsv_symp_long$rsv)
rsv_symp_long$rsv <- ifelse(rsv_symp_long$rsv == 1, 
                            'RSV Cases', rsv_symp_long$rsv)

# Merge on rank
rsv_symp_long_rsv <- rsv_symp_long[rsv_symp_long$rsv == 'RSV Cases', c(2, 3)]
rsv_symp_long_rsv$rank <- rank(rsv_symp_long_rsv$percent, na.last = TRUE)
rsv_symp_long <- left_join(rsv_symp_long, rsv_symp_long_rsv[, c(1, 3)], by = 'symptom')

# Set plot theme
percent_theme <- theme(legend.position = "bottom",
                       axis.text = element_text(size=14),
                       axis.title.y = element_text(size=16),
                       axis.title.x = element_blank(),
                       strip.background = element_blank(),
                       plot.title = element_text(size=18),
                       axis.text.x = element_text(angle = 45, vjust = 1, 
                                                  hjust=1, color = 'black'),
                       legend.title = element_text(size = 16),
                       legend.text = element_text(size = 14))

# Set colors for influenza
flu_colors <- c(
  "Influenza Cases" = "#d7642c",
  "Non-Influenza Cases" = "darkgray")

# Create flight
percent_flu <- ggplot(flu_symp_long, 
                      aes(x = fct_reorder(symptom, desc(rank)), 
                          y = percent)) +
  geom_line(aes(color = flu, lty = flu, group = flu), linewidth = 2, alpha = 0.6) +  
  geom_point(aes(color = flu), size = 4, show.legend = TRUE, alpha = 0.8) +
  theme_minimal() + percent_theme +
  labs(title = "Influenza Syndromic Profile\n") +
  ylab('Proportion of Patients') + xlab('') + expand_limits(y = c(0, 0.5)) +
  scale_color_manual(values = flu_colors)+
  scale_linetype_manual(values = c('solid', 'dashed'))+
  labs(color = "", lty = "") 

# Display
percent_flu

# Set colors for RSV
rsv_colors <- c(
  "RSV Cases" = "#41afaa",
  "Non-RSV Cases" = "darkgray")

# Create flight
percent_rsv <- ggplot(rsv_symp_long, 
                      aes(x = fct_reorder(symptom, desc(rank)), 
                          y = percent)) +
  geom_line(aes(color = rsv, lty = rsv, group = rsv), linewidth = 2, alpha = 0.6) +  
  geom_point(aes(color = rsv), size = 4, show.legend = TRUE, alpha = 0.8) +
  theme_minimal() + percent_theme +
  labs(title = "RSV Syndromic Profile\n") +
  ylab('Proportion of Patients') + xlab('') + expand_limits(y = c(0, 0.5)) +
  scale_color_manual(values = rsv_colors,
                     breaks=c("RSV Cases", "Non-RSV Cases"))+
  scale_linetype_manual(values = c('solid', 'dashed'),
                        breaks=c("RSV Cases", "Non-RSV Cases"))+
  labs(color = "", lty = "") 

# Display
percent_rsv

############################
# CREATE AGE FORREST PLOTS #
############################

#############
# INFLUENZA #
#############

# Add time variables to Part 1 data
flu_dat_p1$month <- month(flu_dat_p1$month_date)
flu_dat_p1$year <- year(flu_dat_p1$month_date)

# Load maps
usa_albers_state <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_state.rds'))   # convert to sf
usa_albers_county <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_county.rds')) # convert to sf

# Load x walk
x_walk <- readRDS('./tmp/county_group_xwalk.rds')

# Collapse county map to county group
usa_albers_county_x_walk <- left_join(x_walk, usa_albers_county[, c(5, 6)], by = c('fips' = 'GEOID'))
usa_albers_county_group <- usa_albers_county_x_walk %>% 
  group_by(county_fips) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()

# Calculate centroid
usa_albers_county_group$centroid <- st_centroid(usa_albers_county_group$geometry)

# Convert to lat/long
usa_albers_county_group$latlon <- st_transform(usa_albers_county_group$centroid, crs = 4326) # EPSG:4326 for WGS84
usa_albers_county_group <- usa_albers_county_group |>
  group_by(county_fips) |>
  mutate(south = ifelse(latlon[[1]][2] <= 30, 1, 0),
         state_fips = substr(county_fips, 1, 2),
         south = ifelse(state_fips == '02', 0, south))

# Load Urban/Rural classification
urban <- read.csv('raw/urban_rural_cat.csv', header = T)
urban <- urban |> 
  mutate(county_fips = ifelse(nchar(as.character(FIPS.code)) == 4,
                              paste0('0', as.character(FIPS.code)),
                              as.character(FIPS.code)),
         urban_code = as.numeric(substr(X2023.Code, 1, 1)))
urban_merge <- left_join(x_walk, urban, by = c('fips' = 'county_fips'))
urban_merge <- urban_merge |> group_by(county_fips) |>
  mutate(urban_code_binary = ifelse(round(mean(urban_code)) < 5, 1, 0)) |>
  distinct(county_fips, urban_code_binary)

# Deal w/ CT new FIPS
urban_merge$urban_code_binary <- ifelse(is.na(urban_merge$urban_code_binary), 0,
                                        urban_merge$urban_code_binary)
urban_merge[urban_merge$county_fips == '09003',]$urban_code_binary <- 1
urban_merge[urban_merge$county_fips == '09013',]$urban_code_binary <- 1
urban_merge[urban_merge$county_fips == '09007',]$urban_code_binary <- 1

# Merge onto flu data
flu_dat_p1 <- left_join(flu_dat_p1, urban_merge, by = 'county_fips')
flu_dat_p1 <- left_join(flu_dat_p1, usa_albers_county_group[, c(1, 5)], by = 'county_fips')

# Save
urban_south_binary <- left_join(urban_merge, usa_albers_county_group[, c(1, 5)], by = 'county_fips')
saveRDS(urban_south_binary, "./tmp/urban_south_binary.rds", )

# Create age group variable
# Part 1
flu_dat_p1$age_group <- factor(as.numeric(flu_dat_p1$age_grp),
                               levels = c(0, 1, 2, 3, 4, 5),
                               labels = c('0-4', '5-12', '13-17', '18-49', '50-64', '65+'))

# Load population data
county_pop <- readRDS('tmp/county_group_pop.rds')

# Create Season Variable
flu_dat_p1$Season <- ifelse(flu_dat_p1$year == 2016, '2016-2017', '')
flu_dat_p1$Season <- ifelse(flu_dat_p1$year == 2017 & flu_dat_p1$month < 9, '2016-2017', flu_dat_p1$Season)
flu_dat_p1$Season <- ifelse(flu_dat_p1$year == 2017 & flu_dat_p1$month > 8, '2017-2018', flu_dat_p1$Season)
flu_dat_p1$Season <- ifelse(flu_dat_p1$year == 2018 & flu_dat_p1$month < 9, '2017-2018', flu_dat_p1$Season)
flu_dat_p1$Season <- ifelse(flu_dat_p1$year == 2018 & flu_dat_p1$month > 8, '2018-2019', flu_dat_p1$Season)
flu_dat_p1$Season <- ifelse(flu_dat_p1$year == 2019 & flu_dat_p1$month < 9, '2018-2019', flu_dat_p1$Season)
flu_dat_p1$Season <- ifelse(flu_dat_p1$year == 2019 & flu_dat_p1$month > 8, '2019-2020', flu_dat_p1$Season)
flu_dat_p1$Season <- ifelse(flu_dat_p1$year == 2020, '2019-2020', flu_dat_p1$Season)

# Restrict to counties with 10,000 or more individuals
flu_dat_p1 <- left_join(flu_dat_p1, county_pop, 
                        by = c('county_fips' = 'county_fips'))
flu_dat_p1 <- flu_dat_p1 |> filter(county_group_pop > 10000)

# Make data individual level (sample 1 million individuals)
# Balance flu cases and non-flu cases
balanced <- flu_dat_p1 |> group_by(flu) |>
  uncount(patient_count_imp) |> sample_n(500000)

# Split the data randomly into testing and training
sample <- sample(c(TRUE, FALSE), nrow(balanced), replace=TRUE, prob=c(0.5, 0.5))
flu_train  <- balanced[sample, ]
flu_test   <- balanced[!sample, ]

## Define candidate models ##

# Fever only model
fev_model <- glm(flu ~ fever,
                 family = binomial(link = "logit"),
                 data = flu_train)
summary(fev_model)

# Save
saveRDS(fev_model, "./tmp/fever_model.rds")

# All symptoms model
symp_model <- glm(flu ~ fever + 
                    cough + 
                    sore_throat +
                    myalgia + 
                    short_breath + 
                    hypoxemia + 
                    nausea_vom + 
                    bronchitis + 
                    chest_pain + 
                    diarrhea + 
                    fatigue + 
                    headache + 
                    congestion + 
                    sneezing,
                  family = binomial(link = "logit"),
                  data = flu_train)
summary(symp_model)

# Save
saveRDS(symp_model, "./tmp/symptom_model.rds")

# Age model
age_model <- glm(flu ~ fever*as.factor(age_group) + 
                   cough*as.factor(age_group) + 
                   sore_throat*as.factor(age_group) +
                   myalgia*as.factor(age_group) + 
                   short_breath*as.factor(age_group) + 
                   hypoxemia*as.factor(age_group) + 
                   nausea_vom*as.factor(age_group) + 
                   bronchitis*as.factor(age_group) + 
                   chest_pain*as.factor(age_group) + 
                   diarrhea*as.factor(age_group) + 
                   fatigue*as.factor(age_group) + 
                   headache*as.factor(age_group) + 
                   congestion*as.factor(age_group) + 
                   sneezing*as.factor(age_group),
                 family = binomial(link = "logit"),
                 data = flu_train)
summary(age_model)

# Save
saveRDS(age_model, "./tmp/age_model.rds")

# Time model
time_model <- glm(flu ~ fever*as.factor(age_group) + 
                    cough*as.factor(age_group) + 
                    sore_throat*as.factor(age_group) +
                    myalgia*as.factor(age_group) + 
                    short_breath*as.factor(age_group) + 
                    hypoxemia*as.factor(age_group) + 
                    nausea_vom*as.factor(age_group) + 
                    bronchitis*as.factor(age_group) + 
                    chest_pain*as.factor(age_group) + 
                    diarrhea*as.factor(age_group) + 
                    fatigue*as.factor(age_group) + 
                    headache*as.factor(age_group) + 
                    congestion*as.factor(age_group) + 
                    sneezing*as.factor(age_group) +
                    as.factor(month),
                  family = binomial(link = "logit"),
                  data = flu_train)
summary(time_model)

# Save
saveRDS(time_model, "./tmp/time_model.rds")

# Geography model
geo_model <- glm(flu ~ fever*as.factor(age_group) + 
                   cough*as.factor(age_group) + 
                   sore_throat*as.factor(age_group) +
                   myalgia*as.factor(age_group) + 
                   short_breath*as.factor(age_group) + 
                   hypoxemia*as.factor(age_group) + 
                   nausea_vom*as.factor(age_group) + 
                   bronchitis*as.factor(age_group) + 
                   chest_pain*as.factor(age_group) + 
                   diarrhea*as.factor(age_group) + 
                   fatigue*as.factor(age_group) + 
                   headache*as.factor(age_group) + 
                   congestion*as.factor(age_group) + 
                   sneezing*as.factor(age_group) +
                   as.factor(month)*urban_code_binary +
                   as.factor(month)*south,
                 family = binomial(link = "logit"),
                 data = flu_train)
summary(geo_model)

# Save
saveRDS(geo_model, "./tmp/geography_model.rds")

# Season model
season_model <- glm(flu ~ fever*as.factor(Season) + 
                   cough*as.factor(Season) + 
                   sore_throat*as.factor(Season) +
                   myalgia*as.factor(Season) + 
                   short_breath*as.factor(Season) + 
                   hypoxemia*as.factor(Season) + 
                   nausea_vom*as.factor(Season) + 
                   bronchitis*as.factor(Season) + 
                   chest_pain*as.factor(Season) + 
                   diarrhea*as.factor(Season) + 
                   fatigue*as.factor(Season) + 
                   headache*as.factor(Season) + 
                   congestion*as.factor(Season) + 
                   sneezing*as.factor(Season) +
                   as.factor(month)*as.factor(Season),
                 family = binomial(link = "logit"),
                 data = flu_train)
summary(season_model)

# Save
saveRDS(season_model, "./tmp/season_model.rds")

# Quick model tests
flu_test$pred <- predict(fev_model, newdata = flu_test, type = "response")
pROC::auc(pROC::roc(flu_test$flu, flu_test$pred))
ggplot(flu_test) + geom_density(aes(x = pred, group = as.factor(flu), color = as.factor(flu)))

flu_test$pred <- predict(symp_model, newdata = flu_test, type = "response")
pROC::auc(pROC::roc(flu_test$flu, flu_test$pred))
ggplot(flu_test) + geom_density(aes(x = pred, group = as.factor(flu), color = as.factor(flu)))

flu_test$pred <- predict(age_model, newdata = flu_test, type = "response")
pROC::auc(pROC::roc(flu_test$flu, flu_test$pred))
ggplot(flu_test) + geom_density(aes(x = pred, group = as.factor(flu), color = as.factor(flu)))

flu_test$pred <- predict(time_model, newdata = flu_test, type = "response")
pROC::auc(pROC::roc(flu_test$flu, flu_test$pred))
ggplot(flu_test) + geom_density(aes(x = pred, group = as.factor(flu), color = as.factor(flu)))

flu_test$pred <- predict(geo_model, newdata = flu_test, type = "response")
pROC::auc(pROC::roc(flu_test$flu, flu_test$pred))
ggplot(flu_test) + geom_density(aes(x = pred, group = as.factor(flu), color = as.factor(flu)))

# Determine the cut off value
cp <- cutpointr(flu_test, pred, flu, 
                method = maximize_metric,
                metric = youden)
summary(cp)

# Overall model with simple time
overall_model <- glm(flu ~ fever + 
                   cough + 
                   sore_throat +
                   myalgia + 
                   short_breath + 
                   hypoxemia + 
                   nausea_vom + 
                   bronchitis + 
                   chest_pain + 
                   diarrhea + 
                   fatigue + 
                   headache + 
                   congestion + 
                   sneezing,
                 family = binomial(link = "logit"),
                 data = flu_train)
summary(overall_model)

# Model by month and age group
# all_model <- glm(flu ~ fever*as.factor(age_group) + 
#                    fever*as.factor(month) +
#                    cough*as.factor(age_group) + 
#                    cough*as.factor(month) +
#                    sore_throat*as.factor(age_group) +
#                    sore_throat*as.factor(month) +
#                    myalgia*as.factor(age_group) + 
#                    myalgia*as.factor(month) +
#                    short_breath*as.factor(age_group) + 
#                    short_breath*as.factor(month) +
#                    hypoxemia*as.factor(age_group) + 
#                    hypoxemia*as.factor(month) +
#                    nausea_vom*as.factor(age_group) + 
#                    nausea_vom*as.factor(month) +
#                    bronchitis*as.factor(age_group) + 
#                    bronchitis*as.factor(month) +
#                    chest_pain*as.factor(age_group) + 
#                    chest_pain*as.factor(month) +
#                    diarrhea*as.factor(age_group) + 
#                    diarrhea*as.factor(month) +
#                    fatigue*as.factor(age_group) + 
#                    fatigue*as.factor(month) +
#                    headache*as.factor(age_group) + 
#                    headache*as.factor(month) +
#                    congestion*as.factor(age_group) + 
#                    congestion*as.factor(month) +
#                    sneezing*as.factor(age_group) +
#                    sneezing*as.factor(month) +
#                    as.factor(month)*south + 
#                    as.factor(month)*urban_code_binary,
#                  family = binomial(link = "logit"),
#                  data = flu_train)
# summary(all_model)

# Models by age
train_model_0_4 <- glm(flu ~ fever + 
                         cough + 
                         sore_throat +
                         myalgia + 
                         short_breath + 
                         hypoxemia + 
                         nausea_vom + 
                         bronchitis + 
                         chest_pain + 
                         diarrhea + 
                         fatigue + 
                         headache + 
                         congestion + 
                         sneezing,
                       family = binomial(link = "logit"),
                       data = flu_train[flu_train$age_group == '0-4',])

train_model_5_12 <- glm(flu ~ fever + 
                          cough + 
                          sore_throat +
                          myalgia + 
                          short_breath + 
                          hypoxemia + 
                          nausea_vom + 
                          bronchitis + 
                          chest_pain + 
                          diarrhea + 
                          fatigue + 
                          headache + 
                          congestion + 
                          sneezing,
                        family = binomial(link = "logit"),
                        data = flu_train[flu_train$age_group == '5-12',])

train_model_13_17 <- glm(flu ~ fever + 
                           cough + 
                           sore_throat +
                           myalgia + 
                           short_breath + 
                           hypoxemia + 
                           nausea_vom + 
                           bronchitis + 
                           chest_pain + 
                           diarrhea + 
                           fatigue + 
                           headache + 
                           congestion + 
                           sneezing,
                         family = binomial(link = "logit"),
                         data = flu_train[flu_train$age_group == '13-17',])

train_model_18_49 <- glm(flu ~ fever + 
                           cough + 
                           sore_throat +
                           myalgia + 
                           short_breath + 
                           hypoxemia + 
                           nausea_vom + 
                           bronchitis + 
                           chest_pain + 
                           diarrhea + 
                           fatigue + 
                           headache + 
                           congestion + 
                           sneezing,
                         family = binomial(link = "logit"),
                         data = flu_train[flu_train$age_group == '18-49',])

train_model_50_64 <- glm(flu ~ fever + 
                           cough + 
                           sore_throat +
                           myalgia + 
                           short_breath + 
                           hypoxemia + 
                           nausea_vom + 
                           bronchitis + 
                           chest_pain + 
                           diarrhea + 
                           fatigue + 
                           headache + 
                           congestion + 
                           sneezing,
                         family = binomial(link = "logit"),
                         data = flu_train[flu_train$age_group == '50-64',])

train_model_65 <- glm(flu ~ fever + 
                        cough + 
                        sore_throat +
                        myalgia + 
                        short_breath + 
                        hypoxemia + 
                        nausea_vom + 
                        bronchitis + 
                        chest_pain + 
                        diarrhea + 
                        fatigue + 
                        headache + 
                        congestion + 
                        sneezing,
                      family = binomial(link = "logit"),
                      data = flu_train[flu_train$age_group == '65+',])

# Models by geography
train_model_south_urban <- glm(flu ~ fever + 
                                 cough + 
                                 sore_throat +
                                 myalgia + 
                                 short_breath + 
                                 hypoxemia + 
                                 nausea_vom + 
                                 bronchitis + 
                                 chest_pain + 
                                 diarrhea + 
                                 fatigue + 
                                 headache + 
                                 congestion + 
                                 sneezing +
                                 as.factor(month),
                       family = binomial(link = "logit"),
                       data = flu_train[flu_train$south == 1 & flu_train$urban_code_binary == 1,])

train_model_south_rural <- glm(flu ~ fever + 
                                 cough + 
                                 sore_throat +
                                 myalgia + 
                                 short_breath + 
                                 hypoxemia + 
                                 nausea_vom + 
                                 bronchitis + 
                                 chest_pain + 
                                 diarrhea + 
                                 fatigue + 
                                 headache + 
                                 congestion + 
                                 sneezing +
                                 as.factor(month),
                               family = binomial(link = "logit"),
                               data = flu_train[flu_train$south == 1 & flu_train$urban_code_binary == 0,])

train_model_north_urban <- glm(flu ~ fever + 
                                 cough + 
                                 sore_throat +
                                 myalgia + 
                                 short_breath + 
                                 hypoxemia + 
                                 nausea_vom + 
                                 bronchitis + 
                                 chest_pain + 
                                 diarrhea + 
                                 fatigue + 
                                 headache + 
                                 congestion + 
                                 sneezing +
                                 as.factor(month),
                               family = binomial(link = "logit"),
                               data = flu_train[flu_train$south == 0 & flu_train$urban_code_binary == 1,])

train_model_north_rural <- glm(flu ~ fever + 
                                 cough + 
                                 sore_throat +
                                 myalgia + 
                                 short_breath + 
                                 hypoxemia + 
                                 nausea_vom + 
                                 bronchitis + 
                                 chest_pain + 
                                 diarrhea + 
                                 fatigue + 
                                 headache + 
                                 congestion + 
                                 sneezing +
                                 as.factor(month),
                               family = binomial(link = "logit"),
                               data = flu_train[flu_train$south == 0 & flu_train$urban_code_binary == 0,])

# Models by season
train_model_16 <- glm(flu ~ fever + 
                        cough + 
                        sore_throat +
                        myalgia + 
                        short_breath + 
                        hypoxemia + 
                        nausea_vom + 
                        bronchitis + 
                        chest_pain + 
                        diarrhea + 
                        fatigue + 
                        headache + 
                        congestion + 
                        sneezing +
                        as.factor(month),
                      family = binomial(link = "logit"),
                      data = flu_train[flu_train$Season == '2016-2017',])

train_model_17 <- glm(flu ~ fever + 
                        cough + 
                        sore_throat +
                        myalgia + 
                        short_breath + 
                        hypoxemia + 
                        nausea_vom + 
                        bronchitis + 
                        chest_pain + 
                        diarrhea + 
                        fatigue + 
                        headache + 
                        congestion + 
                        sneezing +
                        as.factor(month),
                      family = binomial(link = "logit"),
                      data = flu_train[flu_train$Season == '2017-2018',])

train_model_18 <- glm(flu ~ fever + 
                        cough + 
                        sore_throat +
                        myalgia + 
                        short_breath + 
                        hypoxemia + 
                        nausea_vom + 
                        bronchitis + 
                        chest_pain + 
                        diarrhea + 
                        fatigue + 
                        headache + 
                        congestion + 
                        sneezing +
                        as.factor(month),
                      family = binomial(link = "logit"),
                      data = flu_train[flu_train$Season == '2018-2019',])

train_model_19 <- glm(flu ~ fever + 
                        cough + 
                        sore_throat +
                        myalgia + 
                        short_breath + 
                        hypoxemia + 
                        nausea_vom + 
                        bronchitis + 
                        chest_pain + 
                        diarrhea + 
                        fatigue + 
                        headache + 
                        congestion + 
                        sneezing +
                        as.factor(month),
                      family = binomial(link = "logit"),
                      data = flu_train[flu_train$Season == '2019-2020',])

## SET THE OVERALL MODEL ORDER ##

# Create data for the overall model
coefficients <- overall_model$coefficients
se <- summary(overall_model)$coefficients[, 2]
names <- names(overall_model$coefficients)
model_all <- data.frame(names, coefficients, se)
model_all$num <- rank(model_all$coefficients)
model_all$lower <- exp(model_all$coefficients - 1.96*model_all$se)
model_all$upper <- exp(model_all$coefficients + 1.96*model_all$se)
model_all$mean <- exp(model_all$coefficients)

# Edit symptom names
model_all$names <- str_to_title(model_all$names)
model_all$names <- ifelse(model_all$names == 'Sore_throat', 
                          'Sore Throat', model_all$names)
model_all$names <- ifelse(model_all$names == 'Nausea_vom', 
                          'Nausea', model_all$names)
model_all$names <- ifelse(model_all$names == 'Chest_pain', 
                          'Chest Pain', model_all$names)
model_all$names <- ifelse(model_all$names == 'Short_breath', 
                          'Shortness of Breath', model_all$names)

model_all$names <- gsub("As.factor\\(Month\\)", "", model_all$names)

model_all$num <- ifelse(model_all$names == '2', -2, model_all$num)
model_all$num <- ifelse(model_all$names == '3', -3, model_all$num)
model_all$num <- ifelse(model_all$names == '4', -4, model_all$num)
model_all$num <- ifelse(model_all$names == '5', -5, model_all$num)
model_all$num <- ifelse(model_all$names == '6', -6, model_all$num)
model_all$num <- ifelse(model_all$names == '7', -7, model_all$num)
model_all$num <- ifelse(model_all$names == '8', -8, model_all$num)
model_all$num <- ifelse(model_all$names == '9', -9, model_all$num)
model_all$num <- ifelse(model_all$names == '10', -10, model_all$num)
model_all$num <- ifelse(model_all$names == '11', -11, model_all$num)
model_all$num <- ifelse(model_all$names == '12', -12, model_all$num)

# Remove the intercept
#model_all_flu_overall <- model_all[1:15,]
model_all_flu_overall <- model_all |> filter(names != '(Intercept)')


# Save the order
symptom_order_flu <- model_all_flu_overall[, c(1, 4)]

###################
# CREATE AGE PLOT #
###################

coefficients <- train_model_0_4$coefficients
se <- summary(train_model_0_4)$coefficients[, 2]
names <- names(train_model_0_4$coefficients)
model_0_4 <- data.frame(names, coefficients, se)
model_0_4$age_group <- '0-4'
#model_0_4 <- model_0_4[1:15,]

coefficients <- train_model_5_12$coefficients
se <- summary(train_model_5_12)$coefficients[, 2]
names <- names(train_model_5_12$coefficients)
model_5_12 <- data.frame(names, coefficients, se)
model_5_12$age_group <- '5-12'
#model_5_12 <- model_5_12[1:15,]

coefficients <- train_model_13_17$coefficients
se <- summary(train_model_13_17)$coefficients[, 2]
names <- names(train_model_13_17$coefficients)
model_13_17 <- data.frame(names, coefficients, se)
model_13_17$age_group <- '13-17'
#model_13_17 <- model_13_17[1:15,]

coefficients <- train_model_18_49$coefficients
se <- summary(train_model_18_49)$coefficients[, 2]
names <- names(train_model_18_49$coefficients)
model_18_49 <- data.frame(names, coefficients, se)
model_18_49$age_group <- '18-49'
#model_18_49 <- model_18_49[1:15,]

coefficients <- train_model_50_64$coefficients
se <- summary(train_model_50_64)$coefficients[, 2]
names <- names(train_model_50_64$coefficients)
model_50_64 <- data.frame(names, coefficients, se)
model_50_64$age_group <- '50-64'
#model_50_64 <- model_50_64[1:15,]

coefficients <- train_model_65$coefficients
se <- summary(train_model_65)$coefficients[, 2]
names <- names(train_model_65$coefficients)
model_65 <- data.frame(names, coefficients, se)
model_65$age_group <- '65+'
#model_65 <- model_65[1:15,]

# Combine data
models <- rbind(model_0_4, model_5_12, model_13_17, model_18_49, 
                model_50_64, model_65)

# Calculate lower and upper bounds
models$lower <- exp(models$coefficients - 1.96*models$se)
models$upper <- exp(models$coefficients + 1.96*models$se)
models$mean <- exp(models$coefficients)

# Clean up coefficients
models <- models |> filter(names != '(Intercept)')
models$names <- str_to_title(models$names)
models$names <- ifelse(models$names == 'Sore_throat', 
                       'Sore Throat', models$names)
models$names <- ifelse(models$names == 'Nausea_vom', 
                       'Nausea', models$names)
models$names <- ifelse(models$names == 'Chest_pain', 
                       'Chest Pain', models$names)
models$names <- ifelse(models$names == 'Short_breath', 
                       'Shortness of Breath', models$names)

# Create figure
colors <- c('0-4' = '#f3d1c0', '5-12' = '#eaad8f', '13-17' = '#e0885d', 
            '18-49' = '#d7642c', '50-64' = '#b95423', '65+' = '#873e1a')

shapes <- c(' ' = 5)

models$`Age Group Estimate` <- models$age_group
models$`Age Group Estimate` <- factor(models$age_group, levels=c('0-4', '5-12', '13-17', '18-49', '50-64', '65+'))
models <- left_join(models, symptom_order_flu, by = c('names'))

model_all_flu_overall$`Overall Estimate` <- ' '

forrest_theme <- theme(legend.position = "inside",
                       legend.position.inside = c(0.70, 0.20),
                       axis.text = element_text(size=14),
                       axis.title = element_text(size=16),
                       axis.title.y = element_blank(),
                       strip.background = element_blank(),
                       legend.title = element_text(size=16),
                       legend.text = element_text(size=14),
                       plot.title = element_text(size=18, hjust = 0),
                       legend.box="vertical",
                       legend.box.background = element_rect(colour = "black", 
                                                            fill = 'white'))


age_flu_fig <- ggplot(data=models, aes(x=reorder(names, num), y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(aes(group = `Age Group Estimate`, color = `Age Group Estimate`, fill = `Age Group Estimate`),
                  position = position_dodge2(width = 0.4), size = 0.75, alpha = 0.7) + 
  geom_point(data=model_all_flu_overall[1:14,], aes(x=reorder(names, num), y=mean, shape = `Overall Estimate`),
             color = 'black', alpha = 2, size = 3, stroke = 1.2) +
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Odds Ratio (95% CI)") +
  ggtitle("Influenza Regression Results by Age") +
  theme_minimal() +  # use a white background
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(order = 1)) +
  forrest_theme

age_flu_fig

###################
# CREATE AGE PLOT #
###################

coefficients <- age_model$coefficients
se <- summary(age_model)$coefficients[, 2]
names <- names(age_model$coefficients)
model_age <- data.frame(names, coefficients, se)
model_age <- model_age[3:7,]
model_age$names <- substr(model_age$names, 21, nchar(model_age$names))
model_age$num <- rank(model_age$coefficients)

model_age$lower <- exp(model_age$coefficients - 1.96*model_age$se)
model_age$upper <- exp(model_age$coefficients + 1.96*model_age$se)
model_age$mean <- exp(model_age$coefficients)

forrest_other_theme <- theme(legend.position = "none",
                       legend.position.inside = c(0.70, 0.20),
                       axis.text = element_text(size=14),
                       axis.title = element_text(size=16),
                       strip.background = element_blank(),
                       legend.title = element_text(size=16),
                       legend.text = element_text(size=14),
                       plot.title = element_text(size=18, hjust = 0),
                       legend.box="vertical",
                       legend.box.background = element_rect(colour = "black", 
                                                            fill = 'white'))

age_overall_flu_fig <- ggplot(data=model_age, aes(x=reorder(names, -num), y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(aes(color = names), size = 0.75, alpha = 0.7) + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  xlab("Age Group") + ylab("Odds Ratio (95% CI)") + ylim(0, 1.25) +
  ggtitle("Influenza Odds Relative to Age 0-4") +
  scale_color_manual(values = colors) +
  theme_minimal() +  # use a white background
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(order = 1)) +
  forrest_other_theme 

age_overall_flu_fig



# MONTH

coefficients <- train_model_north_rural$coefficients
se <- summary(train_model_north_rural)$coefficients[, 2]
names <- names(train_model_north_rural$coefficients)
model_north_rural <- data.frame(names, coefficients, se)
model_north_rural$`Geography` <- 'North (Rural)'
model_north_rural <- model_north_rural[16:26,]

coefficients <- train_model_south_rural$coefficients
se <- summary(train_model_south_rural)$coefficients[, 2]
names <- names(train_model_south_rural$coefficients)
model_south_rural <- data.frame(names, coefficients, se)
model_south_rural$`Geography` <- 'South (Rural)'
model_south_rural <- model_south_rural[16:26,]

coefficients <- train_model_north_urban$coefficients
se <- summary(train_model_north_urban)$coefficients[, 2]
names <- names(train_model_north_urban$coefficients)
model_north_urban <- data.frame(names, coefficients, se)
model_north_urban$`Geography` <- 'North (Urban)'
model_north_urban <- model_north_urban[16:26,]

coefficients <- train_model_south_urban$coefficients
se <- summary(train_model_south_urban)$coefficients[, 2]
names <- names(train_model_south_urban$coefficients)
model_south_urban <- data.frame(names, coefficients, se)
model_south_urban$`Geography` <- 'South (Urban)'
model_south_urban <- model_south_urban[16:26,]

# Combine data
models <- rbind(model_north_urban, model_south_urban, model_south_rural, model_north_rural)

# Calculate lower and upper bounds
models$lower <- exp(models$coefficients - 1.96*models$se)
models$upper <- exp(models$coefficients + 1.96*models$se)
models$mean <- exp(models$coefficients)



models <- models |> filter(names != '(Intercept)')
models$names <- str_to_title(models$names)
models$names <- ifelse(models$names == 'Sore_throat', 
                       'Sore Throat', models$names)
models$names <- ifelse(models$names == 'Nausea_vom', 
                       'Nausea', models$names)
models$names <- ifelse(models$names == 'Chest_pain', 
                       'Chest Pain', models$names)
models$names <- ifelse(models$names == 'Short_breath', 
                       'Shortness of Breath', models$names)



models$names <- gsub("As.factor\\(Month\\)", "", models$names)


# Clean up coeffitients


# Create figure
coefficients <- time_model$coefficients
se <- summary(time_model)$coefficients[, 2]
names <- names(time_model$coefficients)
time_order <- data.frame(names, coefficients, se)
time_order$mean <- exp(time_order$coefficients)
time_order$lower <- exp(time_order$coefficients - 1.96*time_order$se)
time_order$upper <- exp(time_order$coefficients + 1.96*time_order$se)
time_order$names <- gsub("as.factor\\(month\\)", "", time_order$names)


shapes <- c(' ' = 5)

time_order$`Overall Estimate` <- ' '

time_order <- time_order[21:31,]

colors <- c('South (Urban)' = '#509048', 'North (Urban)' = '#896deb', 
            'South (Rural)' = '#90c058', 'North (Rural)' = '#c2a3fd')

months <- c("Feb","Mar",
              "Apr","May","Jun",
              "Jul","Aug","Sep",
              "Oct","Nov","Dec")

month_num <- as.character(c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))

month <- data.frame(months, month_num)

time_order <- left_join(time_order, month, by = c('names' = 'month_num'))
models <- left_join(models, month, by = c('names' = 'month_num'))


geo_flu_fig <- ggplot(data=models, aes(x=reorder(months, as.numeric(names)), y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(data = models, aes(group = `Geography`, color = `Geography`, fill = `Geography`),
                  position = position_dodge2(width = 0.4), size = 0.75, alpha = 0.7) + 
  geom_point(data=time_order, aes(x=reorder(months, as.numeric(names)), y=mean, shape = `Overall Estimate`),
             color = 'black', alpha = 2, size = 3, stroke = 1.2) +
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  xlab("Month") + ylab("Odds Ratio (95% CI)") +
  ggtitle("Influenza Odds Relative to January by Geography") +
  theme_minimal() +  # use a white background
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = colors) +
  guides(shape = guide_legend(order = 1)) +
  forrest_other_theme +
  theme(plot.title = element_text(size=18, hjust = 0)) +
  theme(legend.position = 'bottom')

geo_flu_fig




## SEASON ##



coefficients <- train_model_16$coefficients
se <- summary(train_model_16)$coefficients[, 2]
names <- names(train_model_16$coefficients)
model_north_rural <- data.frame(names, coefficients, se)
model_north_rural$`Season` <- '16-17'
#model_north_rural <- model_north_rural[21:31,]

coefficients <- train_model_17$coefficients
se <- summary(train_model_17)$coefficients[, 2]
names <- names(train_model_17$coefficients)
model_south_rural <- data.frame(names, coefficients, se)
model_south_rural$`Season` <- '17-18'
#model_south_rural <- model_south_rural[21:31,]

coefficients <- train_model_18$coefficients
se <- summary(train_model_18)$coefficients[, 2]
names <- names(train_model_18$coefficients)
model_north_urban <- data.frame(names, coefficients, se)
model_north_urban$`Season` <- '18-19'
#model_north_urban <- model_north_urban[21:31,]

coefficients <- train_model_19$coefficients
se <- summary(train_model_19)$coefficients[, 2]
names <- names(train_model_19$coefficients)
model_south_urban <- data.frame(names, coefficients, se)
model_south_urban$`Season` <- '19-20'
#model_south_urban <- model_south_urban[21:31,]

# Combine data
models <- rbind(model_north_urban, model_south_urban, model_south_rural, model_north_rural)

# Calculate lower and upper bounds
models$lower <- exp(models$coefficients - 1.96*models$se)
models$upper <- exp(models$coefficients + 1.96*models$se)
models$mean <- exp(models$coefficients)

models <- models |> filter(names != '(Intercept)')
models$names <- str_to_title(models$names)
models$names <- ifelse(models$names == 'Sore_throat', 
                       'Sore Throat', models$names)
models$names <- ifelse(models$names == 'Nausea_vom', 
                       'Nausea', models$names)
models$names <- ifelse(models$names == 'Chest_pain', 
                       'Chest Pain', models$names)
models$names <- ifelse(models$names == 'Short_breath', 
                       'Shortness of Breath', models$names)



models$names <- gsub("As.factor\\(Month\\)", "", models$names)


# Clean up coeffitients


# Create figure


shapes <- c(' ' = 5)

model_all_flu_overall$`Overall Estimate` <- ' '

models <- left_join(models, symptom_order_flu, by = c('names'))


season_flu_fig <- ggplot(data=models, aes(x=reorder(names, num), y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(aes(group = `Season`, color = `Season`, fill = `Season`),
                  position = position_dodge2(width = 0.4), size = 0.75, alpha = 0.7) + 
  geom_point(data=model_all_flu_overall, aes(x=reorder(names, num), y=mean, shape = `Overall Estimate`),
             color = 'black', alpha = 2, size = 3, stroke = 1.2) +
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Odds Ratio (95% CI)") +
  ggtitle("Influenza Regression Results by Season") +
  theme_minimal() +  # use a white background
  scale_shape_manual(values = shapes) +
  guides(shape = guide_legend(order = 1)) +
  forrest_theme

season_flu_fig



season_time_flu_fig <- ggplot(data=models[nchar(models$names) < 3,], aes(x=reorder(names, -num), y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(aes(group = `Season`, color = `Season`, fill = `Season`),
                  position = position_dodge2(width = 0.4), size = 0.75, alpha = 0.7) + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  xlab("Month") + ylab("Odds Ratio (95% CI)") +
  ggtitle("Odds of Influenza Relative to January") +
  theme_minimal() +  # use a white background
  forrest_theme + theme(plot.title = element_text(size=18, hjust = 0)) +
  theme(legend.position = 'bottom')

season_time_flu_fig







###########################
# CREATE MAIN EFFECT PLOT #
###########################

# Models by time
train_model_all <- glm(flu ~ fever + cough + sore_throat +
                         myalgia + short_breath + hypoxemia +
                         nausea_vom + bronchitis + chest_pain + 
                         diarrhea + fatigue + headache + 
                         congestion + sneezing,
                       family = binomial(link = "logit"),
                       data = flu_train)

train_model_age <- glm(flu ~ fever + cough + sore_throat +
                         myalgia + short_breath + hypoxemia +
                         nausea_vom + bronchitis + chest_pain + 
                         diarrhea + fatigue + headache + 
                         congestion + sneezing +
                         age_group,
                       family = binomial(link = "logit"),
                       data = flu_train)

train_model_time <- glm(flu ~ fever + cough + sore_throat +
                          myalgia + short_breath + hypoxemia +
                          nausea_vom + bronchitis + chest_pain + 
                          diarrhea + fatigue + headache + 
                          congestion + sneezing +
                          age_group +
                          as.factor(month),
                       family = binomial(link = "logit"),
                       data = flu_train)

coefficients <- train_model_all$coefficients
se <- summary(train_model_all)$coefficients[, 2]
names <- names(train_model_all$coefficients)
model_all <- data.frame(names, coefficients, se)
model_all$`Model` <- 'All Symptoms'

coefficients <- train_model_age$coefficients
se <- summary(train_model_age)$coefficients[, 2]
names <- names(train_model_age$coefficients)
model_age <- data.frame(names, coefficients, se)
model_age$`Model` <- 'All Symptoms + Age'
model_age <- model_age[1:15,]

coefficients <- train_model_time$coefficients
se <- summary(train_model_time)$coefficients[, 2]
names <- names(train_model_time$coefficients)
model_time <- data.frame(names, coefficients, se)
model_time$`Model` <- 'All Symptoms + Age + Time'
model_time <- model_time[1:15,]

# Combine data
models <- rbind(model_all, model_age, model_time)

# Calculate lower and upper bounds
models$lower <- exp(models$coefficients - 1.96*models$se)
models$upper <- exp(models$coefficients + 1.96*models$se)
models$mean <- exp(models$coefficients)

# Clean up coeffitients
models <- models |> filter(names != '(Intercept)')
models$names <- str_to_title(models$names)
models$names <- ifelse(models$names == 'Sore_throat', 
                       'Sore Throat', models$names)
models$names <- ifelse(models$names == 'Nausea_vom', 
                       'Nausea', models$names)
models$names <- ifelse(models$names == 'Chest_pain', 
                       'Chest Pain', models$names)
models$names <- ifelse(models$names == 'Short_breath', 
                       'Shortness of Breath', models$names)

# Create figure
colors <- c('All Symptoms' = '#d7642c', 'All Symptoms + Age' = '#90c058', 
            'All Symptoms + Age + Time' = '#509048'
)


models <- left_join(models, symptom_order_flu, by = c('names'))

main_flu_fig <- ggplot(data=models, aes(x=reorder(names, num), y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(aes(group = `Model`, color = `Model`, fill = `Model`),
                  position = position_dodge2(width = 0.4), size = 0.75, alpha = 0.7) + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  scale_y_continuous(limits=c(0, 8)) +
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Odds Ratio (95% CI)") +
  ggtitle("Influenza Regression Results by Model") +
  theme_minimal() +  # use a white background
  scale_color_manual(values = colors) +
  forrest_theme + theme(plot.title = element_text(size=18, hjust = 0))

main_flu_fig
 
#######
# RSV #
#######

# Add time variables to Part 1 data
rsv_dat_p1$month <- month(rsv_dat_p1$month_date)
rsv_dat_p1$year <- year(rsv_dat_p1$month_date)

# Create age group variable
# Part 1
rsv_dat_p1$age_group <- factor(as.numeric(rsv_dat_p1$age_grp),
                               levels = c(0, 1, 2, 3, 4),
                               labels = c('0-2', '3-4', '5-17', '18-64', '65+'))

# Load population data
county_pop <- readRDS('tmp/county_group_pop.rds')

# Merge on south and urban info
rsv_dat_p1 <- left_join(rsv_dat_p1, urban_south_binary, by = 'county_fips')

# Restrict to counties with 10,000 or more individuals
rsv_dat_p1 <- left_join(rsv_dat_p1, county_pop, 
                        by = c('county_fips' = 'county_fips'))
rsv_dat_p1 <- rsv_dat_p1 |> filter(county_group_pop > 10000)

# Make data individual level (sample 1 million individuals)
# Balance flu casea and non-flu cases
balanced <- rsv_dat_p1 |> group_by(rsv) |>
  uncount(patient_count_imp) |> sample_n(500000)

# Split the data randomly into testing and training
sample <- sample(c(TRUE, FALSE), nrow(balanced), replace=TRUE, prob=c(0.5, 0.5))
rsv_train  <- balanced[sample, ]
rsv_test   <- balanced[!sample, ]

## Define candidate models ##


# Fever only model
hyp_model <- glm(rsv ~ hypoxemia,
                 family = binomial(link = "logit"),
                 data = rsv_train)
summary(hyp_model)

# Save
saveRDS(hyp_model, "./tmp/hypox_model_rsv.rds")

# All symptoms model
symp_model <- glm(rsv ~ fever + 
                    cough + 
                    sore_throat +
                    short_breath + 
                    hypoxemia + 
                    bronchitis + 
                    loss_appetite + 
                    fatigue + 
                    headache + 
                    congestion + 
                    sneezing,
                  family = binomial(link = "logit"),
                  data = rsv_train)
summary(symp_model)

# Save
saveRDS(symp_model, "./tmp/symptom_model_rsv.rds")

# Age model
age_model <- glm(rsv ~ fever*as.factor(age_group) + 
                   cough*as.factor(age_group) + 
                   sore_throat*as.factor(age_group) +
                   short_breath*as.factor(age_group) + 
                   hypoxemia*as.factor(age_group) + 
                   bronchitis*as.factor(age_group) + 
                   loss_appetite*as.factor(age_group) + 
                   fatigue*as.factor(age_group) + 
                   headache*as.factor(age_group) + 
                   congestion*as.factor(age_group) + 
                   sneezing*as.factor(age_group),
                 family = binomial(link = "logit"),
                 data = rsv_train)
summary(age_model)

# Save
saveRDS(age_model, "./tmp/age_model_rsv.rds")

# Time model
time_model <- glm(rsv ~ fever*as.factor(age_group) + 
                    cough*as.factor(age_group) + 
                    sore_throat*as.factor(age_group) +
                    short_breath*as.factor(age_group) + 
                    hypoxemia*as.factor(age_group) + 
                    bronchitis*as.factor(age_group) + 
                    loss_appetite*as.factor(age_group) + 
                    fatigue*as.factor(age_group) + 
                    headache*as.factor(age_group) + 
                    congestion*as.factor(age_group) + 
                    sneezing*as.factor(age_group) +
                    as.factor(month),
                  family = binomial(link = "logit"),
                  data = rsv_train)
summary(time_model)

# Save
saveRDS(time_model, "./tmp/time_model_rsv.rds")

# Geography model
geo_model <- glm(rsv ~ fever*as.factor(age_group) + 
                   cough*as.factor(age_group) + 
                   sore_throat*as.factor(age_group) +
                   short_breath*as.factor(age_group) + 
                   hypoxemia*as.factor(age_group) + 
                   bronchitis*as.factor(age_group) + 
                   loss_appetite*as.factor(age_group) + 
                   fatigue*as.factor(age_group) + 
                   headache*as.factor(age_group) + 
                   congestion*as.factor(age_group) + 
                   sneezing*as.factor(age_group) +
                   as.factor(month)*urban_code_binary +
                   as.factor(month)*south,
                 family = binomial(link = "logit"),
                 data = rsv_train)
summary(geo_model)

# Save
saveRDS(geo_model, "./tmp/geography_model_rsv.rds")



# Quick model tests
rsv_test$pred <- predict(hyp_model, newdata = rsv_test, type = "response")
pROC::auc(pROC::roc(rsv_test$rsv, rsv_test$pred))
ggplot(rsv_test) + geom_density(aes(x = pred, group = as.factor(rsv), color = as.factor(rsv)))

rsv_test$pred <- predict(symp_model, newdata = rsv_test, type = "response")
pROC::auc(pROC::roc(rsv_test$rsv, rsv_test$pred))
ggplot(rsv_test) + geom_density(aes(x = pred, group = as.factor(rsv), color = as.factor(rsv)))

rsv_test$pred <- predict(age_model, newdata = rsv_test, type = "response")
pROC::auc(pROC::roc(rsv_test$rsv, rsv_test$pred))
ggplot(rsv_test) + geom_density(aes(x = pred, group = as.factor(rsv), color = as.factor(rsv)))

rsv_test$pred <- predict(time_model, newdata = rsv_test, type = "response")
pROC::auc(pROC::roc(rsv_test$rsv, rsv_test$pred))
ggplot(rsv_test) + geom_density(aes(x = pred, group = as.factor(rsv), color = as.factor(rsv)))

rsv_test$pred <- predict(geo_model, newdata = rsv_test, type = "response")
pROC::auc(pROC::roc(rsv_test$rsv, rsv_test$pred))
ggplot(rsv_test) + geom_density(aes(x = pred, group = as.factor(rsv), color = as.factor(rsv)))





# Models by age
train_model_0_2 <- glm(rsv ~ fever + cough + sore_throat +
                         short_breath + hypoxemia + 
                         bronchitis + loss_appetite + 
                         fatigue + headache + congestion + sneezing,
                       family = binomial(link = "logit"),
                       data = rsv_train[rsv_train$age_grp == '0',])

train_model_3_4 <- glm(rsv ~ fever + cough + sore_throat +
                         short_breath + hypoxemia + 
                         bronchitis + loss_appetite + 
                         fatigue + headache + congestion + sneezing,
                       family = binomial(link = "logit"),
                       data = rsv_train[rsv_train$age_grp == '1',])

train_model_5_17 <- glm(rsv ~ fever + cough + sore_throat +
                          short_breath + hypoxemia + 
                          bronchitis + loss_appetite + 
                          fatigue + headache + congestion + sneezing,
                        family = binomial(link = "logit"),
                        data = rsv_train[rsv_train$age_grp == '2',])

train_model_18_64 <- glm(rsv ~ fever + cough + sore_throat +
                           short_breath + hypoxemia + 
                           bronchitis + loss_appetite + 
                           fatigue + headache + congestion + sneezing,
                         family = binomial(link = "logit"),
                         data = rsv_train[rsv_train$age_grp == '3',])

train_model_65 <- glm(rsv ~ fever + cough + sore_throat +
                        short_breath + hypoxemia + 
                        bronchitis + loss_appetite + 
                        fatigue + headache + congestion + sneezing,
                      family = binomial(link = "logit"),
                      data = rsv_train[rsv_train$age_grp == '4',])

# Models by time
train_model_jan <- glm(rsv ~ fever + cough + sore_throat +
                         short_breath + hypoxemia + 
                         bronchitis + loss_appetite + 
                         fatigue + headache + congestion + sneezing,
                       family = binomial(link = "logit"),
                       data = rsv_train[rsv_train$month == 1,])

train_model_apr <- glm(rsv ~ fever + cough + sore_throat +
                         short_breath + hypoxemia + 
                         bronchitis + loss_appetite + 
                         fatigue + headache + congestion + sneezing,
                       family = binomial(link = "logit"),
                       data = rsv_train[rsv_train$month == 4,])

train_model_jul <- glm(rsv ~ fever + cough + sore_throat +
                         short_breath + hypoxemia + 
                         bronchitis + loss_appetite + 
                         fatigue + headache + congestion + sneezing,
                       family = binomial(link = "logit"),
                       data = rsv_train[rsv_train$month == 7,])

train_model_oct <- glm(rsv ~ fever + cough + sore_throat +
                         short_breath + hypoxemia + 
                         bronchitis + loss_appetite + 
                         fatigue + headache + congestion + sneezing,
                       family = binomial(link = "logit"),
                       data = rsv_train[rsv_train$month == 10,])

train_model_south_urban <- glm(rsv ~ fever + cough + sore_throat +
                                 short_breath + hypoxemia + 
                                 bronchitis + loss_appetite + 
                                 fatigue + headache + congestion + sneezing +
                                 as.factor(month),
                               family = binomial(link = "logit"),
                               data = rsv_train[rsv_train$south == 1 & rsv_train$urban_code_binary == 1,])

train_model_south_rural <- glm(rsv ~ fever + cough + sore_throat +
                                 short_breath + hypoxemia + 
                                 bronchitis + loss_appetite + 
                                 fatigue + headache + congestion + sneezing +
                                 as.factor(month),
                               family = binomial(link = "logit"),
                               data = rsv_train[rsv_train$south == 1 & rsv_train$urban_code_binary == 0,])

train_model_north_urban <- glm(rsv ~ fever + cough + sore_throat +
                                 short_breath + hypoxemia + 
                                 bronchitis + loss_appetite + 
                                 fatigue + headache + congestion + sneezing +
                                 as.factor(month),
                               family = binomial(link = "logit"),
                               data = rsv_train[rsv_train$south == 0 & rsv_train$urban_code_binary == 1,])

train_model_north_rural <- glm(rsv ~ fever + cough + sore_throat +
                                 short_breath + hypoxemia + 
                                 bronchitis + loss_appetite + 
                                 fatigue + headache + congestion + sneezing +
                                 as.factor(month),
                               family = binomial(link = "logit"),
                               data = rsv_train[rsv_train$south == 0 & rsv_train$urban_code_binary == 0,])

# Create data for the overall model
coefficients <-symp_model$coefficients
se <- summary(symp_model)$coefficients[, 2]
names <- names(symp_model$coefficients)
model_all <- data.frame(names, coefficients, se)
model_all$num <- rank(model_all$coefficients)
model_all$lower <- exp(model_all$coefficients - 1.96*model_all$se)
model_all$upper <- exp(model_all$coefficients + 1.96*model_all$se)
model_all$mean <- exp(model_all$coefficients)

# Edit symptom names
model_all$names <- str_to_title(model_all$names)
model_all$names <- ifelse(model_all$names == 'Sore_throat', 
                          'Sore Throat', model_all$names)
model_all$names <- ifelse(model_all$names == 'Nausea_vom', 
                          'Nausea', model_all$names)
model_all$names <- ifelse(model_all$names == 'Loss_appetite', 
                                'Loss of Appetite', model_all$names)
model_all$names <- ifelse(model_all$names == 'Chest_pain', 
                          'Chest Pain', model_all$names)
model_all$names <- ifelse(model_all$names == 'Short_breath', 
                          'Shortness of Breath', model_all$names)

# Remove the intercept
model_all_rsv_overall <- model_all |> filter(names != '(Intercept)')

# Save the order
symptom_order_rsv <- model_all_rsv_overall[, c(1, 4)]

###################
# CREATE AGE PLOT #
###################

coefficients <- train_model_0_2$coefficients
se <- summary(train_model_0_2)$coefficients[, 2]
names <- names(train_model_0_2$coefficients)
model_0_2 <- data.frame(names, coefficients, se)
model_0_2$age_group <- '0-2'

coefficients <- train_model_3_4$coefficients
se <- summary(train_model_3_4)$coefficients[, 2]
names <- names(train_model_3_4$coefficients)
model_3_4 <- data.frame(names, coefficients, se)
model_3_4$age_group <- '3-4'

coefficients <- train_model_5_17$coefficients
se <- summary(train_model_5_17)$coefficients[, 2]
names <- names(train_model_5_17$coefficients)
model_5_17 <- data.frame(names, coefficients, se)
model_5_17$age_group <- '5-17'

coefficients <- train_model_18_64$coefficients
se <- summary(train_model_18_64)$coefficients[, 2]
names <- names(train_model_18_64$coefficients)
model_18_64 <- data.frame(names, coefficients, se)
model_18_64$age_group <- '18-64'

coefficients <- train_model_65$coefficients
se <- summary(train_model_65)$coefficients[, 2]
names <- names(train_model_65$coefficients)
model_65 <- data.frame(names, coefficients, se)
model_65$age_group <- '65+  '

# Combine data
models <- rbind(model_0_2, model_3_4, model_5_17, model_18_64, 
                model_65)

# Calculate lower and upper bounds
models$lower <- exp(models$coefficients - 1.96*models$se)
models$upper <- exp(models$coefficients + 1.96*models$se)
models$mean <- exp(models$coefficients)

# Clean up coefficients
models <- models |> filter(names != '(Intercept)')
models$names <- str_to_title(models$names)
models$names <- ifelse(models$names == 'Sore_throat', 
                       'Sore Throat', models$names)
models$names <- ifelse(models$names == 'Nausea_vom', 
                       'Nausea', models$names)
models$names <- ifelse(models$names == 'Chest_pain', 
                       'Chest Pain', models$names)
models$names <- ifelse(models$names == 'Short_breath', 
                       'Shortness of Breath', models$names)
models$names <- ifelse(models$names == 'Loss_appetite', 
                          'Loss of Appetite', models$names)

# Create figure
colors <- c('0-2' = '#9edbd8', '3-4' = '#73cbc7', '5-17' = '#41afaa', 
            '18-64' = '#36928e', '65+  ' = '#2c7672'
)

shapes <- c(' ' = 5)

model_all_rsv_overall$`Overall Estimate` <- ' '

models$`Age Group Estimate` <- models$age_group
models$`Age Group Estimate` <- factor(models$age_group, levels=c('0-2', '3-4', '5-17', '18-64', '65+  '))
models <- left_join(models, symptom_order_rsv, by = c('names'))

age_rsv_fig <- ggplot(data=models, aes(x=reorder(names, num), y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(aes(group = `Age Group Estimate`, color = `Age Group Estimate`, fill = `Age Group Estimate`),
                  position = position_dodge2(width = 0.4), size = 0.75, alpha = 0.7) + 
  geom_point(data=model_all_rsv_overall, aes(x=reorder(names, num), y=mean, shape = `Overall Estimate`),
             color = 'black', alpha = 2, size = 3, stroke = 1.2) +
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Coefficients") + ylab("Odds Ratio (95% CI)") +
  ggtitle("RSV Regression Results by Age") +
  theme_minimal() +  # use a white background
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(order = 1)) +
  forrest_theme + theme(plot.title = element_text(size=18, hjust = 0))

age_rsv_fig


## OTHER AGE PLOT ##

coefficients <- age_model$coefficients
se <- summary(age_model)$coefficients[, 2]
names <- names(age_model$coefficients)
model_age <- data.frame(names, coefficients, se)
model_age <- model_age[3:6,]
model_age$names <- substr(model_age$names, 21, nchar(model_age$names))
model_age$num <- rank(model_age$coefficients)

model_age$lower <- exp(model_age$coefficients - 1.96*model_age$se)
model_age$upper <- exp(model_age$coefficients + 1.96*model_age$se)
model_age$mean <- exp(model_age$coefficients)

model_age[model_age$names == '65+',]$names <- '65+  '
model_age$num <- c(1, 2, 3, 4)

forrest_other_theme <- theme(legend.position = "none",
                             legend.position.inside = c(0.70, 0.20),
                             axis.text = element_text(size=14),
                             axis.title = element_text(size=16),
                             strip.background = element_blank(),
                             legend.title = element_text(size=16),
                             legend.text = element_text(size=14),
                             plot.title = element_text(size=18, hjust = 0),
                             legend.box="vertical",
                             legend.box.background = element_rect(colour = "black", 
                                                                  fill = 'white'))

age_overall_rsv_fig <- ggplot(data=model_age, aes(x=reorder(names, num), y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(aes(color = names), size = 0.75, alpha = 0.7) + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  xlab("Age Group") + ylab("Odds Ratio (95% CI)") + ylim(0, 1) +
  ggtitle("RSV Odds Relative to Age 0-2") +
  scale_color_manual(values = colors) +
  theme_minimal() +  # use a white background
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(order = 1)) +
  forrest_other_theme 

age_overall_rsv_fig


## SEASONALITY PLOT ##

# MONTH

coefficients <- train_model_north_rural$coefficients
se <- summary(train_model_north_rural)$coefficients[, 2]
names <- names(train_model_north_rural$coefficients)
model_north_rural <- data.frame(names, coefficients, se)
model_north_rural$`Geography` <- 'North (Rural)'
model_north_rural <- model_north_rural[13:23,]

coefficients <- train_model_south_rural$coefficients
se <- summary(train_model_south_rural)$coefficients[, 2]
names <- names(train_model_south_rural$coefficients)
model_south_rural <- data.frame(names, coefficients, se)
model_south_rural$`Geography` <- 'South (Rural)'
model_south_rural <- model_south_rural[13:23,]

coefficients <- train_model_north_urban$coefficients
se <- summary(train_model_north_urban)$coefficients[, 2]
names <- names(train_model_north_urban$coefficients)
model_north_urban <- data.frame(names, coefficients, se)
model_north_urban$`Geography` <- 'North (Urban)'
model_north_urban <- model_north_urban[13:23,]

coefficients <- train_model_south_urban$coefficients
se <- summary(train_model_south_urban)$coefficients[, 2]
names <- names(train_model_south_urban$coefficients)
model_south_urban <- data.frame(names, coefficients, se)
model_south_urban$`Geography` <- 'South (Urban)'
model_south_urban <- model_south_urban[13:23,]

# Combine data
models <- rbind(model_north_urban, model_south_urban, model_south_rural, model_north_rural)

# Calculate lower and upper bounds
models$lower <- exp(models$coefficients - 1.96*models$se)
models$upper <- exp(models$coefficients + 1.96*models$se)
models$mean <- exp(models$coefficients)



models <- models |> filter(names != '(Intercept)')
models$names <- str_to_title(models$names)
models$names <- ifelse(models$names == 'Sore_throat', 
                       'Sore Throat', models$names)
models$names <- ifelse(models$names == 'Nausea_vom', 
                       'Nausea', models$names)
models$names <- ifelse(models$names == 'Chest_pain', 
                       'Chest Pain', models$names)
models$names <- ifelse(models$names == 'Short_breath', 
                       'Shortness of Breath', models$names)



models$names <- gsub("As.factor\\(Month\\)", "", models$names)


# Clean up coeffitients


# Create figure
coefficients <- time_model$coefficients
se <- summary(time_model)$coefficients[, 2]
names <- names(time_model$coefficients)
time_order <- data.frame(names, coefficients, se)
time_order$mean <- exp(time_order$coefficients)
time_order$lower <- exp(time_order$coefficients - 1.96*time_order$se)
time_order$upper <- exp(time_order$coefficients + 1.96*time_order$se)
time_order$names <- gsub("as.factor\\(month\\)", "", time_order$names)


shapes <- c(' ' = 5)

time_order$`Overall Estimate` <- ' '

time_order <- time_order[17:27,]

colors <- c('South (Urban)' = '#509048', 'North (Urban)' = '#896deb', 
            'South (Rural)' = '#90c058', 'North (Rural)' = '#c2a3fd')

months <- c("Feb","Mar",
            "Apr","May","Jun",
            "Jul","Aug","Sep",
            "Oct","Nov","Dec")

month_num <- as.character(c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))

month <- data.frame(months, month_num)

time_order <- left_join(time_order, month, by = c('names' = 'month_num'))
models <- left_join(models, month, by = c('names' = 'month_num'))


geo_rsv_fig <- ggplot(data=models, aes(x=reorder(months, as.numeric(names)), y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(data = models, aes(group = `Geography`, color = `Geography`, fill = `Geography`),
                  position = position_dodge2(width = 0.4), size = 0.75, alpha = 0.7) + 
  geom_point(data=time_order, aes(x=reorder(months, as.numeric(names)), y=mean, shape = `Overall Estimate`),
             color = 'black', alpha = 2, size = 3, stroke = 1.2) +
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  xlab("Month") + ylab("Odds Ratio (95% CI)") +
  ggtitle("RSV Odds Relative to January by Geography") +
  theme_minimal() +  # use a white background
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = colors) +
  guides(shape = guide_legend(order = 1)) +
  forrest_other_theme +
  theme(plot.title = element_text(size=18, hjust = 0)) +
  theme(legend.position = 'bottom')

geo_rsv_fig



coefficients <- train_model_all$coefficients
se <- summary(train_model_all)$coefficients[, 2]
names <- names(train_model_all$coefficients)
model_all <- data.frame(names, coefficients, se)
model_all$`Model` <- 'All Symptoms'

coefficients <- train_model_age$coefficients
se <- summary(train_model_age)$coefficients[, 2]
names <- names(train_model_age$coefficients)
model_age <- data.frame(names, coefficients, se)
model_age$`Model` <- 'All Symptoms + Age'
model_age <- model_age[1:12,]

coefficients <- train_model_time$coefficients
se <- summary(train_model_time)$coefficients[, 2]
names <- names(train_model_time$coefficients)
model_time <- data.frame(names, coefficients, se)
model_time$`Model` <- 'All Symptoms + Age + Time'
model_time <- model_time[1:12,]

# Combine data
models <- rbind(model_all, model_age, model_time)

# Calculate lower and upper bounds
models$lower <- exp(models$coefficients - 1.96*models$se)
models$upper <- exp(models$coefficients + 1.96*models$se)
models$mean <- exp(models$coefficients)

# Clean up coeffitients
models <- models |> filter(names != '(Intercept)')
models$names <- str_to_title(models$names)
models$names <- ifelse(models$names == 'Sore_throat', 
                       'Sore Throat', models$names)
models$names <- ifelse(models$names == 'Nausea_vom', 
                       'Nausea', models$names)
models$names <- ifelse(models$names == 'Chest_pain', 
                       'Chest Pain', models$names)
models$names <- ifelse(models$names == 'Short_breath', 
                       'Shortness of Breath', models$names)
models$names <- ifelse(models$names == 'Loss_appetite', 
                       'Loss of Appetite', models$names)

# Create figure
colors <- c('All Symptoms' = '#41afaa', 'All Symptoms + Age' = '#c2a3fd', 
            'All Symptoms + Age + Time' = '#896deb'
)


models <- left_join(models, symptom_order_rsv, by = c('names'))

main_rsv_fig <- ggplot(data=models, aes(x=reorder(names, num), y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(aes(group = `Model`, color = `Model`, fill = `Model`),
                  position = position_dodge2(width = 0.4), size = 0.75, alpha = 0.7) + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Odds Ratio (95% CI)") +
  ggtitle("RSV Regression Results by Model") +
  theme_minimal() +  # use a white background
  scale_color_manual(values = colors) +
  forrest_theme + theme(plot.title = element_text(size=18, hjust = 0))

main_rsv_fig


figure_2_sub_1 <- cowplot::plot_grid(age_overall_flu_fig, geo_flu_fig,
                                     nrow = 2, rel_heights = c(0.6, 1))


figure_2_sub_2


age_overall_rsv_fig





# Create figure, with labels
figure_2 <- cowplot::plot_grid(percent_flu, age_flu_fig,
                               percent_rsv, main_rsv_fig, age_rsv_fig,
                               nrow = 2,
                               labels = c("a", "b", "c", "d", "e", "f"),
                               label_size = 20)

# Save figure 
ggsave('./figs/figure_2.jpg', plot = figure_2, height = 12, width = 20)



















age_rsv_fig <- ggplot(data=models, aes(x=reorder(names, num), y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(aes(group = `Age Group`, color = `Age Group`, fill = `Age Group`),
                  position = position_dodge2(width = 0.4), size = 0.75, alpha = 0.8) + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Odds Ratio (95% CI)") +
  theme_minimal() +  # use a white background
  scale_color_manual(values = colors) + theme(legend.position = 'bottom') +
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 20),
        axis.text = element_text(size=20),
        axis.title = element_text(size=24),
        legend.position = "right",
        legend.box="vertical",
        legend.margin=margin(),
        strip.background = element_blank(),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.key.size = unit(1, 'cm'),
        strip.text = element_text(size = 16),
        plot.title = element_text(size=28, hjust = 0.5))


age_rsv_fig










age_time_model_all <- glm(flu ~ fever + cough + sore_throat +
                            myalgia + short_breath + hypoxemia + nausea_vom + 
                            bronchitis + chest_pain + diarrhea + 
                            fatigue + headache + congestion + sneezing + 
                            as.factor(age_group) + as.factor(month),
                          family = binomial(link = "logit"),
                          data = flu_train)
summary(age_time_model_all)

# Model with all symptoms modified by age
age_model <- glm(flu ~ as.factor(age_group) + fever*as.factor(age_group) + 
                   cough*as.factor(age_group) + sore_throat*as.factor(age_group) +
                   myalgia*as.factor(age_group) + short_breath*as.factor(age_group) + 
                   hypoxemia*as.factor(age_group) +
                   nausea_vom*as.factor(age_group) + bronchitis*as.factor(age_group) + 
                   chest_pain*as.factor(age_group) + 
                   diarrhea*as.factor(age_group) + fatigue*as.factor(age_group) + 
                   headache*as.factor(age_group) + 
                   congestion*as.factor(age_group) + sneezing*as.factor(age_group),
                 family = binomial(link = "logit"),
                 data = flu_train)
summary(age_model)

# Model with all symptoms modified by age in addition to time
time_model <- glm(flu ~ as.factor(age_group) + fever*as.factor(age_group) + cough*as.factor(age_group) + 
                    sore_throat*as.factor(age_group) +
                    myalgia*as.factor(age_group) + short_breath*as.factor(age_group) + 
                    hypoxemia*as.factor(age_group) +
                    nausea_vom*as.factor(age_group) + bronchitis*as.factor(age_group) + 
                    chest_pain*as.factor(age_group) + 
                    diarrhea*as.factor(age_group) + fatigue*as.factor(age_group) + 
                    headache*as.factor(age_group) + 
                    congestion*as.factor(age_group) + sneezing*as.factor(age_group) + as.factor(month),
                  family = binomial(link = "logit"),
                  data = flu_train)
summary(time_model)



figure <- ggplot(data=model_all, aes(x=reorder(names, -num), y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(aes(group = age_group, color = age_group, fill = age_group),
                  position = position_dodge2(width = 0.6), size = 0.75, alpha = 0.8) + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Estimate (95% CI)") +
  theme_minimal() +  # use a white background
  scale_color_manual(values = colors) + theme(legend.position = 'bottom') +
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 20),
        axis.text = element_text(size=20),
        axis.title = element_text(size=24),
        legend.position = "bottom",
        legend.box="vertical",
        legend.margin=margin(),
        strip.background = element_blank(),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.key.size = unit(1, 'cm'),
        strip.text = element_text(size = 16),
        plot.title = element_text(size=28, hjust = 0.5))
