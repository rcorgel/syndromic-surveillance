################################################################################
# File Name: subset_data_p1                                                    #
#                                                                              #
# Purpose:   Subset Influenza, COVID-19, and RSV medical claims data.          #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Subset Influenza data                                          #
#            3. Subset COVID-19 data                                           #
#            3. Subset RSV data                                                #
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

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

######################
# 2. SUBSET FLU DATA #
######################

# Load processed data
flu <- readRDS('./tmp/flu_p1_dat_proc.rds')

# Create year and month variables
flu$year <- year(flu$month_date)
flu$month <- month(flu$month_date)

# Bring data back to season level
flu_16_17 <- flu |>
  filter(month_date > as.Date("2016-08-31")) |>
  filter(month_date < as.Date("2017-09-01"))
flu_17_18 <- flu |>
  filter(month_date > as.Date("2017-08-31")) |>
  filter(month_date < as.Date("2018-09-01"))
flu_18_19 <- flu |>
  filter(month_date > as.Date("2018-08-31")) |>
  filter(month_date < as.Date("2019-09-01"))
flu_19_20 <- flu |>
  filter(month_date > as.Date("2019-08-31")) |>
  filter(month_date < as.Date("2020-09-01"))
# Check that data was subset correctly
(nrow(flu_16_17) + nrow(flu_17_18) + 
         nrow(flu_18_19) + nrow(flu_19_20) == nrow(flu))

# Restrict data to flu season months (September - April)
# flu_16_17 <- flu_16_17 |> filter(!(month > 3 & month < 10))
# flu_17_18 <- flu_17_18 |> filter(!(month > 3 & month < 10))
# flu_18_19 <- flu_18_19 |> filter(!(month > 3 & month < 10))
# flu_19_20 <- flu_19_20 |> filter(!(month > 3 & month < 10))

# Split into cases and non-cases
# 16-17
cases_16_17 <- flu_16_17[flu_16_17$flu == 1, ]
non_cases_16_17 <- flu_16_17[flu_16_17$flu == 0, ]
remove(flu_16_17) # remove season data set
# 17-18
cases_17_18 <- flu_17_18[flu_17_18$flu == 1, ]
non_cases_17_18 <- flu_17_18[flu_17_18$flu == 0, ]
remove(flu_17_18) # remove season data set
# 18-19
cases_18_19 <- flu_18_19[flu_18_19$flu == 1, ]
non_cases_18_19 <- flu_18_19[flu_18_19$flu == 0, ]
remove(flu_18_19) # remove season data set
# 19-20
cases_19_20 <- flu_19_20[flu_19_20$flu == 1, ]
non_cases_19_20 <- flu_19_20[flu_19_20$flu == 0, ]
remove(flu_19_20) # remove season data set
# Check that data was subset correctly
(nrow(cases_16_17) + nrow(non_cases_16_17) + nrow(cases_17_18) + 
    nrow(non_cases_17_18) + nrow(cases_18_19) + nrow(non_cases_18_19) + 
    nrow(cases_19_20) + nrow(non_cases_19_20) == nrow(flu))
remove(flu)

# Expand case data and save
cases_exp_sub_16_17 <- cases_16_17 |> uncount(patient_count_imp) |> 
  group_by(year_month) |> sample_n(10000)
cases_exp_sub_17_18 <- cases_17_18 |> uncount(patient_count_imp) |> 
  group_by(year_month) |> sample_n(10000)
cases_exp_sub_18_19 <- cases_18_19 |> uncount(patient_count_imp) |> 
  group_by(year_month) |> sample_n(10000)
cases_exp_sub_19_20 <- cases_19_20 |> uncount(patient_count_imp) |> 
  group_by(year_month) |> sample_n(10000)

# Save expanded case data
saveRDS(cases_exp_sub_16_17, './tmp/cases_exp_sub_16_17_flu.rds') 
saveRDS(cases_exp_sub_17_18, './tmp/cases_exp_sub_17_18_flu.rds') 
saveRDS(cases_exp_sub_18_19, './tmp/cases_exp_sub_18_19_flu.rds') 
saveRDS(cases_exp_sub_19_20, './tmp/cases_exp_sub_19_20_flu.rds') 
remove(cases_16_17, cases_17_18, cases_18_19, cases_19_20, 
       cases_exp_sub_16_17, cases_exp_sub_17_18, cases_exp_sub_18_19, cases_exp_sub_19_20)

# Expand and subset non-case data
non_cases_exp_sub_16_17 <- non_cases_16_17 |> uncount(patient_count_imp) |> 
                           group_by(year_month) |> sample_n(10000)
non_cases_exp_sub_17_18 <- non_cases_17_18 |> uncount(patient_count_imp) |> 
                           group_by(year_month) |> sample_n(10000)
non_cases_exp_sub_18_19 <- non_cases_18_19 |> uncount(patient_count_imp) |> 
                           group_by(year_month) |> sample_n(10000)
non_cases_exp_sub_19_20 <- non_cases_19_20 |> uncount(patient_count_imp) |> 
                           group_by(year_month) |> sample_n(10000)

# Save expanded non-case data
saveRDS(non_cases_exp_sub_16_17, './tmp/non_cases_exp_sub_16_17_flu.rds') 
saveRDS(non_cases_exp_sub_17_18, './tmp/non_cases_exp_sub_17_18_flu.rds') 
saveRDS(non_cases_exp_sub_18_19, './tmp/non_cases_exp_sub_18_19_flu.rds') 
saveRDS(non_cases_exp_sub_19_20, './tmp/non_cases_exp_sub_19_20_flu.rds') 
rm(list = ls())

########################
# 3. SUBSET COVID DATA #
########################

# Load processed data
covid <- readRDS('./tmp/covid_p1_dat_proc.rds')

# Create year and month variables
covid$year <- year(covid$month_date)
covid$month <- month(covid$month_date)

# Bring data back to season level
covid_wild <- covid |>
  filter(month_date > as.Date("2019-12-31")) |>
  filter(month_date < as.Date("2021-04-01"))
covid_alpha <- covid |>
  filter(month_date > as.Date("2021-03-31")) |>
  filter(month_date < as.Date("2021-07-01"))
covid_delta <- covid |>
  filter(month_date > as.Date("2021-06-30")) |>
  filter(month_date < as.Date("2022-01-01"))
covid_omicron <- covid |>
  filter(month_date > as.Date("2021-12-31")) |>
  filter(month_date < as.Date("2024-01-01"))
# Check that data was subset correctly
(nrow(covid_wild) + nrow(covid_alpha) + 
    nrow(covid_delta) + nrow(covid_omicron) == nrow(covid))

# Split into cases and non-cases
# 16-17
cases_wild <- covid_wild[covid_wild$covid == 1, ]
non_cases_wild <- covid_wild[covid_wild$covid == 0, ]
remove(covid_wild) # remove season data set
# 17-18
cases_alpha <- covid_alpha[covid_alpha$covid == 1, ]
non_cases_alpha <- covid_alpha[covid_alpha$covid == 0, ]
remove(covid_alpha) # remove season data set
# 18-19
cases_delta <- covid_delta[covid_delta$covid == 1, ]
non_cases_delta <- covid_delta[covid_delta$covid == 0, ]
remove(covid_delta) # remove season data set
# 19-20
cases_omicron <- covid_omicron[covid_omicron$covid == 1, ]
non_cases_omicron <- covid_omicron[covid_omicron$covid == 0, ]
remove(covid_omicron) # remove season data set
# Check that data was subset correctly
(nrow(cases_wild) + nrow(non_cases_wild) + nrow(cases_alpha) + 
    nrow(non_cases_alpha) + nrow(cases_delta) + nrow(non_cases_delta) + 
    nrow(cases_omicron) + nrow(non_cases_omicron) == nrow(covid))
remove(covid)

# Expand case data and save
cases_exp_sub_wild <- cases_wild |> uncount(patient_count_imp) |> 
  filter(month_date > as.Date("2020-02-28")) |> # remove pre-March 2020
  group_by(year_month) |> sample_n(9230)
cases_exp_sub_alpha <- cases_alpha |> uncount(patient_count_imp) |> 
  group_by(year_month) |> sample_n(40000)
cases_exp_sub_delta <- cases_delta |> uncount(patient_count_imp) |> 
  group_by(year_month) |> sample_n(20000)
cases_exp_sub_omicron <- cases_omicron |> uncount(patient_count_imp) |> 
  filter(month_date < as.Date("2023-09-01")) |> # remove post-Sept 2023
  group_by(year_month) |> sample_n(6000)

# Save expanded case data
saveRDS(cases_exp_sub_wild, './tmp/cases_exp_sub_wild_covid.rds') 
saveRDS(cases_exp_sub_alpha, './tmp/cases_exp_sub_alpha_covid.rds') 
saveRDS(cases_exp_sub_delta, './tmp/cases_exp_sub_delta_covid.rds') 
saveRDS(cases_exp_sub_omicron, './tmp/cases_exp_sub_omicron_covid.rds') 
remove(cases_wild, cases_alpha, cases_delta, cases_omicron, 
       cases_exp_sub_wild, cases_exp_sub_alpha, cases_exp_sub_delta, cases_exp_sub_omicron)

# Expand and subset non-case data
non_cases_exp_sub_wild <- non_cases_wild |> uncount(patient_count_imp) |> 
  filter(month_date > as.Date("2020-02-28")) |> # remove pre-March 2020
  group_by(year_month) |> sample_n(9230)
remove(non_cases_wild)
non_cases_exp_sub_alpha <- non_cases_alpha |> uncount(patient_count_imp) |> 
  group_by(year_month) |> sample_n(40000)
remove(non_cases_alpha)
non_cases_exp_sub_delta <- non_cases_delta |> uncount(patient_count_imp) |> 
  group_by(year_month) |> sample_n(20000)
remove(non_cases_delta)
# Split omicron into two so it can be processed (R keeps crashing)
non_cases_omicron_1 <- non_cases_omicron |> filter(month_date < as.Date("2022-11-01"))
non_cases_omicron_2 <- non_cases_omicron |> filter(month_date > as.Date("2022-10-31"))
remove(non_cases_omicron)
non_cases_exp_sub_omicron_1 <- non_cases_omicron_1 |> uncount(patient_count_imp) |> 
  group_by(year_month) |> sample_n(6000)
non_cases_exp_sub_omicron_2 <- non_cases_omicron_2 |> uncount(patient_count_imp) |> 
  filter(month_date < as.Date("2023-09-01")) |> # remove post-Sept 2023
  group_by(year_month) |> sample_n(6000)
non_cases_exp_sub_omicron <- rbind(non_cases_exp_sub_omicron_1, non_cases_exp_sub_omicron_2)

# Save expanded non-case data
saveRDS(non_cases_exp_sub_wild, './tmp/non_cases_exp_sub_wild_covid.rds') 
saveRDS(non_cases_exp_sub_alpha, './tmp/non_cases_exp_sub_alpha_covid.rds') 
saveRDS(non_cases_exp_sub_delta, './tmp/non_cases_exp_sub_delta_covid.rds') 
saveRDS(non_cases_exp_sub_omicron, './tmp/non_cases_exp_sub_omicron_covid.rds') 
rm(list = ls())

######################
# 2. SUBSET RSV DATA #
######################

# Load processed data
rsv <- readRDS('./tmp/rsv_p1_dat_proc.rds')

# Create year and month variables
rsv$year <- year(rsv$month_date)
rsv$month <- month(rsv$month_date)

# Bring data back to season level
rsv_16_17 <- rsv |>
  filter(month_date > as.Date("2016-08-31")) |>
  filter(month_date < as.Date("2017-09-01"))
rsv_17_18 <- rsv |>
  filter(month_date > as.Date("2017-08-31")) |>
  filter(month_date < as.Date("2018-09-01"))
rsv_18_19 <- rsv |>
  filter(month_date > as.Date("2018-08-31")) |>
  filter(month_date < as.Date("2019-09-01"))
rsv_19_20 <- rsv |>
  filter(month_date > as.Date("2019-08-31")) |>
  filter(month_date < as.Date("2020-09-01"))
# Check that data was subset correctly
(nrow(rsv_16_17) + nrow(rsv_17_18) + 
    nrow(rsv_18_19) + nrow(rsv_19_20) == nrow(rsv))

# Restrict data to rsv season months (September - April)
# rsv_16_17 <- rsv_16_17 |> filter(!(month > 3 & month < 10))
# rsv_17_18 <- rsv_17_18 |> filter(!(month > 3 & month < 10))
# rsv_18_19 <- rsv_18_19 |> filter(!(month > 3 & month < 10))
# rsv_19_20 <- rsv_19_20 |> filter(!(month > 3 & month < 10))

# Split into cases and non-cases
# 16-17
cases_16_17 <- rsv_16_17[rsv_16_17$rsv == 1, ]
non_cases_16_17 <- rsv_16_17[rsv_16_17$rsv == 0, ]
remove(rsv_16_17) # remove season data set
# 17-18
cases_17_18 <- rsv_17_18[rsv_17_18$rsv == 1, ]
non_cases_17_18 <- rsv_17_18[rsv_17_18$rsv == 0, ]
remove(rsv_17_18) # remove season data set
# 18-19
cases_18_19 <- rsv_18_19[rsv_18_19$rsv == 1, ]
non_cases_18_19 <- rsv_18_19[rsv_18_19$rsv == 0, ]
remove(rsv_18_19) # remove season data set
# 19-20
cases_19_20 <- rsv_19_20[rsv_19_20$rsv == 1, ]
non_cases_19_20 <- rsv_19_20[rsv_19_20$rsv == 0, ]
remove(rsv_19_20) # remove season data set
# Check that data was subset correctly
(nrow(cases_16_17) + nrow(non_cases_16_17) + nrow(cases_17_18) + 
    nrow(non_cases_17_18) + nrow(cases_18_19) + nrow(non_cases_18_19) + 
    nrow(cases_19_20) + nrow(non_cases_19_20) == nrow(rsv))
remove(rsv)

# Expand case data and save
cases_exp_sub_16_17 <- cases_16_17 |> uncount(patient_count_imp) |> 
  group_by(year_month) |> sample_n(1000)
cases_exp_sub_17_18 <- cases_17_18 |> uncount(patient_count_imp) |> 
  group_by(year_month) |> sample_n(1000)
cases_exp_sub_18_19 <- cases_18_19 |> uncount(patient_count_imp) |> 
  group_by(year_month) |> sample_n(1000)
cases_exp_sub_19_20 <- cases_19_20 |> uncount(patient_count_imp) |> 
  group_by(year_month) |> sample_n(1000)

# Save expanded case data
saveRDS(cases_exp_sub_16_17, './tmp/cases_exp_sub_16_17_rsv.rds') 
saveRDS(cases_exp_sub_17_18, './tmp/cases_exp_sub_17_18_rsv.rds') 
saveRDS(cases_exp_sub_18_19, './tmp/cases_exp_sub_18_19_rsv.rds') 
saveRDS(cases_exp_sub_19_20, './tmp/cases_exp_sub_19_20_rsv.rds') 
remove(cases_16_17, cases_17_18, cases_18_19, cases_19_20, 
       cases_exp_sub_16_17, cases_exp_sub_17_18, cases_exp_sub_18_19, cases_exp_sub_19_20)

# Expand and subset non-case data
non_cases_exp_sub_16_17 <- non_cases_16_17 |> uncount(patient_count_imp) |> 
  group_by(year_month) |> sample_n(1000)
non_cases_exp_sub_17_18 <- non_cases_17_18 |> uncount(patient_count_imp) |> 
  group_by(year_month) |> sample_n(1000)
non_cases_exp_sub_18_19 <- non_cases_18_19 |> uncount(patient_count_imp) |> 
  group_by(year_month) |> sample_n(1000)
non_cases_exp_sub_19_20 <- non_cases_19_20 |> uncount(patient_count_imp) |> 
  group_by(year_month) |> sample_n(1000)

# Save expanded non-case data
saveRDS(non_cases_exp_sub_16_17, './tmp/non_cases_exp_sub_16_17_rsv.rds') 
saveRDS(non_cases_exp_sub_17_18, './tmp/non_cases_exp_sub_17_18_rsv.rds') 
saveRDS(non_cases_exp_sub_18_19, './tmp/non_cases_exp_sub_18_19_rsv.rds') 
saveRDS(non_cases_exp_sub_19_20, './tmp/non_cases_exp_sub_19_20_rsv.rds') 
rm(list = ls())

################################################################################
################################################################################
