
####################
# 1. SET-UP SCRIPT #
####################

# Start with a clear environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(cowplot)
library(lubridate)
library(cutpointr)
library(ISOweek)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

# Load models and cut offs
flu_16_17 <- readRDS('./tmp/flu_p2_dat_proc.rds')

flu_16_17_country <- flu_16_17 |> group_by(week_date) |> 
  mutate(`Confirmed Influenza` = sum(flu*patient_count_imp, na.rm = TRUE),
         Fever = sum(fever*patient_count_imp, na.rm = TRUE),
         Cough = sum(cough*patient_count_imp, na.rm = TRUE),
         `Sore Throat` = sum(sore_throat*patient_count_imp, na.rm = TRUE)) |>
  distinct(week_date, `Confirmed Influenza`, Fever, Cough, `Sore Throat`, .keep_all = FALSE)



flu_16_17_country_long <- flu_16_17_country |> pivot_longer(!c(week_date), names_to = "type", values_to = "count")

library(zoo)
flu_16_17_country_long <- flu_16_17_country_long |> group_by(type) |>
  filter(week_date >= as.Date('2016-09-01')) |>
  filter(week_date <= as.Date('2020-08-20')) |>
  mutate(roll_mean = rollmean(count, k = 4, align = 'right', fill = NA))
  

flu_16_17_country_long_filt <- flu_16_17_country_long |>
  filter(type != 'Confirmed Influenza')

ggplot(flu_16_17_country_long_filt, aes(week_date, roll_mean, color = type)) +
  geom_line(size = 1.3, alpha = 0.9) + ylab('Number of Claims') + xlab('Date') +
  ggtitle("Symptom Claims, 2016-2020") + ylim(0, 500000) +
  theme_minimal() + theme(legend.position = 'bottom',
    axis.text.x = element_text(size = 12),
      axis.text.y = element_blank(),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      strip.text = element_text(size=14),
    plot.title = element_text(size=18),
      legend.text=element_text(size=12),
      legend.title=element_text(size=14)) +
  scale_color_manual(name = 'Symptom', values=c("#6BAED6", "#FD8D3C", '#74C476'))

[1] "#EDF8E9" "#BAE4B3" "#74C476" "#31A354" "#006D2C"

remove(flu_16_17)


# Collapse to country
flu_16_17_country <- flu_16_17_county |> group_by(week_date) |> 
  mutate(all_cause_sum = sum(all_cause),
         flu = sum(flu_sum),
         cough = sum(cough_sum),
         sore_throat = sum(sore_throat_sum),
         fever = sum(fever_sum) / all_cause_sum,
         flu_complex = sum(flu_complex_sum) / all_cause_sum) |>
  distinct(week_date, flu, ili, flu_combos, flu_fever, flu_complex, all_cause_sum, .keep_all = FALSE)
