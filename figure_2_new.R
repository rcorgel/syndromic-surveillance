################################################################################
# File Name: figure_2                                                          #
#                                                                              #
# Purpose:   Create figure 2 for the manuscript.                               #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Create figure 2                                                #
#               a. Create profile plots                                        #
#               b. Define candidate models                                     #
#               c. Create age forest plots                                     #
#               d. Create time forest plots                                    #
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
library(cowplot)
library(lme4)
library(car)
library(mgcv)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

######################
# 2. CREATE FIGURE 2 #
######################

# Load flu data
flu_2016 <- readRDS('./tmp/flu_filt_symp_p2_bal_2016.rds')
flu_2017 <- readRDS('./tmp/flu_filt_symp_p2_bal_2017.rds')
flu_2018 <- readRDS('./tmp/flu_filt_symp_p2_bal_2018.rds')
flu_2019 <- readRDS('./tmp/flu_filt_symp_p2_bal_2019.rds')

# Combine flu data to one file
flu <- rbind(flu_2016, flu_2017, flu_2018, flu_2019)
remove(flu_2016, flu_2017, flu_2018, flu_2019)

# Calculate epidemiological week
flu$epi_week <- epiweek(flu$year_week_dt)

# Change reference week to week 26
flu$epi_week_adj <- ((flu$epi_week - 26) %% 52) + 1

# Convert to date
flu$year_week_dt <- as.Date(flu$year_week_dt)

# CHange var to urban_south
flu$urban_south <- ifelse(flu$urban_code_binary == 1 & flu$south == 1, 'Urban South', '')
flu$urban_south <- ifelse(flu$urban_code_binary == 0 & flu$south == 1, 'Rural South', flu$urban_south)
flu$urban_south <- ifelse(flu$urban_code_binary == 1 & flu$south == 0, 'Urban North', flu$urban_south)
flu$urban_south <- ifelse(flu$urban_code_binary == 0 & flu$south == 0, 'Rural North', flu$urban_south)



flu_1 <- flu |> dplyr::filter(flu == 1)
fever_1 <- flu |> dplyr::filter(fever == 1)
p_flu_fever <- sum(fever_1$flu) / nrow(fever_1)
p_fever_flu <- sum(flu_1$fever) / 500000
p_flu <- 0.5
p_fever <- sum(flu$fever) / 1000000
(p_fever_flu * p_flu) / p_fever


# Add nrevss flu data
nrevss_nat <- readRDS('./tmp/nrevss_national.rds')
nrevss_nat$week_date_ALT <- nrevss_nat$week_date + 0
nrevss_state <- readRDS('./tmp/nrevss_state.rds')
nrevss_state$week_date_ALT <- nrevss_state$week_date + 0
flu <- left_join(flu, nrevss_nat[, c(2, 9)], by = c('year_week_dt' = 'week_date_ALT'))
flu <- flu |> rename('nrevss_nat' = 'value')
flu <- left_join(flu, nrevss_state[, c(2, 9, 4)], 
                 by = c('year_week_dt' = 'week_date_ALT',
                        'state_fips' = 'state_fips'))
flu <- flu |> rename('nrevss_state' = 'value')

flu <- flu |> mutate(nrevss_state = ifelse(is.na(nrevss_state), nrevss_nat, nrevss_state))

flu_filt <- flu |> dplyr::filter(!is.na(nrevss_state)) |>
  dplyr::filter(!is.na(nrevss_nat))

# Split the data randomly into testing and training
sample <- sample(c(TRUE, FALSE), nrow(flu_filt), replace=TRUE, prob=c(0.50, 0.50))
flu_train  <- flu_filt[sample, ]
flu_test   <- flu_filt[!sample, ]

# sample <- sample(c(TRUE, FALSE), nrow(flu_sample), replace=TRUE, prob=c(0.50, 0.50))
# flu_train  <- flu_sample[sample, ]
# flu_test   <- flu_sample[!sample, ]

# PLOTS #
# Calculate Percents
flu_train$count <- 1
flu_symp <- flu_train |> 
  group_by(flu) |>
  mutate(total = sum(count),
         fever = sum(fever) / total,
         myalgia = sum(myalgia) / total,
         cough = sum(cough) / total,
         sore_throat = sum(sore_throat) / total, 
         short_breath = sum(short_breath) / total, 
         hypoxemia = sum(hypoxemia) / total, 
         chest_pain = sum(chest_pain) / total, 
         bronchitis = sum(bronchitis) / total,
         nausea_vom = sum(nausea_vom) / total, 
         diarrhea = sum(diarrhea) / total, 
         fatigue = sum(fatigue) / total, 
         headache = sum(headache) / total,
         congestion = sum(congestion) / total, 
         sneezing = sum(sneezing) / total) |>
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
                                'Short of Breath', flu_symp_long$symptom)

# Change case names
flu_symp_long$flu <- ifelse(flu_symp_long$flu == 0, 
                            'Non-Influenza Cases', flu_symp_long$flu)
flu_symp_long$flu <- ifelse(flu_symp_long$flu == 1, 
                            'Influenza Cases', flu_symp_long$flu)

# Merge on rank
flu_symp_long_flu <- flu_symp_long[flu_symp_long$flu == 'Influenza Cases', c(2, 3)]
flu_symp_long_flu$rank <- rank(flu_symp_long_flu$percent, na.last = TRUE)
flu_symp_long <- left_join(flu_symp_long, flu_symp_long_flu[, c(1, 3)], by = 'symptom')

# # Load rsv data
# rsv_dat_p1 <- readRDS('./tmp/rsv_p1_dat_proc_symp.rds')
# 
# # Calculate Percents
# rsv_symp <- rsv_dat_p1 |> 
#   group_by(rsv) |>
#   mutate(fever = sum(fever*patient_count_imp) / sum(patient_count_imp),
#          cough = sum(cough*patient_count_imp) / sum(patient_count_imp),
#          sore_throat = sum(sore_throat*patient_count_imp) / sum(patient_count_imp), 
#          short_breath = sum(short_breath*patient_count_imp) / sum(patient_count_imp), 
#          hypoxemia = sum(hypoxemia*patient_count_imp) / sum(patient_count_imp), 
#          bronchitis = sum(bronchitis*patient_count_imp) / sum(patient_count_imp),
#          loss_appetite = sum(loss_appetite*patient_count_imp) / sum(patient_count_imp), 
#          fatigue = sum(fatigue*patient_count_imp) / sum(patient_count_imp), 
#          headache = sum(headache*patient_count_imp) / sum(patient_count_imp),
#          congestion = sum(congestion*patient_count_imp) / sum(patient_count_imp), 
#          sneezing = sum(sneezing*patient_count_imp) / sum(patient_count_imp)) |>
#   distinct(rsv, fever, cough, sore_throat, 
#            short_breath, hypoxemia, bronchitis,
#            loss_appetite, fatigue, headache,
#            congestion, sneezing)
# 
# # Convert to long
# rsv_symp_long <- rsv_symp |>
#   pivot_longer(!c(rsv), names_to = "symptom", values_to = "percent") |>
#   mutate(symptom = str_to_title(symptom))
# 
# # Change symptom names
# rsv_symp_long$symptom <- ifelse(rsv_symp_long$symptom == 'Sore_throat', 
#                                 'Sore Throat', rsv_symp_long$symptom)
# rsv_symp_long$symptom <- ifelse(rsv_symp_long$symptom == 'Loss_appetite', 
#                                 'Loss of Appetite', rsv_symp_long$symptom)
# rsv_symp_long$symptom <- ifelse(rsv_symp_long$symptom == 'Short_breath', 
#                                 'Shortness of Breath', rsv_symp_long$symptom)
# 
# # Change case names
# rsv_symp_long$rsv <- ifelse(rsv_symp_long$rsv == 0, 
#                             'Non-RSV Cases', rsv_symp_long$rsv)
# rsv_symp_long$rsv <- ifelse(rsv_symp_long$rsv == 1, 
#                             'RSV Cases', rsv_symp_long$rsv)
# 
# # Merge on rank
# rsv_symp_long_rsv <- rsv_symp_long[rsv_symp_long$rsv == 'RSV Cases', c(2, 3)]
# rsv_symp_long_rsv$rank <- rank(rsv_symp_long_rsv$percent, na.last = TRUE)
# rsv_symp_long <- left_join(rsv_symp_long, rsv_symp_long_rsv[, c(1, 3)], by = 'symptom')

# Set plot theme
percent_theme <- theme(legend.position = "inside",
                       legend.position.inside = c(0.845, 0.74),
                       axis.text = element_text(size=16),
                       axis.title.y = element_text(size=18),
                       axis.title.x = element_blank(),
                       strip.background = element_blank(),
                       plot.title = element_text(size=20),
                       axis.text.x = element_text(angle = 35, vjust = 1, 
                                                  hjust=1, color = 'black'),
                       legend.title = element_text(size = 16),
                       legend.text = element_text(size = 14),
                       legend.box.background = element_rect(colour = "black", 
                                                            fill = 'white'))

# Set colors for influenza
flu_colors <- c(
  "Influenza Cases" = "#CC6594",
  "Non-Influenza Cases" = "darkgray")

# Create flight
percent_flu <- ggplot(flu_symp_long, 
                      aes(x = fct_reorder(symptom, desc(rank)), 
                          y = percent)) +
  geom_col(aes(fill = flu, group = flu, color = flu), width = 0.6, alpha = 0.7, position = position_dodge()) +  
  theme_minimal() + percent_theme +
  labs(title = "Observed Syndromic Profiles") +
  ylab('Patient Proportion') + xlab('') + expand_limits(y = c(0, 0.5)) +
  scale_fill_manual('Disease Status', values = flu_colors)+
  scale_color_manual('Disease Status', values = flu_colors)+
  labs(color = "", lty = "") 

# Display
percent_flu

# # Set colors for RSV
# rsv_colors <- c(
#   "RSV Cases" = "#41afaa",
#   "Non-RSV Cases" = "darkgray")
# 
# # Create flight
# percent_rsv <- ggplot(rsv_symp_long, 
#                       aes(x = fct_reorder(symptom, desc(rank)), 
#                           y = percent)) +
#   geom_line(aes(color = rsv, lty = rsv, group = rsv), linewidth = 2, alpha = 0.6) +  
#   geom_point(aes(color = rsv), size = 4, show.legend = TRUE, alpha = 0.8) +
#   theme_minimal() + percent_theme +
#   labs(title = "RSV Syndromic Profile\n") +
#   ylab('Proportion of Patients') + xlab('') + expand_limits(y = c(0, 0.5)) +
#   scale_color_manual(values = rsv_colors,
#                      breaks=c("RSV Cases", "Non-RSV Cases"))+
#   scale_linetype_manual(values = c('solid', 'dashed'),
#                         breaks=c("RSV Cases", "Non-RSV Cases"))+
#   labs(color = "", lty = "") 
# 
# # Display
# percent_rsv

# Fever only model
fev_model <- glm(flu ~ fever,
                 family = binomial(link = "logit"),
                 data = flu_train)
summary(fev_model)
saveRDS(fev_model, "./tmp/fever_model.rds")

# Symptoms model
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
saveRDS(symp_model, "./tmp/symp_model.rds")

# Demography model
demo_model <- glm(flu ~ fever*as.factor(age_grp) + 
                    cough*as.factor(age_grp) + 
                    sore_throat*as.factor(age_grp) +
                    myalgia*as.factor(age_grp) + 
                    short_breath*as.factor(age_grp) + 
                    hypoxemia*as.factor(age_grp) + 
                    nausea_vom*as.factor(age_grp) + 
                    bronchitis*as.factor(age_grp) + 
                    chest_pain*as.factor(age_grp) + 
                    diarrhea*as.factor(age_grp) + 
                    fatigue*as.factor(age_grp) + 
                    headache*as.factor(age_grp) + 
                    congestion*as.factor(age_grp) + 
                    sneezing*as.factor(age_grp) +
                    fever*as.factor(patient_gender_code) + 
                    cough*as.factor(patient_gender_code) + 
                    sore_throat*as.factor(patient_gender_code) +
                    myalgia*as.factor(patient_gender_code) + 
                    short_breath*as.factor(patient_gender_code) + 
                    hypoxemia*as.factor(patient_gender_code) + 
                    nausea_vom*as.factor(patient_gender_code) + 
                    bronchitis*as.factor(patient_gender_code) + 
                    chest_pain*as.factor(patient_gender_code) + 
                    diarrhea*as.factor(patient_gender_code) + 
                    fatigue*as.factor(patient_gender_code) + 
                    headache*as.factor(patient_gender_code) + 
                    congestion*as.factor(patient_gender_code) + 
                    sneezing*as.factor(patient_gender_code) +
                    as.factor(urban_south),
                  family = binomial(link = "logit"),
                  data = flu_train)
summary(demo_model)
saveRDS(demo_model, "./tmp/demo_model.rds")

# Disease model
time_dis_model <- gam(flu ~ fever*as.factor(age_grp) + 
                        cough*as.factor(age_grp) + 
                        sore_throat*as.factor(age_grp) +
                        myalgia*as.factor(age_grp) + 
                        short_breath*as.factor(age_grp) + 
                        hypoxemia*as.factor(age_grp) + 
                        nausea_vom*as.factor(age_grp) + 
                        bronchitis*as.factor(age_grp) + 
                        chest_pain*as.factor(age_grp) + 
                        diarrhea*as.factor(age_grp) + 
                        fatigue*as.factor(age_grp) + 
                        headache*as.factor(age_grp) + 
                        congestion*as.factor(age_grp) + 
                        sneezing*as.factor(age_grp) +
                        fever*as.factor(patient_gender_code) + 
                        cough*as.factor(patient_gender_code) + 
                        sore_throat*as.factor(patient_gender_code) +
                        myalgia*as.factor(patient_gender_code) + 
                        short_breath*as.factor(patient_gender_code) + 
                        hypoxemia*as.factor(patient_gender_code) + 
                        nausea_vom*as.factor(patient_gender_code) + 
                        bronchitis*as.factor(patient_gender_code) + 
                        chest_pain*as.factor(patient_gender_code) + 
                        diarrhea*as.factor(patient_gender_code) + 
                        fatigue*as.factor(patient_gender_code) + 
                        headache*as.factor(patient_gender_code) + 
                        congestion*as.factor(patient_gender_code) + 
                        sneezing*as.factor(patient_gender_code) +
                        s(nrevss_state, bs = "cc", k = 20),
                      family = binomial(link = "logit"),
                      data = flu_train)
summary(time_dis_model)
saveRDS(time_dis_model, "./tmp/time_dis_model.rds")

# Disease Geography model
time_dis_geo_model <- gam(flu ~ fever*as.factor(age_grp) + 
                        cough*as.factor(age_grp) + 
                        sore_throat*as.factor(age_grp) +
                        myalgia*as.factor(age_grp) + 
                        short_breath*as.factor(age_grp) + 
                        hypoxemia*as.factor(age_grp) + 
                        nausea_vom*as.factor(age_grp) + 
                        bronchitis*as.factor(age_grp) + 
                        chest_pain*as.factor(age_grp) + 
                        diarrhea*as.factor(age_grp) + 
                        fatigue*as.factor(age_grp) + 
                        headache*as.factor(age_grp) + 
                        congestion*as.factor(age_grp) + 
                        sneezing*as.factor(age_grp) +
                        fever*as.factor(patient_gender_code) + 
                        cough*as.factor(patient_gender_code) + 
                        sore_throat*as.factor(patient_gender_code) +
                        myalgia*as.factor(patient_gender_code) + 
                        short_breath*as.factor(patient_gender_code) + 
                        hypoxemia*as.factor(patient_gender_code) + 
                        nausea_vom*as.factor(patient_gender_code) + 
                        bronchitis*as.factor(patient_gender_code) + 
                        chest_pain*as.factor(patient_gender_code) + 
                        diarrhea*as.factor(patient_gender_code) + 
                        fatigue*as.factor(patient_gender_code) + 
                        headache*as.factor(patient_gender_code) + 
                        congestion*as.factor(patient_gender_code) + 
                        sneezing*as.factor(patient_gender_code) +
                        as.factor(urban_south) +
                        s(nrevss_state, bs = "cc", k = 20, by = as.factor(urban_south)),
                      family = binomial(link = "logit"),
                      data = flu_train)
summary(time_dis_geo_model)
saveRDS(time_dis_geo_model, "./tmp/time_dis_geo_model.rds")


time_state_model <- gam(flu ~ fever*as.factor(age_grp) + 
                            cough*as.factor(age_grp) + 
                            sore_throat*as.factor(age_grp) +
                            myalgia*as.factor(age_grp) + 
                            short_breath*as.factor(age_grp) + 
                            hypoxemia*as.factor(age_grp) + 
                            nausea_vom*as.factor(age_grp) + 
                            bronchitis*as.factor(age_grp) + 
                            chest_pain*as.factor(age_grp) + 
                            diarrhea*as.factor(age_grp) + 
                            fatigue*as.factor(age_grp) + 
                            headache*as.factor(age_grp) + 
                            congestion*as.factor(age_grp) + 
                            sneezing*as.factor(age_grp) +
                            fever*as.factor(patient_gender_code) + 
                            cough*as.factor(patient_gender_code) + 
                            sore_throat*as.factor(patient_gender_code) +
                            myalgia*as.factor(patient_gender_code) + 
                            short_breath*as.factor(patient_gender_code) + 
                            hypoxemia*as.factor(patient_gender_code) + 
                            nausea_vom*as.factor(patient_gender_code) + 
                            bronchitis*as.factor(patient_gender_code) + 
                            chest_pain*as.factor(patient_gender_code) + 
                            diarrhea*as.factor(patient_gender_code) + 
                            fatigue*as.factor(patient_gender_code) + 
                            headache*as.factor(patient_gender_code) + 
                            congestion*as.factor(patient_gender_code) + 
                            sneezing*as.factor(patient_gender_code) +
                            s(epi_week_adj, bs = "cc", k = 20, by = as.factor(state_fips)),
                          family = binomial(link = "logit"),
                          data = flu_train)
summary(time_state_model)
saveRDS(time_dis_geo_model, "./tmp/time_state_model.rds")

# Briefly test the models
flu_test$pred_fev <- predict(fev_model, newdata = flu_test, type = "response", allow.new.levels = TRUE)
pROC::auc(pROC::roc(flu_test$flu, flu_test$pred_fev))
roc_fev <- pROC::roc(flu_test$flu, flu_test$pred_fev)
ggplot(flu_test) + geom_density(aes(x = pred_fev, group = as.factor(flu), color = as.factor(flu)))

flu_test$pred_symp <- predict(symp_model, newdata = flu_test, type = "response")
pROC::auc(pROC::roc(flu_test$flu, flu_test$pred_symp))
roc_symp <- pROC::roc(flu_test$flu, flu_test$pred_symp)
ggplot(flu_test) + geom_density(aes(x = pred_symp, group = as.factor(flu), color = as.factor(flu)))

flu_test$pred_demo <- predict(demo_model, newdata = flu_test, type = "response")
pROC::auc(pROC::roc(flu_test$flu, flu_test$pred_demo))
roc_demo <- pROC::roc(flu_test$flu, flu_test$pred_demo)
ggplot(flu_test) + geom_density(aes(x = pred_demo, group = as.factor(flu), color = as.factor(flu)))

flu_test$pred_dis <- predict(time_dis_model, newdata = flu_test, type = "response")
pROC::auc(pROC::roc(flu_test$flu, flu_test$pred_dis))
roc_dis <- pROC::roc(flu_test$flu, flu_test$pred_dis)
ggplot(flu_test) + geom_density(aes(x = pred_dis, group = as.factor(flu), color = as.factor(flu)))

flu_test$pred_dis_geo <- predict(time_dis_geo_model, newdata = flu_test, type = "response")
pROC::auc(pROC::roc(flu_test$flu, flu_test$pred_dis_geo))
roc_dis_geo <- pROC::roc(flu_test$flu, flu_test$pred_dis_geo)
ggplot(flu_test) + geom_density(aes(x = pred_dis_geo, group = as.factor(flu), color = as.factor(flu)))



time_dis_geo_model <- readRDS("./tmp/time_dis_geo_model.rds")

library(cutpointr)
cp <- cutpointr(flu_test[flu_test$age_grp == '>=65',], pred_dis_geo, flu, 
                method = maximize_metric, metric = youden)
summary(cp)

cp <- cutpointr(flu_test[flu_test$age_grp == '50-64',], pred_dis_geo, flu, 
                method = maximize_metric, metric = youden)
summary(cp)

cp <- cutpointr(flu_test[flu_test$age_grp == '18-49',], pred_dis_geo, flu, 
                method = maximize_metric, metric = youden)
summary(cp)

cp <- cutpointr(flu_test[flu_test$age_grp == '13-17',], pred_dis_geo, flu, 
                method = maximize_metric, metric = youden)
summary(cp)

cp <- cutpointr(flu_test[flu_test$age_grp == '5-12',], pred_dis_geo, flu, 
                method = maximize_metric, metric = youden)
summary(cp)

cp <- cutpointr(flu_test[flu_test$age_grp == '0-4',], pred_dis_geo, flu, 
                method = maximize_metric, metric = youden)
summary(cp)



# Create ROC Data
roc_data <- data.frame(
  Specificity = c(roc_fev$specificities, roc_symp$specificities, roc_demo$specificities, roc_dis_geo$specificities),
  Sensitivity = c(roc_fev$sensitivities, roc_symp$sensitivities, roc_demo$sensitivities, roc_dis_geo$sensitivities),
  Model = c(
    rep(paste0("Fever Model"), length(roc_fev$specificities)),
    rep(paste0("Symptom Model"), length(roc_symp$specificities)),
    rep(paste0("Demography Model"), length(roc_demo$specificities)),
    rep(paste0("NREVSS Model"), length(roc_dis_geo$specificities))
  )
)

roc_data$Model <- factor(roc_data$Model, levels = c("Fever Model",
                                           "Symptom Model",
                                           "Demography Model",
                                           "NREVSS Model"))

# Plot
roc_plot <- ggplot(roc_data, aes(x = 1 - Specificity, y = Sensitivity, color = Model)) +
  geom_line(linewidth = 3, alpha = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1.5) +
  scale_color_manual(values = c('#CC6594', '#d7642c', '#347DC1', '#45B649')) +
  labs(
    title = "ROC Curves by Model",
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme_minimal() +
  theme(legend.position =  "inside",
        legend.position.inside = c(0.85, 0.12),
        axis.text = element_text(size=14),
        axis.title = element_text(size=16),
        strip.background = element_blank(),
        legend.title = element_text(size=16),
        legend.text = element_text(size=14),
        plot.title = element_text(size=16, hjust = 0),
        legend.box="vertical",
        legend.box.background = element_rect(colour = "black", 
                                             fill = 'white')) +
  coord_equal()

flu_test_test <- flu_test |> ungroup() |>
  mutate(fever = ifelse(fever == 1, 0, fever),
         cough = ifelse(cough == 1, 0, cough),
         sore_throat = ifelse(sore_throat == 1, 0, sore_throat),
         myalgia = ifelse(myalgia == 1, 0,  myalgia),
         short_breath = ifelse(short_breath == 1, 0, short_breath),
         hypoxemia = ifelse(hypoxemia == 1, 0, hypoxemia),
         nausea_vom = ifelse(nausea_vom == 1, 0, nausea_vom),
         bronchitis = ifelse(bronchitis == 1, 0, bronchitis),
         chest_pain = ifelse(chest_pain == 1, 0, chest_pain),
         diarrhea = ifelse(diarrhea == 1, 0, diarrhea),
         fatigue = ifelse(fatigue == 1, 0, fatigue),
         headache = ifelse(headache == 1, 0, headache),
         congestion = ifelse(congestion == 1, 0, congestion),
         sneezing = ifelse(sneezing == 1, 0, sneezing))




plot_data <- plot(time_dis_model, pages = 0)

df <- data.frame(
  nrevss_state = plot_data[[1]]$x,
  fit = plot_data[[1]]$fit,
  se = plot_data[[1]]$se
)
df$lower <- df$fit - 1.96 * df$se
df$upper <- df$fit + 1.96 * df$se

# Convert to odds ratios
df$fit_or <- exp(df$fit)
df$lower_or <- exp(df$lower)
df$upper_or <- exp(df$upper)

test_plot <- ggplot(df, aes(x = nrevss_state, y = fit_or)) +
  geom_line(color = '#45B649', linewidth = 2) +
  geom_ribbon(aes(ymin = lower_or, ymax = upper_or), 
              alpha = 0.2, 
              fill = '#45B649') +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 1.5) +
  labs(
    title = 'Influenza Odds by State Test Positivity',
    x = "Influenza Test Positivity",
    y = "Odds Ratio"
  ) +
  scale_x_continuous(labels = scales::percent_format(scale = 1, accuracy = 1)) +
  theme_minimal() +
  theme(legend.position =  "bottom",
        legend.position.inside = c(0.70, 0.20),
        axis.text = element_text(size=16),
        axis.title = element_text(size=18),
        strip.background = element_blank(),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        plot.title = element_text(size=20, hjust = 0),
        legend.box="vertical",
        legend.box.background = element_rect(colour = "white", 
                                             fill = 'white'))



###############################
# C. CREATE AGE FORREST PLOTS #
###############################

#############
# INFLUENZA #
#############

# Run models trained on each age group
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
                       data = flu_train[flu_train$age_grp == '0-4',])

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
                        data = flu_train[flu_train$age_grp == '5-12',])

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
                         data = flu_train[flu_train$age_grp == '13-17',])

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
                         data = flu_train[flu_train$age_grp == '18-49',])

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
                         data = flu_train[flu_train$age_grp == '50-64',])

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
                      data = flu_train[flu_train$age_grp == '>=65',])

# Create data for the overall model
coefficients <- symp_model$coefficients
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
model_all$names <- ifelse(model_all$names == 'Chest_pain', 
                          'Chest Pain', model_all$names)
model_all$names <- ifelse(model_all$names == 'Short_breath', 
                          'Short of Breath', model_all$names)

# Remove the intercept
model_all_flu_overall <- model_all |> dplyr::filter(names != '(Intercept)')

# Save the order
symptom_order_flu <- model_all_flu_overall[, c(1, 4)]

#########################
# CREATE AGE GROUP PLOT #
#########################

# Create a data frame for each age model
coefficients <- train_model_0_4$coefficients
se <- summary(train_model_0_4)$coefficients[, 2]
names <- names(train_model_0_4$coefficients)
model_0_4 <- data.frame(names, coefficients, se)
model_0_4$age_group <- '0-4'

coefficients <- train_model_5_12$coefficients
se <- summary(train_model_5_12)$coefficients[, 2]
names <- names(train_model_5_12$coefficients)
model_5_12 <- data.frame(names, coefficients, se)
model_5_12$age_group <- '5-12'

coefficients <- train_model_13_17$coefficients
se <- summary(train_model_13_17)$coefficients[, 2]
names <- names(train_model_13_17$coefficients)
model_13_17 <- data.frame(names, coefficients, se)
model_13_17$age_group <- '13-17'

coefficients <- train_model_18_49$coefficients
se <- summary(train_model_18_49)$coefficients[, 2]
names <- names(train_model_18_49$coefficients)
model_18_49 <- data.frame(names, coefficients, se)
model_18_49$age_group <- '18-49'

coefficients <- train_model_50_64$coefficients
se <- summary(train_model_50_64)$coefficients[, 2]
names <- names(train_model_50_64$coefficients)
model_50_64 <- data.frame(names, coefficients, se)
model_50_64$age_group <- '50-64'

coefficients <- train_model_65$coefficients
se <- summary(train_model_65)$coefficients[, 2]
names <- names(train_model_65$coefficients)
model_65 <- data.frame(names, coefficients, se)
model_65$age_group <- '65+'

# Combine data
models <- rbind(model_0_4, model_5_12, model_13_17, model_18_49, 
                model_50_64, model_65)

# Calculate mean, lower, and upper bounds
models$lower <- exp(models$coefficients - 1.96*models$se)
models$upper <- exp(models$coefficients + 1.96*models$se)
models$mean <- exp(models$coefficients)

# Clean up coefficient names
models <- models |> dplyr::filter(names != '(Intercept)')
models$names <- str_to_title(models$names)
models$names <- ifelse(models$names == 'Sore_throat', 
                       'Sore Throat', models$names)
models$names <- ifelse(models$names == 'Nausea_vom', 
                       'Nausea', models$names)
models$names <- ifelse(models$names == 'Chest_pain', 
                       'Chest Pain', models$names)
models$names <- ifelse(models$names == 'Short_breath', 
                       'Short of Breath', models$names)

## CREATE FIGUE ##
# Set colors
colors <- c('0-4' = '#b9d4ec', '5-12' = '#8bb7e0', '13-17' = '#5d9ad3', 
            '18-49' =  '#347DC1', '50-64' = '#285f93', '65+' = '#1b4164')

# Set shape
shapes <- c(' ' = 5)

# Change variable names for plotting
models$`Age Group` <- models$age_group
models$`Age Group` <- factor(models$age_group, levels=c('0-4', '5-12', '13-17', '18-49', '50-64', '65+'))
model_all_flu_overall$`Overall` <- ' '

# Merge on order variable
models <- left_join(models, symptom_order_flu, by = c('names'))

# Set theme
forrest_theme <- theme(legend.position =  "inside",
                       legend.position.inside = c(0.873, 0.237),
                       axis.text = element_text(size=14, color = 'black'),
                       axis.title = element_text(size=16),
                       axis.title.y = element_blank(),
                       strip.background = element_blank(),
                       legend.title = element_text(size=16),
                       legend.text = element_text(size=14),
                       plot.title = element_text(size=18, hjust = 0),
                       legend.box="vertical",
                       legend.box.background = element_rect(colour = "black", 
                                                            fill = 'white'))


# Plot
age_flu_fig <- ggplot(data=models, aes(x=reorder(names, num), y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(aes(group = `Age Group`, color = `Age Group`, fill = `Age Group`),
                  position = position_dodge2(width = 0.5), size = 1, alpha = 0.7) + 
  geom_point(data=model_all_flu_overall[1:14,], aes(x=reorder(names, num), y=mean, shape = `Overall`),
             color = 'black', alpha = 2, size = 4, stroke = 1.2) +
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Odds Ratio") + scale_y_continuous(breaks = seq(0, 8, by = 2), limits = c(-0.01, 9.1)) +
  ggtitle("Influenza Odds by Symptom and Age Group") +
  theme_minimal() +  # use a white background
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(nrow = 6), shape = guide_legend(order = 1)) +
  forrest_theme


# Display
age_flu_fig







# Create figure, with labels
figure_2 <- cowplot::plot_grid(percent_flu, test_plot,
                               nrow = 2,
                               labels = c("a", "c"),
                               label_size = 20)

figure_2_supp <- cowplot::plot_grid(figure_2, age_flu_fig,
                               nrow = 1,
                               labels = c("", "b"),
                               rel_widths = c(1, 0.8),
                               label_size = 20)



# Save figure 
ggsave('./figs/figure_2.jpg', plot = figure_2_supp, height = 8, width = 20)


ggsave('./figs/figure_2_supp.jpg', plot = figure_2_supp, height = 7, width = 20)



################################
# D. CREATE TIME FORREST PLOTS #
################################

#############
# INFLUENZA #
#############

# Define models by geography
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
                                 sin(2*pi * epi_week_adj / 52) +
                                 cos(2*pi * epi_week_adj / 52),
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
                                 sin(2*pi * epi_week_adj / 52) +
                                 cos(2*pi * epi_week_adj / 52),
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
                                 sin(2*pi * epi_week_adj / 52) +
                                 cos(2*pi * epi_week_adj / 52),
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
                                 sin(2*pi * epi_week_adj / 52) +
                                 cos(2*pi * epi_week_adj / 52),
                               family = binomial(link = "logit"),
                               data = flu_train[flu_train$south == 0 & flu_train$urban_code_binary == 0,])

# Create data frames for each model
coefficients <- train_model_north_rural$coefficients
se <- summary(train_model_north_rural)$coefficients[, 2]
names <- names(train_model_north_rural$coefficients)
model_north_rural <- data.frame(names, coefficients, se)
model_north_rural$`Geography` <- 'North (Rural)'
model_north_rural <- model_north_rural[16:17,]

coefficients <- train_model_south_rural$coefficients
se <- summary(train_model_south_rural)$coefficients[, 2]
names <- names(train_model_south_rural$coefficients)
model_south_rural <- data.frame(names, coefficients, se)
model_south_rural$`Geography` <- 'South (Rural)'
model_south_rural <- model_south_rural[16:17,]

coefficients <- train_model_north_urban$coefficients
se <- summary(train_model_north_urban)$coefficients[, 2]
names <- names(train_model_north_urban$coefficients)
model_north_urban <- data.frame(names, coefficients, se)
model_north_urban$`Geography` <- 'North (Urban)'
model_north_urban <- model_north_urban[16:17,]

coefficients <- train_model_south_urban$coefficients
se <- summary(train_model_south_urban)$coefficients[, 2]
names <- names(train_model_south_urban$coefficients)
model_south_urban <- data.frame(names, coefficients, se)
model_south_urban$`Geography` <- 'South (Urban)'
model_south_urban <- model_south_urban[16:17,]

# Combine data
models <- rbind(model_north_urban, model_south_urban, model_south_rural, model_north_rural)

# Calculate mean, lower, and upper bounds
models$lower <- models$coefficients - 1.96*models$se
models$upper <- models$coefficients + 1.96*models$se
models$mean <- models$coefficients

model_dat_n_u <- data.frame(time = seq(1, 52, 1),
                            mean_odds = ((models[1, 2]) * sin(2 * pi * seq(1, 52)/52)) + 
                              ((models[2, 2]) * cos(2 * pi * seq(1, 52)/52)),
                            mean_upper = ((models[1, 6]) * sin(2 * pi * seq(1, 52)/52)) + 
                              ((models[2, 6]) * cos(2 * pi * seq(1, 52)/52)),
                            mean_lower = ((models[1, 5]) * sin(2 * pi * seq(1, 52)/52)) + 
                              ((models[2, 5]) * cos(2 * pi * seq(1, 52)/52)),
                            Geography = 'North (Urban)')

model_dat_s_u <- data.frame(time = seq(1, 52, 1),
                            mean_odds = ((models[3, 2]) * sin(2 * pi * seq(1, 52)/52)) + 
                              ((models[4, 2]) * cos(2 * pi * seq(1, 52)/52)),
                            mean_upper = ((models[3, 6]) * sin(2 * pi * seq(1, 52)/52)) + 
                              ((models[4, 6]) * cos(2 * pi * seq(1, 52)/52)),
                            mean_lower = ((models[3, 5]) * sin(2 * pi * seq(1, 52)/52)) + 
                              ((models[4, 5]) * cos(2 * pi * seq(1, 52)/52)),
                            Geography = 'South (Urban)')

model_dat_s_r <- data.frame(time = seq(1, 52, 1),
                            mean_odds = ((models[5, 2]) * sin(2 * pi * seq(1, 52)/52)) + 
                              ((models[6, 2]) * cos(2 * pi * seq(1, 52)/52)),
                            mean_upper = ((models[5, 6]) * sin(2 * pi * seq(1, 52)/52)) + 
                              ((models[6, 6]) * cos(2 * pi * seq(1, 52)/52)),
                            mean_lower = ((models[5, 5]) * sin(2 * pi * seq(1, 52)/52)) + 
                              ((models[6, 5]) * cos(2 * pi * seq(1, 52)/52)),
                            Geography = 'South (Rural)')

model_dat_n_r <- data.frame(time = seq(1, 52, 1),
                            mean_odds = ((models[7, 2]) * sin(2 * pi * seq(1, 52)/52)) + 
                              ((models[8, 2]) * cos(2 * pi * seq(1, 52)/52)),
                            mean_upper = ((models[7, 6]) * sin(2 * pi * seq(1, 52)/52)) + 
                              ((models[8, 6]) * cos(2 * pi * seq(1, 52)/52)),
                            mean_lower = ((models[7, 5]) * sin(2 * pi * seq(1, 52)/52)) + 
                              ((models[8, 5]) * cos(2 * pi * seq(1, 52)/52)),
                            Geography = 'North (Rural)')


model_dat <- rbind(model_dat_n_u, model_dat_s_u, model_dat_s_r, model_dat_n_r)


# Get overall model results

# Calculate overall results
coefficients <- time_model$coefficients
se <- summary(time_model)$coefficients[, 2]
names <- names(time_model$coefficients)
time_order <- data.frame(names, coefficients, se)
time_order$mean <- time_order$coefficients
time_order$lower <- time_order$coefficients - 1.96*time_order$se
time_order$upper <- time_order$coefficients + 1.96*time_order$se
time_order <- time_order[16:17,]

model_dat_overall <- data.frame(time = seq(1, 52, 1),
                            mean_odds = ((time_order[1, 2]) * sin(2 * pi * seq(1, 52)/52)) + 
                              ((models[2, 2]) * cos(2 * pi * seq(1, 52)/52)),
                            mean_upper = ((time_order[1, 6]) * sin(2 * pi * seq(1, 52)/52)) + 
                              ((models[2, 6]) * cos(2 * pi * seq(1, 52)/52)),
                            mean_lower = ((time_order[1, 5]) * sin(2 * pi * seq(1, 52)/52)) + 
                              ((models[2, 5]) * cos(2 * pi * seq(1, 52)/52)),
                            Overall = ' ')


# Figure out date alignment
flu_time <- flu |> ungroup() |> dplyr::select(c(year_week_dt, epi_week_adj)) |>
  dplyr::filter(year_week_dt > as.Date('2017-06-30') & 
                  year_week_dt < as.Date('2018-07-01')) |>
  distinct(year_week_dt, epi_week_adj)
  
model_dat_overall <- left_join(model_dat_overall, flu_time,
                               by = c('time' = 'epi_week_adj'))

model_dat <- left_join(model_dat, flu_time,
                       by = c('time' = 'epi_week_adj'))

## CREATE FIGUE ##
# Set shape
shapes <- c(' ' = 5)

# Set colors
colors <- c('South (Urban)' = '#45B649', 'North (Urban)' = '#9B59B6', 
            'South (Rural)' = '#b3e2b5', 'North (Rural)' = '#e0cbe8')

# Set theme
forrest_theme <- theme(legend.position =  "right",
                       legend.position.inside = c(0.70, 0.20),
                       axis.text = element_text(size=14),
                       axis.title = element_text(size=16),
                       strip.background = element_blank(),
                       legend.title = element_text(size=14),
                       legend.text = element_text(size=12),
                       plot.title = element_text(size=18, hjust = 0),
                       legend.box="horizontal",
                       legend.box.background = element_rect(colour = "white", 
                                                            fill = 'white'))

# Plot
geo_flu_fig <- ggplot(data=model_dat, aes(x=as.Date(year_week_dt), y=exp(mean_odds), ymin=exp(mean_lower), 
                                          ymax=exp(mean_upper))) +
  geom_line(data = model_dat, aes(group = Geography, color = Geography), linewidth = 1.2, alpha = 0.8) +
  geom_ribbon(aes(x = as.Date(year_week_dt), ymin=exp(mean_lower), ymax=exp(mean_upper), fill = Geography, 
                  group = Geography), alpha = 0.2) +
  geom_point(data = model_dat_overall, aes(x=as.Date(year_week_dt), y=exp(mean_odds), shape = `Overall`), alpha = 0.5, size = 1, stroke = 1) +
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors) +
  scale_shape_manual(values = shapes) +
  xlab("Time") + ylab("Odds Ratio") +
  ggtitle("Influenza Odds Seasonality by Geography") +
  theme_minimal() +  # use a white background
  guides(color = guide_legend(nrow = 4), fill = guide_legend(nrow = 4), shape = guide_legend(order = 1, override.aes = list(size = 3, stroke = 1.2, alpha = 1))) +
  forrest_theme +
  theme(plot.title = element_text(size=18, hjust = 0)) +
  theme(legend.position = 'right', legend.box = "vertical" )  +
  scale_x_date(date_breaks = "2 month", labels = date_format("%b"), limits = c(as.Date('2017-07-01'), as.Date('2018-06-25'))) +
  scale_y_continuous(breaks = seq(0, 8, by = 2)) 

geo_flu_fig



# Display
geo_flu_fig

#######
# RSV #
#######

# Run models trained on each age group
train_model_0_2_rsv <- glm(rsv ~ fever + cough + sore_throat +
                             short_breath + hypoxemia + 
                             bronchitis + loss_appetite + 
                             fatigue + headache + congestion + sneezing,
                           family = binomial(link = "logit"),
                           data = rsv_train[rsv_train$age_grp == '0',])

train_model_3_4_rsv <- glm(rsv ~ fever + cough + sore_throat +
                             short_breath + hypoxemia + 
                             bronchitis + loss_appetite + 
                             fatigue + headache + congestion + sneezing,
                           family = binomial(link = "logit"),
                           data = rsv_train[rsv_train$age_grp == '1',])

train_model_5_17_rsv <- glm(rsv ~ fever + cough + sore_throat +
                              short_breath + hypoxemia + 
                              bronchitis + loss_appetite + 
                              fatigue + headache + congestion + sneezing,
                            family = binomial(link = "logit"),
                            data = rsv_train[rsv_train$age_grp == '2',])

train_model_18_64_rsv <- glm(rsv ~ fever + cough + sore_throat +
                               short_breath + hypoxemia + 
                               bronchitis + loss_appetite + 
                               fatigue + headache + congestion + sneezing,
                             family = binomial(link = "logit"),
                             data = rsv_train[rsv_train$age_grp == '3',])

train_model_65_rsv <- glm(rsv ~ fever + cough + sore_throat +
                            short_breath + hypoxemia + 
                            bronchitis + loss_appetite + 
                            fatigue + headache + congestion + sneezing,
                          family = binomial(link = "logit"),
                          data = rsv_train[rsv_train$age_grp == '4',])

# Create data for the overall model
coefficients <-symp_model_rsv$coefficients
se <- summary(symp_model_rsv)$coefficients[, 2]
names <- names(symp_model_rsv$coefficients)
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

# Create a data frame for each age model
coefficients <- train_model_0_2_rsv$coefficients
se <- summary(train_model_0_2_rsv)$coefficients[, 2]
names <- names(train_model_0_2_rsv$coefficients)
model_0_2 <- data.frame(names, coefficients, se)
model_0_2$age_group <- '0-2'

coefficients <- train_model_3_4_rsv$coefficients
se <- summary(train_model_3_4_rsv)$coefficients[, 2]
names <- names(train_model_3_4_rsv$coefficients)
model_3_4 <- data.frame(names, coefficients, se)
model_3_4$age_group <- '3-4'

coefficients <- train_model_5_17_rsv$coefficients
se <- summary(train_model_5_17_rsv)$coefficients[, 2]
names <- names(train_model_5_17_rsv$coefficients)
model_5_17 <- data.frame(names, coefficients, se)
model_5_17$age_group <- '5-17'

coefficients <- train_model_18_64_rsv$coefficients
se <- summary(train_model_18_64_rsv)$coefficients[, 2]
names <- names(train_model_18_64_rsv$coefficients)
model_18_64 <- data.frame(names, coefficients, se)
model_18_64$age_group <- '18-64'

coefficients <- train_model_65_rsv$coefficients
se <- summary(train_model_65_rsv)$coefficients[, 2]
names <- names(train_model_65_rsv$coefficients)
model_65 <- data.frame(names, coefficients, se)
model_65$age_group <- '65+  '

# Combine data
models_rsv <- rbind(model_0_2, model_3_4, model_5_17, model_18_64, 
                    model_65)

# Calculate mean, lower, and upper bounds
models_rsv$lower <- exp(models_rsv$coefficients - 1.96*models_rsv$se)
models_rsv$upper <- exp(models_rsv$coefficients + 1.96*models_rsv$se)
models_rsv$mean <- exp(models_rsv$coefficients)

# Clean up coefficient names
models_rsv <- models_rsv |> filter(names != '(Intercept)')
models_rsv$names <- str_to_title(models_rsv$names)
models_rsv$names <- ifelse(models_rsv$names == 'Sore_throat', 
                       'Sore Throat', models_rsv$names)
models_rsv$names <- ifelse(models_rsv$names == 'Nausea_vom', 
                       'Nausea', models_rsv$names)
models_rsv$names <- ifelse(models_rsv$names == 'Chest_pain', 
                       'Chest Pain', models_rsv$names)
models_rsv$names <- ifelse(models_rsv$names == 'Short_breath', 
                       'Shortness of Breath', models_rsv$names)
models_rsv$names <- ifelse(models_rsv$names == 'Loss_appetite', 
                       'Loss of Appetite', models_rsv$names)

## CREATE FIGUE ##
# Set colors
colors <- c('0-2' = '#9edbd8', '3-4' = '#73cbc7', '5-17' = '#41afaa', 
            '18-64' = '#36928e', '65+  ' = '#2c7672'
)

# Set shape
shapes <- c(' ' = 5)

# Change variable names for plotting
model_all_rsv_overall$`Overall` <- ' '
models_rsv$`Age Group` <- models_rsv$age_group
models_rsv$`Age Group` <- factor(models_rsv$age_group, levels=c('0-2', '3-4', '5-17', '18-64', '65+  '))

# Merge on order variable
models_rsv <- left_join(models_rsv, symptom_order_rsv, by = c('names'))

# Plot
age_rsv_fig <- ggplot(data=models_rsv, aes(x=reorder(names, num), y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(aes(group = `Age Group`, color = `Age Group`, fill = `Age Group`),
                  position = position_dodge2(width = 0.4), size = 0.75, alpha = 0.7) + 
  geom_point(data=model_all_rsv_overall, aes(x=reorder(names, num), y=mean, shape = `Overall`),
             color = 'black', alpha = 2, size = 3, stroke = 1.2) +
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("Coefficients") + ylab("Odds Ratio") +
  ggtitle("RSV Regression Results by Age") +
  theme_minimal() +  # use a white background
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = colors) +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 1, order = 1)) +
  forrest_theme + theme(plot.title = element_text(size=18, hjust = 0)) + theme(legend.box.just = "left")

# Display 
age_rsv_fig

###############################
# CREATE MAIN EFFECT AGE PLOT #
###############################

# Create a data frame
coefficients <- age_model_rsv$coefficients
se <- summary(age_model_rsv)$coefficients[, 2]
names <- names(age_model_rsv$coefficients)
model_age <- data.frame(names, coefficients, se)
model_age <- model_age[3:6,]
model_age$names <- substr(model_age$names, 21, nchar(model_age$names))
model_age$num <- rank(model_age$coefficients)

# Calculate mean, lower, and upper bounds
model_age$lower <- exp(model_age$coefficients - 1.96*model_age$se)
model_age$upper <- exp(model_age$coefficients + 1.96*model_age$se)
model_age$mean <- exp(model_age$coefficients)

# Reset order
model_age[model_age$names == '65+',]$names <- '65+  '
model_age$num <- c(1, 2, 3, 4)

# Plot
age_overall_rsv_fig <- ggplot(data=model_age, aes(x=reorder(names, num), y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(aes(color = names), size = 0.55, alpha = 0.7) + 
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  xlab("Age Group") + ylab("Odds Ratio") + ylim(0, 1) + scale_y_continuous(breaks = seq(0, 1, by = 0.5)) +
  ggtitle("RSV Odds Relative to Age 0-2") +
  scale_color_manual(values = colors) + 
  theme_minimal() +  # use a white background
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(order = 1)) +
  forrest_other_theme 

# Display
age_overall_rsv_fig

# Combine figures
age_rsv_fig |> 
  ggdraw() +
  draw_plot(
    {
      age_overall_rsv_fig
    },
    x = 0.59, 
    y = 0.266,
    width = 0.4, 
    height = 0.3)



#######
# RSV #
#######

# Define models by geography
train_model_south_urban_rsv <- glm(rsv ~ fever + cough + sore_throat +
                                     short_breath + hypoxemia + 
                                     bronchitis + loss_appetite + 
                                     fatigue + headache + congestion + sneezing +
                                     as.factor(month),
                                   family = binomial(link = "logit"),
                                   data = rsv_train[rsv_train$south == 1 & rsv_train$urban_code_binary == 1,])

train_model_south_rural_rsv <- glm(rsv ~ fever + cough + sore_throat +
                                     short_breath + hypoxemia + 
                                     bronchitis + loss_appetite + 
                                     fatigue + headache + congestion + sneezing +
                                     as.factor(month),
                                   family = binomial(link = "logit"),
                                   data = rsv_train[rsv_train$south == 1 & rsv_train$urban_code_binary == 0,])

train_model_north_urban_rsv <- glm(rsv ~ fever + cough + sore_throat +
                                     short_breath + hypoxemia + 
                                     bronchitis + loss_appetite + 
                                     fatigue + headache + congestion + sneezing +
                                     as.factor(month),
                                   family = binomial(link = "logit"),
                                   data = rsv_train[rsv_train$south == 0 & rsv_train$urban_code_binary == 1,])

train_model_north_rural_rsv <- glm(rsv ~ fever + cough + sore_throat +
                                     short_breath + hypoxemia + 
                                     bronchitis + loss_appetite + 
                                     fatigue + headache + congestion + sneezing +
                                     as.factor(month),
                                   family = binomial(link = "logit"),
                                   data = rsv_train[rsv_train$south == 0 & rsv_train$urban_code_binary == 0,])

# Create data frames for each model
coefficients <- train_model_north_rural_rsv$coefficients
se <- summary(train_model_north_rural_rsv)$coefficients[, 2]
names <- names(train_model_north_rural_rsv$coefficients)
model_north_rural <- data.frame(names, coefficients, se)
model_north_rural$`Geography` <- 'North (Rural)'
model_north_rural <- model_north_rural[13:23,]

coefficients <- train_model_south_rural_rsv$coefficients
se <- summary(train_model_south_rural_rsv)$coefficients[, 2]
names <- names(train_model_south_rural_rsv$coefficients)
model_south_rural <- data.frame(names, coefficients, se)
model_south_rural$`Geography` <- 'South (Rural)'
model_south_rural <- model_south_rural[13:23,]

coefficients <- train_model_north_urban_rsv$coefficients
se <- summary(train_model_north_urban_rsv)$coefficients[, 2]
names <- names(train_model_north_urban_rsv$coefficients)
model_north_urban <- data.frame(names, coefficients, se)
model_north_urban$`Geography` <- 'North (Urban)'
model_north_urban <- model_north_urban[13:23,]

coefficients <- train_model_south_urban_rsv$coefficients
se <- summary(train_model_south_urban_rsv)$coefficients[, 2]
names <- names(train_model_south_urban_rsv$coefficients)
model_south_urban <- data.frame(names, coefficients, se)
model_south_urban$`Geography` <- 'South (Urban)'
model_south_urban <- model_south_urban[13:23,]

# Combine data
models_rsv <- rbind(model_north_urban, model_south_urban, model_south_rural, model_north_rural)

# Calculate mean, lower, and upper bounds
models_rsv$lower <- exp(models_rsv$coefficients - 1.96*models_rsv$se)
models_rsv$upper <- exp(models_rsv$coefficients + 1.96*models_rsv$se)
models_rsv$mean <- exp(models_rsv$coefficients)

# Clean up coefficient names 
models_rsv <- models_rsv |> filter(names != '(Intercept)')
models_rsv$names <- str_to_title(models_rsv$names)
models_rsv$names <- ifelse(models_rsv$names == 'Sore_throat', 
                       'Sore Throat', models_rsv$names)
models_rsv$names <- ifelse(models_rsv$names == 'Nausea_vom', 
                       'Nausea', models_rsv$names)
models_rsv$names <- ifelse(models_rsv$names == 'Chest_pain', 
                       'Chest Pain', models_rsv$names)
models_rsv$names <- ifelse(models_rsv$names == 'Short_breath', 
                       'Shortness of Breath', models_rsv$names)
models_rsv$names <- gsub("As.factor\\(Month\\)", "", models_rsv$names)

# Calculate overall results
coefficients <- time_model_rsv$coefficients
se <- summary(time_model_rsv)$coefficients[, 2]
names <- names(time_model_rsv$coefficients)
time_order <- data.frame(names, coefficients, se)
time_order$mean <- exp(time_order$coefficients)
time_order$lower <- exp(time_order$coefficients - 1.96*time_order$se)
time_order$upper <- exp(time_order$coefficients + 1.96*time_order$se)
time_order$names <- gsub("as.factor\\(month\\)", "", time_order$names)
time_order$`Overall` <- ' '
time_order <- time_order[17:27,]

## CREATE FIGUE ##
# Set shape
shapes <- c(' ' = 5)

# Set colors
colors <- c('South (Urban)' = '#509048', 'North (Urban)' = '#896deb', 
            'South (Rural)' = '#90c058', 'North (Rural)' = '#c2a3fd')

# Change month number to name
months <- c("Feb","Mar",
            "Apr","May","Jun",
            "Jul","Aug","Sep",
            "Oct","Nov","Dec")
month_num <- as.character(c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
month <- data.frame(months, month_num)

# Merge on month name
time_order <- left_join(time_order, month, by = c('names' = 'month_num'))
models_rsv <- left_join(models_rsv, month, by = c('names' = 'month_num'))

# Plot
geo_rsv_fig <- ggplot(data=models_rsv, aes(x=reorder(months, as.numeric(names)), y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(data = models_rsv, aes(group = `Geography`, color = `Geography`, fill = `Geography`),
                  position = position_dodge2(width = 0.4), size = 0.75, alpha = 0.7) + 
  geom_point(data=time_order, aes(x=reorder(months, as.numeric(names)), y=mean, shape = `Overall`),
             color = 'black', alpha = 2, size = 3, stroke = 1.2) +
  geom_hline(yintercept=1, lty=2) +  # add a dotted line at x=1 after flip
  xlab("Month") + ylab("Odds Ratio") +
  ggtitle("RSV Odds Relative to January by Geography") +
  theme_minimal() +  # use a white background
  scale_shape_manual(values = shapes) +
  scale_color_manual(values = colors) +
  guides(color = guide_legend(nrow = 2), shape = guide_legend(order = 1)) +
  forrest_theme +
  theme(plot.title = element_text(size=18, hjust = 0)) +
  theme(legend.position = 'right') + theme(legend.box.just = "left")

# Display
geo_rsv_fig


age_flu_fig |> 
  ggdraw() +
  draw_plot(
    {
      age_overall_flu_fig
    },
    x = 0.59, 
    y = 0.266,
    width = 0.4, 
    height = 0.3)

age_rsv_fig |> 
  ggdraw() +
  draw_plot(
    {
      age_overall_rsv_fig
    },
    x = 0.59, 
    y = 0.266,
    width = 0.4, 
    height = 0.3)


sub_1 <- cowplot::plot_grid(gender_flu_fig, age_flu_fig,
                            nrow = 1,
                            labels = c("b", "c"),
                            label_size = 20)


# Create figure, with labels
figure_2 <- cowplot::plot_grid(percent_flu, sub_1, geo_flu_fig,
                               nrow = 3,
                               labels = c("a", "", "d"),
                               label_size = 20)

figure_2_supp <- cowplot::plot_grid(percent_rsv, age_rsv_fig |> 
                                      ggdraw() +
                                      draw_plot(
                                        {
                                          age_overall_rsv_fig
                                        },
                                        x = 0.529, 
                                        y = 0.244,
                                        width = 0.46, 
                                        height = 0.313), geo_rsv_fig,
                                    nrow = 1,
                                    labels = c("a", "b", "c"),
                                    label_size = 20)

# Save figure 
ggsave('./figs/figure_2.jpg', plot = figure_2, height = 14, width = 14)
ggsave('./figs/figure_2_supp.jpg', plot = figure_2_supp, height = 7, width = 20)

