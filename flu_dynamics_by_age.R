# Start with a clear environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(lubridate)
library(car)
library(plotROC)
library(pROC)
library(cutpointr)
library(zoo)
library(ggcorrplot)

# Set seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

#########################
# LOAD AND PREPARE DATA #
#########################

# Load processed data
# Part 1 full data
flu_p1 <- readRDS('./tmp/flu_p1_dat_proc.rds')

# Part 2 flu data (major counties)
flu_p2 <- readRDS('./tmp/flu_p2_dat_proc_select_counties.rds')

# Add time variables to Part 1 data
flu_p1$month <- month(flu_p1$month_date)
flu_p1$year <- year(flu_p1$month_date)
flu_p1$jan <- ifelse(flu_p1$month == 1, 1, 0)
flu_p1$feb <- ifelse(flu_p1$month == 2, 1, 0)
flu_p1$mar <- ifelse(flu_p1$month == 3, 1, 0)
flu_p1$apr <- ifelse(flu_p1$month == 4, 1, 0)
flu_p1$may <- ifelse(flu_p1$month == 5, 1, 0)
flu_p1$jun <- ifelse(flu_p1$month == 6, 1, 0)
flu_p1$jul <- ifelse(flu_p1$month == 7, 1, 0)
flu_p1$aug <- ifelse(flu_p1$month == 8, 1, 0)
flu_p1$sep <- ifelse(flu_p1$month == 9, 1, 0)
flu_p1$oct <- ifelse(flu_p1$month == 10, 1, 0)
flu_p1$nov <- ifelse(flu_p1$month == 11, 1, 0)
flu_p1$dec <- ifelse(flu_p1$month == 12, 1, 0)

# Add time variables to Part 2 data
flu_p2$month <- month(flu_p2$week_date)
flu_p2$year <- year(flu_p2$week_date)
flu_p2$month_date <- as.Date(paste(flu_p2$month, '01', flu_p2$year, sep = '-'), "%m-%d-%Y")
flu_p2$jan <- ifelse(flu_p2$month == 1, 1, 0)
flu_p2$feb <- ifelse(flu_p2$month == 2, 1, 0)
flu_p2$mar <- ifelse(flu_p2$month == 3, 1, 0)
flu_p2$apr <- ifelse(flu_p2$month == 4, 1, 0)
flu_p2$may <- ifelse(flu_p2$month == 5, 1, 0)
flu_p2$jun <- ifelse(flu_p2$month == 6, 1, 0)
flu_p2$jul <- ifelse(flu_p2$month == 7, 1, 0)
flu_p2$aug <- ifelse(flu_p2$month == 8, 1, 0)
flu_p2$sep <- ifelse(flu_p2$month == 9, 1, 0)
flu_p2$oct <- ifelse(flu_p2$month == 10, 1, 0)
flu_p2$nov <- ifelse(flu_p2$month == 11, 1, 0)
flu_p2$dec <- ifelse(flu_p2$month == 12, 1, 0)

# Create age group variable
# Part 1
flu_p1$age_group <- factor(flu_p1$age_grp,
                           levels = c('0', '1', '2', '3', '4', '5'),
                           labels = c('0-4', '5-12', '13-17', '18-49', '50-64', '65+'))

# Part 2
flu_p2$age_group <- factor(flu_p2$age_grp,
                           levels = c('0', '1', '2', '3', '4', '5'),
                           labels = c('0-4', '5-12', '13-17', '18-49', '50-64', '65+'))

# Restrict data to symptomatic cases only
# Part 1
flu_symp_p1 <- flu_p1 |> dplyr::filter(symp_count > 0)
#remove(flu_p1)
# Part 2
flu_symp_p2 <- flu_p2 |> dplyr::filter(symp_count > 0)
remove(flu_p2)

################
# EXPLORE DATA #
################

# Create some exploratory plots
# Flu dynamics by age
flu_week <- flu_symp_p2 |> 
  group_by(week_date, age_group) |> 
  mutate(flu_sum = sum(flu*patient_count_imp)) |>
  distinct(week_date, age_group, flu_sum, .keep_all = FALSE)

# Select symptom dynamics by age
symptom_week <- flu_symp_p2 |> 
  group_by(week_date, age_group) |> 
  mutate(fever_sum = sum(fever*patient_count_imp),
         cough_sum = sum(cough*patient_count_imp),
         sore_throat_sum = sum(sore_throat*patient_count_imp)) |>
  distinct(week_date, age_group, fever_sum, cough_sum, sore_throat_sum, 
           .keep_all = FALSE)

# Plot flu dynamics by age
ggplot() + 
  geom_line(data = flu_week, 
            aes(x = week_date, y = flu_sum), 
            color = '#C7A939', linewidth = 0.8) + 
  ylab('Flu Count') + xlab('Week') + 
  ggtitle('Flu Dynamics') +
  theme_minimal() + facet_wrap(~age_group)

# Plot symptom dynamics by age
colors <- c("Fever" = "#c36272", "Cough" = "#5a9374", "Sore Throat" = "#347dab")
ggplot() + 
  geom_line(data = symptom_week, 
            aes(x = week_date, y = sore_throat_sum, color = 'Sore Throat'), 
            linewidth = 0.8, alpha = 0.8) + 
  geom_line(data = symptom_week, 
            aes(x = week_date, y = cough_sum, color = 'Cough'), 
            linewidth = 0.8, alpha = 0.8) +
  geom_line(data = symptom_week, 
            aes(x = week_date, y = fever_sum, color = 'Fever'), 
            linewidth = 0.9, alpha = 0.8) +
  ggtitle('Symptom Dynamics') +
  theme_minimal() +
  labs(x = "Week", y = "Symptom Count", color = "Legend") +
  scale_color_manual(values = colors) + theme(legend.position = 'bottom') +
  facet_wrap(~age_group)

# Create correlation matrices
# Reshape long
flu_week_wide <-  spread(flu_week, 
                         key=age_group, 
                         value=flu_sum) 

fever_week_wide <-  spread(symptom_week[, c(1, 2, 3)], 
                         key=age_group, 
                         value=fever_sum) 

cough_week_wide <-  spread(symptom_week[, c(1, 2, 4)], 
                         key=age_group, 
                         value=cough_sum) 

st_week_wide <-  spread(symptom_week[, c(1, 2, 5)], 
                         key=age_group, 
                         value=sore_throat_sum) 

# Flu correlations
corr_flu <- cor(flu_week_wide[, -c(1)])
ggcorrplot(corr_flu, lab = TRUE,
           ggtheme = ggplot2::theme_minimal(),
           colors = c("#c36272", "white", "#347dab")) +
  ggtitle('Flu Age Correlations')

# Fever correlations
corr_fever <- cor(fever_week_wide[, -c(1)])
ggcorrplot(corr_fever, lab = TRUE,
           ggtheme = ggplot2::theme_minimal(),
           colors = c("#c36272", "white", "#347dab")) +
  ggtitle('Fever Age Correlations')

# Cough correlations
corr_cough <- cor(cough_week_wide[, -c(1)])
ggcorrplot(corr_cough, lab = TRUE,
           ggtheme = ggplot2::theme_minimal(),
           colors = c("#c36272", "white", "#347dab")) +
  ggtitle('Cough Age Correlations')

# Sore Throat correlations
st_flu <- cor(st_week_wide[, -c(1)])
ggcorrplot(st_flu, lab = TRUE,
           ggtheme = ggplot2::theme_minimal(),
           colors = c("#c36272", "white", "#347dab")) +
  ggtitle('Sore Throat Age Correlations')

# Percentage of flu cases with symptom by age group
flu_descriptives_symp <- flu_p1 |> mutate(
  fever = ifelse(fever == 1, patient_count_imp, 0),
  myalgia = ifelse( myalgia == 1, patient_count_imp, 0),
  short_breath = ifelse(short_breath == 1, patient_count_imp, 0),
  hypoxemia = ifelse(hypoxemia == 1, patient_count_imp, 0),
  chest_pain = ifelse(chest_pain == 1, patient_count_imp, 0),
  cough = ifelse(cough == 1, patient_count_imp, 0),
  bronchitis = ifelse(bronchitis == 1, patient_count_imp, 0),
  nausea_vom = ifelse(nausea_vom == 1, patient_count_imp, 0),
  sore_throat = ifelse(sore_throat == 1, patient_count_imp, 0),
  diarrhea = ifelse(diarrhea == 1, patient_count_imp, 0),
  fatigue = ifelse(fatigue == 1, patient_count_imp, 0),
  headache = ifelse(headache == 1, patient_count_imp, 0),
  congestion = ifelse(congestion == 1, patient_count_imp, 0),
  sneezing = ifelse(sneezing == 1, patient_count_imp, 0)) |>
  group_by(flu, age_group) |> mutate(Fever = sum(fever),
                               Myalgia = sum(myalgia),
                               `Short Breath` = sum(short_breath),
                               Hypoxemia = sum(hypoxemia),
                               `Chest Pain` = sum(chest_pain),
                               Cough = sum(cough),
                               Bronchitis = sum(bronchitis),
                               `Nausea/Vomit` = sum(nausea_vom),
                               `Sore Throat` = sum(sore_throat),
                               Diarrhea = sum(diarrhea),
                               Fatigue = sum(fatigue),
                               Headache = sum(headache),
                               Congestion = sum(congestion),
                               Sneezing = sum(sneezing),
                               patient_sum = sum(patient_count_imp)) |>
  distinct(flu, age_group, Fever, Myalgia, `Short Breath`, Hypoxemia, `Chest Pain`,
           Cough, Bronchitis, `Nausea/Vomit`, `Sore Throat`, Diarrhea, Fatigue,
           Headache, Congestion, Sneezing, patient_sum, .keep_all = FALSE) |>
  pivot_longer(cols=c('Fever', 'Myalgia', 'Short Breath', 'Hypoxemia', 'Chest Pain',
                      'Cough', 'Bronchitis', 'Nausea/Vomit', 'Sore Throat', 'Diarrhea', 'Fatigue',
                      'Headache', 'Congestion', 'Sneezing'),
               names_to='Symptom',
               values_to='Count') |>
  mutate(Percent = Count/patient_sum)

ggplot(flu_descriptives_symp[flu_descriptives_symp$flu == 1,], aes(x = Symptom, y = Percent)) +
  geom_bar(aes(), fill = '#C7A939', stat = "identity") + ggtitle('Symptom Percent: Flu Cases') +
  theme_minimal() + xlab('') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_wrap(~age_group)

ggplot(flu_descriptives_symp[flu_descriptives_symp$flu == 0,], aes(x = Symptom, y = Percent)) +
  geom_bar(aes(), fill = '#9086ba', stat = "identity") + ggtitle('Symptom Percent: Non-Cases') +
  theme_minimal() + xlab('') + ylim(0, 0.4) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  facet_wrap(~age_group)

##############
# TRAIN DATA #
##############

# Split up counties into test and train
counties <- flu_symp_p1 |> group_by(county_fips) |>
  distinct(county_fips)
sample_counties <- sample(c(TRUE, FALSE), nrow(counties), replace=TRUE, prob=c(0.5, 0.5))
counties_train  <- counties[sample_counties, ]
counties_train$train <- 1

# Split the data by merging on train counties
flu_symp_p1 <- left_join(flu_symp_p1, counties_train, by = c('county_fips' = 'county_fips'))
flu_symp_p1_train <- flu_symp_p1 |> filter(train == 1)
flu_symp_p1_test <- flu_symp_p1 |> filter(is.na(train))

# Sample from the full Part 1 data
# Unbalanced data
unbalanced <- flu_symp_p1_train |> group_by(age_group) |> 
  uncount(patient_count_imp) |> sample_n(20000)

# Balanced Data
# cases <- flu_symp_p1_train |> filter(flu == 1) |> group_by(age_group) |>
#   uncount(patient_count_imp) |> sample_n(10000)
# non_cases <- flu_symp_p1_train |> filter(flu == 0) |> group_by(age_group) |>
#   uncount(patient_count_imp) |> sample_n(10000)
# balanced <- rbind(cases, non_cases)

# Train the model by age
train_model_0_4 <- glm(flu ~ fever + cough + sore_throat +
                     myalgia + short_breath + hypoxemia +
                     nausea_vom + bronchitis + chest_pain + 
                     diarrhea + fatigue + headache + 
                     congestion + sneezing + jan + feb + mar + 
                     apr + may + jun + jul + aug + sep + oct + nov + dec,
                   family = binomial(link = "logit"),
                   data = unbalanced[unbalanced$age_group == '0-4',])

train_model_5_12 <- glm(flu ~ fever + cough + sore_throat +
                         myalgia + short_breath + hypoxemia +
                         nausea_vom + bronchitis + chest_pain + 
                         diarrhea + fatigue + headache + 
                         congestion + sneezing + jan + feb + mar + 
                         apr + may + jun + jul + aug + sep + oct + nov + dec,
                       family = binomial(link = "logit"),
                       data = unbalanced[unbalanced$age_group == '5-12',])

train_model_13_17 <- glm(flu ~ fever + cough + sore_throat +
                         myalgia + short_breath + hypoxemia +
                         nausea_vom + bronchitis + chest_pain + 
                         diarrhea + fatigue + headache + 
                         congestion + sneezing + jan + feb + mar + 
                         apr + may + jun + jul + aug + sep + oct + nov + dec,
                       family = binomial(link = "logit"),
                       data = unbalanced[unbalanced$age_group == '13-17',])

train_model_18_49 <- glm(flu ~ fever + cough + sore_throat +
                         myalgia + short_breath + hypoxemia +
                         nausea_vom + bronchitis + chest_pain + 
                         diarrhea + fatigue + headache + 
                         congestion + sneezing + jan + feb + mar + 
                         apr + may + jun + jul + aug + sep + oct + nov + dec,
                       family = binomial(link = "logit"),
                       data = unbalanced[unbalanced$age_group == '18-49',])

train_model_50_64 <- glm(flu ~ fever + cough + sore_throat +
                         myalgia + short_breath + hypoxemia +
                         nausea_vom + bronchitis + chest_pain + 
                         diarrhea + fatigue + headache + 
                         congestion + sneezing + jan + feb + mar + 
                         apr + may + jun + jul + aug + sep + oct + nov + dec,
                       family = binomial(link = "logit"),
                       data = unbalanced[unbalanced$age_group == '50-64',])

train_model_65 <- glm(flu ~ fever + cough + sore_throat +
                           myalgia + short_breath + hypoxemia +
                           nausea_vom + bronchitis + chest_pain + 
                           diarrhea + fatigue + headache + 
                           congestion + sneezing + jan + feb + mar + 
                           apr + may + jun + jul + aug + sep + oct + nov + dec,
                         family = binomial(link = "logit"),
                         data = unbalanced[unbalanced$age_group == '65+',])

train_model <- glm(flu ~ fever*as.factor(age_grp) + cough*as.factor(age_grp) + 
                     sore_throat*as.factor(age_grp) +
                     myalgia*as.factor(age_grp) + short_breath*as.factor(age_grp) + 
                     hypoxemia*as.factor(age_grp) +
                     nausea_vom *as.factor(age_grp)+ bronchitis*as.factor(age_grp) + 
                     chest_pain*as.factor(age_grp) + 
                     diarrhea*as.factor(age_grp) + fatigue*as.factor(age_grp) + 
                     headache*as.factor(age_grp) + congestion*as.factor(age_grp) + 
                     sneezing*as.factor(age_grp) + jan + feb + mar + 
                     apr + may + jun + jul + aug + sep + oct + nov + dec,
                   family = binomial(link = "logit"),
                   data = unbalanced)

train_model_time <- glm(flu ~ jan + feb + mar + 
                     apr + may + jun + jul + aug + sep + oct + nov + dec,
                   family = binomial(link = "logit"),
                   data = unbalanced)

train_model_fever <- glm(flu ~ fever,
                        family = binomial(link = "logit"),
                        data = unbalanced)

train_model_symptom <- glm(flu ~ 0 + fever + cough + sore_throat +
                             myalgia + short_breath + hypoxemia +
                             nausea_vom + bronchitis + chest_pain + 
                             diarrhea + fatigue + headache + 
                             congestion + sneezing,
                         family = binomial(link = "logit"),
                         data = unbalanced)

train_model_fever_time <- glm(flu ~ fever + jan + feb + mar + 
                                apr + may + jun + jul + aug + sep + oct + nov + dec,
                         family = binomial(link = "logit"),
                         data = unbalanced)

train_model_symptom_time <- glm(flu ~ fever + cough + sore_throat +
                             + myalgia + short_breath + hypoxemia +
                             nausea_vom + bronchitis + chest_pain + 
                             diarrhea + fatigue + headache + 
                             congestion + sneezing + jan + feb + mar + 
                               apr + may + jun + jul + aug + sep + oct + nov + dec,
                           family = binomial(link = "logit"),
                           data = unbalanced)

train_model_symptom_age <- glm(flu ~ fever*as.factor(age_grp) + cough*as.factor(age_grp) + 
                                 sore_throat*as.factor(age_grp) +
                                 + myalgia*as.factor(age_grp) + short_breath*as.factor(age_grp) + 
                                 hypoxemia*as.factor(age_grp) +
                                 nausea_vom *as.factor(age_grp)+ bronchitis*as.factor(age_grp) + 
                                 chest_pain*as.factor(age_grp) + 
                                 diarrhea*as.factor(age_grp) + fatigue*as.factor(age_grp) + 
                                 headache*as.factor(age_grp) + congestion*as.factor(age_grp) + 
                                 sneezing*as.factor(age_grp),
                                family = binomial(link = "logit"),
                                data = unbalanced)

################################
# DISPLAY ESTIMATE DIFFERENCES #
################################

# Make data sets with coefficients and se
coefficients <- train_model_0_4$coefficients
se <- summary(train_model_0_4)$coefficients[, 2]
names <- names(train_model_0_4$coefficients)
model_0_4 <- data.frame(names, coefficients, se)
model_0_4$age_group <- '0-4'
model_0_4$num <- seq(1, 26, 1)
summary(train_model_0_4)
coefficients <- train_model_5_12$coefficients
se <- summary(train_model_5_12)$coefficients[, 2]
names <- names(train_model_5_12$coefficients)
model_5_12 <- data.frame(names, coefficients, se)
model_5_12$age_group <- '5-12'
model_5_12$num <- seq(1, 26, 1)

coefficients <- train_model_13_17$coefficients
se <- summary(train_model_13_17)$coefficients[, 2]
names <- names(train_model_13_17$coefficients)
model_13_17 <- data.frame(names, coefficients, se)
model_13_17$age_group <- '13-17'
model_13_17$num <- seq(1, 26, 1)

coefficients <- train_model_18_49$coefficients
se <- summary(train_model_18_49)$coefficients[, 2]
names <- names(train_model_18_49$coefficients)
model_18_49 <- data.frame(names, coefficients, se)
model_18_49$age_group <- '18-49'
model_18_49$num <- seq(1, 26, 1)

coefficients <- train_model_50_64$coefficients
se <- summary(train_model_50_64)$coefficients[, 2]
names <- names(train_model_50_64$coefficients)
model_50_64 <- data.frame(names, coefficients, se)
model_50_64$age_group <- '50-64'
model_50_64$num <- seq(1, 26, 1)

coefficients <- train_model_65$coefficients
se <- summary(train_model_65)$coefficients[, 2]
names <- names(train_model_65$coefficients)
model_65 <- data.frame(names, coefficients, se)
model_65$age_group <- '65+'
model_65$num <- seq(1, 26, 1)

# Combine data
models <- rbind(model_0_4, model_5_12, model_13_17, model_18_49, 
                model_50_64, model_65)

# Calculate lower and upper bounds
models$lower <- exp(models$coefficients - 1.96*models$se)
models$upper <- exp(models$coefficients + 1.96*models$se)
models$mean <- exp(models$coefficients)

# Create figure

colors <- c('0-4' = '#c36272', '5-12' = '#C7A939', '13-17' = '#5a9374', 
            '18-49' = '#347dab', '50-64' = '#9086ba', '65+' = 'darkgrey'
)

figure <- ggplot(data=models, aes(x=reorder(names, -num), y=mean, ymin=lower, ymax=upper)) +
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

ggsave('./figs/coeffs_by_age.jpg', plot = figure, height = 15, width = 8)

################
# PREDICT DATA #
################

# Predict testing data
# Expand the testing data
flu_symp_p1_test_expand <- flu_symp_p1_test |> group_by(age_group) |>
  uncount(patient_count_imp) |> sample_n(20000)

# Separate out the age groups
flu_symp_p1_test_0_4 <- 
  flu_symp_p1_test_expand[flu_symp_p1_test_expand$age_group == '0-4',]
flu_symp_p1_test_5_12 <- 
  flu_symp_p1_test_expand[flu_symp_p1_test_expand$age_group == '5-12',]
flu_symp_p1_test_13_17 <- 
  flu_symp_p1_test_expand[flu_symp_p1_test_expand$age_group == '13-17',]
flu_symp_p1_test_18_49 <- 
  flu_symp_p1_test_expand[flu_symp_p1_test_expand$age_group == '18-49',]
flu_symp_p1_test_50_64 <- 
  flu_symp_p1_test_expand[flu_symp_p1_test_expand$age_group == '50-64',]
flu_symp_p1_test_65 <- 
  flu_symp_p1_test_expand[flu_symp_p1_test_expand$age_group == '65+',]

flu_symp_p1_test_expand_test <- flu_symp_p1_test_expand |> sample_n(20000)

# Predict
flu_symp_p1_test_0_4$pred <- predict(train_model_symptom_age, 
                                        newdata = flu_symp_p1_test_0_4, 
                                        type = "response")
flu_symp_p1_test_5_12$pred <- predict(train_model_symptom_age, 
                                        newdata = flu_symp_p1_test_5_12, 
                                        type = "response")
flu_symp_p1_test_13_17$pred <- predict(train_model_symptom_age, 
                                        newdata = flu_symp_p1_test_13_17, 
                                        type = "response")
flu_symp_p1_test_18_49$pred <- predict(train_model_symptom_age, 
                                        newdata = flu_symp_p1_test_18_49, 
                                        type = "response")
flu_symp_p1_test_50_64$pred <- predict(train_model_symptom_age, 
                                        newdata = flu_symp_p1_test_50_64, 
                                        type = "response")
flu_symp_p1_test_65$pred <- predict(train_model_symptom_age, 
                                        newdata = flu_symp_p1_test_65, 
                                        type = "response")
flu_symp_p1_test_expand_test$pred_symptom_age_time <- predict(train_model, 
                                 newdata = flu_symp_p1_test_expand_test, 
                                 type = "response")
flu_symp_p1_test_expand_test$pred_time <- predict(train_model_time, 
                                             newdata = flu_symp_p1_test_expand_test, 
                                             type = "response")

flu_symp_p1_test_expand_test$pred_fever <- predict(train_model_fever, 
                                                  newdata = flu_symp_p1_test_expand_test, 
                                                  type = "response")

flu_symp_p1_test_expand_test$pred_symp <- predict(train_model_symptom, 
                                                   newdata = flu_symp_p1_test_expand_test, 
                                                   type = "response")
flu_symp_p1_test_expand_test$pred_fever_time <- predict(train_model_fever_time, 
                                                   newdata = flu_symp_p1_test_expand_test, 
                                                   type = "response")
flu_symp_p1_test_expand_test$pred_symp_time <- predict(train_model_symptom_time, 
                                                  newdata = flu_symp_p1_test_expand_test, 
                                                  type = "response")
flu_symp_p1_test_expand_test$pred_symp_age <- predict(train_model_symptom_age, 
                                                       newdata = flu_symp_p1_test_expand_test, 
                                                       type = "response")

pROC::auc(pROC::roc(flu_symp_p1_test_0_4$flu, flu_symp_p1_test_0_4$pred))
pROC::auc(pROC::roc(flu_symp_p1_test_5_12$flu, flu_symp_p1_test_5_12$pred))
pROC::auc(pROC::roc(flu_symp_p1_test_13_17$flu, flu_symp_p1_test_13_17$pred))
pROC::auc(pROC::roc(flu_symp_p1_test_18_49$flu, flu_symp_p1_test_18_49$pred))
pROC::auc(pROC::roc(flu_symp_p1_test_50_64$flu, flu_symp_p1_test_50_64$pred))
pROC::auc(pROC::roc(flu_symp_p1_test_65$flu, flu_symp_p1_test_65$pred))

pROC::auc(pROC::roc(flu_symp_p1_test_expand_test$flu, flu_symp_p1_test_expand_test$pred_symptom_age_time))
pROC::auc(pROC::roc(flu_symp_p1_test_expand_test$flu, flu_symp_p1_test_expand_test$pred_time))
pROC::auc(pROC::roc(flu_symp_p1_test_expand_test$flu, flu_symp_p1_test_expand_test$pred_fever))
pROC::auc(pROC::roc(flu_symp_p1_test_expand_test$flu, flu_symp_p1_test_expand_test$pred_symp))
pROC::auc(pROC::roc(flu_symp_p1_test_expand_test$flu, flu_symp_p1_test_expand_test$pred_fever_time))
pROC::auc(pROC::roc(flu_symp_p1_test_expand_test$flu, flu_symp_p1_test_expand_test$pred_symp_time))
pROC::auc(pROC::roc(flu_symp_p1_test_expand_test$flu, flu_symp_p1_test_expand_test$pred_symp_age))
summary(train_model_time)

ggplot(flu_symp_p1_test_expand_test) + geom_density(aes(x = pred_symp_time, group = flu, color = flu))


flu_symp_p1_test_65$pred_symp <- predict(train_model_symptom, 
                                     newdata = flu_symp_p1_test_65, 
                                     type = "response")
flu_symp_p1_test_65$pred_symp_age <- predict(train_model_symptom_age, 
                                          newdata = flu_symp_p1_test_65, 
                                          type = "response")
pROC::auc(pROC::roc(flu_symp_p1_test_65$flu, flu_symp_p1_test_65$pred_symp))
pROC::auc(pROC::roc(flu_symp_p1_test_65$flu, flu_symp_p1_test_65$pred_symp_age))

# Plot ROC
basicplot <- ggplot(flu_symp_p1_test_expand_test, aes(d = flu, m = pred_time)) + 
  geom_roc(n.cuts=40,labels=FALSE) + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  theme_minimal() + ylab('True Positive Rate') + xlab('False Positive Rate') +
  theme_minimal() + theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size=14),
    plot.title = element_text(size=18),
    legend.text=element_text(size=12),
    legend.title=element_text(size=14))
basicplot 








flu_descriptives_symp <- unbalanced |>
  group_by(flu) |> mutate(Fever = sum(fever),
                          Myalgia = sum(myalgia),
                          `Short Breath` = sum(short_breath),
                          Hypoxemia = sum(hypoxemia),
                          `Chest Pain` = sum(chest_pain),
                          Cough = sum(cough),
                          Bronchitis = sum(bronchitis),
                          `Nausea/Vomit` = sum(nausea_vom),
                          `Sore Throat` = sum(sore_throat),
                          Diarrhea = sum(diarrhea),
                          Fatigue = sum(fatigue),
                          Headache = sum(headache),
                          Congestion = sum(congestion),
                          Sneezing = sum(sneezing),
                          count = 1,
                          patient_sum = sum(count)) |>
  distinct(flu, Fever, Myalgia, `Short Breath`, Hypoxemia, `Chest Pain`,
           Cough, Bronchitis, `Nausea/Vomit`, `Sore Throat`, Diarrhea, Fatigue,
           Headache, Congestion, Sneezing, patient_sum, .keep_all = FALSE) |>
  pivot_longer(cols=c('Fever', 'Myalgia', 'Short Breath', 'Hypoxemia', 'Chest Pain',
                      'Cough', 'Bronchitis', 'Nausea/Vomit', 'Sore Throat', 'Diarrhea', 'Fatigue',
                      'Headache', 'Congestion', 'Sneezing'),
               names_to='Symptom',
               values_to='Count') |>
  mutate(Percent = Count/patient_sum)

ggplot(flu_descriptives_symp[flu_descriptives_symp$flu == 1,], aes(x = Symptom, y = Percent)) +
  geom_bar(aes(), fill = '#C7A939', stat = "identity") + ggtitle('Symptom Percent: Flu Cases') +
  theme_minimal() + xlab('') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(flu_descriptives_symp[flu_descriptives_symp$flu == 0,], aes(x = Symptom, y = Percent)) +
  geom_bar(aes(), fill = '#9086ba', stat = "identity") + ggtitle('Symptom Percent: Non-Cases') +
  theme_minimal() + xlab('') + ylim(0, 0.4) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
