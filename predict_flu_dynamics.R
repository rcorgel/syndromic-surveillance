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

# Set seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

# Load processed data
# Part 1 full data
flu <- readRDS('./tmp/flu_p1_dat_proc.rds')

# Part 2 major counties
flu_p2 <- readRDS('./tmp/flu_p2_dat_proc_select_counties.rds')

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

# Add month vars to Part 1 data
flu$month <- month(flu$month_date)
flu$jan <- ifelse(flu$month == 1, 1, 0)
flu$feb <- ifelse(flu$month == 2, 1, 0)
flu$mar <- ifelse(flu$month == 3, 1, 0)
flu$apr <- ifelse(flu$month == 4, 1, 0)
flu$may <- ifelse(flu$month == 5, 1, 0)
flu$jun <- ifelse(flu$month == 6, 1, 0)
flu$jul <- ifelse(flu$month == 7, 1, 0)
flu$aug <- ifelse(flu$month == 8, 1, 0)
flu$sep <- ifelse(flu$month == 9, 1, 0)
flu$oct <- ifelse(flu$month == 10, 1, 0)
flu$nov <- ifelse(flu$month == 11, 1, 0)
flu$dec <- ifelse(flu$month == 12, 1, 0)

# Create lagged fever, sore throat, and cough percent variables for Part 2
lagged_vars <- flu_p2 |> group_by(county_fips, month_date) |>
  mutate(fever_perc_lag = sum(fever*patient_count_imp) / sum(patient_count_imp),
         cough_perc_lag = sum(cough*patient_count_imp) / sum(patient_count_imp),
         sore_throat_perc_lag = sum(sore_throat*patient_count_imp) / sum(patient_count_imp)) |>
  distinct(county_fips, month_date, fever_perc_lag, cough_perc_lag, sore_throat_perc_lag, .keep_all = FALSE) |>
  mutate(lag_month_date = month_date %m+% months(1)) |> ungroup() |>
  select(-c(month_date))

# Merge on lagged variables for Part 2
flu_p2 <- left_join(flu_p2, lagged_vars, by = c('county_fips' = 'county_fips',
                                                'month_date' = 'lag_month_date'))

# Remove asymptomatic cases for disease and non-diseased from part 1 and 2
flu_symp <- flu |> dplyr::filter(symp_count > 0)


# PLOTS FOR PRESENTATION
flu_week <- flu_p2 |> filter(county_fips == '36061') |> 
  group_by(week_date) |> mutate(flu_sum = sum(flu*patient_count_imp)) |>
  distinct(week_date, flu_sum, .keep_all = FALSE)

symptom_week <- flu_p2 |> filter(county_fips == '36061') |> 
  group_by(week_date) |> 
  mutate(fever_sum = sum(fever*patient_count_imp),
         cough_sum = sum(cough*patient_count_imp),
         sore_throat_sum = sum(sore_throat*patient_count_imp)) |>
  distinct(week_date, fever_sum, cough_sum, sore_throat_sum, .keep_all = FALSE)
  
ggplot() + 
  geom_line(data = flu_week, 
            aes(x = week_date, y = flu_sum), color = '#C7A939', linewidth = 1.0) + 
  ylim(0, 2000) + ylab('Flu Count') + xlab('Week') + ggtitle('Flu Dynamics: New York, NY') +
  theme_minimal()

colors <- c("Fever" = "#c36272", "Cough" = "#5a9374", "Sore Throat" = "#347dab")
ggplot() + 
  geom_line(data = symptom_week, 
            aes(x = week_date, y = sore_throat_sum, color = 'Sore Throat'), linewidth = 1.0, alpha = 0.9) + 
  geom_line(data = symptom_week, 
            aes(x = week_date, y = cough_sum, color = 'Cough'), linewidth = 1.0, alpha = 0.9) +
  geom_line(data = symptom_week, 
            aes(x = week_date, y = fever_sum, color = 'Fever'), linewidth = 1.0, alpha = 0.9) +
  ylim(0, 3700) + ggtitle('Symptom Dynamics: New York, NY') +
  theme_minimal() +
  labs(x = "Week",
       y = "Symptom Count",
       color = "Legend") +
  scale_color_manual(values = colors) + theme(legend.position = 'bottom')

# Collapse data to claim count
ppl_2018 <- flu |> filter(year(month_date) == 2018) |>
  group_by(county_fips) |> mutate(ppl_count = sum(patient_count_imp)) |>
  distinct(county_fips, ppl_count, .keep_all = FALSE)

# Create county to county group x-walk
x_walk <- ppl_2018 |> separate(county_fips, sep = "_", remove = FALSE, 
                                      into = c('fips_1', 'fips_2', 'fips_3', 'fips_4',
                                               'fips_5', 'fips_6', 'fips_7', 'fips_8')) |>
  pivot_longer(cols=c('fips_1', 'fips_2', 'fips_3', 'fips_4', 'fips_5', 'fips_6', 
                      'fips_7', 'fips_8'),
               names_to='fips_num',
               values_to='fips') |>
  filter(!is.na(fips))
x_walk$state_fips <- substr(x_walk$county_fips, 1, 2)
# Remove PR and unknown (99999)
x_walk <- x_walk |> filter(state_fips != '72') |>
  filter(county_fips != '99999')

###################
# 3. EXPLORE DATA #
###################

# Load maps
library(sf)
usa_albers_state <- st_as_sf(readRDS('./tmp/usa_albers_state.rds'))   # convert to sf
usa_albers_county <- st_as_sf(readRDS('./tmp/usa_albers_county.rds')) # convert to sf

# Collapse county map to county group
usa_albers_county_x_walk <- left_join(x_walk, usa_albers_county[, c(5, 6)], by = c('fips' = 'GEOID'))
usa_albers_county_group <- usa_albers_county_x_walk %>% 
  group_by(county_fips, ppl_count) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()

# Map
map <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_group), aes(fill = log(ppl_count), group = county_fips), color= 'black', linewidth = 0.15) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.4) +
  scale_fill_continuous('Log(Patients)', high = 'black', low = 'white') +
  theme_void() + ggtitle('Log Number of Patients, 2018') + theme(legend.position = 'right',
                                     plot.title = element_text(size = 14),
                                     panel.border = element_rect(fill=NA, linewidth = 0.8, color = 'white'),
                                     legend.title = element_text(size = 12),
                                     legend.text = element_text(size = 12),
                                     legend.key.size = unit(0.5, 'cm')) 

map

ppl_2018 <- flu |> 
  group_by(county_fips) |> mutate(ppl_count = sum(patient_count_imp)) |>
  distinct(county_fips, ppl_count, .keep_all = FALSE)

flu_descriptives_symp <- flu |> mutate(
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
                          patient_sum = sum(patient_count_imp)) |>
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
  theme_minimal() + xlab('') + ylim(0, 0.3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#####################

flu_p2_symp <- flu_p2 |>
  filter(symp_count > 0) |>
  mutate(state_fips = substr(county_fips, 1, 2))

# Split up counties into test and train
counties <- flu_symp |> group_by(county_fips) |>
  distinct(county_fips)

sample_counties <- sample(c(TRUE, FALSE), nrow(counties), replace=TRUE, prob=c(0.5, 0.5))
counties_train  <- counties[sample_counties, ]
counties_test   <- counties[!sample_counties, ]
counties_train$train <- 1

# Split data
flu_symp <- left_join(flu_symp, counties_train, by = c('county_fips' = 'county_fips'))
flu_symp_train <- flu_symp |> filter(train == 1)
flu_symp_test <- flu_symp |>
  filter(is.na(train))

# flu_county_select <- flu_symp |> 
#   dplyr::filter(county_fips == "17031" | county_fips == "53033" | 
#                   county_fips == "06037" | county_fips == "04013" |
#                   county_fips == "48113" | county_fips == "36061" |
#                   county_fips == "11001" | county_fips == "13121" |
#                   county_fips == "29510")

# Bring data back to season level
flu_16_17 <- flu_symp |>
  filter(month_date > as.Date("2016-08-31")) |>
  filter(month_date < as.Date("2017-09-01"))

# Add season
flu_16_17$season <- ifelse(flu_16_17$month < 9 & flu_16_17$month > 3, 0, 1)

# Sample from the full part 1 data
as_is <- flu_symp_train |> uncount(patient_count_imp) |> 
  sample_n(100000)

cases <- flu_symp_train |> filter(flu == 1) |> uncount(patient_count_imp) |> 
  sample_n(50000)

non_cases <- flu_symp_train |> filter(flu == 0) |> uncount(patient_count_imp) |> 
  sample_n(50000)

balanced <- rbind(cases, non_cases)

# Split into training and testing
sample <- sample(c(TRUE, FALSE), nrow(balanced), replace=TRUE, prob=c(0.5, 0.5))
flu_train  <- balanced[sample, ]
flu_test   <- balanced[!sample, ]

# Train the model
train_model <- glm(flu ~ fever + cough + sore_throat +
                     + myalgia + short_breath + hypoxemia +
                     nausea_vom + bronchitis + chest_pain + 
                     diarrhea + fatigue + headache + 
                     congestion + sneezing + 0 + cough_perc_lag +
                     fever_perc_lag + sore_throat_perc_lag,
                   family = binomial(link = "logit"),
                   data = balanced)

summary(train_model)
vif(train_model)

# Predict testing data
flu_test$pred <- predict(train_model, newdata = flu_test, type = "response")

pROC::auc(pROC::roc(flu_test$flu, flu_test$pred))


# Plot ROC
basicplot <- ggplot(flu_test, aes(d = flu, m = pred)) + 
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

flu_p2_symp$season <- ifelse(flu_p2_symp$month < 9 & flu_p2_symp$month > 3, 0, 1)

# 






  
# Reproduce on Part 2 data
flu_p2_symp$pred <- predict(train_model, newdata = flu_p2_symp, type = "response", se.fit = TRUE)$fit
flu_p2_symp$se <- predict(train_model, newdata = flu_p2_symp, type = "response", se.fit = TRUE)$se.fit
flu_p2_symp$upr <- flu_p2_symp$pred + (1.96 * flu_p2_symp$se)
flu_p2_symp$lwr <- flu_p2_symp$pred - (1.96 * flu_p2_symp$se)

p2_week_county_sum <- flu_p2_symp |> group_by(week_date, county_fips) |> 
  mutate(sum_flu = sum(flu),
         sum_pred = sum(pred),
         sum_up = sum(upr),
         sum_low = sum(lwr)) |>
  distinct(week_date, county_fips, sum_flu, sum_pred, sum_up, sum_low, .keep_all = FALSE) |>
  ungroup() |>
  group_by(county_fips) |>
  mutate(roll_flu = rollmean(x = sum_flu, k = 3, fill = NA, align = "right"),
         roll_pred = rollmean(x = sum_pred, k = 3, fill = NA, align = "right"),
         roll_up = rollmean(x = sum_up, k = 3, fill = NA, align = "right"),
         roll_low = rollmean(x = sum_low, k = 3, fill = NA, align = "right")) 

ggplot(p2_week_county_sum, aes(x=week_date, y=roll_flu)) + 
  geom_line(aes(y = roll_flu), color = '#9086ba') + 
  geom_line(aes(y = roll_pred), color = '#5a9374') +
  geom_ribbon(aes(x=week_date, ymin=roll_low, ymax=roll_up), fill="#5a9374", alpha=0.8) +
  facet_wrap(vars(county_fips), scales = "free") +
  theme_minimal() + ylab('Number of Cases') + xlab('Time (weeks)')


test <- p2_week_county_sum |>
  filter(week_date > as.Date("2016-08-31")) |>
  filter(week_date < as.Date("2017-09-01")) |>
  group_by(county_fips) |>
  mutate(pred = sum(sum_pred)) |> distinct(county_fips, pred)




