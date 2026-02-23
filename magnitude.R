
rm(list = ls())


# Load libraries
library(tidyverse)
library(lubridate)
library(car)
library(plotROC)
library(pROC)
library(cutpointr)
library(zoo)
library(ISOweek)

# Set seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')


fev_model <- readRDS("./tmp/fever_model.rds")
all_model <- readRDS("./tmp/symptoms_model.rds")
age_model <- readRDS("./tmp/age_model.rds")
time_model <- readRDS("./tmp/time_model.rds")

# Load part 2 data
flu_proc_symp <- readRDS('./tmp/flu_p2_dat_proc_symp.rds')

# Create empty lists to fill
flu_nat_filt <- list()

for (i in 1:4) {
  # Print the progress and set season data
  print(i)
  flu <- as.data.frame(flu_proc_symp[i])
  
  # add month and age variables
  flu$month <- month(flu$week_date)
  flu$age_group <- factor(flu$age_grp,
                          levels = c('0', '1', '2', '3', '4', '5'),
                          labels = c('0-4', '5-12', '13-17', '18-49', '50-64', '65+'))
  
  # Restrict data to the influenza season
  flu <- flu |> dplyr::filter(month > 9 | month < 6)
  
  # Collapse to the national level
  flu_nat <- flu |> group_by(flu, fever, cough, sore_throat, myalgia, 
                             short_breath, hypoxemia, nausea_vom, 
                             bronchitis, chest_pain, diarrhea,
                             fatigue, headache, congestion, sneezing,
                             age_group, month) |>
    mutate(patient_sum = sum(patient_count_imp)) |>
    distinct(flu, fever, cough, sore_throat, myalgia, 
             short_breath, hypoxemia, nausea_vom, 
             bronchitis, chest_pain, diarrhea,
             fatigue, headache, congestion, sneezing,
             age_group, month, patient_sum)
  
  # Fill lists
  flu_nat_filt[[i]] <- as.data.frame(flu_nat)
  remove(flu, flu_nat)
}

flu_16_17 <- flu_nat_filt[[4]]
flu_16_17$pred_fev <- predict(fev_model, 
                              newdata = flu_16_17, 
                              type = "response")
flu_16_17$pred_all <- predict(all_model, 
                              newdata = flu_16_17, 
                              type = "response")
flu_16_17$pred_age <- predict(age_model, 
                              newdata = flu_16_17, 
                              type = "response")
flu_16_17$pred_time <- predict(time_model, 
                               newdata = flu_16_17, 
                               type = "response")
flu_16_17$pred_fev_count <- flu_16_17$pred_fev * flu_16_17$patient_sum
flu_16_17$pred_all_count <- flu_16_17$pred_all * flu_16_17$patient_sum
flu_16_17$pred_age_count <- flu_16_17$pred_age * flu_16_17$patient_sum
flu_16_17$pred_time_count <- flu_16_17$pred_time * flu_16_17$patient_sum

flu_16_17$pred_fev_count <- ifelse(flu_16_17$pred_fev > 0.5, 1, 0) * flu_16_17$patient_sum
flu_16_17$pred_all_count <- ifelse(flu_16_17$pred_all > 0.5, 1, 0) * flu_16_17$patient_sum
flu_16_17$pred_age_count <- ifelse(flu_16_17$pred_age > 0.42, 1, 0) * flu_16_17$patient_sum
flu_16_17$pred_time_count <- ifelse(flu_16_17$pred_time > 0.5, 1, 0) * flu_16_17$patient_sum


all_cause <- read_csv('raw/county_season_ac_v2.csv')
all_cause_week <- read_csv('raw/county_week_ac_v2.csv')



all_cause$patient_count_imp <- all_cause$all_cause
all_cause$patient_count_imp[all_cause$patient_count_imp == '<=5'] <- 
  sample(1:5, length(all_cause$patient_count_imp[all_cause$patient_count_imp == '<=5']), 
         replace= TRUE)
# Convert imputed count to numeric
all_cause$patient_count_imp <- as.numeric(all_cause$patient_count_imp)

all_cause_16_17 <- all_cause |> dplyr::filter(season == '2019-2020') |>
  mutate(state_fips = substr(county_fips, 1, 2)) |>
  filter(!(state_fips %in% c("72","99","")))


all_cause_week$patient_count_imp <- all_cause_week$all_cause
all_cause_week$patient_count_imp[all_cause_week$patient_count_imp == '<=5'] <- 
  sample(1:5, length(all_cause_week$patient_count_imp[all_cause_week$patient_count_imp == '<=5']), 
         replace= TRUE)
# Convert imputed count to numeric
all_cause_week$patient_count_imp <- as.numeric(all_cause_week$patient_count_imp)

all_cause_16_17_week <- all_cause_week |> 
  mutate(state_fips = substr(county_fips, 1, 2)) |>
  filter(!(state_fips %in% c("72","99",""))) |>
  # Convert week date to date
  mutate(week_date = sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", year_week)) |>
  mutate(week_date = ISOweek2date(week_date)) |>
  mutate(year = year(week_date),
         month = month(week_date))


  
all_cause_16_17_week$season <- ifelse(all_cause_16_17_week$year == 2016 & all_cause_16_17_week$month > 8, '2016-2017', '')
all_cause_16_17_week$season <- ifelse(all_cause_16_17_week$year == 2017 & all_cause_16_17_week$month < 9, '2016-2017', all_cause_16_17_week$season)
all_cause_16_17_week$season <- ifelse(all_cause_16_17_week$year == 2017 & all_cause_16_17_week$month > 8, '2017-2018', all_cause_16_17_week$season)
all_cause_16_17_week$season <- ifelse(all_cause_16_17_week$year == 2018 & all_cause_16_17_week$month < 9, '2017-2018', all_cause_16_17_week$season)
all_cause_16_17_week$season <- ifelse(all_cause_16_17_week$year == 2018 & all_cause_16_17_week$month > 8, '2018-2019', all_cause_16_17_week$season)
all_cause_16_17_week$season <- ifelse(all_cause_16_17_week$year == 2019 & all_cause_16_17_week$month < 9, '2018-2019', all_cause_16_17_week$season)
all_cause_16_17_week$season <- ifelse(all_cause_16_17_week$year == 2019 & all_cause_16_17_week$month > 8, '2019-2020', all_cause_16_17_week$season)
all_cause_16_17_week$season <- ifelse(all_cause_16_17_week$year == 2020 & all_cause_16_17_week$month < 9, '2019-2020', all_cause_16_17_week$season)

all_cause_16_17_week <- all_cause_16_17_week |> dplyr::filter(season == '2019-2020')

sum(all_cause_16_17_week$patient_count_imp) / sum(all_cause_16_17$patient_count_imp)

all_cause_16_17_week <- all_cause_16_17_week |>
  dplyr::filter(month > 9 | month < 6)

# ALL CAUSE RATIO
# sum(all_cause_16_17_week$patient_count_imp) / sum(all_cause_16_17$patient_count_imp)

# SYMPTOMATIC RATIO
# flu_symp <- as.data.frame(flu_proc_symp[1])
# flu_symp$season <- ifelse(year(flu_symp$week_date) == 2016 & month(flu_symp$week_date) > 8, '2016-2017', '')
# flu_symp$season <- ifelse(year(flu_symp$week_date) == 2017 & month(flu_symp$week_date) < 9, '2016-2017', flu_symp$season)
# flu_symp_16_17 <- flu_symp |> dplyr::filter(season == '2016-2017')
# 
# sum(flu_symp_16_17$patient_count_imp)
# 
# flu_p1 <- readRDS('./tmp/flu_p1_dat_proc.rds')
# flu_symp_p1 <- flu_p1 |> dplyr::filter(symp_count > 0)
# flu_symp_p1$season <- ifelse(year(flu_symp_p1$month_date) == 2016 & month(flu_symp_p1$month_date) > 8, '2016-2017', '')
# flu_symp_p1$season <- ifelse(year(flu_symp_p1$month_date) == 2017 & month(flu_symp_p1$month_date) < 9, '2016-2017', flu_symp_p1$season)
# flu_symp_16_17_p1 <- flu_symp_p1 |> dplyr::filter(season == '2016-2017')
# 
# sum(flu_symp_16_17$patient_count_imp * flu_symp_16_17$flu) / sum(flu_symp_16_17_p1$patient_count_imp * flu_symp_16_17_p1$flu)
# 
# sum(flu_16_17$pred_age_count) / sum(flu_16_17$flu * flu_16_17$patient_sum) 

percent <- (sum(flu_16_17$flu * flu_16_17$patient_sum)  * runif(10000, min = 1/0.58, max = 1/0.42)) / 
  sum(all_cause_16_17_week$patient_count_imp)
flu <- data.frame(Cat = 'FLU', percent)
percent <- (sum(flu_16_17$pred_fev_count) * runif(10000, min = 1/0.58, max = 1/0.42)) / 
  sum(all_cause_16_17_week$patient_count_imp)
fever <- data.frame(Cat = 'FEVER', percent)
percent <- (sum(flu_16_17$pred_all_count) * runif(10000, min = 1/0.58, max = 1/0.42)) / 
    sum(all_cause_16_17_week$patient_count_imp)
all <- data.frame(Cat = 'ALL', percent)
percent <- (sum(flu_16_17$pred_age_count) * runif(10000, min = 1/0.58, max = 1/0.42)) / 
  sum(all_cause_16_17_week$patient_count_imp)
age <- data.frame(Cat = 'ALL + AGE', percent)
percent <- (sum(flu_16_17$pred_time_count) * runif(10000, min = 1/0.58, max = 1/0.42)) / 
  sum(all_cause_16_17_week$patient_count_imp)
time <- data.frame(Cat = 'ALL + AGE + TIME', percent)

mag <- rbind(flu, fever, all, age, time)

mag_2019 <- mag |> group_by(Cat) |>
  mutate(min = min(percent), max = max(percent)) |>
  distinct(Cat, min, max)

mag_2016$Season <- '2016 - 2017'
mag_2017$Season <- '2017 - 2018'
mag_2018$Season <- '2018 - 2019'
mag_2019$Season <- '2019 - 2020'




mag_all <- rbind(mag_2016, mag_2017, mag_2018, mag_2019)


ggplot() + 
  geom_errorbar(data=mag_all, aes(x=Season, ymin=min, ymax=max, group = Cat, color = Cat), width = 0.1, linewidth = 1, 
                position = position_dodge2(width = )) +
  coord_flip() + xlab('') 



library(ggridges)
ggplot(mag, aes(x = percent, y = Cat, fill = Cat)) + geom_density_ridges(scale = 1)


# CDC ESTIMATES
season <- c('2011-2012','2012-2013', '2013-2014', '2014-2015',
            '2015-2016', '2016-2017', '2017-2018', '2018-2019',
            '2019-2020')

population <- c(311.6, 313.9, 316.1, 318.4, 320.7, 323.1, 325.1,
                326.8, 328.3)

cdc_flu <- c(9.3, 34, 30, 30, 24, 29, 41, 29, 36)

cdc_dat <- data.frame(season, population, cdc_flu)

cdc_dat$prev <- cdc_dat$cdc_flu / cdc_dat$population
cdc_dat$x <- "CDC Estimated Prevalence"

range <- data.frame(lab = 'CDC Estimated Prevalence', high = 0.10, low = 0.025)


ggplot() + 
  geom_errorbar(data=range, aes(x=lab, ymin=low, ymax=high), width = 0.1, linewidth = 1) +
  coord_flip() + xlab('') +
  geom_point(data=cdc_dat, aes(x = x, y = prev, color = season), size = 4, alpha = 0.8)
  

  d=data.frame(drink=c("coffee","tea","water"), mean=c(3,6,2), lower=c(2.6,5.6,1.8), upper=c(3.5,6.3,2.8))




# Create empty lists to fill
flu_nat_filt_county <- list()

for (i in 1:4) {
  # Print the progress and set season data
  print(i)
  flu <- as.data.frame(flu_proc_symp[i])
  
  # add month and age variables
  flu$month <- month(flu$week_date)
  flu$age_group <- factor(flu$age_grp,
                          levels = c('0', '1', '2', '3', '4', '5'),
                          labels = c('0-4', '5-12', '13-17', '18-49', '50-64', '65+'))
  
  # Restrict data to the influenza season
  flu <- flu |> dplyr::filter(month > 9 | month < 6)
  
  # Collapse to the national level
  flu_nat <- flu |> group_by(county_fips, week_date, flu, fever, cough, sore_throat, myalgia, 
                             short_breath, hypoxemia, nausea_vom, 
                             bronchitis, chest_pain, diarrhea,
                             fatigue, headache, congestion, sneezing,
                             age_group, month) |>
    mutate(patient_sum = sum(patient_count_imp)) |>
    distinct(county_fips, flu, fever, cough, sore_throat, myalgia, 
             short_breath, hypoxemia, nausea_vom, 
             bronchitis, chest_pain, diarrhea,
             fatigue, headache, congestion, sneezing,
             age_group, month, patient_sum)
  
  # Fill lists
  flu_nat_filt_county[[i]] <- as.data.frame(flu_nat)
  remove(flu, flu_nat)
}



flu <- flu_nat_filt_county[[4]]

flu$pred_fev <- predict(fev_model, 
                              newdata = flu, 
                              type = "response")
flu$pred_all <- predict(all_model, 
                              newdata = flu, 
                              type = "response")
flu$pred_age <- predict(age_model, 
                              newdata = flu, 
                              type = "response")

flu$pred_age_count <- ifelse(flu$pred_age > 0.42, 1, 0) * flu$patient_sum


flu_county <- flu |> group_by(county_fips) |>
  mutate(pred_sum = sum(pred_age_count),
         flu_sum = sum(flu * patient_sum)) |>
  distinct(county_fips, pred_sum, flu_sum)


all_cause_16_17_sum <- all_cause_16_17_week |> group_by(county_fips) |>
  mutate(cause_sum = sum(patient_count_imp)) |>
  distinct(county_fips,cause_sum)

flu_county <- left_join(flu_county, all_cause_16_17_sum, by = c('county_fips'))

flu_county$inci <- (flu_county$flu_sum * 2) / flu_county$cause_sum
flu_county$inci_pred <- (flu_county$pred_sum * 2) / flu_county$cause_sum


ggplot(flu_county) + geom_density(aes(x = inci_pred))

library(sf)
# Load maps
usa_albers_state <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_state.rds'))   # convert to sf
usa_albers_county <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_county.rds')) # convert to sf

x_walk <- readRDS('./tmp/county_group_xwalk.rds')

# Collapse county map to county group
usa_albers_county_x_walk <- left_join(x_walk, usa_albers_county[, c(5, 6)], by = c('fips' = 'GEOID'))
usa_albers_county_group <- usa_albers_county_x_walk %>% 
  group_by(county_fips) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()

usa_albers_county_group <- left_join(usa_albers_county_group, flu_county,
                                     by = c('county_fips' = 'county_fips'))


# Load population data
county_pop <- readRDS('tmp/county_group_pop.rds')

# Restrict to counties with 10,000 or more individuals
usa_albers_county_group <- left_join(usa_albers_county_group, county_pop, 
                         by = c('county_fips' = 'county_fips'))
usa_albers_county_group_filt <- usa_albers_county_group |> filter(county_group_pop > 20000)


# Create maps
map <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_group), aes(fill = inci_pred, group = county_fips), color= 'black', linewidth = 0.10) +
  scale_fill_gradient('Prev \n', low = "white", high = "black") + ggtitle('Flu Prev, 2019-2020') + 
  theme_void() + theme(legend.position = 'bottom',
                       plot.title = element_text(size = 16, hjust = 0.5),
                       panel.border = element_rect(fill=NA, linewidth = 0.8, color = 'white'),
                       legend.title = element_text(size = 12),
                       legend.text = element_text(size = 12),
                       legend.key.height = unit(0.4, 'cm'),
                       legend.key.width = unit(1.5, "cm")) 
map


ggplot(flu_county) + geom_point(aes(x = inci_pred, y = inci))
cor(flu_county$inci_pred, flu_county$inci)

# Create empty lists to fill
flu_nat_filt_county_week <- list()

for (i in 1:4) {
  # Print the progress and set season data
  print(i)
  flu <- as.data.frame(flu_proc_symp[i])
  
  # add month and age variables
  flu$month <- month(flu$week_date)
  flu$age_group <- factor(flu$age_grp,
                          levels = c('0', '1', '2', '3', '4', '5'),
                          labels = c('0-4', '5-12', '13-17', '18-49', '50-64', '65+'))
  
  # Restrict data to the influenza season
  flu <- flu |> dplyr::filter(month > 9 | month < 6)
  
  # Collapse to the national level
  flu_nat <- flu |> group_by(state_fips, week_date, flu, fever, cough, sore_throat, myalgia, 
                             short_breath, hypoxemia, nausea_vom, 
                             bronchitis, chest_pain, diarrhea,
                             fatigue, headache, congestion, sneezing,
                             age_group, month) |>
    mutate(patient_sum = sum(patient_count_imp)) |>
    distinct(state_fips, flu, fever, cough, sore_throat, myalgia, 
             short_breath, hypoxemia, nausea_vom, 
             bronchitis, chest_pain, diarrhea,
             fatigue, headache, congestion, sneezing,
             age_group, month, patient_sum)
  
  # Fill lists
  flu_nat_filt_county_week[[i]] <- as.data.frame(flu_nat)
  remove(flu, flu_nat)
}

flu <- flu_nat_filt_county_week[[1]]

flu$pred_fev <- predict(fev_model, 
                        newdata = flu, 
                        type = "response")
flu$pred_all <- predict(all_model, 
                        newdata = flu, 
                        type = "response")
flu$pred_age <- predict(age_model, 
                        newdata = flu, 
                        type = "response")
flu$pred_time <- predict(time_model, 
                         newdata = flu, 
                         type = "response")

flu$pred_age_count <- ifelse(flu$pred_age > 0.42, 1, 0) * flu$patient_sum

flu_state <- flu |> group_by(state_fips, week_date) |>
  mutate(pred_sum = sum(pred_age_count)) |>
  distinct(state_fips, week_date, pred_sum)


all_cause_16_17_week_state <- all_cause_16_17_week |> group_by(state_fips, week_date) |>
  mutate(cause_sum = sum(patient_count_imp)) |>
  distinct(state_fips, week_date, cause_sum)

flu_state <- left_join(flu_state, all_cause_16_17_week_state, by = c('state_fips', 'week_date'))

flu_state$inci <- flu_state$pred_sum / flu_state$cause_sum

ggplot(flu_state) + geom_line(aes(y = pred_sum , x = week_date), color = 'black') +
  # geom_line(aes(y = flu , x = week_date), color = 'red') +
  facet_wrap(vars(state_fips), scales = "free")

ili_state <- read.csv('./raw/ILINet_State.csv', skip = 1, header = TRUE, sep = ',')

ili_state$WEEK_2 <- ifelse(nchar(ili_state$WEEK) == 1, paste0('0', ili_state$WEEK), ili_state$WEEK)

ili_state$year_week <- paste(ili_state$YEAR, ili_state$WEEK_2, sep = '-')

ili_state$week_date = sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", ili_state$year_week)
ili_state$week_date = ISOweek2date(ili_state$week_date)

library(usmap)

ili_state$state_fips <- fips(ili_state$REGION)
ili_state[, c(6, 18, 19)]
flu_state <- left_join(flu_state, ili_state[, c(6, 18, 19)], by = c('state_fips', 'week_date'))
flu_state$ili_rate <- as.numeric(flu_state$X.UNWEIGHTED.ILI)
flu_state <- flu_state |> group_by(state_fips) |>
  mutate(ili_norm = (ili_rate - mean(ili_rate)) / sd(ili_rate),
         pred_norm = (pred_sum - mean(pred_sum)) / sd(pred_sum)) 

ggplot(flu_state) + geom_line(aes(y = ili_norm , x = week_date), color = 'black') +
  geom_line(aes(y = pred_norm , x = week_date), color = 'red') +
  # geom_line(aes(y = flu , x = week_date), color = 'red') +
  facet_wrap(vars(state_fips), scales = "free")
