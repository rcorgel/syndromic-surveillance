################################################################################
# File Name: create_intro_figure                                               #
#                                                                              #
# Purpose:   Analyze influenza from medical claims symptoms.                   #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load data                                                      #
#            3. Analyze data                                                   #
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
library(cowplot)
library(lubridate)
library(ISOweek)
library(sf)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1CgyyYhIWFVG2aMITI-HBpKamvwwpQ_ED/CHC_data/')

################
# 2. LOAD DATA #
################

# Load all cause data
all_cause_week <- read.csv('flu_cases/2024-07-18/county_week_ac_v3.csv')
all_cause_season <- read.csv('flu_cases/2024-07-18/county_season_ac_v3.csv')

# Load population data
population <- read.csv('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/raw/co-est2019-alldata.csv')
population$STATE <- ifelse(nchar(population$STATE) == 1, 
                           paste0("0", population$STATE), population$STATE)
population$COUNTY <- ifelse(nchar(population$COUNTY) == 1, 
                           paste0("00", population$COUNTY), population$COUNTY)
population$COUNTY <- ifelse(nchar(population$COUNTY) == 2, 
                            paste0("0", population$COUNTY), population$COUNTY)
population$fips <- paste0(population$STATE, population$COUNTY)
population_2019 <- population[, c('fips', 'POPESTIMATE2019')] 
  
# Load disease data
covid <- read.csv('covid_cases/2024-07-18/county_week_covid_v3.csv')
flu <- read.csv('flu_cases/2024-07-18/county_week_flu_v3.csv')
rsv <- read.csv('rsv_cases/2024-07-26/county_week_rsv_v3.csv')

# Load symptom data
flu_16_17 <- read.csv('syndromic_surveillance/2024-07-26/flu_2016-09_2017-08_p2.csv')
flu_17_18 <- read.csv('syndromic_surveillance/2024-07-26/flu_2017-09_2018-08_p2.csv')
flu_18_19 <- read.csv('syndromic_surveillance/2024-07-31/flu_2018-09_2019-08_p3.csv')
flu_19_20 <- read.csv('syndromic_surveillance/2024-07-26/flu_2019-09_2020-08_p2.csv')
covid_wild <- read.csv('syndromic_surveillance/2024-07-26/covid_wild_p2.csv')
covid_alpha <- read.csv('syndromic_surveillance/2024-07-26/covid_alpha_p2.csv')
covid_delta <- read.csv('syndromic_surveillance/2024-07-26/covid_delta_p2.csv')
covid_omicron <- read.csv('syndromic_surveillance/2024-07-26/covid_omicron_p2.csv')

# CREATE MAP

# Create county to county group x-walk
x_walk <- flu_19_20 |> separate(county_fips, sep = "_", remove = FALSE, 
                                into = c('fips_1', 'fips_2', 'fips_3', 'fips_4',
                                         'fips_5', 'fips_6', 'fips_7', 'fips_8')) |>
  pivot_longer(cols=c('fips_1', 'fips_2', 'fips_3', 'fips_4', 'fips_5', 'fips_6', 
                      'fips_7', 'fips_8'),
               names_to='fips_num',
               values_to='fips') |>
  filter(!is.na(fips)) |> select(c(county_fips, fips))
x_walk$state_fips <- substr(x_walk$county_fips, 1, 2)
# Remove PR and unknown (99999)
x_walk <- x_walk |> filter(state_fips != '72') |>
  filter(county_fips != '99999')

# Load maps
usa_albers_state <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_state.rds'))   # convert to sf
usa_albers_county <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_county.rds')) # convert to sf

# Collapse county map to county group
usa_albers_county_x_walk <- left_join(x_walk, usa_albers_county[, c(5, 6)], by = c('fips' = 'GEOID'))
usa_albers_county_group <- usa_albers_county_x_walk %>% 
  group_by(county_fips) %>%
  summarise(geometry = sf::st_union(geometry)) %>%
  ungroup()

# Collapse population data to county group
x_walk_unique <- x_walk |> distinct(fips, county_fips)
population_2019 <- left_join(population_2019, x_walk_unique, by = c('fips' = 'fips'))
population_2019_group <- population_2019 %>% 
  group_by(county_fips) %>% mutate(population_group = sum(POPESTIMATE2019)) |>
  distinct(county_fips, population_group)

usa_albers_county_group <- left_join(usa_albers_county_group, population_2019_group,
                                     by = c('county_fips' = 'county_fips'))

all_cause_season_2019 <- all_cause_season |> filter(season == '2019-2020')

usa_albers_county_group <- left_join(usa_albers_county_group, all_cause_season_2019,
                                     by = c('county_fips' = 'county_fips'))
usa_albers_county_group$Proportion <- usa_albers_county_group$all_cause / usa_albers_county_group$population_group.y
# Create maps
map <- ggplot() +
  geom_sf(data = st_as_sf(usa_albers_county_group), aes(fill = Proportion, group = county_fips), color= 'black', linewidth = 0.10) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), fill = '#FFFFFF00', color= 'black', linewidth = 0.3) +
  scale_fill_gradient('Patients / Population  \n', low = "white", high = "black", limits = c(0, 1.7)) + ggtitle('Claims Population Coverage, 2019-2020') + 
  theme_void() + theme(legend.position = 'bottom',
                       plot.title = element_text(size = 16, hjust = 0.5),
                       panel.border = element_rect(fill=NA, linewidth = 0.8, color = 'white'),
                       legend.title = element_text(size = 12),
                       legend.text = element_text(size = 12),
                       legend.key.height = unit(0.4, 'cm'),
                       legend.key.width = unit(1.5, "cm")) 
map


# Create Disease time series
flu  <- mutate(flu, conf_flu = ifelse(is.na(conf_flu), 0, conf_flu))
flu$conf_flu[flu$conf_flu == '<=5'] <- 
  sample(1:5, length(flu$conf_flu[flu$conf_flu == '<=5']), replace= TRUE)
flu$Virus <- 'Influenza'
flu <- rename(flu, 'count' = 'conf_flu')

covid  <- mutate(covid, covid = ifelse(is.na(covid), 0, covid))
covid$covid[covid$covid == '<=5'] <- 
  sample(1:5, length(covid$covid[covid$covid == '<=5']), replace= TRUE)
covid$Virus <- 'COVID-19'
covid <- rename(covid, 'count' = 'covid')

rsv  <- mutate(rsv, rsv = ifelse(is.na(rsv), 0, rsv))
rsv$rsv[rsv$rsv == '<=5'] <- 
  sample(1:5, length(rsv$rsv[rsv$rsv == '<=5']), replace= TRUE)
rsv$Virus <- 'RSV'
rsv <- rename(rsv, 'count' = 'rsv')

viruses <- rbind(flu, covid, rsv)

viruses$state_fips <- substr(viruses$county_fips, 1, 2)
# Remove PR and unknown (99999)
viruses <- viruses |> filter(state_fips != '72') |>
  filter(county_fips != '99999')
library(ISOweek)
viruses$week_date <- sub("(\\d{4}-)(\\d{2})", "\\1W\\2-1", viruses$year_week)
viruses$week_date <- ISOweek2date(viruses$week_date)

viruses_us <- viruses |> group_by(week_date, Virus) |>
  mutate(count_sum = sum(as.numeric(count))) |>
  distinct(week_date, Virus, count_sum) |> ungroup() |>
  group_by(Virus) |>
  mutate(average = mean(count_sum),
         sd = sd(count_sum),
         z = (count_sum - average) / sd)


ggplot(viruses_us, aes(x = week_date, y = z)) + 
  geom_line(aes(color = Virus)) + theme_minimal() +
  theme(legend.position = 'bottom') + ylab('Virus Activity') + 
  xlab('Week') + ggtitle('Influenza, COVID-19, and RSV Activity, 2016-2024')


flu_symp <- rbind(flu_16_17, flu_17_18, flu_18_19, flu_19_20)
flu_symp$patient_count_imp <- flu_symp$patient_count
flu_symp$patient_count_imp[flu_symp$patient_count_imp == '<=5'] <- 
  sample(1:5, length(flu_symp$patient_count_imp[flu_symp$patient_count_imp == '<=5']), replace= TRUE)
flu_symp$patient_count_imp <- as.numeric(flu_symp$patient_count_imp)
symptoms_count <- flu_symp |> group_by()
