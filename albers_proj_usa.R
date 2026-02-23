################################################################################
# File Name: albers_proj_usa                                                   #
#                                                                              #
# Purpose:   Create map of the United States with an Albers projection.        #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load data                                                      #
#            3. Project map                                                    #
#            4. Save maps                                                      #
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
library(splitstackshape)
library(sf)
#library(albersusa)
#library(maptools)
library(sp)
library(mapproj)
#library(rgeos)
#library(rgdal)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

################
# 2. LOAD DATA #
################

# Load county-level shape files
cbsa <- read_sf(dsn = './raw/tl_2025_us_cbsa/', 
                  layer = 'tl_2025_us_cbsa')

# Load county-level shape files
county <- read_sf(dsn = './raw/cb_2021_us_county_20m/', 
                  layer = 'cb_2021_us_county_20m')

# Load state-level shape file
state <- read_sf(dsn = './raw/cb_2021_us_state_20m/', 
                  layer = 'cb_2021_us_state_20m')

##################
# 3. PROJECT MAP #
##################

##################
# CBSA-level map #
##################

# Project data
cbsa_trans <- st_transform(cbsa, crs = st_crs("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"))
cbsa_trans <- as(cbsa_trans, 'Spatial')

# Find Alasks and Hawaii
cbsa_trans$AK <- ifelse(grepl(", AK", cbsa_trans$NAME, ignore.case = TRUE) == TRUE, 1, 0)
cbsa_trans$HI <- ifelse(grepl(", HI", cbsa_trans$NAME, ignore.case = TRUE) == TRUE, 1, 0)
cbsa_trans$PR <- ifelse(grepl(", PR", cbsa_trans$NAME, ignore.case = TRUE) == TRUE, 1, 0)
cbsa_trans$AS <- ifelse(grepl(", AS", cbsa_trans$NAME, ignore.case = TRUE) == TRUE, 1, 0)
cbsa_trans$GU <- ifelse(grepl(", GU", cbsa_trans$NAME, ignore.case = TRUE) == TRUE, 1, 0)
cbsa_trans$MH <- ifelse(grepl(", MH", cbsa_trans$NAME, ignore.case = TRUE) == TRUE, 1, 0)
cbsa_trans$FM <- ifelse(grepl(", FM", cbsa_trans$NAME, ignore.case = TRUE) == TRUE, 1, 0)
cbsa_trans$MP <- ifelse(grepl(", MP", cbsa_trans$NAME, ignore.case = TRUE) == TRUE, 1, 0)
cbsa_trans$PW <- ifelse(grepl(", PW", cbsa_trans$NAME, ignore.case = TRUE) == TRUE, 1, 0)
cbsa_trans$VI <- ifelse(grepl(", VI", cbsa_trans$NAME, ignore.case = TRUE) == TRUE, 1, 0)

# Alaska
alaska <- cbsa_trans[cbsa_trans$AK==1,]
alaska <- elide(alaska, rotate=-50)
alaska <- elide(alaska, scale=max(apply(bbox(alaska), 1, diff)) / 2.3)
alaska <- elide(alaska, shift=c(-2300000, -2500000))
proj4string(alaska) <- proj4string(cbsa_trans)

# Hawaii
hawaii <- cbsa_trans[cbsa_trans$HI==1,]
hawaii <- elide(hawaii, rotate=-35)
hawaii <- elide(hawaii, shift=c(5400000, -1400000))
proj4string(hawaii) <- proj4string(cbsa_trans)

# Add Alaska and Hawaii to the full map
cbsa_trans_sf <- st_as_sf(cbsa_trans)
cbsa_trans_sf <- cbsa_trans_sf |> dplyr::filter(AK != 1) |> 
  dplyr::filter(HI != 1) |>
  dplyr::filter(PR != 1) |>
  dplyr::filter(AS != 1) |>
  dplyr::filter(GU != 1) |>
  dplyr::filter(MH != 1) |>
  dplyr::filter(FM != 1) |>
  dplyr::filter(MP != 1) |>
  dplyr::filter(PW != 1) |>
  dplyr::filter(VI != 1) 
alaska_sf <- st_as_sf(alaska)
hawaii_sf <- st_as_sf(hawaii)
cbsa_trans_sf <- rbind(cbsa_trans_sf, alaska_sf, hawaii_sf)

# Examine county-level data
ggplot(data = cbsa_trans_sf) +
  geom_sf() + theme_void()

####################
# County-level map #
####################

# Project data
county_trans <- st_transform(county, crs = st_crs("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"))
county_trans <- as(county_trans, 'Spatial')

# Alaska
alaska <- county_trans[county_trans$STATEFP=="02",]
alaska <- elide(alaska, rotate=-50)
alaska <- elide(alaska, scale=max(apply(bbox(alaska), 1, diff)) / 2.3)
alaska <- elide(alaska, shift=c(-2300000, -2500000))
proj4string(alaska) <- proj4string(county_trans)

# Hawaii
hawaii <- county_trans[county_trans$STATEFP=="15",]
hawaii <- elide(hawaii, rotate=-35)
hawaii <- elide(hawaii, shift=c(5400000, -1400000))
proj4string(hawaii) <- proj4string(county_trans)

# PR
pr <- county_trans[county_trans$STATEFP=="72",]
pr <- elide(pr, rotate=20)
pr <- elide(pr, scale=max(apply(bbox(pr), 1, diff)) *2)
pr <- elide(pr, shift=c(-2200000, 170000))
proj4string(pr) <- proj4string(county_trans)

# Add Alaska and Hawaii to the full map
county_trans <- county_trans[!county_trans$STATEFP %in% c("02", "15", "72"),]
county_trans <- rbind(county_trans, alaska, hawaii)

# Convert to sf object
county_trans_sf <- st_as_sf(county_trans)

# Examine county-level data
ggplot(data = county_trans_sf) +
  geom_sf() + theme_void()

####################
# State-level map #
####################

# Project data
state_trans <- st_transform(state, crs = st_crs("+proj=laea +lat_0=45 +lon_0=-100 +x_0=0 +y_0=0 +a=6370997 +b=6370997 +units=m +no_defs"))
state_trans <- as(state_trans, 'Spatial')

# Alaska
alaska <- state_trans[state_trans$STATEFP=="02",]
alaska <- elide(alaska, rotate=-50)
alaska <- elide(alaska, scale=max(apply(bbox(alaska), 1, diff)) / 2.3)
alaska <- elide(alaska, shift=c(-2300000, -2500000))
proj4string(alaska) <- proj4string(state_trans)

# Hawaii
hawaii <- state_trans[state_trans$STATEFP=="15",]
hawaii <- elide(hawaii, rotate=-35)
hawaii <- elide(hawaii, shift=c(5400000, -1400000))
proj4string(hawaii) <- proj4string(state_trans)

# PR
pr <- state_trans[state_trans$STATEFP=="72",]
pr <- elide(pr, rotate=20)
pr <- elide(pr, scale=max(apply(bbox(pr), 1, diff)) *2)
pr <- elide(pr, shift=c(-2200000, 170000))
proj4string(pr) <- proj4string(state_trans)

# Add Alaska and Hawaii to the full map
state_trans <- state_trans[!state_trans$STATEFP %in% c("02", "15", "72"),]
state_trans <- rbind(state_trans, alaska, hawaii)

# Convert to sf object
state_trans_sf <- st_as_sf(state_trans)

# Examine state-level data
ggplot(data = state_trans_sf) +
  geom_sf() + theme_void()

################
# 4. SAVE MAPS #
################

saveRDS(state_trans_sf, "./tmp/usa_albers_state.rds")
saveRDS(county_trans_sf, "./tmp/usa_albers_county.rds")
saveRDS(cbsa_trans_sf, "./tmp/usa_albers_cbsa.rds")

################################################################################
################################################################################
