####################
# 1. SET-UP SCRIPT #
####################

# Start with a clear environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(lubridate)
library(rstan)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

###############
# 2. ANALYSIS #
###############

# Load data
flu <- readRDS('./tmp/flu_p1_dat_proc.rds')

# Restrict data to a single county
flu_la <- flu |> filter(county_fips == '06037') |> uncount(patient_count_imp)

# Restrict data to flu cases (for conditional probability)
flu_la_flu <- flu_la |> filter(flu == 1)

# Calculate observed probabilities
N <- length(flu_la$flu)
p.flu <- sum(flu_la$flu) / N
p.fever <- sum(flu_la$fever) / N
p.fever.flu <- sum(flu_la_flu$fever) / length(flu_la_flu$flu)

(p.fever.flu*p.flu) / p.fever


# Create data
example_dat <- list(N = N, 
                    p_flu = p.flu,
                    p_fever = p.fever,
                    p_fever_flu = p.fever.flu)

# Set stan options
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Fit stan model
fit <- stan(file = 'syndromic-surveillance/basic_model.stan', data = example_dat,
            iter = 10000)



pairs(fit)

print(fit)
plot(fit)
pairs(fit)
