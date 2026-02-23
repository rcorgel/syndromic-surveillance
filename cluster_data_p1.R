################################################################################
# File Name: cluster_data_p1                                                   #
#                                                                              #
# Purpose:   Cluster Influenza, COVID-19, and RSV medical claims data.         #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Cluster Influenza data                                         #
#            3. Cluster COVID-19 data                                          #
#            3. Cluster RSV data                                               #
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
library(dbscan)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/')

#######################
# 2. CLUSTER FLU DATA #
#######################

# Load data
flu <- readRDS('./tmp/flu_p1_dat_proc.rds')

# Collapse data to symptoms and disease
flu_cluster <- flu |> 
  # Keep only flu cases and drop asymptomatic cases 
  filter(symp_count != 0) |>
  filter(flu == 1) |>
  group_by(flu, fever, myalgia, cough, sore_throat, 
           short_breath, hypoxemia, chest_pain, bronchitis,
           nausea_vom, diarrhea, fatigue, headache,
           congestion, sneezing) |>
  mutate(weight = sum(patient_count_imp)) |>
  distinct(flu, fever, myalgia, cough, sore_throat, 
           short_breath, hypoxemia, chest_pain, bronchitis,
           nausea_vom, diarrhea, fatigue, headache,
           congestion, sneezing, weight)

clusters <- hdbscan(as.matrix(flu_cluster[, c(2:15)]), minPts = 5)

# Conduct cluster analysis
clusters <- hclust(dist(flu_cluster[, c(2:15)]), members = flu_cluster$weight)
plot(clusters)

# Determine cuts
clusterCut <- cutree(clusters, 8)

# Merge on cluster
flu_cluster$cluster <- clusterCut

# Summarize clusters
cluster <- flu_cluster |> 
  group_by(cluster) |>
  mutate(flu = sum(flu*weight) / sum(weight),
         fever = sum(fever*weight) / sum(weight),
         myalgia = sum(myalgia*weight) / sum(weight),
         cough = sum(cough*weight) / sum(weight),
         sore_throat = sum(sore_throat*weight) / sum(weight), 
         short_breath = sum(short_breath*weight) / sum(weight), 
         hypoxemia = sum(hypoxemia*weight) / sum(weight), 
         chest_pain = sum(chest_pain*weight) / sum(weight), 
         bronchitis = sum(bronchitis*weight) / sum(weight),
         nausea_vom = sum(nausea_vom*weight) / sum(weight), 
         diarrhea = sum(diarrhea*weight) / sum(weight), 
         fatigue = sum(fatigue*weight) / sum(weight), 
         headache = sum(headache*weight) / sum(weight),
         congestion = sum(congestion*weight) / sum(weight), 
         sneezing = sum(sneezing*weight) / sum(weight),
         group_size = sum(weight)) |>
  distinct(flu, fever, myalgia, cough, sore_throat, 
           short_breath, hypoxemia, chest_pain, bronchitis,
           nausea_vom, diarrhea, fatigue, headache,
           congestion, sneezing, cluster, group_size)

library(factoextra)
library(NbClust)
library(cluster)

# Elbow method
fviz_nbclust(flu_cluster[, c(2:15)], pam, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

#########################
# 3. CLUSTER COVID DATA #
#########################


#######################
# 4. CLUSTER RSV DATA #
#######################
