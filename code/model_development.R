### Title: Model Development
### Author: Rebecca Howard
### Date: 09/13/2022

# Libraries
library(mgcv)
library(maps)
library(mapdata)
library(here)

# Load data
yoy_hake <- readRDS(here('data', 'yoy_hake.rdata'))

# Functions

# Hake ----
# Use models selected during model exploration
hake_total <- gam(catch + 1 ~ factor(year) + s(longitude, latitude) + s(bottom_depth, k = 4) +
                    s(julian) + s(temperature, k = 4) + s(salinity, k = 4) + s(longitude, latitude, by = NPGO_pos),
                  family = tw(link = "log"),
                  method = 'REML',
                  data = yoy_hake)
hake_lower <- gam(lower_cpue + 1 ~ factor(year) + s(longitude, latitude) + s(bottom_depth, k = 4) + 
                    s(julian) + s(temperature, k = 4) + s(salinity, k = 4) + s(longitude, latitude, by = NPGO_pos),
                  family = tw(link = "log"),
                  method = 'REML',
                  data = yoy_hake)
hake_upper <- gam(upper_cpue + 1 ~ factor(year) + s(longitude, latitude) + s(bottom_depth, k = 4) + 
                    s(julian) + s(temperature, k = 4) + s(salinity, k = 4) + s(longitude, latitude, by = NPGO_pos),
                  family = tw(link = "log"),
                  method = 'REML',
                  data = yoy_hake)
summary(hake_total)
summary(hake_lower)
summary(hake_upper)

