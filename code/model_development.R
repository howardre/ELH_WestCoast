### Title: Model Development
### Author: Rebecca Howard
### Date: 09/13/2022

# Libraries
library(mgcv)
library(maps)
library(mapdata)
library(here)
library(tidyr)
library(purrr)
library(ggplot2)

# Load data
yoy_hake <- readRDS(here('data', 'yoy_hake.rdata')) %>% 
  drop_na(temperature, salinity, bottom_depth)  %>% 
  dplyr::filter(catch < 2500)

ggplot(yoy_hake) +
  geom_point(aes(catch, bottom_depth),
             alpha = 0.1,
             size = 2)

# Functions

# Hake ----
# Leave one group out cross validation
hake_formula <- formula(catch + 1 ~ s(year, bs = "re") + s(longitude, latitude) + s(bottom_depth, k = 4) +
                       s(julian) + s(temperature, k = 4) + s(salinity, k = 4) + s(longitude, latitude, by = NPGO_pos))

# Run GAMs with each year left out
hake_gams <- lapply(unique(yoy_hake$year), function(x) {
  output <- gam(hake_formula,
                family = tw(link = "log"),
                method = 'REML',
                data = yoy_hake[yoy_hake$year != x, ])
})

hake_data <- split(yoy_hake, yoy_hake$year)

# Get predictions
hake_results <- hake_data
for(i in seq_along(hake_gams)){
  for(j in seq_along(hake_data)){
  hake_results[[j]]$pred <- predict(hake_gams[[i]],
                               newdata = hake_data[[j]],
                               type = "response",                                   
                               exclude = "s(year)")
  }}


# Use models selected during model exploration
hake_total <- gam(hake_formula,
                  family = tw(link = "log"),
                  method = 'REML',
                  data = yoy_hake)
hake_lower <- gam(hake_formula,
                  family = tw(link = "log"),
                  method = 'REML',
                  data = yoy_hake)
hake_upper <- gam(hake_formula,
                  family = tw(link = "log"),
                  method = 'REML',
                  data = yoy_hake)
summary(hake_total)
summary(hake_lower)
summary(hake_upper)

