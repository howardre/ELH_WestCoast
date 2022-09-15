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

# Load data
yoy_hake <- readRDS(here('data', 'yoy_hake.rdata')) %>% 
  drop_na(temperature, salinity, bottom_depth) %>%
  filter(catch < 711)

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

hake_cv <- lapply(hake_gams, function(x) {
  exp(predict(x, 
              newdata = hake_data, 
              type = "link"))
})
  
  
  
# Test prediction method (getting very large predictions that don't make sense)
yoy_hake <- dplyr::filter(yoy_hake, catch < 711) # seems like the issue is with outliers?
wo2000_gam <- gam(hake_formula,
                  family = tw(link = "log"),
                  method = "REML",
                  data = yoy_hake[yoy_hake$year != 2000, ])

data_2000 <- dplyr::filter(yoy_hake, year == 2000)

data_2000$pred <- predict(wo2000_gam,
                          newdata = data_2000,
                          type = "response",
                          exclude = "s(year)")

sqrt(mean(data_2000$catch - data_2000$pred, na.rm = T)^2) 

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

