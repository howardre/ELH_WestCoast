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
library(Metrics)
library(dplyr)

# Load data
yoy_hake <- readRDS(here('data', 'yoy_hake.rdata')) %>% 
  drop_na(temperature, salinity, bottom_depth)  %>% 
  filter(catch < 2500)

ctds <- readRDS(here('data', 'ctd_means.rdata'))

ctd_means <- ctds %>%
  group_by(year) %>%
  summarise(across(temperature, mean, na.rm = TRUE))

# Look at outliers
ggplot(yoy_hake) +
  geom_point(aes(catch, bottom_depth),
             alpha = 0.1,
             size = 2)

# Functions

# Hake ----
# Aggregate model
hake_formula <- formula(catch + 1 ~ s(year, bs = "re") + s(longitude, latitude) + s(bottom_depth, k = 4) +
                          s(julian) + s(temperature, k = 4) + s(salinity, k = 4) + s(longitude, latitude, by = NPGO_pos))

# Use models selected during model exploration
hake_total <- gam(hake_formula,
                  family = tw(link = "log"),
                  method = 'REML',
                  data = yoy_hake)
summary(hake_total)

# Leave one group out cross validation
# Run GAMs with each year left out
# Leave out one year, run model on remaining data
hake_gams <- lapply(unique(yoy_hake$year), function(x) {
  output <- gam(hake_formula,
                family = tw(link = "log"),
                method = 'REML',
                data = yoy_hake[yoy_hake$year != x, ])
})

hake_data <- split(yoy_hake, yoy_hake$year)

# Get predictions
# Predict on the left out year's data
hake_results <- hake_data
for(i in seq_along(hake_gams)){
  for(j in seq_along(hake_data)){
  hake_results[[j]]$pred <- predict(hake_gams[[i]],
                               newdata = hake_data[[j]],
                               type = "response",                                   
                               exclude = "s(year)")
  }}

# Calculate RMSE
# Get values for each year and overall value
hake_RMSE <- lapply(hake_results, function(x) {
  rmse(x$catch, x$pred)
})

range(yoy_hake$catch)
mean(unlist(hake_RMSE))

hake_error <- as.data.frame(do.call(rbind, hake_RMSE))
hake_error$year <- rownames(hake_error)
rownames(hake_error) <- NULL
colnames(hake_error)[1] <- "RMSE"

# Plot the RMSE for each year
hake_error$temperature <- ctd_means$temperature[match(hake_error$year, ctd_means$year)]

ggplot(hake_error) +
  geom_line(aes(year, RMSE),
             size = 1,
             group = 1) 

ggplot(hake_error) +
  geom_line(aes(year, temperature),
            size = 1,
            group = 1) 


# Size explicit
# Small
hake_small_formula <- formula(lower_cpue + 1 ~ s(year, bs = "re") + s(longitude, latitude) + s(bottom_depth, k = 4) +
                          s(julian) + s(temperature, k = 4) + s(salinity, k = 4) + s(longitude, latitude, by = NPGO_pos))

# Use models selected during model exploration
hake_small <- gam(hake_small_formula,
                  family = tw(link = "log"),
                  method = 'REML',
                  data = yoy_hake)
summary(hake_total)

# Leave one group out cross validation
hake_small_gams <- lapply(unique(yoy_hake$year), function(x) {
  output <- gam(hake_small_formula,
                family = tw(link = "log"),
                method = 'REML',
                data = yoy_hake[yoy_hake$year != x, ])
})

# Get predictions
# Predict on the left out year's data
hake_small_results <- hake_data
for(i in seq_along(hake_small_gams)){
  for(j in seq_along(hake_data)){
    hake_small_results[[j]]$pred_small <- predict(hake_small_gams[[i]],
                                                  newdata = hake_data[[j]],
                                                  type = "response",
                                                  exclude = "s(year)")
  }}

# Calculate RMSE
# Get values for each year and overall value
hake_small_RMSE <- lapply(hake_small_results, function(x) {
  rmse(x$lower_cpue, x$pred_small)
})

mean(unlist(hake_small_RMSE))

hake_small_error <- as.data.frame(do.call(rbind, hake_small_RMSE))
hake_small_error$year <- rownames(hake_small_error)
rownames(hake_small_error) <- NULL
colnames(hake_small_error)[1] <- "RMSE"

# Plot the RMSE for each year
ggplot(hake_small_error) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) # something weird going on in 2004

# Large
hake_large_formula <- formula(upper_cpue + 1 ~ s(year, bs = "re") + s(longitude, latitude) + s(bottom_depth, k = 4) +
                                s(julian) + s(temperature, k = 4) + s(salinity, k = 4) + s(longitude, latitude, by = NPGO_pos))

# Use models selected during model exploration
hake_large <- gam(hake_large_formula,
                  family = tw(link = "log"),
                  method = 'REML',
                  data = yoy_hake)
summary(hake_total)

# Leave one group out cross validation
hake_large_gams <- lapply(unique(yoy_hake$year), function(x) {
  output <- gam(hake_large_formula,
                family = tw(link = "log"),
                method = 'REML',
                data = yoy_hake[yoy_hake$year != x, ])
})

# Get predictions
# Predict on the left out year's data
hake_large_results <- hake_data
for(i in seq_along(hake_large_gams)){
  for(j in seq_along(hake_data)){
    hake_large_results[[j]]$pred_large <- predict(hake_large_gams[[i]],
                                                  newdata = hake_data[[j]],
                                                  type = "response",
                                                  exclude = "s(year)")
  }}

# Calculate RMSE
# Get values for each year and overall value
hake_large_RMSE <- lapply(hake_large_results, function(x) {
  rmse(x$upper_cpue, x$pred_large)
})

mean(unlist(hake_large_RMSE))

hake_large_error <- as.data.frame(do.call(rbind, hake_large_RMSE))
hake_large_error$year <- rownames(hake_large_error)
rownames(hake_large_error) <- NULL
colnames(hake_large_error)[1] <- "RMSE"

# Plot the RMSE for each year
ggplot(hake_large_error) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

# Add the predictions for small and large together
hake_combined_results <- map2(hake_small_results, 
                              hake_large_results,
                              ~ .x %>% bind_cols(.y %>% select(pred_large)))

hake_combined_preds <- lapply(hake_combined_results, function(x) {
  mutate(x, pred = pred_large + pred_small)    
})

hake_combined_RMSE <- lapply(hake_combined_preds, function(x) {
  rmse(x$catch, x$pred)
})

mean(unlist(hake_combined_RMSE))

hake_combined_error <- as.data.frame(do.call(rbind, hake_combined_RMSE))
hake_combined_error$year <- rownames(hake_combined_error)
rownames(hake_combined_error) <- NULL
colnames(hake_combined_error)[1] <- "RMSE"

# Plot the RMSE for each year
ggplot(hake_combined_error) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 