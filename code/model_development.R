### Title: Model Development
### Author: Rebecca Howard
### Date: 09/13/2022

# Libraries ----
library(mgcv)
library(maps)
library(mapdata)
library(here)
library(purrr)
library(ggplot2)
library(Metrics)
library(dplyr)
library(fields)
library(colorspace)

# Functions ----
source(here('code/functions', 'vis_gam_COLORS.R'))
read_data <- function(file){
  yoy <- readRDS(here('data', file)) %>% 
    tidyr::drop_na(temperature, salinity, bottom_depth) %>%
    tidyr::drop_na(bottom_depth) %>%
    filter(catch < 2500)
  yoy$lncpue <- log(yoy$catch + 1)
  return(yoy)
}

contour_col <- rgb(0, 0, 255, max = 255, alpha = 0, names = "white")
jet.colors <- colorRampPalette(c(sequential_hcl(15, palette = "Mint")))

location_plot <- function(gam, species_subset, yaxis, title, value) {
  myvis_gam(gam,
            view = c('lon', 'lat'),
            too.far = 0.07,
            plot.type = 'contour',
            contour.col = contour_col,
            color = "jet" ,
            type = 'link',
            xlim = c(-125.7, -116.5),
            ylim = range(species_subset$lat, na.rm = TRUE) + c(-.4, .5),
            family = "serif",
            xlab = "lon",
            ylab = yaxis,
            main = title,
            cex.lab = 2.5,
            cex.axis =  2.5,
            cex.main = 3)
  symbols(species_subset$lon,
          species_subset$lat,
          circle = value,
          inches = 0.2,
          add = T,
          bg = alpha('dimgray', 0.4),
          fg = alpha('black', 0.1))  
  maps::map('worldHires',
            add = T,
            col = 'antiquewhite4',
            fill = T)
  image.plot(legend.only = T,
             col = jet.colors(100),
             legend.shrink = 0.2,
             smallplot = c(.28, .31, .11, .24),
             legend.cex = 1.3,
             axis.args = list(cex.axis = 1.8,
                              family = "serif"),
             legend.width = 0.8,
             legend.mar = 6,
             zlim = c(min(gam$linear.predictors), 
                      max(gam$linear.predictors)),
             legend.args = list("CPUE",
                                side = 2, cex = 1.4,
                                family = "serif"))
}

variable_coefficient <- function(gam, data, variable){
  preds <- predict(gam, type = 'terms', se.fit = T)
  pred_slope <- preds[[1]][, 7] / variable
  pred_slope_se <- 1.96 * preds[[2]][, 7]
  pred_slope_up <- (preds[[1]][, 7] + pred_slope_se) / variable
  pred_slope_low <- (preds[[1]][, 7] - pred_slope_se) / variable
  sign_slope_pos <- (1:length(pred_slope))[pred_slope_low > 0]
  sign_slope_neg <- (1:length(pred_slope))[pred_slope_up < 0]
  return(list(sign_slope_neg, sign_slope_pos, pred_slope))
}

plot_var_coef <- function(my_gam, species_subset, predictions){
  par(mar = c(6.4, 7.2, .5, 0.6) + 0.1,
      oma = c(1, 1, 1, 1),
      mgp = c(5, 2, 0))
  myvis_gam(my_gam,
            view = c('lon', 'lat'),
            too.far = 0.07,
            plot.type = 'contour',
            contour.col = contour_col,
            color = "jet" ,
            type = 'link',
            xlim = c(-125.7,-116.5),
            ylim = range(species_subset$lat, na.rm = TRUE) + c(-.4, .5),
            family = "serif",
            xlab = "lon",
            ylab = "lat",
            main = " ",
            cex.lab = 2.5,
            cex.axis =  2.5)
  symbols(species_subset$lon[predictions[[2]]],
          species_subset$lat[predictions[[2]]],
          circle = predictions[[3]][predictions[[2]]],
          inches = 0.12,
          add = T,
          bg = alpha('darkred', 0.4),
          fg = alpha('black', 0.08))
  symbols(species_subset$lon[predictions[[1]]],
          species_subset$lat[predictions[[1]]],
          circle = (-1) * predictions[[3]][predictions[[1]]],
          inches = 0.12,
          add = T,
          bg = alpha('navy', 0.4),
          fg = alpha('black', 0.08))
  map("worldHires",
      fill = T,
      col = "wheat4",
      add = T)
  image.plot(legend.only = T,
             col = jet.colors(100),
             legend.shrink = 0.2,
             smallplot = c(.28, .31, .11, .24),
             legend.cex = 1.3,
             axis.args = list(cex.axis = 1.8,
                              family = "serif"),
             legend.width = 0.8,
             legend.mar = 6,
             zlim = c(min(my_gam$linear.predictors),
                      max(my_gam$linear.predictors)),
             legend.args = list("log(cpue+1)",
                                side = 2, cex = 2,
                                family = "serif"))
}

# Load data ----
yoy_hake <- read_data('yoy_hake.Rdata')
yoy_anchovy <- read_data('yoy_anch.Rdata')
yoy_widow <- read_data('yoy_widw.Rdata')
yoy_shortbelly <- read_data('yoy_sbly.Rdata')
yoy_sdab <- read_data('yoy_dab.Rdata')

ctds <- readRDS(here('data', 'ctd_means.rdata'))

ctd_means <- ctds %>%
  group_by(year) %>%
  summarise(across(temperature, mean, na.rm = TRUE))

# Hake ----
# Aggregate model
# Use models selected during model exploration
hake_total <- gam(lncpue ~ factor(year) + 
                    s(lon, lat) + 
                    s(bottom_depth, k = 4) +
                    s(jday) + 
                    s(temperature, k = 4) +
                    s(salinity, k = 4) +
                    s(lon, lat, by = NPGO_pos),
                  family = tw(link = "log"),
                  method = "REML",
                  data = yoy_hake)
summary(hake_total)

test_values <- predict(hake_total, type = "response")

year_effect_hake <- predict.gam(hake_total,
                                type = "terms")[, 1]

yoy_hake$y_catch <- yoy_hake$catch + year_effect_hake
yoy_hake$y_catch_adj <- yoy_hake$y_catch + abs(min(yoy_hake$y_catch)) # issue with negative values, doesn't work otherwise

# Leave one group out cross validation
# Run GAMs with each year left out
# Leave out one year, run model on remaining data
hake_formula <- formula(y_catch_adj + 1 ~ s(lon, lat) + 
                          s(bottom_depth, k = 4) +
                          s(jday) + 
                          s(temperature, k = 4) +
                          s(salinity, k = 4) +
                          s(lon, lat, by = NPGO_pos)) # Note no factor(year), added into response

hake_gams <- lapply(unique(yoy_hake$year), function(x) {
  output <- gam(hake_formula,
                family = tw(link = "log"),
                method = 'REML',
                data = yoy_hake[yoy_hake$year != x, ])
}) # Gives the list of GAMs with each year left out

# Create list with data from each year
hake_data <- split(yoy_hake, yoy_hake$year)

# Get predictions
# Predict on the left out year's data
hake_results <- hake_data
for(i in seq_along(hake_gams)){
  for(j in seq_along(hake_data)){
  hake_results[[j]]$pred <- predict(hake_gams[[i]],
                               newdata = hake_results[[j]],
                               type = "response")
  }} # Predicts on the left out year

# Calculate RMSE
# Get values for each year and overall value
hake_RMSE <- lapply(hake_results, function(x) {
  rmse(x$catch, x$pred)
})

range(yoy_hake$catch)
mean(unlist(hake_RMSE)) # 185

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
hake_small_formula <- formula(y_small_adj + 1 ~ s(lon, lat) + 
                                s(bottom_depth, k = 4) +
                                s(jday) + 
                                s(temperature, k = 4) +
                                s(salinity, k = 4) +
                                s(lon, lat, by = NPGO_pos))

# Use models selected during model exploration
hake_small <- gam(small + 1 ~ factor(year) +
                    s(lon, lat) + 
                    s(bottom_depth, k = 4) +
                    s(jday) + 
                    s(temperature, k = 4) +
                    s(salinity, k = 4) +
                    s(lon, lat, by = NPGO_pos),
                  family = tw(link = "log"),
                  method = 'REML',
                  data = yoy_hake)
summary(hake_small)

year_small_hake <- predict(hake_small, 
                            type = "terms")[, 1]

yoy_hake$y_small <- yoy_hake$small + year_small_hake
yoy_hake$y_small_adj <- yoy_hake$y_small + abs(min(yoy_hake$y_small))

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
                                                  newdata = hake_small_results[[j]],
                                                  type = "response")
  }}

# Calculate RMSE
# Get values for each year and overall value
hake_small_RMSE <- lapply(hake_small_results, function(x) {
  rmse(x$small, x$pred_small)
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
            group = 1) 

# Large
hake_large_formula <- formula(y_large_adj + 1 ~ s(lon, lat) + 
                                s(bottom_depth, k = 4) +
                                s(jday) + 
                                s(temperature, k = 4) +
                                s(salinity, k = 4) +
                                s(lon, lat, by = NPGO_pos))

# Use models selected during model exploration
hake_large <- gam(large + 1 ~ factor(year) +
                    s(lon, lat) + 
                    s(bottom_depth, k = 4) +
                    s(jday) + 
                    s(temperature, k = 4) +
                    s(salinity, k = 4) +
                    s(lon, lat, by = NPGO_pos),
                  family = tw(link = "log"),
                  method = 'REML',
                  data = yoy_hake)
summary(hake_large)

year_large_hake <- predict(hake_large, 
                           type = "terms")[, 1]

yoy_hake$y_large <- yoy_hake$large + year_large_hake
yoy_hake$y_large_adj <- yoy_hake$y_large + abs(min(yoy_hake$y_large))

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
                                                  newdata = hake_large_results[[j]],
                                                  type = "response")
  }}

# Calculate RMSE
# Get values for each year and overall value
hake_large_RMSE <- lapply(hake_large_results, function(x) {
  rmse(x$large, x$pred_large)
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

# Add the large and small together
hake_combined_results <- map2(hake_large_results, 
                              hake_small_results, 
                              ~left_join(.x, .y))
hake_added_results <- lapply(hake_combined_results, function(x){
  rmse(x$catch, x$pred_small + x$pred_large)
})

mean(unlist(hake_added_results)) # 182

# Maps
par(mfrow = c(1, 3),
    mar = c(6.4, 7.2, 2.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
location_plot(hake_total, yoy_hake, "lat", "All Sizes", yoy_hake$catch)
location_plot(hake_large, yoy_hake, " ", "Small Sizes", yoy_hake$large)
location_plot(hake_small, yoy_hake, " ", "Large Sizes", yoy_hake$small)
dev.copy(jpeg, here('results/RREAS_preliminary', 'hake_distributions.jpg'), 
         height = 15, width = 20, units = 'in', res = 200)
dev.off()


# Anchovy ----
# Aggregate model
anchovy_formula <- formula(y_catch_adj + 1 ~ s(lon, lat) + 
                          s(bottom_depth, k = 4) +
                          s(jday) + 
                          s(temperature, k = 4) +
                          s(salinity, k = 4) +
                          s(lon, lat, by = NPGO_pos))

# Use models selected during model exploration
anchovy_total <- gam(catch + 1 ~ factor(year) + 
                    s(lon, lat) + 
                    s(bottom_depth, k = 4) +
                    s(jday) + 
                    s(temperature, k = 4) +
                    s(salinity, k = 4) +
                    s(lon, lat, by = NPGO_pos),
                  family = tw(link = "log"),
                  method = 'REML',
                  data = yoy_anchovy)
summary(anchovy_total)

year_effect_anchovy <- predict(anchovy_total, 
                            type = "terms")[, 1]

yoy_anchovy$y_catch <- yoy_anchovy$catch + year_effect_anchovy
yoy_anchovy$y_catch_adj <- yoy_anchovy$y_catch + abs(min(yoy_anchovy$y_catch))

# Leave one group out cross validation
# Run GAMs with each year left out
# Leave out one year, run model on remaining data
anchovy_gams <- lapply(unique(yoy_anchovy$year_f), function(x) {
  output <- gam(anchovy_formula,
                family = tw(link = "log"),
                method = 'REML',
                data = yoy_anchovy[yoy_anchovy$year_f != x, ])
})

anchovy_data <- split(yoy_anchovy, yoy_anchovy$year)

# Get predictions
# Predict on the left out year's data
anchovy_results <- anchovy_data
for(i in seq_along(anchovy_gams)){
  for(j in seq_along(anchovy_data)){
    anchovy_results[[j]]$pred <- predict(anchovy_gams[[i]],
                                      newdata = anchovy_results[[j]],
                                      type = "response")
  }}

# Calculate RMSE
# Get values for each year and overall value
anchovy_RMSE <- lapply(anchovy_results, function(x) {
  rmse(x$catch, x$pred)
})

range(yoy_anchovy$catch)
mean(unlist(anchovy_RMSE)) # 370

anchovy_error <- as.data.frame(do.call(rbind, anchovy_RMSE))
anchovy_error$year <- rownames(anchovy_error)
rownames(anchovy_error) <- NULL
colnames(anchovy_error)[1] <- "RMSE"


# Plot the RMSE for each year
anchovy_error$temperature <- ctd_means$temperature[match(anchovy_error$year, ctd_means$year)]

ggplot(anchovy_error) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(anchovy_error) +
  geom_line(aes(year, temperature),
            size = 1,
            group = 1) 


# Size explicit
# Small
anchovy_small_formula <- formula(y_small_adj + 1 ~ s(lon, lat) + 
                                s(bottom_depth, k = 4) +
                                s(jday) + 
                                s(temperature, k = 4) +
                                s(salinity, k = 4) +
                                s(lon, lat, by = NPGO_pos))

# Use models selected during model exploration
anchovy_small <- gam(small + 1 ~ factor(year) +
                    s(lon, lat) + 
                    s(bottom_depth, k = 4) +
                    s(jday) + 
                    s(temperature, k = 4) +
                    s(salinity, k = 4) +
                    s(lon, lat, by = NPGO_pos),
                  family = tw(link = "log"),
                  method = 'REML',
                  data = yoy_anchovy)
summary(anchovy_small)

year_small_anchovy <- predict(anchovy_small, 
                           type = "terms")[, 1]

yoy_anchovy$y_small <- yoy_anchovy$small + year_small_anchovy
yoy_anchovy$y_small_adj <- yoy_anchovy$y_small + abs(min(yoy_anchovy$y_small))

# Leave one group out cross validation
anchovy_small_gams <- lapply(unique(yoy_anchovy$year), function(x) {
  output <- gam(anchovy_small_formula,
                family = tw(link = "log"),
                method = 'REML',
                data = yoy_anchovy[yoy_anchovy$year != x, ])
})

# Get predictions
# Predict on the left out year's data
anchovy_small_results <- anchovy_data
for(i in seq_along(anchovy_small_gams)){
  for(j in seq_along(anchovy_data)){
    anchovy_small_results[[j]]$pred_small <- predict(anchovy_small_gams[[i]],
                                                  newdata = anchovy_small_results[[j]],
                                                  type = "response")
  }}

# Calculate RMSE
# Get values for each year and overall value
anchovy_small_RMSE <- lapply(anchovy_small_results, function(x) {
  rmse(x$small, x$pred_small)
})

mean(unlist(anchovy_small_RMSE))

anchovy_small_error <- as.data.frame(do.call(rbind, anchovy_small_RMSE))
anchovy_small_error$year <- rownames(anchovy_small_error)
rownames(anchovy_small_error) <- NULL
colnames(anchovy_small_error)[1] <- "RMSE"

# Plot the RMSE for each year
ggplot(anchovy_small_error) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

# Large
anchovy_large_formula <- formula(y_large_adj + 1 ~ s(lon, lat) + 
                                s(bottom_depth, k = 4) +
                                s(jday) + 
                                s(temperature, k = 4) +
                                s(salinity, k = 4) +
                                s(lon, lat, by = NPGO_pos))

# Use models selected during model exploration
anchovy_large <- gam(large + 1 ~ factor(year) +
                    s(lon, lat) + 
                    s(bottom_depth, k = 4) +
                    s(jday) + 
                    s(temperature, k = 4) +
                    s(salinity, k = 4) +
                    s(lon, lat, by = NPGO_pos),
                  family = tw(link = "log"),
                  method = 'REML',
                  data = yoy_anchovy)
summary(anchovy_large)

year_large_anchovy <- predict(anchovy_large, 
                           type = "terms")[, 1]

yoy_anchovy$y_large <- yoy_anchovy$large + year_large_anchovy
yoy_anchovy$y_large_adj <- yoy_anchovy$y_large + abs(min(yoy_anchovy$y_large))

# Leave one group out cross validation
anchovy_large_gams <- lapply(unique(yoy_anchovy$year), function(x) {
  output <- gam(anchovy_large_formula,
                family = tw(link = "log"),
                method = 'REML',
                data = yoy_anchovy[yoy_anchovy$year != x, ])
})

# Get predictions
# Predict on the left out year's data
anchovy_large_results <- anchovy_data
for(i in seq_along(anchovy_large_gams)){
  for(j in seq_along(anchovy_data)){
    anchovy_large_results[[j]]$pred_large <- predict(anchovy_large_gams[[i]],
                                                  newdata = anchovy_large_results[[j]],
                                                  type = "response")
  }}

# Calculate RMSE
# Get values for each year and overall value
anchovy_large_RMSE <- lapply(anchovy_large_results, function(x) {
  rmse(x$large, x$pred_large)
})

mean(unlist(anchovy_large_RMSE))

anchovy_large_error <- as.data.frame(do.call(rbind, anchovy_large_RMSE))
anchovy_large_error$year <- rownames(anchovy_large_error)
rownames(anchovy_large_error) <- NULL
colnames(anchovy_large_error)[1] <- "RMSE"

# Plot the RMSE for each year
ggplot(anchovy_large_error) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

# Add the large and small together
anchovy_combined_results <- map2(anchovy_large_results, 
                              anchovy_small_results, 
                              ~left_join(.x, .y))
anchovy_added_results <- lapply(anchovy_combined_results, function(x){
  rmse(x$catch, x$pred_small + x$pred_large)
})

mean(unlist(anchovy_added_results)) # 554

# Widow Rockfish ----
# Aggregate model
widow_formula <- formula(y_catch_adj + 1 ~ s(lon, lat) + 
                             s(bottom_depth, k = 4) +
                             s(jday) + 
                             s(temperature, k = 4) +
                             s(salinity, k = 4) +
                             s(lon, lat, by = NPGO_pos))

# Use models selected during model exploration
widow_total <- gam(catch + 1 ~ factor(year) + 
                       s(lon, lat) + 
                       s(bottom_depth, k = 4) +
                       s(jday) + 
                       s(temperature, k = 4) +
                       s(salinity, k = 4) +
                       s(lon, lat, by = NPGO_pos),
                     family = tw(link = "log"),
                     method = 'REML',
                     data = yoy_widow)
summary(widow_total)

year_effect_widow <- predict(widow_total, 
                               type = "terms")[, 1]

yoy_widow$y_catch <- yoy_widow$catch + year_effect_widow
yoy_widow$y_catch_adj <- yoy_widow$y_catch + abs(min(yoy_widow$y_catch))

# Leave one group out cross validation
# Run GAMs with each year left out
# Leave out one year, run model on remaining data
widow_gams <- lapply(unique(yoy_widow$year_f), function(x) {
  output <- gam(widow_formula,
                family = tw(link = "log"),
                method = 'REML',
                data = yoy_widow[yoy_widow$year_f != x, ])
})

widow_data <- split(yoy_widow, yoy_widow$year)

# Get predictions
# Predict on the left out year's data
widow_results <- widow_data
for(i in seq_along(widow_gams)){
  for(j in seq_along(widow_data)){
    widow_results[[j]]$pred <- predict(widow_gams[[i]],
                                         newdata = widow_results[[j]],
                                         type = "response")
  }}

# Calculate RMSE
# Get values for each year and overall value
widow_RMSE <- lapply(widow_results, function(x) {
  rmse(x$catch, x$pred)
})

range(yoy_widow$catch)
mean(unlist(widow_RMSE)) # 23

widow_error <- as.data.frame(do.call(rbind, widow_RMSE))
widow_error$year <- rownames(widow_error)
rownames(widow_error) <- NULL
colnames(widow_error)[1] <- "RMSE"


# Plot the RMSE for each year
widow_error$temperature <- ctd_means$temperature[match(widow_error$year, ctd_means$year)]

ggplot(widow_error) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(widow_error) +
  geom_line(aes(year, temperature),
            size = 1,
            group = 1) 


# Size explicit
# Small
widow_small_formula <- formula(y_small_adj + 1 ~ s(lon, lat) + 
                                   s(bottom_depth, k = 4) +
                                   s(jday) + 
                                   s(temperature, k = 4) +
                                   s(salinity, k = 4) +
                                   s(lon, lat, by = NPGO_pos))

# Use models selected during model exploration
widow_small <- gam(small + 1 ~ factor(year) +
                       s(lon, lat) + 
                       s(bottom_depth, k = 4) +
                       s(jday) + 
                       s(temperature, k = 4) +
                       s(salinity, k = 4) +
                       s(lon, lat, by = NPGO_pos),
                     family = tw(link = "log"),
                     method = 'REML',
                     data = yoy_widow)
summary(widow_small)

year_small_widow <- predict(widow_small, 
                              type = "terms")[, 1]

yoy_widow$y_small <- yoy_widow$small + year_small_widow
yoy_widow$y_small_adj <- yoy_widow$y_small + abs(min(yoy_widow$y_small))

# Leave one group out cross validation
widow_small_gams <- lapply(unique(yoy_widow$year), function(x) {
  output <- gam(widow_small_formula,
                family = tw(link = "log"),
                method = 'REML',
                data = yoy_widow[yoy_widow$year != x, ])
})

# Get predictions
# Predict on the left out year's data
widow_small_results <- widow_data
for(i in seq_along(widow_small_gams)){
  for(j in seq_along(widow_data)){
    widow_small_results[[j]]$pred_small <- predict(widow_small_gams[[i]],
                                                     newdata = widow_small_results[[j]],
                                                     type = "response")
  }}

# Calculate RMSE
# Get values for each year and overall value
widow_small_RMSE <- lapply(widow_small_results, function(x) {
  rmse(x$small, x$pred_small)
})

mean(unlist(widow_small_RMSE))

widow_small_error <- as.data.frame(do.call(rbind, widow_small_RMSE))
widow_small_error$year <- rownames(widow_small_error)
rownames(widow_small_error) <- NULL
colnames(widow_small_error)[1] <- "RMSE"

# Plot the RMSE for each year
ggplot(widow_small_error) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

# Large
widow_large_formula <- formula(y_large_adj + 1 ~ s(lon, lat) + 
                                   s(bottom_depth, k = 4) +
                                   s(jday) + 
                                   s(temperature, k = 4) +
                                   s(salinity, k = 4) +
                                   s(lon, lat, by = NPGO_pos))

# Use models selected during model exploration
widow_large <- gam(large + 1 ~ factor(year) +
                       s(lon, lat) + 
                       s(bottom_depth, k = 4) +
                       s(jday) + 
                       s(temperature, k = 4) +
                       s(salinity, k = 4) +
                       s(lon, lat, by = NPGO_pos),
                     family = tw(link = "log"),
                     method = 'REML',
                     data = yoy_widow)
summary(widow_large)

year_large_widow <- predict(widow_large, 
                              type = "terms")[, 1]

yoy_widow$y_large <- yoy_widow$large + year_large_widow
yoy_widow$y_large_adj <- yoy_widow$y_large + abs(min(yoy_widow$y_large))

# Leave one group out cross validation
widow_large_gams <- lapply(unique(yoy_widow$year), function(x) {
  output <- gam(widow_large_formula,
                family = tw(link = "log"),
                method = 'REML',
                data = yoy_widow[yoy_widow$year != x, ])
})

# Get predictions
# Predict on the left out year's data
widow_large_results <- widow_data
for(i in seq_along(widow_large_gams)){
  for(j in seq_along(widow_data)){
    widow_large_results[[j]]$pred_large <- predict(widow_large_gams[[i]],
                                                     newdata = widow_large_results[[j]],
                                                     type = "response")
  }}

# Calculate RMSE
# Get values for each year and overall value
widow_large_RMSE <- lapply(widow_large_results, function(x) {
  rmse(x$large, x$pred_large)
})

mean(unlist(widow_large_RMSE))

widow_large_error <- as.data.frame(do.call(rbind, widow_large_RMSE))
widow_large_error$year <- rownames(widow_large_error)
rownames(widow_large_error) <- NULL
colnames(widow_large_error)[1] <- "RMSE"

# Plot the RMSE for each year
ggplot(widow_large_error) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

# Add the large and small together
widow_combined_results <- map2(widow_large_results, 
                                 widow_small_results, 
                                 ~left_join(.x, .y))
widow_added_results <- lapply(widow_combined_results, function(x){
  rmse(x$catch, x$pred_small + x$pred_large)
})

mean(unlist(widow_added_results)) # 22

# Shortbelly Rockfish ----
# Aggregate model
shortbelly_formula <- formula(y_catch_adj + 1 ~ s(lon, lat) + 
                           s(bottom_depth, k = 4) +
                           s(jday) + 
                           s(temperature, k = 4) +
                           s(salinity, k = 4) +
                           s(lon, lat, by = NPGO_pos))

# Use models selected during model exploration
shortbelly_total <- gam(catch + 1 ~ factor(year) + 
                     s(lon, lat) + 
                     s(bottom_depth, k = 4) +
                     s(jday) + 
                     s(temperature, k = 4) +
                     s(salinity, k = 4) +
                     s(lon, lat, by = NPGO_pos),
                   family = tw(link = "log"),
                   method = 'REML',
                   data = yoy_shortbelly)
summary(shortbelly_total)

year_effect_shortbelly <- predict(shortbelly_total, 
                             type = "terms")[, 1]

yoy_shortbelly$y_catch <- yoy_shortbelly$catch + year_effect_shortbelly
yoy_shortbelly$y_catch_adj <- yoy_shortbelly$y_catch + abs(min(yoy_shortbelly$y_catch))

# Leave one group out cross validation
# Run GAMs with each year left out
# Leave out one year, run model on remaining data
shortbelly_gams <- lapply(unique(yoy_shortbelly$year_f), function(x) {
  output <- gam(shortbelly_formula,
                family = tw(link = "log"),
                method = 'REML',
                data = yoy_shortbelly[yoy_shortbelly$year_f != x, ])
})

shortbelly_data <- split(yoy_shortbelly, yoy_shortbelly$year)

# Get predictions
# Predict on the left out year's data
shortbelly_results <- shortbelly_data
for(i in seq_along(shortbelly_gams)){
  for(j in seq_along(shortbelly_data)){
    shortbelly_results[[j]]$pred <- predict(shortbelly_gams[[i]],
                                       newdata = shortbelly_results[[j]],
                                       type = "response")
  }}

# Calculate RMSE
# Get values for each year and overall value
shortbelly_RMSE <- lapply(shortbelly_results, function(x) {
  rmse(x$catch, x$pred)
})

range(yoy_shortbelly$catch)
mean(unlist(shortbelly_RMSE)) # 114

shortbelly_error <- as.data.frame(do.call(rbind, shortbelly_RMSE))
shortbelly_error$year <- rownames(shortbelly_error)
rownames(shortbelly_error) <- NULL
colnames(shortbelly_error)[1] <- "RMSE"


# Plot the RMSE for each year
shortbelly_error$temperature <- ctd_means$temperature[match(shortbelly_error$year, ctd_means$year)]

ggplot(shortbelly_error) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(shortbelly_error) +
  geom_line(aes(year, temperature),
            size = 1,
            group = 1) 


# Size explicit
# Small
shortbelly_small_formula <- formula(y_small_adj + 1 ~ s(lon, lat) + 
                                 s(bottom_depth, k = 4) +
                                 s(jday) + 
                                 s(temperature, k = 4) +
                                 s(salinity, k = 4) +
                                 s(lon, lat, by = NPGO_pos))

# Use models selected during model exploration
shortbelly_small <- gam(small + 1 ~ factor(year) +
                     s(lon, lat) + 
                     s(bottom_depth, k = 4) +
                     s(jday) + 
                     s(temperature, k = 4) +
                     s(salinity, k = 4) +
                     s(lon, lat, by = NPGO_pos),
                   family = tw(link = "log"),
                   method = 'REML',
                   data = yoy_shortbelly)
summary(shortbelly_small)

year_small_shortbelly <- predict(shortbelly_small, 
                            type = "terms")[, 1]

yoy_shortbelly$y_small <- yoy_shortbelly$small + year_small_shortbelly
yoy_shortbelly$y_small_adj <- yoy_shortbelly$y_small + abs(min(yoy_shortbelly$y_small))

# Leave one group out cross validation
shortbelly_small_gams <- lapply(unique(yoy_shortbelly$year), function(x) {
  output <- gam(shortbelly_small_formula,
                family = tw(link = "log"),
                method = 'REML',
                data = yoy_shortbelly[yoy_shortbelly$year != x, ])
})

# Get predictions
# Predict on the left out year's data
shortbelly_small_results <- shortbelly_data
for(i in seq_along(shortbelly_small_gams)){
  for(j in seq_along(shortbelly_data)){
    shortbelly_small_results[[j]]$pred_small <- predict(shortbelly_small_gams[[i]],
                                                   newdata = shortbelly_small_results[[j]],
                                                   type = "response")
  }}

# Calculate RMSE
# Get values for each year and overall value
shortbelly_small_RMSE <- lapply(shortbelly_small_results, function(x) {
  rmse(x$small, x$pred_small)
})

mean(unlist(shortbelly_small_RMSE))

shortbelly_small_error <- as.data.frame(do.call(rbind, shortbelly_small_RMSE))
shortbelly_small_error$year <- rownames(shortbelly_small_error)
rownames(shortbelly_small_error) <- NULL
colnames(shortbelly_small_error)[1] <- "RMSE"

# Plot the RMSE for each year
ggplot(shortbelly_small_error) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

# Large
shortbelly_large_formula <- formula(y_large_adj + 1 ~ s(lon, lat) + 
                                 s(bottom_depth, k = 4) +
                                 s(jday) + 
                                 s(temperature, k = 4) +
                                 s(salinity, k = 4) +
                                 s(lon, lat, by = NPGO_pos))

# Use models selected during model exploration
shortbelly_large <- gam(large + 1 ~ factor(year) +
                     s(lon, lat) + 
                     s(bottom_depth, k = 4) +
                     s(jday) + 
                     s(temperature, k = 4) +
                     s(salinity, k = 4) +
                     s(lon, lat, by = NPGO_pos),
                   family = tw(link = "log"),
                   method = 'REML',
                   data = yoy_shortbelly)
summary(shortbelly_large)

year_large_shortbelly <- predict(shortbelly_large, 
                            type = "terms")[, 1]

yoy_shortbelly$y_large <- yoy_shortbelly$large + year_large_shortbelly
yoy_shortbelly$y_large_adj <- yoy_shortbelly$y_large + abs(min(yoy_shortbelly$y_large))

# Leave one group out cross validation
shortbelly_large_gams <- lapply(unique(yoy_shortbelly$year), function(x) {
  output <- gam(shortbelly_large_formula,
                family = tw(link = "log"),
                method = 'REML',
                data = yoy_shortbelly[yoy_shortbelly$year != x, ])
})

# Get predictions
# Predict on the left out year's data
shortbelly_large_results <- shortbelly_data
for(i in seq_along(shortbelly_large_gams)){
  for(j in seq_along(shortbelly_data)){
    shortbelly_large_results[[j]]$pred_large <- predict(shortbelly_large_gams[[i]],
                                                   newdata = shortbelly_large_results[[j]],
                                                   type = "response")
  }}

# Calculate RMSE
# Get values for each year and overall value
shortbelly_large_RMSE <- lapply(shortbelly_large_results, function(x) {
  rmse(x$large, x$pred_large)
})

mean(unlist(shortbelly_large_RMSE))

shortbelly_large_error <- as.data.frame(do.call(rbind, shortbelly_large_RMSE))
shortbelly_large_error$year <- rownames(shortbelly_large_error)
rownames(shortbelly_large_error) <- NULL
colnames(shortbelly_large_error)[1] <- "RMSE"

# Plot the RMSE for each year
ggplot(shortbelly_large_error) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

# Add the large and small together
shortbelly_combined_results <- map2(shortbelly_large_results, 
                               shortbelly_small_results, 
                               ~left_join(.x, .y))
shortbelly_added_results <- lapply(shortbelly_combined_results, function(x){
  rmse(x$catch, x$pred_small + x$pred_large)
})

mean(unlist(shortbelly_added_results)) # 106

# Pacific Sanddab ----
# Aggregate model
sdab_formula <- formula(y_catch_adj + 1 ~ s(lon, lat) + 
                                s(bottom_depth, k = 4) +
                                s(jday) + 
                                s(temperature, k = 4) +
                                s(salinity, k = 4) +
                                s(lon, lat, by = NPGO_pos))

# Use models selected during model exploration
sdab_total <- gam(catch + 1 ~ factor(year) + 
                          s(lon, lat) + 
                          s(bottom_depth, k = 4) +
                          s(jday) + 
                          s(temperature, k = 4) +
                          s(salinity, k = 4) +
                          s(lon, lat, by = NPGO_pos),
                        family = tw(link = "log"),
                        method = 'REML',
                        data = yoy_sdab)
summary(sdab_total)

year_effect_sdab <- predict(sdab_total, 
                                  type = "terms")[, 1]

yoy_sdab$y_catch <- yoy_sdab$catch + year_effect_sdab
yoy_sdab$y_catch_adj <- yoy_sdab$y_catch + abs(min(yoy_sdab$y_catch))

# Leave one group out cross validation
# Run GAMs with each year left out
# Leave out one year, run model on remaining data
sdab_gams <- lapply(unique(yoy_sdab$year_f), function(x) {
  output <- gam(sdab_formula,
                family = tw(link = "log"),
                method = 'REML',
                data = yoy_sdab[yoy_sdab$year_f != x, ])
})

sdab_data <- split(yoy_sdab, yoy_sdab$year)

# Get predictions
# Predict on the left out year's data
sdab_results <- sdab_data
for(i in seq_along(sdab_gams)){
  for(j in seq_along(sdab_data)){
    sdab_results[[j]]$pred <- predict(sdab_gams[[i]],
                                            newdata = sdab_results[[j]],
                                            type = "response")
  }}

# Calculate RMSE
# Get values for each year and overall value
sdab_RMSE <- lapply(sdab_results, function(x) {
  rmse(x$catch, x$pred)
})

range(yoy_sdab$catch)
mean(unlist(sdab_RMSE)) # 64

sdab_error <- as.data.frame(do.call(rbind, sdab_RMSE))
sdab_error$year <- rownames(sdab_error)
rownames(sdab_error) <- NULL
colnames(sdab_error)[1] <- "RMSE"


# Plot the RMSE for each year
sdab_error$temperature <- ctd_means$temperature[match(sdab_error$year, ctd_means$year)]

ggplot(sdab_error) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(sdab_error) +
  geom_line(aes(year, temperature),
            size = 1,
            group = 1) 


# Size explicit
# Small
sdab_small_formula <- formula(y_small_adj + 1 ~ s(lon, lat) + 
                                      s(bottom_depth, k = 4) +
                                      s(jday) + 
                                      s(temperature, k = 4) +
                                      s(salinity, k = 4) +
                                      s(lon, lat, by = NPGO_pos))

# Use models selected during model exploration
sdab_small <- gam(small + 1 ~ factor(year) +
                          s(lon, lat) + 
                          s(bottom_depth, k = 4) +
                          s(jday) + 
                          s(temperature, k = 4) +
                          s(salinity, k = 4) +
                          s(lon, lat, by = NPGO_pos),
                        family = tw(link = "log"),
                        method = 'REML',
                        data = yoy_sdab)
summary(sdab_small)

year_small_sdab <- predict(sdab_small, 
                                 type = "terms")[, 1]

yoy_sdab$y_small <- yoy_sdab$small + year_small_sdab
yoy_sdab$y_small_adj <- yoy_sdab$y_small + abs(min(yoy_sdab$y_small))

# Leave one group out cross validation
sdab_small_gams <- lapply(unique(yoy_sdab$year), function(x) {
  output <- gam(sdab_small_formula,
                family = tw(link = "log"),
                method = 'REML',
                data = yoy_sdab[yoy_sdab$year != x, ])
})

# Get predictions
# Predict on the left out year's data
sdab_small_results <- sdab_data
for(i in seq_along(sdab_small_gams)){
  for(j in seq_along(sdab_data)){
    sdab_small_results[[j]]$pred_small <- predict(sdab_small_gams[[i]],
                                                        newdata = sdab_small_results[[j]],
                                                        type = "response")
  }}

# Calculate RMSE
# Get values for each year and overall value
sdab_small_RMSE <- lapply(sdab_small_results, function(x) {
  rmse(x$small, x$pred_small)
})

mean(unlist(sdab_small_RMSE))

sdab_small_error <- as.data.frame(do.call(rbind, sdab_small_RMSE))
sdab_small_error$year <- rownames(sdab_small_error)
rownames(sdab_small_error) <- NULL
colnames(sdab_small_error)[1] <- "RMSE"

# Plot the RMSE for each year
ggplot(sdab_small_error) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

# Large
sdab_large_formula <- formula(y_large_adj + 1 ~ s(lon, lat) + 
                                      s(bottom_depth, k = 4) +
                                      s(jday) + 
                                      s(temperature, k = 4) +
                                      s(salinity, k = 4) +
                                      s(lon, lat, by = NPGO_pos))

# Use models selected during model exploration
sdab_large <- gam(large + 1 ~ factor(year) +
                          s(lon, lat) + 
                          s(bottom_depth, k = 4) +
                          s(jday) + 
                          s(temperature, k = 4) +
                          s(salinity, k = 4) +
                          s(lon, lat, by = NPGO_pos),
                        family = tw(link = "log"),
                        method = 'REML',
                        data = yoy_sdab)
summary(sdab_large)

year_large_sdab <- predict(sdab_large, 
                                 type = "terms")[, 1]

yoy_sdab$y_large <- yoy_sdab$large + year_large_sdab
yoy_sdab$y_large_adj <- yoy_sdab$y_large + abs(min(yoy_sdab$y_large))

# Leave one group out cross validation
sdab_large_gams <- lapply(unique(yoy_sdab$year), function(x) {
  output <- gam(sdab_large_formula,
                family = tw(link = "log"),
                method = 'REML',
                data = yoy_sdab[yoy_sdab$year != x, ])
})

# Get predictions
# Predict on the left out year's data
sdab_large_results <- sdab_data
for(i in seq_along(sdab_large_gams)){
  for(j in seq_along(sdab_data)){
    sdab_large_results[[j]]$pred_large <- predict(sdab_large_gams[[i]],
                                                        newdata = sdab_large_results[[j]],
                                                        type = "response")
  }}

# Calculate RMSE
# Get values for each year and overall value
sdab_large_RMSE <- lapply(sdab_large_results, function(x) {
  rmse(x$large, x$pred_large)
})

mean(unlist(sdab_large_RMSE))

sdab_large_error <- as.data.frame(do.call(rbind, sdab_large_RMSE))
sdab_large_error$year <- rownames(sdab_large_error)
rownames(sdab_large_error) <- NULL
colnames(sdab_large_error)[1] <- "RMSE"

# Plot the RMSE for each year
ggplot(sdab_large_error) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

# Add the large and small together
sdab_combined_results <- map2(sdab_large_results, 
                                    sdab_small_results, 
                                    ~left_join(.x, .y))
sdab_added_results <- lapply(sdab_combined_results, function(x){
  rmse(x$catch, x$pred_small + x$pred_large)
})

mean(unlist(sdab_added_results)) # 54