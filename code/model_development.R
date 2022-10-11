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

# Load data ----
yoy_hake <- readRDS(here('data', 'yoy_hake.Rdata')) %>% 
#  tidyr::drop_na(temperature, salinity, bottom_depth) %>%
  tidyr::drop_na(bottom_depth) %>%
  filter(catch < 2500)
yoy_hake$year_f <- as.factor(yoy_hake$year)

ctds <- readRDS(here('data', 'ctd_means.rdata'))

ctd_means <- ctds %>%
  group_by(year) %>%
  summarise(across(temperature, mean, na.rm = TRUE))

# Functions ----
source(here('code/functions', 'vis_gam_COLORS.R'))
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
             legend.args = list("log(cpue+1)",
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

# Hake ----
# Aggregate model
hake_formula <- formula(y_catch_adj + 1 ~ s(lon, lat) + 
                          s(bottom_depth, k = 4) +
                          s(jday) + 
                          s(lon, lat, by = NPGO_pos))

# Use models selected during model exploration
hake_total <- gam(catch + 1 ~ factor(year) + 
                    s(lon, lat) + 
                    s(bottom_depth, k = 4) +
                    s(jday) + 
                    s(lon, lat, by = NPGO_pos),
                  family = tw(link = "log"),
                  method = 'REML',
                  data = yoy_hake)
summary(hake_total)

year_effect_hake <- predict(hake_total, 
                            type = "terms")[, 1]

yoy_hake$y_catch <- yoy_hake$catch + year_effect_hake
yoy_hake$y_catch_adj <- yoy_hake$y_catch + abs(min(yoy_hake$y_catch))

# Leave one group out cross validation
# Run GAMs with each year left out
# Leave out one year, run model on remaining data
hake_gams <- lapply(unique(yoy_hake$year_f), function(x) {
  output <- gam(hake_formula,
                family = tw(link = "log"),
                method = 'REML',
                data = yoy_hake[yoy_hake$year_f != x, ])
})

hake_data <- split(yoy_hake, yoy_hake$year)

# Get predictions
# Predict on the left out year's data
hake_results <- hake_data
for(i in seq_along(hake_gams)){
  for(j in seq_along(hake_data)){
  hake_results[[j]]$pred <- predict(hake_gams[[i]],
                               newdata = hake_results[[j]],
                               type = "response")
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

ggplot(hake_error_woy) +
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
                                s(lon, lat, by = NPGO_pos))

# Use models selected during model exploration
hake_small <- gam(small + 1 ~ factor(year) +
                    s(lon, lat) + 
                    s(bottom_depth, k = 4) +
                    s(jday) + 
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
                                s(lon, lat, by = NPGO_pos))

# Use models selected during model exploration
hake_large <- gam(large + 1 ~ factor(year) +
                    s(lon, lat) + 
                    s(bottom_depth, k = 4) +
                    s(jday) + 
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

# Maps
par(mfrow = c(1, 3),
    mar = c(6.4, 7.2, 2.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
location_plot(hake_total, yoy_hake, "lat", "All Sizes", yoy_hake$catch)
location_plot(hake_large, yoy_hake, " ", "Small Sizes", yoy_hake$upper_cpue)
location_plot(hake_small, yoy_hake, " ", "Large Sizes", yoy_hake$lower_cpue)
dev.copy(jpeg, here('results/RREAS_preliminary', 'hake_distributions.jpg'), 
         height = 15, width = 20, units = 'in', res = 200)
dev.off()

