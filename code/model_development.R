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
library(ggthemes)

# Functions ----
source(here('code/functions', 'vis_gam_COLORS.R'))
read_data <- function(file){
  yoy <- readRDS(here('data', file)) %>% 
    tidyr::drop_na(temperature, salinity, bottom_depth) %>%
    tidyr::drop_na(bottom_depth) %>%
    filter(catch < 2500)
  yoy$year_f <- as.factor(yoy$year)
  yoy$catch1 <- yoy$catch + 1
  yoy$small_catch1 <- yoy$small + 1
  yoy$large_catch1 <- yoy$large + 1
  return(yoy)
}

contour_col <- rgb(0, 0, 255, max = 255, alpha = 0, names = "white")
jet.colors <- colorRampPalette(c(sequential_hcl(15, palette = "BurgYl")))

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

plot_var_coef <- function(my_gam, species_subset, predictions, yaxis, size){
  myvis_gam(my_gam,
            view = c('lon', 'lat'),
            too.far = 0.07,
            plot.type = 'contour',
            contour.col = contour_col,
            color = "jet" ,
            type = 'link',
            xlim = c(-125.7, -116.5),
            ylim = range(species_subset$lat, na.rm = TRUE) + c(-.4, .5),
            family = "Lato",
            xlab = "Longitude",
            ylab = yaxis,
            main = size,
            cex.lab = 7.5,
            cex.axis =  7,
            cex.main = 8)
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
  maps::map("state",
            boundary = FALSE,
            fill = TRUE,
            col = "wheat4",
            add = TRUE)
  text(x = state_labels$lon, 
       y = state_labels$lat,
       state_labels$name, 
       pos = 1,
       col = "black",
       cex = 6,
       family = "Lato")
  image.plot(legend.only = T,
             col = jet.colors(100),
             legend.shrink = 0.2,
             smallplot = c(.24, .29, .08, .21),
             legend.cex = 4,
             axis.args = list(cex.axis = 5,
                              family = "Lato"),
             legend.width = 0.8,
             legend.mar = 6,
             zlim = c(min(my_gam$linear.predictors),
                      max(my_gam$linear.predictors)),
             legend.args = list("log(cpue+1)",
                                side = 2, 
                                cex = 4,
                                family = "Lato",
                                line = 1.5))
}

plot_var_coef2 <- function(my_gam, species_subset, predictions, yaxis, size){
  myvis_gam(my_gam,
            view = c('lon', 'lat'),
            too.far = 0.07,
            plot.type = 'contour',
            contour.col = contour_col,
            color = "jet" ,
            type = 'link',
            xlim = c(-125.7, -116.5),
            ylim = range(species_subset$lat, na.rm = TRUE) + c(-.4, .5),
            family = "Lato",
            xlab = "Longitude",
            ylab = yaxis,
            main = size,
            cex.lab = 7.5,
            cex.axis =  7,
            cex.main = 8)
  # symbols(species_subset$lon[predictions[[2]]],
  #         species_subset$lat[predictions[[2]]],
  #         circle = predictions[[3]][predictions[[2]]],
  #         inches = 0.12,
  #         add = T,
  #         bg = alpha('darkred', 0.4),
  #         fg = alpha('black', 0.08))
  symbols(species_subset$lon[predictions[[1]]],
          species_subset$lat[predictions[[1]]],
          circle = (-1) * predictions[[3]][predictions[[1]]],
          inches = 0.12,
          add = T,
          bg = alpha('navy', 0.4),
          fg = alpha('black', 0.08))
  maps::map("state",
            boundary = FALSE,
            fill = TRUE,
            col = "wheat4",
            add = TRUE)
  text(x = state_labels$lon, 
       y = state_labels$lat,
       state_labels$name, 
       pos = 1,
       col = "black",
       cex = 6,
       family = "Lato")
  image.plot(legend.only = T,
             col = jet.colors(100),
             legend.shrink = 0.2,
             smallplot = c(.24, .29, .08, .21),
             legend.cex = 4,
             axis.args = list(cex.axis = 5,
                              family = "Lato"),
             legend.width = 0.8,
             legend.mar = 6,
             zlim = c(min(my_gam$linear.predictors),
                      max(my_gam$linear.predictors)),
             legend.args = list("log(cpue+1)",
                                side = 2, 
                                cex = 4,
                                family = "Lato",
                                line = 1.5))
}

plot_var_coef3 <- function(my_gam, species_subset, predictions, yaxis, size){
  myvis_gam(my_gam,
            view = c('lon', 'lat'),
            too.far = 0.07,
            plot.type = 'contour',
            contour.col = contour_col,
            color = "jet" ,
            type = 'link',
            xlim = c(-125.7, -116.5),
            ylim = range(species_subset$lat, na.rm = TRUE) + c(-.4, .5),
            family = "Lato",
            xlab = "Longitude",
            ylab = yaxis,
            main = size,
            cex.lab = 7.5,
            cex.axis =  7,
            cex.main = 8)
  symbols(species_subset$lon[predictions[[2]]],
          species_subset$lat[predictions[[2]]],
          circle = predictions[[3]][predictions[[2]]],
          inches = 0.12,
          add = T,
          bg = alpha('darkred', 0.4),
          fg = alpha('black', 0.08))
  # symbols(species_subset$lon[predictions[[1]]],
  #         species_subset$lat[predictions[[1]]],
  #         circle = (-1) * predictions[[3]][predictions[[1]]],
  #         inches = 0.12,
  #         add = T,
  #         bg = alpha('navy', 0.4),
  #         fg = alpha('black', 0.08))
  maps::map("state",
            boundary = FALSE,
            fill = TRUE,
            col = "wheat4",
            add = TRUE)
  text(x = state_labels$lon, 
       y = state_labels$lat,
       state_labels$name, 
       pos = 1,
       col = "black",
       cex = 6,
       family = "Lato")
  image.plot(legend.only = T,
             col = jet.colors(100),
             legend.shrink = 0.2,
             smallplot = c(.24, .29, .08, .21),
             legend.cex = 4,
             axis.args = list(cex.axis = 5,
                              family = "Lato"),
             legend.width = 0.8,
             legend.mar = 6,
             zlim = c(min(my_gam$linear.predictors),
                      max(my_gam$linear.predictors)),
             legend.args = list("log(cpue+1)",
                                side = 2, 
                                cex = 4,
                                family = "Lato",
                                line = 1.5))
}

RMSE_calc <- function(results, data){
  RMSE <- lapply(results, function(x) {
    rmse(x$catch1, x$pred)
  })
  
  range(data$catch1)
  avg <- mean(unlist(RMSE)) # 185
  
  error <- as.data.frame(do.call(rbind, RMSE))
  error$year <- rownames(error)
  rownames(error) <- NULL
  colnames(error)[1] <- "RMSE"
  return(list(avg, error))
}

RMSE_calc_small <- function(results, data){
  RMSE <- lapply(results, function(x) {
    rmse(x$small_catch1, x$pred_small)
  })
  
  range(data$small_catch1)
  avg <- mean(unlist(RMSE)) # 185
  
  error <- as.data.frame(do.call(rbind, RMSE))
  error$year <- rownames(error)
  rownames(error) <- NULL
  colnames(error)[1] <- "RMSE"
  return(list(avg, error))
}

RMSE_calc_large <- function(results, data){
  RMSE <- lapply(results, function(x) {
    rmse(x$large_catch1, x$pred_large)
  })
  
  range(data$large_catch1)
  avg <- mean(unlist(RMSE)) # 185
  
  error <- as.data.frame(do.call(rbind, RMSE))
  error$year <- rownames(error)
  rownames(error) <- NULL
  colnames(error)[1] <- "RMSE"
  return(list(avg, error))
}

year_adjust <- function(gam, data){
  year_pred <- predict.gam(gam,
                           type = "terms",
                           terms = "year_f")
  year_effect <- gam$family$linkinv(year_pred)
}

LOYO_validation <- function(data, formula){
  hake_gams <- lapply(unique(data$year_f), function(x) {
    output <- gam(formula,
                  family = tw(link = "log"),
                  method = 'REML',
                  data = data[data$year_f != x,])
  }) # Gives the list of GAMs with each year left out
  return(hake_gams)
}

LOYO_preds <- function(gam_list, data, results){
  for(i in seq_along(gam_list)){
    for(j in seq_along(data)){
      results[[j]]$pred <- predict(gam_list[[i]],
                                   newdata = results[[j]],
                                   type = "response")
    }} # Predicts on the left out year
  return(results)
}

LOYO_preds_small <- function(gam_list, data, results){
  for(i in seq_along(gam_list)){
    for(j in seq_along(data)){
      results[[j]]$pred_small <- predict(gam_list[[i]],
                                   newdata = results[[j]],
                                   type = "response")
    }} # Predicts on the left out year
  return(results)
}

LOYO_preds_large <- function(gam_list, data, results){
  for(i in seq_along(gam_list)){
    for(j in seq_along(data)){
      results[[j]]$pred_large <- predict(gam_list[[i]],
                                         newdata = results[[j]],
                                         type = "response")
    }} # Predicts on the left out year
  return(results)
}

# Load data ----
yoy_hake <- read_data('yoy_hake.Rdata')
yoy_anchovy <- read_data('yoy_anch.Rdata')
yoy_widow <- read_data('yoy_widw.Rdata')
yoy_widow <- filter(yoy_widow, catch < 2000) # two large hauls in 2016 caused huge errors
yoy_shortbelly <- read_data('yoy_sbly.Rdata')
yoy_sdab <- read_data('yoy_dab.Rdata')

ctds <- readRDS(here('data', 'ctd_means.rdata'))

ctd_means <- ctds %>%
  group_by(year) %>%
  summarise(across(temperature, mean, na.rm = TRUE))

state_labels <- data.frame(name = c("Washington", "Oregon", "California"),
                           lat = c(46.4, 44.0, 37.0),
                           lon = c(-121.0, -121.0, -120.0))

sysfonts::font_add_google("Lato")
showtext::showtext_auto()

# Hake ----
# Aggregate model
# Use models selected during model exploration
hake_total <- gam(catch1 ~ year_f + 
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

hake_year_effect <- year_adjust(hake_total, yoy_hake)
yoy_hake$y_catch <- yoy_hake$catch1 + hake_year_effect[, 1]

# Leave one group out cross validation
# Run GAMs with each year left out
# Leave out one year, run model on remaining data
hake_formula <- formula(y_catch ~ s(lon, lat) + 
                          s(bottom_depth, k = 4) +
                          s(jday) + 
                          s(temperature, k = 4) +
                          s(salinity, k = 4) +
                          s(lon, lat, by = NPGO_pos)) # Note no factor(year), added into response

hake_gams <- LOYO_validation(yoy_hake, hake_formula)

# Create list with data from each year
hake_data <- split(yoy_hake, yoy_hake$year_f)

# Get predictions
# Predict on the left out year's data
hake_results <- hake_data
hake_results <- LOYO_preds(hake_gams, hake_data, hake_results)

# Calculate RMSE
# Get values for each year and overall value
hake_error <- RMSE_calc(hake_results, yoy_hake)
mean(hake_error[[2]]$RMSE)

# Plot the RMSE for each year
hake_error[[2]]$temperature <- ctd_means$temperature[match(hake_error[[2]]$year, ctd_means$year)]

ggplot(hake_error[[2]]) +
  geom_line(aes(year, RMSE),
             size = 1,
             group = 1) 

ggplot(hake_error[[2]]) +
  geom_line(aes(year, temperature),
            size = 1,
            group = 1) 


# Size explicit
# Small
# Use models selected during model exploration
hake_small <- gam(small_catch1 ~ year_f +
                    s(lon, lat) +
                    s(bottom_depth, k = 4) +
                    s(jday) +
                    s(temperature, k = 4) +
                    s(salinity, k = 4) +
                    s(lon, lat, by = NPGO_pos),
                  family = tw(link = "log"),
                  method = "REML",
                  data = yoy_hake)
summary(hake_small)

hake_year_small <- year_adjust(hake_small, yoy_hake)
yoy_hake$y_small_catch <- yoy_hake$small_catch1 + hake_year_small[, 1]

# Leave one group out cross validation
# Run GAMs with each year left out
# Leave out one year, run model on remaining data
hake_small_formula <- formula(y_small_catch ~ s(lon, lat) +
                                s(bottom_depth, k = 4) +
                                s(jday) +
                                s(temperature, k = 4) +
                                s(salinity, k = 4) +
                                s(lon, lat, by = NPGO_pos)) # Note no year factor, added into response

hake_small_gams <- LOYO_validation(yoy_hake, hake_small_formula)

# Get predictions
# Predict on the left out year's data
hake_small_results <- hake_data
hake_small_results <- LOYO_preds_small(hake_small_gams, hake_data, hake_small_results)

# Calculate RMSE
# Get values for each year and overall value
hake_small_error <- RMSE_calc_small(hake_small_results, yoy_hake)

# Plot the RMSE for each year
hake_small_error[[2]]$temperature <- ctd_means$temperature[match(hake_small_error[[2]]$year, ctd_means$year)]

ggplot(hake_small_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(hake_small_error[[2]]) +
  geom_line(aes(year, temperature),
            size = 1,
            group = 1)

# Large
# Use models selected during model exploration
hake_large <- gam(large_catch1 ~ year_f +
                    s(lon, lat) +
                    s(bottom_depth, k = 4) +
                    s(jday) +
                    s(temperature, k = 4) +
                    s(salinity, k = 4) +
                    s(lon, lat, by = NPGO_pos),
                  family = tw(link = "log"),
                  method = "REML",
                  data = yoy_hake)
summary(hake_large)

hake_year_large <- year_adjust(hake_large, yoy_hake)
yoy_hake$y_large_catch <- yoy_hake$large_catch1 + hake_year_large[, 1]

# Leave one group out cross validation
# Run GAMs with each year left out
# Leave out one year, run model on remaining data
hake_large_formula <- formula(y_large_catch ~ s(lon, lat) +
                                s(bottom_depth, k = 4) +
                                s(jday) +
                                s(temperature, k = 4) +
                                s(salinity, k = 4) +
                                s(lon, lat, by = NPGO_pos)) # Note no year factor, added into response

hake_large_gams <- LOYO_validation(yoy_hake, hake_large_formula)

# Get predictions
# Predict on the left out year's data
hake_large_results <- hake_data
hake_large_results <- LOYO_preds_large(hake_large_gams, hake_data, hake_large_results)

# Calculate RMSE
# Get values for each year and overall value
hake_large_error <- RMSE_calc_large(hake_large_results, yoy_hake)

# Plot the RMSE for each year
hake_large_error[[2]]$temperature <- ctd_means$temperature[match(hake_large_error[[2]]$year, ctd_means$year)]

ggplot(hake_large_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(hake_large_error[[2]]) +
  geom_line(aes(year, temperature),
            size = 1,
            group = 1)

# Add the large and small together
hake_combined_results <- map2(hake_large_results, 
                              hake_small_results, 
                              ~left_join(.x, .y))

hake_added_results <- lapply(hake_combined_results, function(x){
  rmse(x$catch1, x$pred_small + x$pred_large)
})

mean(unlist(hake_added_results)) # 182

hake_combined_df <- data.frame(year = names(hake_added_results), 
                               RMSE = unlist(hake_added_results))


ggplot(hake_combined_df) +
  geom_line(aes(year, RMSE),
            size = 1.3,
            group = 1,
            color = "maroon4") +
  labs(x = "Year",
       y = "RMSE",
       title = "Yearly Error for Pacific Hake") +
  scale_x_discrete(breaks = c(1987, 1997, 2007, 2017)) +
  theme_tufte() +
  theme(axis.title = element_text(size = 26,
                                  color = "maroon4",
                                  family = "Lato"),
        axis.text = element_text(size = 22,
                                 family = "Lato"),
        plot.title = element_text(size = 28, 
                                  face = "bold",
                                  family = "Lato"),
        axis.line = element_line(color = "black"))

# Maps
# General distributions
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

# Variable coefficient plots
pred_hake_all <- variable_coefficient(hake_total, yoy_hake, yoy_hake$NPGO_pos)
pred_hake_small <- variable_coefficient(hake_small, yoy_hake, yoy_hake$NPGO_pos)
pred_hake_large <- variable_coefficient(hake_large, yoy_hake, yoy_hake$NPGO_pos)

windows()
par(mfrow = c(1, 3),
    mar = c(6.4, 7.2, 2.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
plot_var_coef(hake_total, yoy_hake, pred_hake_all, "Latitude", "All Sizes")
plot_var_coef(hake_small, yoy_hake, pred_hake_small, "", "Small Sizes")
plot_var_coef2(hake_large, yoy_hake, pred_hake_large, "", "Large Sizes")
dev.copy(jpeg, here('results/RREAS_preliminary', 'yoy_hake_var_coef.jpg'), 
         height = 15, width = 20, units = 'in', res = 200)
dev.off()

# Anchovy ----
# Aggregate model
# Use models selected during model exploration
anchovy_total <- gam(catch1 ~ year_f + 
                    s(lon, lat) + 
                    s(bottom_depth, k = 4) +
                    s(jday) + 
                    s(temperature, k = 4) +
                    s(salinity, k = 4) +
                    s(lon, lat, by = NPGO_pos),
                  family = tw(link = "log"),
                  method = "REML",
                  data = yoy_anchovy)
summary(anchovy_total)

anchovy_year_effect <- year_adjust(anchovy_total, yoy_anchovy)
yoy_anchovy$y_catch <- yoy_anchovy$catch1 + anchovy_year_effect[, 1]

# Leave one group out cross validation
# Run GAMs with each year left out
# Leave out one year, run model on remaining data
anchovy_formula <- formula(y_catch ~ s(lon, lat) + 
                          s(bottom_depth, k = 4) +
                          s(jday) + 
                          s(temperature, k = 4) +
                          s(salinity, k = 4) +
                          s(lon, lat, by = NPGO_pos)) # Note no factor(year), added into response

anchovy_gams <- LOYO_validation(yoy_anchovy, anchovy_formula)

# Create list with data from each year
anchovy_data <- split(yoy_anchovy, yoy_anchovy$year_f)

# Get predictions
# Predict on the left out year's data
anchovy_results <- anchovy_data
anchovy_results <- LOYO_preds(anchovy_gams, anchovy_data, anchovy_results)

# Calculate RMSE
# Get values for each year and overall value
anchovy_error <- RMSE_calc(anchovy_results, yoy_anchovy)
mean(anchovy_error[[2]]$RMSE) #371

# Plot the RMSE for each year
anchovy_error[[2]]$temperature <- ctd_means$temperature[match(anchovy_error[[2]]$year, ctd_means$year)]

ggplot(anchovy_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(anchovy_error[[2]]) +
  geom_line(aes(year, temperature),
            size = 1,
            group = 1) 


# Size explicit
# Small
# Use models selected during model exploration
anchovy_small <- gam(small_catch1 ~ year_f +
                    s(lon, lat) +
                    s(bottom_depth, k = 4) +
                    s(jday) +
                    s(temperature, k = 4) +
                    s(salinity, k = 4) +
                    s(lon, lat, by = NPGO_pos),
                  family = tw(link = "log"),
                  method = "REML",
                  data = yoy_anchovy)
summary(anchovy_small)

anchovy_year_small <- year_adjust(anchovy_small, yoy_anchovy)
yoy_anchovy$y_small_catch <- yoy_anchovy$small_catch1 + anchovy_year_small[, 1]

# Leave one group out cross validation
# Run GAMs with each year left out
# Leave out one year, run model on remaining data
anchovy_small_formula <- formula(y_small_catch ~ s(lon, lat) +
                                s(bottom_depth, k = 4) +
                                s(jday) +
                                s(temperature, k = 4) +
                                s(salinity, k = 4) +
                                s(lon, lat, by = NPGO_pos)) # Note no year factor, added into response

anchovy_small_gams <- LOYO_validation(yoy_anchovy, anchovy_small_formula)

# Get predictions
# Predict on the left out year's data
anchovy_small_results <- anchovy_data
anchovy_small_results <- LOYO_preds_small(anchovy_small_gams, anchovy_data, anchovy_small_results)

# Calculate RMSE
# Get values for each year and overall value
anchovy_small_error <- RMSE_calc_small(anchovy_small_results, yoy_anchovy)

# Plot the RMSE for each year
anchovy_small_error[[2]]$temperature <- ctd_means$temperature[match(anchovy_small_error[[2]]$year, ctd_means$year)]

ggplot(anchovy_small_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(anchovy_small_error[[2]]) +
  geom_line(aes(year, temperature),
            size = 1,
            group = 1)

# Large
# Use models selected during model exploration
anchovy_large <- gam(large_catch1 ~ year_f +
                    s(lon, lat) +
                    s(bottom_depth, k = 4) +
                    s(jday) +
                    s(temperature, k = 4) +
                    s(salinity, k = 4) +
                    s(lon, lat, by = NPGO_pos),
                  family = tw(link = "log"),
                  method = "REML",
                  data = yoy_anchovy)
summary(anchovy_large)

anchovy_year_large <- year_adjust(anchovy_large, yoy_anchovy)
yoy_anchovy$y_large_catch <- yoy_anchovy$large_catch1 + anchovy_year_large[, 1]

# Leave one group out cross validation
# Run GAMs with each year left out
# Leave out one year, run model on remaining data
anchovy_large_formula <- formula(y_large_catch ~ s(lon, lat) +
                                s(bottom_depth, k = 4) +
                                s(jday) +
                                s(temperature, k = 4) +
                                s(salinity, k = 4) +
                                s(lon, lat, by = NPGO_pos)) # Note no year factor, added into response

anchovy_large_gams <- LOYO_validation(yoy_anchovy, anchovy_large_formula)

# Get predictions
# Predict on the left out year's data
anchovy_large_results <- anchovy_data
anchovy_large_results <- LOYO_preds_large(anchovy_large_gams, anchovy_data, anchovy_large_results)

# Calculate RMSE
# Get values for each year and overall value
anchovy_large_error <- RMSE_calc_large(anchovy_large_results, yoy_anchovy)

anchovy_combined_df <- data.frame(year = names(anchovy_added_results), 
                               RMSE = unlist(anchovy_added_results))


ggplot(anchovy_combined_df) +
  geom_line(aes(year, RMSE),
            size = 1.3,
            group = 1,
            color = "maroon4") +
  labs(x = "Year",
       y = "RMSE",
       title = "Yearly RMSE for Northern Anchovy") +
  theme_tufte() +
  theme(axis.title = element_text(size = 24,
                                  color = "maroon4"),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 36),
        axis.line = element_line(color = "black"))

# Plot the RMSE for each year
anchovy_large_error[[2]]$temperature <- ctd_means$temperature[match(anchovy_large_error[[2]]$year, ctd_means$year)]

ggplot(anchovy_large_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(anchovy_large_error[[2]]) +
  geom_line(aes(year, temperature),
            size = 1,
            group = 1)

# Add the large and small together
anchovy_combined_results <- map2(anchovy_large_results, 
                              anchovy_small_results, 
                              ~left_join(.x, .y))

anchovy_added_results <- lapply(anchovy_combined_results, function(x){
  rmse(x$catch1, x$pred_small + x$pred_large)
})

mean(unlist(anchovy_added_results)) # 321

anchovy_combined_df <- data.frame(year = names(anchovy_added_results), 
                               RMSE = unlist(anchovy_added_results))

ggplot(anchovy_combined_df) +
  geom_line(aes(year, RMSE),
            size = 1.3,
            group = 1,
            color = "maroon4") +
  labs(x = "Year",
       y = "RMSE",
       title = "Yearly RMSE for Northern Anchovy") +
  theme_tufte() +
  theme(axis.title = element_text(size = 24,
                                  color = "maroon4"),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 36),
        axis.line = element_line(color = "black"))

# Maps
par(mfrow = c(1, 3),
    mar = c(6.4, 7.2, 2.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
location_plot(anchovy_total, yoy_anchovy, "lat", "All Sizes", yoy_anchovy$catch)
location_plot(anchovy_large, yoy_anchovy, " ", "Small Sizes", yoy_anchovy$large)
location_plot(anchovy_small, yoy_anchovy, " ", "Large Sizes", yoy_anchovy$small)
dev.copy(jpeg, here('results/RREAS_preliminary', 'anchovy_distributions.jpg'), 
         height = 15, width = 20, units = 'in', res = 200)
dev.off()

# Variable coefficient plots
pred_anchovy_all <- variable_coefficient(anchovy_total, yoy_anchovy, yoy_anchovy$NPGO_pos)
pred_anchovy_small <- variable_coefficient(anchovy_small, yoy_anchovy, yoy_anchovy$NPGO_pos)
pred_anchovy_large <- variable_coefficient(anchovy_large, yoy_anchovy, yoy_anchovy$NPGO_pos)

windows()
par(mfrow = c(1, 3),
    mar = c(6.4, 7.2, 2.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
plot_var_coef3(anchovy_total, yoy_anchovy, pred_anchovy_all, "Latitude", "All Sizes")
plot_var_coef3(anchovy_small, yoy_anchovy, pred_anchovy_small, "", "Small Sizes")
plot_var_coef3(anchovy_large, yoy_anchovy, pred_anchovy_large, "", "Large Sizes")
dev.copy(jpeg, here('results/RREAS_preliminary', 'yoy_anchovy_var_coef.jpg'), 
         height = 15, width = 20, units = 'in', res = 200)
dev.off()

# Widow Rockfish ----
# Aggregate model
# Use models selected during model exploration
widow_total <- gam(catch1 ~ year_f + 
                    s(lon, lat) + 
                    s(bottom_depth, k = 4) +
                    s(jday) + 
                    s(temperature, k = 4) +
                    s(salinity, k = 4) +
                    s(lon, lat, by = ONI_pos),
                  family = tw(link = "log"),
                  method = "REML",
                  data = yoy_widow)
summary(widow_total)

widow_year_effect <- year_adjust(widow_total, yoy_widow)
yoy_widow$y_catch <- yoy_widow$catch1 + widow_year_effect[, 1]

# Leave one group out cross validation
# Run GAMs with each year left out
# Leave out one year, run model on remaining data
widow_formula <- formula(y_catch ~ s(lon, lat) + 
                          s(bottom_depth, k = 4) +
                          s(jday) + 
                          s(temperature, k = 4) +
                          s(salinity, k = 4) +
                          s(lon, lat, by = ONI_pos)) # Note no factor(year), added into response

widow_gams <- LOYO_validation(yoy_widow, widow_formula)

# Create list with data from each year
widow_data <- split(yoy_widow, yoy_widow$year_f)

# Get predictions
# Predict on the left out year's data
widow_results <- widow_data
widow_results <- LOYO_preds(widow_gams, widow_data, widow_results)

# Calculate RMSE
# Get values for each year and overall value
widow_error <- RMSE_calc(widow_results, yoy_widow)
mean(widow_error[[2]]$RMSE) #485

# Plot the RMSE for each year
widow_error[[2]]$temperature <- ctd_means$temperature[match(widow_error[[2]]$year, ctd_means$year)]

ggplot(widow_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(widow_error[[2]]) +
  geom_line(aes(year, temperature),
            size = 1,
            group = 1) 


# Size explicit
# Small
# Use models selected during model exploration
widow_small <- gam(small_catch1 ~ year_f +
                    s(lon, lat) +
                    s(bottom_depth, k = 4) +
                    s(jday) +
                    s(temperature, k = 4) +
                    s(salinity, k = 4) +
                    s(lon, lat, by = ONI_pos),
                  family = tw(link = "log"),
                  method = "REML",
                  data = yoy_widow)
summary(widow_small)

widow_year_small <- year_adjust(widow_small, yoy_widow)
yoy_widow$y_small_catch <- yoy_widow$small_catch1 + widow_year_small[, 1]

# Leave one group out cross validation
# Run GAMs with each year left out
# Leave out one year, run model on remaining data
widow_small_formula <- formula(y_small_catch ~ s(lon, lat) +
                                s(bottom_depth, k = 4) +
                                s(jday) +
                                s(temperature, k = 4) +
                                s(salinity, k = 4) +
                                s(lon, lat, by = ONI_pos)) # Note no year factor, added into response

widow_small_gams <- LOYO_validation(yoy_widow, widow_small_formula)

# Get predictions
# Predict on the left out year's data
widow_small_results <- widow_data
widow_small_results <- LOYO_preds_small(widow_small_gams, widow_data, widow_small_results)

# Calculate RMSE
# Get values for each year and overall value
widow_small_error <- RMSE_calc_small(widow_small_results, yoy_widow)

# Plot the RMSE for each year
widow_small_error[[2]]$temperature <- ctd_means$temperature[match(widow_small_error[[2]]$year, ctd_means$year)]

ggplot(widow_small_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(widow_small_error[[2]]) +
  geom_line(aes(year, temperature),
            size = 1,
            group = 1)

# Large
# Use models selected during model exploration
widow_large <- gam(large_catch1 ~ year_f +
                    s(lon, lat) +
                    s(bottom_depth, k = 4) +
                    s(jday) +
                    s(temperature, k = 4) +
                    s(salinity, k = 4) +
                    s(lon, lat, by = ONI_pos),
                  family = tw(link = "log"),
                  method = "REML",
                  data = yoy_widow)
summary(widow_large)

widow_year_large <- year_adjust(widow_large, yoy_widow)
yoy_widow$y_large_catch <- yoy_widow$large_catch1 + widow_year_large[, 1]

# Leave one group out cross validation
# Run GAMs with each year left out
# Leave out one year, run model on remaining data
widow_large_formula <- formula(y_large_catch ~ s(lon, lat) +
                                s(bottom_depth, k = 4) +
                                s(jday) +
                                s(temperature, k = 4) +
                                s(salinity, k = 4) +
                                s(lon, lat, by = ONI_pos)) # Note no year factor, added into response

widow_large_gams <- LOYO_validation(yoy_widow, widow_large_formula)

# Get predictions
# Predict on the left out year's data
widow_large_results <- widow_data
widow_large_results <- LOYO_preds_large(widow_large_gams, widow_data, widow_large_results)

# Calculate RMSE
# Get values for each year and overall value
widow_large_error <- RMSE_calc_large(widow_large_results, yoy_widow)

# Plot the RMSE for each year
widow_large_error[[2]]$temperature <- ctd_means$temperature[match(widow_large_error[[2]]$year, ctd_means$year)]

ggplot(widow_large_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(widow_large_error[[2]]) +
  geom_line(aes(year, temperature),
            size = 1,
            group = 1)

# Add the large and small together
widow_combined_results <- map2(widow_large_results, 
                              widow_small_results, 
                              ~left_join(.x, .y))

widow_added_results <- lapply(widow_combined_results, function(x){
  rmse(x$catch1, x$pred_small + x$pred_large)
})

mean(unlist(widow_added_results)) # 25
# When including the outliers, splitting up the data into size bins drastically improves the RMSE
# When aggregated with the outliers, the RMSE is extremely high

widow_combined_df <- data.frame(year = names(widow_added_results), 
                               RMSE = unlist(widow_added_results))

ggplot(widow_combined_df) +
  geom_line(aes(year, RMSE),
            size = 1.3,
            group = 1,
            color = "maroon4") +
  labs(x = "Year",
       y = "RMSE",
       title = "Yearly RMSE for Widow Rockfish") +
  scale_x_discrete(breaks = c(1987, 1997, 2007, 2017)) +
  theme_tufte() +
  theme(axis.title = element_text(size = 24,
                                  color = "maroon4"),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 36),
        axis.line = element_line(color = "black"))

# Maps
par(mfrow = c(1, 3),
    mar = c(6.4, 7.2, 2.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
location_plot(widow_total, yoy_widow, "lat", "All Sizes", yoy_widow$catch)
location_plot(widow_large, yoy_widow, " ", "Small Sizes", yoy_widow$large)
location_plot(widow_small, yoy_widow, " ", "Large Sizes", yoy_widow$small)
dev.copy(jpeg, here('results/RREAS_preliminary', 'widow_distributions.jpg'), 
         height = 15, width = 20, units = 'in', res = 200)
dev.off()

# Variable coefficient plots
pred_widow_all <- variable_coefficient(widow_total, yoy_widow, yoy_widow$ONI_pos)
pred_widow_small <- variable_coefficient(widow_small, yoy_widow, yoy_widow$ONI_pos)
pred_widow_large <- variable_coefficient(widow_large, yoy_widow, yoy_widow$ONI_pos)

windows()
par(mfrow = c(1, 3),
    mar = c(6.4, 7.2, 2.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
plot_var_coef(widow_total, yoy_widow, pred_widow_all, "Latitude", "All Sizes")
plot_var_coef(widow_small, yoy_widow, pred_widow_small, "", "Small Sizes")
plot_var_coef3(widow_large, yoy_widow, pred_widow_large, "", "Large Sizes")
dev.copy(jpeg, here('results/RREAS_preliminary', 'yoy_widow_var_coef.jpg'), 
         height = 15, width = 20, units = 'in', res = 200)
dev.off()


# Shortbelly Rockfish ----
# Aggregate model
# Use models selected during model exploration
shortbelly_total <- gam(catch1 ~ year_f + 
                    s(lon, lat) + 
                    s(bottom_depth, k = 4) +
                    s(jday) + 
                    s(temperature, k = 4) +
                    s(salinity, k = 4) +
                    s(lon, lat, by = PDO_pos),
                  family = tw(link = "log"),
                  method = "REML",
                  data = yoy_shortbelly)
summary(shortbelly_total)

shortbelly_year_effect <- year_adjust(shortbelly_total, yoy_shortbelly)
yoy_shortbelly$y_catch <- yoy_shortbelly$catch1 + shortbelly_year_effect[, 1]

# Leave one group out cross validation
# Run GAMs with each year left out
# Leave out one year, run model on remaining data
shortbelly_formula <- formula(y_catch ~ s(lon, lat) + 
                          s(bottom_depth, k = 4) +
                          s(jday) + 
                          s(temperature, k = 4) +
                          s(salinity, k = 4) +
                          s(lon, lat, by = PDO_pos)) # Note no factor(year), added into response

shortbelly_gams <- LOYO_validation(yoy_shortbelly, shortbelly_formula)

# Create list with data from each year
shortbelly_data <- split(yoy_shortbelly, yoy_shortbelly$year_f)

# Get predictions
# Predict on the left out year's data
shortbelly_results <- shortbelly_data
shortbelly_results <- LOYO_preds(shortbelly_gams, shortbelly_data, shortbelly_results)

# Calculate RMSE
# Get values for each year and overall value
shortbelly_error <- RMSE_calc(shortbelly_results, yoy_shortbelly)
mean(shortbelly_error[[2]]$RMSE) #111

# Plot the RMSE for each year
shortbelly_error[[2]]$temperature <- ctd_means$temperature[match(shortbelly_error[[2]]$year, ctd_means$year)]

ggplot(shortbelly_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(shortbelly_error[[2]]) +
  geom_line(aes(year, temperature),
            size = 1,
            group = 1) 


# Size explicit
# Small
# Use models selected during model exploration
shortbelly_small <- gam(small_catch1 ~ year_f +
                    s(lon, lat) +
                    s(bottom_depth, k = 4) +
                    s(jday) +
                    s(temperature, k = 4) +
                    s(salinity, k = 4) +
                    s(lon, lat, by = PDO_pos),
                  family = tw(link = "log"),
                  method = "REML",
                  data = yoy_shortbelly)
summary(shortbelly_small)

shortbelly_year_small <- year_adjust(shortbelly_small, yoy_shortbelly)
yoy_shortbelly$y_small_catch <- yoy_shortbelly$small_catch1 + shortbelly_year_small[, 1]

# Leave one group out cross validation
# Run GAMs with each year left out
# Leave out one year, run model on remaining data
shortbelly_small_formula <- formula(y_small_catch ~ s(lon, lat) +
                                s(bottom_depth, k = 4) +
                                s(jday) +
                                s(temperature, k = 4) +
                                s(salinity, k = 4) +
                                s(lon, lat, by = PDO_pos)) # Note no year factor, added into response

shortbelly_small_gams <- LOYO_validation(yoy_shortbelly, shortbelly_small_formula)

# Get predictions
# Predict on the left out year's data
shortbelly_small_results <- shortbelly_data
shortbelly_small_results <- LOYO_preds_small(shortbelly_small_gams, shortbelly_data, shortbelly_small_results)

# Calculate RMSE
# Get values for each year and overall value
shortbelly_small_error <- RMSE_calc_small(shortbelly_small_results, yoy_shortbelly)

# Plot the RMSE for each year
shortbelly_small_error[[2]]$temperature <- ctd_means$temperature[match(shortbelly_small_error[[2]]$year, ctd_means$year)]

ggplot(shortbelly_small_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(shortbelly_small_error[[2]]) +
  geom_line(aes(year, temperature),
            size = 1,
            group = 1)

# Large
# Use models selected during model exploration
shortbelly_large <- gam(large_catch1 ~ year_f +
                    s(lon, lat) +
                    s(bottom_depth, k = 4) +
                    s(jday) +
                    s(temperature, k = 4) +
                    s(salinity, k = 4) +
                    s(lon, lat, by = PDO_pos),
                  family = tw(link = "log"),
                  method = "REML",
                  data = yoy_shortbelly)
summary(shortbelly_large)

shortbelly_year_large <- year_adjust(shortbelly_large, yoy_shortbelly)
yoy_shortbelly$y_large_catch <- yoy_shortbelly$large_catch1 + shortbelly_year_large[, 1]

# Leave one group out cross validation
# Run GAMs with each year left out
# Leave out one year, run model on remaining data
shortbelly_large_formula <- formula(y_large_catch ~ s(lon, lat) +
                                s(bottom_depth, k = 4) +
                                s(jday) +
                                s(temperature, k = 4) +
                                s(salinity, k = 4) +
                                s(lon, lat, by = PDO_pos)) # Note no year factor, added into response

shortbelly_large_gams <- LOYO_validation(yoy_shortbelly, shortbelly_large_formula)

# Get predictions
# Predict on the left out year's data
shortbelly_large_results <- shortbelly_data
shortbelly_large_results <- LOYO_preds_large(shortbelly_large_gams, shortbelly_data, shortbelly_large_results)

# Calculate RMSE
# Get values for each year and overall value
shortbelly_large_error <- RMSE_calc_large(shortbelly_large_results, yoy_shortbelly)

# Plot the RMSE for each year
shortbelly_large_error[[2]]$temperature <- ctd_means$temperature[match(shortbelly_large_error[[2]]$year, ctd_means$year)]

ggplot(shortbelly_large_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(shortbelly_large_error[[2]]) +
  geom_line(aes(year, temperature),
            size = 1,
            group = 1)

# Add the large and small together
shortbelly_combined_results <- map2(shortbelly_large_results, 
                              shortbelly_small_results, 
                              ~left_join(.x, .y))

shortbelly_added_results <- lapply(shortbelly_combined_results, function(x){
  rmse(x$catch1, x$pred_small + x$pred_large)
})

mean(unlist(shortbelly_added_results)) # 110

shortbelly_combined_df <- data.frame(year = names(shortbelly_added_results), 
                               RMSE = unlist(shortbelly_added_results))


ggplot(shortbelly_combined_df) +
  geom_line(aes(year, RMSE),
            size = 1.3,
            group = 1,
            color = "maroon4") +
  labs(x = "Year",
       y = "RMSE",
       title = "Yearly Error for Shortbelly Rockfish") +
  scale_x_discrete(breaks = c(1987, 1997, 2007, 2017)) +
  theme_tufte() +
  theme(axis.title = element_text(size = 26,
                                  color = "maroon4",
                                  family = "Lato"),
        axis.text = element_text(size = 22,
                                 family = "Lato"),
        plot.title = element_text(size = 28, 
                                  face = "bold",
                                  family = "Lato"),
        axis.line = element_line(color = "black"))

# Maps
par(mfrow = c(1, 3),
    mar = c(6.4, 7.2, 2.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
location_plot(shortbelly_total, yoy_shortbelly, "lat", "All Sizes", yoy_shortbelly$catch)
location_plot(shortbelly_large, yoy_shortbelly, " ", "Small Sizes", yoy_shortbelly$large)
location_plot(shortbelly_small, yoy_shortbelly, " ", "Large Sizes", yoy_shortbelly$small)
dev.copy(jpeg, here('results/RREAS_preliminary', 'shortbelly_distributions.jpg'), 
         height = 15, width = 20, units = 'in', res = 200)
dev.off()

# Variable coefficient plots
pred_shortbelly_all <- variable_coefficient(shortbelly_total, yoy_shortbelly, yoy_shortbelly$PDO_pos)
pred_shortbelly_small <- variable_coefficient(shortbelly_small, yoy_shortbelly, yoy_shortbelly$PDO_pos)
pred_shortbelly_large <- variable_coefficient(shortbelly_large, yoy_shortbelly, yoy_shortbelly$PDO_pos)

windows()
par(mfrow = c(1, 3),
    mar = c(6.4, 7.2, 2.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
plot_var_coef(shortbelly_total, yoy_shortbelly, pred_shortbelly_all, "Latitude", "All Sizes")
plot_var_coef2(shortbelly_small, yoy_shortbelly, pred_shortbelly_small, "", "Small Sizes")
plot_var_coef3(shortbelly_large, yoy_shortbelly, pred_shortbelly_large, "", "Large Sizes")
dev.copy(jpeg, here('results/RREAS_preliminary', 'yoy_shortbelly_var_coef.jpg'), 
         height = 15, width = 20, units = 'in', res = 200)
dev.off()

# Pacific Sanddab ----
# Aggregate model
# Use models selected during model exploration
sdab_total <- gam(catch1 ~ year_f + 
                    s(lon, lat) + 
                    s(bottom_depth, k = 4) +
                    s(jday) + 
                    s(temperature, k = 4) +
                    s(salinity, k = 4) +
                    s(lon, lat, by = NPGO_pos),
                  family = tw(link = "log"),
                  method = "REML",
                  data = yoy_sdab)
summary(sdab_total)

sdab_year_effect <- year_adjust(sdab_total, yoy_sdab)
yoy_sdab$y_catch <- yoy_sdab$catch1 + sdab_year_effect[, 1]

# Leave one group out cross validation
# Run GAMs with each year left out
# Leave out one year, run model on remaining data
sdab_formula <- formula(y_catch ~ s(lon, lat) + 
                          s(bottom_depth, k = 4) +
                          s(jday) + 
                          s(temperature, k = 4) +
                          s(salinity, k = 4) +
                          s(lon, lat, by = NPGO_pos)) # Note no factor(year), added into response

sdab_gams <- LOYO_validation(yoy_sdab, sdab_formula)

# Create list with data from each year
sdab_data <- split(yoy_sdab, yoy_sdab$year_f)

# Get predictions
# Predict on the left out year's data
sdab_results <- sdab_data
sdab_results <- LOYO_preds(sdab_gams, sdab_data, sdab_results)

# Calculate RMSE
# Get values for each year and overall value
sdab_error <- RMSE_calc(sdab_results, yoy_sdab)
mean(sdab_error[[2]]$RMSE) #65

# Plot the RMSE for each year
sdab_error[[2]]$temperature <- ctd_means$temperature[match(sdab_error[[2]]$year, ctd_means$year)]

ggplot(sdab_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1.3,
            group = 1,
            color = "darkslategray4") +
  labs(x = "Year",
       y = "RMSE",
       title = "Aggregate Pacific Sanddab") +
  scale_x_discrete(breaks = c(1987, 1997, 2007, 2017)) +
  theme_tufte() +
  theme(axis.title = element_text(size = 24,
                                  color = "darkslategray4"),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 36),
        axis.line = element_line(color = "black"))

ggplot(sdab_error[[2]]) +
  geom_line(aes(year, temperature),
            size = 1,
            group = 1) 


# Size explicit
# Small
# Use models selected during model exploration
sdab_small <- gam(small_catch1 ~ year_f +
                    s(lon, lat) +
                    s(bottom_depth, k = 4) +
                    s(jday) +
                    s(temperature, k = 4) +
                    s(salinity, k = 4) +
                    s(lon, lat, by = NPGO_pos),
                  family = tw(link = "log"),
                  method = "REML",
                  data = yoy_sdab)
summary(sdab_small)

sdab_year_small <- year_adjust(sdab_small, yoy_sdab)
yoy_sdab$y_small_catch <- yoy_sdab$small_catch1 + sdab_year_small[, 1]

# Leave one group out cross validation
# Run GAMs with each year left out
# Leave out one year, run model on remaining data
sdab_small_formula <- formula(y_small_catch ~ s(lon, lat) +
                                s(bottom_depth, k = 4) +
                                s(jday) +
                                s(temperature, k = 4) +
                                s(salinity, k = 4) +
                                s(lon, lat, by = NPGO_pos)) # Note no year factor, added into response

sdab_small_gams <- LOYO_validation(yoy_sdab, sdab_small_formula)

# Get predictions
# Predict on the left out year's data
sdab_small_results <- sdab_data
sdab_small_results <- LOYO_preds_small(sdab_small_gams, sdab_data, sdab_small_results)

# Calculate RMSE
# Get values for each year and overall value
sdab_small_error <- RMSE_calc_small(sdab_small_results, yoy_sdab)

# Plot the RMSE for each year
sdab_small_error[[2]]$temperature <- ctd_means$temperature[match(sdab_small_error[[2]]$year, ctd_means$year)]

ggplot(sdab_small_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(sdab_small_error[[2]]) +
  geom_line(aes(year, temperature),
            size = 1,
            group = 1)

# Large
# Use models selected during model exploration
sdab_large <- gam(large_catch1 ~ year_f +
                    s(lon, lat) +
                    s(bottom_depth, k = 4) +
                    s(jday) +
                    s(temperature, k = 4) +
                    s(salinity, k = 4) +
                    s(lon, lat, by = NPGO_pos),
                  family = tw(link = "log"),
                  method = "REML",
                  data = yoy_sdab)
summary(sdab_large)

sdab_year_large <- year_adjust(sdab_large, yoy_sdab)
yoy_sdab$y_large_catch <- yoy_sdab$large_catch1 + sdab_year_large[, 1]

# Leave one group out cross validation
# Run GAMs with each year left out
# Leave out one year, run model on remaining data
sdab_large_formula <- formula(y_large_catch ~ s(lon, lat) +
                                s(bottom_depth, k = 4) +
                                s(jday) +
                                s(temperature, k = 4) +
                                s(salinity, k = 4) +
                                s(lon, lat, by = NPGO_pos)) # Note no year factor, added into response

sdab_large_gams <- LOYO_validation(yoy_sdab, sdab_large_formula)

# Get predictions
# Predict on the left out year's data
sdab_large_results <- sdab_data
sdab_large_results <- LOYO_preds_large(sdab_large_gams, sdab_data, sdab_large_results)

# Calculate RMSE
# Get values for each year and overall value
sdab_large_error <- RMSE_calc_large(sdab_large_results, yoy_sdab)

# Plot the RMSE for each year
sdab_large_error[[2]]$temperature <- ctd_means$temperature[match(sdab_large_error[[2]]$year, ctd_means$year)]

ggplot(sdab_large_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(sdab_large_error[[2]]) +
  geom_line(aes(year, temperature),
            size = 1,
            group = 1)

# Add the large and small together
sdab_combined_results <- map2(sdab_large_results, 
                              sdab_small_results, 
                              ~left_join(.x, .y))

sdab_added_results <- lapply(sdab_combined_results, function(x){
  rmse(x$catch1, x$pred_small + x$pred_large)
})

mean(unlist(sdab_added_results)) # 54

sdab_combined_df <- data.frame(year = names(sdab_added_results), 
                               RMSE = unlist(sdab_added_results))

ggplot(sdab_combined_df) +
  geom_line(aes(year, RMSE),
            size = 1.3,
            group = 1,
            color = "maroon4") +
  labs(x = "Year",
       y = "RMSE",
       title = "Yearly RMSE for Pacific Sanddab") +
  scale_x_discrete(breaks = c(1987, 1997, 2007, 2017)) +
  theme_tufte() +
  theme(axis.title = element_text(size = 24,
                                  color = "maroon4"),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 36),
        axis.line = element_line(color = "black"))

# Maps
par(mfrow = c(1, 3),
    mar = c(6.4, 7.2, 2.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
location_plot(sdab_total, yoy_sdab, "lat", "All Sizes", yoy_sdab$catch)
location_plot(sdab_large, yoy_sdab, " ", "Small Sizes", yoy_sdab$large)
location_plot(sdab_small, yoy_sdab, " ", "Large Sizes", yoy_sdab$small)
dev.copy(jpeg, here('results/RREAS_preliminary', 'sdab_distributions.jpg'), 
         height = 15, width = 20, units = 'in', res = 200)
dev.off()

# Variable coefficient plots
pred_sdab_all <- variable_coefficient(sdab_total, yoy_sdab, yoy_sdab$NPGO_pos)
pred_sdab_small <- variable_coefficient(sdab_small, yoy_sdab, yoy_sdab$NPGO_pos)
pred_sdab_large <- variable_coefficient(sdab_large, yoy_sdab, yoy_sdab$NPGO_pos)

windows()
par(mfrow = c(1, 3),
    mar = c(6.4, 7.2, 2.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
plot_var_coef2(sdab_total, yoy_sdab, pred_sdab_all, "Latitude", "All Sizes")
plot_var_coef(sdab_small, yoy_sdab, pred_sdab_small, "", "Small Sizes")
plot_var_coef3(sdab_large, yoy_sdab, pred_sdab_large, "", "Large Sizes")
dev.copy(jpeg, here('results/RREAS_preliminary', 'yoy_sdab_var_coef.jpg'), 
         height = 15, width = 20, units = 'in', res = 200)
dev.off()