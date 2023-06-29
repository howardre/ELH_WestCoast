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
    tidyr::drop_na(roms_temperature, roms_salinity, roms_ssh, bottom_depth, year, jday, lat, lon) %>%
    filter(catch < 2500 &
             year != 2020 &
             year < 2019 &
             year > 2010) %>%
    mutate(catch1 = catch + 1,
           small_catch1 = small + 1,
           large_catch1 = large + 1,
           year_f = as.factor(year),
           ssh_pos = year_ssh + abs(min(year_ssh)) + 10)
  yoy <- yoy[!(yoy$small == 0 & yoy$large == 0 & yoy$catch > 0), ]
  return(yoy)
}

contour_col <- rgb(0, 0, 255, max = 255, alpha = 0, names = "white")
jet.colors <- colorRampPalette(c(sequential_hcl(15, palette = "BluGrn")))

location_plot <- function(gam, species_subset, yaxis, title, value) {
  myvis_gam(gam,
            view = c('lon', 'lat'),
            too.far = 0.07,
            plot.type = 'contour',
            contour.col = contour_col,
            color = "jet" ,
            type = 'response',
            xlim = c(-125.7, -116.5),
            ylim = range(species_subset$lat, na.rm = TRUE) + c(-.4, .5),
            family = "serif",
            xlab = "Longitude",
            ylab = yaxis,
            main = title,
            cex.lab = 4.5,
            cex.axis =  4.2,
            cex.main = 4.8)
  symbols(species_subset$lon,
          species_subset$lat,
          circle = value,
          inches = 0.3,
          add = T,
          bg = alpha('violetred3', 0.4),
          fg = alpha('violetred4', 0.1))  
  maps::map('worldHires',
            add = T,
            col = 'antiquewhite4',
            fill = T)
  image.plot(legend.only = T,
             col = jet.colors(100),
             legend.shrink = 0.2,
             smallplot = c(.24, .29, .08, .21),
             legend.cex = 2.5,
             axis.args = list(cex.axis = 3.5,
                              family = "serif"),
             legend.width = 0.8,
             legend.mar = 6,
             zlim = c(min(gam$linear.predictors),
                      max(gam$linear.predictors)),
             legend.args = list("log(CPUE+1)",
                                side = 2, 
                                cex = 2.5,
                                family = "serif",
                                line = 1.5))
}

variable_coefficient <- function(gam, data, variable, number){
  preds <- predict(gam, type = 'terms', se.fit = T)
  pred_slope <- preds[[1]][, number] / variable # change number depending on if terms are removed
  pred_slope_se <- 1.96 * preds[[2]][, number]
  pred_slope_up <- (preds[[1]][, number] + pred_slope_se) / variable
  pred_slope_low <- (preds[[1]][, number] - pred_slope_se) / variable
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
            type = 'response',
            xlim = c(-125.7, -116.5),
            ylim = range(species_subset$lat, na.rm = TRUE) + c(-.4, .5),
            family = "serif",
            xlab = "Longitude",
            ylab = yaxis,
            main = size,
            cex.lab = 4.5,
            cex.axis =  4.2,
            cex.main = 4.8)
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
          bg = alpha('goldenrod', 0.4),
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
       cex = 4,
       family = "serif")
  image.plot(legend.only = T,
             col = jet.colors(100),
             legend.shrink = 0.2,
             smallplot = c(.24, .29, .08, .21),
             legend.cex = 2.5,
             axis.args = list(cex.axis = 3.5,
                              family = "serif"),
             legend.width = 0.8,
             legend.mar = 6,
             zlim = c(min(my_gam$linear.predictors),
                      max(my_gam$linear.predictors)),
             legend.args = list("log(CPUE+1)",
                                side = 2, 
                                cex = 2.5,
                                family = "serif",
                                line = 1.5))
}

plot_var_coef2 <- function(my_gam, species_subset, predictions, yaxis, size){
  myvis_gam(my_gam,
            view = c('lon', 'lat'),
            too.far = 0.07,
            plot.type = 'contour',
            contour.col = contour_col,
            color = "jet" ,
            type = 'response',
            xlim = c(-125.7, -116.5),
            ylim = range(species_subset$lat, na.rm = TRUE) + c(-.4, .5),
            family = "serif",
            xlab = "Longitude",
            ylab = yaxis,
            main = size,
            cex.lab = 4.5,
            cex.axis =  4.2,
            cex.main = 4.8)
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
          bg = alpha('goldenrod', 0.4),
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
       cex = 4,
       family = "serif")
  image.plot(legend.only = T,
             col = jet.colors(100),
             legend.shrink = 0.2,
             smallplot = c(.24, .29, .08, .21),
             legend.cex = 2.5,
             axis.args = list(cex.axis = 3.5,
                              family = "serif"),
             legend.width = 0.8,
             legend.mar = 6,
             zlim = c(min(my_gam$linear.predictors),
                      max(my_gam$linear.predictors)),
             legend.args = list("log(CPUE+1)",
                                side = 2, 
                                cex = 2.5,
                                family = "serif",
                                line = 1.5))
}

plot_var_coef3 <- function(my_gam, species_subset, predictions, yaxis, size){
  myvis_gam(my_gam,
            view = c('lon', 'lat'),
            too.far = 0.07,
            plot.type = 'contour',
            contour.col = contour_col,
            color = "jet" ,
            type = 'response',
            xlim = c(-125.7, -116.5),
            ylim = range(species_subset$lat, na.rm = TRUE) + c(-.4, .5),
            family = "serif",
            xlab = "Longitude",
            ylab = yaxis,
            main = size,
            cex.lab = 4.5,
            cex.axis =  4.2,
            cex.main = 4.8)
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
       cex = 4,
       family = "serif")
  image.plot(legend.only = T,
             col = jet.colors(100),
             legend.shrink = 0.2,
             smallplot = c(.24, .29, .08, .21),
             legend.cex = 2.5,
             axis.args = list(cex.axis = 3.5,
                              family = "serif"),
             legend.width = 0.8,
             legend.mar = 6,
             zlim = c(min(my_gam$linear.predictors),
                      max(my_gam$linear.predictors)),
             legend.args = list("log(CPUE+1)",
                                side = 2, 
                                cex = 2.5,
                                family = "serif",
                                line = 1.5))
}

RMSE_calc <- function(results, data){
  RMSE <- lapply(results, function(x) {
    rmse(x$catch1, x$pred)
  })
  
  range(data$catch1)
  avg <- mean(unlist(RMSE)) 
  
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
  avg <- mean(unlist(RMSE)) 
  
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
  avg <- mean(unlist(RMSE)) 
  
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
  gams <- lapply(unique(data$year_f), function(x) {
    output <- gam(formula,
                  family = tw(link = "log"),
                  method = 'REML',
                  data = data[data$year_f != x,])
  }) # Gives the list of GAMs with each year left out
  return(gams)
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

RMSE_plot <- function(data, title){
  ggplot(data) +
    geom_line(aes(year, RMSE),
              size = 1.3,
              group = 1,
              color = "maroon4") +
    labs(x = "Year",
         y = "RMSE",
         title = title) +
    # scale_x_discrete(breaks = c(1987, 1997, 2007, 2017)) +
    theme_tufte() +
    theme(axis.title = element_text(size = 26,
                                    color = "maroon4",
                                    family = "serif"),
          axis.text = element_text(size = 22,
                                   family = "serif"),
          plot.title = element_text(size = 28, 
                                    face = "bold",
                                    family = "serif"),
          axis.line = element_line(color = "black"))
}

plot_variable <- function(gam, covariate, bounds, variable, ylabel, yvalues){
  plot(gam,
       pages = 0,
       select = covariate, 
       shade = T,
       shade.col = alpha('deepskyblue4', 0.4),
       ylim = bounds,
       xlab = variable,
       ylab = ylabel,
       yaxt = yvalues,
       seWithMean = T,
       scale = 0,
       cex.axis = 6.5,
       cex.lab = 6.5,
       family = "serif",
       lwd = 2.5)
}


# Load data ----
# See function for modifications made
# Currently only have hindcast up to 2018, so data is filtered to this
yoy_hake <- read_data('yoy_hake.Rdata') 
yoy_anchovy <- read_data('yoy_anch.Rdata') 
yoy_anchovy <- filter(yoy_anchovy, year > 2013 & jday < 164)
yoy_widow <- read_data('yoy_widw.Rdata')
yoy_widow <- filter(yoy_widow, catch < 2000) # two large hauls in 2016 caused huge errors
yoy_shortbelly <- read_data('yoy_sbly.Rdata') 
yoy_sdab <- read_data('yoy_dab.Rdata') 

state_labels <- data.frame(name = c("Washington", "Oregon", "California"),
                           lat = c(47, 44.0, 37.0),
                           lon = c(-121.0, -121.0, -120.0))

# Hake ----
# Aggregate model
# Use models selected during model exploration
hake_total <- gam(catch1 ~ year_f +
                    s(lon, lat) +
                    s(bottom_depth, k = 5) +
                    s(jday) +
                    s(roms_temperature, k = 5) +
                    s(roms_salinity, k = 5) + 
                    s(roms_ssh, k = 5) +
                    s(lon, lat, by = ssh_pos),
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
                          s(bottom_depth, k = 5) +
                          s(jday) +
                          s(roms_temperature, k = 5) +
                          s(roms_salinity, k = 5) +
                          s(roms_ssh, k = 5) +
                          s(lon, lat, by = ssh_pos)) # Note no factor(year), added into response

hake_gams <- LOYO_validation(yoy_hake, hake_formula)

# Create list with data from each year
hake_data <- split(yoy_hake, yoy_hake$year)

# Get predictions
# Predict on the left out year's data
hake_results <- hake_data
hake_results <- LOYO_preds(hake_gams, hake_data, hake_results)

# Calculate RMSE
# Get values for each year and overall value
hake_error <- RMSE_calc(hake_results, yoy_hake)
mean(hake_error[[2]]$RMSE) # 196

# Plot the RMSE for each year
hake_error[[2]]$roms_temperature <- yoy_hake$roms_temperature[match(hake_error[[2]]$year, yoy_hake$year)]

ggplot(hake_error[[2]]) +
  geom_line(aes(year, RMSE),
             linewidth = 1,
             group = 1) 

ggplot(hake_error[[2]]) +
  geom_line(aes(year, roms_temperature),
            linewidth = 1,
            group = 1) 


# Size explicit
# Small
# Use models selected during model exploration
hake_small <- gam(small_catch1 ~ year_f +
                    s(lon, lat) +
                    s(bottom_depth, k = 5) +
                    s(jday) +
                    s(roms_temperature, k = 5) +
                    s(roms_salinity, k = 5) +
                    s(roms_ssh, k = 5) +
                    s(lon, lat, by = ssh_pos),
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
                                s(bottom_depth, k = 5) +
                                s(jday) +
                                s(roms_temperature, k = 5) +
                                s(roms_salinity, k = 5) +
                                s(roms_ssh, k = 5) +
                                s(lon, lat, by = ssh_pos)) # Note no year factor, added into response

hake_small_gams <- LOYO_validation(yoy_hake, hake_small_formula)

# Get predictions
# Predict on the left out year's data
hake_small_results <- hake_data
hake_small_results <- LOYO_preds_small(hake_small_gams, hake_data, hake_small_results)

# Calculate RMSE
# Get values for each year and overall value
hake_small_error <- RMSE_calc_small(hake_small_results, yoy_hake)

# Plot the RMSE for each year
hake_small_error[[2]]$roms_salinity <- yoy_hake$roms_salinity[match(hake_small_error[[2]]$year, yoy_hake$year)]

ggplot(hake_small_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(hake_small_error[[2]]) +
  geom_line(aes(year, roms_salinity),
            size = 1,
            group = 1)

# Large
# Use models selected during model exploration
hake_large <- gam(large_catch1 ~ year_f +
                    s(lon, lat) +
                    s(bottom_depth, k = 5) +
                    s(jday) +
                    s(roms_temperature, k = 5) +
                    s(roms_salinity, k = 5) +
                    s(roms_ssh, k = 5) +
                    s(lon, lat, by = ssh_pos),
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
                                s(bottom_depth, k = 5) +
                                s(jday) +
                                s(roms_temperature, k = 5) +
                                s(roms_salinity, k = 5) +
                                s(roms_ssh, k = 5) +
                                s(lon, lat, by = ssh_pos)) # Note no year factor, added into response

hake_large_gams <- LOYO_validation(yoy_hake, hake_large_formula)

# Get predictions
# Predict on the left out year's data
hake_large_results <- hake_data
hake_large_results <- LOYO_preds_large(hake_large_gams, hake_data, hake_large_results)

# Calculate RMSE
# Get values for each year and overall value
hake_large_error <- RMSE_calc_large(hake_large_results, yoy_hake)

# Plot the RMSE for each year
hake_large_error[[2]]$roms_ssh <- yoy_hake$roms_ssh[match(hake_large_error[[2]]$year, yoy_hake$year)]

ggplot(hake_large_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(hake_large_error[[2]]) +
  geom_line(aes(year, roms_ssh),
            size = 1,
            group = 1)

# Add the large and small together
hake_combined_results <- map2(hake_large_results, 
                              hake_small_results, 
                              ~left_join(.x, .y))

hake_added_results <- lapply(hake_combined_results, function(x){
  rmse(x$catch1, x$pred_small + x$pred_large)
})

mean(unlist(hake_added_results)) # 235

hake_combined_df <- data.frame(year = names(hake_added_results), 
                               RMSE = unlist(hake_added_results))

windows(width = 10,
        height = 8)
RMSE_plot(hake_combined_df, 
          "Yearly Error for Pacific Hake")
dev.copy(jpeg, 
         here('results/hindcast_output/yoy_hake', 
              'hake_explicit_RMSE.jpg'), 
         height = 8, 
         width = 10, 
         units = 'in',
         res = 200)
dev.off()

# Partial dependence plots
# Aggregate model
tiff(here('results/hindcast_output/yoy_hake',
          'hake_partial_dependence.jpg'),
     units = "in",
     width = 56,
     height = 12,
     res = 200)
par(mfrow = c(1, 5),
    mar = c(11, 15, .5, 0.6) + 0.1,
    oma = c(3, 1, 1, 1),
    mgp = c(9, 4, 0))
plot_variable(hake_total,
              covariate = 2,
              bounds = c(-4.5, 6),
              "Depth",
              "Effect on Species Abundance",
              "s")
plot_variable(hake_total,
              covariate = 4,
              bounds = c(-4.5, 6),
              "Temperature",
              " ",
              "n")
plot_variable(hake_total,
              covariate = 5,
              bounds = c(-4.5, 6),
              "Salinity",
              " ",
              "n")
plot_variable(hake_total,
              covariate = 6,
              bounds = c(-4.5, 6),
              "Sea Surface Height",
              " ",
              "n")
plot_variable(hake_total,
              covariate = 3,
              bounds = c(-4.5, 6),
              "Day of Year",
              " ",
              "n")
dev.off()

# Small model
tiff(here('results/hindcast_output/yoy_hake',
          'hake_partial_dependence_small.jpg'),
     units = "in",
     width = 56,
     height = 12,
     res = 200)
par(mfrow = c(1, 5),
    mar = c(11, 15, .5, 0.6) + 0.1,
    oma = c(3, 1, 1, 1),
    mgp = c(9, 4, 0))
plot_variable(hake_small,
              covariate = 2,
              bounds = c(-2, 4),
              "Depth",
              "Effect on Species Abundance",
              "s")
plot_variable(hake_small,
              covariate = 4,
              bounds = c(-2, 4),
              "Temperature",
              " ",
              "n")
plot_variable(hake_small,
              covariate = 5,
              bounds = c(-2, 4),
              "Salinity",
              " ",
              "n")
plot_variable(hake_small,
              covariate = 6,
              bounds = c(-4, 2),
              "Sea Surface Height",
              " ",
              "n")
plot_variable(hake_small,
              covariate = 3,
              bounds = c(-2, 4),
              "Day of Year",
              " ",
              "n")
dev.off()

# Large model
tiff(here('results/hindcast_output/yoy_hake',
          'hake_partial_dependence_large.jpg'),
     units = "in",
     width = 56,
     height = 12,
     res = 200)
par(mfrow = c(1, 5),
    mar = c(11, 15, .5, 0.6) + 0.1,
    oma = c(3, 1, 1, 1),
    mgp = c(9, 4, 0))
plot_variable(hake_large,
              covariate = 2,
              bounds = c(-7.5, 7.5),
              "Depth",
              "Effect on Species Abundance",
              "s")
plot_variable(hake_large,
              covariate = 4,
              bounds = c(-7.5, 7.5),
              "Temperature",
              " ",
              "n")
plot_variable(hake_large,
              covariate = 5,
              bounds = c(-7.5, 7.5),
              "Salinity",
              " ",
              "n")
plot_variable(hake_large,
              covariate = 6,
              bounds = c(-7.5, 7.5),
              "Sea Surface Height",
              " ",
              "n")
plot_variable(hake_large,
              covariate = 3,
              bounds = c(-7.5, 7.5),
              "Day of Year",
              " ",
              "n")
dev.off()

# Maps
# General distributions
par(mfrow = c(1, 3),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
location_plot(hake_total, yoy_hake, "Latitude", "All Sizes", log(yoy_hake$catch1))
location_plot(hake_large, yoy_hake, " ", "Small Sizes (7-35 mm)", log(yoy_hake$large_catch1))
location_plot(hake_small, yoy_hake, " ", "Large Sizes (36-134 mm)", log(yoy_hake$small_catch1))
dev.copy(jpeg, here('results/hindcast_output/yoy_hake', 
                    'hake_distributions.jpg'), 
         height = 15, 
         width = 20, 
         units = 'in', 
         res = 200)
dev.off()

# Variable coefficient plots
pred_hake_all <- variable_coefficient(hake_total, yoy_hake, yoy_hake$ssh_pos, 8)
pred_hake_small <- variable_coefficient(hake_small, yoy_hake, yoy_hake$ssh_pos, 8)
pred_hake_large <- variable_coefficient(hake_large, yoy_hake, yoy_hake$ssh_pos, 8)

windows()
par(mfrow = c(1, 3),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
plot_var_coef3(hake_total, yoy_hake, pred_hake_all, "Latitude", "All Sizes")
plot_var_coef(hake_small, yoy_hake, pred_hake_small, "", "Small Sizes (7-35 mm)")
plot_var_coef(hake_large, yoy_hake, pred_hake_large, "", "Large Sizes (36-134 mm)")
dev.copy(jpeg, here('results/hindcast_output/yoy_hake', 
                    'yoy_hake_var_coef.jpg'), 
         height = 15, 
         width = 20, 
         units = 'in', 
         res = 200)
dev.off()

# Anchovy ----
# Aggregate model
# Use models selected during model exploration
anchovy_total <- gam(catch1 ~ year_f +
                       s(lon, lat) +
                       s(bottom_depth, k = 5) +
                       s(jday) +
                       s(roms_temperature, k = 5) + 
                       # s(roms_salinity, k = 5) +
                       s(roms_ssh, k = 5) +
                       s(lon, lat, by = ssh_pos),
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
                             s(bottom_depth, k = 5) +
                             s(jday) +
                             s(roms_temperature, k = 5) + 
                             # s(roms_salinity, k = 5) +
                             s(roms_ssh, k = 5) +
                             s(lon, lat, by = ssh_pos)) # Note no factor(year), added into response

anchovy_gams <- LOYO_validation(yoy_anchovy, anchovy_formula)

# Create list with data from each year
anchovy_data <- split(yoy_anchovy, yoy_anchovy$year)

# Get predictions
# Predict on the left out year's data
anchovy_results <- anchovy_data
anchovy_results <- LOYO_preds(anchovy_gams, anchovy_data, anchovy_results)

# Calculate RMSE
# Get values for each year and overall value
anchovy_error <- RMSE_calc(anchovy_results, yoy_anchovy)
mean(anchovy_error[[2]]$RMSE) # 330

# Plot the RMSE for each year
anchovy_error[[2]]$roms_temperature <- yoy_anchovy$roms_temperature[match(anchovy_error[[2]]$year, yoy_anchovy$year)]

ggplot(anchovy_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(anchovy_error[[2]]) +
  geom_line(aes(year, roms_temperature),
            size = 1,
            group = 1) 


# Size explicit
# Small
# Use models selected during model exploration
anchovy_small <- gam(small_catch1 ~ year_f +
                       s(lon, lat) +
                       s(bottom_depth, k = 5) +
                       s(jday) +
                       s(roms_temperature, k = 5) +
                       s(roms_salinity, k = 5) +
                       s(roms_ssh, k = 5) +
                       s(lon, lat, by = ssh_pos),
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
                                   s(bottom_depth, k = 5) +
                                   s(jday) +
                                   s(roms_temperature, k = 5) +
                                   s(roms_salinity, k = 5) +
                                   s(roms_ssh, k = 5) +
                                   s(lon, lat, by = ssh_pos)) # Note no year factor, added into response

anchovy_small_gams <- LOYO_validation(yoy_anchovy, anchovy_small_formula)

# Get predictions
# Predict on the left out year's data
anchovy_small_results <- anchovy_data
anchovy_small_results <- LOYO_preds_small(anchovy_small_gams, anchovy_data, anchovy_small_results)

# Calculate RMSE
# Get values for each year and overall value
anchovy_small_error <- RMSE_calc_small(anchovy_small_results, yoy_anchovy)

# Plot the RMSE for each year
anchovy_small_error[[2]]$roms_salinity <- yoy_anchovy$roms_salinity[match(anchovy_small_error[[2]]$year, yoy_anchovy$year)]

ggplot(anchovy_small_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(anchovy_small_error[[2]]) +
  geom_line(aes(year, roms_salinity),
            size = 1,
            group = 1)

# Large
# Use models selected during model exploration
anchovy_large <- gam(large_catch1 ~ year_f +
                       s(lon, lat) +
                       s(bottom_depth, k = 5) +
                       s(jday) +
                       # s(roms_temperature, k = 5) +
                       s(roms_salinity, k = 5) +
                       s(roms_ssh, k = 5) +
                       s(lon, lat, by = ssh_pos),
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
                                   s(bottom_depth, k = 5) +
                                   s(jday) +
                                   # s(roms_temperature, k = 5) +
                                   s(roms_salinity, k = 5) +
                                   s(roms_ssh, k = 5) +
                                   s(lon, lat, by = ssh_pos)) # Note no year factor, added into response

anchovy_large_gams <- LOYO_validation(yoy_anchovy, anchovy_large_formula)

# Get predictions
# Predict on the left out year's data
anchovy_large_results <- anchovy_data
anchovy_large_results <- LOYO_preds_large(anchovy_large_gams, anchovy_data, anchovy_large_results)

# Calculate RMSE
# Get values for each year and overall value
anchovy_large_error <- RMSE_calc_large(anchovy_large_results, yoy_anchovy)

# Plot the RMSE for each year
anchovy_large_error[[2]]$roms_ssh <- yoy_anchovy$roms_ssh[match(anchovy_large_error[[2]]$year, yoy_anchovy$year)]

ggplot(anchovy_large_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(anchovy_large_error[[2]]) +
  geom_line(aes(year, roms_ssh),
            size = 1,
            group = 1)

# Add the large and small together
anchovy_combined_results <- map2(anchovy_large_results, 
                              anchovy_small_results, 
                              ~left_join(.x, .y))

anchovy_added_results <- lapply(anchovy_combined_results, function(x){
  rmse(x$catch1, x$pred_small + x$pred_large)
})

mean(unlist(anchovy_added_results)) # 298

anchovy_combined_df <- data.frame(year = names(anchovy_added_results), 
                                  RMSE = unlist(anchovy_added_results))

windows(width = 10,
        height = 8)
RMSE_plot(anchovy_combined_df, 
          "Yearly Error for Northern Anchovy")
dev.copy(jpeg, 
         here('results/hindcast_output/yoy_anchovy', 
              'anchovy_explicit_RMSE.jpg'), 
         height = 8, 
         width = 10, 
         units = 'in',
         res = 200)
dev.off()

# Partial dependence plots
# Aggregate model
tiff(here('results/hindcast_output/yoy_anchovy',
          'anchovy_partial_dependence.jpg'),
     units = "in",
     width = 56,
     height = 12,
     res = 200)
par(mfrow = c(1, 5),
    mar = c(11, 15, .5, 0.6) + 0.1,
    oma = c(3, 1, 1, 1),
    mgp = c(9, 4, 0))
plot_variable(anchovy_total,
              covariate = 2,
              bounds = c(-3, 6),
              "Depth",
              "Effect on Species Abundance",
              "s")
plot_variable(anchovy_total,
              covariate = 4,
              bounds = c(-3, 6),
              "Temperature",
              " ",
              "n")
plot(0, xaxt = "n", yaxt = "n", bty = "n", pch = "", ylab = "", xlab = "")
plot_variable(anchovy_total,
              covariate = 5,
              bounds = c(-3, 6),
              "Sea Surface Height",
              " ",
              "n")
plot_variable(anchovy_total,
              covariate = 3,
              bounds = c(-3, 6),
              "Day of Year",
              " ",
              "n")
dev.off()

# Small model
tiff(here('results/hindcast_output/yoy_anchovy',
          'anchovy_partial_dependence_small.jpg'),
     units = "in",
     width = 56,
     height = 12,
     res = 200)
par(mfrow = c(1, 5),
    mar = c(11, 15, .5, 0.6) + 0.1,
    oma = c(3, 1, 1, 1),
    mgp = c(9, 4, 0))
plot_variable(anchovy_small,
              covariate = 2,
              bounds = c(-4.5, 5),
              "Depth",
              "Effect on Species Abundance",
              "s")
plot_variable(anchovy_small,
              covariate = 4,
              bounds = c(-4.5, 5),
              "Temperature",
              " ",
              "n")
plot_variable(anchovy_small,
              covariate = 5,
              bounds = c(-4.5, 5),
              "Salinity",
              " ",
              "n")
plot_variable(anchovy_small,
              covariate = 6,
              bounds = c(-4.5, 5),
              "Sea Surface Height",
              " ",
              "n")
plot_variable(anchovy_small,
              covariate = 3,
              bounds = c(-4.5, 5),
              "Day of Year",
              " ",
              "n")
dev.off()

# Large model
tiff(here('results/hindcast_output/yoy_anchovy',
          'anchovy_partial_dependence_large.jpg'),
     units = "in",
     width = 56,
     height = 12,
     res = 200)
par(mfrow = c(1, 5),
    mar = c(11, 15, .5, 0.6) + 0.1,
    oma = c(3, 1, 1, 1),
    mgp = c(9, 4, 0))
plot_variable(anchovy_large,
              covariate = 2,
              bounds = c(-3, 6.2),
              "Depth",
              "Effect on Species Abundance",
              "s")
plot(0, xaxt = "n", yaxt = "n", bty = "n", pch = "", ylab = "", xlab = "")
plot_variable(anchovy_large,
              covariate = 4,
              bounds = c(-3, 6.2),
              "Salinity",
              " ",
              "n")
plot_variable(anchovy_large,
              covariate = 5,
              bounds = c(-3, 6.2),
              "Sea Surface Height",
              " ",
              "n")
plot_variable(anchovy_large,
              covariate = 3,
              bounds = c(-3, 6.2),
              "Day of Year",
              " ",
              "n")
dev.off()

# Maps
# General distributions
par(mfrow = c(1, 3),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
location_plot(anchovy_total, yoy_anchovy, "Latitude", "All Sizes", log(yoy_anchovy$catch1))
location_plot(anchovy_large, yoy_anchovy, " ", "Small Sizes (15-35 mm)", log(yoy_anchovy$large_catch1))
location_plot(anchovy_small, yoy_anchovy, " ", "Large Sizes (36-92 mm)", log(yoy_anchovy$small_catch1))
dev.copy(jpeg, here('results/hindcast_output/yoy_anchovy', 
                    'anchovy_distributions.jpg'), 
         height = 15, 
         width = 20, 
         units = 'in', 
         res = 200)
dev.off()

# Variable coefficient plots
pred_anchovy_all <- variable_coefficient(anchovy_total, yoy_anchovy, yoy_anchovy$ssh_pos, 7)
pred_anchovy_small <- variable_coefficient(anchovy_small, yoy_anchovy, yoy_anchovy$ssh_pos, 8)
pred_anchovy_large <- variable_coefficient(anchovy_large, yoy_anchovy, yoy_anchovy$ssh_pos, 7)

windows()
par(mfrow = c(1, 3),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
plot_var_coef(anchovy_total, yoy_anchovy, pred_anchovy_all, "Latitude", "All Sizes")
plot_var_coef(anchovy_small, yoy_anchovy, pred_anchovy_small, "", "Small Sizes (15-35 mm)")
plot_var_coef(anchovy_large, yoy_anchovy, pred_anchovy_large, "", "Large Sizes (36-92 mm)")
dev.copy(jpeg, here('results/hindcast_output/yoy_anchovy', 
                    'yoy_anchovy_var_coef.jpg'), 
         height = 15, 
         width = 20, 
         units = 'in', 
         res = 200)
dev.off()

# Widow Rockfish ----
# Aggregate model
# Use models selected during model exploration
widow_total <- gam(catch1 ~ year_f +
                     s(lon, lat) +
                     s(bottom_depth, k = 5) +
                     s(jday) +
                     s(roms_temperature, k = 5) + 
                     # s(roms_salinity, k = 5) + 
                     s(roms_ssh, k = 5) +
                     s(lon, lat, by = ssh_pos),
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
                           s(bottom_depth, k = 5) +
                           s(jday) +
                           s(roms_temperature, k = 5) +
                           # s(roms_salinity, k = 5) +
                           s(roms_ssh, k = 5) +
                           s(lon, lat, by = ssh_pos)) # Note no factor(year), added into response

widow_gams <- LOYO_validation(yoy_widow, widow_formula)

# Create list with data from each year
widow_data <- split(yoy_widow, yoy_widow$year)

# Get predictions
# Predict on the left out year's data
widow_results <- widow_data
widow_results <- LOYO_preds(widow_gams, widow_data, widow_results)

# Calculate RMSE
# Get values for each year and overall value
widow_error <- RMSE_calc(widow_results, yoy_widow)
mean(widow_error[[2]]$RMSE) # 29

# Plot the RMSE for each year
widow_error[[2]]$roms_temperature <- yoy_widow$roms_temperature[match(widow_error[[2]]$year, yoy_widow$year)]

ggplot(widow_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(widow_error[[2]]) +
  geom_line(aes(year, roms_temperature),
            size = 1,
            group = 1) 


# Size explicit
# Small
# Use models selected during model exploration
widow_small <- gam(small_catch1 ~ year_f +
                     s(lon, lat) +
                     s(bottom_depth, k = 5) +
                     s(jday) +
                     # s(roms_temperature, k = 5) + 
                     s(roms_salinity, k = 5) +
                     s(roms_ssh, k = 5) +
                     s(lon, lat, by = ssh_pos),
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
                                 s(bottom_depth, k = 5) +
                                 s(jday) +
                                 # s(roms_temperature, k = 5) +
                                 s(roms_salinity, k = 5) +
                                 s(roms_ssh, k = 5) +
                                 s(lon, lat,  by = ssh_pos)) # Note no year factor, added into response

widow_small_gams <- LOYO_validation(yoy_widow, widow_small_formula)

# Get predictions
# Predict on the left out year's data
widow_small_results <- widow_data
widow_small_results <- LOYO_preds_small(widow_small_gams, widow_data, widow_small_results)

# Calculate RMSE
# Get values for each year and overall value
widow_small_error <- RMSE_calc_small(widow_small_results, yoy_widow)

# Plot the RMSE for each year
widow_small_error[[2]]$roms_salinity <- yoy_widow$roms_salinity[match(widow_small_error[[2]]$year, yoy_widow$year)]

ggplot(widow_small_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(widow_small_error[[2]]) +
  geom_line(aes(year, roms_salinity),
            size = 1,
            group = 1)

# Large
# Use models selected during model exploration
widow_large <- gam(large_catch1 ~ year_f +
                     s(lon, lat) +
                     s(bottom_depth, k = 5) +
                     s(jday) +
                     s(roms_temperature, k = 5) +
                     s(roms_salinity, k = 5) +
                     s(roms_ssh, k = 5) +
                     s(lon, lat, by = ssh_pos),
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
                                 s(bottom_depth, k = 5) +
                                 s(jday) +
                                 s(roms_temperature, k = 5) +
                                 s(roms_salinity, k = 5) +
                                 s(roms_ssh, k = 5) +
                                 s(lon, lat,  by = ssh_pos)) # Note no year factor, added into response

widow_large_gams <- LOYO_validation(yoy_widow, widow_large_formula)

# Get predictions
# Predict on the left out year's data
widow_large_results <- widow_data
widow_large_results <- LOYO_preds_large(widow_large_gams, widow_data, widow_large_results)

# Calculate RMSE
# Get values for each year and overall value
widow_large_error <- RMSE_calc_large(widow_large_results, yoy_widow)

# Plot the RMSE for each year
widow_large_error[[2]]$roms_ssh <- yoy_widow$roms_ssh[match(widow_large_error[[2]]$year, yoy_widow$year)]

ggplot(widow_large_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(widow_large_error[[2]]) +
  geom_line(aes(year, roms_ssh),
            size = 1,
            group = 1)

# Add the large and small together
widow_combined_results <- map2(widow_large_results, 
                              widow_small_results, 
                              ~left_join(.x, .y))

widow_added_results <- lapply(widow_combined_results, function(x){
  rmse(x$catch1, x$pred_small + x$pred_large)
})

mean(unlist(widow_added_results)) # 29

widow_combined_df <- data.frame(year = names(widow_added_results), 
                               RMSE = unlist(widow_added_results))

windows(width = 10,
        height = 8)
RMSE_plot(widow_combined_df, 
          "Yearly Error for Widow Rockfish")
dev.copy(jpeg, 
         here('results/hindcast_output/yoy_widow', 
              'widow_explicit_RMSE.jpg'), 
         height = 8, 
         width = 10, 
         units = 'in',
         res = 200)
dev.off()

# Partial dependence plots
# Aggregate model
tiff(here('results/hindcast_output/yoy_widow',
          'widow_partial_dependence.jpg'),
     units = "in",
     width = 42,
     height = 12,
     res = 200)
par(mfrow = c(1, 4),
    mar = c(11, 15, .5, 0.6) + 0.1,
    oma = c(3, 1, 1, 1),
    mgp = c(9, 4, 0))
plot_variable(widow_total,
              covariate = 2,
              bounds = c(-1.5, 1.2),
              "Depth",
              "Effect on Species Abundance",
              "s")
plot_variable(widow_total,
              covariate = 4,
              bounds = c(-1.5, 1.2),
              "Temperature",
              " ",
              "n")
plot_variable(widow_total,
              covariate = 5,
              bounds = c(-1.5, 1.2),
              "Sea Surface Height",
              " ",
              "n")
plot_variable(widow_total,
              covariate = 3,
              bounds = c(-1.5, 1.2),
              "Day of Year",
              " ",
              "n")
dev.off()

# Small model
tiff(here('results/hindcast_output/yoy_widow',
          'widow_partial_dependence_small.jpg'),
     units = "in",
     width = 42,
     height = 12,
     res = 200)
par(mfrow = c(1, 4),
    mar = c(11, 15, .5, 0.6) + 0.1,
    oma = c(3, 1, 1, 1),
    mgp = c(9, 4, 0))
plot_variable(widow_small,
              covariate = 2,
              bounds = c(-2.5, 1),
              "Depth",
              "Effect on Species Abundance",
              "s")
plot_variable(widow_small,
              covariate = 4,
              bounds = c(-2.5, 1),
              "Temperature",
              " ",
              "n")
plot_variable(widow_small,
              covariate = 5,
              bounds = c(-2.5, 1),
              "Sea Surface Height",
              " ",
              "n")
plot_variable(widow_small,
              covariate = 3,
              bounds = c(-2.5, 1),
              "Day of Year",
              " ",
              "n")
dev.off()

# Large model
tiff(here('results/hindcast_output/yoy_widow',
          'widow_partial_dependence_large.jpg'),
     units = "in",
     width = 56,
     height = 12,
     res = 200)
par(mfrow = c(1, 5),
    mar = c(11, 15, .5, 0.6) + 0.1,
    oma = c(3, 1, 1, 1),
    mgp = c(9, 4, 0))
plot_variable(widow_large,
              covariate = 2,
              bounds = c(-1, 1.5),
              "Depth",
              "Effect on Species Abundance",
              "s")
plot_variable(widow_large,
              covariate = 4,
              bounds = c(-1, 1.5),
              "Temperature",
              " ",
              "n")
plot_variable(widow_large,
              covariate = 5,
              bounds = c(-1, 1.5),
              "Salinity",
              " ",
              "n")
plot_variable(widow_large,
              covariate = 6,
              bounds = c(-1, 1.5),
              "Sea Surface Height",
              " ",
              "n")
plot_variable(widow_large,
              covariate = 3,
              bounds = c(-1, 1.5),
              "Day of Year",
              " ",
              "n")
dev.off()

# Maps
# General distributions
par(mfrow = c(1, 3),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
location_plot(widow_total, yoy_widow, "Latitude", "All Sizes", log(yoy_widow$catch1))
location_plot(widow_large, yoy_widow, " ", "Small Sizes (11-45 mm)", log(yoy_widow$large_catch1))
location_plot(widow_small, yoy_widow, " ", "Large Sizes (46-73 mm)", log(yoy_widow$small_catch1))
dev.copy(jpeg, here('results/hindcast_output/yoy_widow', 
                    'widow_distributions.jpg'), 
         height = 15, 
         width = 20, 
         units = 'in', 
         res = 200)
dev.off()

# Variable coefficient plots
pred_widow_all <- variable_coefficient(widow_total, yoy_widow, yoy_widow$ssh_pos, 7)
pred_widow_small <- variable_coefficient(widow_small, yoy_widow, yoy_widow$ssh_pos, 7)
pred_widow_large <- variable_coefficient(widow_large, yoy_widow, yoy_widow$ssh_pos, 8)

windows()
par(mfrow = c(1, 3),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
plot_var_coef(widow_total, yoy_widow, pred_widow_all, "Latitude", "All Sizes")
plot_var_coef3(widow_small, yoy_widow, pred_widow_small, "", "Small Sizes (11-45 mm)")
plot_var_coef(widow_large, yoy_widow, pred_widow_large, "", "Large Sizes (46-73 mm)")
dev.copy(jpeg, here('results/hindcast_output/yoy_widow', 
                    'yoy_widow_var_coef.jpg'), 
         height = 15, 
         width = 20, 
         units = 'in', 
         res = 200)
dev.off()

# Shortbelly Rockfish ----
# Aggregate model
# Use models selected during model exploration
shortbelly_total <- gam(catch1 ~ year_f +
                          s(lon, lat) +
                          s(bottom_depth, k = 5) +
                          s(jday) +
                          s(roms_temperature, k = 5) +
                          s(roms_salinity, k = 5) +
                          # s(roms_ssh, k = 5) + # removal reduced AIC
                          s(lon, lat, by = ssh_pos),
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
                                s(bottom_depth, k = 5) +
                                s(jday) +
                                s(roms_temperature, k = 5) +
                                s(roms_salinity, k = 5) +
                                # s(roms_ssh, k = 5) +
                                s(lon, lat, by = ssh_pos)) # Note no factor(year), added into response

shortbelly_gams <- LOYO_validation(yoy_shortbelly, shortbelly_formula)

# Create list with data from each year
shortbelly_data <- split(yoy_shortbelly, yoy_shortbelly$year)

# Get predictions
# Predict on the left out year's data
shortbelly_results <- shortbelly_data
shortbelly_results <- LOYO_preds(shortbelly_gams, shortbelly_data, shortbelly_results)

# Calculate RMSE
# Get values for each year and overall value
shortbelly_error <- RMSE_calc(shortbelly_results, yoy_shortbelly)
mean(shortbelly_error[[2]]$RMSE) # 157

# Plot the RMSE for each year
shortbelly_error[[2]]$roms_temperature <- yoy_shortbelly$roms_temperature[match(shortbelly_error[[2]]$year, yoy_shortbelly$year)]

ggplot(shortbelly_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(shortbelly_error[[2]]) +
  geom_line(aes(year, roms_temperature),
            size = 1,
            group = 1) 


# Size explicit
# Small
# Use models selected during model exploration
shortbelly_small <- gam(small_catch1 ~ year_f +
                          s(lon, lat) +
                          s(bottom_depth, k = 5) +
                          s(jday) +
                          s(roms_temperature, k = 5) +
                          s(roms_salinity, k = 5) +
                          s(roms_ssh, k = 5) +
                          s(lon, lat, by = ssh_pos),
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
                                      s(bottom_depth, k = 5) +
                                      s(jday) +
                                      s(roms_temperature, k = 5) +
                                      s(roms_salinity, k = 5) +
                                      s(roms_ssh, k = 5) +
                                      s(lon, lat, by = ssh_pos)) # Note no year factor, added into response

shortbelly_small_gams <- LOYO_validation(yoy_shortbelly, shortbelly_small_formula)

# Get predictions
# Predict on the left out year's data
shortbelly_small_results <- shortbelly_data
shortbelly_small_results <- LOYO_preds_small(shortbelly_small_gams, shortbelly_data, shortbelly_small_results)

# Calculate RMSE
# Get values for each year and overall value
shortbelly_small_error <- RMSE_calc_small(shortbelly_small_results, yoy_shortbelly)

# Plot the RMSE for each year
shortbelly_small_error[[2]]$roms_salinity <- yoy_shortbelly$roms_salinity[match(shortbelly_small_error[[2]]$year, yoy_shortbelly$year)]

ggplot(shortbelly_small_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(shortbelly_small_error[[2]]) +
  geom_line(aes(year, roms_salinity),
            size = 1,
            group = 1)

# Large
# Use models selected during model exploration
shortbelly_large <- gam(large_catch1 ~ year_f +
                          s(lon, lat) +
                          s(bottom_depth, k = 5) +
                          s(jday) +
                          s(roms_temperature, k = 5) +
                          # s(roms_salinity, k = 5) +
                          s(roms_ssh, k = 5) +
                          s(lon, lat, by = ssh_pos),
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
                                      s(bottom_depth, k = 5) +
                                      s(jday) +
                                      s(roms_temperature, k = 5) +
                                      # s(roms_salinity, k = 5) +
                                      s(roms_ssh, k = 5) +
                                      s(lon, lat, by = ssh_pos)) # Note no year factor, added into response

shortbelly_large_gams <- LOYO_validation(yoy_shortbelly, shortbelly_large_formula)

# Get predictions
# Predict on the left out year's data
shortbelly_large_results <- shortbelly_data
shortbelly_large_results <- LOYO_preds_large(shortbelly_large_gams, shortbelly_data, shortbelly_large_results)

# Calculate RMSE
# Get values for each year and overall value
shortbelly_large_error <- RMSE_calc_large(shortbelly_large_results, yoy_shortbelly)

# Plot the RMSE for each year
shortbelly_large_error[[2]]$roms_ssh <- yoy_shortbelly$roms_ssh[match(shortbelly_large_error[[2]]$year, yoy_shortbelly$year)]

ggplot(shortbelly_large_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(shortbelly_large_error[[2]]) +
  geom_line(aes(year, roms_ssh),
            size = 1,
            group = 1)

# Add the large and small together
shortbelly_combined_results <- map2(shortbelly_large_results, 
                                    shortbelly_small_results, 
                                    ~left_join(.x, .y))

shortbelly_added_results <- lapply(shortbelly_combined_results, function(x){
  rmse(x$catch1, x$pred_small + x$pred_large)
})

mean(unlist(shortbelly_added_results)) # 154

shortbelly_combined_df <- data.frame(year = names(shortbelly_added_results), 
                                     RMSE = unlist(shortbelly_added_results))

windows(width = 10,
        height = 8)
RMSE_plot(shortbelly_combined_df, 
          "Yearly Error for Shortbelly Rockfish")
dev.copy(jpeg, 
         here('results/hindcast_output/yoy_shortbelly', 
              'shortbelly_explicit_RMSE.jpg'), 
         height = 8, 
         width = 10, 
         units = 'in',
         res = 200)
dev.off()

# Partial dependence plots
# Aggregate model
tiff(here('results/hindcast_output/yoy_shortbelly',
          'shortbelly_partial_dependence.jpg'),
     units = "in",
     width = 56,
     height = 12,
     res = 200)
par(mfrow = c(1, 5),
    mar = c(11, 15, .5, 0.6) + 0.1,
    oma = c(3, 1, 1, 1),
    mgp = c(9, 4, 0))
plot_variable(shortbelly_total,
              covariate = 2,
              bounds = c(-3.5, 2.5),
              "Depth",
              "Effect on Species Abundance",
              "s")
plot_variable(shortbelly_total,
              covariate = 4,
              bounds = c(-3.5, 2.5),
              "Temperature",
              " ",
              "n")
plot_variable(shortbelly_total,
              covariate = 5,
              bounds = c(-3.5, 2.5),
              "Salinity",
              " ",
              "n")
plot_variable(shortbelly_total,
              covariate = 6,
              bounds = c(-3.5, 2.5),
              "Sea Surface Height",
              " ",
              "n")
plot_variable(shortbelly_total,
              covariate = 3,
              bounds = c(-3.5, 2.5),
              "Day of Year",
              " ",
              "n")
dev.off()

# Small model
tiff(here('results/hindcast_output/yoy_shortbelly',
          'shortbelly_partial_dependence_small.jpg'),
     units = "in",
     width = 56,
     height = 12,
     res = 150)
par(mfrow = c(1, 5),
    mar = c(11, 15, .5, 0.6) + 0.1,
    oma = c(3, 1, 1, 1),
    mgp = c(9, 4, 0))
plot_variable(shortbelly_small,
              covariate = 2,
              bounds = c(-4.5, 2.5),
              "Depth",
              "Effect on Species Abundance",
              "s")
plot_variable(shortbelly_small,
              covariate = 4,
              bounds = c(-4.5, 2.5),
              "Temperature",
              " ",
              "n")
plot_variable(shortbelly_small,
              covariate = 5,
              bounds = c(-4.5, 2.5),
              "Salinity",
              " ",
              "n")
plot_variable(shortbelly_small,
              covariate = 6,
              bounds = c(-4.5, 2.5),
              "Sea Surface Height",
              " ",
              "n")
plot_variable(shortbelly_small,
              covariate = 3,
              bounds = c(-4.5, 2.5),
              "Day of Year",
              " ",
              "n")
dev.off()

# Large model
tiff(here('results/hindcast_output/yoy_shortbelly',
          'shortbelly_partial_dependence_large.jpg'),
     units = "in",
     width = 42,
     height = 12,
     res = 200)
par(mfrow = c(1, 4),
    mar = c(11, 15, .5, 0.6) + 0.1,
    oma = c(3, 1, 1, 1),
    mgp = c(9, 4, 0))
plot_variable(shortbelly_large,
              covariate = 2,
              bounds = c(-4.5, 1),
              "Depth",
              "Effect on Species Abundance",
              "s")
plot_variable(shortbelly_large,
              covariate = 4,
              bounds = c(-4.5, 1),
              "Salinity",
              " ",
              "n")
plot_variable(shortbelly_large,
              covariate = 5,
              bounds = c(-4.5, 1),
              "Sea Surface Height",
              " ",
              "n")
plot_variable(shortbelly_large,
              covariate = 3,
              bounds = c(-4.5, 1),
              "Day of Year",
              " ",
              "n")
dev.off()

# Maps
# General distributions
par(mfrow = c(1, 3),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
location_plot(shortbelly_total, yoy_shortbelly, "Latitude", "All Sizes", log(yoy_shortbelly$catch1))
location_plot(shortbelly_large, yoy_shortbelly, " ", "Small Sizes (8-35 mm)", log(yoy_shortbelly$large_catch1))
location_plot(shortbelly_small, yoy_shortbelly, " ", "Large Sizes (36-85 mm)", log(yoy_shortbelly$small_catch1))
dev.copy(jpeg, here('results/hindcast_output/yoy_shortbelly', 
                    'shortbelly_distributions.jpg'), 
         height = 15, 
         width = 20, 
         units = 'in', 
         res = 200)
dev.off()

# Variable coefficient plots
pred_shortbelly_all <- variable_coefficient(shortbelly_total, yoy_shortbelly, yoy_shortbelly$ssh_pos, 7)
pred_shortbelly_small <- variable_coefficient(shortbelly_small, yoy_shortbelly, yoy_shortbelly$ssh_pos, 6)
pred_shortbelly_large <- variable_coefficient(shortbelly_large, yoy_shortbelly, yoy_shortbelly$ssh_pos, 6)

windows()
par(mfrow = c(1, 3),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
plot_var_coef3(shortbelly_total, yoy_shortbelly, pred_shortbelly_all, "Latitude", "All Sizes")
plot_var_coef3(shortbelly_small, yoy_shortbelly, pred_shortbelly_small, "", "Small Sizes (8-35 mm)")
plot_var_coef(shortbelly_large, yoy_shortbelly, pred_shortbelly_large, "", "Large Sizes (36-85 mm)")
dev.copy(jpeg, here('results/hindcast_output/yoy_shortbelly', 
                    'yoy_shortbelly_var_coef.jpg'), 
         height = 15, 
         width = 20, 
         units = 'in', 
         res = 200)
dev.off()

# Pacific Sanddab ----
# Aggregate model
# Use models selected during model exploration
sdab_total <- gam(catch1 ~ year_f + 
                    s(lon, lat) + 
                    s(bottom_depth, k = 5) +
                    s(jday) + 
                    # s(roms_temperature, k = 5) +
                    s(roms_salinity, k = 5) +
                    s(roms_ssh, k = 5) +
                    s(lon, lat, by = ssh_pos),
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
                          s(bottom_depth, k = 5) +
                          s(jday) + 
                          # s(roms_temperature, k = 5) +
                          s(roms_salinity, k = 5) +
                          s(roms_ssh, k = 5) +
                          s(lon, lat, by = ssh_pos)) # Note no factor(year), added into response

sdab_gams <- LOYO_validation(yoy_sdab, sdab_formula)

# Create list with data from each year
sdab_data <- split(yoy_sdab, yoy_sdab$year)

# Get predictions
# Predict on the left out year's data
sdab_results <- sdab_data
sdab_results <- LOYO_preds(sdab_gams, sdab_data, sdab_results)

# Calculate RMSE
# Get values for each year and overall value
sdab_error <- RMSE_calc(sdab_results, yoy_sdab)
mean(sdab_error[[2]]$RMSE) # 134

# Plot the RMSE for each year
sdab_error[[2]]$roms_temperature <- yoy_sdab$roms_temperature[match(sdab_error[[2]]$year, yoy_sdab$year)]

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
  geom_line(aes(year, roms_temperature),
            size = 1,
            group = 1) 


# Size explicit
# Small
# Use models selected during model exploration
sdab_small <- gam(small_catch1 ~ year_f +
                    s(lon, lat) +
                    s(bottom_depth, k = 5) +
                    s(jday) +
                    s(roms_temperature, k = 5) +
                    s(roms_salinity, k = 5) +
                    s(roms_ssh, k = 5) +
                    s(lon, lat, by = ssh_pos),
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
                                s(bottom_depth, k = 5) +
                                s(jday) +
                                s(roms_temperature, k = 5) +
                                s(roms_salinity, k = 5) +
                                s(roms_ssh, k = 5) +
                                s(lon, lat, by = ssh_pos)) # Note no year factor, added into response

sdab_small_gams <- LOYO_validation(yoy_sdab, sdab_small_formula)

# Get predictions
# Predict on the left out year's data
sdab_small_results <- sdab_data
sdab_small_results <- LOYO_preds_small(sdab_small_gams, sdab_data, sdab_small_results)

# Calculate RMSE
# Get values for each year and overall value
sdab_small_error <- RMSE_calc_small(sdab_small_results, yoy_sdab)

# Plot the RMSE for each year
sdab_small_error[[2]]$roms_salinity <- yoy_sdab$roms_salinity[match(sdab_small_error[[2]]$year, yoy_sdab$year)]

ggplot(sdab_small_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(sdab_small_error[[2]]) +
  geom_line(aes(year, roms_salinity),
            size = 1,
            group = 1)

# Large
# Use models selected during model exploration
sdab_large <- gam(large_catch1 ~ year_f +
                    s(lon, lat) +
                    s(bottom_depth, k = 5) +
                    s(jday) +
                    # s(roms_temperature, k = 5) +
                    # s(roms_salinity, k = 5) +
                    s(roms_ssh, k = 5) +
                    s(lon, lat, by = ssh_pos),
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
                                s(bottom_depth, k = 5) +
                                s(jday) +
                                # s(roms_temperature, k = 5) +
                                # s(roms_salinity, k = 5) +
                                s(roms_ssh, k = 5) +
                                s(lon, lat, by = ssh_pos)) # Note no year factor, added into response

sdab_large_gams <- LOYO_validation(yoy_sdab, sdab_large_formula)

# Get predictions
# Predict on the left out year's data
sdab_large_results <- sdab_data
sdab_large_results <- LOYO_preds_large(sdab_large_gams, sdab_data, sdab_large_results)

# Calculate RMSE
# Get values for each year and overall value
sdab_large_error <- RMSE_calc_large(sdab_large_results, yoy_sdab)

# Plot the RMSE for each year
sdab_large_error[[2]]$roms_ssh <- yoy_sdab$roms_ssh[match(sdab_large_error[[2]]$year, yoy_sdab$year)]

ggplot(sdab_large_error[[2]]) +
  geom_line(aes(year, RMSE),
            size = 1,
            group = 1) 

ggplot(sdab_large_error[[2]]) +
  geom_line(aes(year, roms_ssh),
            size = 1,
            group = 1)

# Add the large and small together
sdab_combined_results <- map2(sdab_large_results, 
                              sdab_small_results, 
                              ~left_join(.x, .y))

sdab_added_results <- lapply(sdab_combined_results, function(x){
  rmse(x$catch1, x$pred_small + x$pred_large)
})

mean(unlist(sdab_added_results)) # 169

sdab_combined_df <- data.frame(year = names(sdab_added_results), 
                               RMSE = unlist(sdab_added_results))

windows(width = 10,
        height = 8)
RMSE_plot(sdab_combined_df, 
          "Yearly Error for Pacific Sanddab")
dev.copy(jpeg, 
         here('results/hindcast_output/yoy_sanddab', 
              'sanddab_explicit_RMSE.jpg'), 
         height = 8, 
         width = 10, 
         units = 'in',
         res = 200)
dev.off()

# Partial dependence plots
# Aggregate model
tiff(here('results/hindcast_output/yoy_sanddab',
          'sanddab_partial_dependence.jpg'),
     units = "in",
     width = 56,
     height = 12,
     res = 200)
par(mfrow = c(1, 5),
    mar = c(11, 15, .5, 0.6) + 0.1,
    oma = c(3, 1, 1, 1),
    mgp = c(9, 4, 0))
plot_variable(sdab_total,
              covariate = 2,
              bounds = c(-2, 1.5),
              "Depth",
              "Effect on Species Abundance",
              "s")
plot_variable(sdab_total,
              covariate = 4,
              bounds = c(-2, 1.5),
              "Temperature",
              " ",
              "n")
plot_variable(sdab_total,
              covariate = 5,
              bounds = c(-2, 1.5),
              "Salinity",
              " ",
              "n")
plot_variable(sdab_total,
              covariate = 6,
              bounds = c(-2, 1.5),
              "Sea Surface Height",
              " ",
              "n")
plot_variable(sdab_total,
              covariate = 3,
              bounds = c(-2, 1.5),
              "Day of Year",
              " ",
              "n")
dev.off()

# Small model
tiff(here('results/hindcast_output/yoy_sanddab',
          'sanddab_partial_dependence_small.jpg'),
     units = "in",
     width = 56,
     height = 12,
     res = 150)
par(mfrow = c(1, 5),
    mar = c(11, 15, .5, 0.6) + 0.1,
    oma = c(3, 1, 1, 1),
    mgp = c(9, 4, 0))
plot_variable(sdab_small,
              covariate = 2,
              bounds = c(-3, 1.5),
              "Depth",
              "Effect on Species Abundance",
              "s")
plot_variable(sdab_small,
              covariate = 4,
              bounds = c(-3, 1.5),
              "Temperature",
              " ",
              "n")
plot_variable(sdab_small,
              covariate = 5,
              bounds = c(-3, 1.5),
              "Salinity",
              " ",
              "n")
plot_variable(sdab_small,
              covariate = 6,
              bounds = c(-3, 1.5),
              "Sea Surface Height",
              " ",
              "n")
plot_variable(sdab_small,
              covariate = 3,
              bounds = c(-3, 1.5),
              "Day of Year",
              " ",
              "n")
dev.off()

# Large model
tiff(here('results/hindcast_output/yoy_sanddab',
          'sanddab_partial_dependence_large.jpg'),
     units = "in",
     width = 42,
     height = 12,
     res = 200)
par(mfrow = c(1, 4),
    mar = c(11, 15, .5, 0.6) + 0.1,
    oma = c(3, 1, 1, 1),
    mgp = c(9, 4, 0))
plot_variable(sdab_large,
              covariate = 2,
              bounds = c(-2.2, 2.7),
              "Depth",
              "Effect on Species Abundance",
              "s")
plot_variable(sdab_large,
              covariate = 4,
              bounds = c(-2.2, 2.7),
              "Temperature",
              " ",
              "n")
plot_variable(sdab_large,
              covariate = 5,
              bounds = c(-2.2, 2.7),
              "Salinity",
              " ",
              "n")
plot_variable(sdab_large,
              covariate = 3,
              bounds = c(-2.2, 2.7),
              "Day of Year",
              " ",
              "n")
dev.off()

# Maps
# General distributions
par(mfrow = c(1, 3),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
location_plot(sdab_total, yoy_sdab, "Latitude", "All Sizes", log(yoy_sdab$catch1))
location_plot(sdab_large, yoy_sdab, " ", "Small Sizes (11-30 mm)", log(yoy_sdab$large_catch1))
location_plot(sdab_small, yoy_sdab, " ", "Large Sizes (31-82 mm)", log(yoy_sdab$small_catch1))
dev.copy(jpeg, here('results/hindcast_output/yoy_sanddab', 
                    'sanddab_distributions.jpg'), 
         height = 15, 
         width = 20, 
         units = 'in', 
         res = 200)
dev.off()

# Variable coefficient plots
pred_sdab_all <- variable_coefficient(sdab_total, yoy_sdab, yoy_sdab$ssh_pos, 7) 
pred_sdab_small <- variable_coefficient(sdab_small, yoy_sdab, yoy_sdab$ssh_pos, 8)
pred_sdab_large <- variable_coefficient(sdab_large, yoy_sdab, yoy_sdab$ssh_pos, 6)

windows()
par(mfrow = c(1, 3),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0))
plot_var_coef3(sdab_total, yoy_sdab, pred_sdab_all, "Latitude", "All Sizes")
plot_var_coef(sdab_small, yoy_sdab, pred_sdab_small, "", "Small Sizes (11-30 mm)")
plot_var_coef3(sdab_large, yoy_sdab, pred_sdab_large, "", "Large Sizes (31-82 mm)")
dev.copy(jpeg, here('results/hindcast_output/yoy_sanddab', 
                    'yoy_sanddab_var_coef.jpg'), 
         height = 15, 
         width = 20, 
         units = 'in', 
         res = 200)
dev.off()
