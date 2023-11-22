### Title: sdmTMB Model Development
### Author: Rebecca Howard
### Date: 06/30/2023

# Load libraries ----
library(readxl)
library(tibble)
library(here)
library(ggplot2)
library(dplyr)
library(sdmTMB) # do not update TMB past 1.9.4
library(visreg)
library(raster)
library(sf)
library(colorspace)
library(mapdata)
library(fields)
library(ggpubr)

# Load data and functions ----
# Functions
source(here('code/functions', 'vis_gam_COLORS.R'))
source(here('code/functions', 'distance_function.R'))
source(here('code/functions', 'sdmTMB_grid.R'))
source(here('code/functions', 'sdmTMB_map.R'))
source(here('code/functions', 'read_data.R'))

# Functions
# More info here: https://pbs-assess.github.io/sdmTMB/index.html
LOYO_validation <- function(df){
  models <- lapply(unique(df$year), function(x) {
    the_mesh <- make_mesh(df[df$year != x, ], 
                          xy_cols = c("X", "Y"), 
                          cutoff = 10,
                          seed = 1993)
    output <- sdmTMB(catch ~ 0 + as.factor(year) +
                       s(bottom_depth, k = 5) +
                       s(roms_temperature, k = 5) +
                       s(roms_salinity, k = 5) +
                       s(jday),
                     mesh = the_mesh,
                     time = "year",
                     family = tweedie(link = "log"),
                     data = df[df$year != x, ])
  }) # Gives the list of models with each year left out
  return(models)
}

sdmTMB_formula <- function(df, mesh){
  sdmTMB(catch ~ 0 +
           s(bottom_depth, k = 5) +
           s(roms_temperature, k = 5) +
           s(roms_salinity, k = 5) +
           s(ssh_anom, k = 5) +
           s(jday),
       # spatial_varying = ~ 0 + ssh_pos, # change to new variables
         data = df,
         mesh = mesh,
         time = "year",
         spatial = "on",
         family = tweedie(link = "log"),
         spatiotemporal = "ar1",
         control = sdmTMBcontrol(newton_loops = 1))
}

sdmTMB_small <- function(df, mesh){
  sdmTMB(small ~ 0 +
           s(bottom_depth, k = 5) +
           s(roms_temperature, k = 5) +
           s(roms_salinity, k = 5) +
           s(ssh_anom, k = 5) +
           s(jday),
#         spatial_varying = ~ 0 + ssh_pos, 
         data = df,
         mesh = mesh,
         time = "year",
         spatial = "on",
         family = tweedie(link = "log"),
         spatiotemporal = "ar1",
         control = sdmTMBcontrol(newton_loops = 1))
}

sdmTMB_large <- function(df, mesh){
  sdmTMB(large ~ 0 +
           s(bottom_depth, k = 5) +
           s(roms_temperature, k = 5) +
           s(roms_salinity, k = 5) +
           s(ssh_anom, k = 5) +
           s(jday),
       #  spatial_varying = ~ 0 + ssh_pos, 
         data = df,
         mesh = mesh,
         time = "year",
         spatial = "on",
         family = tweedie(link = "log"),
         spatiotemporal = "ar1",
         control = sdmTMBcontrol(newton_loops = 1))
}

# Data
# Need to fix length issues (allocating all to large sizes if no measurements)
yoy_hake <- read_data('yoy_hake.Rdata') 
yoy_anchovy <- read_data('yoy_anch.Rdata')
yoy_anchovy <- filter(yoy_anchovy, jday < 164)
yoy_widow <- read_data('yoy_widw.Rdata')
yoy_widow <- filter(yoy_widow, catch < 2000) # two large hauls in 2016 caused huge errors
yoy_shortbelly <- read_data('yoy_sbly.Rdata')
yoy_sdab <- read_data('yoy_dab.Rdata')

# Pacific Hake ----
# Make mesh object with matrices
yoy_hake_mesh <- make_mesh(yoy_hake, 
                           xy_cols = c("X", "Y"), 
                           n_knots = 200, # increase back to 200
                           type = "cutoff_search",
                           seed = 1993)
plot(yoy_hake_mesh)

# Fit models
# Currently having issues with spatially varying term
# May improve with the new variables?
hake_model_small <- sdmTMB_small(yoy_hake,
                                 yoy_hake_mesh) # currently takes 6 minutes
hake_model_large <- sdmTMB_large(yoy_hake,
                                 yoy_hake_mesh)

sanity(hake_model_small) # sigma_z is the SD of the spatially varying coefficient field
tidy(hake_model_small,
     effect = "ran_pars", 
     conf.int = T)
sanity(hake_model_large)
tidy(hake_model_large,
     effect = "ran_pars", 
     conf.int = T)

# Leave-one-year-out cross validation
hake_results <- LOYO_validation(yoy_hake) # Test to see if all converge

# Plot covariates
individual_plot <- function(output, variable, xlab){
  ggplot(output$fit, aes(x = variable, y = visregFit)) +
    geom_line(color = "black",
              linewidth = 1,
              show.legend = FALSE) +
    geom_ribbon(aes(ymin = visregLwr, 
                    ymax = visregUpr,
                    fill = "coral2"),
                alpha = 0.5,
                show.legend = FALSE) +
    labs(x = xlab,
         y = "Abundance Anomalies") 
}

plot_variables <- function(model, data){
  depth <- visreg(model,
                  data = data,
                  xvar = "bottom_depth",
                  plot = FALSE)
  depth_plot <- ggplot(depth$fit, aes(x = bottom_depth, y = visregFit)) +
    geom_line(color = "black",
              linewidth = 1,
              show.legend = FALSE) +
    geom_ribbon(aes(ymin = visregLwr, 
                    ymax = visregUpr,
                    fill = "coral2"),
                alpha = 0.5,
                show.legend = FALSE) +
    labs(x = 'Depth (m)',
         y = "Abundance Anomalies") 

  temp <- visreg(model,
                 data = data,
                 xvar = "roms_temperature",
                 plot = FALSE)
  temp_plot <- ggplot(temp$fit, aes(x = roms_temperature, y = visregFit)) +
    geom_line(color = "black",
              linewidth = 1,
              show.legend = FALSE) +
    geom_ribbon(aes(ymin = visregLwr, 
                    ymax = visregUpr,
                    fill = "coral2"),
                alpha = 0.5,
                show.legend = FALSE) +
    labs(x = 'Temperature (\u00B0C)',
         y = "Abundance Anomalies") 

  salt <- visreg(model,
                 data = data,
                 xvar = "roms_salinity",
                 plot = FALSE)
  salt_plot <- ggplot(salt$fit, aes(x = roms_salinity, y = visregFit)) +
    geom_line(color = "black",
              linewidth = 1,
              show.legend = FALSE) +
    geom_ribbon(aes(ymin = visregLwr, 
                    ymax = visregUpr,
                    fill = "coral2"),
                alpha = 0.5,
                show.legend = FALSE) +
    labs(x = 'Salinity',
         y = "Abundance Anomalies") 

  ssh <- visreg(model,
                 data = data,
                 xvar = "ssh_anom",
                 plot = FALSE)
  ssh_plot <- ggplot(ssh$fit, aes(x = ssh_anom, y = visregFit)) +
    geom_line(color = "black",
              linewidth = 1,
              show.legend = FALSE) +
    geom_ribbon(aes(ymin = visregLwr, 
                    ymax = visregUpr,
                    fill = "coral2"),
                alpha = 0.5,
                show.legend = FALSE) +
    labs(x = 'Sea Surface Height',
         y = "Abundance Anomalies") 

  doy <- visreg(model,
                data = data,
                xvar = "jday",
                plot = FALSE)
  doy_plot <- ggplot(doy$fit, aes(x = jday, y = visregFit)) +
    geom_line(color = "black",
              linewidth = 1,
              show.legend = FALSE) +
    geom_ribbon(aes(ymin = visregLwr, 
                    ymax = visregUpr,
                    fill = "coral2"),
                alpha = 0.5,
                show.legend = FALSE) +
    labs(x = 'Day of Year',
         y = "Abundance Anomalies") 

  ggarrange(depth_plot, temp_plot, salt_plot, ssh_plot, doy_plot) +
    theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text = element_text(family = "serif", size = 18),
          axis.title = element_text(family = "serif", size = 22),
          axis.text.x = element_text(angle = 0, vjust = 0.7))
}

plot_variables(hake_model, yoy_hake)
plot_variables(hake_model_small, yoy_hake)
plot_variables(hake_model_large, yoy_hake)


# Predict and plot
nlat = 40
nlon = 60
latd = seq(min(yoy_hake$lat), max(yoy_hake$lat), length.out = nlat)
lond = seq(min(yoy_hake$lon), max(yoy_hake$lon), length.out = nlon)

hake_pred <- sdmTMB_grid(yoy_hake, hake_model)
hake_pred_small <- sdmTMB_grid(yoy_hake, hake_model_small)
hake_pred_large <- sdmTMB_grid(yoy_hake, hake_model_large)


windows(height = 15, width = 24)
par(mfrow = c(1, 3),
    mar = c(6.4, 7.2, 1.6, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
sdmTMB_map(yoy_hake, hake_pred)
sdmTMB_map(yoy_hake, hake_pred_small)
sdmTMB_map(yoy_hake, hake_pred_large)