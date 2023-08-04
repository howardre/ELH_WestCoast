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
  sdmTMB(catch ~ 0 + as.factor(year) +
           s(bottom_depth, k = 5) +
           s(roms_temperature, k = 5) +
           s(roms_salinity, k = 5) +
           s(jday),
         spatial_varying = ~ 0 + ssh_pos, # Not sure this is the right way to do this with SSH
         data = df,
         mesh = mesh,
         time = "year",
         spatial = "on",
         family = tweedie(link = "log"),
         spatiotemporal = "ar1",
         control = sdmTMBcontrol(newton_loops = 1))
}

sdmTMB_small <- function(df, mesh){
  sdmTMB(small ~ 0 + as.factor(year) +
           s(bottom_depth, k = 5) +
           s(roms_temperature, k = 5) +
           s(roms_salinity, k = 5) +
           s(jday),
         spatial_varying = ~ 0 + ssh_pos, 
         data = df,
         mesh = mesh,
         time = "year",
         spatial = "on",
         family = tweedie(link = "log"),
         spatiotemporal = "ar1",
         control = sdmTMBcontrol(newton_loops = 1))
}

sdmTMB_large <- function(df, mesh){
  sdmTMB(large ~ 0 + as.factor(year) +
           s(bottom_depth, k = 5) +
           s(roms_temperature, k = 5) +
           s(roms_salinity, k = 5) +
           s(jday),
         spatial_varying = ~ 0 + ssh_pos, 
         data = df,
         mesh = mesh,
         time = "year",
         spatial = "on",
         family = tweedie(link = "log"),
         spatiotemporal = "ar1",
         control = sdmTMBcontrol(newton_loops = 1))
}

# Data
yoy_hake <- read_data('yoy_hake.Rdata') 
# yoy_anchovy <- read_data('yoy_anch.Rdata') 
# yoy_anchovy <- filter(yoy_anchovy, year > 2013 & jday < 164)
# yoy_widow <- read_data('yoy_widw.Rdata')
# yoy_widow <- filter(yoy_widow, catch < 2000) # two large hauls in 2016 caused huge errors
# yoy_shortbelly <- read_data('yoy_sbly.Rdata') 
# yoy_sdab <- read_data('yoy_dab.Rdata') 

# Pacific Hake ----
# Make mesh object with matrices
yoy_hake_mesh <- make_mesh(yoy_hake, 
                           xy_cols = c("X", "Y"), 
                           n_knots = 200,
                           type = "cutoff_search",
                           seed = 1993)
plot(yoy_hake_mesh)

# Fit models
hake_model <- sdmTMB_formula(yoy_hake, 
                             yoy_hake_mesh)
hake_model_small <- sdmTMB_small(yoy_hake,
                                 yoy_hake_mesh)
hake_model_large <- sdmTMB_large(yoy_hake,
                                 yoy_hake_mesh)

sanity(hake_model)
tidy(hake_model)
tidy(hake_model,
     effect = "ran_pars", 
     conf.int = T)
sanity(hake_model_small)
tidy(hake_model_small)
tidy(hake_model_small,
     effect = "ran_pars", 
     conf.int = T)
sanity(hake_model_large)
tidy(hake_model_large)
tidy(hake_model_large,
     effect = "ran_pars", 
     conf.int = T)

# Leave-one-year-out cross validation
hake_results <- LOYO_validation(yoy_hake) # Test to see if all converge

# Plot covariates
visreg(hake_model, 
       xvar = "bottom_depth", 
       scale = "response")

# Predict and plot
hake_pred <- sdmTMB_grid(yoy_hake, hake_model)
sdmTMB_map(yoy_hake, hake_pred)

hake_pred_small <- sdmTMB_grid(yoy_hake, hake_model_small)
sdmTMB_map(yoy_hake, hake_pred_small)

hake_pred_large <- sdmTMB_grid(yoy_hake, hake_model_large)
sdmTMB_map(yoy_hake, hake_pred_large)