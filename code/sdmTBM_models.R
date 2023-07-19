### Title: sdmTMB Model Development
### Author: Rebecca Howard
### Date: 06/30/2023

# Load libraries ----
library(readxl)
library(tibble)
library(here)
library(ggplot2)
library(dplyr)
library(sdmTMB)
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
source(here('code/functions', 'sdmTBM_grid.R'))
source(here('code/functions', 'sdmTBM_map.R'))
source(here('code/functions', 'read_data.R'))

# Data
yoy_hake <- read_data('yoy_hake.Rdata') 
# yoy_anchovy <- read_data('yoy_anch.Rdata') 
# yoy_anchovy <- filter(yoy_anchovy, year > 2013 & jday < 164)
# yoy_widow <- read_data('yoy_widw.Rdata')
# yoy_widow <- filter(yoy_widow, catch < 2000) # two large hauls in 2016 caused huge errors
# yoy_shortbelly <- read_data('yoy_sbly.Rdata') 
# yoy_sdab <- read_data('yoy_dab.Rdata') 

# Models


# Make mesh object with matrices
yoy_hake_mesh <- make_mesh(yoy_hake, 
                           xy_cols = c("X", "Y"), 
                           cutoff = 10,
                           seed = 1993)
plot(yoy_hake_mesh)

# Fit model
# More info here: https://pbs-assess.github.io/sdmTMB/index.html
# Can add k-fold cross validation
hake_model <- sdmTMB(catch ~ 0 + as.factor(year) +
                       s(bottom_depth, k = 5) +
                       s(roms_temperature, k = 5) +
                       s(roms_salinity, k = 5) +
                       s(jday) ,
                     # spatial_varying = ~ 0 + ssh_pos, # Not sure this is the right way to do this with SSH
                     data = yoy_hake,
                     mesh = yoy_hake_mesh,
                     time = "year",
                     family = tweedie(link = "log"),
                     spatiotemporal = "rw")
sanity(hake_model)

tidy(hake_model)
tidy(hake_model,
     effect = "ran_pars", 
     conf.int = T)

# Plot covariates
visreg(hake_model, 
       xvar = "bottom_depth", 
       scale = "response")

# Predict and plot
hake_pred <- sdmTMB_grid(yoy_hake, hake_model)
sdmTMB_map(yoy_hake, hake_pred)
