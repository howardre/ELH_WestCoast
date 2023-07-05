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

# Load data and functions ----
# Functions
source(here('code/functions', 'vis_gam_COLORS.R'))
read_data <- function(file){
  yoy <- readRDS(here('data', file)) %>% 
    tidyr::drop_na(roms_temperature, roms_salinity, roms_ssh, bottom_depth, year, jday, lat, lon) %>%
    filter(catch < 2500 &
             year < 2020 &
             lat < 42) %>%
    mutate(catch1 = catch + 1,
           small_catch1 = small + 1,
           large_catch1 = large + 1,
           year_f = as.factor(year),
           ssh_pos = year_ssh + abs(min(year_ssh)) + 10)
  yoy <- yoy[!(yoy$small == 0 & yoy$large == 0 & yoy$catch > 0), ]
  yoy_utm <- add_utm_columns(yoy, c("lon", "lat")) # add UTM coordinates
  return(yoy_utm)
}

# Data
yoy_hake <- read_data('yoy_hake.Rdata') 
yoy_anchovy <- read_data('yoy_anch.Rdata') 
yoy_anchovy <- filter(yoy_anchovy, year > 2013 & jday < 164)
yoy_widow <- read_data('yoy_widw.Rdata')
yoy_widow <- filter(yoy_widow, catch < 2000) # two large hauls in 2016 caused huge errors
yoy_shortbelly <- read_data('yoy_sbly.Rdata') 
yoy_sdab <- read_data('yoy_dab.Rdata') 

# Make mesh object with matrices
yoy_hake_mesh <- make_mesh(yoy_hake, 
                           xy_cols = c("X", "Y"), 
                           cutoff = 10,
                           seed = 1993)
plot(yoy_hake_mesh)

# Fit model
# Can add k-fold cross validation
hake_model <- sdmTMB(catch ~ 0 + as.factor(year) +
                       s(bottom_depth, k = 5) +
                       s(roms_temperature, k = 5) +
                       s(roms_salinity, k = 5) +
                       s(jday),
                     spatial_varying = ~ ssh_pos,
                     data = yoy_hake,
                     mesh = yoy_hake_mesh,
                     time = "year",
                     family = tweedie(link = "log"),
                     spatiotemporal = "off",
                     control = sdmTMBcontrol(newton_loops = 1))
sanity(hake_model)

# Plot
hake_model
tidy(hake_model)
tidy(hake_model,
     effect = "ran_pars", 
     conf.int = T)

visreg(hake_model, xvar = "bottom_depth", scale = "response")

# Variable coefficient
hake_vc_grid <- replicate_df(yoy_hake, "ssh_pos", unique(yoy_hake$ssh_pos))
hake_vc_grid$scaled_ssh <- (hake_vc_grid$ssh_pos - mean(yoy_hake$ssh_pos)) / sd(yoy_hake$ssh_pos)
hake_pred <- as.data.frame(predict(hake_model,
                     newdata = hake_vc_grid))

plot_map_raster <- function(dat, column = est) {
  ggplot(dat, aes(X, Y, fill = {{ column }})) +
    geom_tile() +
    facet_wrap(~ ssh_pos) +
    coord_fixed() +
    scale_fill_viridis_c()
}

plot_map_raster(filter(hake_pred, ssh_pos == 10), zeta_s_ssh_pos)
