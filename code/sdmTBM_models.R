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
library(visreg) # still waiting on v2.7.1 to fix SVC issue
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

sdmTMB_small <- function(df, mesh){
  sdmTMB(small ~ 0 + 
           s(bottom_depth, k = 5) +
           s(roms_temperature, k = 5) +
           s(roms_salinity, k = 5) +
           s(ssh_anom, k = 5) +
           s(jday, k = 15),
         spatial_varying = ~ 0 + ssh_annual_scaled,
         data = df,
         mesh = mesh,
         time = "year",
         spatial = "on",
         family = tweedie(link = "log"),
         spatiotemporal = "iid",
         control = sdmTMBcontrol(newton_loops = 1,
                                 nlminb_loops = 2))
}

sdmTMB_large <- function(df, mesh){
  sdmTMB(large ~ 0 + 
           s(bottom_depth, k = 5) +
           s(roms_temperature, k = 5) +
           s(roms_salinity, k = 5) +
           s(ssh_anom, k = 5) +
           s(jday, k = 15),
         spatial_varying = ~ 0 + ssh_annual_scaled,
         data = df,
         mesh = mesh,
         time = "year",
         spatial = "on",
         family = tweedie(link = "log"),
         spatiotemporal = "iid",
       control = sdmTMBcontrol(newton_loops = 1,
                               nlminb_loops = 2))
}

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
                  plot = FALSE,
                  cond = list(ssh_annual_scaled = 1)) # currently required with SVCs
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
         y = "Abundance Anomalies") +
    theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text = element_text(family = "serif", size = 18),
          axis.title = element_text(family = "serif", size = 22),
          axis.text.x = element_text(angle = 0, vjust = 0.7)) 
  
  temp <- visreg(model,
                 data = data,
                 xvar = "roms_temperature",
                 plot = FALSE,
                 cond = list(ssh_annual_scaled = 1))
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
         y = "Abundance Anomalies") +
    theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text = element_text(family = "serif", size = 18),
          axis.title = element_text(family = "serif", size = 22),
          axis.text.x = element_text(angle = 0, vjust = 0.7)) 
  
  salt <- visreg(model,
                 data = data,
                 xvar = "roms_salinity",
                 plot = FALSE,
                 cond = list(ssh_annual_scaled = 1))
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
         y = "Abundance Anomalies") +
    theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text = element_text(family = "serif", size = 18),
          axis.title = element_text(family = "serif", size = 22),
          axis.text.x = element_text(angle = 0, vjust = 0.7)) 
  
  ssh <- visreg(model,
                data = data,
                xvar = "ssh_anom",
                plot = FALSE,
                cond = list(ssh_annual_scaled = 1))
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
         y = "Abundance Anomalies") +
    theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text = element_text(family = "serif", size = 18),
          axis.title = element_text(family = "serif", size = 22),
          axis.text.x = element_text(angle = 0, vjust = 0.7)) 
  
  doy <- visreg(model,
                data = data,
                xvar = "jday",
                plot = FALSE,
                cond = list(ssh_annual_scaled = 1))
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
         y = "Abundance Anomalies") +
    theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text = element_text(family = "serif", size = 18),
          axis.title = element_text(family = "serif", size = 22),
          axis.text.x = element_text(angle = 0, vjust = 0.7))
  
  ggarrange(depth_plot, temp_plot, salt_plot, ssh_plot, doy_plot)
}

# Data
# May have to filter salinity and depth due to outliers?
yoy_hake <- read_data('yoy_hake.Rdata') 
yoy_anchovy <- read_data('yoy_anch.Rdata')
yoy_anchovy <- filter(yoy_anchovy, jday < 164 & survey == "RREAS")
yoy_widow <- read_data('yoy_widw.Rdata')
yoy_widow <- filter(yoy_widow, catch < 2000) # two large hauls in 2016 caused huge errors
yoy_shortbelly <- read_data('yoy_sbly.Rdata')
yoy_sdab <- read_data('yoy_dab.Rdata')

# Pacific Hake ----
# Make mesh object with matrices
yoy_hake_mesh <- make_mesh(yoy_hake, 
                           xy_cols = c("X", "Y"), 
                           n_knots = 200,
                           type = "cutoff_search",
                           seed = 2024)
plot(yoy_hake_mesh)

# Fit models
hake_model_small <- sdmTMB_small(yoy_hake,
                                 yoy_hake_mesh) 
hake_model_large <- sdmTMB_large(yoy_hake,
                                 yoy_hake_mesh) # requires model simplification

sanity(hake_model_small) # sigma_z is the SD of the spatially varying coefficient field
tidy(hake_model_small, # no std error reported when using log link
     effect = "ran_pars", 
     conf.int = TRUE)
sanity(hake_model_large) 
tidy(hake_model_large,
     effect = "ran_pars", 
     conf.int = TRUE)
hake_model_small
hake_model_large

# Plot covariates
plot_variables(hake_model_small, yoy_hake)
plot_variables(hake_model_large, yoy_hake)

# Predict and plot
nlat = 40
nlon = 60
latd = seq(min(yoy_hake$lat), max(yoy_hake$lat), length.out = nlat)
lond = seq(min(yoy_hake$lon), max(yoy_hake$lon), length.out = nlon)

hake_pred_small <- sdmTMB_grid(yoy_hake, hake_model_small)
hake_pred_large <- sdmTMB_grid(yoy_hake, hake_model_large)

windows(height = 15, width = 20)
par(mfrow = c(1, 2),
    mar = c(6.4, 7.2, 1.6, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
sdmTMB_map(yoy_hake, hake_pred_small)
sdmTMB_map(yoy_hake, hake_pred_large)


# Northern Anchovy ----
# Make mesh object with matrices
yoy_anchovy_mesh <- make_mesh(yoy_anchovy,
                              xy_cols = c("X", "Y"), 
                              n_knots = 200,
                              type = "cutoff_search",
                              seed = 2024)
plot(yoy_anchovy_mesh)

# Fit models
anchovy_model_small <- sdmTMB_small(yoy_anchovy,
                                    yoy_anchovy_mesh) # this SVC not a good fit
anchovy_model_large <- sdmTMB_large(yoy_anchovy,
                                    yoy_anchovy_mesh) # needs simplification

sanity(anchovy_model_small) 
tidy(anchovy_model_small, 
     effect = "ran_pars", 
     conf.int = T)
sanity(anchovy_model_large)
tidy(anchovy_model_large,
     effect = "ran_pars", 
     conf.int = T)
anchovy_model_small
anchovy_model_large

# Plot covariates
plot_variables(anchovy_model_small, yoy_anchovy)
plot_variables(anchovy_model_large, yoy_anchovy)

# Predict and plot
nlat = 40
nlon = 60
latd = seq(min(yoy_anchovy$lat), max(yoy_anchovy$lat), length.out = nlat)
lond = seq(min(yoy_anchovy$lon), max(yoy_anchovy$lon), length.out = nlon)

anchovy_pred_small <- sdmTMB_grid(yoy_anchovy, anchovy_model_small)
anchovy_pred_large <- sdmTMB_grid(yoy_anchovy, anchovy_model_large)

windows(height = 15, width = 20)
par(mfrow = c(1, 2),
    mar = c(6.4, 7.2, 1.6, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
sdmTMB_map(yoy_anchovy, anchovy_pred_small)
sdmTMB_map(yoy_anchovy, anchovy_pred_large)


# Pacific Sanddab ----
# Make mesh object with matrices
yoy_sdab_mesh <- make_mesh(yoy_sdab,
                           xy_cols = c("X", "Y"),
                           n_knots = 200,
                           type = "cutoff_search",
                           seed =  2024)
plot(yoy_sdab_mesh)

# Fit models
sdab_model_small <- sdmTMB_small(yoy_sdab,
                                 yoy_sdab_mesh) # this SVC not a good fit
sdab_model_large <- sdmTMB_large(yoy_sdab,
                                 yoy_sdab_mesh) # needs simplification

sanity(sdab_model_small) 
tidy(sdab_model_small, 
     effect = "ran_pars", 
     conf.int = T)
sanity(sdab_model_large)
tidy(sdab_model_large,
     effect = "ran_pars", 
     conf.int = T)
sdab_model_small
sdab_model_large

# Plot covariates
plot_variables(sdab_model_small, yoy_sdab)
plot_variables(sdab_model_large, yoy_sdab)

# Predict and plot
nlat = 40
nlon = 60
latd = seq(min(yoy_sdab$lat), max(yoy_sdab$lat), length.out = nlat)
lond = seq(min(yoy_sdab$lon), max(yoy_sdab$lon), length.out = nlon)

sdab_pred_small <- sdmTMB_grid(yoy_sdab, sdab_model_small)
sdab_pred_large <- sdmTMB_grid(yoy_sdab, sdab_model_large)

windows(height = 15, width = 20)
par(mfrow = c(1, 2),
    mar = c(6.4, 7.2, 1.6, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
sdmTMB_map(yoy_sdab, sdab_pred_small)
sdmTMB_map(yoy_sdab, sdab_pred_large)


# Shortbelly Rockfish ----
# Make mesh object with matrices
yoy_shortbelly_mesh <- make_mesh(yoy_shortbelly,
                              xy_cols = c("X", "Y"), 
                              n_knots = 200,
                              type = "cutoff_search",
                              seed = 2024)
plot(yoy_shortbelly_mesh)

# Fit models
shortbelly_model_small <- sdmTMB_small(yoy_shortbelly,
                                    yoy_shortbelly_mesh) # lots of problems
shortbelly_model_large <- sdmTMB_large(yoy_shortbelly,
                                    yoy_shortbelly_mesh) # SVC may not be good fit

sanity(shortbelly_model_small) 
tidy(shortbelly_model_small, 
     effect = "ran_pars", 
     conf.int = T)
sanity(shortbelly_model_large)
tidy(shortbelly_model_large,
     effect = "ran_pars", 
     conf.int = T)
shortbelly_model_small
shortbelly_model_large

# Plot covariates
plot_variables(shortbelly_model_small, yoy_shortbelly)
plot_variables(shortbelly_model_large, yoy_shortbelly)

# Predict and plot
nlat = 40
nlon = 60
latd = seq(min(yoy_shortbelly$lat), max(yoy_shortbelly$lat), length.out = nlat)
lond = seq(min(yoy_shortbelly$lon), max(yoy_shortbelly$lon), length.out = nlon)

shortbelly_pred_small <- sdmTMB_grid(yoy_shortbelly, shortbelly_model_small)
shortbelly_pred_large <- sdmTMB_grid(yoy_shortbelly, shortbelly_model_large)

windows(height = 15, width = 20)
par(mfrow = c(1, 2),
    mar = c(6.4, 7.2, 1.6, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
sdmTMB_map(yoy_shortbelly, shortbelly_pred_small)
sdmTMB_map(yoy_shortbelly, shortbelly_pred_large)


# Widow Rockfish ----
# Make mesh object with matrices
yoy_widow_mesh <- make_mesh(yoy_widow,
                              xy_cols = c("X", "Y"), 
                              n_knots = 200,
                              type = "cutoff_search",
                              seed = 2024)
plot(yoy_widow_mesh)

# Fit models
widow_model_small <- sdmTMB_small(yoy_widow,
                                    yoy_widow_mesh) # didn't converge
widow_model_large <- sdmTMB_large(yoy_widow,
                                    yoy_widow_mesh) # bad

sanity(widow_model_small) 
tidy(widow_model_small, 
     effect = "ran_pars", 
     conf.int = T)
sanity(widow_model_large)
tidy(widow_model_large,
     effect = "ran_pars", 
     conf.int = T)
widow_model_small
widow_model_large

# Plot covariates
plot_variables(widow_model_small, yoy_widow)
plot_variables(widow_model_large, yoy_widow)

# Predict and plot
nlat = 40
nlon = 60
latd = seq(min(yoy_widow$lat), max(yoy_widow$lat), length.out = nlat)
lond = seq(min(yoy_widow$lon), max(yoy_widow$lon), length.out = nlon)

widow_pred_small <- sdmTMB_grid(yoy_widow, widow_model_small)
widow_pred_large <- sdmTMB_grid(yoy_widow, widow_model_large)

windows(height = 15, width = 20)
par(mfrow = c(1, 2),
    mar = c(6.4, 7.2, 1.6, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
sdmTMB_map(yoy_widow, widow_pred_small)
sdmTMB_map(yoy_widow, widow_pred_large)
