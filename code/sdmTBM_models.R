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
library(RColorBrewer)
library(mapdata)
library(fields)
library(ggpubr)
library(future)

# Load data and functions ----
# Functions
source(here('code/functions', 'vis_gam_COLORS.R'))
source(here('code/functions', 'distance_function.R'))
source(here('code/functions', 'sdmTMB_grid.R'))
source(here('code/functions', 'sdmTMB_map.R'))
source(here('code/functions', 'read_data.R'))
source(here('code/functions', 'sdmTMB_select.R')) # cannot be spatiotemporal & spatial


# Functions
# More info here: https://pbs-assess.github.io/sdmTMB/index.html
individual_plot <- function(model, data, variable, xlab){
  object <- visreg(model,
                   data = data,
                   xvar = variable,
                   plot = FALSE,
                   cond = list(ssh_annual_scaled = 1)) # currently required with SVCs
  object_plot <- ggplot(object$fit, aes(x = !!sym(variable), y = visregFit)) + # !!sym() removes quotes
    geom_line(color = "black",
              linewidth = 1,
              show.legend = FALSE) +
    geom_ribbon(aes(ymin = visregLwr, 
                    ymax = visregUpr,
                    fill = "coral2"),
                alpha = 0.5,
                show.legend = FALSE) +
    labs(x = xlab,
         y = "Abundance Anomalies") +
    theme_classic() +
    theme(axis.ticks = element_blank(),
          axis.text = element_text(family = "serif", size = 38),
          axis.title = element_text(family = "serif", size = 42),
          axis.text.x = element_text(angle = 0, vjust = 0.7),
          plot.margin = margin(2, 2, 2, 2, "cm")) 
  return(object_plot)
} # use for individual variables

plot_variables <- function(model, data){
  temp <- visreg(model,
                 data = data,
                 xvar = "sst_scaled",
                 plot = FALSE)
  temp_plot <- ggplot(temp$fit, aes(x = sst_scaled, y = visregFit)) +
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
          axis.text = element_text(family = "serif", size = 38),
          axis.title = element_text(family = "serif", size = 42),
          axis.text.x = element_text(angle = 0, vjust = 0.7),
          plot.margin = margin(2, 2, 2, 2, "cm")) 
  
  salt <- visreg(model,
                 data = data,
                 xvar = "sss_scaled",
                 plot = FALSE)
  salt_plot <- ggplot(salt$fit, aes(x = sss_scaled, y = visregFit)) +
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
          axis.text = element_text(family = "serif", size = 38),
          axis.title = element_text(family = "serif", size = 42),
          axis.text.x = element_text(angle = 0, vjust = 0.7),
          plot.margin = margin(2, 2, 2, 2, "cm")) 
  doy <- visreg(model,
                data = data,
                xvar = "jday_scaled",
                plot = FALSE)
  doy_plot <- ggplot(doy$fit, aes(x = jday_scaled, y = visregFit)) +
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
          axis.text = element_text(family = "serif", size = 38),
          axis.title = element_text(family = "serif", size = 42),
          axis.text.x = element_text(angle = 0, vjust = 0.7),
          plot.margin = margin(2, 2, 2, 2, "cm"))
  
  ggarrange(temp_plot, salt_plot, doy_plot, ncol = 3, nrow = 1)
} # works only if all variables retained

# Data
# May have to filter salinity and depth due to outliers?
yoy_hake <- read_data('yoy_hake.Rdata') 
yoy_anchovy <- read_data('yoy_anch.Rdata')
yoy_widow <- filter(read_data('yoy_widw.Rdata')) # two large hauls in 2016 caused huge errors
yoy_shortbelly <- read_data('yoy_sbly.Rdata')
yoy_sdab <- filter(read_data('yoy_dab.Rdata'))
yoy_squid <- filter(read_data('yoy_squid.Rdata')) # no lengths before 2004

state_labels <- data.frame(name = c("Washington", "Oregon", "California"),
                           lat = c(47, 44.0, 37.0),
                           lon = c(-121.0, -121.0, -120.0))
nlat = 40
nlon = 60

# Sandbox ----
the_mesh <- make_mesh(yoy_widow,
                      xy_cols = c("X", "Y"),
                      cutoff = 20)
sdm_test <- sdmTMB(large ~  0 + sst_scaled +
                     sss_scaled +
                     vgeo,
                   spatial_varying = ~ 0 + vgeo,
                   data = yoy_widow,
                   mesh = the_mesh,
                   spatial = "off",
                   time = "year",
                   family = tweedie(link = "log"),
                   spatiotemporal = "ar1",
                   control = sdmTMBcontrol(newton_loops = 1,
                                           nlminb_loops = 2))
sanity(sdm_test)
tidy(sdm_test, 
     conf.int = TRUE)

# Pacific Hake ----
# Make mesh object with matrices
yoy_hake_mesh <- make_mesh(yoy_hake, 
                           xy_cols = c("X", "Y"),
                           cutoff = 18)
plot(yoy_hake_mesh) 

# Select models
hake_model_small <- sdmTMB_select_small(yoy_hake, yoy_hake_mesh) # time elapsed: 1:20
hake_model_small[[2]] # spice_iso26
hake_model_large <- sdmTMB_select_large(yoy_hake, yoy_hake_mesh) 
hake_model_large[[2]] # spice_iso26

sanity(hake_model_small[[2]]) # sigma_z is the SD of the spatially varying coefficient field
tidy(hake_model_small[[2]], # no std error reported when using log link
     conf.int = TRUE)
sanity(hake_model_large[[2]]) 
tidy(hake_model_large[[2]],
     conf.int = TRUE)

# Plot covariates
tiff(here('results/hindcast_output/yoy_hake',
          'hake_partial_dependence_small_sdmtmb.jpg'),
     units = "in",
     width = 38,
     height = 12,
     res = 200)
plot_variables(hake_model_small[[2]], yoy_hake)
dev.off()

tiff(here('results/hindcast_output/yoy_hake',
          'hake_partial_dependence_large_sdmtmb.jpg'),
     units = "in",
     width = 38,
     height = 12,
     res = 200)
plot_variables(hake_model_large[[2]], yoy_hake)
dev.off()

# Predict and plot
latd = seq(min(yoy_hake$lat), max(yoy_hake$lat), length.out = nlat)
lond = seq(min(yoy_hake$lon), max(yoy_hake$lon), length.out = nlon)

hake_pred_small <- sdmTMB_grid(yoy_hake, hake_model_small[[2]])
hake_pred_large <- sdmTMB_grid(yoy_hake, hake_model_large[[2]])

# Overall predictions
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
sdmTMB_map(yoy_hake, hake_pred_small, "Small (15-35 mm)", "Latitude")
sdmTMB_map(yoy_hake, hake_pred_large, "Large (36-81 mm)", " ")
dev.copy(jpeg, here('results/hindcast_output/yoy_hake', 
                    'hake_distributions_sdmtmb.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# SVC maps
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
sdmTMB_SVC(yoy_hake, hake_pred_small, "Small (15-35 mm)", "Latitude")
sdmTMB_SVC(yoy_hake, hake_pred_large, "Large (36-81 mm)", " ")
dev.copy(jpeg, here('results/hindcast_output/yoy_hake', 
                    'hake_SVC_sdmtmb.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Cross validation
plan(multisession)
clust <- as.numeric(as.factor(yoy_hake$year)) # year clusters, may not work with spatiotemporal "on"
hake_small_cv <- sdmTMB_cv(small ~ 0 + 
                             s(bottom_depth, k = 5) +
                             s(roms_temperature, k = 5) +
                             s(roms_salinity, k = 5) +
                             s(ssh_anom, k = 5) +
                             s(jday, k = 15),
                           spatial_varying = ~ 0 + ssh_annual_scaled,
                           data = yoy_hake,
                           mesh = yoy_hake_mesh,
                           time = "year",
                           family = tweedie(link = "log"),
                           spatiotemporal = "iid",
                           control = sdmTMBcontrol(newton_loops = 1,
                                                   nlminb_loops = 2),
                           fold_ids = clust,
                           k_folds = length(unique(clust)))

plot(hake_small_cv$fold_elpd) # higher values better
hake_small_cv$elpd

plot(hake_small_cv$fold_loglik) # higher values better
hake_small_cv$sum_loglik

clust <- kmeans(yoy_hake[, c("X", "Y")], k = 50)$cluster # spatial clustering
hake_large_cv <- sdmTMB_cv(large ~ 0 + 
                             s(bottom_depth, k = 5) +
                             s(roms_temperature, k = 5) +
                             s(roms_salinity, k = 5) +
                             s(ssh_anom, k = 5) +
                             s(jday, k = 15),
                           spatial_varying = ~ 0 + ssh_annual_scaled,
                           data = yoy_hake,
                           mesh = yoy_hake_mesh,
                           time = "year",
                           family = tweedie(link = "log"),
                           spatiotemporal = "iid",
                           control = sdmTMBcontrol(newton_loops = 1,
                                                   nlminb_loops = 2),
                           parallel = TRUE, 
                           fold_ids = clust,
                           k_folds = length(unique(clust))) # not working

# Northern Anchovy ----
# Make mesh object with matrices
yoy_anchovy_mesh <- make_mesh(yoy_anchovy,
                              xy_cols = c("X", "Y"),
                              cutoff = 20)
plot(yoy_anchovy_mesh) 

# Select models
anchovy_model_small <- sdmTMB_select_small(yoy_anchovy, yoy_anchovy_mesh) 
anchovy_model_small[[2]] # spice_iso26
anchovy_model_large <- sdmTMB_select_large(yoy_anchovy, yoy_anchovy_mesh) 
anchovy_model_large[[2]] # none

sanity(anchovy_model_small[[2]]) 
tidy(anchovy_model_small[[2]], 
     conf.int = TRUE)
sanity(anchovy_model_large[[2]]) 
tidy(anchovy_model_large[[2]],
     conf.int = TRUE)

# Plot covariates
tiff(here('results/hindcast_output/yoy_anchovy',
          'anchovy_partial_dependence_small_sdmtmb.jpg'),
     units = "in",
     width = 56,
     height = 12,
     res = 200)
plot_variables(anchovy_model_small, yoy_anchovy)
dev.off()

tiff(here('results/hindcast_output/yoy_anchovy',
          'anchovy_partial_dependence_large_sdmtmb.jpg'),
     units = "in",
     width = 56,
     height = 12,
     res = 200)
plot_variables(anchovy_model_large, yoy_anchovy)
dev.off()


# Predict and plot
latd = seq(min(yoy_anchovy$lat), max(yoy_anchovy$lat), length.out = nlat)
lond = seq(min(yoy_anchovy$lon), max(yoy_anchovy$lon), length.out = nlon)

anchovy_pred_small <- sdmTMB_grid(yoy_anchovy, anchovy_model_small)
anchovy_pred_large <- sdmTMB_grid(yoy_anchovy, anchovy_model_large)

# Overall predictions
windows(height = 15, width = 18)
par(mfrow = c(2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
sdmTMB_map(yoy_anchovy, anchovy_pred_small, "Small (15-35 mm)", "Latitude")
sdmTMB_map(yoy_anchovy, anchovy_pred_large, "Large (36-92 mm)", " ")
dev.copy(jpeg, here('results/hindcast_output/yoy_anchovy', 
                    'anchovy_distributions_sdmtmb.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# SVC maps
windows(height = 15, width = 18)
par(mfrow = c(2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
sdmTMB_SVC(yoy_anchovy, anchovy_pred_small, "Small (15-35 mm)", "Latitude")
sdmTMB_SVC(yoy_anchovy, anchovy_pred_large, "Large (36-92 mm)", " ")
dev.copy(jpeg, here('results/hindcast_output/yoy_anchovy', 
                    'anchovy_SVC_sdmtmb.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Pacific Sanddab ----
# Make mesh object with matrices
yoy_sdab_mesh <- make_mesh(yoy_sdab,
                           xy_cols = c("X", "Y"),
                           cutoff = 20)
plot(yoy_sdab_mesh)

# Select models
sdab_model_small <- sdmTMB_select_small(yoy_sdab, yoy_sdab_mesh) 
sdab_model_small[[2]] # u_vint_100m
sdab_model_large <- sdmTMB_select_large(yoy_sdab, yoy_sdab_mesh) 
sdab_model_large[[2]] # u_vint_50m

sanity(sdab_model_small[[2]]) 
tidy(sdab_model_small[[2]], 
     effect = "ran_pars", 
     conf.int = TRUE)
sanity(sdab_model_large[[2]]) 
tidy(sdab_model_large[[2]],
     effect = "ran_pars", 
     conf.int = TRUE)

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
par(mfrow = c(2),
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
                              cutoff = 18)
                                                                     
plot(yoy_shortbelly_mesh)

# Select models
shortbelly_model_small <- sdmTMB_select_small(yoy_shortbelly, yoy_shortbelly_mesh) 
shortbelly_model_small[[2]] # spice_iso26
shortbelly_model_large <- sdmTMB_select_large(yoy_shortbelly, yoy_shortbelly_mesh) 
shortbelly_model_large[[2]] 

sanity(shortbelly_model_small[[2]]) 
tidy(shortbelly_model_small[[2]], 
     effect = "ran_pars", 
     conf.int = TRUE)
sanity(shortbelly_model_large[[2]]) 
tidy(shortbelly_model_large[[2]],
     effect = "ran_pars", 
     conf.int = TRUE)

# Plot covariates
tiff(here('results/hindcast_output/yoy_shortbelly',
          'shortbelly_partial_dependence_small_sdmtmb.jpg'),
     units = "in",
     width = 38,
     height = 12,
     res = 200)
plot_variables(shortbelly_model_small[[2]], yoy_shortbelly)
dev.off()

tiff(here('results/hindcast_output/yoy_shortbelly',
          'shortbelly_partial_dependence_large_sdmtmb.jpg'),
     units = "in",
     width = 38,
     height = 12,
     res = 200)
plot_variables(shortbelly_model_large[[2]], yoy_shortbelly)
dev.off()

# Predict and plot
nlat = 40
nlon = 60
latd = seq(min(yoy_shortbelly$lat), max(yoy_shortbelly$lat), length.out = nlat)
lond = seq(min(yoy_shortbelly$lon), max(yoy_shortbelly$lon), length.out = nlon)

shortbelly_pred_small <- sdmTMB_grid(yoy_shortbelly, shortbelly_model_small)
shortbelly_pred_large <- sdmTMB_grid(yoy_shortbelly, shortbelly_model_large)

windows(height = 15, width = 20)
par(mfrow = c(2),
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
                            cutoff = 20)
plot(yoy_widow_mesh)

# Select models
widow_model_small <- sdmTMB_select_small(yoy_widow, yoy_widow_mesh) 
widow_model_small[[2]]
widow_model_large <- sdmTMB_select_large(yoy_widow, yoy_widow_mesh) 
widow_model_large[[2]] # v_cu

sanity(widow_model_small[[2]]) 
tidy(widow_model_small[[2]], 
     effect = "ran_pars", 
     conf.int = TRUE)
sanity(widow_model_large[[2]]) # works with cutoff = 20
tidy(widow_model_large[[2]],
     effect = "ran_pars", 
     conf.int = TRUE)

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
par(mfrow = c(2),
    mar = c(6.4, 7.2, 1.6, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
sdmTMB_map(yoy_widow, widow_pred_small)
sdmTMB_map(yoy_widow, widow_pred_large)


# Market Squid ----
# Make mesh object with matrices
yoy_squid_mesh <- make_mesh(yoy_squid,
                            xy_cols = c("X", "Y"),
                            cutoff = 20)
plot(yoy_squid_mesh)

# Select models
squid_model_small <- sdmTMB_select_small(yoy_squid, yoy_squid_mesh) 
squid_model_small[[2]] 
squid_model_large <- sdmTMB_select_large(yoy_squid, yoy_squid_mesh) 
squid_model_large[[2]]

sanity(squid_model_small[[2]]) 
tidy(squid_model_small[[2]], 
     effect = "ran_pars", 
     conf.int = TRUE)
sanity(squid_model_large[[2]]) 
tidy(squid_model_large[[2]],
     effect = "ran_pars", 
     conf.int = TRUE)
