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
library(scales)

# Load data and functions ----
# Functions
source(here('code/functions', 'vis_gam_COLORS.R'))
source(here('code/functions', 'distance_function.R'))
source(here('code/functions', 'sdmTMB_grid.R'))
source(here('code/functions', 'sdmTMB_map.R'))
source(here('code/functions', 'sdmTMB_select_cv.R'))
source(here('code/functions', 'read_data.R'))
source(here('code/functions', 'plot_variables.R'))
source(here('code/functions', 'sdmTMB_select.R')) # cannot be spatiotemporal & spatial

# Data
# Removing years with just the core area and no PWCC
yoy_hake <- filter(read_data('yoy_hake.Rdata'), year > 2002) 
yoy_anchovy <- filter(read_data('yoy_anch.Rdata'), latitude < 42, year > 2013)
yoy_widow <- filter(read_data('yoy_widw.Rdata'), year > 2000)
yoy_shortbelly <- filter(read_data('yoy_sbly.Rdata'), year > 2000)
yoy_sdab <- filter(read_data('yoy_dab.Rdata'), latitude > 36, year > 2012) # doesn't work otherwise
yoy_squid <- read_data('yoy_squid.Rdata')

state_labels <- data.frame(name = c("Washington", "Oregon", "California"),
                           lat = c(47, 44.0, 37.0),
                           lon = c(-121.0, -121.0, -120.0))
nlat = 40
nlon = 60
extra_years <- c(2020:2100)
nsim <- 500
bubble_color <- colorRampPalette(c(sequential_hcl(15, palette = "Viridis")))


# Pacific Hake ----
# Make mesh object with matrices
yoy_hake_mesh <- make_mesh(yoy_hake, 
                           xy_cols = c("X", "Y"),
                           cutoff = 18)
plot(yoy_hake_mesh) 

# Select models
# Calculate deviance explained compared to null model
hake_model_small <- sdmTMB_select_small(yoy_hake, yoy_hake_mesh) 
hake_small_stat <- calc_stat_small(hake_model_small, yoy_hake_mesh, yoy_hake)
saveRDS(hake_model_small, here('data', 'hake_models_small'))

hake_model_large <- sdmTMB_select_large(yoy_hake, yoy_hake_mesh) 
hake_large_stat <- calc_stat_large(hake_model_large, yoy_hake_mesh, yoy_hake)
saveRDS(hake_model_large, here('data', 'hake_models_large'))

# Load models
hake_model_small <- readRDS(here('data', 'hake_models_small'))
hake_model_large <- readRDS(here('data', 'hake_models_large'))

# Get residuals
hake_data <- hake_model_small$sdm_iso26$data
hake_data$small_resid <- residuals(hake_model_small$sdm_v_cu)
hake_data$large_resid <- residuals(hake_model_large$sdm_iso26)

# Normal QQ plots
windows(height = 8, width = 15)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
qqnorm(hake_data$small_resid, main = "Small Sizes Q-Q Plot")
qqline(hake_data$small_resid)

qqnorm(hake_data$large_resid, main = "Large Sizes Q-Q Plot")
qqline(hake_data$large_resid)

dev.copy(jpeg, here('results/hindcast_output/yoy_hake', 
                    'hake_qq.jpg'), 
         height = 8, 
         width = 15, 
         units = 'in', 
         res = 200)
dev.off()

# Spatial residuals
ggplot(hake_data, 
       aes(X, Y, col = small_resid)) +
  scale_color_gradient2() +
  geom_point() +
  coord_fixed()

ggplot(hake_data, 
       aes(X, Y, col = large_resid)) +
  scale_color_gradient2() +
  geom_point() +
  coord_fixed()

# Use 5-fold cross validation to calculate log likelihood
# Tried LFO, LOYO, 10%, and 10-fold - none worked
hake_train <- filter(yoy_hake, year < 2017)
hake_test <- filter(yoy_hake, year >= 2017)
hake_train_mesh <- make_mesh(hake_train, 
                             xy_cols = c("X", "Y"),
                             cutoff = 18)
test_years <- c(2017:2019)
hake_small_v_cu <- sdmTMB_cv(small ~ 0 + v_cu +
                                s(jday_scaled, k = 3) +
                                s(sst_scaled, k = 3) +
                               s(sss_scaled, k = 3),
                             spatial_varying = ~ 0 + v_cu,
                             data = yoy_hake,
                             mesh = yoy_hake_mesh,
                             family = tweedie(link = "log"),
                             lfo = TRUE,
                             lfo_forecast = 2,
                             lfo_validations = 1,
                             time = "year")

hake_small_v_cu$fold_loglik # fold log-likelihood
hake_small_v_cu$sum_loglik # total log-likelihood


# Get correlation coefficient
hake_test_pred <- predict(hake_small_v_cu,
                          newdata = hake_test,
                          type = "response")
small_sp <- cor.test(test$larvalcatchper10m2, 
                     test$small_pred, 
                     method = 'spearman',
                     exact = FALSE)

plot(yoy_hake_mesh) 


hake_model_small_cv <- sdmTMB_cv_small(yoy_hake, yoy_hake_mesh) 
saveRDS(hake_model_small_cv, here('data', 'hake_models_small_cv'))

hake_model_large_cv <- sdmTMB_cv_large(yoy_hake, yoy_hake_mesh) 
saveRDS(hake_model_large_cv, here('data', 'hake_models_large_cv'))

# Plot covariates
tiff(here('results/hindcast_output/yoy_hake',
          'hake_partial_dependence_small_sdmtmb.jpg'),
     units = "in",
     width = 38,
     height = 12,
     res = 200)
plot_variables(hake_model_small$sdm_iso26, hake_data)
dev.off()

tiff(here('results/hindcast_output/yoy_hake',
          'hake_partial_dependence_large_sdmtmb.jpg'),
     units = "in",
     width = 38,
     height = 12,
     res = 200)
plot_variables(hake_model_large$sdm_iso26, hake_data)
dev.off()

# Predict and plot
latd = seq(min(yoy_hake$lat), max(yoy_hake$lat), length.out = nlat)
lond = seq(min(yoy_hake$lon), max(yoy_hake$lon), length.out = nlon)

hake_pred_small <- sdmTMB_grid(yoy_hake, hake_model_small$sdm_iso26)
hake_pred_small$zeta_s_depth_iso26[hake_pred_small$dist > 60000] <- NA 
hake_pred_large <- sdmTMB_grid(yoy_hake, hake_model_large$sdm_iso26)
hake_pred_large$zeta_s_depth_iso26[hake_pred_large$dist > 60000] <- NA 

# Niche overlap
# Schoener's D seems to make the most sense (Carroll et al., 2019)
hake_pred_small$large_scaled <- hake_pred_large$preds_scaled
hake_pred_small$p_small <- hake_pred_small$preds_scaled / sum(hake_pred_small$preds_scaled, na.rm = T)
hake_pred_small$p_large <- hake_pred_small$large_scaled / sum(hake_pred_small$large_scaled, na.rm = T)
hake_schoener <- 1 - 0.5 * sum(abs(hake_pred_small$p_small - hake_pred_small$p_large), na.rm = TRUE)
hake_schoener

# Overall predictions
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
sdmTMB_map(yoy_hake, hake_pred_small, "Small (15-35 mm)", "Latitude \u00B0N")
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
sdmTMB_SVC(yoy_hake, hake_pred_small, "Small (15-35 mm)", "Latitude \u00B0N", 
           hake_pred_small$zeta_s_depth_iso26, "26 kg/m\u00B3 Isopycnal \n Depth Effect")
sdmTMB_SVC(yoy_hake, hake_pred_large, "Large (36-81 mm)", " ",
           hake_pred_large$zeta_s_depth_iso26, "26 kg/m\u00B3 Isopycnal \n Depth Effect")
dev.copy(jpeg, here('results/hindcast_output/yoy_hake', 
                    'hake_SVC_sdmtmb.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Alternate SVC calculation
zeta_s <- predict(hake_model_small$sdm_vgeo, 
                  newdata = yoy_hake, 
                  nsim = 200, 
                  sims_var = 'zeta_s')
sims <- spread_sims(hake_model_small$sdm_vgeo, nsim = 200)
combined <- sims$vgeo + t(zeta_s)
yoy_hake$vgeo_effect <- apply(combined, 2, median)
yoy_hake$vgeo_effect_lwr <- apply(combined, 2, quantile, probs = 0.1)
yoy_hake$vgeo_effect_upr <- apply(combined, 2, quantile, probs = 0.9)

ggplot(yoy_hake, aes(X, Y)) + 
  geom_point(aes(color = vgeo_effect_upr)) +
  scale_color_viridis()

# Northern Anchovy ----
# Make mesh object with matrices
yoy_anchovy_mesh <- make_mesh(yoy_anchovy, 
                           xy_cols = c("X", "Y"),
                           cutoff = 15)
plot(yoy_anchovy_mesh) 

# Select models
# Calculate deviance explained compared to null model
anchovy_model_small <- sdmTMB_select_small(yoy_anchovy, yoy_anchovy_mesh) 
anchovy_small_stat <- calc_stat_small(anchovy_model_small, yoy_anchovy_mesh, yoy_anchovy)
saveRDS(anchovy_model_small, here('data', 'anchovy_models_small'))

anchovy_model_large <- sdmTMB_select_large(yoy_anchovy, yoy_anchovy_mesh) 
anchovy_large_stat <- calc_stat_large(anchovy_model_large, yoy_anchovy_mesh, yoy_anchovy)
saveRDS(anchovy_model_large, here('data', 'anchovy_models_large'))

# Load models
anchovy_model_small <- readRDS(here('data', 'anchovy_models_small'))
anchovy_model_large <- readRDS(here('data', 'anchovy_models_large'))

# Get residuals
anchovy_data <- anchovy_model_small$sdm_uvint100m$data
anchovy_data$small_resid <- residuals(anchovy_model_small$sdm_uvint100m)
anchovy_data$large_resid <- residuals(anchovy_model_large$sdm_uvint100m)

# Normal QQ plots
windows(height = 8, width = 15)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
qqnorm(anchovy_data$small_resid, main = "Small Sizes Q-Q Plot")
qqline(anchovy_data$small_resid)

qqnorm(anchovy_data$large_resid, main = "Large Sizes Q-Q Plot")
qqline(anchovy_data$large_resid)

dev.copy(jpeg, here('results/hindcast_output/yoy_anchovy', 
                    'anchovy_qq.jpg'), 
         height = 8, 
         width = 15, 
         units = 'in', 
         res = 200)
dev.off()

# Spatial residuals
ggplot(anchovy_data, 
       aes(X, Y, col = small_resid)) +
  scale_color_gradient2() +
  geom_point() +
  coord_fixed()

ggplot(anchovy_data, 
       aes(X, Y, col = large_resid)) +
  scale_color_gradient2() +
  geom_point() +
  coord_fixed()

# Cross validation
anchovy_model_small_cv <- sdmTMB_cv_small(yoy_anchovy, yoy_anchovy_mesh) 
saveRDS(anchovy_model_small_cv, here('data', 'anchovy_models_small_cv'))

anchovy_model_large_cv <- sdmTMB_cv_large(yoy_anchovy, yoy_anchovy_mesh) 
saveRDS(anchovy_model_large_cv, here('data', 'anchovy_models_large_cv'))

# Plot covariates
tiff(here('results/hindcast_output/yoy_anchovy',
          'anchovy_partial_dependence_small_sdmtmb.jpg'),
     units = "in",
     width = 38,
     height = 12,
     res = 200)
plot_variables(anchovy_model_small$sdm_vgeo, anchovy_data)
dev.off()

tiff(here('results/hindcast_output/yoy_anchovy',
          'anchovy_partial_dependence_large_sdmtmb.jpg'),
     units = "in",
     width = 38,
     height = 12,
     res = 200)
plot_variables(anchovy_model_large$sdm_uvint100m, anchovy_data)
dev.off()


# Predict and plot
latd = seq(min(yoy_anchovy$lat), max(yoy_anchovy$lat), length.out = nlat)
lond = seq(min(yoy_anchovy$lon), max(yoy_anchovy$lon), length.out = nlon)

anchovy_pred_small <- sdmTMB_grid(yoy_anchovy, anchovy_model_small$sdm_uvint100m)
anchovy_pred_small$zeta_s_u_vint_100m[anchovy_pred_small$dist > 60000] <- NA 
anchovy_pred_large <- sdmTMB_grid(yoy_anchovy, anchovy_model_large$sdm_uvint100m)
anchovy_pred_large$zeta_s_u_vint_100m[anchovy_pred_large$dist > 60000] <- NA 

# Niche overlap
# Schoener's D seems to make the most sense (Carroll et al., 2019)
anchovy_pred_small$large_scaled <- anchovy_pred_large$preds_scaled
anchovy_pred_small$p_small <- anchovy_pred_small$preds_scaled / sum(anchovy_pred_small$preds_scaled, na.rm = T)
anchovy_pred_small$p_large <- anchovy_pred_small$large_scaled / sum(anchovy_pred_small$large_scaled, na.rm = T)
anchovy_schoener <- 1 - 0.5 * sum(abs(anchovy_pred_small$p_small - anchovy_pred_small$p_large), na.rm = TRUE)
anchovy_schoener

# Overall predictions
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
sdmTMB_map(yoy_anchovy, anchovy_pred_small, "Small (21-35 mm)", "Latitude \u00B0N")
sdmTMB_map(yoy_anchovy, anchovy_pred_large, "Large (36-85 mm)", " ")
dev.copy(jpeg, here('results/hindcast_output/yoy_anchovy', 
                    'anchovy_distributions_sdmtmb.jpg'), 
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
sdmTMB_SVC(yoy_anchovy, anchovy_pred_small, "Small (21-35 mm)", "Latitude \u00B0N", 
           anchovy_pred_small$zeta_s_u_vint_100m, "Eastward u vertically \n integrated 0-100m Effect")
sdmTMB_SVC(yoy_anchovy, anchovy_pred_large, "Large (36-85 mm)", " ",
           anchovy_pred_large$zeta_s_u_vint_100m, "Eastward u vertically \n integrated 0-100m Effect")
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
                           cutoff = 15)
plot(yoy_sdab_mesh) 

# Select models
# Calculate deviance explained compared to null model
sdab_model_small <- sdmTMB_select_small(yoy_sdab, yoy_sdab_mesh) 
sdab_small_stat <- calc_stat_small(sdab_model_small, yoy_sdab_mesh, yoy_sdab)
saveRDS(sdab_model_small, here('data', 'sdab_models_small'))

sdab_model_large <- sdmTMB_select_large(yoy_sdab, yoy_sdab_mesh) 
sdab_large_stat <- calc_stat_large(sdab_model_large, yoy_sdab_mesh, yoy_sdab)
saveRDS(sdab_model_large, here('data', 'sdab_models_large'))

# Load models
sdab_model_small <- readRDS(here('data', 'sdab_models_small'))
sdab_model_large <- readRDS(here('data', 'sdab_models_large'))

# Get residuals
sdab_data <- sdab_model_small$sdm_uvint50m$data
sdab_data$small_resid <- residuals(sdab_model_small$sdm_uvint50m)
sdab_data$large_resid <- residuals(sdab_model_large$sdm_uvint50m)

# Normal QQ plots
windows(height = 8, width = 15)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
qqnorm(sdab_data$small_resid, main = "Small Sizes Q-Q Plot")
qqline(sdab_data$small_resid)

qqnorm(sdab_data$large_resid, main = "Large Sizes Q-Q Plot")
qqline(sdab_data$large_resid)

dev.copy(jpeg, here('results/hindcast_output/yoy_sanddab', 
                    'sdab_qq.jpg'), 
         height = 8, 
         width = 15, 
         units = 'in', 
         res = 200)
dev.off()

# Spatial residuals
ggplot(sdab_data, 
       aes(X, Y, col = small_resid)) +
  scale_color_gradient2() +
  geom_point() +
  coord_fixed()

ggplot(sdab_data, 
       aes(X, Y, col = large_resid)) +
  scale_color_gradient2() +
  geom_point() +
  coord_fixed()

# Cross validation
sdab_model_small_cv <- sdmTMB_cv_small(yoy_sdab, yoy_sdab_mesh) 
saveRDS(sdab_model_small_cv, here('data', 'sdab_models_small_cv'))

sdab_model_large_cv <- sdmTMB_cv_large(yoy_sdab, yoy_sdab_mesh) 
saveRDS(sdab_model_large_cv, here('data', 'sdab_models_large_cv'))

# Plot covariates
tiff(here('results/hindcast_output/yoy_sanddab',
          'sdab_partial_dependence_small_sdmtmb.jpg'),
     units = "in",
     width = 38,
     height = 12,
     res = 200)
plot_variables(sdab_model_small$sdm_uvint50m, sdab_data)
dev.off()

tiff(here('results/hindcast_output/yoy_sanddab',
          'sdab_partial_dependence_large_sdmtmb.jpg'),
     units = "in",
     width = 38,
     height = 12,
     res = 200)
plot_variables(sdab_model_large$sdm_uvint100m, sdab_data)
dev.off()


# Predict and plot
latd = seq(min(yoy_sdab$lat), max(yoy_sdab$lat), length.out = nlat)
lond = seq(min(yoy_sdab$lon), max(yoy_sdab$lon), length.out = nlon)

sdab_pred_small <- sdmTMB_grid(yoy_sdab, sdab_model_small$sdm_uvint50m)
sdab_pred_small$zeta_s_u_vint_50m[sdab_pred_small$dist > 60000] <- NA 
sdab_pred_large <- sdmTMB_grid(yoy_sdab, sdab_model_large$sdm_uvint50m)
sdab_pred_large$zeta_s_u_vint_50m[sdab_pred_large$dist > 60000] <- NA 

# Niche overlap
# Schoener's D seems to make the most sense (Carroll et al., 2019)
sdab_pred_small$large_scaled <- sdab_pred_large$preds_scaled
sdab_pred_small$p_small <- sdab_pred_small$preds_scaled / sum(sdab_pred_small$preds_scaled, na.rm = T)
sdab_pred_small$p_large <- sdab_pred_small$large_scaled / sum(sdab_pred_small$large_scaled, na.rm = T)
sdab_schoener <- 1 - 0.5 * sum(abs(sdab_pred_small$p_small - sdab_pred_small$p_large), na.rm = TRUE)
sdab_schoener

# Overall predictions
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
sdab_map(yoy_sdab, sdab_pred_small, "Small (16-25 mm)", "Latitude \u00B0N")
sdab_map(yoy_sdab, sdab_pred_large, "Large (26-55 mm)", " ")
dev.copy(jpeg, here('results/hindcast_output/yoy_sanddab', 
                    'sdab_distributions_sdmtmb.jpg'), 
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
sdab_SVC(yoy_sdab, sdab_pred_small, "Small (16-25 mm)", "Latitude \u00B0N", 
         sdab_pred_small$zeta_s_u_vint_50m, "Eastward u Vertically \n Integrated 0-100m Effect")
sdab_SVC(yoy_sdab, sdab_pred_large, "Large (26-55 mm)", " ",
         sdab_pred_large$zeta_s_u_vint_50m, "Eastward u Vertically \n Integrated 0-50m Effect")
dev.copy(jpeg, here('results/hindcast_output/yoy_sanddab', 
                    'sdab_SVC_sdmtmb.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()


# Shortbelly Rockfish ----
# Make mesh object with matrices
yoy_shortbelly_mesh <- make_mesh(yoy_shortbelly, 
                                 xy_cols = c("X", "Y"),
                                 cutoff = 15)
plot(yoy_shortbelly_mesh) 

# Select models
# Calculate deviance explained compared to null model
shortbelly_model_small <- sdmTMB_select_small(yoy_shortbelly, yoy_shortbelly_mesh) 
shortbelly_small_stat <- calc_stat_small(shortbelly_model_small, yoy_shortbelly_mesh, yoy_shortbelly)
saveRDS(shortbelly_model_small, here('data', 'shortbelly_models_small'))

shortbelly_model_large <- sdmTMB_select_large(yoy_shortbelly, yoy_shortbelly_mesh) 
shortbelly_large_stat <- calc_stat_large(shortbelly_model_large, yoy_shortbelly_mesh, yoy_shortbelly)
saveRDS(shortbelly_model_large, here('data', 'shortbelly_models_large'))

# Load models
shortbelly_model_small <- readRDS(here('data', 'shortbelly_models_small'))
shortbelly_model_large <- readRDS(here('data', 'shortbelly_models_large'))

# Get residuals
shortbelly_data <- shortbelly_model_small$sdm_iso26$data
shortbelly_data$small_resid <- residuals(shortbelly_model_small$sdm_iso26)
shortbelly_data$large_resid <- residuals(shortbelly_model_large$sdm_vgeo)

# Normal QQ plots
windows(height = 8, width = 15)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
qqnorm(shortbelly_data$small_resid, main = "Small Sizes Q-Q Plot")
qqline(shortbelly_data$small_resid)

qqnorm(shortbelly_data$large_resid, main = "Large Sizes Q-Q Plot")
qqline(shortbelly_data$large_resid)

dev.copy(jpeg, here('results/hindcast_output/yoy_shortbelly', 
                    'shortbelly_qq.jpg'), 
         height = 8, 
         width = 15, 
         units = 'in', 
         res = 200)
dev.off()

# Spatial residuals
ggplot(shortbelly_data, 
       aes(X, Y, col = small_resid)) +
  scale_color_gradient2() +
  geom_point() +
  coord_fixed()

ggplot(shortbelly_data, 
       aes(X, Y, col = large_resid)) +
  scale_color_gradient2() +
  geom_point() +
  coord_fixed()

# Cross validation
shortbelly_model_small_cv <- sdmTMB_cv_small(yoy_shortbelly, yoy_shortbelly_mesh) 
saveRDS(shortbelly_model_small_cv, here('data', 'shortbelly_models_small_cv'))

shortbelly_model_large_cv <- sdmTMB_cv_large(yoy_shortbelly, yoy_shortbelly_mesh) 
saveRDS(shortbelly_model_large_cv, here('data', 'shortbelly_models_large_cv'))

# Plot covariates
tiff(here('results/hindcast_output/yoy_shortbelly',
          'shortbelly_partial_dependence_small_sdmtmb.jpg'),
     units = "in",
     width = 38,
     height = 12,
     res = 200)
plot_variables(shortbelly_model_small$sdm_iso26, shortbelly_data)
dev.off()

tiff(here('results/hindcast_output/yoy_shortbelly',
          'shortbelly_partial_dependence_large_sdmtmb.jpg'),
     units = "in",
     width = 38,
     height = 12,
     res = 200)
plot_variables(shortbelly_model_large$sdm_vgeo, shortbelly_data)
dev.off()


# Predict and plot
latd = seq(min(yoy_shortbelly$lat), max(yoy_shortbelly$lat), length.out = nlat)
lond = seq(min(yoy_shortbelly$lon), max(yoy_shortbelly$lon), length.out = nlon)

shortbelly_pred_small <- sdmTMB_grid(yoy_shortbelly, shortbelly_model_small$sdm_iso26)
shortbelly_pred_small$zeta_s_depth_iso26[shortbelly_pred_small$dist > 60000] <- NA 
shortbelly_pred_large <- sdmTMB_grid(yoy_shortbelly, shortbelly_model_large$sdm_vgeo)
shortbelly_pred_large$zeta_s_vgeo[shortbelly_pred_large$dist > 60000] <- NA 

# Niche overlap
# Schoener's D seems to make the most sense (Carroll et al., 2019)
shortbelly_pred_small$large_scaled <- shortbelly_pred_large$preds_scaled
shortbelly_pred_small$p_small <- shortbelly_pred_small$preds_scaled / sum(shortbelly_pred_small$preds_scaled, na.rm = T)
shortbelly_pred_small$p_large <- shortbelly_pred_small$large_scaled / sum(shortbelly_pred_small$large_scaled, na.rm = T)
shortbelly_schoener <- 1 - 0.5 * sum(abs(shortbelly_pred_small$p_small - shortbelly_pred_small$p_large), na.rm = TRUE)
shortbelly_schoener

# Overall predictions
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
sdmTMB_map(yoy_shortbelly, shortbelly_pred_small, "Small (11-35 mm)", "Latitude \u00B0N")
sdmTMB_map(yoy_shortbelly, shortbelly_pred_large, "Large (36-78 mm)", " ")
dev.copy(jpeg, here('results/hindcast_output/yoy_shortbelly', 
                    'shortbelly_distributions_sdmtmb.jpg'), 
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
sdmTMB_SVC(yoy_shortbelly, shortbelly_pred_small, "Small (11-35 mm)", "Latitude \u00B0N", 
           shortbelly_pred_small$zeta_s_depth_iso26, "26 kg/m\u00B3 Isopycnal \n Depth Effect")
sdmTMB_SVC(yoy_shortbelly, shortbelly_pred_large, "Large (36-78 mm)", " ",
           shortbelly_pred_large$zeta_s_vgeo, "Geostrophic Current \n Effect")
dev.copy(jpeg, here('results/hindcast_output/yoy_shortbelly', 
                    'shortbelly_SVC_sdmtmb.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()


# Widow Rockfish ----
# Make mesh object with matrices
yoy_widow_mesh <- make_mesh(yoy_widow, 
                           xy_cols = c("X", "Y"),
                           cutoff = 15)
plot(yoy_widow_mesh) 

# Select models
# Calculate deviance explained compared to null model
widow_model_small <- sdmTMB_select_small(yoy_widow, yoy_widow_mesh) 
widow_small_stat <- calc_stat_small(widow_model_small, yoy_widow_mesh, yoy_widow)
saveRDS(widow_model_small, here('data', 'widow_models_small'))

widow_model_large <- sdmTMB_select_large(yoy_widow, yoy_widow_mesh) 
widow_large_stat <- calc_stat_large(widow_model_large, yoy_widow_mesh, yoy_widow)
saveRDS(widow_model_large, here('data', 'widow_models_large'))

# Load models
widow_model_small <- readRDS(here('data', 'widow_models_small'))
widow_model_large <- readRDS(here('data', 'widow_models_large'))

# Get residuals
widow_data <- widow_model_small$sdm_vmax_cu$data
widow_data$small_resid <- residuals(widow_model_small$sdm_vmax_cu)
widow_data$large_resid <- residuals(widow_model_large$sdm_spice)

# Normal QQ plots
windows(height = 8, width = 15)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
qqnorm(widow_data$small_resid, main = "Small Sizes Q-Q Plot")
qqline(widow_data$small_resid)

qqnorm(widow_data$large_resid, main = "Large Sizes Q-Q Plot")
qqline(widow_data$large_resid)

dev.copy(jpeg, here('results/hindcast_output/yoy_widow', 
                    'widow_qq.jpg'), 
         height = 8, 
         width = 15, 
         units = 'in', 
         res = 200)
dev.off()

# Spatial residuals
ggplot(widow_data, 
       aes(X, Y, col = small_resid)) +
  scale_color_gradient2() +
  geom_point() +
  coord_fixed()

ggplot(widow_data, 
       aes(X, Y, col = large_resid)) +
  scale_color_gradient2() +
  geom_point() +
  coord_fixed()

# Cross validation
widow_model_small_cv <- sdmTMB_cv_small(yoy_widow, yoy_widow_mesh) 
saveRDS(widow_model_small_cv, here('data', 'widow_models_small_cv'))

widow_model_large_cv <- sdmTMB_cv_large(yoy_widow, yoy_widow_mesh) 
saveRDS(widow_model_large_cv, here('data', 'widow_models_large_cv'))

# Plot covariates
tiff(here('results/hindcast_output/yoy_widow',
          'widow_partial_dependence_small_sdmtmb.jpg'),
     units = "in",
     width = 38,
     height = 12,
     res = 200)
plot_variables(widow_model_small$sdm_iso26, widow_data)
dev.off()

tiff(here('results/hindcast_output/yoy_widow',
          'widow_partial_dependence_large_sdmtmb.jpg'),
     units = "in",
     width = 38,
     height = 12,
     res = 200)
plot_variables(widow_model_large$sdm_spice, widow_data)
dev.off()


# Predict and plot
latd = seq(min(yoy_widow$lat), max(yoy_widow$lat), length.out = nlat)
lond = seq(min(yoy_widow$lon), max(yoy_widow$lon), length.out = nlon)

widow_pred_small <- sdmTMB_grid(yoy_widow, widow_model_small$sdm_vmax_cu)
widow_pred_small$zeta_s_vmax_cu[widow_pred_small$dist > 60000] <- NA 
widow_pred_large <- sdmTMB_grid(yoy_widow, widow_model_large$sdm_spice)
widow_pred_large$zeta_s_spice_iso26[widow_pred_large$dist > 60000] <- NA 

# Niche overlap
# Schoener's D seems to make the most sense (Carroll et al., 2019)
widow_pred_small$large_scaled <- widow_pred_large$preds_scaled
widow_pred_small$p_small <- widow_pred_small$preds_scaled / sum(widow_pred_small$preds_scaled, na.rm = T)
widow_pred_small$p_large <- widow_pred_small$large_scaled / sum(widow_pred_small$large_scaled, na.rm = T)
widow_schoener <- 1 - 0.5 * sum(abs(widow_pred_small$p_small - widow_pred_small$p_large), na.rm = TRUE)
widow_schoener

# Overall predictions
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
sdmTMB_map(yoy_widow, widow_pred_small, "Small (17-32 mm)", "Latitude \u00B0N")
sdmTMB_map(yoy_widow, widow_pred_large, "Large (33-64 mm)", " ")
dev.copy(jpeg, here('results/hindcast_output/yoy_widow', 
                    'widow_distributions_sdmtmb.jpg'), 
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
sdmTMB_SVC(yoy_widow, widow_pred_small, "Small (17-32 mm)", "Latitude \u00B0N", 
           widow_pred_small$zeta_s_vmax_cu, "Max Velocity CA \n Undercurrent Effect ")
sdmTMB_SVC(yoy_widow, widow_pred_large, "Large (33-64 mm)", " ",
           widow_pred_large$zeta_s_spice_iso26, "Spiciness Effect")
dev.copy(jpeg, here('results/hindcast_output/yoy_widow', 
                    'widow_SVC_sdmtmb.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()


# Market Squid ----
# Make mesh object with matrices
yoy_squid_mesh <- make_mesh(yoy_squid, 
                            xy_cols = c("X", "Y"),
                            cutoff = 15)
plot(yoy_squid_mesh) 

# Select models
# Calculate deviance explained compared to null model
squid_model_small <- sdmTMB_select_small(yoy_squid, yoy_squid_mesh) 
squid_small_stat <- calc_stat_small(squid_model_small, yoy_squid_mesh, yoy_squid)
saveRDS(squid_model_small, here('data', 'squid_models_small'))

squid_model_large <- sdmTMB_select_large(yoy_squid, yoy_squid_mesh) 
squid_large_stat <- calc_stat_large(squid_model_large, yoy_squid_mesh, yoy_squid)
saveRDS(squid_model_large, here('data', 'squid_models_large'))

# Load models
squid_model_small <- readRDS(here('data', 'squid_models_small'))
squid_model_large <- readRDS(here('data', 'squid_models_large'))

# Get residuals
squid_data <- squid_model_small$sdm_iso26$data
squid_data$small_resid <- residuals(squid_model_small$sdm_iso26)
squid_data$large_resid <- residuals(squid_model_large$sdm_spice)

# Normal QQ plots
windows(height = 8, width = 15)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
qqnorm(squid_data$small_resid, main = "Small Sizes Q-Q Plot")
qqline(squid_data$small_resid)

qqnorm(squid_data$large_resid, main = "Large Sizes Q-Q Plot")
qqline(squid_data$large_resid)

dev.copy(jpeg, here('results/hindcast_output/yoy_squid', 
                    'squid_qq.jpg'), 
         height = 8, 
         width = 15, 
         units = 'in', 
         res = 200)
dev.off()

# Spatial residuals
ggplot(squid_data, 
       aes(X, Y, col = small_resid)) +
  scale_color_gradient2() +
  geom_point() +
  coord_fixed()

ggplot(squid_data, 
       aes(X, Y, col = large_resid)) +
  scale_color_gradient2() +
  geom_point() +
  coord_fixed()

# Cross validation
squid_model_small_cv <- sdmTMB_cv_small(yoy_squid, yoy_squid_mesh) 
saveRDS(squid_model_small_cv, here('data', 'squid_models_small_cv'))

squid_model_large_cv <- sdmTMB_cv_large(yoy_squid, yoy_squid_mesh) 
saveRDS(squid_model_large_cv, here('data', 'squid_models_large_cv'))

# Plot covariates
tiff(here('results/hindcast_output/yoy_squid',
          'squid_partial_dependence_small_sdmtmb.jpg'),
     units = "in",
     width = 38,
     height = 12,
     res = 200)
plot_variables(squid_model_small$sdm_iso26, squid_data)
dev.off()

tiff(here('results/hindcast_output/yoy_squid',
          'squid_partial_dependence_large_sdmtmb.jpg'),
     units = "in",
     width = 38,
     height = 12,
     res = 200)
plot_variables(squid_model_large$sdm_spice, squid_data)
dev.off()


# Predict and plot
latd = seq(min(yoy_squid$lat), max(yoy_squid$lat), length.out = nlat)
lond = seq(min(yoy_squid$lon), max(yoy_squid$lon), length.out = nlon)

squid_pred_small <- sdmTMB_grid(yoy_squid, squid_model_small$sdm_iso26)
squid_pred_small$zeta_s_depth_iso26[squid_pred_small$dist > 60000] <- NA 
squid_pred_large <- sdmTMB_grid(yoy_squid, squid_model_large$sdm_spice)
squid_pred_large$zeta_s_spice_iso26[squid_pred_large$dist > 60000] <- NA 

# Niche overlap
# Schoener's D seems to make the most sense (Carroll et al., 2019)
squid_pred_small$large_scaled <- squid_pred_large$preds_scaled
squid_pred_small$p_small <- squid_pred_small$preds_scaled / sum(squid_pred_small$preds_scaled, na.rm = T)
squid_pred_small$p_large <- squid_pred_small$large_scaled / sum(squid_pred_small$large_scaled, na.rm = T)
squid_schoener <- 1 - 0.5 * sum(abs(squid_pred_small$p_small - squid_pred_small$p_large), na.rm = TRUE)
squid_schoener

# Overall predictions
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
sdmTMB_map(yoy_squid, squid_pred_small, "Small (11-70 mm)", "Latitude \u00B0N")
sdmTMB_map(yoy_squid, squid_pred_large, "Large (71-139 mm)", " ")
dev.copy(jpeg, here('results/hindcast_output/yoy_squid', 
                    'squid_distributions_sdmtmb.jpg'), 
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
sdmTMB_SVC(yoy_squid, squid_pred_small, "Small (11-70 mm)", "Latitude \u00B0N", 
           squid_pred_small$zeta_s_depth_iso26, "26 kg/m\u00B3 Isopycnal \n Depth Effect")
sdmTMB_SVC(yoy_squid, squid_pred_large, "Large (71-139 mm)", " ",
           squid_pred_large$zeta_s_spice_iso26, "Spiciness Effect")
dev.copy(jpeg, here('results/hindcast_output/yoy_squid', 
                    'squid_SVC_sdmtmb.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()


# # Sandbox ----
# the_mesh <- make_mesh(yoy_hake,
#                       xy_cols = c("X", "Y"),
#                       cutoff = 15)
# sdm_test <- sdmTMB(small ~  0 +
#                      s(jday_scaled, k = 3) +
#                      s(sst_scaled, k = 3) +
#                      s(sss_scaled, k = 3) +
#                      v_cu,
#                    spatial_varying = ~ 0 + v_cu,
#                    extra_time = extra_years,
#                    data = yoy_hake,
#                    mesh = the_mesh,
#                    spatial = "on",
#                    time = "year",
#                    family = tweedie(link = "log"),
#                    spatiotemporal = "off",
#                    control = sdmTMBcontrol(newton_loops = 1,
#                                            nlminb_loops = 2))
# sanity(sdm_test)
# tidy(sdm_test,
#      conf.int = TRUE)
# 
# 
# yoy_hake <- mutate(yoy_hake, fold = ifelse(year < 2013, 1, 2))
# 
# library(caret)
# clust <- createFolds(factor(yoy_hake$small), k = 3, list = FALSE)
# 
# cv_test <- try(sdmTMB_cv(small ~ 0 + vgeo +
#                        s(jday_scaled, k = 3) +
#                        s(sst_scaled, k = 3) +
#                        s(sss_scaled, k = 3),
#                      spatial_varying = ~ 0 + vgeo,
#                      data = yoy_hake,
#                      mesh = the_mesh,
#                      spatial = "off",
#                      spatiotemporal = "off",
#                      family = tweedie(link = "log"),
#                      k_folds = length(unique(clust)),
#                      fold_ids = clust,
#                      parallel = TRUE))