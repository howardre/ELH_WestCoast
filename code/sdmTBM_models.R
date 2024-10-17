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
library(fmesher)

# Load data and functions ----
# Functions
source(here('code/functions', 'vis_gam_COLORS.R'))
source(here('code/functions', 'distance_function.R'))
source(here('code/functions', 'sdmTMB_select_cv.R'))
source(here('code/functions', 'read_data.R'))
source(here('code/functions', 'plot_variables.R'))
source(here('code/functions', 'sdmTMB_select.R')) # cannot be spatiotemporal & spatial
source(here('code/functions', 'sdmTMB_map.R'))
source(here('code/functions', 'sdmTMB_grid.R'))

options(future.globals.maxSize = 8000 * 1024^2) # required for CV

# New mesh function
get_mesh <- function(data) {
  # Coastline mesh
  # Calculate max edge. The spatial range is about 1/3 of the study area;
  # the max.edge is about 1/5 of that. 
  spatial_range = diff(range(data$Y)) / 3
  max_edge = spatial_range / 5
  
  mesh <- make_mesh(data, 
                    xy_cols = c("X", "Y"), 
                    fmesher_func = fmesher::fm_mesh_2d_inla,
                    max_edge = c(1, 2) * max_edge, # inner and outer max triangle lengths
                    offset = c(1, 2) * max_edge,  # inner and outer border widths
                    cutoff = max_edge / 5 # shortest allowed distance between points
  )
  
  return(mesh)
}

# Data
# Removing years with just the core area and no PWCC
yoy_hake <- filter(read_data('yoy_hake.Rdata'), year > 2002) 
yoy_anchovy <- filter(read_data('yoy_anch.Rdata'), latitude < 42, year > 2013)
yoy_widow <- filter(read_data('yoy_widw.Rdata'), year > 2000)
yoy_shortbelly <- filter(read_data('yoy_sbly.Rdata'), year > 2000)
yoy_sdab <- filter(read_data('yoy_dab.Rdata'), year > 2012) # doesn't work otherwise
yoy_squid <- read_data('yoy_squid.Rdata')

nep_avgs <- readRDS(here('data', 'nep_avgs.Rdata'))
nep_avgs$large_grp <- findInterval(nep_avgs$latitude,
                                   c(32, 36, 40, 44, 48))
nep_avgs$small_grp <- findInterval(nep_avgs$latitude,
                                   c(32, 34, 36, 38, 40, 42, 44, 48))
nep_large <- nep_avgs %>% 
  group_by(years, large_grp) %>%
  summarize(across(c(vgeo, v_cu, vmax_cu), mean))
nep_small <- nep_avgs %>%
  group_by(years, small_grp) %>%
  summarize(across(c(u_vint_50m, u_vint_100m, depth_iso26, spice_iso26), mean))

state_labels <- data.frame(name = c("Washington", "Oregon", "California"),
                           lat = c(47, 44.0, 37.0),
                           lon = c(-121.0, -121.0, -120.0))
nlat = 40
nlon = 60
extra_years <- c(2020:2100)
nsim <- 500

# Pacific Hake ----
# Make mesh object with matrices
yoy_hake_mesh <- get_mesh(yoy_hake)
plot(yoy_hake_mesh) 

# Select models
# Calculate deviance explained compared to null model
hake_model_small_select <- sdmTMB_select(yoy_hake, yoy_hake_mesh, "small") 
hake_small_stat <- calc_stat(hake_model_small_select, yoy_hake_mesh, yoy_hake, "small")
saveRDS(hake_model_small_select, here('data', 'hake_models_small'))

hake_model_large_select <- sdmTMB_select(yoy_hake, yoy_hake_mesh, "large") 
hake_large_stat <- calc_stat(hake_model_large_select, yoy_hake_mesh, yoy_hake, "large")
saveRDS(hake_model_large_select, here('data', 'hake_models_large'))

# Cross validation
hake_model_small_cv <- sdmTMB_compare(yoy_hake, yoy_hake_mesh, "small")
hake_model_large_cv <- sdmTMB_compare(yoy_hake, yoy_hake_mesh, "large")

hake_small_best <- hake_model_small_cv[[which.max(sapply(1:length(hake_model_small_cv), 
                                                         function(x) (hake_model_small_cv[[x]]$sum_loglik)))]]
hake_large_best <- hake_model_large_cv[[which.max(sapply(1:length(hake_model_large_cv), 
                                                         function(x) (hake_model_large_cv[[x]]$sum_loglik)))]]

hake_small_best$models[[1]]$formula
hake_large_best$models[[1]]$formula

saveRDS(hake_model_small_cv, here('data', 'hake_models_small_cv'))
saveRDS(hake_model_large_cv, here('data', 'hake_models_large_cv'))

# Load models
hake_model_small <- sdmTMB(small ~ s(jday_scaled, k = 3) +
                             s(sst_scaled, k = 3) +
                             u_vint_50m - 1,
                           extra_time = extra_years,
                           spatial_varying = ~ 0 + u_vint_50m,
                           data = yoy_hake,
                           mesh = yoy_hake_mesh,
                           spatial = "on",
                           time = "year",
                           family = tweedie(link = "log"),
                           spatiotemporal = "off",
                           control = sdmTMBcontrol(newton_loops = 1,
                                                   nlminb_loops = 2))
hake_model_large <- sdmTMB(large ~ s(jday_scaled, k = 3) +
                             s(sst_scaled, k = 3) +
                             vmax_cu - 1,
                           extra_time = extra_years,
                           spatial_varying = ~ 0 + vmax_cu,
                           data = yoy_hake,
                           mesh = yoy_hake_mesh,
                           spatial = "on",
                           time = "year",
                           family = tweedie(link = "log"),
                           spatiotemporal = "off",
                           control = sdmTMBcontrol(newton_loops = 1,
                                                   nlminb_loops = 2))

# Error checks
hake_small_stat
rownames(hake_small_stat)[which.max(hake_small_stat$log_likelihood)]
sanity(hake_model_small)
tidy(hake_model_small, 
     conf.int = TRUE,
     conf.level = 0.95) 

hake_large_stat
rownames(hake_large_stat)[which.max(hake_large_stat$log_likelihood)]
sanity(hake_model_large)
tidy(hake_model_large, 
     conf.int = TRUE,
     conf.level = 0.99) # 0.01

# Get residuals
yoy_hake$small_resid <- residuals(hake_model_small,
                                   type = "mle-mvn")
yoy_hake$large_resid <- residuals(hake_model_large,
                                   type = "mle-mvn")

# Normal QQ plots
windows(height = 8, width = 15)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
qqnorm(yoy_hake$small_resid, main = "Small Sizes Q-Q Plot")
qqline(yoy_hake$small_resid)

qqnorm(yoy_hake$large_resid, main = "Large Sizes Q-Q Plot")
qqline(yoy_hake$large_resid)

dev.copy(jpeg, here('results/hindcast_output/yoy_hake', 
                    'hake_qq.jpg'), 
         height = 8, 
         width = 15, 
         units = 'in', 
         res = 200)
dev.off()

# Spatial residuals
ggplot(yoy_hake, 
       aes(X, Y, col = small_resid)) +
  scale_color_gradient2() +
  geom_point() +
  coord_fixed()

ggplot(yoy_hake, 
       aes(X, Y, col = large_resid)) +
  scale_color_gradient2() +
  geom_point() +
  coord_fixed()

# Plot covariates
tiff(here('results/hindcast_output/yoy_hake',
          'hake_partial_dependence_small_sdmtmb.jpg'),
     units = "in",
     width = 30,
     height = 12,
     res = 200)
plot_variables(hake_model_small, yoy_hake)
dev.off()

tiff(here('results/hindcast_output/yoy_hake',
          'hake_partial_dependence_large_sdmtmb.jpg'),
     units = "in",
     width = 30,
     height = 12,
     res = 200)
plot_variables(hake_model_large, yoy_hake)
dev.off()

# Get correlation coefficient
hake_small_pred <- predict(hake_model_small,
                           newdata = yoy_hake,
                           type = "response")
small_sp <- cor.test(hake_small_pred$small, 
                     hake_small_pred$est, 
                     method = 'spearman',
                     exact = FALSE)
small_sp


# Alternate SVC calculation
# zeta_s <- predict(hake_model_small$sdm_vgeo, 
#                   newdata = yoy_hake, 
#                   nsim = 200, 
#                   sims_var = 'zeta_s')
# sims <- spread_sims(hake_model_small$sdm_vgeo, nsim = 200)
# combined <- sims$vgeo + t(zeta_s)
# yoy_hake$vgeo_effect <- apply(combined, 2, median)
# yoy_hake$vgeo_effect_lwr <- apply(combined, 2, quantile, probs = 0.1)
# yoy_hake$vgeo_effect_upr <- apply(combined, 2, quantile, probs = 0.9)
# 
# ggplot(yoy_hake, aes(X, Y)) + 
#   geom_point(aes(color = vgeo_effect_upr)) +
#   scale_color_viridis()

# Northern Anchovy ----
# Make mesh object with matrices
yoy_anchovy_mesh <- get_mesh(yoy_anchovy)
plot(yoy_anchovy_mesh) 

# Select models
# Calculate deviance explained compared to null model
anchovy_model_small_select <- sdmTMB_select(yoy_anchovy, yoy_anchovy_mesh, "small") 
anchovy_small_stat <- calc_stat(anchovy_model_small_select, yoy_anchovy_mesh, yoy_anchovy, "small")
saveRDS(anchovy_model_small_select, here('data', 'anchovy_models_small'))

anchovy_model_large_select <- sdmTMB_select(yoy_anchovy, yoy_anchovy_mesh, "large") 
anchovy_large_stat <- calc_stat(anchovy_model_large_select, yoy_anchovy_mesh, yoy_anchovy, "large")
saveRDS(anchovy_model_large_select, here('data', 'anchovy_models_large'))

# Cross validation
anchovy_model_small_cv <- sdmTMB_compare(yoy_anchovy, yoy_anchovy_mesh, "small")
anchovy_model_large_cv <- sdmTMB_compare(yoy_anchovy, yoy_anchovy_mesh, "large")

anchovy_small_best <- anchovy_model_small_cv[[which.max(sapply(1:length(anchovy_model_small_cv), 
                                                         function(x) (anchovy_model_small_cv[[x]]$sum_loglik)))]]
anchovy_large_best <- anchovy_model_large_cv[[which.max(sapply(1:length(anchovy_model_large_cv), 
                                                         function(x) (anchovy_model_large_cv[[x]]$sum_loglik)))]]

anchovy_small_best$models[[1]]$formula
anchovy_large_best$models[[1]]$formula

saveRDS(anchovy_model_small_cv, here('data', 'anchovy_models_small_cv'))
saveRDS(anchovy_model_large_cv, here('data', 'anchovy_models_large_cv'))

# Load models
anchovy_model_small <- sdmTMB(small ~ s(jday_scaled, k = 3) +
                                s(sst_scaled, k = 3) +
                                u_vint_100m - 1,
                           extra_time = extra_years,
                           spatial_varying = ~ 0 + u_vint_100m,
                           data = yoy_anchovy,
                           mesh = yoy_anchovy_mesh,
                           spatial = "on",
                           time = "year",
                           family = tweedie(link = "log"),
                           spatiotemporal = "off",
                           control = sdmTMBcontrol(newton_loops = 1,
                                                   nlminb_loops = 2))
anchovy_model_large <- sdmTMB(large ~ s(jday_scaled, k = 3) +
                             s(sst_scaled, k = 3) +
                             vmax_cu - 1,
                           extra_time = extra_years,
                           spatial_varying = ~ 0 + vmax_cu,
                           data = yoy_anchovy,
                           mesh = yoy_anchovy_mesh,
                           spatial = "on",
                           time = "year",
                           family = tweedie(link = "log"),
                           spatiotemporal = "off",
                           control = sdmTMBcontrol(newton_loops = 1,
                                                   nlminb_loops = 2))

# Error checks
anchovy_small_stat
rownames(anchovy_small_stat)[which.max(anchovy_small_stat$log_likelihood)]
sanity(anchovy_model_small)
tidy(anchovy_model_small, 
     conf.int = TRUE,
     conf.level = 0.95) # 0.05

anchovy_large_stat
rownames(anchovy_large_stat)[which.max(anchovy_large_stat$log_likelihood)]
sanity(anchovy_model_large)
tidy(anchovy_model_large, 
     conf.int = TRUE,
     conf.level = 0.95) 

# Get residuals
yoy_anchovy$small_resid <- residuals(anchovy_model_small,
                                     type = "mle-mvn")
yoy_anchovy$large_resid <- residuals(anchovy_model_large,
                                     type = "mle-mvn")

# Normal QQ plots
windows(height = 8, width = 15)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
qqnorm(yoy_anchovy$small_resid, main = "Small Sizes Q-Q Plot")
qqline(yoy_anchovy$small_resid)

qqnorm(yoy_anchovy$large_resid, main = "Large Sizes Q-Q Plot")
qqline(yoy_anchovy$large_resid)

dev.copy(jpeg, here('results/hindcast_output/yoy_anchovy', 
                    'anchovy_qq.jpg'), 
         height = 8, 
         width = 15, 
         units = 'in', 
         res = 200)
dev.off()

# Spatial residuals
ggplot(yoy_anchovy, 
       aes(X, Y, col = small_resid)) +
  scale_color_gradient2() +
  geom_point() +
  coord_fixed()

ggplot(yoy_anchovy, 
       aes(X, Y, col = large_resid)) +
  scale_color_gradient2() +
  geom_point() +
  coord_fixed()

# Plot covariates
tiff(here('results/hindcast_output/yoy_anchovy',
          'anchovy_partial_dependence_small_sdmtmb.jpg'),
     units = "in",
     width = 30,
     height = 12,
     res = 200)
plot_variables(anchovy_model_small, yoy_anchovy)
dev.off()

tiff(here('results/hindcast_output/yoy_anchovy',
          'anchovy_partial_dependence_large_sdmtmb.jpg'),
     units = "in",
     width = 30,
     height = 12,
     res = 200)
plot_variables(anchovy_model_large, yoy_anchovy)
dev.off()

# Get correlation coefficient
anchovy_small_pred <- predict(anchovy_model_small,
                           newdata = yoy_anchovy,
                           type = "response")
small_sp <- cor.test(anchovy_small_pred$small, 
                     anchovy_small_pred$est, 
                     method = 'spearman',
                     exact = FALSE)
small_sp



# Pacific Sanddab ----
# Make mesh object with matrices
yoy_sdab_mesh <- get_mesh(yoy_sdab)
plot(yoy_sdab_mesh) 

# Select models
# Calculate deviance explained compared to null model
sdab_model_small_select <- sdmTMB_select(yoy_sdab, yoy_sdab_mesh, "small") 
sdab_small_stat <- calc_stat(sdab_model_small_select, yoy_sdab_mesh, yoy_sdab, "small")
saveRDS(sdab_model_small_select, here('data', 'sdab_models_small'))

sdab_model_large_select <- sdmTMB_select(yoy_sdab, yoy_sdab_mesh, "large") 
sdab_large_stat <- calc_stat(sdab_model_large_select, yoy_sdab_mesh, yoy_sdab, "large")
saveRDS(sdab_model_large_select, here('data', 'sdab_models_large'))

# Cross validation
sdab_model_small_cv <- sdmTMB_compare(yoy_sdab, yoy_sdab_mesh, "small")
sdab_model_large_cv <- sdmTMB_compare(yoy_sdab, yoy_sdab_mesh, "large")

sdab_small_best <- sdab_model_small_cv[[which.max(sapply(1:length(sdab_model_small_cv), 
                                                         function(x) (sdab_model_small_cv[[x]]$sum_loglik)))]]
sdab_large_best <- sdab_model_large_cv[[which.max(sapply(1:length(sdab_model_large_cv), 
                                                         function(x) (sdab_model_large_cv[[x]]$sum_loglik)))]]

sdab_small_best$models[[1]]$formula
sdab_large_best$models[[1]]$formula

saveRDS(sdab_model_small_cv, here('data', 'sdab_models_small_cv'))
saveRDS(sdab_model_large_cv, here('data', 'sdab_models_large_cv'))

# Load models
sdab_model_small <- sdmTMB(small ~ s(jday_scaled, k = 3) +
                             s(sst_scaled, k = 3) +
                             u_vint_100m - 1,
                           extra_time = extra_years,
                           spatial_varying = ~ 0 + u_vint_100m,
                           data = yoy_sdab,
                           mesh = yoy_sdab_mesh,
                           spatial = "on",
                           time = "year",
                           family = tweedie(link = "log"),
                           spatiotemporal = "off",
                           control = sdmTMBcontrol(newton_loops = 1,
                                                   nlminb_loops = 2))
sdab_model_large <- sdmTMB(large ~ s(jday_scaled, k = 3) +
                             s(sst_scaled, k = 3) +
                             spice_iso26 - 1,
                           extra_time = extra_years,
                           spatial_varying = ~ 0 + spice_iso26,
                           data = yoy_sdab,
                           mesh = yoy_sdab_mesh,
                           spatial = "on",
                           time = "year",
                           family = tweedie(link = "log"),
                           spatiotemporal = "off",
                           control = sdmTMBcontrol(newton_loops = 1,
                                                   nlminb_loops = 2))

# Error checks
sdab_small_stat
rownames(sdab_small_stat)[which.max(sdab_small_stat$log_likelihood)]
sanity(sdab_model_small)
tidy(sdab_model_small, 
     conf.int = TRUE,
     conf.level = 0.99) # 0.01

sdab_large_stat
rownames(sdab_large_stat)[which.max(sdab_large_stat$log_likelihood)]
sanity(sdab_model_large)
tidy(sdab_model_large, 
     conf.int = TRUE,
     conf.level = 0.99) # 0.01

# Get residuals
yoy_sdab$small_resid <- residuals(sdab_model_small,
                                  type = "mle-mvn")
yoy_sdab$large_resid <- residuals(sdab_model_large,
                                  type = "mle-mvn")

# Normal QQ plots
windows(height = 8, width = 15)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
qqnorm(yoy_sdab$small_resid, main = "Small Sizes Q-Q Plot")
qqline(yoy_sdab$small_resid)

qqnorm(yoy_sdab$large_resid, main = "Large Sizes Q-Q Plot")
qqline(yoy_sdab$large_resid)

dev.copy(jpeg, here('results/hindcast_output/yoy_sanddab', 
                    'sdab_qq.jpg'), 
         height = 8, 
         width = 15, 
         units = 'in', 
         res = 200)
dev.off()

# Spatial residuals
ggplot(yoy_sdab, 
       aes(X, Y, col = small_resid)) +
  scale_color_gradient2() +
  geom_point() +
  coord_fixed()

ggplot(yoy_sdab, 
       aes(X, Y, col = large_resid)) +
  scale_color_gradient2() +
  geom_point() +
  coord_fixed()

# Plot covariates
tiff(here('results/hindcast_output/yoy_sanddab',
          'sdab_partial_dependence_small_sdmtmb.jpg'),
     units = "in",
     width = 30,
     height = 12,
     res = 200)
plot_variables(sdab_model_small, yoy_sdab)
dev.off()

tiff(here('results/hindcast_output/yoy_sanddab',
          'sdab_partial_dependence_large_sdmtmb.jpg'),
     units = "in",
     width = 30,
     height = 12,
     res = 200)
plot_variables(sdab_model_large, yoy_sdab)
dev.off()

# Get correlation coefficient
sdab_small_pred <- predict(sdab_model_small,
                           newdata = yoy_sdab,
                           type = "response")
small_sp <- cor.test(sdab_small_pred$small, 
                     sdab_small_pred$est, 
                     method = 'spearman',
                     exact = FALSE)
small_sp


# Shortbelly Rockfish ----
# Make mesh object with matrices
yoy_shortbelly_mesh <- get_mesh(yoy_shortbelly)
plot(yoy_shortbelly_mesh) 

# Select models
# Calculate deviance explained compared to null model
shortbelly_model_small_select <- sdmTMB_select(yoy_shortbelly, yoy_shortbelly_mesh, "small") 
shortbelly_small_stat <- calc_stat(shortbelly_model_small_select, yoy_shortbelly_mesh, yoy_shortbelly, "small")
saveRDS(shortbelly_model_small_select, here('data', 'shortbelly_models_small'))

shortbelly_model_large_select <- sdmTMB_select(yoy_shortbelly, yoy_shortbelly_mesh, "large") 
shortbelly_large_stat <- calc_stat(shortbelly_model_large_select, yoy_shortbelly_mesh, yoy_shortbelly, "large")
saveRDS(shortbelly_model_large_select, here('data', 'shortbelly_models_large'))

# Cross validation
shortbelly_model_small_cv <- sdmTMB_compare(yoy_shortbelly, yoy_shortbelly_mesh, "small")
shortbelly_model_large_cv <- sdmTMB_compare(yoy_shortbelly, yoy_shortbelly_mesh, "large")

shortbelly_small_best <- shortbelly_model_small_cv[[which.max(sapply(1:length(shortbelly_model_small_cv), 
                                                         function(x) (shortbelly_model_small_cv[[x]]$sum_loglik)))]]
shortbelly_large_best <- shortbelly_model_large_cv[[which.max(sapply(1:length(shortbelly_model_large_cv), 
                                                         function(x) (shortbelly_model_large_cv[[x]]$sum_loglik)))]]

shortbelly_small_best$models[[1]]$formula
shortbelly_large_best$models[[1]]$formula

saveRDS(shortbelly_model_small_cv, here('data', 'shortbelly_models_small_cv'))
saveRDS(shortbelly_model_large_cv, here('data', 'shortbelly_models_large_cv'))

# Load models
shortbelly_model_small <- sdmTMB(small ~ s(jday_scaled, k = 3) +
                             s(sst_scaled, k = 3) +
                             depth_iso26 - 1,
                           extra_time = extra_years,
                           spatial_varying = ~ 0 + depth_iso26,
                           data = yoy_shortbelly,
                           mesh = yoy_shortbelly_mesh,
                           spatial = "on",
                           time = "year",
                           family = tweedie(link = "log"),
                           spatiotemporal = "off",
                           control = sdmTMBcontrol(newton_loops = 1,
                                                   nlminb_loops = 2))
shortbelly_model_large <- sdmTMB(large ~ s(jday_scaled, k = 3) +
                             s(sst_scaled, k = 3) +
                             spice_iso26 - 1,
                           extra_time = extra_years,
                           spatial_varying = ~ 0 + spice_iso26,
                           data = yoy_shortbelly,
                           mesh = yoy_shortbelly_mesh,
                           spatial = "on",
                           time = "year",
                           family = tweedie(link = "log"),
                           spatiotemporal = "off",
                           control = sdmTMBcontrol(newton_loops = 1,
                                                   nlminb_loops = 2))

# Error checks
shortbelly_small_stat
rownames(shortbelly_small_stat)[which.max(shortbelly_small_stat$log_likelihood)]
sanity(shortbelly_model_small)
tidy(shortbelly_model_small, 
     conf.int = TRUE,
     conf.level = 0.95) 

shortbelly_large_stat
rownames(shortbelly_large_stat)[which.max(shortbelly_large_stat$log_likelihood)]
sanity(shortbelly_model_large)
tidy(shortbelly_model_large, 
     conf.int = TRUE,
     conf.level = 0.99) # 0.01

# Get residuals
yoy_shortbelly$small_resid <- residuals(shortbelly_model_small,
                                  type = "mle-mvn")
yoy_shortbelly$large_resid <- residuals(shortbelly_model_large,
                                  type = "mle-mvn")

# Normal QQ plots
windows(height = 8, width = 15)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
qqnorm(yoy_shortbelly$small_resid, main = "Small Sizes Q-Q Plot")
qqline(yoy_shortbelly$small_resid)

qqnorm(yoy_shortbelly$large_resid, main = "Large Sizes Q-Q Plot")
qqline(yoy_shortbelly$large_resid)

dev.copy(jpeg, here('results/hindcast_output/yoy_shortbelly', 
                    'shortbelly_qq.jpg'), 
         height = 8, 
         width = 15, 
         units = 'in', 
         res = 200)
dev.off()

# Spatial residuals
ggplot(yoy_shortbelly, 
       aes(X, Y, col = small_resid)) +
  scale_color_gradient2() +
  geom_point() +
  coord_fixed()

ggplot(yoy_shortbelly, 
       aes(X, Y, col = large_resid)) +
  scale_color_gradient2() +
  geom_point() +
  coord_fixed()

# Plot covariates
tiff(here('results/hindcast_output/yoy_shortbelly',
          'shortbelly_partial_dependence_small_sdmtmb.jpg'),
     units = "in",
     width = 30,
     height = 12,
     res = 200)
plot_variables(shortbelly_model_small, yoy_shortbelly)
dev.off()

tiff(here('results/hindcast_output/yoy_shortbelly',
          'shortbelly_partial_dependence_large_sdmtmb.jpg'),
     units = "in",
     width = 30,
     height = 12,
     res = 200)
plot_variables(shortbelly_model_large, yoy_shortbelly)
dev.off()

# Get correlation coefficient
shortbelly_small_pred <- predict(shortbelly_model_small,
                           newdata = yoy_shortbelly,
                           type = "response")
small_sp <- cor.test(shortbelly_small_pred$small, 
                     shortbelly_small_pred$est, 
                     method = 'spearman',
                     exact = FALSE)
small_sp


# Widow Rockfish ----
# Make mesh object with matrices
yoy_widow_mesh <- get_mesh(yoy_widow)
plot(yoy_widow_mesh) 

# Select models
# Calculate deviance explained compared to null model
widow_model_small_select <- sdmTMB_select(yoy_widow, yoy_widow_mesh, "small") 
widow_small_stat <- calc_stat(widow_model_small_select, yoy_widow_mesh, yoy_widow, "small")
saveRDS(widow_model_small_select, here('data', 'widow_models_small'))

widow_model_large_select <- sdmTMB_select(yoy_widow, yoy_widow_mesh, "large") 
widow_large_stat <- calc_stat(widow_model_large_select, yoy_widow_mesh, yoy_widow, "large")
saveRDS(widow_model_large_select, here('data', 'widow_models_large'))

# Cross validation
widow_model_small_cv <- sdmTMB_compare(yoy_widow, yoy_widow_mesh, "small")
widow_model_large_cv <- sdmTMB_compare(yoy_widow, yoy_widow_mesh, "large")

widow_small_best <- widow_model_small_cv[[which.max(sapply(1:length(widow_model_small_cv), 
                                                         function(x) (widow_model_small_cv[[x]]$sum_loglik)))]]
widow_large_best <- widow_model_large_cv[[which.max(sapply(1:length(widow_model_large_cv), 
                                                         function(x) (widow_model_large_cv[[x]]$sum_loglik)))]]

widow_small_best$models[[1]]$formula
widow_large_best$models[[1]]$formula

saveRDS(widow_model_small_cv, here('data', 'widow_models_small_cv'))
saveRDS(widow_model_large_cv, here('data', 'widow_models_large_cv'))

# Load models
widow_model_small <- sdmTMB(small ~ s(jday_scaled, k = 3) +
                             s(sst_scaled, k = 3) +
                             v_cu - 1,
                           extra_time = extra_years,
                           spatial_varying = ~ 0 + v_cu,
                           data = yoy_widow,
                           mesh = yoy_widow_mesh,
                           spatial = "on",
                           time = "year",
                           family = tweedie(link = "log"),
                           spatiotemporal = "off",
                           control = sdmTMBcontrol(newton_loops = 1,
                                                   nlminb_loops = 2))
widow_model_large <- sdmTMB(large ~ s(jday_scaled, k = 3) +
                             s(sst_scaled, k = 3) +
                             v_cu - 1,
                           extra_time = extra_years,
                           spatial_varying = ~ 0 + v_cu,
                           data = yoy_widow,
                           mesh = yoy_widow_mesh,
                           spatial = "on",
                           time = "year",
                           family = tweedie(link = "log"),
                           spatiotemporal = "off",
                           control = sdmTMBcontrol(newton_loops = 1,
                                                   nlminb_loops = 2))

# Error checks
widow_small_stat
rownames(widow_small_stat)[which.max(widow_small_stat$log_likelihood)]
sanity(widow_model_small)
tidy(widow_model_small, 
     conf.int = TRUE,
     conf.level = 0.95) # 0.01

widow_large_stat
rownames(widow_large_stat)[which.max(widow_large_stat$log_likelihood)]
sanity(widow_model_large)
tidy(widow_model_large, 
     conf.int = TRUE,
     conf.level = 0.95) # 0.01

# Get residuals
yoy_widow$small_resid <- residuals(widow_model_small,
                                  type = "mle-mvn")
yoy_widow$large_resid <- residuals(widow_model_large,
                                  type = "mle-mvn")

# Normal QQ plots
windows(height = 8, width = 15)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
qqnorm(yoy_widow$small_resid, main = "Small Sizes Q-Q Plot")
qqline(yoy_widow$small_resid)

qqnorm(yoy_widow$large_resid, main = "Large Sizes Q-Q Plot")
qqline(yoy_widow$large_resid)

dev.copy(jpeg, here('results/hindcast_output/yoy_widow', 
                    'widow_qq.jpg'), 
         height = 8, 
         width = 15, 
         units = 'in', 
         res = 200)
dev.off()

# Spatial residuals
ggplot(yoy_widow, 
       aes(X, Y, col = small_resid)) +
  scale_color_gradient2() +
  geom_point() +
  coord_fixed()

ggplot(yoy_widow, 
       aes(X, Y, col = large_resid)) +
  scale_color_gradient2() +
  geom_point() +
  coord_fixed()

# Plot covariates
tiff(here('results/hindcast_output/yoy_widow',
          'widow_partial_dependence_small_sdmtmb.jpg'),
     units = "in",
     width = 30,
     height = 12,
     res = 200)
plot_variables(widow_model_small, yoy_widow)
dev.off()

tiff(here('results/hindcast_output/yoy_widow',
          'widow_partial_dependence_large_sdmtmb.jpg'),
     units = "in",
     width = 30,
     height = 12,
     res = 200)
plot_variables(widow_model_large, yoy_widow)
dev.off()

# Get correlation coefficient
widow_small_pred <- predict(widow_model_small,
                           newdata = yoy_widow,
                           type = "response")
small_sp <- cor.test(widow_small_pred$small, 
                     widow_small_pred$est, 
                     method = 'spearman',
                     exact = FALSE)
small_sp


# Market Squid ----
# Make mesh object with matrices
yoy_squid_mesh <- get_mesh(yoy_squid)
plot(yoy_squid_mesh) 

# Select models
# Calculate deviance explained compared to null model
squid_model_small_select <- sdmTMB_select(yoy_squid, yoy_squid_mesh, "small") 
squid_small_stat <- calc_stat(squid_model_small_select, yoy_squid_mesh, yoy_squid, "small")
saveRDS(squid_model_small_select, here('data', 'squid_models_small'))

squid_model_large_select <- sdmTMB_select(yoy_squid, yoy_squid_mesh, "large") 
squid_large_stat <- calc_stat(squid_model_large_select, yoy_squid_mesh, yoy_squid, "large")
saveRDS(squid_model_large_select, here('data', 'squid_models_large'))

# Cross validation
squid_model_small_cv <- sdmTMB_compare(yoy_squid, yoy_squid_mesh, "small")
squid_model_large_cv <- sdmTMB_compare(yoy_squid, yoy_squid_mesh, "large")

squid_small_best <- squid_model_small_cv[[which.max(sapply(1:length(squid_model_small_cv), 
                                                         function(x) (squid_model_small_cv[[x]]$sum_loglik)))]]
squid_large_best <- squid_model_large_cv[[which.max(sapply(1:length(squid_model_large_cv), 
                                                         function(x) (squid_model_large_cv[[x]]$sum_loglik)))]]

squid_small_best$models[[1]]$formula
squid_large_best$models[[1]]$formula

saveRDS(squid_model_small_cv, here('data', 'squid_models_small_cv'))
saveRDS(squid_model_large_cv, here('data', 'squid_models_large_cv'))

# Load models
squid_model_small <- sdmTMB(small ~ s(jday_scaled, k = 3) +
                             s(sst_scaled, k = 3) +
                             u_vint_100m - 1,
                           extra_time = extra_years,
                           spatial_varying = ~ 0 + u_vint_100m,
                           data = yoy_squid,
                           mesh = yoy_squid_mesh,
                           spatial = "on",
                           time = "year",
                           family = tweedie(link = "log"),
                           spatiotemporal = "off",
                           control = sdmTMBcontrol(newton_loops = 1,
                                                   nlminb_loops = 2))
squid_model_large <- sdmTMB(large ~ s(jday_scaled, k = 3) +
                             s(sst_scaled, k = 3) - 1,
                           extra_time = extra_years,
                           data = yoy_squid,
                           mesh = yoy_squid_mesh,
                           spatial = "on",
                           time = "year",
                           family = tweedie(link = "log"),
                           spatiotemporal = "off",
                           control = sdmTMBcontrol(newton_loops = 1,
                                                   nlminb_loops = 2))

# Error checks
squid_small_stat
rownames(squid_small_stat)[which.max(squid_small_stat$log_likelihood)]
sanity(squid_model_small)
tidy(squid_model_small, 
     conf.int = TRUE,
     conf.level = 0.95) 

squid_large_stat
rownames(squid_large_stat)[which.max(squid_large_stat$log_likelihood)]
sanity(squid_model_large)
tidy(squid_model_large, 
     conf.int = TRUE,
     conf.level = 0.99) 

# Get residuals
yoy_squid$small_resid <- residuals(squid_model_small,
                                  type = "mle-mvn")
yoy_squid$large_resid <- residuals(squid_model_large,
                                  type = "mle-mvn")

# Normal QQ plots
windows(height = 8, width = 15)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
qqnorm(yoy_squid$small_resid, main = "Small Sizes Q-Q Plot")
qqline(yoy_squid$small_resid)

qqnorm(yoy_squid$large_resid, main = "Large Sizes Q-Q Plot")
qqline(yoy_squid$large_resid)

dev.copy(jpeg, here('results/hindcast_output/yoy_squid', 
                    'squid_qq.jpg'), 
         height = 8, 
         width = 15, 
         units = 'in', 
         res = 200)
dev.off()

# Spatial residuals
ggplot(yoy_squid, 
       aes(X, Y, col = small_resid)) +
  scale_color_gradient2() +
  geom_point() +
  coord_fixed()

ggplot(yoy_squid, 
       aes(X, Y, col = large_resid)) +
  scale_color_gradient2() +
  geom_point() +
  coord_fixed()

# Plot covariates
tiff(here('results/hindcast_output/yoy_squid',
          'squid_partial_dependence_small_sdmtmb.jpg'),
     units = "in",
     width = 30,
     height = 12,
     res = 200)
plot_variables(squid_model_small, yoy_squid)
dev.off()

tiff(here('results/hindcast_output/yoy_squid',
          'squid_partial_dependence_large_sdmtmb.jpg'),
     units = "in",
     width = 30,
     height = 12,
     res = 200)
plot_variables(squid_model_large, yoy_squid)
dev.off()

# Get correlation coefficient
squid_small_pred <- predict(squid_model_small,
                           newdata = yoy_squid,
                           type = "response")
small_sp <- cor.test(squid_small_pred$small, 
                     squid_small_pred$est, 
                     method = 'spearman',
                     exact = FALSE)
small_sp

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