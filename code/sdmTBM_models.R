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
yoy_hake_mesh <- make_mesh(yoy_hake, 
                           xy_cols = c("X", "Y"),
                           cutoff = 18)
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
                             depth_iso26 - 1,
                           extra_time = extra_years,
                           spatial_varying = ~ 0 + depth_iso26,
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
                             depth_iso26 - 1,
                           extra_time = extra_years,
                           spatial_varying = ~ 0 + depth_iso26,
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
                             s(sst_scaled, k = 3),
                           extra_time = extra_years,
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
                             depth_iso26 - 1,
                           extra_time = extra_years,
                           spatial_varying = ~ 0 + depth_iso26,
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
                             depth_iso26 - 1,
                           extra_time = extra_years,
                           spatial_varying = ~ 0 + depth_iso26,
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
     conf.level = 0.95) 

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

dev.copy(jpeg, here('results/hindcast_output/yoy_sdab', 
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
tiff(here('results/hindcast_output/yoy_sdab',
          'sdab_partial_dependence_small_sdmtmb.jpg'),
     units = "in",
     width = 30,
     height = 12,
     res = 200)
plot_variables(sdab_model_small, yoy_sdab)
dev.off()

tiff(here('results/hindcast_output/yoy_sdab',
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
yoy_shortbelly_mesh <- make_mesh(yoy_shortbelly, 
                                 xy_cols = c("X", "Y"),
                                 cutoff = 18)
plot(yoy_shortbelly_mesh) 

# Select models
# Calculate deviance explained compared to null model
shortbelly_model_small <- sdmTMB_select(yoy_shortbelly, yoy_shortbelly_mesh, "small") 
shortbelly_small_stat <- calc_stat(shortbelly_model_small, yoy_shortbelly_mesh, yoy_shortbelly, "small")
saveRDS(shortbelly_model_small, here('data', 'shortbelly_models_small'))

shortbelly_model_large <- sdmTMB_select(yoy_shortbelly, yoy_shortbelly_mesh, "large") 
shortbelly_large_stat <- calc_stat(shortbelly_model_large, yoy_shortbelly_mesh, yoy_shortbelly, "large")
saveRDS(shortbelly_model_large, here('data', 'shortbelly_models_large'))

# Load models
shortbelly_model_small <- readRDS(here('data', 'shortbelly_models_small'))
shortbelly_model_large <- readRDS(here('data', 'shortbelly_models_large'))

# Error checks
shortbelly_small_stat
rownames(shortbelly_small_stat)[which.max(shortbelly_small_stat$log_likelihood)]
sanity(shortbelly_model_small$sdm_iso26)
tidy(shortbelly_model_small$sdm_iso26, 
     conf.int = TRUE,
     conf.level = 0.90) # no

shortbelly_large_stat
rownames(shortbelly_large_stat)[which.max(shortbelly_large_stat$log_likelihood)]
sanity(shortbelly_model_large$sdm_vgeo)
tidy(shortbelly_model_large$sdm_vgeo, 
     conf.int = TRUE,
     conf.level = 0.90) # 0.001

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
shortbelly_model_small_cv <- sdmTMB_compare(yoy_shortbelly, yoy_shortbelly_mesh, "small")
shortbelly_model_large_cv <- sdmTMB_compare(yoy_shortbelly, yoy_shortbelly_mesh, "large")

shortbelly_small_best <- shortbelly_model_small_cv[[which.max(sapply(1:length(shortbelly_model_small_cv), 
                                                         function(x) (shortbelly_model_small_cv[[x]]$sum_loglik)))]]
shortbelly_large_best <- shortbelly_model_large_cv[[which.max(sapply(1:length(shortbelly_model_large_cv), 
                                                         function(x) (shortbelly_model_large_cv[[x]]$sum_loglik)))]]

saveRDS(shortbelly_model_small_cv, here('data', 'shortbelly_models_small_cv'))
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


# Widow Rockfish ----
# Make mesh object with matrices
yoy_widow_mesh <- make_mesh(yoy_widow,
                            xy_cols = c("X", "Y"),
                            cutoff = 18)
plot(yoy_widow_mesh) 

# Select models
# Calculate deviance explained compared to null model
widow_model_small <- sdmTMB_select(yoy_widow, yoy_widow_mesh, "small") 
widow_small_stat <- calc_stat(widow_model_small, yoy_widow_mesh, yoy_widow, "small")
saveRDS(widow_model_small, here('data', 'widow_models_small'))

widow_model_large <- sdmTMB_select(yoy_widow, yoy_widow_mesh, "large") 
widow_large_stat <- calc_stat(widow_model_large, yoy_widow_mesh, yoy_widow, "large")
saveRDS(widow_model_large, here('data', 'widow_models_large'))

# Load models
widow_model_small <- readRDS(here('data', 'widow_models_small'))
widow_model_large <- readRDS(here('data', 'widow_models_large'))

# Error checks
widow_small_stat
rownames(widow_small_stat)[which.max(widow_small_stat$log_likelihood)]
sanity(widow_model_small$sdm_vmax_cu)
tidy(widow_model_small$sdm_vmax_cu, 
     conf.int = TRUE,
     conf.level = 0.9) # no

widow_large_stat
rownames(widow_large_stat)[which.max(widow_large_stat$log_likelihood)]
sanity(widow_model_large$sdm_spice)
tidy(widow_model_large$sdm_spice, 
     conf.int = TRUE,
     conf.level = 0.9) # no

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
widow_model_small_cv <- sdmTMB_compare(yoy_widow, yoy_widow_mesh, "small")
widow_model_large_cv <- sdmTMB_compare(yoy_widow, yoy_widow_mesh, "large")

widow_small_best <- widow_model_small_cv[[which.max(sapply(1:length(widow_model_small_cv), 
                                                         function(x) (widow_model_small_cv[[x]]$sum_loglik)))]]
widow_large_best <- widow_model_large_cv[[which.max(sapply(1:length(widow_model_large_cv), 
                                                         function(x) (widow_model_large_cv[[x]]$sum_loglik)))]]

saveRDS(widow_model_small_cv, here('data', 'widow_models_small_cv'))
saveRDS(widow_model_large_cv, here('data', 'widow_models_large_cv'))


# Plot covariates
tiff(here('results/hindcast_output/yoy_widow',
          'widow_partial_dependence_small_sdmtmb.jpg'),
     units = "in",
     width = 38,
     height = 12,
     res = 200)
plot_variables(widow_model_small$sdm_vmax_cu, widow_data)
dev.off()

tiff(here('results/hindcast_output/yoy_widow',
          'widow_partial_dependence_large_sdmtmb.jpg'),
     units = "in",
     width = 38,
     height = 12,
     res = 200)
plot_variables(widow_model_large$sdm_spice, widow_data)
dev.off()


# Market Squid ----
# Make mesh object with matrices
yoy_squid_mesh <- make_mesh(yoy_squid, 
                            xy_cols = c("X", "Y"),
                            cutoff = 18)
plot(yoy_squid_mesh) 

# Select models
# Calculate deviance explained compared to null model
squid_model_small <- sdmTMB_select(yoy_squid, yoy_squid_mesh, "small") 
squid_small_stat <- calc_stat(squid_model_small, yoy_squid_mesh, yoy_squid, "small")
saveRDS(squid_model_small, here('data', 'squid_models_small'))

squid_model_large <- sdmTMB_select(yoy_squid, yoy_squid_mesh, "large") 
squid_large_stat <- calc_stat(squid_model_large, yoy_squid_mesh, yoy_squid, "large")
saveRDS(squid_model_large, here('data', 'squid_models_large'))

# Load models
squid_model_small <- readRDS(here('data', 'squid_models_small'))
squid_model_large <- readRDS(here('data', 'squid_models_large'))

# Error checks
squid_small_stat
rownames(squid_small_stat)[which.max(squid_small_stat$log_likelihood)]
sanity(squid_model_small$sdm_iso26)
tidy(squid_model_small$sdm_iso26, 
     conf.int = TRUE,
     conf.level = 0.999) # 0.001

squid_large_stat
rownames(squid_large_stat)[which.max(squid_large_stat$log_likelihood)]
sanity(squid_model_large$sdm_spice)
tidy(squid_model_large$sdm_spice, 
     conf.int = TRUE,
     conf.level = 0.99) # 0.001

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
squid_model_small_cv <- sdmTMB_compare(yoy_squid, yoy_squid_mesh, "small")
squid_model_large_cv <- sdmTMB_compare(yoy_squid, yoy_squid_mesh, "large")

squid_small_best <- squid_model_small_cv[[which.max(sapply(1:length(squid_model_small_cv), 
                                                         function(x) (squid_model_small_cv[[x]]$sum_loglik)))]]
squid_large_best <- squid_model_large_cv[[which.max(sapply(1:length(squid_model_large_cv), 
                                                         function(x) (squid_model_large_cv[[x]]$sum_loglik)))]]

saveRDS(squid_model_small_cv, here('data', 'squid_models_small_cv'))
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