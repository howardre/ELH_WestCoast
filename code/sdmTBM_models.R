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
yoy_anchovy <- filter(read_data('yoy_anch.Rdata'), latitude < 42)
yoy_widow <- read_data('yoy_widw.Rdata') 
yoy_shortbelly <- read_data('yoy_sbly.Rdata')
yoy_sdab <- filter(read_data('yoy_dab.Rdata'), latitude > 36)
yoy_squid <- read_data('yoy_squid.Rdata')

state_labels <- data.frame(name = c("Washington", "Oregon", "California"),
                           lat = c(47, 44.0, 37.0),
                           lon = c(-121.0, -121.0, -120.0))
nlat = 40
nlon = 60
extra_years <- c(2020:2099)
nsim <- 500
bubble_color <- colorRampPalette(c(sequential_hcl(15, palette = "Viridis")))

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

# Pacific Hake ----
# Make mesh object with matrices
yoy_hake_mesh <- make_mesh(yoy_hake, 
                           xy_cols = c("X", "Y"),
                           cutoff = 15)
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
hake_data <- hake_model_small$sdm_vgeo$data
hake_data$small_resid <- residuals(hake_model_small$sdm_vgeo)
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
plot_variables(hake_model_small$sdm_vgeo, hake_data)
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

hake_pred_small <- sdmTMB_grid(yoy_hake, hake_model_small$sdm_vgeo)
hake_pred_small$zeta_s_vgeo[hake_pred_small$dist > 60000] <- NA 
hake_pred_large <- sdmTMB_grid(yoy_hake, hake_model_large$sdm_iso26)
hake_pred_large$zeta_s_depth_iso26[hake_pred_large$dist > 60000] <- NA 

# Niche overlap
# Schoener's D seems to make the most sense (Carroll et al., 2019)
hake_pred_small$large_scaled <- hake_pred_large$preds_scaled
hake_pred_small$p_small <- hake_pred_small$preds_scaled / sum(hake_pred_small$preds_scaled, na.rm = T)
hake_pred_small$p_large <- hake_pred_small$large_scaled / sum(hake_pred_small$large_scaled, na.rm = T)
hake_schoener <- 1 - 0.5 * sum(abs(hake_pred_small$p_small - hake_pred_small$p_large), na.rm = TRUE)

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

# Spatial maps
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
sdmTMB_SVC(yoy_hake, hake_pred_small, "Small (15-35 mm)", "Latitude \u00B0N", 
           hake_pred_small$zeta_s_vgeo, "Geostrophic Current")
sdmTMB_SVC(yoy_hake, hake_pred_large, "Large (36-81 mm)", " ",
           hake_pred_large$zeta_s_depth_iso26, "26 kg/m\u00B3 Isopycnal \n Depth")
dev.copy(jpeg, here('results/hindcast_output/yoy_hake', 
                    'hake_SVC_sdmtmb.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# SVC terms
# Create predictions for original data
yoy_hake_pred_small <- predict(hake_model_small$sdm_vgeo, newdata = yoy_hake)
yoy_hake_pred_large <- predict(hake_model_large$sdm_iso26, newdata = yoy_hake)
colnames(yoy_hake_pred_small)[which(names(yoy_hake_pred_small) == "est")] <- "est_small"
colnames(yoy_hake_pred_large)[which(names(yoy_hake_pred_large) == "est")] <- "est_large"
yoy_hake_preds <- merge(yoy_hake_pred_small,
                        yoy_hake_pred_large[, c("jday", "latitude", "longitude", "year", "est_large")],
                        by = c("jday", "latitude", "longitude", "year"),
                        all.x = TRUE)

# Small
hake_small_tidy <- tidy(hake_model_small$sdm_vgeo, conf.int = TRUE)
hake_small_zeta_sim <- predict(hake_model_small$sdm_vgeo, 
                               nsim = nsim,
                               sims_var = "zeta_s")
hake_small_sims <- spread_sims(hake_model_small$sdm_vgeo,
                               nsim = nsim)
hake_small_beta <- hake_small_sims$vgeo
hake_small_combined <- hake_small_beta + t(hake_small_zeta_sim)

# Large
hake_large_tidy <- tidy(hake_model_large$sdm_iso26, conf.int = TRUE)
hake_large_zeta_sim <- predict(hake_model_large$sdm_iso26, 
                               nsim = nsim,
                               sims_var = "zeta_s")
hake_large_sims <- spread_sims(hake_model_large$sdm_iso26,
                               nsim = nsim)
hake_large_beta <- hake_large_sims$depth_iso26
hake_large_combined <- hake_large_beta + t(hake_large_zeta_sim)

yoy_hake_preds$small_zeta_sim <- apply(t(hake_small_zeta_sim), 2, median)
yoy_hake_preds$small_vgeo_sim <- apply(hake_small_combined, 2, median)
yoy_hake_preds$small_vgeo_lwr <- apply(hake_small_combined, 2, quantile, probs = 0.10)
yoy_hake_preds$small_vgeo_upr <- apply(hake_small_combined, 2, quantile, probs = 0.90)
yoy_hake_preds$large_zeta_sim <- apply(t(hake_large_zeta_sim), 2, median)
yoy_hake_preds$large_iso26_sim <- apply(hake_large_combined, 2, median)
yoy_hake_preds$large_iso26_lwr <- apply(hake_large_combined, 2, quantile, probs = 0.10)
yoy_hake_preds$large_iso26_upr <- apply(hake_large_combined, 2, quantile, probs = 0.90)

small_hake_effect <- hake_small_tidy[hake_small_tidy$term == "vgeo", 2]
large_hake_effect <- hake_large_tidy[hake_large_tidy$term == "depth_iso26", 2]

yoy_hake_mean <- yoy_hake_preds %>%
  group_by(station) %>%
  summarise(year = max(year),
            latitude = mean(latitude),
            longitude = mean(longitude),
            median_est_count_small = median(est_small),
            median_est_count_large = median(est_large),
            small_effect = as.numeric(hake_small_tidy[hake_small_tidy$term == "vgeo", 2]) + 
              hake_small_zeta_sim$zeta_s_vgeo,
            small_vgeo_sim = median(small_vgeo_sim),
            small_vgeo_lwr = median(small_vgeo_lwr),
            small_vgeo_upr = median(small_vgeo_upr),
            large_effect = as.numeric(hake_large_tidy[hake_large_tidy$term == "depth_iso26", 2]) + 
              mean(hake_large_zeta_sim$zeta_s_vgeo),
            large_iso26_sim = median(large_iso26_sim),
            large_iso26_lwr = median(large_iso26_lwr),
            large_iso26_upr = median(large_iso26_upr))

hake_small_cut_pts <- seq(0, max(yoy_hake_mean$small_effect), by = 0.3)
hake_small_cuts <- with(yoy_hake_mean, cut(small_effect + abs(small_effect) + 0.02, 
                                            hake_small_cut_pts))
hake_large_cut_pts <- seq(0, max(yoy_hake_mean$large_effect), by = 0.3)
hake_large_cuts <- with(yoy_hake_mean, cut(large_effect + abs(large_effect) + 0.08, 
                                            hake_large_cut_pts))


windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
image(lond,
      latd,
      t(matrix(hake_pred_small$est,
               nrow = length(latd),
               ncol = length(lond),
               byrow = T)),
      xlim = c(-126, -116),
      ylim = range(yoy_hake_mean$latitude, na.rm = TRUE) + c(-.4, .5),
      axes = FALSE,
      xlab = "",
      ylab = "")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "mintcream")
par(new = TRUE)
symbols(yoy_hake_mean$longitude,
        yoy_hake_mean$latitude,
        circle = yoy_hake_mean$median_est_count_small + 
          abs(min(yoy_hake_mean$median_est_count_small)) + 0.01,
        inches = 0.2,
        add = T,
        fg = bubble_color(10)[hake_small_cuts],
        bg = alpha(bubble_color(10)[hake_small_cuts], 0.6))  
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
     cex = 2.6,
     family = "serif")

image(lond,
      latd,
      t(matrix(hake_pred_large$est,
               nrow = length(latd),
               ncol = length(lond),
               byrow = T)),
      xlim = c(-126, -116),
      ylim = range(yoy_hake_mean$latitude, na.rm = TRUE) + c(-.4, .5),
      axes = FALSE,
      xlab = "",
      ylab = "")
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "mintcream")
par(new = TRUE)
symbols(yoy_hake_mean$longitude,
        yoy_hake_mean$latitude,
        circle = yoy_hake_mean$median_est_count_large + 
          abs(min(yoy_hake_mean$median_est_count_large)) + 0.01,
        inches = 0.2,
        add = T,
        fg = bubble_color(10)[hake_large_cuts],
        bg = alpha(bubble_color(10)[hake_large_cuts], 0.6))  
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
     cex = 2.6,
     family = "serif")

dev.copy(jpeg, here('results/hindcast_output/yoy_hake', 
                    'hake_SVC_sdmtmb.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()
  

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
anchovy_data <- anchovy_model_small$sdm_vgeo$data
anchovy_data$small_resid <- residuals(anchovy_model_small$sdm_vgeo)
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

anchovy_pred_small <- sdmTMB_grid(yoy_anchovy, anchovy_model_small$sdm_vgeo)
anchovy_pred_small$zeta_s_vgeo[anchovy_pred_small$dist > 60000] <- NA 
anchovy_pred_large <- sdmTMB_grid(yoy_anchovy, anchovy_model_large$sdm_uvint100m)
anchovy_pred_large$zeta_s_u_vint_100m[anchovy_pred_large$dist > 60000] <- NA 

# Niche overlap
# Schoener's D seems to make the most sense (Carroll et al., 2019)
anchovy_pred_small$large_scaled <- anchovy_pred_large$preds_scaled
anchovy_pred_small$p_small <- anchovy_pred_small$preds_scaled / sum(anchovy_pred_small$preds_scaled, na.rm = T)
anchovy_pred_small$p_large <- anchovy_pred_small$large_scaled / sum(anchovy_pred_small$large_scaled, na.rm = T)
anchovy_schoener <- 1 - 0.5 * sum(abs(anchovy_pred_small$p_small - anchovy_pred_small$p_large), na.rm = TRUE)

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
           anchovy_pred_small$zeta_s_vgeo, "Geostrophic Current")
sdmTMB_SVC(yoy_anchovy, anchovy_pred_large, "Large (36-85 mm)", " ",
           anchovy_pred_large$zeta_s_u_vint_100m, "Eastward u vertically \n integrated 0-100 m")
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
sdab_data <- sdab_model_small$sdm_v_cu$data
sdab_data$small_resid <- residuals(sdab_model_small$sdm_uvint50m)
sdab_data$large_resid <- residuals(sdab_model_large$sdm_uvint100m)

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
sdab_pred_small$zeta_s_u_vint_50m[sdab_pred_small$dist > 50000] <- NA 
sdab_pred_large <- sdmTMB_grid(yoy_sdab, sdab_model_large$sdm_uvint100m)
sdab_pred_large$zeta_s_u_vint_100m[sdab_pred_large$dist > 50000] <- NA 

# Niche overlap
# Schoener's D seems to make the most sense (Carroll et al., 2019)
sdab_pred_small$large_scaled <- sdab_pred_large$preds_scaled
sdab_pred_small$p_small <- sdab_pred_small$preds_scaled / sum(sdab_pred_small$preds_scaled, na.rm = T)
sdab_pred_small$p_large <- sdab_pred_small$large_scaled / sum(sdab_pred_small$large_scaled, na.rm = T)
sdab_schoener <- 1 - 0.5 * sum(abs(sdab_pred_small$p_small - sdab_pred_small$p_large), na.rm = TRUE)

# Overall predictions
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
sdmTMB_map(yoy_sdab, sdab_pred_small, "Small (16-25 mm)", "Latitude \u00B0N")
sdmTMB_map(yoy_sdab, sdab_pred_large, "Large (26-55 mm)", " ")
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
sdmTMB_SVC(yoy_sdab, sdab_pred_small, "Small (16-25 mm)", "Latitude \u00B0N", 
           sdab_pred_small$zeta_s_u_vint_50m, "Eastward u vertically \n integrated 0-50 m")
sdmTMB_SVC(yoy_sdab, sdab_pred_large, "Large (26-55 mm)", " ",
           sdab_pred_large$zeta_s_u_vint_100m, "Eastward u vertically \n integrated 0-100 m")
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
shortbelly_pred_small$zeta_s_depth_iso26[shortbelly_pred_small$dist > 50000] <- NA 
shortbelly_pred_large <- sdmTMB_grid(yoy_shortbelly, shortbelly_model_large$sdm_vgeo)
shortbelly_pred_large$zeta_s_vgeo[shortbelly_pred_large$dist > 50000] <- NA 

# Niche overlap
# Schoener's D seems to make the most sense (Carroll et al., 2019)
shortbelly_pred_small$large_scaled <- shortbelly_pred_large$preds_scaled
shortbelly_pred_small$p_small <- shortbelly_pred_small$preds_scaled / sum(shortbelly_pred_small$preds_scaled, na.rm = T)
shortbelly_pred_small$p_large <- shortbelly_pred_small$large_scaled / sum(shortbelly_pred_small$large_scaled, na.rm = T)
shortbelly_schoener <- 1 - 0.5 * sum(abs(shortbelly_pred_small$p_small - shortbelly_pred_small$p_large), na.rm = TRUE)

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
           shortbelly_pred_small$zeta_s_depth_iso26, "26 kg/m\u00B3 Isopycnal \n Depth")
sdmTMB_SVC(yoy_shortbelly, shortbelly_pred_large, "Large (36-78 mm)", " ",
           shortbelly_pred_large$zeta_s_vgeo, "Geostrophic Current")
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
widow_data <- widow_model_small$sdm_iso26$data
widow_data$small_resid <- residuals(widow_model_small$sdm_iso26)
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

widow_pred_small <- sdmTMB_grid(yoy_widow, widow_model_small$sdm_iso26)
widow_pred_small$zeta_s_depth_iso26[widow_pred_small$dist > 50000] <- NA 
widow_pred_large <- sdmTMB_grid(yoy_widow, widow_model_large$sdm_spice)
widow_pred_large$zeta_s_spice_iso26[widow_pred_large$dist > 50000] <- NA 

# Niche overlap
# Schoener's D seems to make the most sense (Carroll et al., 2019)
widow_pred_small$large_scaled <- widow_pred_large$preds_scaled
widow_pred_small$p_small <- widow_pred_small$preds_scaled / sum(widow_pred_small$preds_scaled, na.rm = T)
widow_pred_small$p_large <- widow_pred_small$large_scaled / sum(widow_pred_small$large_scaled, na.rm = T)
widow_schoener <- 1 - 0.5 * sum(abs(widow_pred_small$p_small - widow_pred_small$p_large), na.rm = TRUE)

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
           widow_pred_small$zeta_s_depth_iso26, "26 kg/m\u00B3 Isopycnal \n Depth")
sdmTMB_SVC(yoy_widow, widow_pred_large, "Large (33-64 mm)", " ",
           widow_pred_large$zeta_s_spice_iso26, "Spiciness")
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
squid_pred_small$zeta_s_depth_iso26[squid_pred_small$dist > 50000] <- NA 
squid_pred_large <- sdmTMB_grid(yoy_squid, squid_model_large$sdm_spice)
squid_pred_large$zeta_s_spice_iso26[squid_pred_large$dist > 50000] <- NA 

# Niche overlap
# Schoener's D seems to make the most sense (Carroll et al., 2019)
squid_pred_small$large_scaled <- squid_pred_large$preds_scaled
squid_pred_small$p_small <- squid_pred_small$preds_scaled / sum(squid_pred_small$preds_scaled, na.rm = T)
squid_pred_small$p_large <- squid_pred_small$large_scaled / sum(squid_pred_small$large_scaled, na.rm = T)
squid_schoener <- 1 - 0.5 * sum(abs(squid_pred_small$p_small - squid_pred_small$p_large), na.rm = TRUE)

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
           squid_pred_small$zeta_s_depth_iso26, "26 kg/m\u00B3 Isopycnal \n Depth")
sdmTMB_SVC(yoy_squid, squid_pred_large, "Large (71-139 mm)", " ",
           squid_pred_large$zeta_s_spice_iso26, "Spiciness")
dev.copy(jpeg, here('results/hindcast_output/yoy_squid', 
                    'squid_SVC_sdmtmb.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()