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

extra_years <- c(2020:2099)

withJavaLogging = function(expr, silentSuccess=FALSE, stopIsFatal=TRUE) {
  hasFailed = FALSE
  messages = list()
  warnings = list()
  logger = function(obj) {
    # Change behaviour based on type of message
    level = sapply(class(obj), switch, debug="DEBUG", message="INFO", warning="WARN", caughtError = "ERROR",
                   error=if (stopIsFatal) "FATAL" else "ERROR", "")
    level = c(level[level != ""], "ERROR")[1]
    simpleMessage = switch(level, DEBUG=,INFO=TRUE, FALSE)
    quashable = switch(level, DEBUG=,INFO=,WARN=TRUE, FALSE)
    
    # Format message
    time  = format(Sys.time(), "%Y-%m-%d %H:%M:%OS3")
    txt   = conditionMessage(obj)
    if (!simpleMessage) txt = paste(txt, "\n", sep="")
    msg = paste(time, level, txt, sep=" ")
    calls = sys.calls()
    calls = calls[1:length(calls)-1]
    trace = limitedLabels(c(calls, attr(obj, "calls")))
    if (!simpleMessage && length(trace) > 0) {
      trace = trace[length(trace):1]
      msg = paste(msg, "  ", paste("at", trace, collapse="\n  "), "\n", sep="")
    }
    
    # Output message
    if (silentSuccess && !hasFailed && quashable) {
      messages <<- append(messages, msg)
      if (level == "WARN") warnings <<- append(warnings, msg)
    } else {
      if (silentSuccess && !hasFailed) {
        cat(paste(messages, collapse=""))
        hasFailed <<- TRUE
      }
      cat(msg)
    }
    
    # Muffle any redundant output of the same message
    optionalRestart = function(r) { res = findRestart(r); if (!is.null(res)) invokeRestart(res) }
    optionalRestart("muffleMessage")
    optionalRestart("muffleWarning")
  }
  vexpr = withCallingHandlers(withVisible(expr),
                              debug=logger, message=logger, warning=logger, caughtError=logger, error=logger)
  if (silentSuccess && !hasFailed) {
    cat(paste(warnings, collapse=""))
  }
  if (vexpr$visible) vexpr$value else invisible(vexpr$value)
}

# Sandbox ----
the_mesh <- make_mesh(yoy_hake,
                      xy_cols = c("X", "Y"),
                      cutoff = 15)
sdm_test <- sdmTMB(small ~  0 +
                     s(jday_scaled, k = 3) +
                     s(sst_scaled, k = 3) +
                     s(sss_scaled, k = 3) +
                     v_cu,
                   spatial_varying = ~ 0 + v_cu,
                   extra_time = extra_years,
                   data = yoy_hake,
                   mesh = the_mesh,
                   spatial = "on",
                   time = "year",
                   family = tweedie(link = "log"),
                   spatiotemporal = "off",
                   control = sdmTMBcontrol(newton_loops = 1,
                                           nlminb_loops = 2))
sanity(sdm_test)
tidy(sdm_test,
     conf.int = TRUE)


yoy_hake <- mutate(yoy_hake, fold = ifelse(year < 2013, 1, 2))


cv_test <- try(sdmTMB_cv(small ~ 0 + v_cu +
                            s(jday_scaled, k = 3) +
                            s(sst_scaled, k = 3) +
                            s(sss_scaled, k = 3),
                          spatial_varying = ~ 0 + v_cu,
                          data = yoy_hake,
                          mesh = the_mesh,
                          spatiotemporal = "off",
                          family = tweedie(link = "log"),
                          k_folds = 2,
                          fold_ids = yoy_hake$fold,
                          parallel = TRUE))
  log_lik[[i]] <- sum(models[[i]]$data$cv_loglik[which(models[[i]]$data$year == test_years[i])])


cv_test <- try(sdmTMB_cv(small ~ 0 + vmax_cu +
                       s(jday_scaled, k = 3) +
                       s(sst_scaled, k = 3) +
                       s(sss_scaled, k = 3),
                     spatial_varying = ~ 0 + vmax_cu,
                     data = yoy_hake,
                     mesh = the_mesh,
                     spatiotemporal = "off",
                     family = tweedie(link = "log"),
                     k_folds = 5,
                     parallel = TRUE))

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

# Get residuals
hake_data <- hake_model_small$sdm_vgeo$data
hake_data$small_resid <- residuals(hake_model_small$sdm_vgeo)
hake_data$large_resid <- residuals(hake_model_large$sdm_iso26)

# Normal QQ plots
qqnorm(hake_data$small_resid)
qqline(hake_data$small_resid)

qqnorm(hake_data$large_resid)
qqline(hake_data$large_resid)

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
hake_pred_small$zeta_s_u_vint_100m[hake_pred_small$dist > 50000] <- NA 
hake_pred_large <- sdmTMB_grid(yoy_hake, hake_model_large[[2]])
hake_pred_large$zeta_s_spice_iso26[hake_pred_large$dist > 50000] <- NA 


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
sdmTMB_SVC(yoy_hake, hake_pred_small, "Small (15-35 mm)", "Latitude \u00B0N", 
           hake_pred_small$zeta_s_u_vint_100m, "u (0-100m)")
sdmTMB_SVC(yoy_hake, hake_pred_large, "Large (36-81 mm)", " ",
           hake_pred_large$zeta_s_spice_iso26, "Spice")
dev.copy(jpeg, here('results/hindcast_output/yoy_hake', 
                    'hake_SVC_sdmtmb.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Cross validation
plan(multisession)
clust <- as.numeric(as.factor(unique(yoy_hake$year))) # year clusters, may not work with spatiotemporal "on"
hake_small_cv <- sdmTMB_cv(small ~ 0 + spice_iso26 +
                             jday_scaled +
                             sst_scaled +
                             sss_scaled,
                           spatial_varying = ~ 0 + spice_iso26,
                           data = yoy_hake,
                           mesh = yoy_hake_mesh,
                           spatial = "off",
                           time = "year",
                           family = tweedie(link = "log"),
                           spatiotemporal = "iid",
                           control = sdmTMBcontrol(newton_loops = 1,
                                                   nlminb_loops = 2),
                           fold_ids = clust,
                           k_folds = length(unique(clust)))

plot(hake_small_cv$fold_loglik) # higher values better
hake_small_cv$sum_loglik

clust <- kmeans(yoy_hake[, c("X", "Y")], k = 50)$cluster # spatial clustering
hake_large_cv <- sdmTMB_cv(large ~ 0 + spice_iso26 +
                             jday_scaled +
                             sst_scaled +
                             sss_scaled,
                           spatial_varying = ~ 0 + spice_iso26,
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

# Get residuals
anchovy_data <- anchovy_model_small$sdm_vgeo$data
anchovy_data$small_resid <- residuals(anchovy_model_small$sdm_iso26)
anchovy_data$large_resid <- residuals(anchovy_model_large$sdm_vmax_cu)

# Normal QQ plots
qqnorm(anchovy_data$small_resid)
qqline(anchovy_data$small_resid)

qqnorm(anchovy_data$large_resid)
qqline(anchovy_data$large_resid)

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

anchovy_model_small_cv <- sdmTMB_cv_small(yoy_anchovy, yoy_anchovy_mesh) 
saveRDS(anchovy_model_small_cv, here('data', 'anchovy_models_small_cv'))

anchovy_model_large_cv <- sdmTMB_cv_large(yoy_anchovy, yoy_anchovy_mesh) 
saveRDS(anchovy_model_large_cv, here('data', 'anchovy_models_large_cv'))

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

# Get residuals
sdab_data <- sdab_model_small$sdm_vgeo$data
sdab_data$small_resid <- residuals(sdab_model_small$sdm_v_cu)
sdab_data$large_resid <- residuals(sdab_model_large$sdm_spice)

# Normal QQ plots
qqnorm(sdab_data$small_resid)
qqline(sdab_data$small_resid)

qqnorm(sdab_data$large_resid)
qqline(sdab_data$large_resid)

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

sdab_model_small_cv <- sdmTMB_cv_small(yoy_sdab, yoy_sdab_mesh) 
saveRDS(sdab_model_small_cv, here('data', 'sdab_models_small_cv'))

sdab_model_large_cv <- sdmTMB_cv_large(yoy_sdab, yoy_sdab_mesh) 
saveRDS(sdab_model_large_cv, here('data', 'sdab_models_large_cv'))

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

# Get residuals
shortbelly_data <- shortbelly_model_small$sdm_vgeo$data
shortbelly_data$small_resid <- residuals(shortbelly_model_small$sdm_iso26)
shortbelly_data$large_resid <- residuals(shortbelly_model_large$sdm_vgeo)

# Normal QQ plots
qqnorm(shortbelly_data$small_resid)
qqline(shortbelly_data$small_resid)

qqnorm(shortbelly_data$large_resid)
qqline(shortbelly_data$large_resid)

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

# Get residuals
widow_data <- widow_model_small$sdm_vgeo$data
widow_data$small_resid <- residuals(widow_model_small$sdm_iso26)
widow_data$large_resid <- residuals(widow_model_large$sdm_spice)

# Normal QQ plots
qqnorm(widow_data$small_resid)
qqline(widow_data$small_resid)

qqnorm(widow_data$large_resid)
qqline(widow_data$large_resid)

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

widow_model_small_cv <- sdmTMB_cv_small(yoy_widow, yoy_widow_mesh) 
saveRDS(widow_model_small_cv, here('data', 'widow_models_small_cv'))

widow_model_large_cv <- sdmTMB_cv_large(yoy_widow, yoy_widow_mesh) 
saveRDS(widow_model_large_cv, here('data', 'widow_models_large_cv'))

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

# Get residuals
squid_data <- squid_model_small$sdm_vgeo$data
squid_data$small_resid <- residuals(squid_model_small$sdm_vgeo)
squid_data$large_resid <- residuals(squid_model_large$sdm_iso26)

# Normal QQ plots
qqnorm(squid_data$small_resid)
qqline(squid_data$small_resid)

qqnorm(squid_data$large_resid)
qqline(squid_data$large_resid)

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
