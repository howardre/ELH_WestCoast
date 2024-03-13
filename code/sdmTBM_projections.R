### Title: sdmTMB Projecting
### Author: Rebecca Howard
### Date: 03/12/2024

# Load libraries ----
library(tibble)
library(here)
library(ggplot2)
library(dplyr)
library(sdmTMB) 
library(visreg)
library(raster)
library(sf)
library(colorspace)
library(RColorBrewer)
library(mapdata)
library(fields)
library(ggpubr)
library(future)
library(scales)
library(tidync)

# Load data and functions ----
# Functions
source(here('code/functions', 'vis_gam_COLORS.R'))
source(here('code/functions', 'read_data.R'))
source(here('code/functions', 'distance_function.R'))
source(here('code/functions', 'projection_functions.R'))
source(here('code/functions', 'map_project.R'))

# Data
roms_ss <- readRDS(here('data', 'nep_ipsl.Rdata'))
roms_means <- readRDS(here('data', 'nep_ipsl_means.Rdata'))

state_labels <- data.frame(name = c("Washington", "Oregon", "California"),
                           lat = c(47, 44.0, 37.0),
                           lon = c(-121.0, -121.0, -120.0))

extra_years <- c(2020:2100)

### Pacific Hake ---------------------------------------------------------------------------------------------------------------------------------
yoy_hake <- read_data('yoy_hake.Rdata') 
hake_mesh <- make_mesh(data, 
                       xy_cols = c("X", "Y"),
                       cutoff = 15)

set.seed(1993)
hake_small <- sdmTMB(small ~ 0 + vgeo +
                          s(jday_scaled, k = 3) +
                          s(sst_scaled, k = 3) +
                          s(sss_scaled, k = 3),
                        spatial_varying = ~ 0 + vgeo,
                        extra_time = extra_years,
                        data = yoy_hake,
                        mesh = hake_mesh,
                        spatial = "on",
                        time = "year",
                        family = tweedie(link = "log"),
                        spatiotemporal = "off",
                        control = sdmTMBcontrol(newton_loops = 1,
                                                nlminb_loops = 2))
set.seed(1993)
hake_large <- sdmTMB(large ~ 0 + depth_iso26 +
                          s(jday_scaled, k = 3) +
                          s(sst_scaled, k = 3) +
                          s(sss_scaled, k = 3),
                        spatial_varying = ~ 0 + depth_iso26,
                        extra_time = extra_years,
                        data = yoy_hake,
                        mesh = hake_mesh,
                        spatial = "on",
                        time = "year",
                        family = tweedie(link = "log"),
                        spatiotemporal = "off",
                        control = sdmTMBcontrol(newton_loops = 1,
                                                nlminb_loops = 2))

#### IPSL
##### 2020-2040------------------------------------------------------------------------------------------------------------------------------------
hake_ipsl1 <- sdm_cells(yoy_hake, hake_small, hake_large,
                        roms_means, roms_ss, 2020:2040)
saveRDS(hake_ipsl1, file = here("data", "hake_ipsl1.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(hake_ipsl1[[1]], "Small (15-35 mm)", "Latitude \u00B0N")
map_project(hake_ipsl1[[2]], "Large (36-81 mm)", "")
mtext("Pacific Hake 2020-2040", 
      side = 2, 
      line = 50, 
      outer = TRUE, 
      cex = 7,
      at = 0.89)
dev.copy(jpeg, here('results/forecast_output/yoy_hake', 
                    'hake_ipsl1_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

##### 2050-2070------------------------------------------------------------------------------------------------------------------------------------
hake_ipsl2 <- sdm_cells(yoy_hake, hake_small, hake_large,
                        roms_means, roms_ss, 2050:2070)
saveRDS(hake_ipsl2, file = here("data", "hake_ipsl2.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(hake_ipsl2[[1]], "Small (15-35 mm)", "Latitude \u00B0N")
map_project(hake_ipsl2[[2]], "Large (36-81 mm)", "")
dev.copy(jpeg, here('results/forecast_output/yoy_hake', 
                    'hake_ipsl2_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

##### 2080-2100------------------------------------------------------------------------------------------------------------------------------------
hake_ipsl3 <- sdm_cells(yoy_hake, hake_small, hake_large,
                        roms_means, roms_ss, 2080:2100)
saveRDS(hake_ipsl3, file = here("data", "hake_ipsl3.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(hake_ipsl3[[1]], "Small (15-35 mm)", "Latitude \u00B0N")
map_project(hake_ipsl3[[2]], "Large (36-81 mm)", "")
dev.copy(jpeg, here('results/forecast_output/yoy_hake', 
                    'hake_ipsl3_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()


# Northern Anchovy ----
yoy_anchovy <- filter(read_data('yoy_anch.Rdata'), latitude < 42)
anchovy_mesh <- make_mesh(data,
                          xy_cols = c("X", "Y"),
                          cutoff = 15)

set.seed(1993)
anchovy_small <- sdmTMB(small ~ 0 + vgeo +
                          s(jday_scaled, k = 3) +
                          s(sst_scaled, k = 3) +
                          s(sss_scaled, k = 3),
                        spatial_varying = ~ 0 + vgeo,
                        extra_time = extra_years,
                        data = yoy_anchovy,
                        mesh = anchovy_mesh,
                        spatial = "on",
                        time = "year",
                        family = tweedie(link = "log"),
                        spatiotemporal = "off",
                        control = sdmTMBcontrol(newton_loops = 1,
                                                nlminb_loops = 2))
set.seed(1993)
anchovy_large <- sdmTMB(large ~ 0 + u_vint_100m +
                          s(jday_scaled, k = 3) +
                          s(sst_scaled, k = 3) +
                          s(sss_scaled, k = 3),
                        spatial_varying = ~ 0 + u_vint_100m,
                        extra_time = extra_years,
                        data = yoy_anchovy,
                        mesh = anchovy_mesh,
                        spatial = "on",
                        time = "year",
                        family = tweedie(link = "log"),
                        spatiotemporal = "off",
                        control = sdmTMBcontrol(newton_loops = 1,
                                                nlminb_loops = 2))

# Pacific Sanddab ----
yoy_sdab <- filter(read_data('yoy_dab.Rdata'), latitude > 36)
sdab_mesh <- make_mesh(data,
                       xy_cols = c("X", "Y"),
                       cutoff = 15)

set.seed(1993)
sdab_small <- sdmTMB(small ~ 0 + u_vint_50m +
                          s(jday_scaled, k = 3) +
                          s(sst_scaled, k = 3) +
                          s(sss_scaled, k = 3),
                        spatial_varying = ~ 0 + u_vint_50m,
                        extra_time = extra_years,
                        data = yoy_sdab,
                        mesh = sdab_mesh,
                        spatial = "on",
                        time = "year",
                        family = tweedie(link = "log"),
                        spatiotemporal = "off",
                        control = sdmTMBcontrol(newton_loops = 1,
                                                nlminb_loops = 2))
set.seed(1993)
sdab_large <- sdmTMB(large ~ 0 + u_vint_100m +
                          s(jday_scaled, k = 3) +
                          s(sst_scaled, k = 3) +
                          s(sss_scaled, k = 3),
                        spatial_varying = ~ 0 + u_vint_100m,
                        extra_time = extra_years,
                        data = yoy_sdab,
                        mesh = sdab_mesh,
                        spatial = "on",
                        time = "year",
                        family = tweedie(link = "log"),
                        spatiotemporal = "off",
                        control = sdmTMBcontrol(newton_loops = 1,
                                                nlminb_loops = 2))

# Widow Rockfish ----
yoy_widow <- read_data('yoy_widw.Rdata') 
widow_mesh <- make_mesh(data,
                          xy_cols = c("X", "Y"),
                          cutoff = 15)

set.seed(1993)
widow_small <- sdmTMB(small ~ 0 + depth_iso26 +
                          s(jday_scaled, k = 3) +
                          s(sst_scaled, k = 3) +
                          s(sss_scaled, k = 3),
                        spatial_varying = ~ 0 + depth_iso26,
                        extra_time = extra_years,
                        data = yoy_widow,
                        mesh = widow_mesh,
                        spatial = "on",
                        time = "year",
                        family = tweedie(link = "log"),
                        spatiotemporal = "off",
                        control = sdmTMBcontrol(newton_loops = 1,
                                                nlminb_loops = 2))
set.seed(1993)
widow_large <- sdmTMB(large ~ 0 + spice_iso26 +
                          s(jday_scaled, k = 3) +
                          s(sst_scaled, k = 3) +
                          s(sss_scaled, k = 3),
                        spatial_varying = ~ 0 + spice_iso26,
                        extra_time = extra_years,
                        data = yoy_widow,
                        mesh = widow_mesh,
                        spatial = "on",
                        time = "year",
                        family = tweedie(link = "log"),
                        spatiotemporal = "off",
                        control = sdmTMBcontrol(newton_loops = 1,
                                                nlminb_loops = 2))

# Shortbelly Rockfish ----
yoy_shortbelly <- read_data('yoy_sbly.Rdata')
shortbelly_mesh <- make_mesh(data,
                        xy_cols = c("X", "Y"),
                        cutoff = 15)

set.seed(1993)
shortbelly_small <- sdmTMB(small ~ 0 + depth_iso26 +
                        s(jday_scaled, k = 3) +
                        s(sst_scaled, k = 3) +
                        s(sss_scaled, k = 3),
                      spatial_varying = ~ 0 + depth_iso26,
                      extra_time = extra_years,
                      data = yoy_shortbelly,
                      mesh = shortbelly_mesh,
                      spatial = "on",
                      time = "year",
                      family = tweedie(link = "log"),
                      spatiotemporal = "off",
                      control = sdmTMBcontrol(newton_loops = 1,
                                              nlminb_loops = 2))
set.seed(1993)
shortbelly_large <- sdmTMB(large ~ 0 + vgeo +
                        s(jday_scaled, k = 3) +
                        s(sst_scaled, k = 3) +
                        s(sss_scaled, k = 3),
                      spatial_varying = ~ 0 + vgeo,
                      extra_time = extra_years,
                      data = yoy_shortbelly,
                      mesh = shortbelly_mesh,
                      spatial = "on",
                      time = "year",
                      family = tweedie(link = "log"),
                      spatiotemporal = "off",
                      control = sdmTMBcontrol(newton_loops = 1,
                                              nlminb_loops = 2))

# Market Squid ----
yoy_squid <- read_data('yoy_squid.Rdata')
squid_mesh <- make_mesh(data,
                        xy_cols = c("X", "Y"),
                        cutoff = 15)

set.seed(1993)
squid_small <- sdmTMB(small ~ 0 + depth_iso26 +
                        s(jday_scaled, k = 3) +
                        s(sst_scaled, k = 3) +
                        s(sss_scaled, k = 3),
                      spatial_varying = ~ 0 + depth_iso26,
                      extra_time = extra_years,
                      data = yoy_squid,
                      mesh = squid_mesh,
                      spatial = "on",
                      time = "year",
                      family = tweedie(link = "log"),
                      spatiotemporal = "off",
                      control = sdmTMBcontrol(newton_loops = 1,
                                              nlminb_loops = 2))
set.seed(1993)
squid_large <- sdmTMB(large ~ 0 + spice_iso26 +
                        s(jday_scaled, k = 3) +
                        s(sst_scaled, k = 3) +
                        s(sss_scaled, k = 3),
                      spatial_varying = ~ 0 + spice_iso26,
                      extra_time = extra_years,
                      data = yoy_squid,
                      mesh = squid_mesh,
                      spatial = "on",
                      time = "year",
                      family = tweedie(link = "log"),
                      spatiotemporal = "off",
                      control = sdmTMBcontrol(newton_loops = 1,
                                              nlminb_loops = 2))