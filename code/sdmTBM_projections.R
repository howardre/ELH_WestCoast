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
hindcast_ss <- readRDS(here('data', 'nep_final.Rdata'))
hindcast_means <- readRDS(here('data', 'nep_avgs.Rdata'))

state_labels <- data.frame(name = c("Washington", "Oregon", "California"),
                           lat = c(47, 44.0, 37.0),
                           lon = c(-121.0, -121.0, -120.0))

extra_years <- c(2020:2100)

### Pacific Hake ---------------------------------------------------------------------------------------------------------------------------------
yoy_hake <- read_data('yoy_hake.Rdata') 
hake_mesh <- make_mesh(yoy_hake, 
                       xy_cols = c("X", "Y"),
                       cutoff = 15)

set.seed(1993)
hake_small <- sdmTMB(small ~ 0 + depth_iso26 +
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

#### Hindcast--------------------------------------------------------------------------------------------------------------------------------------
hake_hindcast <- sdm_cells(yoy_hake, hake_small, hake_large,
                           hindcast_means, hindcast_ss, 1995:2019)

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(hake_hindcast[[1]], "Small (15-35 mm)", "Latitude \u00B0N")
map_project(hake_hindcast[[2]], "Large (36-81 mm)", "")
dev.copy(jpeg, here('results/forecast_output/yoy_hake', 
                    'hake_hindcast_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(hake_hindcast[[3]]))

#### IPSL------------------------------------------------------------------------------------------------------------------------------------------
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
dev.copy(jpeg, here('results/forecast_output/yoy_hake', 
                    'hake_ipsl1_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(hake_ipsl1[[3]]))

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

# Overlap
mean(as.numeric(hake_ipsl2[[3]]))

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

# Overlap
mean(as.numeric(hake_ipsl3[[3]]))

# Remove objects
rm(hake_hindcast, hake_ipsl1, hake_ipsl2,
   hake_ipsl3, yoy_hake, hake_mesh, hake_large,
   hake_small)


# Northern Anchovy --------------------------------------------------------------------------------------------------------------------------------
yoy_anchovy <- filter(read_data('yoy_anch.Rdata'), latitude < 42, year > 2013)
anchovy_mesh <- make_mesh(yoy_anchovy,
                          xy_cols = c("X", "Y"),
                          cutoff = 15)

set.seed(1993)
anchovy_small <- sdmTMB(small ~ 0 + u_vint_100m +
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

#### Hindcast--------------------------------------------------------------------------------------------------------------------------------------
anchovy_hindcast <- sdm_cells(yoy_anchovy, anchovy_small, anchovy_large,
                              hindcast_means, hindcast_ss, 2014:2019)

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(anchovy_hindcast[[1]], "Small (15-35 mm)", "Latitude \u00B0N")
map_project(anchovy_hindcast[[2]], "Large (36-85 mm)", "")
dev.copy(jpeg, here('results/forecast_output/yoy_anchovy', 
                    'anchovy_hindcast_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(anchovy_hindcast[[3]]))

#### IPSL------------------------------------------------------------------------------------------------------------------------------------------
##### 2020-2040------------------------------------------------------------------------------------------------------------------------------------
anchovy_ipsl1 <- sdm_cells(yoy_anchovy, anchovy_small, anchovy_large,
                        roms_means, roms_ss, 2020:2040)
saveRDS(anchovy_ipsl1, file = here("data", "anchovy_ipsl1.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(anchovy_ipsl1[[1]], "Small (15-35 mm)", "Latitude \u00B0N")
map_project(anchovy_ipsl1[[2]], "Large (36-85 mm)", "")
dev.copy(jpeg, here('results/forecast_output/yoy_anchovy', 
                    'anchovy_ipsl1_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(anchovy_ipsl1[[3]]))

##### 2050-2070------------------------------------------------------------------------------------------------------------------------------------
anchovy_ipsl2 <- sdm_cells(yoy_anchovy, anchovy_small, anchovy_large,
                        roms_means, roms_ss, 2050:2070)
saveRDS(anchovy_ipsl2, file = here("data", "anchovy_ipsl2.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(anchovy_ipsl2[[1]], "Small (15-35 mm)", "Latitude \u00B0N")
map_project(anchovy_ipsl2[[2]], "Large (36-85 mm)", "")
dev.copy(jpeg, here('results/forecast_output/yoy_anchovy', 
                    'anchovy_ipsl2_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(anchovy_ipsl2[[3]]))

##### 2080-2100------------------------------------------------------------------------------------------------------------------------------------
anchovy_ipsl3 <- sdm_cells(yoy_anchovy, anchovy_small, anchovy_large,
                        roms_means, roms_ss, 2080:2100)
saveRDS(anchovy_ipsl3, file = here("data", "anchovy_ipsl3.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(anchovy_ipsl3[[1]], "Small (15-35 mm)", "Latitude \u00B0N")
map_project(anchovy_ipsl3[[2]], "Large (36-81 mm)", "")
dev.copy(jpeg, here('results/forecast_output/yoy_anchovy', 
                    'anchovy_ipsl3_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(anchovy_ipsl3[[3]]))

# Remove objects
rm(anchovy_hindcast, anchovy_ipsl1, anchovy_ipsl2,
   anchovy_ipsl3, yoy_anchovy, anchovy_mesh, anchovy_large,
   anchovy_small)


# Pacific Sanddab ---------------------------------------------------------------------------------------------------------------------------------
yoy_sdab <- filter(read_data('yoy_dab.Rdata'), latitude > 36, year > 2012)
sdab_mesh <- make_mesh(yoy_sdab,
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
sdab_large <- sdmTMB(large ~ 0 + u_vint_50m +
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

#### Hindcast--------------------------------------------------------------------------------------------------------------------------------------
sdab_hindcast <- sdm_cells(yoy_sdab, sdab_small, sdab_large,
                              hindcast_means, hindcast_ss, 2013:2019)

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(sdab_hindcast[[1]], "Small (16-25 mm)", "Latitude \u00B0N")
map_project(sdab_hindcast[[2]], "Large (26-55 mm)", "")
dev.copy(jpeg, here('results/forecast_output/yoy_sanddab', 
                    'sdab_hindcast_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(sdab_hindcast[[3]]))

#### IPSL------------------------------------------------------------------------------------------------------------------------------------------
##### 2020-2040------------------------------------------------------------------------------------------------------------------------------------
sdab_ipsl1 <- sdm_cells(yoy_sdab, sdab_small, sdab_large,
                           roms_means, roms_ss, 2020:2040)
saveRDS(sdab_ipsl1, file = here("data", "sdab_ipsl1.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(sdab_ipsl1[[1]], "Small (16-25 mm)", "Latitude \u00B0N")
map_project(sdab_ipsl1[[2]], "Large (26-55 mm)", "")
dev.copy(jpeg, here('results/forecast_output/yoy_sanddab', 
                    'sdab_ipsl1_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(sdab_ipsl1[[3]]))

##### 2050-2070------------------------------------------------------------------------------------------------------------------------------------
sdab_ipsl2 <- sdm_cells(yoy_sdab, sdab_small, sdab_large,
                           roms_means, roms_ss, 2050:2070)
saveRDS(sdab_ipsl2, file = here("data", "sdab_ipsl2.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(sdab_ipsl2[[1]], "Small (16-25 mm)", "Latitude \u00B0N")
map_project(sdab_ipsl2[[2]], "Large (26-55 mm)", "")
dev.copy(jpeg, here('results/forecast_output/yoy_sanddab', 
                    'sdab_ipsl2_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(sdab_ipsl2[[3]]))

##### 2080-2100------------------------------------------------------------------------------------------------------------------------------------
sdab_ipsl3 <- sdm_cells(yoy_sdab, sdab_small, sdab_large,
                           roms_means, roms_ss, 2080:2100)
saveRDS(sdab_ipsl3, file = here("data", "sdab_ipsl3.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(sdab_ipsl3[[1]], "Small (16-25 mm)", "Latitude \u00B0N")
map_project(sdab_ipsl3[[2]], "Large (26-55 mm)", "")
dev.copy(jpeg, here('results/forecast_output/yoy_sanddab', 
                    'sdab_ipsl3_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(sdab_ipsl3[[3]]))
<<<<<<< HEAD

# Remove objects
rm(sdab_hindcast, sdab_ipsl1, sdab_ipsl2,
   sdab_ipsl3, yoy_sdab, sdab_mesh, sdab_large,
   sdab_small)
=======
>>>>>>> d9a5ea490fddc154f4fbd7af6e55e5274449e9dd


# Shortbelly Rockfish ---------------------------------------------------------------------------------------------------------------------------------
yoy_shortbelly <- read_data('yoy_sbly.Rdata')
shortbelly_mesh <- make_mesh(yoy_shortbelly,
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

#### Hindcast--------------------------------------------------------------------------------------------------------------------------------------
shortbelly_hindcast <- sdm_cells(yoy_shortbelly, shortbelly_small, shortbelly_large,
                                 hindcast_means, hindcast_ss, 1995:2019)

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(shortbelly_hindcast[[1]], "Small (11-35 mm)", "Latitude \u00B0N")
map_project(shortbelly_hindcast[[2]], "Large (36-78 mm)", "")
dev.copy(jpeg, here('results/forecast_output/yoy_shortbelly', 
                    'shortbelly_hindcast_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(shortbelly_hindcast[[3]]))

#### IPSL------------------------------------------------------------------------------------------------------------------------------------------
##### 2020-2040------------------------------------------------------------------------------------------------------------------------------------
shortbelly_ipsl1 <- sdm_cells(yoy_shortbelly, shortbelly_small, shortbelly_large,
                              roms_means, roms_ss, 2020:2040)
saveRDS(shortbelly_ipsl1, file = here("data", "shortbelly_ipsl1.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(shortbelly_ipsl1[[1]], "Small (11-35 mm)", "Latitude \u00B0N")
map_project(shortbelly_ipsl1[[2]], "Large (36-78 mm)", "")
dev.copy(jpeg, here('results/forecast_output/yoy_shortbelly', 
                    'shortbelly_ipsl1_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(shortbelly_ipsl1[[3]]))

##### 2050-2070------------------------------------------------------------------------------------------------------------------------------------
shortbelly_ipsl2 <- sdm_cells(yoy_shortbelly, shortbelly_small, shortbelly_large,
                              roms_means, roms_ss, 2050:2070)
saveRDS(shortbelly_ipsl2, file = here("data", "shortbelly_ipsl2.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(shortbelly_ipsl2[[1]], "Small (11-35 mm)", "Latitude \u00B0N")
map_project(shortbelly_ipsl2[[2]], "Large (36-78 mm)", "")
dev.copy(jpeg, here('results/forecast_output/yoy_shortbelly', 
                    'shortbelly_ipsl2_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(shortbelly_ipsl2[[3]]))

##### 2080-2100------------------------------------------------------------------------------------------------------------------------------------
shortbelly_ipsl3 <- sdm_cells(yoy_shortbelly, shortbelly_small, shortbelly_large,
                              roms_means, roms_ss, 2080:2100)
saveRDS(shortbelly_ipsl3, file = here("data", "shortbelly_ipsl3.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(shortbelly_ipsl3[[1]], "Small (11-35 mm)", "Latitude \u00B0N")
map_project(shortbelly_ipsl3[[2]], "Large (36-78 mm)", "")
dev.copy(jpeg, here('results/forecast_output/yoy_shortbelly', 
                    'shortbelly_ipsl3_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(shortbelly_ipsl3[[3]]))

<<<<<<< HEAD
# Remove objects
rm(shortbelly_hindcast, shortbelly_ipsl1, shortbelly_ipsl2,
   shortbelly_ipsl3, yoy_shortbelly, shortbelly_mesh, shortbelly_large,
   shortbelly_small)

=======
>>>>>>> d9a5ea490fddc154f4fbd7af6e55e5274449e9dd

# Widow Rockfish ---------------------------------------------------------------------------------------------------------------------------------
yoy_widow <- read_data('yoy_widw.Rdata') 
widow_mesh <- make_mesh(yoy_widow,
                        xy_cols = c("X", "Y"),
                        cutoff = 15)

set.seed(1993)
widow_small <- sdmTMB(small ~ 0 + vmax_cu +
                          s(jday_scaled, k = 3) +
                          s(sst_scaled, k = 3) +
                          s(sss_scaled, k = 3),
                        spatial_varying = ~ 0 + vmax_cu,
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

#### Hindcast--------------------------------------------------------------------------------------------------------------------------------------
widow_hindcast <- sdm_cells(yoy_widow, widow_small, widow_large,
                           hindcast_means, hindcast_ss, 1995:2019)

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(widow_hindcast[[1]], "Small (17-32 mm)", "Latitude \u00B0N")
map_project(widow_hindcast[[2]], "Large (33-64 mm)", "")
dev.copy(jpeg, here('results/forecast_output/yoy_widow', 
                    'widow_hindcast_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(widow_hindcast[[3]]))

#### IPSL------------------------------------------------------------------------------------------------------------------------------------------
##### 2020-2040------------------------------------------------------------------------------------------------------------------------------------
widow_ipsl1 <- sdm_cells(yoy_widow, widow_small, widow_large,
                        roms_means, roms_ss, 2020:2040)
saveRDS(widow_ipsl1, file = here("data", "widow_ipsl1.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(widow_ipsl1[[1]], "Small (17-32 mm)", "Latitude \u00B0N")
map_project(widow_ipsl1[[2]], "Large (33-64 mm)", "")
dev.copy(jpeg, here('results/forecast_output/yoy_widow', 
                    'widow_ipsl1_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(widow_ipsl1[[3]]))

##### 2050-2070------------------------------------------------------------------------------------------------------------------------------------
widow_ipsl2 <- sdm_cells(yoy_widow, widow_small, widow_large,
                        roms_means, roms_ss, 2050:2070)
saveRDS(widow_ipsl2, file = here("data", "widow_ipsl2.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(widow_ipsl2[[1]], "Small (17-32 mm)", "Latitude \u00B0N")
map_project(widow_ipsl2[[2]], "Large (33-64 mm)", "")
dev.copy(jpeg, here('results/forecast_output/yoy_widow', 
                    'widow_ipsl2_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(widow_ipsl2[[3]]))

##### 2080-2100------------------------------------------------------------------------------------------------------------------------------------
widow_ipsl3 <- sdm_cells(yoy_widow, widow_small, widow_large,
                        roms_means, roms_ss, 2080:2100)
saveRDS(widow_ipsl3, file = here("data", "widow_ipsl3.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(widow_ipsl3[[1]], "Small (17-32 mm)", "Latitude \u00B0N")
map_project(widow_ipsl3[[2]], "Large (33-64 mm)", "")
dev.copy(jpeg, here('results/forecast_output/yoy_widow', 
                    'widow_ipsl3_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(widow_ipsl3[[3]]))

<<<<<<< HEAD
# Remove objects
rm(widow_hindcast, widow_ipsl1, widow_ipsl2,
   widow_ipsl3, yoy_widow, widow_mesh, widow_large,
   widow_small)

=======
>>>>>>> d9a5ea490fddc154f4fbd7af6e55e5274449e9dd

# Market Squid ------------------------------------------------------------------------------------------------------------------------------------
yoy_squid <- read_data('yoy_squid.Rdata')
squid_mesh <- make_mesh(yoy_squid,
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

#### Hindcast--------------------------------------------------------------------------------------------------------------------------------------
squid_hindcast <- sdm_cells(yoy_squid, squid_small, squid_large,
                                 hindcast_means, hindcast_ss, 2004:2019)

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(squid_hindcast[[1]], "Small (11-70 mm)", "Latitude \u00B0N")
map_project(squid_hindcast[[2]], "Large (71-139 mm)", "")
dev.copy(jpeg, here('results/forecast_output/yoy_squid', 
                    'squid_hindcast_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(squid_hindcast[[3]]))

#### IPSL------------------------------------------------------------------------------------------------------------------------------------------
##### 2020-2040------------------------------------------------------------------------------------------------------------------------------------
squid_ipsl1 <- sdm_cells(yoy_squid, squid_small, squid_large,
                              roms_means, roms_ss, 2020:2040)
saveRDS(squid_ipsl1, file = here("data", "squid_ipsl1.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(squid_ipsl1[[1]], "Small (11-70 mm)", "Latitude \u00B0N")
map_project(squid_ipsl1[[2]], "Large (71-139 mm)", "")
dev.copy(jpeg, here('results/forecast_output/yoy_squid', 
                    'squid_ipsl1_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(squid_ipsl1[[3]]))

##### 2050-2070------------------------------------------------------------------------------------------------------------------------------------
squid_ipsl2 <- sdm_cells(yoy_squid, squid_small, squid_large,
                              roms_means, roms_ss, 2050:2070)
saveRDS(squid_ipsl2, file = here("data", "squid_ipsl2.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(squid_ipsl2[[1]], "Small (11-70 mm)", "Latitude \u00B0N")
map_project(squid_ipsl2[[2]], "Large (71-139 mm)", "")
dev.copy(jpeg, here('results/forecast_output/yoy_squid', 
                    'squid_ipsl2_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(squid_ipsl2[[3]]))

##### 2080-2100------------------------------------------------------------------------------------------------------------------------------------
squid_ipsl3 <- sdm_cells(yoy_squid, squid_small, squid_large,
                              roms_means, roms_ss, 2080:2100)
saveRDS(squid_ipsl3, file = here("data", "squid_ipsl3.rds"))

# Plot
windows(height = 15, width = 18)
par(mfrow = c(1, 2),
    mar = c(6.6, 7.6, 3.5, 0.6) + 0.1,
    oma = c(1, 1, 1, 1),
    mgp = c(5, 2, 0),
    family = "serif")
map_project(squid_ipsl3[[1]], "Small (11-70 mm)", "Latitude \u00B0N")
map_project(squid_ipsl3[[2]], "Large (71-139 mm)", "")
dev.copy(jpeg, here('results/forecast_output/yoy_squid', 
                    'squid_ipsl3_plot.jpg'), 
         height = 15, 
         width = 16, 
         units = 'in', 
         res = 200)
dev.off()

# Overlap
mean(as.numeric(squid_ipsl3[[3]]))
<<<<<<< HEAD

# Remove objects
rm(squid_hindcast, squid_ipsl1, squid_ipsl2,
   squid_ipsl3, yoy_squid, squid_mesh, squid_large,
   squid_small)
=======
>>>>>>> d9a5ea490fddc154f4fbd7af6e55e5274449e9dd
